#pragma once

#include <complex>
#include <cmath>
#include <cstddef>
#include <functional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// These globals are used by existing headers but are not declared there.
// Declare them here to make this header self-contained.
extern bool use_pseudo_quadratic_element;
extern bool use_true_quadratic_element;
extern double simulation_time;
extern double coupling_tol;

#include "BEM_setBoundaryTypes.hpp"
#include "BEM_solveBVP.hpp"

// -----------------------------------------------------------------------------
// Linear Frequency-Domain BEM (Rankine source, boundary-explicit)
//
// Scope:
// - Solves Laplace BVP with mixed boundary conditions on a triangulated boundary.
// - Supports: Neumann (given phi_n), Dirichlet (given phi), Robin on free-surface:
//     phi_n = kappa * phi, where kappa may be complex (sponge zone).
//
// Design choice (important for corners):
// - Robin faces are treated as "Dirichlet-type" for indexing (unknown is phi_n),
//   so that intersections with Neumann faces become CORNER points and get the same
//   continuity constraint handling as the time-domain code.
// -----------------------------------------------------------------------------

namespace bem_frequency_domain {

using Complex = std::complex<double>;
using Id = std::tuple<networkPoint *, networkFace *>;

enum class FaceBC {
  Neumann,   // phi_n is prescribed, unknown is phi (standard)
  Dirichlet, // phi is prescribed, unknown is phi_n (standard)
  Robin      // phi_n = kappa * phi (free-surface); handled as Dirichlet-type (unknown phi_n)
};

struct LinearFSBC {
  double omega = 0.0;   // [rad/s]
  double gravity = 9.81;
  // Sponge coefficient mu(x) >= 0 (0 means no sponge). If null, treated as 0.
  std::function<double(const networkPoint &)> sponge_mu;

  Complex kappa_at(const networkPoint &p) const {
    const double mu = sponge_mu ? sponge_mu(p) : 0.0;
    const Complex s(mu, -omega); // (-i*omega + mu)
    if (gravity == 0.0)
      throw std::runtime_error("LinearFSBC: gravity must be non-zero");
    // Linear free-surface: (d/dt + mu)^2 phi + g * phi_n = 0 -> phi_n = -(s^2/g) * phi.
    return -(s * s) / gravity;
  }
};

struct BoundaryData {
  // Face classification.
  std::function<FaceBC(const networkFace &)> face_bc;

  // Neumann value: prescribed phi_n on Neumann faces (may depend on face; face can be nullptr for single-ID nodes).
  std::function<Complex(const networkPoint &, const networkFace *)> neumann_phin;

  // Dirichlet value: prescribed phi on Dirichlet faces (phi is pointwise unique in this codebase).
  std::function<Complex(const networkPoint &)> dirichlet_phi;
};

struct Solution {
  std::size_t n = 0;
  std::vector<Id> id_by_index;

  // Unknown vector u (size n): for Dirichlet-type IDs -> phi_n, for Neumann IDs -> phi.
  std::vector<Complex> u;

  // Reconstructed boundary values (size n, aligned with id_by_index).
  std::vector<Complex> phi;
  std::vector<Complex> phin;

  // Residual of the full BIE (not the modified corner rows).
  std::vector<Complex> bie_residual;
  double bie_residual_l2 = 0.0;
};

// Fortran LAPACK (complex double). Assumes symbols with trailing underscore.
extern "C" void zgetrf_(const int *m, const int *n, Complex *a, const int *lda, int *ipiv, int *info);
extern "C" void zgetrs_(const char *trans, const int *n, const int *nrhs, const Complex *a, const int *lda, const int *ipiv, Complex *b, const int *ldb, int *info);

struct lapack_zlu {
  int n = 0;
  int lda = 0;
  int info = 0;
  std::vector<int> ipiv;
  std::vector<Complex> a_col_major;

  lapack_zlu() = default;

  lapack_zlu(int n_in, std::vector<Complex> a_in_col_major) : n(n_in), lda(n_in), ipiv(static_cast<std::size_t>(n_in)), a_col_major(std::move(a_in_col_major)) {
    if (n <= 0)
      throw std::runtime_error("lapack_zlu: n must be > 0");
    if (static_cast<int>(a_col_major.size()) != n * n)
      throw std::runtime_error("lapack_zlu: matrix size mismatch");
    zgetrf_(&n, &n, a_col_major.data(), &lda, ipiv.data(), &info);
    if (info != 0)
      throw std::runtime_error("lapack_zlu: zgetrf_ failed (info=" + std::to_string(info) + ")");
  }

  void solve_in_place(std::vector<Complex> &b) const {
    if (static_cast<int>(b.size()) != n)
      throw std::runtime_error("lapack_zlu::solve_in_place: RHS size mismatch");
    int nrhs = 1;
    int ldb = n;
    int info_solve = 0;
    const char trans = 'N';
    zgetrs_(&trans, &n, &nrhs, a_col_major.data(), &lda, ipiv.data(), b.data(), &ldb, &info_solve);
    if (info_solve != 0)
      throw std::runtime_error("lapack_zlu: zgetrs_ failed (info=" + std::to_string(info_solve) + ")");
  }
};

inline double l2_norm(const std::vector<Complex> &v) {
  long double sum = 0.0L;
  for (const auto &x : v) {
    const long double a = static_cast<long double>(std::abs(x));
    sum += a * a;
  }
  return std::sqrt(static_cast<double>(sum));
}

struct BoundaryStateGuard {
  struct FaceState {
    networkFace *f = nullptr;
    bool Dirichlet = false;
    bool Neumann = false;
    bool isLinearElement = true;
    bool isPseudoQuadraticElement = false;
    bool isTrueQuadraticElement = false;
  };
  struct LineState {
    networkLine *l = nullptr;
    bool Dirichlet = false;
    bool Neumann = false;
    bool CORNER = false;
  };
  struct PointState {
    networkPoint *p = nullptr;
    bool Dirichlet = false;
    bool Neumann = false;
    bool CORNER = false;
    bool isMultipleNode = false;
    std::unordered_map<networkFace *, int> f2Index;
  };

  std::vector<FaceState> faces;
  std::vector<LineState> lines;
  std::vector<PointState> points;

  explicit BoundaryStateGuard(const std::vector<Network *> &nets) {
    std::unordered_set<networkFace *> uniq_faces;
    std::unordered_set<networkLine *> uniq_lines;
    std::unordered_set<networkPoint *> uniq_points;
    for (auto *net : nets) {
      if (!net)
        continue;
      for (auto *f : net->getBoundaryFaces())
        uniq_faces.emplace(f);
      for (auto *l : net->getLines())
        uniq_lines.emplace(l);
      for (auto *p : net->getPoints())
        uniq_points.emplace(p);
    }

    faces.reserve(uniq_faces.size());
    for (auto *f : uniq_faces) {
      faces.push_back(FaceState{f, f->Dirichlet, f->Neumann, f->isLinearElement, f->isPseudoQuadraticElement, f->isTrueQuadraticElement});
    }
    lines.reserve(uniq_lines.size());
    for (auto *l : uniq_lines) {
      lines.push_back(LineState{l, l->Dirichlet, l->Neumann, l->CORNER});
    }
    points.reserve(uniq_points.size());
    for (auto *p : uniq_points) {
      points.push_back(PointState{p, p->Dirichlet, p->Neumann, p->CORNER, p->isMultipleNode, p->f2Index});
    }
  }

  void restore() const {
    for (const auto &s : faces) {
      if (!s.f)
        continue;
      s.f->Dirichlet = s.Dirichlet;
      s.f->Neumann = s.Neumann;
      s.f->isLinearElement = s.isLinearElement;
      s.f->isPseudoQuadraticElement = s.isPseudoQuadraticElement;
      s.f->isTrueQuadraticElement = s.isTrueQuadraticElement;
    }
    for (const auto &s : lines) {
      if (!s.l)
        continue;
      s.l->Dirichlet = s.Dirichlet;
      s.l->Neumann = s.Neumann;
      s.l->CORNER = s.CORNER;
    }
    for (const auto &s : points) {
      if (!s.p)
        continue;
      s.p->Dirichlet = s.Dirichlet;
      s.p->Neumann = s.Neumann;
      s.p->CORNER = s.CORNER;
      s.p->isMultipleNode = s.isMultipleNode;
      s.p->f2Index = s.f2Index;
    }
  }
};

inline void apply_face_bc_to_mesh(const std::vector<Network *> &nets, const BoundaryData &bc, std::unordered_set<networkFace *> &robin_faces_out) {
  robin_faces_out.clear();
  if (!bc.face_bc)
    throw std::runtime_error("apply_face_bc_to_mesh: face_bc is not set");

  for (auto *net : nets) {
    if (!net)
      continue;
    net->setGeometricPropertiesForce();
    for (auto *f : net->getBoundaryFaces()) {
      const FaceBC t = bc.face_bc(*f);
      if (t == FaceBC::Neumann) {
        f->Neumann = true;
        f->Dirichlet = false;
      } else {
        // Dirichlet or Robin are "Dirichlet-type" for indexing (unknown is phi_n).
        f->Dirichlet = true;
        f->Neumann = false;
        if (t == FaceBC::Robin)
          robin_faces_out.emplace(f);
      }
    }

    for (auto *l : net->getLines()) {
      const auto &fs = l->getBoundaryFaces();
      l->Neumann = std::ranges::all_of(fs, [](const auto *f) { return f->Neumann; });
      l->Dirichlet = std::ranges::all_of(fs, [](const auto *f) { return f->Dirichlet; });
      l->CORNER = (!l->Neumann && !l->Dirichlet);
    }

    for (auto *p : net->getPoints()) {
      const auto &fs = p->getBoundaryFaces();
      p->Neumann = std::ranges::all_of(fs, [](const auto *f) { return f->Neumann; });
      p->Dirichlet = std::ranges::all_of(fs, [](const auto *f) { return f->Dirichlet; });
      p->CORNER = (!p->Neumann && !p->Dirichlet);
      setIsMultipleNode(p);
    }

    for (auto *f : net->getBoundaryFaces()) {
      if (use_true_quadratic_element) {
        f->isTrueQuadraticElement = true;
        f->isPseudoQuadraticElement = false;
        f->isLinearElement = false;
      } else {
        f->isTrueQuadraticElement = false;
        f->isPseudoQuadraticElement = use_pseudo_quadratic_element && f->Dirichlet;
        f->isLinearElement = !f->isPseudoQuadraticElement;
      }
    }
  }
}

inline Solution solve_linear_bvp(const std::vector<Network *> &nets, const LinearFSBC &fsbc, const BoundaryData &bc) {
  if (nets.empty())
    throw std::runtime_error("solve_linear_bvp: empty networks");

  BoundaryStateGuard guard(nets);

  std::unordered_set<networkFace *> robin_faces;
  apply_face_bc_to_mesh(nets, bc, robin_faces);

  // Build ID indices (this fills p->f2Index).
  const std::size_t n = setNodeFaceIndices(nets);
  if (n == 0)
    throw std::runtime_error("solve_linear_bvp: matrix size is 0");

  // Prepare BIE coefficients (dense, geometry-only).
  BEM_BVP bvp(nets);
  bvp.matrix_size = static_cast<int>(n);
  bvp.setIGIGn();

  std::vector<Id> id_by_index(n);
  for (auto *net : nets) {
    for (auto *p : net->getBoundaryPoints()) {
      for (const auto &[f, i] : p->f2Index) {
        if (i < 0 || static_cast<std::size_t>(i) >= n)
          continue;
        id_by_index[static_cast<std::size_t>(i)] = {p, f};
      }
    }
  }

  // Classify Dirichlet IDs into Robin-vs-true-Dirichlet using point adjacency.
  std::unordered_map<networkPoint *, bool> point_is_robin;
  point_is_robin.reserve(n);
  for (const auto &[p, f] : id_by_index) {
    (void)f;
    if (!p)
      continue;
    if (point_is_robin.contains(p))
      continue;
    bool is_robin = false;
    for (auto *adj : p->getBoundaryFaces()) {
      if (robin_faces.contains(adj)) {
        is_robin = true;
        break;
      }
    }
    point_is_robin.emplace(p, is_robin);
  }

  // Scaling for constraint rows (mimic generateBIEMatrix() idea).
  double max_value = 1.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double d = std::abs(bvp.IGIGn[i][i][0]);
    if (d > max_value)
      max_value = d;
  }

  // Assemble A*u=b (complex), column-major.
  std::vector<Complex> A(n * n, Complex{0.0, 0.0});
  std::vector<Complex> b(n, Complex{0.0, 0.0});

  auto aidx = [&](std::size_t row, std::size_t col) -> std::size_t { return row + col * n; };

  auto get_dirichlet_phi = [&](const networkPoint *p) -> Complex {
    if (!bc.dirichlet_phi)
      return Complex{0.0, 0.0};
    return bc.dirichlet_phi(*p);
  };
  auto get_neumann_phin = [&](const networkPoint *p, const networkFace *f) -> Complex {
    if (!bc.neumann_phin)
      return Complex{0.0, 0.0};
    return bc.neumann_phin(*p, f);
  };

  // Main assembly loop.
  for (std::size_t i = 0; i < n; ++i) {
    auto [p_row, f_row] = id_by_index[i];
    if (!p_row)
      throw std::runtime_error("solve_linear_bvp: id_by_index has null point");

    const bool row_is_corner_neumann = (p_row->CORNER && isNeumannID_BEM(p_row, f_row));
    if (row_is_corner_neumann) {
      // Replace BIE row by continuity constraint: phi(neumann-id) == phi(dirichlet-id).
      A[aidx(i, i)] = max_value;

      const int j_dir_i = pf2Index(p_row, nullptr);
      if (j_dir_i < 0 || static_cast<std::size_t>(j_dir_i) >= n)
        throw std::runtime_error("solve_linear_bvp: missing Dirichlet ID at CORNER point");
      const std::size_t j_dir = static_cast<std::size_t>(j_dir_i);

      const bool robin = point_is_robin[p_row];
      if (robin) {
        const Complex kappa = fsbc.kappa_at(*p_row);
        if (std::abs(kappa) < 1e-15)
          throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin CORNER point");
        A[aidx(i, j_dir)] += -max_value / kappa; // phi_dir = phin/kappa
        b[i] = Complex{0.0, 0.0};
      } else {
        b[i] = max_value * get_dirichlet_phi(p_row);
      }
      continue;
    }

    // Standard BIE row: sum IG*phin = sum IGn*phi (with substitutions for Robin).
    for (std::size_t j = 0; j < n; ++j) {
      const double IG = bvp.IGIGn[i][j][0];
      const double IGn = bvp.IGIGn[i][j][1];

      auto [p_col, f_col] = id_by_index[j];
      if (!p_col)
        throw std::runtime_error("solve_linear_bvp: id_by_index has null point (col)");

      if (isDirichletID_BEM(p_col, f_col)) {
        // Unknown is phi_n (u_j).
        const bool robin = point_is_robin[p_col];
        if (robin) {
          const Complex kappa = fsbc.kappa_at(*p_col);
          if (std::abs(kappa) < 1e-15)
            throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin point");
          A[aidx(i, j)] += Complex{IG, 0.0} - Complex{IGn, 0.0} / kappa; // phi = phin/kappa
        } else {
          A[aidx(i, j)] += Complex{IG, 0.0};
          b[i] += Complex{IGn, 0.0} * get_dirichlet_phi(p_col);
        }
      } else {
        // Neumann ID: unknown is phi (u_j), phi_n is known.
        A[aidx(i, j)] += Complex{-IGn, 0.0};
        b[i] += -Complex{IG, 0.0} * get_neumann_phin(p_col, f_col);
      }
    }
  }

  // Solve
  lapack_zlu lu(static_cast<int>(n), A);
  lu.solve_in_place(b); // b becomes the solution u

  Solution sol;
  sol.n = n;
  sol.id_by_index = std::move(id_by_index);
  sol.u = b;
  sol.phi.assign(n, Complex{0.0, 0.0});
  sol.phin.assign(n, Complex{0.0, 0.0});

  // Reconstruct (phi, phin) on each ID.
  for (std::size_t j = 0; j < n; ++j) {
    auto [p_col, f_col] = sol.id_by_index[j];
    if (isDirichletID_BEM(p_col, f_col)) {
      sol.phin[j] = sol.u[j];
      const bool robin = point_is_robin[p_col];
      if (robin) {
        const Complex kappa = fsbc.kappa_at(*p_col);
        sol.phi[j] = sol.u[j] / kappa;
      } else {
        sol.phi[j] = get_dirichlet_phi(p_col);
      }
    } else {
      sol.phi[j] = sol.u[j];
      sol.phin[j] = get_neumann_phin(p_col, f_col);
    }
  }

  // Full BIE residual check (for diagnostics).
  sol.bie_residual.assign(n, Complex{0.0, 0.0});
  for (std::size_t i = 0; i < n; ++i) {
    Complex r{0.0, 0.0};
    for (std::size_t j = 0; j < n; ++j) {
      const double IG = bvp.IGIGn[i][j][0];
      const double IGn = bvp.IGIGn[i][j][1];
      r += Complex{IG, 0.0} * sol.phin[j] - Complex{IGn, 0.0} * sol.phi[j];
    }
    sol.bie_residual[i] = r;
  }
  sol.bie_residual_l2 = l2_norm(sol.bie_residual);

  // Restore original boundary flags/indexing.
  guard.restore();
  return sol;
}

} // namespace bem_frequency_domain
