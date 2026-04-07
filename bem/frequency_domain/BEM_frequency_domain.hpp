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
extern bool use_quadratic_linear_hybrid;
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

// DOF identifier: either a point or a line (edge midpoint), paired with a face.
struct DofId {
  networkPoint* point = nullptr;
  networkLine* line = nullptr;
  networkFace* face = nullptr;

  bool is_point() const { return point != nullptr; }
  bool is_line() const { return line != nullptr; }
};

// Legacy alias
using Id = std::tuple<networkPoint*, networkFace*>;

enum class FaceBC {
  Neumann,   // phi_n is prescribed, unknown is phi (standard)
  Dirichlet, // phi is prescribed, unknown is phi_n (standard)
  Robin      // phi_n = kappa * phi (free-surface); handled as Dirichlet-type (unknown phi_n)
};

struct LinearFSBC {
  double omega = 0.0; // [rad/s]
  double gravity = 9.81;
  // Sponge coefficient mu(x) >= 0 (0 means no sponge). If null, treated as 0.
  std::function<double(const networkPoint&)> sponge_mu;

  Complex kappa_at(const networkPoint& p) const {
    const double mu = sponge_mu ? sponge_mu(p) : 0.0;
    const Complex s(mu, -omega); // (-i*omega + mu)
    if (gravity == 0.0)
      throw std::runtime_error("LinearFSBC: gravity must be non-zero");
    // Linear free-surface: (d/dt + mu)^2 phi + g * phi_n = 0 -> phi_n = -(s^2/g) * phi.
    return -(s * s) / gravity;
  }

  // Midpoint evaluation for edge nodes (true quadratic)
  Complex kappa_at_midpoint(const networkLine& l) const {
    auto [p0, p1] = l.getPoints();
    // Evaluate sponge at midpoint coordinate
    const double mu0 = sponge_mu ? sponge_mu(*p0) : 0.0;
    const double mu1 = sponge_mu ? sponge_mu(*p1) : 0.0;
    const double mu = 0.5 * (mu0 + mu1);
    const Complex s(mu, -omega);
    if (gravity == 0.0)
      throw std::runtime_error("LinearFSBC: gravity must be non-zero");
    return -(s * s) / gravity;
  }
};

struct BoundaryData {
  // Face classification.
  std::function<FaceBC(const networkFace&)> face_bc;

  // Point (vertex node) callbacks.
  std::function<Complex(const networkPoint&, const networkFace*)> neumann_phin;
  std::function<Complex(const networkPoint&)> dirichlet_phi;

  // Line (edge node) callbacks for true_quadratic.
  // Evaluated at midpoint position with edge normal. If null, endpoint average is used.
  std::function<Complex(const networkLine&, const networkFace*)> neumann_phin_line;
  std::function<Complex(const networkLine&)> dirichlet_phi_line;
};

struct Solution {
  std::size_t n = 0;
  std::vector<DofId> id_by_index;

  // Unknown vector u (size n): for Dirichlet-type IDs -> phi_n, for Neumann IDs -> phi.
  std::vector<Complex> u;

  // Reconstructed boundary values (size n, aligned with id_by_index).
  std::vector<Complex> phi;
  std::vector<Complex> phin;

  // Residual of the full BIE (not the modified corner rows).
  std::vector<Complex> bie_residual;
  double bie_residual_l2 = 0.0;

  // Legacy access: get point/face tuple (returns nullptr for line DOFs)
  Id get_point_id(std::size_t i) const {
    const auto& d = id_by_index[i];
    return {d.point, d.face};
  }
};

// Fortran LAPACK (complex double). Assumes symbols with trailing underscore.
extern "C" void zgetrf_(const int* m, const int* n, Complex* a, const int* lda, int* ipiv, int* info);
extern "C" void zgetrs_(const char* trans, const int* n, const int* nrhs, const Complex* a, const int* lda, const int* ipiv, Complex* b, const int* ldb, int* info);

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

  void solve_in_place(std::vector<Complex>& b) const {
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

inline double l2_norm(const std::vector<Complex>& v) {
  long double sum = 0.0L;
  for (const auto& x : v) {
    const long double a = static_cast<long double>(std::abs(x));
    sum += a * a;
  }
  return std::sqrt(static_cast<double>(sum));
}

struct BoundaryStateGuard {
  struct FaceState {
    networkFace* f = nullptr;
    bool isLinearElement = true;
    bool isPseudoQuadraticElement = false;
    bool isTrueQuadraticElement = false;
    Network* penetratedBody = nullptr;
  };
  struct LineState {
    networkLine* l = nullptr;
    bool Dirichlet = false;
    bool Neumann = false;
    bool CORNER = false;
    bool isMultipleNode = false;
    std::unordered_map<networkFace*, NodeFaceState> dofs;
    int midpoint_index = -1;
  };
  struct PointState {
    networkPoint* p = nullptr;
    bool Dirichlet = false;
    bool Neumann = false;
    bool CORNER = false;
    bool isMultipleNode = false;
    std::unordered_map<networkFace*, NodeFaceState> dofs;
  };

  std::vector<FaceState> faces;
  std::vector<LineState> lines;
  std::vector<PointState> points;

  explicit BoundaryStateGuard(const std::vector<Network*>& nets) {
    std::unordered_set<networkFace*> uniq_faces;
    std::unordered_set<networkLine*> uniq_lines;
    std::unordered_set<networkPoint*> uniq_points;
    for (auto* net : nets) {
      if (!net)
        continue;
      for (auto* f : net->getBoundaryFaces())
        uniq_faces.emplace(f);
      for (auto* l : net->getLines())
        uniq_lines.emplace(l);
      for (auto* p : net->getPoints())
        uniq_points.emplace(p);
    }

    faces.reserve(uniq_faces.size());
    for (auto* f : uniq_faces)
      faces.push_back({f, f->isLinearElement, f->isPseudoQuadraticElement, f->isTrueQuadraticElement, f->penetratedBody});
    lines.reserve(uniq_lines.size());
    for (auto* l : uniq_lines)
      lines.push_back({l, l->Dirichlet, l->Neumann, l->CORNER, l->isMultipleNode, l->dofs, l->midpoint_index});
    points.reserve(uniq_points.size());
    for (auto* p : uniq_points)
      points.push_back({p, p->Dirichlet, p->Neumann, p->CORNER, p->isMultipleNode, p->dofs});
  }

  void restore() const {
    for (const auto& s : faces) {
      if (!s.f)
        continue;
      s.f->isLinearElement = s.isLinearElement;
      s.f->isPseudoQuadraticElement = s.isPseudoQuadraticElement;
      s.f->isTrueQuadraticElement = s.isTrueQuadraticElement;
      s.f->penetratedBody = s.penetratedBody;
    }
    for (const auto& s : lines) {
      if (!s.l)
        continue;
      s.l->Dirichlet = s.Dirichlet;
      s.l->Neumann = s.Neumann;
      s.l->CORNER = s.CORNER;
      s.l->isMultipleNode = s.isMultipleNode;
      s.l->dofs = s.dofs;
      s.l->midpoint_index = s.midpoint_index;
    }
    for (const auto& s : points) {
      if (!s.p)
        continue;
      s.p->Dirichlet = s.Dirichlet;
      s.p->Neumann = s.Neumann;
      s.p->CORNER = s.CORNER;
      s.p->isMultipleNode = s.isMultipleNode;
      s.p->dofs = s.dofs;
    }
  }
};

// Build face BC map and populate robin_faces_out. Returns the map for reuse.
inline std::unordered_map<const networkFace*, FaceBC>
build_face_bc_map(const std::vector<Network*>& nets, const BoundaryData& bc,
                  std::unordered_set<networkFace*>& robin_faces_out) {
  robin_faces_out.clear();
  if (!bc.face_bc)
    throw std::runtime_error("build_face_bc_map: face_bc is not set");
  std::unordered_map<const networkFace*, FaceBC> face_bc_map;
  for (auto* net : nets) {
    if (!net)
      continue;
    net->setGeometricPropertiesForce();
    for (auto* f : net->getBoundaryFaces()) {
      const FaceBC t = bc.face_bc(*f);
      face_bc_map[f] = t;
      if (t == FaceBC::Robin)
        robin_faces_out.emplace(f);
    }
  }
  return face_bc_map;
}

// Apply frequency-domain BC to the node-face state system.
// This is the single entry point for BC setup in the frequency domain.
// It populates dofs (source of truth) and derives summary flags from them.
inline void apply_frequency_domain_bc(
    const std::vector<Network*>& nets,
    const std::unordered_map<const networkFace*, FaceBC>& face_bc_map) {

  auto set_dofs_for_entity = [&](auto* entity) {
    for (auto* f : entity->getBoundaryFaces()) {
      auto it = face_bc_map.find(f);
      if (it == face_bc_map.end())
        continue;
      auto& state = entity->dof(f);
      state.detached_by_pressure = false;
      if (it->second == FaceBC::Neumann)
        state.contact_opponent_faces = {const_cast<networkFace*>(it->first)};
      else
        state.contact_opponent_faces.clear();
    }
  };

  for (auto* net : nets) {
    if (!net)
      continue;

    // Clear penetratedBody to avoid stale fallback in getNodeFaceBoundaryType
    for (auto* f : net->getBoundaryFaces())
      f->penetratedBody = nullptr;

    // Create dofs entries via shared infrastructure
    initializeNodeFaceStates(net);

    // Populate dofs contact_opponent_faces based on face BC
    for (auto* p : net->getBoundaryPoints())
      set_dofs_for_entity(p);
    if (use_true_quadratic_element) {
      for (auto* l : net->getBoundaryLines())
        set_dofs_for_entity(l);
    }

    // Derive summary flags from dofs (dofs is the source of truth)
    for (auto* l : net->getLines()) {
      const auto& fs = l->getBoundaryFaces();
      l->Neumann = std::ranges::all_of(fs, [&](const auto* f) { return isNeumannBoundaryState(l, f); });
      l->Dirichlet = std::ranges::all_of(fs, [&](const auto* f) { return isDirichletBoundaryState(l, f); });
      l->CORNER = (!l->Neumann && !l->Dirichlet);
    }
    for (auto* p : net->getPoints()) {
      const auto& fs = p->getBoundaryFaces();
      p->Neumann = std::ranges::all_of(fs, [&](const auto* f) { return isNeumannBoundaryState(p, f); });
      p->Dirichlet = std::ranges::all_of(fs, [&](const auto* f) { return isDirichletBoundaryState(p, f); });
      p->CORNER = (!p->Neumann && !p->Dirichlet);
      setMultipleNode(p);
    }
    if (use_true_quadratic_element) {
      for (auto* l : net->getBoundaryLines())
        setMultipleNode(l);
    }

    // Element type assignment
    for (auto* f : net->getBoundaryFaces()) {
      if (use_true_quadratic_element) {
        f->isTrueQuadraticElement = true;
        f->isPseudoQuadraticElement = false;
        f->isLinearElement = false;
      } else {
        f->isTrueQuadraticElement = false;
        auto it = face_bc_map.find(f);
        const bool face_is_dirichlet_type = it != face_bc_map.end() && it->second != FaceBC::Neumann;
        f->isPseudoQuadraticElement = use_pseudo_quadratic_element && face_is_dirichlet_type;
        f->isLinearElement = !f->isPseudoQuadraticElement;
      }
    }
  }
}

inline Solution solve_linear_bvp(const std::vector<Network*>& nets, const LinearFSBC& fsbc, const BoundaryData& bc) {
  if (nets.empty())
    throw std::runtime_error("solve_linear_bvp: empty networks");

  BoundaryStateGuard guard(nets);

  std::unordered_set<networkFace*> robin_faces;
  auto face_bc_map = build_face_bc_map(nets, bc, robin_faces);
  apply_frequency_domain_bc(nets, face_bc_map);

  // Build DOF indices (fills entity->dofs[face].index).
  const std::size_t n = setNodeFaceIndices(nets);
  if (n == 0)
    throw std::runtime_error("solve_linear_bvp: matrix size is 0");

  // Prepare BIE coefficients (dense, geometry-only).
  BEM_BVP bvp(nets);
  bvp.matrix_size = static_cast<int>(n);
  bvp.setIGIGn();

  std::vector<DofId> id_by_index(n);
  for (auto* net : nets) {
    for (auto* p : net->getBoundaryPoints()) {
      for (const auto& [f, state] : p->dofs) {
        if (state.index < 0 || static_cast<std::size_t>(state.index) >= n)
          continue;
        id_by_index[static_cast<std::size_t>(state.index)] = {p, nullptr, f};
      }
    }
    if (use_true_quadratic_element) {
      for (auto* l : net->getBoundaryLines()) {
        for (const auto& [f, state] : l->dofs) {
          if (state.index < 0 || static_cast<std::size_t>(state.index) >= n)
            continue;
          id_by_index[static_cast<std::size_t>(state.index)] = {nullptr, l, f};
        }
      }
    }
  }

  // Helper: check if a DOF's entity is adjacent to any Robin face
  auto is_entity_robin = [&](const DofId& d) -> bool {
    auto check = [&](const auto* entity) -> bool {
      for (auto* adj : entity->getBoundaryFaces())
        if (robin_faces.contains(adj))
          return true;
      return false;
    };
    if (d.is_point())
      return check(d.point);
    if (d.is_line())
      return check(d.line);
    return false;
  };

  // Helper: get kappa at DOF location
  auto kappa_at_dof = [&](const DofId& d) -> Complex {
    if (d.is_point())
      return fsbc.kappa_at(*d.point);
    if (d.is_line())
      return fsbc.kappa_at_midpoint(*d.line);
    throw std::runtime_error("kappa_at_dof: empty DofId");
  };

  // Helper: is this DOF Neumann or Dirichlet-type?
  auto dof_is_dirichlet_key = [](const DofId& d) -> bool {
    if (d.is_point())
      return isDirichletBieDofKey(d.point, d.face);
    if (d.is_line())
      return isDirichletBieDofKey(d.line, d.face);
    return false;
  };
  auto dof_is_neumann_key = [](const DofId& d) -> bool {
    if (d.is_point())
      return isNeumannBieDofKey(d.point, d.face);
    if (d.is_line())
      return isNeumannBieDofKey(d.line, d.face);
    return false;
  };
  auto dof_is_corner = [](const DofId& d) -> bool {
    if (d.is_point())
      return d.point->CORNER;
    if (d.is_line())
      return d.line->CORNER;
    return false;
  };

  // Helper: get BC values for a DOF
  auto get_dirichlet_phi = [&](const DofId& d) -> Complex {
    if (d.is_point() && bc.dirichlet_phi)
      return bc.dirichlet_phi(*d.point);
    if (d.is_line() && bc.dirichlet_phi_line)
      return bc.dirichlet_phi_line(*d.line);
    if (d.is_line() && bc.dirichlet_phi) {
      // Fallback: average of endpoint values
      auto [p0, p1] = d.line->getPoints();
      return 0.5 * (bc.dirichlet_phi(*p0) + bc.dirichlet_phi(*p1));
    }
    return Complex{0.0, 0.0};
  };
  auto get_neumann_phin = [&](const DofId& d) -> Complex {
    if (d.is_point() && bc.neumann_phin)
      return bc.neumann_phin(*d.point, d.face);
    if (d.is_line() && bc.neumann_phin_line)
      return bc.neumann_phin_line(*d.line, d.face);
    if (d.is_line() && bc.neumann_phin) {
      // Fallback: midpoint evaluation using adjacent face normals
      // (consistent with GMRES path in main_freq_domain.cpp)
      auto [p0, p1] = d.line->getPoints();
      // Use endpoint average for the callback value, but note:
      // for accurate midpoint evaluation, set neumann_phin_line explicitly
      return 0.5 * (bc.neumann_phin(*p0, d.face) + bc.neumann_phin(*p1, d.face));
    }
    return Complex{0.0, 0.0};
  };

  // Get Dirichlet index for CORNER continuity constraint
  auto get_dirichlet_index = [&](const DofId& d) -> int {
    if (d.is_point())
      return pf2Index(d.point, nullptr);
    if (d.is_line())
      return lf2Index(d.line, nullptr);
    return -1;
  };

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

  // Main assembly loop.
  for (std::size_t i = 0; i < n; ++i) {
    const auto& d_row = id_by_index[i];
    if (!d_row.is_point() && !d_row.is_line())
      throw std::runtime_error("solve_linear_bvp: id_by_index has empty DofId");

    const bool row_is_corner_neumann = dof_is_corner(d_row) && dof_is_neumann_key(d_row);
    if (row_is_corner_neumann) {
      // Replace BIE row by continuity constraint: phi(neumann-id) == phi(dirichlet-id).
      A[aidx(i, i)] = max_value;

      const int j_dir_i = get_dirichlet_index(d_row);
      if (j_dir_i < 0 || static_cast<std::size_t>(j_dir_i) >= n)
        throw std::runtime_error("solve_linear_bvp: missing Dirichlet ID at CORNER entity");
      const std::size_t j_dir = static_cast<std::size_t>(j_dir_i);

      const bool robin = is_entity_robin(d_row);
      if (robin) {
        const Complex kappa = kappa_at_dof(d_row);
        if (std::abs(kappa) < 1e-15)
          throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin CORNER");
        A[aidx(i, j_dir)] += -max_value / kappa;
        b[i] = Complex{0.0, 0.0};
      } else {
        b[i] = max_value * get_dirichlet_phi(d_row);
      }
      continue;
    }

    // Standard BIE row: sum IG*phin = sum IGn*phi (with substitutions for Robin).
    for (std::size_t j = 0; j < n; ++j) {
      const double IG = bvp.IGIGn[i][j][0];
      const double IGn = bvp.IGIGn[i][j][1];

      const auto& d_col = id_by_index[j];
      if (!d_col.is_point() && !d_col.is_line())
        throw std::runtime_error("solve_linear_bvp: id_by_index has empty DofId (col)");

      if (dof_is_dirichlet_key(d_col)) {
        // Unknown is phi_n (u_j).
        const bool robin = is_entity_robin(d_col);
        if (robin) {
          const Complex kappa = kappa_at_dof(d_col);
          if (std::abs(kappa) < 1e-15)
            throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin point");
          A[aidx(i, j)] += Complex{IG, 0.0} - Complex{IGn, 0.0} / kappa;
        } else {
          A[aidx(i, j)] += Complex{IG, 0.0};
          b[i] += Complex{IGn, 0.0} * get_dirichlet_phi(d_col);
        }
      } else {
        // Neumann ID: unknown is phi (u_j), phi_n is known.
        A[aidx(i, j)] += Complex{-IGn, 0.0};
        b[i] += -Complex{IG, 0.0} * get_neumann_phin(d_col);
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
    const auto& d = sol.id_by_index[j];
    if (dof_is_dirichlet_key(d)) {
      sol.phin[j] = sol.u[j];
      const bool robin = is_entity_robin(d);
      if (robin) {
        const Complex kappa = kappa_at_dof(d);
        sol.phi[j] = sol.u[j] / kappa;
      } else {
        sol.phi[j] = get_dirichlet_phi(d);
      }
    } else {
      sol.phi[j] = sol.u[j];
      sol.phin[j] = get_neumann_phin(d);
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
