#pragma once

#include <complex>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
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
// - Robin faces can be indexed either as "Dirichlet-type" (unknown is phi_n)
//   or as "Neumann-type" (unknown is phi). The default keeps the legacy phin
//   path; the phi policy is used by the frequency-domain open-boundary probe
//   because it restores BIE consistency for the current Sommerfeld tests.
// -----------------------------------------------------------------------------

namespace bem_frequency_domain {

using Complex = std::complex<double>;
using Id = std::tuple<BEM_DOF_Base *, networkFace *>;

struct CondensedInterfaceExpr {
  bool eliminated = false;
  std::size_t dir_index = 0;
  Complex scale{0.0, 0.0};
  Complex constant{0.0, 0.0};
};

enum class FaceBC {
  Neumann,   // phi_n is prescribed, unknown is phi (standard)
  Dirichlet, // phi is prescribed, unknown is phi_n (standard)
  Robin      // phi_n = kappa * phi (free-surface)
};

enum class RobinUnknownPolicy {
  Phin, // unknown is phi_n, reconstruct phi=(phi_n-rhs)/kappa
  Phi   // unknown is phi,   reconstruct phi_n=kappa*phi+rhs
};

inline const char *robin_unknown_policy_name(RobinUnknownPolicy policy) {
  switch (policy) {
  case RobinUnknownPolicy::Phin:
    return "phin";
  case RobinUnknownPolicy::Phi:
    return "phi";
  }
  return "unknown";
}

struct LinearFSBC {
  double omega = 0.0;   // [rad/s]
  double gravity = 9.81;
  double kappa_scale = 1.0; // diagnostic multiplier for the linear free-surface coefficient
  // Sponge coefficient mu(x) >= 0 (0 means no sponge). If null, treated as 0.
  std::function<double(const BEM_DOF_Base &)> sponge_mu;
  // Optional affine Robin offset r(x):
  //   phi_n = kappa(x) * phi + r(x)
  // This lets a relaxation zone nudge the solution toward a known target field
  // without damping that target itself. If null, r=0 and the legacy path is kept.
  std::function<Complex(const BEM_DOF_Base &)> robin_rhs;

  Complex kappa_at(const BEM_DOF_Base &node) const {
    const double mu = sponge_mu ? sponge_mu(node) : 0.0;
    const Complex s(mu, -omega); // (-i*omega + mu)
    if (gravity == 0.0)
      throw std::runtime_error("LinearFSBC: gravity must be non-zero");
    // Linear free-surface: (d/dt + mu)^2 phi + g * phi_n = 0 -> phi_n = -(s^2/g) * phi.
    return kappa_scale * (-(s * s) / gravity);
  }

  Complex robin_rhs_at(const BEM_DOF_Base &node) const {
    return robin_rhs ? robin_rhs(node) : Complex{0.0, 0.0};
  }
};

struct BoundaryData {
  // Face classification.
  std::function<FaceBC(const networkFace &)> face_bc;

  // Neumann value: prescribed phi_n on Neumann faces (may depend on face; face can be nullptr for single-ID nodes).
  std::function<Complex(const BEM_DOF_Base &, const networkFace *)> neumann_phin;

  // Dirichlet value: prescribed phi on Dirichlet faces (phi is pointwise unique in this codebase).
  std::function<Complex(const BEM_DOF_Base &)> dirichlet_phi;

  // Optional face-aware Robin coefficient overrides:
  //   phi_n = kappa(node, face) * phi + rhs(node, face).
  // If null, LinearFSBC::kappa_at()/robin_rhs_at() are used.  This is kept
  // optional so the standard free-surface path remains node-based.
  std::function<Complex(const BEM_DOF_Base &, const networkFace *)> robin_kappa;
  std::function<Complex(const BEM_DOF_Base &, const networkFace *)> robin_rhs_override;

  // Robin unknown policy.  The legacy bool is kept as an alias for existing
  // probes; new callers should set robin_unknown_policy explicitly.
  RobinUnknownPolicy robin_unknown_policy = RobinUnknownPolicy::Phin;

  // Legacy alias:
  //   false: unknown is phi_n, phi = (phi_n-rhs) / kappa (default)
  //   true : unknown is phi,   phi_n = kappa * phi + rhs
  bool robin_unknown_phi = false;

  bool robin_unknown_is_phi() const {
    return robin_unknown_phi || robin_unknown_policy == RobinUnknownPolicy::Phi;
  }

  // Diagnostic/controlled discretization hook. When true for a Neumann-capable
  // entity, Neumann DOFs are kept per face even if the local geometry is smooth.
  std::function<bool(const BEM_DOF_Base &)> force_neumann_multiple;

  // Diagnostic hook. When true for a BCInterface Neumann row, keep the normal
  // BIE row instead of replacing it with the phi-continuity constraint.
  std::function<bool(const BEM_DOF_Base &, const networkFace *)> use_bie_row_for_interface;

  // Diagnostic mode: eliminate BCInterface continuity constraints algebraically
  // before solving, instead of keeping the strong constraint rows in the LU
  // system. The full solution is reconstructed on the original IDs afterward.
  bool eliminate_interface_constraints = false;

  // Diagnostic mode for BCInterface continuity constraint row scaling.
  // Empty/default preserves the existing max-diagonal scale.  Other supported
  // values are "unit", "area", and "lumped_area".
  std::string interface_row_scale_mode;

  // Diagnostic probe: symmetrize only the algebraic pair between a
  // BCInterface constraint row and its active Dirichlet/Robin column.  This is
  // not a physical default; it is used to test whether raw reciprocity is
  // dominated by interface/Robin adjoint pairing.
  bool interface_robin_adjoint_probe = false;

  // Optional diagnostic hook called after IG/IGn assembly and ID indexing.
  // It must not mutate the mesh or the BVP.
  std::function<void(const BEM_BVP &, const std::vector<Id> &)> debug_bie_operator;

  // Optional experimental hook called after IG/IGn assembly and ID indexing,
  // before assembling the linear system. This is for diagnostic operator
  // transformations only; the default physical path leaves it empty.
  std::function<void(BEM_BVP &, const std::vector<Id> &)> transform_bie_operator;

  // Optional diagnostic hook called after assembling the complex linear
  // system A*u=b and before factorization. Matrix A is column-major.
  std::function<void(const std::vector<Complex> &,
                     const std::vector<Complex> &,
                     const std::vector<Id> &,
                     const std::unordered_set<networkFace *> &,
                     const std::unordered_map<BEM_DOF_Base *, bool> &)> debug_linear_system;

  // Optional diagnostic hook called after factorization/solve while the dense
  // matrix is still available. Arguments are A, original RHS, solution u, IDs,
  // Robin faces, and node-Robin map.
  std::function<void(const std::vector<Complex> &,
                     const std::vector<Complex> &,
                     const std::vector<Complex> &,
                     const std::vector<Id> &,
                     const std::unordered_set<networkFace *> &,
                     const std::unordered_map<BEM_DOF_Base *, bool> &)> debug_solved_linear_system;

  // Optional diagnostic hook for the matrix-free/FMM path.  When set,
  // solve_linear_bvp also assembles the dense operator as a reference, evaluates
  // both operators on the same solved vector, and passes the dense matrix/RHS,
  // FMM matvec/RHS, solution vector, and ID metadata to the caller.  This is
  // expensive and intended for small/linear diagnostic runs only.
  std::function<void(const std::vector<Complex> &,
                     const std::vector<Complex> &,
                     const std::function<std::vector<Complex>(const std::vector<Complex> &)> &,
                     const std::vector<Complex> &,
                     const std::vector<Complex> &,
                     const std::vector<Id> &,
                     const std::unordered_set<networkFace *> &,
                     const std::unordered_map<BEM_DOF_Base *, bool> &)> debug_fmm_matvec_compare;

  // Optional diagnostic hook for the algebraically condensed system used when
  // eliminate_interface_constraints=true.  Matrices are column-major.  The
  // retained_index vector maps original IDs to condensed indices; eliminated
  // entries store u_i = scale * u_dir + constant.
  std::function<void(const std::vector<Complex> &,
                     const std::vector<Complex> &,
                     const std::vector<Complex> &,
                     const std::vector<Id> &,
                     const std::vector<Id> &,
                     const std::vector<int> &,
                     const std::vector<CondensedInterfaceExpr> &,
                     const std::unordered_set<networkFace *> &,
                     const std::unordered_map<BEM_DOF_Base *, bool> &)> debug_condensed_solved_linear_system;

  // Diagnostic: compute A^T*u - A*u for the assembled algebraic system.
  // This is expensive for dense LU and is intended only for small debug runs.
  bool compute_adjoint_residual = false;

  // Dense complex linear solver for the assembled frequency-domain system.
  // "lu" is the existing LAPACK direct path. "gmres" uses the same dense
  // assembled matrix but avoids LU factorization; it is a bridge toward a
  // future matrix-free complex GMRES/FMM path.
  std::string linear_solver = "lu";
  double gmres_tol = 1e-9;
  int gmres_max_iter = 500;
  int gmres_restart = 100;
  std::string gmres_preconditioner = "diagonal";
  bool fmm_coordinate_scaling = true;
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

  // Row layout captured before boundary flags are restored after the solve.
  std::vector<unsigned char> interface_constraint_row;

  // Algebraic adjoint residual A^T*u - A*u.  For A*u=b, this measures the
  // non-symmetry seen by the solved vector under the bilinear transpose form.
  std::vector<Complex> adjoint_residual;
  double adjoint_residual_l2 = 0.0;

  std::string linear_solver;
  int linear_solver_iterations = 0;
  double linear_solver_residual_l2 = 0.0;
  double linear_solver_relative_residual_l2 = 0.0;
  bool linear_solver_converged = false;
  std::string linear_solver_effective_preconditioner;

  std::string fmm_setup_source;
  bool fmm_coordinate_scaling = false;
  bool fmm_morton_reindex = false;
  bool fmm_reused_sources = false;
  bool fmm_reused_static = false;
  std::size_t fmm_targets = 0;
  std::size_t fmm_sources = 0;
  std::size_t fmm_total_near_terms = 0;
  double fmm_mean_near_terms_per_target = 0.0;
  std::size_t fmm_max_near_terms_per_target = 0;
  std::size_t fmm_vertex_targets = 0;
  std::size_t fmm_vertex_total_near_terms = 0;
  double fmm_vertex_mean_near_terms_per_target = 0.0;
  std::size_t fmm_vertex_max_near_terms_per_target = 0;
  std::size_t fmm_midpoint_targets = 0;
  std::size_t fmm_midpoint_total_near_terms = 0;
  double fmm_midpoint_mean_near_terms_per_target = 0.0;
  std::size_t fmm_midpoint_max_near_terms_per_target = 0;
  double fmm_max_source_offset = 0.0;
  double fmm_mean_source_offset = 0.0;
  double fmm_p95_source_offset = 0.0;
  double fmm_coordinate_scale_factor = 1.0;
  std::vector<std::pair<int, double>> fmm_top_source_offsets;
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

inline Complex complex_dot_conj(const std::vector<Complex> &a, const std::vector<Complex> &b) {
  Complex out{0.0, 0.0};
  const std::size_t n = std::min(a.size(), b.size());
  for (std::size_t i = 0; i < n; ++i)
    out += std::conj(a[i]) * b[i];
  return out;
}

inline std::vector<Complex> dense_matvec_col_major(const std::vector<Complex> &A, const std::vector<Complex> &x, std::size_t n) {
  std::vector<Complex> y(n, Complex{0.0, 0.0});
  for (std::size_t col = 0; col < n; ++col) {
    const Complex xc = x[col];
    if (std::abs(xc) == 0.0)
      continue;
    const Complex *a_col = A.data() + col * n;
    for (std::size_t row = 0; row < n; ++row)
      y[row] += a_col[row] * xc;
  }
  return y;
}

inline std::vector<Complex> solve_small_dense(std::vector<Complex> A, std::vector<Complex> b, std::size_t n) {
  auto idx = [n](std::size_t row, std::size_t col) { return row + col * n; };
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t pivot = k;
    double pivot_abs = std::abs(A[idx(k, k)]);
    for (std::size_t r = k + 1; r < n; ++r) {
      const double v = std::abs(A[idx(r, k)]);
      if (v > pivot_abs) {
        pivot = r;
        pivot_abs = v;
      }
    }
    if (!(pivot_abs > 0.0))
      throw std::runtime_error("solve_small_dense: singular normal equation");
    if (pivot != k) {
      for (std::size_t c = k; c < n; ++c)
        std::swap(A[idx(k, c)], A[idx(pivot, c)]);
      std::swap(b[k], b[pivot]);
    }
    const Complex diag = A[idx(k, k)];
    for (std::size_t r = k + 1; r < n; ++r) {
      const Complex factor = A[idx(r, k)] / diag;
      if (std::abs(factor) == 0.0)
        continue;
      A[idx(r, k)] = Complex{0.0, 0.0};
      for (std::size_t c = k + 1; c < n; ++c)
        A[idx(r, c)] -= factor * A[idx(k, c)];
      b[r] -= factor * b[k];
    }
  }
  std::vector<Complex> x(n, Complex{0.0, 0.0});
  for (std::size_t rr = 0; rr < n; ++rr) {
    const std::size_t i = n - 1 - rr;
    Complex sum = b[i];
    for (std::size_t c = i + 1; c < n; ++c)
      sum -= A[idx(i, c)] * x[c];
    x[i] = sum / A[idx(i, i)];
  }
  return x;
}

inline std::vector<Complex> least_squares_normal_eq(const std::vector<std::vector<Complex>> &h,
                                                    const std::vector<Complex> &g,
                                                    std::size_t rows,
                                                    std::size_t cols) {
  std::vector<Complex> normal(cols * cols, Complex{0.0, 0.0});
  std::vector<Complex> rhs(cols, Complex{0.0, 0.0});
  auto nidx = [cols](std::size_t row, std::size_t col) { return row + col * cols; };
  for (std::size_t i = 0; i < cols; ++i) {
    for (std::size_t r = 0; r < rows; ++r)
      rhs[i] += std::conj(h[r][i]) * g[r];
    for (std::size_t j = 0; j < cols; ++j) {
      Complex v{0.0, 0.0};
      for (std::size_t r = 0; r < rows; ++r)
        v += std::conj(h[r][i]) * h[r][j];
      normal[nidx(i, j)] = v;
    }
  }
  return solve_small_dense(std::move(normal), std::move(rhs), cols);
}

struct DenseGMRESResult {
  std::vector<Complex> x;
  int iterations = 0;
  double residual_l2 = 0.0;
  bool converged = false;
};

inline DenseGMRESResult solve_complex_gmres_operator(const std::function<std::vector<Complex>(const std::vector<Complex> &)> &matvec,
                                                     const std::vector<Complex> &b,
                                                     int restart_in,
                                                     int max_iter_in,
                                                     double tol_in,
                                                     const std::vector<Complex> &left_scale_in = {}) {
  const std::size_t n = b.size();
  if (n == 0)
    return {};
  if (!matvec)
    throw std::runtime_error("solve_complex_gmres_operator: empty matvec");
  const int restart = std::max(1, std::min<int>(restart_in > 0 ? restart_in : 50, static_cast<int>(n)));
  const int max_iter = std::max(restart, max_iter_in > 0 ? max_iter_in : restart);
  const double tol = (std::isfinite(tol_in) && tol_in > 0.0) ? tol_in : 1e-9;

  std::vector<Complex> left_scale = left_scale_in;
  if (left_scale.empty())
    left_scale.assign(n, Complex{1.0, 0.0});
  if (left_scale.size() != n)
    throw std::runtime_error("solve_complex_gmres_operator: left_scale size mismatch");

  DenseGMRESResult result;
  result.x.assign(n, Complex{0.0, 0.0});

  std::vector<Complex> b_pre = b;
  for (std::size_t i = 0; i < n; ++i)
    b_pre[i] *= left_scale[i];
  const double norm_b_pre = std::max(l2_norm(b_pre), 1e-300);

  auto matvec_pre = [&](const std::vector<Complex> &x) {
    auto y = matvec(x);
    if (y.size() != n)
      throw std::runtime_error("solve_complex_gmres_operator: matvec size mismatch");
    for (std::size_t i = 0; i < n; ++i)
      y[i] *= left_scale[i];
    return y;
  };
  auto residual_pre = [&]() {
    auto ax = matvec_pre(result.x);
    for (std::size_t i = 0; i < n; ++i)
      ax[i] = b_pre[i] - ax[i];
    return ax;
  };
  auto residual_true_l2 = [&](const std::vector<Complex> &x) {
    auto ax = matvec(x);
    if (ax.size() != n)
      throw std::runtime_error("solve_complex_gmres_operator: matvec size mismatch");
    for (std::size_t i = 0; i < n; ++i)
      ax[i] = b[i] - ax[i];
    return l2_norm(ax);
  };

  std::vector<Complex> r = residual_pre();
  double residual_pre_l2 = l2_norm(r);
  result.residual_l2 = residual_true_l2(result.x);
  if (residual_pre_l2 <= tol * norm_b_pre) {
    result.converged = true;
    return result;
  }

  while (result.iterations < max_iter) {
    const double beta = l2_norm(r);
    if (!(beta > 0.0) || !std::isfinite(beta))
      break;
    std::vector<std::vector<Complex>> v(static_cast<std::size_t>(restart) + 1, std::vector<Complex>(n, Complex{0.0, 0.0}));
    for (std::size_t i = 0; i < n; ++i)
      v[0][i] = r[i] / beta;
    std::vector<std::vector<Complex>> h(static_cast<std::size_t>(restart) + 1, std::vector<Complex>(static_cast<std::size_t>(restart), Complex{0.0, 0.0}));
    std::vector<Complex> g(static_cast<std::size_t>(restart) + 1, Complex{0.0, 0.0});
    g[0] = Complex{beta, 0.0};

    int inner_done = 0;
    std::vector<Complex> best_x = result.x;
    double best_res = result.residual_l2;

    for (int j = 0; j < restart && result.iterations < max_iter; ++j) {
      std::vector<Complex> w = matvec_pre(v[static_cast<std::size_t>(j)]);
      for (int i = 0; i <= j; ++i) {
        h[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = complex_dot_conj(v[static_cast<std::size_t>(i)], w);
        const Complex hij = h[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        for (std::size_t k = 0; k < n; ++k)
          w[k] -= hij * v[static_cast<std::size_t>(i)][k];
      }
      const double h_next = l2_norm(w);
      h[static_cast<std::size_t>(j + 1)][static_cast<std::size_t>(j)] = Complex{h_next, 0.0};
      if (h_next > 0.0) {
        for (std::size_t k = 0; k < n; ++k)
          v[static_cast<std::size_t>(j + 1)][k] = w[k] / h_next;
      }
      ++result.iterations;
      inner_done = j + 1;

      const auto y = least_squares_normal_eq(h, g, static_cast<std::size_t>(inner_done) + 1, static_cast<std::size_t>(inner_done));
      std::vector<Complex> candidate = result.x;
      for (int col = 0; col < inner_done; ++col) {
        for (std::size_t row = 0; row < n; ++row)
          candidate[row] += v[static_cast<std::size_t>(col)][row] * y[static_cast<std::size_t>(col)];
      }
      auto ax = matvec(candidate);
      if (ax.size() != n)
        throw std::runtime_error("solve_complex_gmres_operator: matvec size mismatch");
      for (std::size_t row = 0; row < n; ++row)
        ax[row] = b[row] - ax[row];
      const double candidate_res = l2_norm(ax);
      auto ax_pre = matvec_pre(candidate);
      for (std::size_t row = 0; row < n; ++row)
        ax_pre[row] = b_pre[row] - ax_pre[row];
      const double candidate_pre_res = l2_norm(ax_pre);
      if (candidate_res < best_res) {
        best_res = candidate_res;
        best_x = std::move(candidate);
      }
      if (candidate_pre_res <= tol * norm_b_pre) {
        result.x = std::move(best_x);
        result.residual_l2 = candidate_res;
        result.converged = true;
        return result;
      }
      if (h_next == 0.0)
        break;
    }

    result.x = std::move(best_x);
    result.residual_l2 = best_res;
    r = residual_pre();
    residual_pre_l2 = l2_norm(r);
    if (residual_pre_l2 <= tol * norm_b_pre) {
      result.converged = true;
      return result;
    }
    if (inner_done == 0)
      break;
  }

  return result;
}

inline DenseGMRESResult solve_dense_complex_gmres(const std::vector<Complex> &A,
                                                  const std::vector<Complex> &b,
                                                  int restart_in,
                                                  int max_iter_in,
                                                  double tol_in,
                                                  const std::string &preconditioner_in = "diagonal") {
  const std::size_t n = b.size();
  if (n == 0)
    return {};
  const int restart = std::max(1, std::min<int>(restart_in > 0 ? restart_in : 50, static_cast<int>(n)));
  const int max_iter = std::max(restart, max_iter_in > 0 ? max_iter_in : restart);
  const double tol = (std::isfinite(tol_in) && tol_in > 0.0) ? tol_in : 1e-9;

  DenseGMRESResult result;
  result.x.assign(n, Complex{0.0, 0.0});
  auto aidx = [n](std::size_t row, std::size_t col) { return row + col * n; };
  const std::string preconditioner = preconditioner_in.empty() ? "diagonal" : preconditioner_in;
  std::vector<Complex> left_scale(n, Complex{1.0, 0.0});
  if (preconditioner == "diagonal") {
    for (std::size_t i = 0; i < n; ++i) {
      const Complex d = A[aidx(i, i)];
      if (std::abs(d) > 1e-300)
        left_scale[i] = Complex{1.0, 0.0} / d;
    }
  } else if (preconditioner == "row_norm") {
    std::vector<double> row_norm2(n, 0.0);
    for (std::size_t col = 0; col < n; ++col) {
      for (std::size_t row = 0; row < n; ++row)
        row_norm2[row] += std::norm(A[aidx(row, col)]);
    }
    for (std::size_t i = 0; i < n; ++i) {
      const double row_norm = std::sqrt(row_norm2[i]);
      if (row_norm > 1e-300)
        left_scale[i] = Complex{1.0 / row_norm, 0.0};
    }
  } else if (preconditioner == "none") {
    // Preserve unit left scaling.
  } else {
    throw std::runtime_error("solve_dense_complex_gmres: unsupported preconditioner=" + preconditioner);
  }

  return solve_complex_gmres_operator(
      [&](const std::vector<Complex> &x) { return dense_matvec_col_major(A, x, n); },
      b,
      restart,
      max_iter,
      tol,
      left_scale);
}

struct BoundaryStateGuard {
  struct FaceState {
    networkFace *f = nullptr;
    bool isLinearElement = true;
    bool isPseudoQuadraticElement = false;
    bool isTrueQuadraticElement = false;
    Network *penetratedBody = nullptr;
  };
  struct LineState {
    networkLine *l = nullptr;
    bool Dirichlet = false;
    bool Neumann = false;
    bool BCInterface = false;
    bool isMultipleNode = false;
    std::unordered_map<networkFace *, NodeFaceState> dofs;
    Tdd phiphin{};
    Tdd phiphin_t{};
  };
  struct PointState {
    networkPoint *p = nullptr;
    bool Dirichlet = false;
    bool Neumann = false;
    bool BCInterface = false;
    bool isMultipleNode = false;
    std::unordered_map<networkFace *, NodeFaceState> dofs;
    Tdd phiphin{};
    Tdd phiphin_t{};
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
      faces.push_back(FaceState{f, f->isLinearElement, f->isPseudoQuadraticElement, f->isTrueQuadraticElement, f->penetratedBody});
    }
    lines.reserve(uniq_lines.size());
    for (auto *l : uniq_lines) {
      lines.push_back(LineState{l, l->Dirichlet, l->Neumann, l->BCInterface, l->isMultipleNode, l->dofs, l->phiphin, l->phiphin_t});
    }
    points.reserve(uniq_points.size());
    for (auto *p : uniq_points) {
      points.push_back(PointState{p, p->Dirichlet, p->Neumann, p->BCInterface, p->isMultipleNode, p->dofs, p->phiphin, p->phiphin_t});
    }
  }

  void restore() const {
    for (const auto &s : faces) {
      if (!s.f)
        continue;
      s.f->isLinearElement = s.isLinearElement;
      s.f->isPseudoQuadraticElement = s.isPseudoQuadraticElement;
      s.f->isTrueQuadraticElement = s.isTrueQuadraticElement;
      s.f->penetratedBody = s.penetratedBody;
    }
    for (const auto &s : lines) {
      if (!s.l)
        continue;
      s.l->Dirichlet = s.Dirichlet;
      s.l->Neumann = s.Neumann;
      s.l->BCInterface = s.BCInterface;
      s.l->isMultipleNode = s.isMultipleNode;
      s.l->dofs = s.dofs;
      s.l->phiphin = s.phiphin;
      s.l->phiphin_t = s.phiphin_t;
    }
    for (const auto &s : points) {
      if (!s.p)
        continue;
      s.p->Dirichlet = s.Dirichlet;
      s.p->Neumann = s.Neumann;
      s.p->BCInterface = s.BCInterface;
      s.p->isMultipleNode = s.isMultipleNode;
      s.p->dofs = s.dofs;
      s.p->phiphin = s.phiphin;
      s.p->phiphin_t = s.phiphin_t;
    }
  }
};

inline void apply_face_bc_to_mesh(const std::vector<Network *> &nets, const BoundaryData &bc, std::unordered_set<networkFace *> &robin_faces_out) {
  robin_faces_out.clear();
  if (!bc.face_bc)
    throw std::runtime_error("apply_face_bc_to_mesh: face_bc is not set");

  std::unordered_map<networkFace *, FaceBC> face_bc_map;

  for (auto *net : nets) {
    if (!net)
      continue;
    ensureLineXmidInitialized(net);
    net->setGeometricPropertiesForce();
    for (auto *f : net->getBoundaryFaces()) {
      const FaceBC t = bc.face_bc(*f);
      face_bc_map[f] = t;
      if (t == FaceBC::Robin)
        robin_faces_out.emplace(f);
    }

    auto apply_node_face_state = [&](auto *node) {
      for (auto *f : node->getBoundaryFaces()) {
        auto it = face_bc_map.find(f);
        auto &d = node->dof(f);
        d.contact_opponent_faces.clear();
        if (it != face_bc_map.end() && (it->second == FaceBC::Neumann || (it->second == FaceBC::Robin && bc.robin_unknown_is_phi())))
          d.contact_opponent_faces.push_back(f);
        d.detached_by_pressure = false;
      }
    };

    for (auto *l : net->getBoundaryLines()) {
      apply_node_face_state(l);
      const auto summary = getEdgeSummarizedBoundaryType(l);
      l->Neumann = (summary == SummarizedNodeBoundaryType::Neumann);
      l->Dirichlet = (summary == SummarizedNodeBoundaryType::Dirichlet);
      l->BCInterface = (summary == SummarizedNodeBoundaryType::Interface);
      setMultipleNode(l);
      if (bc.force_neumann_multiple && hasAnyNeumannBoundaryState(l) && bc.force_neumann_multiple(*l))
        l->isMultipleNode = true;
    }

    for (auto *p : net->getBoundaryPoints()) {
      apply_node_face_state(p);
      const auto summary = getVertexSummarizedBoundaryType(p);
      p->Neumann = (summary == SummarizedNodeBoundaryType::Neumann);
      p->Dirichlet = (summary == SummarizedNodeBoundaryType::Dirichlet);
      p->BCInterface = (summary == SummarizedNodeBoundaryType::Interface);
      setMultipleNode(p);
      if (bc.force_neumann_multiple && hasAnyNeumannBoundaryState(p) && bc.force_neumann_multiple(*p))
        p->isMultipleNode = true;
    }

    for (auto *f : net->getBoundaryFaces()) {
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

inline Solution solve_linear_bvp(const std::vector<Network *> &nets, const LinearFSBC &fsbc, const BoundaryData &bc) {
  if (nets.empty())
    throw std::runtime_error("solve_linear_bvp: empty networks");

  BoundaryStateGuard guard(nets);

  std::unordered_set<networkFace *> robin_faces;
  apply_face_bc_to_mesh(nets, bc, robin_faces);

  // Build ID indices (this fills BEM_DOF_Base::dofs[*].index).
  const std::size_t n = setNodeFaceIndices(nets, [&](const BEM_DOF_Base *node) {
    return node && bc.force_neumann_multiple && bc.force_neumann_multiple(*node);
  });
  if (n == 0)
    throw std::runtime_error("solve_linear_bvp: matrix size is 0");

  // Frequency-domain Robin coefficients depend on physical coordinates through
  // the sponge profile.  The shared time-domain FMM setup may scale coordinates,
  // so cache kappa before entering the FMM session.
  std::unordered_map<const BEM_DOF_Base *, Complex> physical_kappa_cache;
  std::unordered_map<const BEM_DOF_Base *, Complex> physical_robin_rhs_cache;
  auto cache_physical_robin_values = [&](const BEM_DOF_Base *node) {
    if (!node || physical_kappa_cache.contains(node))
      return;
    bool touches_robin = false;
    for (auto *f : node->getBoundaryFaces()) {
      if (robin_faces.contains(f)) {
        touches_robin = true;
        break;
      }
    }
    if (touches_robin) {
      physical_kappa_cache.emplace(node, fsbc.kappa_at(*node));
      physical_robin_rhs_cache.emplace(node, fsbc.robin_rhs_at(*node));
    }
  };
  for (auto *net : nets) {
    if (!net)
      continue;
    for (auto *p : net->getBoundaryPoints())
      cache_physical_robin_values(p);
    for (auto *l : net->getBoundaryLines())
      cache_physical_robin_values(l);
  }

  // Prepare BIE coefficients (dense, geometry-only).
  BEM_BVP bvp(nets);
  bvp.matrix_size = static_cast<int>(n);
  std::unique_ptr<BEM_BVP::FMMMatvecSession> fmm_session;
  const bool use_matrix_free_fmm = (bc.linear_solver == "fmm_gmres");
  const bool build_dense_fmm_reference = use_matrix_free_fmm && static_cast<bool>(bc.debug_fmm_matvec_compare);

  auto collect_id_by_index = [&]() {
    std::vector<Id> ids(n);
    auto collect_active_ids = [&](auto *node) {
      for (const auto &[f, d] : node->dofs) {
        if (d.index < 0 || static_cast<std::size_t>(d.index) >= n)
          continue;
        ids[static_cast<std::size_t>(d.index)] = {node, f};
      }
    };
    for (auto *net : nets) {
      for (auto *p : net->getBoundaryPoints())
        collect_active_ids(p);
      for (auto *l : net->getBoundaryLines())
        collect_active_ids(l);
    }
    return ids;
  };

  std::vector<Id> pre_fmm_id_by_index;
  struct GlobalPreconditionerGuard {
    std::string saved;
    bool active = false;
    ~GlobalPreconditionerGuard() {
      if (active)
        preconditioner_type = saved;
    }
  } preconditioner_guard{preconditioner_type, use_matrix_free_fmm};
  if (use_matrix_free_fmm)
    preconditioner_type = "NONE";
  if (use_matrix_free_fmm) {
    if (bc.eliminate_interface_constraints)
      throw std::runtime_error("solve_linear_bvp: fmm_gmres does not support eliminate_interface_constraints yet");
    if (bc.interface_robin_adjoint_probe)
      throw std::runtime_error("solve_linear_bvp: fmm_gmres does not support interface_robin_adjoint_probe yet");
    if (bc.transform_bie_operator || bc.debug_bie_operator || bc.debug_linear_system || bc.debug_solved_linear_system || bc.debug_condensed_solved_linear_system || bc.compute_adjoint_residual)
      throw std::runtime_error("solve_linear_bvp: fmm_gmres currently supports solve/residual only, not dense operator diagnostics");
    if (build_dense_fmm_reference) {
      pre_fmm_id_by_index = collect_id_by_index();
      bvp.setIGIGn();
    }
    fmm_session = bvp.prepareFMMMatvecSession(bc.fmm_coordinate_scaling, true, true, "time_domain_shared");
  } else {
    bvp.setIGIGn();
  }

  std::vector<Id> id_by_index = collect_id_by_index();
  if (bc.transform_bie_operator)
    bc.transform_bie_operator(bvp, id_by_index);
  if (bc.debug_bie_operator)
    bc.debug_bie_operator(bvp, id_by_index);

  // Classify Dirichlet IDs into Robin-vs-true-Dirichlet using point adjacency.
  std::unordered_map<BEM_DOF_Base *, bool> node_is_robin;
  node_is_robin.reserve(n);
  for (const auto &[p, f] : id_by_index) {
    (void)f;
    if (!p)
      continue;
    if (node_is_robin.contains(p))
      continue;
    bool is_robin = false;
    for (auto *adj : p->getBoundaryFaces()) {
      if (robin_faces.contains(adj)) {
        is_robin = true;
        break;
      }
    }
    node_is_robin.emplace(p, is_robin);
  }

  // Scaling for constraint rows (mimic generateBIEMatrix() idea).
  double max_value = 1.0;
  if (use_matrix_free_fmm) {
    for (const double d0 : bvp.diag_coeffs) {
      const double d = std::abs(d0);
      if (d > max_value)
        max_value = d;
    }
  } else {
    for (std::size_t i = 0; i < n; ++i) {
      const double d = std::abs(bvp.IGIGn[i][i][0]);
      if (d > max_value)
        max_value = d;
    }
  }

  auto interface_row_scale = [&](const BEM_DOF_Base *p, const networkFace *f) -> double {
    if (bc.interface_row_scale_mode.empty() || bc.interface_row_scale_mode == "default")
      return max_value;
    if (bc.interface_row_scale_mode == "unit")
      return 1.0;
    if (bc.interface_row_scale_mode == "area")
      return f ? std::max(f->area, 1e-300) : max_value;
    if (bc.interface_row_scale_mode == "lumped_area") {
      if (!p)
        return max_value;
      long double a = 0.0L;
      for (auto *adj : p->getBoundaryFaces()) {
        if (adj)
          a += static_cast<long double>(adj->area) / 3.0L;
      }
      return std::max(static_cast<double>(a), 1e-300);
    }
    throw std::runtime_error("solve_linear_bvp: unsupported interface_row_scale_mode=" + bc.interface_row_scale_mode);
  };

  // Assemble A*u=b (complex), column-major unless the FMM path supplies A as
  // a matrix-free operator.
  std::vector<Complex> A;
  if (!use_matrix_free_fmm)
    A.assign(n * n, Complex{0.0, 0.0});
  std::vector<Complex> b(n, Complex{0.0, 0.0});

  auto aidx = [&](std::size_t row, std::size_t col) -> std::size_t { return row + col * n; };

  auto get_dirichlet_phi = [&](const BEM_DOF_Base *p) -> Complex {
    if (!bc.dirichlet_phi)
      return Complex{0.0, 0.0};
    return bc.dirichlet_phi(*p);
  };
  auto get_neumann_phin = [&](const BEM_DOF_Base *p, const networkFace *f) -> Complex {
    if (!bc.neumann_phin)
      return Complex{0.0, 0.0};
    return bc.neumann_phin(*p, f);
  };
  auto is_robin_id = [&](const BEM_DOF_Base *p, const networkFace *f) -> bool {
    if (f)
      return robin_faces.contains(const_cast<networkFace *>(f));
    if (!p)
      return false;
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base *>(p));
    return it != node_is_robin.end() && it->second;
  };

  auto get_kappa = [&](const BEM_DOF_Base *p, const networkFace *f) -> Complex {
    if (!p)
      throw std::runtime_error("solve_linear_bvp: kappa requested for null node");
    if (bc.robin_kappa)
      return bc.robin_kappa(*p, f);
    const auto it = physical_kappa_cache.find(p);
    if (it != physical_kappa_cache.end())
      return it->second;
    return fsbc.kappa_at(*p);
  };
  auto get_robin_rhs = [&](const BEM_DOF_Base *p, const networkFace *f) -> Complex {
    if (!p)
      throw std::runtime_error("solve_linear_bvp: Robin RHS requested for null node");
    if (bc.robin_rhs_override)
      return bc.robin_rhs_override(*p, f);
    const auto it = physical_robin_rhs_cache.find(p);
    if (it != physical_robin_rhs_cache.end())
      return it->second;
    return fsbc.robin_rhs_at(*p);
  };

  auto is_interface_constraint_row = [&](const BEM_DOF_Base *p_row, const networkFace *f_row) -> bool {
    return p_row && p_row->BCInterface && isNeumannBieDofKey(p_row, const_cast<networkFace *>(f_row)) &&
           !(bc.use_bie_row_for_interface && bc.use_bie_row_for_interface(*p_row, const_cast<networkFace *>(f_row)));
  };

  std::vector<std::size_t> coeff_index_by_current(n);
  for (std::size_t i = 0; i < n; ++i)
    coeff_index_by_current[i] = i;
  if (build_dense_fmm_reference) {
    struct LocalIdHash {
      std::size_t operator()(const Id &id) const noexcept {
        const auto *p = std::get<0>(id);
        const auto *f = std::get<1>(id);
        return std::hash<const void *>{}(p) ^ (std::hash<const void *>{}(f) << 1);
      }
    };
    struct LocalIdEq {
      bool operator()(const Id &a, const Id &b) const noexcept {
        return std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b);
      }
    };
    std::unordered_map<Id, std::size_t, LocalIdHash, LocalIdEq> old_index_by_id;
    old_index_by_id.reserve(pre_fmm_id_by_index.size());
    for (std::size_t i = 0; i < pre_fmm_id_by_index.size(); ++i)
      old_index_by_id.emplace(pre_fmm_id_by_index[i], i);
    for (std::size_t i = 0; i < n; ++i) {
      const auto it = old_index_by_id.find(id_by_index[i]);
      if (it == old_index_by_id.end())
        throw std::runtime_error("solve_linear_bvp: failed to map FMM-reindexed ID back to dense reference ID");
      coeff_index_by_current[i] = it->second;
    }
  }

  auto assemble_dense_system = [&](const std::vector<std::size_t> &coeff_index) {
    std::vector<Complex> A_out(n * n, Complex{0.0, 0.0});
    std::vector<Complex> b_out(n, Complex{0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i) {
      auto [p_row, f_row] = id_by_index[i];
      if (!p_row)
        throw std::runtime_error("solve_linear_bvp: id_by_index has null point");

      const bool row_is_corner_neumann = is_interface_constraint_row(p_row, f_row);
      if (row_is_corner_neumann) {
        // Replace BIE row by continuity constraint: phi(neumann-id) == phi(dirichlet-id).
        const double constraint_scale = interface_row_scale(p_row, f_row);
        A_out[aidx(i, i)] = constraint_scale;

        const auto *d_dir = p_row->findActiveBieDof(nullptr);
        if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
          throw std::runtime_error("solve_linear_bvp: missing Dirichlet ID at BCInterface point");
        const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);

        const bool robin = node_is_robin[p_row];
        if (robin) {
          const Complex kappa = get_kappa(p_row, f_row);
          if (std::abs(kappa) < 1e-15)
            throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin BCInterface point");
          const Complex r = get_robin_rhs(p_row, f_row);
          A_out[aidx(i, j_dir)] += -constraint_scale / kappa; // phi_dir = (phin-r)/kappa
          b_out[i] = -constraint_scale * r / kappa;
        } else {
          b_out[i] = constraint_scale * get_dirichlet_phi(p_row);
        }
        continue;
      }

      // Standard BIE row: sum IG*phin = sum IGn*phi (with substitutions for Robin).
      const std::size_t ci = coeff_index[i];
      for (std::size_t j = 0; j < n; ++j) {
        const std::size_t cj = coeff_index[j];
        const double IG = bvp.IGIGn[ci][cj][0];
        const double IGn = bvp.IGIGn[ci][cj][1];

        auto [p_col, f_col] = id_by_index[j];
        if (!p_col)
          throw std::runtime_error("solve_linear_bvp: id_by_index has null point (col)");

        if (isDirichletBieDofKey(p_col, f_col)) {
          // Unknown is phi_n (u_j).
          const bool robin = is_robin_id(p_col, f_col);
          if (robin) {
            if (bc.robin_unknown_is_phi())
              throw std::runtime_error("solve_linear_bvp: Robin Dirichlet-type ID appears while robin_unknown_phi=true");
            const Complex kappa = get_kappa(p_col, f_col);
            if (std::abs(kappa) < 1e-15)
              throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin point");
            const Complex r = get_robin_rhs(p_col, f_col);
            A_out[aidx(i, j)] += Complex{IG, 0.0} - Complex{IGn, 0.0} / kappa; // phi = (phin-r)/kappa
            b_out[i] += -Complex{IGn, 0.0} * r / kappa;
          } else {
            A_out[aidx(i, j)] += Complex{IG, 0.0};
            b_out[i] += Complex{IGn, 0.0} * get_dirichlet_phi(p_col);
          }
        } else {
          // Neumann ID: unknown is phi (u_j), phi_n is known.
          const bool robin = bc.robin_unknown_is_phi() && is_robin_id(p_col, f_col);
          if (robin) {
            const Complex kappa = get_kappa(p_col, f_col);
            const Complex r = get_robin_rhs(p_col, f_col);
            A_out[aidx(i, j)] += Complex{IG, 0.0} * kappa - Complex{IGn, 0.0}; // phi_n = kappa * phi
            b_out[i] += -Complex{IG, 0.0} * r;
          } else {
            A_out[aidx(i, j)] += Complex{-IGn, 0.0};
            b_out[i] += -Complex{IG, 0.0} * get_neumann_phin(p_col, f_col);
          }
        }
      }
    }
    return std::pair<std::vector<Complex>, std::vector<Complex>>{std::move(A_out), std::move(b_out)};
  };

  std::vector<Complex> dense_fmm_reference_A;
  std::vector<Complex> dense_fmm_reference_rhs;
  if (!use_matrix_free_fmm || build_dense_fmm_reference) {
    auto assembled = assemble_dense_system(coeff_index_by_current);
    if (use_matrix_free_fmm) {
      dense_fmm_reference_A = std::move(assembled.first);
      dense_fmm_reference_rhs = std::move(assembled.second);
    } else {
      A = std::move(assembled.first);
      b = std::move(assembled.second);
    }
  }

  if (bc.debug_linear_system)
    bc.debug_linear_system(A, b, id_by_index, robin_faces, node_is_robin);

  if (bc.interface_robin_adjoint_probe) {
    for (std::size_t i = 0; i < n; ++i) {
      auto [p_row, f_row] = id_by_index[i];
      const bool row_is_corner_neumann =
          (p_row && p_row->BCInterface && isNeumannBieDofKey(p_row, f_row) &&
           !(bc.use_bie_row_for_interface && bc.use_bie_row_for_interface(*p_row, f_row)));
      if (!row_is_corner_neumann)
        continue;
      const auto *d_dir = p_row->findActiveBieDof(nullptr);
      if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
        continue;
      const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);
      if (!is_robin_id(p_row, nullptr))
        continue;
      const Complex sym = 0.5 * (A[aidx(i, j_dir)] + A[aidx(j_dir, i)]);
      A[aidx(i, j_dir)] = sym;
      A[aidx(j_dir, i)] = sym;
    }
  }

  const double dense_rhs_l2 = use_matrix_free_fmm ? 0.0 : l2_norm(b);
  const bool keep_original_rhs = bc.compute_adjoint_residual || static_cast<bool>(bc.debug_solved_linear_system);
  const std::vector<Complex> rhs_original = keep_original_rhs ? b : std::vector<Complex>{};

  std::vector<Complex> u_full;
  int solver_iterations = 0;
  double solver_residual_l2 = 0.0;
  double solver_relative_residual_l2 = 0.0;
  bool solver_converged = true;
  std::string solver_effective_preconditioner = (bc.linear_solver == "lu" || bc.linear_solver.empty()) ? "none" : bc.gmres_preconditioner;
  std::function<std::vector<Complex>(const std::vector<Complex> &)> matrix_free_matvec;
  std::vector<Complex> matrix_free_rhs;

  if (use_matrix_free_fmm) {
    auto reconstruct_fields = [&](const std::vector<Complex> &u, const bool unknown_part) {
      std::vector<Complex> phi(n, Complex{0.0, 0.0});
      std::vector<Complex> phin(n, Complex{0.0, 0.0});
      for (std::size_t j = 0; j < n; ++j) {
        auto [p_col, f_col] = id_by_index[j];
        if (!p_col)
          throw std::runtime_error("solve_linear_bvp: id_by_index has null point (fmm)");
        if (isDirichletBieDofKey(p_col, f_col)) {
          const bool robin = is_robin_id(p_col, f_col);
          if (unknown_part) {
            phin[j] = u[j];
            if (robin) {
              if (bc.robin_unknown_is_phi())
                throw std::runtime_error("solve_linear_bvp: fmm_gmres Robin Dirichlet-type ID with robin_unknown_phi=true");
              const Complex kappa = get_kappa(p_col, f_col);
              if (std::abs(kappa) < 1e-15)
                throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin point (fmm)");
              phi[j] = u[j] / kappa;
            }
          } else if (robin) {
            const Complex kappa = get_kappa(p_col, f_col);
            if (std::abs(kappa) < 1e-15)
              throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin point (fmm known)");
            phi[j] = -get_robin_rhs(p_col, f_col) / kappa;
          } else if (!robin) {
            phi[j] = get_dirichlet_phi(p_col);
          }
        } else {
          const bool robin = bc.robin_unknown_is_phi() && is_robin_id(p_col, f_col);
          if (unknown_part) {
            phi[j] = u[j];
            if (robin) {
              const Complex kappa = get_kappa(p_col, f_col);
              phin[j] = kappa * u[j];
            }
          } else if (robin) {
            phin[j] = get_robin_rhs(p_col, f_col);
          } else if (!robin) {
            phin[j] = get_neumann_phin(p_col, f_col);
          }
        }
      }
      return std::pair<std::vector<Complex>, std::vector<Complex>>{std::move(phi), std::move(phin)};
    };

    const double fmm_phin_scale =
        (bvp.last_fmm_matvec_stats.coordinate_scaling && bvp.last_fmm_matvec_stats.coordinate_scale_factor > 1e-10)
            ? bvp.last_fmm_matvec_stats.coordinate_scale_factor
            : 1.0;
    auto set_fmm_component = [&](const std::vector<Complex> &phi, const std::vector<Complex> &phin, const bool imag_part) {
      for (std::size_t j = 0; j < n; ++j) {
        auto [p_col, f_col] = id_by_index[j];
        auto &d = p_col->dof(const_cast<networkFace *>(f_col));
        d.phi_FMM = imag_part ? phi[j].imag() : phi[j].real();
        d.phin_FMM = fmm_phin_scale * (imag_part ? phin[j].imag() : phin[j].real());
        if (d.index >= 0 && static_cast<std::size_t>(d.index) < bvp.cache_phi_val_D_by_index.size()) {
          bvp.cache_phi_val_D_by_index[static_cast<std::size_t>(d.index)] = d.phi_FMM;
          bvp.cache_phin_val_D_by_index[static_cast<std::size_t>(d.index)] = d.phin_FMM;
        }
      }
    };

    auto eval_bie_fields = [&](const std::vector<Complex> &phi, const std::vector<Complex> &phin) {
      set_fmm_component(phi, phin, false);
      const auto yr = bvp.compute_Ax_minus_b(false);
      set_fmm_component(phi, phin, true);
      const auto yi = bvp.compute_Ax_minus_b(false);
      if (yr.size() != n || yi.size() != n)
        throw std::runtime_error("solve_linear_bvp: fmm matvec size mismatch");
      std::vector<Complex> y(n, Complex{0.0, 0.0});
      for (std::size_t i = 0; i < n; ++i)
        y[i] = Complex{yr[i], yi[i]};
      return y;
    };

    matrix_free_matvec = [&](const std::vector<Complex> &u) {
      if (u.size() != n)
        throw std::runtime_error("solve_linear_bvp: fmm_gmres vector size mismatch");
      auto [phi, phin] = reconstruct_fields(u, true);
      auto y = eval_bie_fields(phi, phin);
      for (std::size_t i = 0; i < n; ++i) {
        auto [p_row, f_row] = id_by_index[i];
        if (!is_interface_constraint_row(p_row, f_row))
          continue;
        const double constraint_scale = interface_row_scale(p_row, f_row);
        y[i] = constraint_scale * u[i];
        const auto *d_dir = p_row->findActiveBieDof(nullptr);
        if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
          throw std::runtime_error("solve_linear_bvp: missing Dirichlet ID at BCInterface point (fmm)");
        const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);
        if (node_is_robin[p_row]) {
          const Complex kappa = get_kappa(p_row, f_row);
          if (std::abs(kappa) < 1e-15)
            throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin BCInterface point (fmm)");
          y[i] += (-constraint_scale / kappa) * u[j_dir];
        }
      }
      return y;
    };

    {
      std::vector<Complex> zero(n, Complex{0.0, 0.0});
      auto [phi_known, phin_known] = reconstruct_fields(zero, false);
      matrix_free_rhs = eval_bie_fields(phi_known, phin_known);
      for (auto &v : matrix_free_rhs)
        v = -v;
      for (std::size_t i = 0; i < n; ++i) {
        auto [p_row, f_row] = id_by_index[i];
        if (!is_interface_constraint_row(p_row, f_row))
          continue;
        const double constraint_scale = interface_row_scale(p_row, f_row);
        if (node_is_robin[p_row]) {
          const Complex kappa = get_kappa(p_row, f_row);
          if (std::abs(kappa) < 1e-15)
            throw std::runtime_error("solve_linear_bvp: |kappa| too small at a Robin BCInterface point (fmm rhs)");
          matrix_free_rhs[i] = -constraint_scale * get_robin_rhs(p_row, f_row) / kappa;
        } else {
          matrix_free_rhs[i] = constraint_scale * get_dirichlet_phi(p_row);
        }
      }
    }
    const double matrix_free_rhs_l2 = l2_norm(matrix_free_rhs);

	    std::vector<Complex> left_scale(n, Complex{1.0, 0.0});
	    std::string effective_preconditioner = "none";
	    if (bc.gmres_preconditioner == "diagonal" || bc.gmres_preconditioner == "fmm_near_diagonal") {
	      if (bc.gmres_preconditioner == "diagonal" && build_dense_fmm_reference) {
	        for (std::size_t i = 0; i < n; ++i) {
	          const Complex d = dense_fmm_reference_A[aidx(i, i)];
	          if (std::abs(d) > 1e-300)
	            left_scale[i] = Complex{1.0, 0.0} / d;
	        }
	        effective_preconditioner = "dense_reference_diagonal";
	      } else {
	        auto near_diagonal_for = [&](const BEM_DOF_Base *node, const std::size_t row) {
	          const auto *target = dynamic_cast<const target4FMM *>(node);
	          if (!target)
	            return std::array<double, 2>{0.0, 0.0};
	          double IG = 0.0;
	          double IGn_without_jump = 0.0;
	          for (std::size_t k = 0; k < target->near_indices.size(); ++k) {
	            if (target->near_indices[k] != static_cast<int32_t>(row))
	              continue;
	            IG += target->near_weights_phi[k];
	            IGn_without_jump += target->near_weights_phin[k];
	          }
	          return std::array<double, 2>{IG, IGn_without_jump};
	        };
	        for (std::size_t i = 0; i < n; ++i) {
	          auto [p_row, f_row] = id_by_index[i];
	          Complex diag{1.0, 0.0};
	          if (is_interface_constraint_row(p_row, f_row)) {
	            diag = Complex{interface_row_scale(p_row, f_row), 0.0};
	          } else {
	            const auto [IG, IGn_without_jump] = near_diagonal_for(p_row, i);
	            const double jump = (i < bvp.diag_coeffs.size()) ? bvp.diag_coeffs[i] : 0.0;
	            const double IGn = IGn_without_jump + jump;
	            if (isDirichletBieDofKey(p_row, f_row)) {
	              if (is_robin_id(p_row, f_row)) {
	                const Complex kappa = get_kappa(p_row, f_row);
	                diag = Complex{IG, 0.0} - Complex{IGn, 0.0} / kappa;
	              } else {
	                diag = Complex{IG, 0.0};
	              }
	            } else {
              if (bc.robin_unknown_is_phi() && is_robin_id(p_row, f_row)) {
	                const Complex kappa = get_kappa(p_row, f_row);
	                diag = kappa * Complex{IG, 0.0} - Complex{IGn, 0.0};
	              } else {
	                diag = Complex{-IGn, 0.0};
	              }
	            }
	          }
	          if (std::abs(diag) > 1e-300)
	            left_scale[i] = Complex{1.0, 0.0} / diag;
	        }
	        effective_preconditioner = "fmm_near_diagonal";
	      }
	    } else if (bc.gmres_preconditioner == "row_norm") {
	      constexpr int probe_count = 8;
	      std::vector<double> row_norm2(n, 0.0);
      for (int probe = 0; probe < probe_count; ++probe) {
        std::vector<Complex> x_probe(n, Complex{0.0, 0.0});
        for (std::size_t i = 0; i < n; ++i) {
          // Deterministic Rademacher probes.  For random +/-1 vectors,
          // E |(A r)_i|^2 equals the squared row norm; this approximates a
          // left row scaling without assembling A.
          const std::uint64_t h = (static_cast<std::uint64_t>(i) + 0x9e3779b97f4a7c15ULL) ^
                                  (static_cast<std::uint64_t>(probe) * 0xbf58476d1ce4e5b9ULL);
          x_probe[i] = (static_cast<int>((h ^ (h >> 31)) & 1ULL) ? Complex{1.0, 0.0} : Complex{-1.0, 0.0});
        }
        const auto y_probe = matrix_free_matvec(x_probe);
        if (y_probe.size() != n)
          throw std::runtime_error("solve_linear_bvp: fmm_gmres row_norm probe size mismatch");
        for (std::size_t i = 0; i < n; ++i)
          row_norm2[i] += std::norm(y_probe[i]);
      }
      for (std::size_t i = 0; i < n; ++i) {
        const double estimate = std::sqrt(row_norm2[i] / static_cast<double>(probe_count));
        if (estimate > 1e-300 && std::isfinite(estimate))
          left_scale[i] = Complex{1.0 / estimate, 0.0};
      }
      effective_preconditioner = "row_norm_estimate";
    } else if (!bc.gmres_preconditioner.empty() && bc.gmres_preconditioner != "none") {
      std::cout << "[complex_fmm_gmres] preconditioner=" << bc.gmres_preconditioner
                << " is not available for matrix-free v1; using none" << std::endl;
    }
    const auto gm = solve_complex_gmres_operator(matrix_free_matvec,
                                                matrix_free_rhs,
                                                bc.gmres_restart,
                                                bc.gmres_max_iter,
                                                bc.gmres_tol,
                                                left_scale);
    solver_iterations = gm.iterations;
    solver_residual_l2 = gm.residual_l2;
    solver_relative_residual_l2 = (matrix_free_rhs_l2 > 1e-300) ? solver_residual_l2 / matrix_free_rhs_l2 : solver_residual_l2;
    solver_converged = gm.converged;
	    solver_effective_preconditioner = effective_preconditioner;
	    std::cout << "[complex_fmm_gmres] n=" << n
	              << " preconditioner=" << effective_preconditioner
	              << " iter=" << gm.iterations
	              << " residual=" << gm.residual_l2
	              << " relative_residual=" << solver_relative_residual_l2
	              << " converged=" << (gm.converged ? 1 : 0) << std::endl;
	    u_full = gm.x;
	    if (bc.debug_fmm_matvec_compare) {
	      if (dense_fmm_reference_A.empty() || dense_fmm_reference_rhs.empty())
	        throw std::runtime_error("solve_linear_bvp: missing dense reference for fmm matvec comparison");
	      bc.debug_fmm_matvec_compare(dense_fmm_reference_A,
	                                  dense_fmm_reference_rhs,
	                                  matrix_free_matvec,
	                                  matrix_free_rhs,
	                                  u_full,
	                                  id_by_index,
	                                  robin_faces,
	                                  node_is_robin);
	    }
	  } else if (bc.eliminate_interface_constraints) {
    std::vector<CondensedInterfaceExpr> eliminated(n);
    std::vector<int> retained_index(n, -1);
    std::vector<Id> retained_ids;
    retained_ids.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
      auto [p_row, f_row] = id_by_index[i];
      const bool row_is_corner_neumann =
          (p_row && p_row->BCInterface && isNeumannBieDofKey(p_row, f_row) &&
           !(bc.use_bie_row_for_interface && bc.use_bie_row_for_interface(*p_row, f_row)));
      if (row_is_corner_neumann) {
        const auto *d_dir = p_row->findActiveBieDof(nullptr);
        if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
          throw std::runtime_error("solve_linear_bvp: missing Dirichlet ID at eliminated BCInterface point");
        const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);
        const bool robin = node_is_robin[p_row];
        eliminated[i].eliminated = true;
        eliminated[i].dir_index = j_dir;
        if (robin) {
          const Complex kappa = get_kappa(p_row, f_row);
          if (std::abs(kappa) < 1e-15)
            throw std::runtime_error("solve_linear_bvp: |kappa| too small at an eliminated Robin BCInterface point");
          eliminated[i].scale = Complex{1.0, 0.0} / kappa;
          eliminated[i].constant = -get_robin_rhs(p_row, f_row) / kappa;
        } else {
          eliminated[i].constant = get_dirichlet_phi(p_row);
        }
        continue;
      }
      retained_index[i] = static_cast<int>(retained_ids.size());
      retained_ids.push_back(id_by_index[i]);
    }

    const std::size_t m = retained_ids.size();
    if (m == 0)
      throw std::runtime_error("solve_linear_bvp: condensed matrix size is 0");
    for (std::size_t i = 0; i < n; ++i) {
      if (!eliminated[i].eliminated)
        continue;
      if (retained_index[eliminated[i].dir_index] < 0)
        throw std::runtime_error("solve_linear_bvp: eliminated interface constraint points to another eliminated ID");
    }

    std::vector<Complex> A_condensed(m * m, Complex{0.0, 0.0});
    std::vector<Complex> b_condensed(m, Complex{0.0, 0.0});
    auto cidx = [&](std::size_t row, std::size_t col) -> std::size_t { return row + col * m; };

    for (std::size_t r = 0; r < n; ++r) {
      const int rr_i = retained_index[r];
      if (rr_i < 0)
        continue;
      const std::size_t rr = static_cast<std::size_t>(rr_i);
      Complex rhs = b[r];
      for (std::size_t c = 0; c < n; ++c) {
        const Complex value = A[aidx(r, c)];
        if (std::abs(value) == 0.0)
          continue;
        if (eliminated[c].eliminated) {
          const auto &expr = eliminated[c];
          rhs -= value * expr.constant;
          if (std::abs(expr.scale) > 0.0) {
            const std::size_t cc = static_cast<std::size_t>(retained_index[expr.dir_index]);
            A_condensed[cidx(rr, cc)] += value * expr.scale;
          }
        } else {
          const std::size_t cc = static_cast<std::size_t>(retained_index[c]);
          A_condensed[cidx(rr, cc)] += value;
        }
      }
      b_condensed[rr] = rhs;
    }

    const std::vector<Complex> rhs_condensed_original =
        bc.debug_condensed_solved_linear_system ? b_condensed : std::vector<Complex>{};
    if (bc.linear_solver == "gmres") {
      const double rhs_l2 = l2_norm(b_condensed);
      const auto gm = solve_dense_complex_gmres(A_condensed,
                                                b_condensed,
                                                bc.gmres_restart,
                                                bc.gmres_max_iter,
                                                bc.gmres_tol,
                                                bc.gmres_preconditioner);
      b_condensed = gm.x;
      solver_iterations = gm.iterations;
      solver_residual_l2 = gm.residual_l2;
      solver_relative_residual_l2 = (rhs_l2 > 1e-300) ? solver_residual_l2 / rhs_l2 : solver_residual_l2;
      solver_converged = gm.converged;
      solver_effective_preconditioner = bc.gmres_preconditioner;
      std::cout << "[complex_gmres] condensed n=" << m
                << " preconditioner=" << bc.gmres_preconditioner
                << " iter=" << gm.iterations
                << " residual=" << gm.residual_l2
                << " relative_residual=" << solver_relative_residual_l2
                << " converged=" << (gm.converged ? 1 : 0) << std::endl;
    } else if (bc.linear_solver == "lu" || bc.linear_solver.empty()) {
      lapack_zlu lu(static_cast<int>(m), A_condensed);
      lu.solve_in_place(b_condensed);
      solver_relative_residual_l2 = 0.0;
      solver_effective_preconditioner = "none";
    } else {
      throw std::runtime_error("solve_linear_bvp: unsupported linear_solver=" + bc.linear_solver);
    }

    if (bc.debug_condensed_solved_linear_system) {
      bc.debug_condensed_solved_linear_system(A_condensed,
                                             rhs_condensed_original,
                                             b_condensed,
                                             id_by_index,
                                             retained_ids,
                                             retained_index,
                                             eliminated,
                                             robin_faces,
                                             node_is_robin);
    }

    u_full.assign(n, Complex{0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i) {
      const int ii = retained_index[i];
      if (ii >= 0)
        u_full[i] = b_condensed[static_cast<std::size_t>(ii)];
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (!eliminated[i].eliminated)
        continue;
      const auto &expr = eliminated[i];
      u_full[i] = expr.scale * u_full[expr.dir_index] + expr.constant;
    }
  } else {
    if (bc.linear_solver == "gmres") {
      const auto gm = solve_dense_complex_gmres(A,
                                                b,
                                                bc.gmres_restart,
                                                bc.gmres_max_iter,
                                                bc.gmres_tol,
                                                bc.gmres_preconditioner);
      solver_iterations = gm.iterations;
      solver_residual_l2 = gm.residual_l2;
      solver_relative_residual_l2 = (dense_rhs_l2 > 1e-300) ? solver_residual_l2 / dense_rhs_l2 : solver_residual_l2;
      solver_converged = gm.converged;
      solver_effective_preconditioner = bc.gmres_preconditioner;
      std::cout << "[complex_gmres] n=" << n
                << " preconditioner=" << bc.gmres_preconditioner
                << " iter=" << gm.iterations
                << " residual=" << gm.residual_l2
                << " relative_residual=" << solver_relative_residual_l2
                << " converged=" << (gm.converged ? 1 : 0) << std::endl;
      u_full = gm.x;
    } else if (bc.linear_solver == "lu" || bc.linear_solver.empty()) {
      lapack_zlu lu(static_cast<int>(n), A);
      lu.solve_in_place(b); // b becomes the solution u
      u_full = b;
      solver_relative_residual_l2 = 0.0;
      solver_effective_preconditioner = "none";
    } else {
      throw std::runtime_error("solve_linear_bvp: unsupported linear_solver=" + bc.linear_solver);
    }
  }

  Solution sol;
  sol.n = n;
  sol.id_by_index = id_by_index;
  sol.u = std::move(u_full);
  sol.interface_constraint_row.assign(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    auto [p_row, f_row] = id_by_index[i];
    sol.interface_constraint_row[i] = is_interface_constraint_row(p_row, f_row) ? 1 : 0;
  }
  sol.linear_solver = (bc.linear_solver.empty() ? "lu" : bc.linear_solver);
  sol.linear_solver_iterations = solver_iterations;
  sol.linear_solver_residual_l2 = solver_residual_l2;
  sol.linear_solver_relative_residual_l2 = solver_relative_residual_l2;
  sol.linear_solver_converged = solver_converged;
  sol.linear_solver_effective_preconditioner = solver_effective_preconditioner;
  if (use_matrix_free_fmm) {
    const auto &stats = bvp.last_fmm_matvec_stats;
    sol.fmm_setup_source = stats.setup_source;
    sol.fmm_coordinate_scaling = stats.coordinate_scaling;
    sol.fmm_morton_reindex = stats.morton_reindex;
    sol.fmm_reused_sources = stats.reused_sources;
    sol.fmm_reused_static = stats.reused_static_fmm;
    sol.fmm_targets = stats.targets;
    sol.fmm_sources = stats.sources;
    sol.fmm_total_near_terms = stats.total_near_terms;
    sol.fmm_mean_near_terms_per_target = stats.mean_near_terms_per_target;
    sol.fmm_max_near_terms_per_target = stats.max_near_terms_per_target;
    sol.fmm_vertex_targets = stats.vertex_targets;
    sol.fmm_vertex_total_near_terms = stats.vertex_total_near_terms;
    sol.fmm_vertex_mean_near_terms_per_target = stats.vertex_mean_near_terms_per_target;
    sol.fmm_vertex_max_near_terms_per_target = stats.vertex_max_near_terms_per_target;
    sol.fmm_midpoint_targets = stats.midpoint_targets;
    sol.fmm_midpoint_total_near_terms = stats.midpoint_total_near_terms;
    sol.fmm_midpoint_mean_near_terms_per_target = stats.midpoint_mean_near_terms_per_target;
    sol.fmm_midpoint_max_near_terms_per_target = stats.midpoint_max_near_terms_per_target;
    sol.fmm_max_source_offset = stats.max_source_offset;
    sol.fmm_mean_source_offset = stats.mean_source_offset;
    sol.fmm_p95_source_offset = stats.p95_source_offset;
    sol.fmm_coordinate_scale_factor = stats.coordinate_scale_factor;
    sol.fmm_top_source_offsets = stats.top_source_offsets;
  }
  if (use_matrix_free_fmm && fmm_session) {
    // Matrix-free matvec evaluates known Neumann data while coordinates are
    // scaled and applies the matching FMM phin scale internally.  The public
    // Solution fields, however, must be reconstructed in physical coordinates.
    fmm_session.reset();
  }
  sol.phi.assign(n, Complex{0.0, 0.0});
  sol.phin.assign(n, Complex{0.0, 0.0});

  // Reconstruct (phi, phin) on each ID.
  for (std::size_t j = 0; j < n; ++j) {
    auto [p_col, f_col] = sol.id_by_index[j];
    if (isDirichletBieDofKey(p_col, f_col)) {
      sol.phin[j] = sol.u[j];
      const bool robin = is_robin_id(p_col, f_col);
      if (robin) {
        if (bc.robin_unknown_is_phi())
          throw std::runtime_error("solve_linear_bvp: Robin Dirichlet-type reconstruction with robin_unknown_phi=true");
        const Complex kappa = get_kappa(p_col, f_col);
        sol.phi[j] = (sol.u[j] - get_robin_rhs(p_col, f_col)) / kappa;
      } else {
        sol.phi[j] = get_dirichlet_phi(p_col);
      }
    } else {
      sol.phi[j] = sol.u[j];
      if (bc.robin_unknown_is_phi() && is_robin_id(p_col, f_col)) {
        const Complex kappa = get_kappa(p_col, f_col);
        sol.phin[j] = kappa * sol.u[j] + get_robin_rhs(p_col, f_col);
      } else {
        sol.phin[j] = get_neumann_phin(p_col, f_col);
      }
    }
  }

  // Full BIE residual check (for diagnostics).
  sol.bie_residual.assign(n, Complex{0.0, 0.0});
  if (use_matrix_free_fmm) {
    // The matrix-free GMRES result already reports the true unpreconditioned
    // ||b - A u|| from the same FMM matvec.  Avoid a second post-solve FMM
    // pass here; it is redundant and keeps v1 focused on the solver path.
    sol.bie_residual_l2 = solver_residual_l2;
  } else {
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
  }

  if (bc.compute_adjoint_residual) {
    sol.adjoint_residual.assign(n, Complex{0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i) {
      Complex atu{0.0, 0.0};
      for (std::size_t j = 0; j < n; ++j)
        atu += A[aidx(j, i)] * sol.u[j];
      sol.adjoint_residual[i] = atu - rhs_original[i];
    }
    sol.adjoint_residual_l2 = l2_norm(sol.adjoint_residual);
  }

  if (bc.debug_solved_linear_system)
    bc.debug_solved_linear_system(A, rhs_original, sol.u, sol.id_by_index, robin_faces, node_is_robin);

  // Restore original boundary flags/indexing.
  guard.restore();
  return sol;
}

} // namespace bem_frequency_domain
