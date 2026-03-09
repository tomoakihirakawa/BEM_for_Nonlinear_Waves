#ifndef TanakaSolitaryWave_HPP
#define TanakaSolitaryWave_HPP

#include <vector>
#include <cmath>
#include <array>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <cassert>

/*
 * TanakaSolitaryWave: Fully nonlinear solitary wave solution
 * using Tanaka's (1986) iterative method.
 *
 * Reference: M. Tanaka, "The stability of solitary waves,"
 *            Physics of Fluids 29(3), 650-655, 1986.
 *
 * Algorithm summary (Section II of the paper):
 *   1. Work in wave frame: h=1 (depth), c=1 (phase speed)
 *   2. Complex velocity potential W = Φ + iΨ in strip 0 < Ψ < 1
 *   3. Ω = ln(dW/dz) = τ - iθ, where q = exp(τ), θ = velocity angle
 *   4. Cauchy integral eq (3) relates τ and θ on free surface (Ψ=1)
 *   5. Bernoulli eq (1): dq³/dΦ = -(3/F²) sin θ
 *   6. Iterate: assume τ → solve for θ → compute F² → update q → τ = ln q
 *   7. Variable transform Φ = αγ + γ^m concentrates mesh near crest
 */

struct TanakaSolitaryWave {

   // ===== Physical parameters (dimensional) =====
   double depth = 0;
   double wave_height = 0;
   double phase_speed = 0;
   double x_offset = 0;
   double bottom_z = 0;
   double L = 0; // effective wavelength (for absorber gamma)

   // ===== Non-dimensional solver output =====
   double F2 = 0;
   double qc = 0;

   // ===== Tabulated solution on uniform x-grid (non-dim, |ξ|/h) =====
   std::vector<double> x_tab;   // ≥ 0
   std::vector<double> eta_tab; // surface elevation (even)
   std::vector<double> phi_tab; // lab-frame potential on surface (stored for ξ ≥ 0, odd extension)
   std::vector<double> us_tab;  // lab-frame u on surface = c(1 - q cos θ) / c (even)
   std::vector<double> ws_tab;  // |w| on surface = c|q sin θ| / c (even, ≥ 0)
   double dx_tab = 0;
   int n_tab = 0;
   double phi_asymptotic = 0; // asymptotic value of phi_nd as |ξ| → ∞

   // =========================================================================
   //  Public interface (compatible with WaterWaveTheory pattern)
   // =========================================================================

   /**
    * Solve for the solitary wave.
    * @param H_over_h  Wave height / depth ratio (0 < H/h < ~0.83)
    * @param h         Water depth [m]
    * @param bz        Bottom z-coordinate [m]
    * @param x0        Initial crest x-position [m]
    * @param N_solver  Mesh points for Tanaka solver (default 100)
    */
   void solve(double H_over_h, double h, double bz = 0, double x0 = 0, int N_solver = 100) {
      if (H_over_h <= 0 || H_over_h >= 0.84)
         throw std::runtime_error("TanakaSolitaryWave: H/h must be in (0, ~0.83)");
      depth = h;
      bottom_z = bz;
      x_offset = x0;

      solve_for_qc(H_over_h, N_solver);

      phase_speed = std::sqrt(F2 * _GRAVITY_ * depth);
      wave_height = H_over_h * depth;

      // Effective wavelength: width where η > 0.01 * η_max
      double threshold = 0.01;
      L = 0;
      for (int i = 0; i < n_tab; ++i) {
         if (eta_tab[i] > threshold * eta_tab[0])
            L = x_tab[i] * depth * 2.0;
      }
      if (L < depth) L = 10.0 * depth;

      std::cout << "TanakaSolitaryWave solved: H/h=" << H_over_h
                << " F²=" << F2 << " c=" << phase_speed
                << " L_eff=" << L << " m"
                << " phi_asymptotic=" << phi_asymptotic * phase_speed * depth
                << std::endl;
   }

   /** Free surface elevation (returns z-coordinate). */
   double eta(const std::array<double, 3> &X, double t) const {
      double xi_nd = compute_xi_nd(X, t);
      double eta_nd = interp_even(std::abs(xi_nd));
      return eta_nd * depth + depth + bottom_z;
   }

   /** Velocity potential at (X, t). */
   double phi(const std::array<double, 3> &X, double t) const {
      double xi_nd = compute_xi_nd(X, t);
      double phi_nd = interp_phi_odd(xi_nd);
      return phi_nd * phase_speed * depth;
   }

   /** Velocity (u, v, w) at (X, t). */
   std::array<double, 3> gradPhi(const std::array<double, 3> &X, double t) const {
      double xi_nd = compute_xi_nd(X, t);
      double abs_xi = std::abs(xi_nd);
      double y_nd = (X[2] - bottom_z) / depth;

      double u_s = interp_tab(abs_xi, us_tab);
      double w_s = interp_tab(abs_xi, ws_tab);
      double eta_nd = interp_tab(abs_xi, eta_tab);

      double total_d = 1.0 + eta_nd;
      double y_ratio = (total_d > 1e-10) ? std::clamp(y_nd / total_d, 0.0, 1.0) : 0.0;

      double w_sign = (xi_nd > 0) ? 1.0 : ((xi_nd < 0) ? -1.0 : 0.0);

      return {u_s * phase_speed,
              0.0,
              w_s * w_sign * y_ratio * phase_speed};
   }

   /** Time derivative of velocity: ∂(∇φ)/∂t = -c ∂(∇φ)/∂x */
   std::array<double, 3> gradPhi_t(const std::array<double, 3> &X, double t) const {
      double eps = depth * 1e-5;
      auto Xp = X, Xm = X;
      Xp[0] += eps;
      Xm[0] -= eps;
      auto vp = gradPhi(Xp, t);
      auto vm = gradPhi(Xm, t);
      return {-phase_speed * (vp[0] - vm[0]) / (2 * eps),
              0.0,
              -phase_speed * (vp[2] - vm[2]) / (2 * eps)};
   }

private:
   // =========================================================================
   //  Coordinate helpers
   // =========================================================================

   double compute_xi_nd(const std::array<double, 3> &X, double t) const {
      return (X[0] - x_offset - phase_speed * t) / depth;
   }

   // =========================================================================
   //  Interpolation on uniform x_tab grid
   // =========================================================================

   double interp_tab(double abs_x, const std::vector<double> &tab) const {
      if (abs_x <= 0) return tab[0];
      if (dx_tab <= 0) return 0;
      double fi = abs_x / dx_tab;
      int i = static_cast<int>(fi);
      if (i >= n_tab - 1) return 0;
      double frac = fi - i;
      return tab[i] * (1.0 - frac) + tab[i + 1] * frac;
   }

   double interp_even(double abs_x) const {
      return interp_tab(abs_x, eta_tab);
   }

   double interp_phi_odd(double xi_nd) const {
      double abs_xi = std::abs(xi_nd);
      double val;
      if (n_tab > 0 && abs_xi >= x_tab[n_tab - 1])
         val = phi_asymptotic; // use asymptotic value beyond table
      else
         val = interp_tab(abs_xi, phi_tab);
      return (xi_nd >= 0) ? val : -val;
   }

   // =========================================================================
   //  Lagrange interpolation (8-point) for shifted mesh points
   // =========================================================================

   static double lagrange8(const std::vector<double> &f, int N, double idx_frac) {
      int i0 = static_cast<int>(std::floor(idx_frac)) - 3;
      double s = idx_frac - std::floor(idx_frac) + 3.0; // local coord, s ∈ [3, 4)
      double result = 0;
      for (int k = 0; k < 8; ++k) {
         int ii = i0 + k;
         double fval = 0;
         if (ii >= 0 && ii < N) fval = f[ii];
         double wt = 1.0;
         for (int j = 0; j < 8; ++j) {
            if (j != k) wt *= (s - j) / (k - j);
         }
         result += fval * wt;
      }
      return result;
   }

   // =========================================================================
   //  Gauss elimination for dense linear system Ax = b
   // =========================================================================

   static void solve_dense(std::vector<std::vector<double>> &A,
                           std::vector<double> &b, int n) {
      for (int col = 0; col < n; ++col) {
         // Partial pivoting
         int pivot = col;
         for (int row = col + 1; row < n; ++row)
            if (std::abs(A[row][col]) > std::abs(A[pivot][col])) pivot = row;
         if (pivot != col) {
            std::swap(A[col], A[pivot]);
            std::swap(b[col], b[pivot]);
         }
         double diag = A[col][col];
         if (std::abs(diag) < 1e-30)
            throw std::runtime_error("Tanaka solver: singular matrix");
         for (int row = col + 1; row < n; ++row) {
            double factor = A[row][col] / diag;
            for (int j = col + 1; j < n; ++j) A[row][j] -= factor * A[col][j];
            b[row] -= factor * b[col];
         }
      }
      // Back substitution
      for (int row = n - 1; row >= 0; --row) {
         for (int j = row + 1; j < n; ++j) b[row] -= A[row][j] * b[j];
         b[row] /= A[row][row];
      }
   }

   // =========================================================================
   //  Core Tanaka solver
   // =========================================================================

   struct SolverResult {
      std::vector<double> Phi;     // potential values at mesh points (j=0..N-1)
      std::vector<double> tau;     // ln(q)
      std::vector<double> theta;   // velocity angle
      std::vector<double> x_phys;  // physical x from crest
      std::vector<double> y_phys;  // physical y (non-dim)
      double F2;
      int N;
   };

   SolverResult run_tanaka(double qc_in, int N) {
      // Solver parameters
      const double alpha = 0.01 + 0.1 * qc_in;
      const int m_exp = 5;
      const double gamma_max = std::pow(40.0, 1.0 / m_exp); // Φ_max ≈ 40
      const double Delta = gamma_max / (N - 1);

      // Build mesh
      std::vector<double> gamma_pts(N), Phi_pts(N), dPhi_dg(N);
      for (int j = 0; j < N; ++j) {
         double g = j * Delta;
         gamma_pts[j] = g;
         Phi_pts[j] = alpha * g + std::pow(g, m_exp);
         dPhi_dg[j] = alpha + m_exp * std::pow(g, m_exp - 1);
      }

      // Quadrature weights (trapezoidal rule in γ)
      std::vector<double> w(N);
      for (int j = 0; j < N; ++j) {
         w[j] = dPhi_dg[j] * Delta;
         if (j == 0 || j == N - 1) w[j] *= 0.5;
      }

      // Shifted mesh points for PV integral
      int N_shift = N - 1;
      std::vector<double> gamma_s(N_shift), Phi_s(N_shift), dPhi_s(N_shift), w_s(N_shift);
      for (int k = 0; k < N_shift; ++k) {
         double g = (k + 0.5) * Delta;
         gamma_s[k] = g;
         Phi_s[k] = alpha * g + std::pow(g, m_exp);
         dPhi_s[k] = alpha + m_exp * std::pow(g, m_exp - 1);
         w_s[k] = dPhi_s[k] * Delta;
      }

      // Initial guess: τ decays like sech²
      double tau_c = std::log(qc_in);
      std::vector<double> tau(N), theta(N, 0.0);
      double decay_scale = 5.0;
      for (int j = 0; j < N; ++j) {
         double Phi_j = Phi_pts[j];
         double sech = 1.0 / std::cosh(Phi_j / decay_scale);
         tau[j] = tau_c * sech * sech;
      }

      double F2_val = 1.0 + 0.5; // initial guess

      // ===== Outer iteration =====
      for (int outer = 0; outer < 40; ++outer) {
         // Step 1: Solve integral equation (3) for θ given τ
         compute_theta(Phi_pts, tau, theta, w, Phi_s, w_s, gamma_pts, Delta, N, N_shift);

         // Step 2: Compute F² from eq (5)
         double S = 0;
         for (int j = 0; j < N; ++j)
            S += std::sin(theta[j]) * w[j];
         double F2_new = -3.0 * S / (1.0 - qc_in * qc_in * qc_in);

         // Step 3: Update q from eq (4)
         std::vector<double> tau_new(N);
         tau_new[0] = std::log(qc_in);
         double integral = 0;
         for (int j = 1; j < N; ++j) {
            // Trapezoidal integration of sin θ from 0 to Φ_j
            integral += 0.5 * (std::sin(theta[j - 1]) * w[j - 1] + std::sin(theta[j]) * w[j]);
            // Actually: cumulative trap from j-1 to j
            // But weights already include dΦ, so:
            // ∫_0^{Φ_j} sin θ dΦ ≈ Σ_{k=0}^{j} sin(θ_k) w_k  (with endpoint correction)
            // Let me use a proper cumulative sum
         }
         // Redo: proper cumulative trapezoidal integration
         std::vector<double> cum_sin(N, 0);
         for (int j = 1; j < N; ++j) {
            // Integral from Φ_{j-1} to Φ_j using trapezoidal rule
            double dPhi = Phi_pts[j] - Phi_pts[j - 1];
            cum_sin[j] = cum_sin[j - 1] + 0.5 * (std::sin(theta[j - 1]) + std::sin(theta[j])) * dPhi;
         }
         for (int j = 0; j < N; ++j) {
            double q3 = qc_in * qc_in * qc_in - (3.0 / F2_new) * cum_sin[j];
            double q_val = std::cbrt(std::max(q3, 1e-30));
            tau_new[j] = std::log(q_val);
         }

         // Convergence check
         double dF2 = std::abs(F2_new - F2_val);
         F2_val = F2_new;
         tau = tau_new;

         if (dF2 < 1e-12 && outer > 3) break;
      }

      // ===== Compute physical profile =====
      SolverResult res;
      res.N = N;
      res.F2 = F2_val;
      res.Phi = Phi_pts;
      res.tau = tau;
      res.theta = theta;
      res.x_phys.resize(N);
      res.y_phys.resize(N);

      // Integrate x and y from crest (Φ=0)
      res.x_phys[0] = 0;
      res.y_phys[0] = 0; // will add y_crest later
      for (int j = 1; j < N; ++j) {
         double dPhi = Phi_pts[j] - Phi_pts[j - 1];
         double e_tau_0 = std::exp(-tau[j - 1]);
         double e_tau_1 = std::exp(-tau[j]);
         res.x_phys[j] = res.x_phys[j - 1] + 0.5 * dPhi *
                              (e_tau_0 * std::cos(theta[j - 1]) + e_tau_1 * std::cos(theta[j]));
         res.y_phys[j] = res.y_phys[j - 1] + 0.5 * dPhi *
                              (e_tau_0 * std::sin(theta[j - 1]) + e_tau_1 * std::sin(theta[j]));
      }
      // y_crest: at Φ→∞, y should be 1.0 (still water)
      // y(∞) = y_crest_raw + y_phys[last] = 1.0
      // But y_phys[last] includes the drop from crest to far field
      // Actually y_phys[j] = ∫_0^{Φ_j} exp(-τ) sin θ dΦ (cumulative from crest)
      // y(Φ_j) = y_crest + y_phys[j]
      // y(∞) = 1.0, so y_crest = 1.0 - y_phys[last]
      double y_crest = 1.0 - res.y_phys[N - 1];
      for (int j = 0; j < N; ++j)
         res.y_phys[j] = y_crest + res.y_phys[j]; // now y_phys = actual y coordinate

      return res;
   }

   // =========================================================================
   //  Integral equation solver: compute θ given τ
   // =========================================================================

   void compute_theta(const std::vector<double> &Phi,
                      const std::vector<double> &tau,
                      std::vector<double> &theta,
                      const std::vector<double> &w,
                      const std::vector<double> &Phi_s,
                      const std::vector<double> &w_s,
                      const std::vector<double> &gamma_pts,
                      double Delta, int N, int N_shift) {
      // Integral equation (using symmetry, Φ > 0):
      // θ_i + (2/π) Σ_j θ_j K_θ(Φ_j, Φ_i) w_j
      //   = (1/π) Σ_j τ_j K_τ(Φ_j, Φ_i) w_j
      //     - (1/π) Σ_k τ(Φ'_k) × 2Φ_i/(Φ'_k² - Φ_i²) × w'_k
      //
      // K_θ(φ, Φ) = 1/((φ-Φ)²+4) - 1/((φ+Φ)²+4)
      // K_τ(φ, Φ) = (φ-Φ)/((φ-Φ)²+4) - (φ+Φ)/((φ+Φ)²+4)

      theta[0] = 0; // symmetry
      int M = N - 1; // unknowns: θ_1, ..., θ_{N-1}

      // Interpolate τ at shifted points
      std::vector<double> tau_shift(N_shift);
      for (int k = 0; k < N_shift; ++k) {
         tau_shift[k] = lagrange8(tau, N, k + 0.5);
      }

      // Build matrix A and RHS b
      std::vector<std::vector<double>> A(M, std::vector<double>(M, 0));
      std::vector<double> b(M, 0);

      for (int ii = 0; ii < M; ++ii) {
         int i = ii + 1; // Φ_i = Phi[i]
         double Pi = Phi[i];

         // Identity
         A[ii][ii] = 1.0;

         // Kernel terms for θ (j = 1, ..., N-1)
         for (int jj = 0; jj < M; ++jj) {
            int j = jj + 1;
            double Pj = Phi[j];
            double d1 = Pj - Pi, d2 = Pj + Pi;
            double K_theta = 1.0 / (d1 * d1 + 4.0) - 1.0 / (d2 * d2 + 4.0);
            A[ii][jj] += (2.0 / M_PI) * K_theta * w[j];
         }

         // RHS: τ kernel (j = 0, ..., N-1)
         double rhs_tau = 0;
         for (int j = 0; j < N; ++j) {
            double Pj = Phi[j];
            double d1 = Pj - Pi, d2 = Pj + Pi;
            double K_tau = d1 / (d1 * d1 + 4.0) - d2 / (d2 * d2 + 4.0);
            rhs_tau += tau[j] * K_tau * w[j];
         }
         rhs_tau /= M_PI;

         // RHS: PV integral using shifted points
         double rhs_pv = 0;
         for (int k = 0; k < N_shift; ++k) {
            double Pk = Phi_s[k];
            double denom = Pk * Pk - Pi * Pi;
            if (std::abs(denom) > 1e-30) {
               rhs_pv += tau_shift[k] * 2.0 * Pi / denom * w_s[k];
            }
         }
         rhs_pv /= M_PI;

         b[ii] = rhs_tau - rhs_pv;
      }

      // Solve linear system
      solve_dense(A, b, M);

      // Store result
      for (int ii = 0; ii < M; ++ii)
         theta[ii + 1] = b[ii];
   }

   // =========================================================================
   //  Bisection to find qc for desired H/h
   // =========================================================================

   void solve_for_qc(double H_over_h, int N_solver) {
      // H/h = η_crest = F²(1 - qc²)/2  approximately
      // H/h is monotonically decreasing with qc

      double qc_lo = 0.02;
      double qc_hi = 0.98;

      SolverResult best_result;
      double best_qc = 0.5;

      for (int iter = 0; iter < 50; ++iter) {
         double qc_mid = (qc_lo + qc_hi) / 2.0;
         auto result = run_tanaka(qc_mid, N_solver);

         // Actual wave height from profile
         double eta_crest = result.y_phys[0] - 1.0;

         if (std::abs(eta_crest - H_over_h) < 1e-8) {
            best_qc = qc_mid;
            best_result = result;
            break;
         }

         if (eta_crest > H_over_h) {
            qc_lo = qc_mid; // wave too high → increase qc (weaker wave)
         } else {
            qc_hi = qc_mid; // wave too low → decrease qc (stronger wave)
         }

         best_qc = qc_mid;
         best_result = result;
      }

      qc = best_qc;
      F2 = best_result.F2;
      build_tables(best_result);
   }

   // =========================================================================
   //  Build uniform x-grid tables from solver result
   // =========================================================================

   void build_tables(const SolverResult &res) {
      // Resample to uniform x-grid
      double x_max = res.x_phys[res.N - 1];
      if (x_max < 1.0) x_max = 50.0;

      n_tab = 2000;
      dx_tab = x_max / (n_tab - 1);
      x_tab.resize(n_tab);
      eta_tab.resize(n_tab);
      phi_tab.resize(n_tab);
      us_tab.resize(n_tab);
      ws_tab.resize(n_tab);

      // Precompute on the Tanaka mesh: lab-frame quantities
      int N = res.N;
      std::vector<double> eta_raw(N), phi_raw(N), us_raw(N), ws_raw(N);
      for (int j = 0; j < N; ++j) {
         eta_raw[j] = res.y_phys[j] - 1.0; // surface elevation
         double q = std::exp(res.tau[j]);
         double costh = std::cos(res.theta[j]);
         double sinth = std::sin(res.theta[j]);
         us_raw[j] = 1.0 - q * costh; // lab-frame u / c (right-moving)
         ws_raw[j] = std::abs(q * sinth); // |w| / c on surface
      }

      // Lab-frame potential on surface (right-moving):
      // φ_lab(ξ) = ξ - Φ_wave(ξ) for ξ > 0 (anti-symmetric)
      // Here ξ = x_phys[j] and Φ_wave = Phi[j]
      for (int j = 0; j < N; ++j)
         phi_raw[j] = res.x_phys[j] - res.Phi[j];

      // Asymptotic value: phi_raw → constant (≠0) as ξ → ∞
      // This is the Stokes drift contribution.
      phi_asymptotic = phi_raw[N - 1];

      // Resample using linear interpolation from non-uniform x_phys to uniform x_tab
      for (int i = 0; i < n_tab; ++i) {
         double x = i * dx_tab;
         x_tab[i] = x;

         // Find interval in x_phys (which is monotonically increasing)
         if (x <= 0) {
            eta_tab[i] = eta_raw[0];
            phi_tab[i] = phi_raw[0];
            us_tab[i] = us_raw[0];
            ws_tab[i] = ws_raw[0];
            continue;
         }
         if (x >= res.x_phys[N - 1]) {
            eta_tab[i] = 0;
            phi_tab[i] = phi_asymptotic; // non-zero asymptotic value
            us_tab[i] = 0;
            ws_tab[i] = 0;
            continue;
         }

         // Binary search
         int lo = 0, hi = N - 1;
         while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (res.x_phys[mid] <= x)
               lo = mid;
            else
               hi = mid;
         }
         double t = (x - res.x_phys[lo]) / (res.x_phys[hi] - res.x_phys[lo]);
         eta_tab[i] = eta_raw[lo] + t * (eta_raw[hi] - eta_raw[lo]);
         phi_tab[i] = phi_raw[lo] + t * (phi_raw[hi] - phi_raw[lo]);
         us_tab[i] = us_raw[lo] + t * (us_raw[hi] - us_raw[lo]);
         ws_tab[i] = ws_raw[lo] + t * (ws_raw[hi] - ws_raw[lo]);
      }
   }
};

#endif
