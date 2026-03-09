#pragma once

#include "Network.hpp"          // Networkクラスの定義(Node, Face等)が必要
#include "basic.hpp"            // Tddd, Norm, Cross, Dot などの基本的な関数定義を想定
#include "integrationOfODE.hpp" // RungeKutta
#include <algorithm>            // std::remove_if
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// M_PIの定義ガード
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief 3次元渦粒子 (Vortex Particle)
 * 支配方程式 (ラグランジュ記述):dx/dt = u
 * d(alpha)/dt = (alpha . nabla) u  (ストレッチング) + 粘性項(PSEなどで別途計算)
 */
struct VortexParticle {
  std::array<double, 3> x;     //!< 粒子位置
  std::array<double, 3> alpha; //!< 渦強度ベクトル (alpha = omega * volume)
  double sigma;                //!< コア半径 (正則化パラメータ)　ｙ
  double volume;               //!< 粒子体積

  // 時間発展計算用の一時変数
  std::array<double, 3> u_omega_VPM = {0., 0., 0.};     //!< 粒子による誘起速度
  std::array<double, 3> u_potential_BEM = {0., 0., 0.}; //!< ポテンシャル流れによる誘起速度
  std::array<double, 3> u_total = {0., 0., 0.};         //!< 粒子速度
  std::array<double, 3> d_alpha_dt = {0., 0., 0.};      //!< ストレッチングによる渦強度の変化率
  RungeKutta<std::array<double, 3>> RK_x;
  RungeKutta<std::array<double, 3>> RK_alpha;

  VortexParticle(const std::array<double, 3> &x_in, const std::array<double, 3> &alpha_in, double sigma_in, double vol_in) : x(x_in), alpha(alpha_in), sigma(sigma_in), volume(vol_in) {}
};

/**
 * @brief 渦粒子法 管理クラス
 */
class VortexMethod {
public:
  enum class StretchingScheme {
    Standard,  // dα/dt = (α·∇)u = (∇u) α
    Transpose, // dα/dt = (∇u)^T α   (Cottet & Koumoutsakos Eq. 3.1.7)
  };

  // PSE (diffusion) kernel correction level:
  // - None: original uncorrected PSE
  // - Gradient: enforce 1st moments (cancel gradient contamination)
  // - Curvature: enforce 1st + 2nd moments (match curvature / Laplacian moments)
  enum class PSECorrectionMode {
    None,
    Gradient,
    Curvature,
  };

private:
  double nu; // basic_constants.hppで定義されていると仮定
  StretchingScheme stretching_scheme = StretchingScheme::Standard;
  PSECorrectionMode pse_correction_mode = PSECorrectionMode::None;
  // For the current (Gaussian) PSE kernel in this file, the continuous second moment satisfies:
  //   ∫ η(r;σ) r_i^2 dV = 1   (i = x,y,z)
  // so the Laplacian (curvature) reproduction target becomes 2.0 (to cancel Taylor's 1/2).
  static constexpr double pse_kernel_second_moment = 1.0;
  // Only used when `pse_correction_mode == Curvature`.
  // For Laplacian (2nd-derivative) reproduction, the natural target is 2.0
  // (to cancel the Taylor-expansion 1/2 factor in the quadratic term).
  double pse_correction_second_moment_target = 2.0 * pse_kernel_second_moment;
  double pse_correction_regularization_rel = 1e-12; // diagonal Tikhonov regularization (relative to matrix scale)
  double pse_correction_max_dimless_coeff = 1e3;    // safety clamp for ill-conditioned neighborhoods
  std::vector<VortexParticle> particles;
  std::function<std::array<double, 3>(const std::array<double, 3> &)> potentialField;

  template <std::size_t N> static bool solveLinearSystemGaussian(std::array<std::array<double, N>, N> A, std::array<double, N> b, std::array<double, N> &x) {
    double max_abs = 0.0;
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        max_abs = std::max(max_abs, std::abs(A[i][j]));
    if (!(max_abs > 0.0) || !std::isfinite(max_abs))
      return false;
    const double eps = max_abs * 1e-14 + 1e-30;

    for (std::size_t k = 0; k < N; ++k) {
      std::size_t pivot = k;
      double pivot_abs = std::abs(A[k][k]);
      for (std::size_t i = k + 1; i < N; ++i) {
        const double v = std::abs(A[i][k]);
        if (v > pivot_abs) {
          pivot_abs = v;
          pivot = i;
        }
      }
      if (!(pivot_abs > eps))
        return false;

      if (pivot != k) {
        std::swap(A[pivot], A[k]);
        std::swap(b[pivot], b[k]);
      }

      const double akk = A[k][k];
      for (std::size_t i = k + 1; i < N; ++i) {
        const double factor = A[i][k] / akk;
        if (factor == 0.0)
          continue;
        A[i][k] = 0.0;
        for (std::size_t j = k + 1; j < N; ++j)
          A[i][j] -= factor * A[k][j];
        b[i] -= factor * b[k];
      }
    }

    for (std::size_t ii = 0; ii < N; ++ii) {
      const std::size_t i = N - 1 - ii;
      double sum = b[i];
      for (std::size_t j = i + 1; j < N; ++j)
        sum -= A[i][j] * x[j];
      const double aii = A[i][i];
      if (!(std::abs(aii) > eps))
        return false;
      x[i] = sum / aii;
      if (!std::isfinite(x[i]))
        return false;
    }

    return true;
  }

  template <std::size_t N> static bool solveLinearSystemLapackLU(const std::array<std::array<double, N>, N> &A, const std::array<double, N> &b, std::array<double, N> &x) {
    // Reuse the project's LAPACK wrapper (basic_linear_systems.hpp is included via basic.hpp).
    // It may throw on singular/illegal inputs; keep it contained (OpenMP-safe) by returning false.
    try {
      lapack_lu lu(A, x, b); // solves A.x = b at construction
      for (double v : x)
        if (!std::isfinite(v))
          return false;
      return true;
    } catch (...) {
      return false;
    }
  }

  template <std::size_t K> static std::array<double, K> pseBasis(const std::array<double, 3> &r) {
    const double x = r[0];
    const double y = r[1];
    const double z = r[2];
    if constexpr (K == 3) {
      return {x, y, z};
    } else {
      return {x, y, z, x * x, y * y, z * z, x * y, x * z, y * z};
    }
  }

  template <std::size_t K> bool computePSECorrectionCoeffs(const VortexParticle &pi, std::size_t i, const std::vector<VortexParticle> &particles, std::array<double, K> &coeff_out, double second_moment_target, double cutoff_ratio, double regularization_rel, double max_dimless_coeff) const {
    static constexpr double volume_eps = 1e-12;
    static constexpr double sigma_eps = 1e-12;
    static const double pi32 = std::pow(M_PI, 1.5);

    std::array<double, K> base{};
    std::array<std::array<double, K>, K> M{};
    std::size_t neighbor_count = 0;

    for (std::size_t j = 0; j < particles.size(); ++j) {
      if (j == i)
        continue;
      const auto &pj = particles[j];
      if (!(pj.volume > volume_eps))
        continue;

      const double sigma = 0.5 * (pi.sigma + pj.sigma);
      if (!(sigma > sigma_eps))
        continue;
      const double sigma2 = sigma * sigma;
      const double cutoff2 = (cutoff_ratio * cutoff_ratio) * sigma2;

      const auto r_vec = pj.x - pi.x; // r_ij (from i to j)
      const double r2 = Dot(r_vec, r_vec);
      if (r2 > cutoff2)
        continue;

      const double sigma5 = sigma2 * sigma2 * sigma;
      const double eta_prefactor = 2.0 / (pi32 * sigma5);
      const double eta = eta_prefactor * std::exp(-r2 / sigma2);
      const double w0 = pj.volume * eta;
      if (!(w0 > 0.0))
        continue;

      const auto phi = pseBasis<K>(r_vec);
      for (std::size_t a = 0; a < K; ++a)
        base[a] += w0 * phi[a];

      for (std::size_t a = 0; a < K; ++a) {
        for (std::size_t b = a; b < K; ++b)
          M[a][b] += w0 * phi[a] * phi[b];
      }
      ++neighbor_count;
    }

    if (neighbor_count < K)
      return false;

    for (std::size_t a = 0; a < K; ++a)
      for (std::size_t b = 0; b < a; ++b)
        M[a][b] = M[b][a];

    double max_diag = 0.0;
    for (std::size_t a = 0; a < K; ++a)
      max_diag = std::max(max_diag, std::abs(M[a][a]));
    const double lambda = (regularization_rel > 0.0) ? (regularization_rel * (max_diag > 0.0 ? max_diag : 1.0)) : 0.0;
    if (lambda > 0.0) {
      for (std::size_t a = 0; a < K; ++a)
        M[a][a] += lambda;
    }

    std::array<double, K> rhs{};
    for (std::size_t a = 0; a < K; ++a) {
      double target = 0.0;
      if constexpr (K == 9) {
        if (a == 3 || a == 4 || a == 5)
          target = second_moment_target;
      }
      rhs[a] = target - base[a]; // because p(r)=1 + Σ c_n φ_n
    }

    std::array<double, K> c{};
    if (!solveLinearSystemLapackLU<K>(M, rhs, c)) {
      if (!solveLinearSystemGaussian<K>(M, rhs, c))
        return false;
    }

    const double h = std::max(std::abs(pi.sigma), sigma_eps);
    double max_dimless = 0.0;
    for (std::size_t a = 0; a < K; ++a) {
      const double scale = (a < 3) ? h : (h * h);
      max_dimless = std::max(max_dimless, std::abs(c[a]) * scale);
    }
    if (!std::isfinite(max_dimless) || !(max_dimless <= max_dimless_coeff))
      return false;

    coeff_out = c;
    return true;
  }

  struct WallFluxStats {
    std::size_t faces_total = 0;
    std::size_t faces_absorbed = 0;
    std::size_t faces_shed = 0;
    std::size_t added = 0;
    // Sum of applied wall vorticity flux (alpha = omega*volume) per call.
    // `absorb` distributes into existing particles; `shed` creates new particles.
    Tddd sum_alpha_flux_absorbed = {0., 0., 0.};
    Tddd sum_alpha_flux_shed = {0., 0., 0.};
    double min_sigma_face = std::numeric_limits<double>::infinity();
    double max_sigma_face = 0.0;
    double min_sigma = std::numeric_limits<double>::infinity();
    double max_sigma = 0.0;
    double min_sigma_search = std::numeric_limits<double>::infinity();
    double max_sigma_search = 0.0;
  };

	  WallFluxStats injectWallVorticityFluxPSE_core(Network *boundaryNet, double dt, std::size_t min_absorb_receivers, double min_absorb_total_weight, double sigma_factor, bool allow_shed) {
	    WallFluxStats stats;
	    if (!boundaryNet)
	      return stats;
	    boundaryNet->setGeometricPropertiesForce();

    static constexpr double area_eps = 1e-12;
    static constexpr double alpha_eps = 1e-10;
    static constexpr double sigma_eps = 1e-12;
    static constexpr double dt_eps = 1e-16;
    static constexpr double length_eps = 1e-14;

    if (!(sigma_factor > 0.0) || !std::isfinite(sigma_factor))
      sigma_factor = 1.5;

    for (const auto &face : boundaryNet->getBoundaryFaces()) {
      if (!face->Neumann)
        continue;

      ++stats.faces_total;
      auto [p0, p1, p2] = face->getPoints();
      const double area = face->area;
      if (!(area > area_eps))
        continue;

      Network *contactBody = face->penetratedBody;
      if (!contactBody) {
        auto findBody = [&](networkPoint *p) -> Network * {
          if (p->penetratedBody)
            return p->penetratedBody;
          auto cf = p->getNearestContactFace(face);
          return cf ? cf->getNetwork() : nullptr;
        };
        contactBody = findBody(p0);
        if (!contactBody)
          contactBody = findBody(p1);
        if (!contactBody)
          contactBody = findBody(p2);
      }
      if (contactBody) {
        bool shed_enabled = true;
        contactBody->inputJSON.find("shed_vortices", [&](const auto &val) {
          if (!val.empty() && (val[0] == "false" || val[0] == "False" || val[0] == "FALSE"))
            shed_enabled = false;
        });
        if (!shed_enabled)
          continue;
      }
      if (!contactBody)
        continue;

      const Tddd normal = face->normal;
      const Tddd S_vec = normal * area;
      const Tddd center = (p0->X + p1->X + p2->X) / 3.0;

      auto get_u_diff = [&](networkPoint *p) {
        const Tddd u_potential_BEM = p->u_potential_BEM;
        const Tddd u_omega_VPM = computeVelocity(p->X);
        const Tddd u_body = contactBody->velocityRigidBody(p->X);
        const Tddd u_total = u_potential_BEM + u_omega_VPM;
        return u_total - u_body;
      };

	      const Tddd u_diff_avg = (get_u_diff(p0) + get_u_diff(p1) + get_u_diff(p2)) / 3.0;
	      // Convert a vorticity *flux rate* into an alpha increment for this time-step.
	      // Without the dt factor, the total injected circulation scales with the number of steps
	      // (i.e., does not converge as dt -> 0).
	      const Tddd alpha_flux = Cross(u_diff_avg, S_vec) * dt;
	      if (!(Norm(alpha_flux) > alpha_eps))
	        continue;

      const double l01 = Norm(p0->X - p1->X);
      const double l12 = Norm(p1->X - p2->X);
      const double l20 = Norm(p2->X - p0->X);
      const double mean_edge_length = (l01 + l12 + l20) / 3.0;
      if (!(mean_edge_length > length_eps) || !std::isfinite(mean_edge_length))
        continue;

      // sigma_face: 基本的な渦粒子間の距離の推定値．正三角形要素の中心間距離に等しい
      const double sigma_face = std::sqrt(3.0) * mean_edge_length / 3.0;
      if (!(sigma_face > sigma_eps) || !std::isfinite(sigma_face))
        continue;
      stats.min_sigma_face = std::min(stats.min_sigma_face, sigma_face);
      stats.max_sigma_face = std::max(stats.max_sigma_face, sigma_face);

      // sigma: 固定されたコアサイズ（新しい渦粒子でもこの値）およそ[1~2]*渦粒子間
      const double sigma = sigma_factor * sigma_face;
      if (!(sigma > sigma_eps) || !std::isfinite(sigma))
        continue;
      stats.min_sigma = std::min(stats.min_sigma, sigma);
      stats.max_sigma = std::max(stats.max_sigma, sigma);

      // sigma_heat: 物理的な拡散長
      double sigma_heat = 0.0;
      if (nu > 0.0 && dt > dt_eps) {
        const double sigma_heat2 = 4.0 * nu * dt;
        if (sigma_heat2 > sigma_eps * sigma_eps)
          sigma_heat = std::sqrt(sigma_heat2);
      }

      // Used only during injection:
      // - sigma_search: range for deciding absorption (distribute into existing near-wall particles).
      //   (We use `min_absorb_total_weight` to robustly decide absorb vs shed, even if the neighborhood is sparse.)
      // - sigma_search also defines the injection distance when shedding new particles.
      //   (No separate sigma_inject variable to avoid confusion.)
      const double sigma_search = std::max(sigma, sigma_heat);
      if (!(sigma_search > sigma_eps) || !std::isfinite(sigma_search))
        continue;
      stats.min_sigma_search = std::min(stats.min_sigma_search, sigma_search);
      stats.max_sigma_search = std::max(stats.max_sigma_search, sigma_search);

      // Query / shed point: place at a fixed wall-normal distance equal to the core size `sigma`,
      // so the injection location is determined by `sigma_factor` (via sigma = sigma_factor*sigma_face).
      // Search uses `sigma_search = max(sigma, sigma_heat)` to robustly decide absorb vs shed.
      const Tddd x_query = center - normal * sigma;

      // Absorption search cutoff: use sigma_search as the search radius.
      const double cutoff2 = sigma_search * sigma_search;

      double sum_w = 0.0;
      std::vector<std::pair<std::size_t, double>> receivers;
      receivers.reserve(16);
      for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto &pi = particles[i];
        const auto r = pi.x - x_query;
        const double r2 = Dot(r, r);
        if (r2 > cutoff2)
          continue;
        const double w = std::max(0.0, pi.volume) * std::exp(-r2 / (sigma_search * sigma_search));
        if (!(w > 0.0))
          continue;
        receivers.emplace_back(i, w);
        sum_w += w;
      }

      const bool can_absorb = (receivers.size() >= min_absorb_receivers) && (sum_w > min_absorb_total_weight);
      if (can_absorb) {
        for (const auto &[i, w] : receivers) {
          const double s = w / sum_w;
          particles[i].alpha = particles[i].alpha + alpha_flux * s;
        }
        stats.sum_alpha_flux_absorbed = stats.sum_alpha_flux_absorbed + alpha_flux;
        ++stats.faces_absorbed;
        continue;
      }

      if (!allow_shed)
        continue;

      const Tddd x_shed = x_query;
      const double volume = area * sigma;
      this->addParticle(x_shed, alpha_flux, sigma, volume);
      stats.sum_alpha_flux_shed = stats.sum_alpha_flux_shed + alpha_flux;
      ++stats.faces_shed;
      ++stats.added;
    }

    return stats;
  }

public:
  VortexMethod(double nu_in = _WATER_NU_10deg_) : nu(nu_in) {}
  ~VortexMethod() {}

  void setStretchingScheme(StretchingScheme scheme) { stretching_scheme = scheme; }
  StretchingScheme getStretchingScheme() const { return stretching_scheme; }

  void setPSECorrectionMode(PSECorrectionMode mode) {
    pse_correction_mode = mode;
    // Auto-select defaults from the current kernel + chosen correction mode.
    if (pse_correction_mode == PSECorrectionMode::Curvature)
      pse_correction_second_moment_target = 2.0 * pse_kernel_second_moment;
  }
  PSECorrectionMode getPSECorrectionMode() const { return pse_correction_mode; }
  void setPSECorrectionSecondMomentTarget(double second_moment_target) {
    if (!(second_moment_target > 0.0) || !std::isfinite(second_moment_target))
      return;
    pse_correction_second_moment_target = second_moment_target;
  }
  double getPSECorrectionSecondMomentTarget() const { return pse_correction_second_moment_target; }

  void addParticle(const std::array<double, 3> &x, const std::array<double, 3> &alpha, double sigma, double volume) { particles.emplace_back(x, alpha, sigma, volume); }

  const std::vector<VortexParticle> &getParticles() const { return particles; }
  std::vector<VortexParticle> &getParticles() { return particles; }

  struct SuggestedTimeStep {
    double dt = std::numeric_limits<double>::infinity();
    double dt_strain = std::numeric_limits<double>::infinity();    // accuracy constraint (convection)
    double dt_diffusion = std::numeric_limits<double>::infinity(); // stability constraint (PSE, explicit)
    double dt_move = std::numeric_limits<double>::infinity();      // overlap/accuracy heuristic
    double max_strain_rate = 0.0;
    double max_u = 0.0;
    double min_sigma = std::numeric_limits<double>::infinity();
    std::size_t particle_count = 0;
  };

  SuggestedTimeStep suggestTimeStep(double max_dt, double C_strain = 0.2, double C_diffusion = 0.25, double C_move = 0.5) const {
    SuggestedTimeStep out;
    out.particle_count = particles.size();
    if (particles.empty()) {
      out.dt = max_dt;
      return out;
    }

    constexpr double eps = 1e-16;
    for (const auto &p : particles) {
      out.max_u = std::max(out.max_u, Norm(p.u_total));
      out.min_sigma = std::min(out.min_sigma, p.sigma);
      const double a = Norm(p.alpha);
      const double da = Norm(p.d_alpha_dt);
      if (a > eps) {
        const double strain = da / a; // ||dα/dt|| / ||α|| ≈ O(||∇u||)
        out.max_strain_rate = std::max(out.max_strain_rate, strain);
      }
    }

    if (out.max_strain_rate > eps)
      out.dt_strain = C_strain / out.max_strain_rate;

    if (nu > 0.0 && std::isfinite(out.min_sigma) && out.min_sigma > eps)
      out.dt_diffusion = C_diffusion * (out.min_sigma * out.min_sigma) / nu;

    if (out.max_u > eps && std::isfinite(out.min_sigma) && out.min_sigma > eps)
      out.dt_move = C_move * out.min_sigma / out.max_u;

    out.dt = std::min({max_dt, out.dt_strain, out.dt_diffusion, out.dt_move});
    return out;
  }

  /**
   * @brief 条件に合致する粒子を削除する
   * @param predicate trueを返すとその粒子は削除される
   */
  void removeParticles(const std::function<bool(const VortexParticle &)> &predicate) {
    auto it = std::remove_if(particles.begin(), particles.end(), predicate);
    particles.erase(it, particles.end());
  }

  /**
   * @brief 全粒子に対して処理を適用する（座標の修正などに使用）
   */
  void processParticles(const std::function<void(VortexParticle &)> &processor) {
    for (auto &p : particles)
      processor(p);
  }

  /**
   * @brief Biot-Savartの法則 (Rosenhead Moore カーネル)
   * * u(x) = (1/4pi) * sum_j [ (alpha_j x r_ij) / (r_ij^2 + sigma_j^2)^(3/2) ]
   */

  std::array<double, 3> computeVelocity(const std::array<double, 3> &target_x) const {
    std::array<double, 3> u = {0.0, 0.0, 0.0};
    for (const auto &p : particles) {
      // Rosenheadカーネル (代数カーネル)
      // 分母: (r^2 + sigma^2)^(3/2)
      std::array<double, 3> r_vec = target_x - p.x;
      u += CrossDouble(p.alpha, r_vec) / (4.0 * M_PI * std::pow(Dot(r_vec, r_vec) + p.sigma * p.sigma, 1.5));
    }
    return u;
  }

  /**
   * @brief ストレッチング項の計算
   * Standard:  dα_i/dt = (α_i·∇)u = (∇u) α_i
   * Transpose: dα_i/dt = (∇u)^T α_i   (Cottet & Koumoutsakos Eq. 3.1.7)
   */
  void computeStretching() {
// ストレッチング計算はO(N^2)なので並列化推奨
	#pragma omp parallel for
	    for (size_t i = 0; i < particles.size(); ++i) {
	      std::array<double, 3> stretching = {0.0, 0.0, 0.0};

	      for (size_t j = 0; j < particles.size(); ++j) {
	        if (i == j)
	          continue;

	        const auto &p_j = particles[j];
	        const auto &p_i = particles[i];

        std::array<double, 3> r_vec = p_i.x - p_j.x; // r_ij
        double r2 = Dot(r_vec, r_vec);
        // 正則化パラメータ sigma はソース側(j)のものを使うのが一般的だが、
        // 対称性を保つために (sigma_i^2 + sigma_j^2)/2 などを使う手法もある。
        // ここではシンプルに p_j.sigma を使用。
        double R2 = r2 + p_j.sigma * p_j.sigma;

        double R = std::sqrt(R2);
        double R3 = R * R2;   // R^3 = (r^2+σ^2)^{3/2}
        double R5 = R3 * R2;  // R^5 = (r^2+σ^2)^{5/2}

        double factor1 = 1.0 / (4.0 * M_PI * R3);
        double factor2 = 3.0 / (4.0 * M_PI * R5);

        if (stretching_scheme == StretchingScheme::Transpose) {
          // (∇u)^T α = ∇(u·α) for constant α.
          // For Rosenhead kernel u = (α_j × r)/ (4π R^3):
          // (∇u)^T α_i = (1/4π)[ (α_i×α_j)/R^3 - 3 ((α_i×α_j)·r) r / R^5 ].
          const auto cross_ai_aj = Cross(p_i.alpha, p_j.alpha);
          const auto term1 = cross_ai_aj * factor1;
          const double dot_cross_r = Dot(cross_ai_aj, r_vec);
          const auto term2 = r_vec * (factor2 * dot_cross_r);
          stretching = stretching + (term1 - term2);
        } else {
          // Standard: (α_i·∇)u (directional derivative)
          // dα_i/dt = (1/4π) * sum_j [ (α_j×α_i)/R^3 - 3 (α_i·r)(α_j×r)/R^5 ].
          const auto term1 = Cross(p_j.alpha, p_i.alpha) * factor1;
          const double dot_ai_r = Dot(p_i.alpha, r_vec);
          const auto cross_aj_r = Cross(p_j.alpha, r_vec);
          const auto term2 = cross_aj_r * (factor2 * dot_ai_r);
          stretching = stretching + (term1 - term2);
        }
      }
      particles[i].d_alpha_dt = stretching;
    }
  }

  /**
   * @brief 粘性拡散項の計算 (PSE: Potential Singularity Expansion)
   * * d(alpha_i)/dt = nu * nabla^2 (alpha_i)
   * * ここでは、他の粒子との相互作用による拡散項を計算する。
   * @param nu 粘性係数
   */

  void computeDiffusionPSE() {
    if (nu <= 0.0 || particles.size() < 2)
      return;

    static constexpr double volume_eps = 1e-12;
    static constexpr double cutoff_ratio = 4.0; // 4σ 以遠は無視
    static const double pi32 = std::pow(M_PI, 1.5);
    // With the current Gaussian η prefactor in this file, the continuous 2nd moment is:
    //   ∫ η(r;σ) r_i^2 dV = 1
    // so the Laplacian prefactor becomes 2 to cancel Taylor's 1/2. Curvature correction enforces that
    // moment target directly (≈2), so we only need this explicit factor for None/Gradient modes.
    static constexpr double laplacian_scale_uncorrected = 2.0;

    const std::size_t n_perm = particles.size();

    if (pse_correction_mode == PSECorrectionMode::None) {
#pragma omp parallel for
      for (std::size_t i = 0; i < n_perm; ++i) {
        auto &pi = particles[i];
        if (pi.volume <= volume_eps)
          continue;
        const auto omega_i = pi.alpha * (1.0 / pi.volume);
        std::array<double, 3> diffusion_i = {0.0, 0.0, 0.0};

        for (std::size_t j = 0; j < n_perm; ++j) {
          if (i == j)
            continue;

          const auto &pj = particles[j];
          if (pj.volume <= volume_eps)
            continue;
          const auto omega_j = pj.alpha * (1.0 / pj.volume);

          const double sigma = 0.5 * (pi.sigma + pj.sigma);
          if (sigma <= 1e-12)
            continue;
          const double sigma2 = sigma * sigma;
          const double cutoff2 = cutoff_ratio * cutoff_ratio * sigma2;

          const auto r_vec = pi.x - pj.x;
          const double r2 = Dot(r_vec, r_vec);
          if (r2 > cutoff2)
            continue;

          const double sigma5 = sigma2 * sigma2 * sigma;
          const double eta_prefactor = 2.0 / (pi32 * sigma5);
          const double eta = eta_prefactor * std::exp(-r2 / sigma2);

          const auto omega_diff = omega_j - omega_i;
          const double weight = laplacian_scale_uncorrected * nu * pi.volume * pj.volume * eta;
          diffusion_i = diffusion_i + omega_diff * weight;
        }
        pi.d_alpha_dt = pi.d_alpha_dt + diffusion_i;
      }
      return;
    }

    const bool use_curvature = (pse_correction_mode == PSECorrectionMode::Curvature);
    const double laplacian_scale = use_curvature ? 1.0 : laplacian_scale_uncorrected;

    // Precompute per-particle correction coefficients (only for permanent particles).
    // To preserve the antisymmetry of the exchange flux (and thus conservation),
    // the pairwise correction factor is symmetrized as: f_ij = 0.5*(p_i(r_ij) + p_j(r_ji)).
    std::vector<std::array<double, 9>> corr(n_perm);
#pragma omp parallel for
    for (std::size_t i = 0; i < n_perm; ++i) {
      corr[i].fill(0.0);
      const auto &pi = particles[i];
      if (pi.volume <= volume_eps)
        continue;

      if (pse_correction_mode == PSECorrectionMode::Gradient) {
        std::array<double, 3> c3{};
        if (computePSECorrectionCoeffs<3>(pi, i, particles, c3, 0.0, cutoff_ratio, pse_correction_regularization_rel, pse_correction_max_dimless_coeff)) {
          corr[i][0] = c3[0];
          corr[i][1] = c3[1];
          corr[i][2] = c3[2];
        }
      } else if (pse_correction_mode == PSECorrectionMode::Curvature) {
        std::array<double, 9> c9{};
        if (computePSECorrectionCoeffs<9>(pi, i, particles, c9, pse_correction_second_moment_target, cutoff_ratio, pse_correction_regularization_rel, pse_correction_max_dimless_coeff)) {
          corr[i] = c9;
        }
      }
    }

#pragma omp parallel for
    for (std::size_t i = 0; i < n_perm; ++i) {
      auto &pi = particles[i];
      if (pi.volume <= volume_eps)
        continue;
      const auto omega_i = pi.alpha * (1.0 / pi.volume);
      std::array<double, 3> diffusion_i = {0.0, 0.0, 0.0};

      for (std::size_t j = 0; j < n_perm; ++j) {
        if (i == j)
          continue;
        const auto &pj = particles[j];
        if (pj.volume <= volume_eps)
          continue;
        const auto omega_j = pj.alpha * (1.0 / pj.volume);

        const double sigma = 0.5 * (pi.sigma + pj.sigma);
        if (sigma <= 1e-12)
          continue;
        const double sigma2 = sigma * sigma;
        const double cutoff2 = cutoff_ratio * cutoff_ratio * sigma2;

        const auto r_ij = pj.x - pi.x; // r_ij (from i to j)
        const double r2 = Dot(r_ij, r_ij);
        if (r2 > cutoff2)
          continue;

        const double sigma5 = sigma2 * sigma2 * sigma;
        const double eta_prefactor = 2.0 / (pi32 * sigma5);
        const double eta = eta_prefactor * std::exp(-r2 / sigma2);

        const double x = r_ij[0];
        const double y = r_ij[1];
        const double z = r_ij[2];
        const double x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const double xy = x * y;
        const double xz = x * z;
        const double yz = y * z;

        const auto &ci = corr[i];
        double p_i = 1.0 + ci[0] * x + ci[1] * y + ci[2] * z;
        if (use_curvature)
          p_i += ci[3] * x2 + ci[4] * y2 + ci[5] * z2 + ci[6] * xy + ci[7] * xz + ci[8] * yz;

        const auto &cj = corr[j];
        // r_ji = -r_ij: linear terms flip sign, quadratic terms don't.
        double p_j = 1.0 + cj[0] * (-x) + cj[1] * (-y) + cj[2] * (-z);
        if (use_curvature)
          p_j += cj[3] * x2 + cj[4] * y2 + cj[5] * z2 + cj[6] * xy + cj[7] * xz + cj[8] * yz;

        const double f_sym = 0.5 * (p_i + p_j);
        const double weight = laplacian_scale * nu * pi.volume * pj.volume * eta * f_sym;
        diffusion_i = diffusion_i + (omega_j - omega_i) * weight;
      }

      pi.d_alpha_dt = pi.d_alpha_dt + diffusion_i;
    }
  }

  /**
   * @brief PSE "integral technique" style wall vorticity injection:
   * First try to distribute the wall vorticity flux into existing near-wall particles.
   * Only if there are no receivers (or insufficient overlap), create ("shed") new particles.
   *
   * This targets the coupling point between diffusion (PSE) and boundary conditions (no-slip via slip velocity).
   *
   * Notes:
   * - This updates `particles[i].alpha` directly (i.e., an accumulated increment for the current step),
   *   which is consistent with the existing explicit boundary injection style in this codebase.
   */
  // Convenience overload with explicit defaults (kept minimal to avoid silently ignoring JSON-configured values in callers).
	  void injectWallVorticityFluxPSE(Network *boundaryNet, double dt) { this->injectWallVorticityFluxPSE(boundaryNet, dt, 1, 1e-14, 1.5, true); }

	  void injectWallVorticityFluxPSE(Network *boundaryNet, double dt, std::size_t min_absorb_receivers, double min_absorb_total_weight, double sigma_factor, bool allow_shed) {
	    const auto stats = injectWallVorticityFluxPSE_core(boundaryNet, dt, min_absorb_receivers, min_absorb_total_weight, sigma_factor, allow_shed);
	    if (stats.faces_absorbed > 0 || stats.faces_shed > 0) {
        Tddd sum_alpha = {0., 0., 0.};
        for (const auto &p : particles)
          sum_alpha = sum_alpha + p.alpha;
        const Tddd sum_alpha_flux = stats.sum_alpha_flux_absorbed + stats.sum_alpha_flux_shed;

		      std::cout << "[VPM:wall_bc] mode=absorbPSE faces=" << stats.faces_total << " absorbed=" << stats.faces_absorbed << " shed=" << stats.faces_shed << " added=" << stats.added;
		      if (std::isfinite(stats.min_sigma_face))
		        std::cout << " sigma_face[min,max]=" << stats.min_sigma_face << "," << stats.max_sigma_face;
	      if (std::isfinite(stats.min_sigma))
	        std::cout << " sigma[min,max]=" << stats.min_sigma << "," << stats.max_sigma;
	      if (std::isfinite(stats.min_sigma_search))
	        std::cout << " sigma_search[min,max]=" << stats.min_sigma_search << "," << stats.max_sigma_search;
	      std::cout << " particles=" << particles.size();
        std::cout << " sum_alpha_z=" << sum_alpha[2];
        std::cout << " sum_alpha_flux_z=" << sum_alpha_flux[2];
        std::cout << " sum_alpha_flux_absorb_z=" << stats.sum_alpha_flux_absorbed[2];
        std::cout << " sum_alpha_flux_shed_z=" << stats.sum_alpha_flux_shed[2];
        std::cout << std::endl;
	    }
	  }

  void set_u_potential_BEM(std::function<std::array<double, 3>(const std::array<double, 3> &)> potentialField) {
    this->potentialField = potentialField;
    if (!this->potentialField)
      return;
    for (auto &p : particles) {
      p.u_potential_BEM = this->potentialField(p.x);
    }
  }

  // --- RK Integration Helpers ---

  void initializeRK(double dt, double t0, int steps = 4) {
    for (auto &p : particles) {
      p.RK_x.initialize(dt, t0, p.x, steps);
      p.RK_alpha.initialize(dt, t0, p.alpha, steps);
    }
  }

  void calcVelocityAndDerivatives() {
    if (potentialField) {
#pragma omp parallel for
      for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].u_potential_BEM = potentialField(particles[i].x);
      }
    }

#pragma omp parallel for
    for (size_t i = 0; i < particles.size(); ++i) {
      particles[i].u_omega_VPM = computeVelocity(particles[i].x);
      particles[i].u_total = particles[i].u_omega_VPM + particles[i].u_potential_BEM;
    }

    // NOTE: Diffusion is intentionally NOT included here.
    // We advance advection + stretching with RK, and apply PSE diffusion once per time-step
    // with an explicit (Euler) step (operator splitting) in `applyDiffusionEuler(...)`.
    computeStretching(); // overwrites d_alpha_dt
  }

  // Explicit diffusion step applied once per time-step (operator splitting).
  // This updates alpha as: alpha <- alpha + dt * PSE(alpha).
  // `d_alpha_dt` is restored to the pre-diffusion value so that dt heuristics based on
  // stretching (e.g., `suggestTimeStep`) are not polluted by the diffusion-rate values.
  void applyDiffusionEuler(double dt) {
    if (!(dt > 0.0))
      return;
    if (nu <= 0.0 || particles.size() < 2)
      return;

    std::vector<std::array<double, 3>> backup_dadt;
    backup_dadt.reserve(particles.size());
    for (const auto &p : particles)
      backup_dadt.push_back(p.d_alpha_dt);

    for (auto &p : particles)
      p.d_alpha_dt = {0.0, 0.0, 0.0};

    computeDiffusionPSE(); // adds diffusion to d_alpha_dt

    for (auto &p : particles)
      p.alpha = p.alpha + p.d_alpha_dt * dt;

    for (std::size_t i = 0; i < particles.size(); ++i)
      particles[i].d_alpha_dt = backup_dadt[i];
  }

  void calcVelocityAndStretching() { calcVelocityAndDerivatives(); }

  // Explicit Euler step for convection (x) and stretching (alpha).
  // Requires `potentialField` to be set (e.g. via `set_u_potential_BEM(...)`).
  void applyAdvectionStretchingEuler(double dt) {
    if (!(dt > 0.0))
      return;
    if (particles.empty())
      return;

    calcVelocityAndStretching(); // sets u_total and d_alpha_dt (stretching only)

    for (auto &p : particles) {
      p.x = p.x + p.u_total * dt;
      p.alpha = p.alpha + p.d_alpha_dt * dt;
    }
  }

  void pushRK() {
    for (auto &p : particles) {
      p.RK_x.push(p.u_total);
      p.RK_alpha.push(p.d_alpha_dt);
      p.x = p.RK_x.getX();
      p.alpha = p.RK_alpha.getX();
    }
  }
};
