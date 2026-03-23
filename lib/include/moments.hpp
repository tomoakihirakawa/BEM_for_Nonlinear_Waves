#pragma once

// Extracted Moments struct from lib_spatial_partitioning.hpp (nested Buckets<T,N>::Moments)
// This header is preparatory: after integration, remove the nested struct in lib_spatial_partitioning.hpp
// and include this header instead. Until that removal, DO NOT include this file in a TU that also
// includes the old nested definition to avoid ODR / duplicate symbol issues.

// Forward declaration for runtime configuration

#include <array>
#include <cctype>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
#if defined(_OPENMP)
struct OmpScheduleSpec {
  omp_sched_t kind;
  int chunk;
};

inline OmpScheduleSpec parse_omp_schedule(const char* env, omp_sched_t fallback_kind, int fallback_chunk) {
  if (env == nullptr || *env == '\0')
    return {fallback_kind, fallback_chunk};
  std::string s(env);
  for (char& ch : s)
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  const auto comma = s.find(',');
  const std::string kind_str = (comma == std::string::npos) ? s : s.substr(0, comma);
  int chunk = fallback_chunk;
  if (comma != std::string::npos) {
    const char* chunk_str = s.c_str() + comma + 1;
    char* endp = nullptr;
    const long v = std::strtol(chunk_str, &endp, 10);
    if (endp != chunk_str && v > 0)
      chunk = static_cast<int>(v);
  }

  if (kind_str == "static")
    return {omp_sched_static, chunk};
  if (kind_str == "dynamic" || kind_str == "dyn")
    return {omp_sched_dynamic, chunk};
  if (kind_str == "guided")
    return {omp_sched_guided, chunk};
  if (kind_str == "auto")
    return {omp_sched_auto, chunk};
  return {fallback_kind, fallback_chunk};
}
#endif
} // namespace

namespace {
// Real-field optimization (m<0 conjugate symmetry)
//
// If the boundary data is real (time-domain), the multipole/local coefficients satisfy
//   M_{n,-m} = conj(M_{n,m})   (for this codebase's Rnm / l2p basis)
// and similarly for the local expansion.
//
// In that case we can compute/update only m>=0 coefficients and fill m<0 by conjugation.
// This roughly halves the work in M2M/M2L/L2L for real-field runs.
//
// IMPORTANT:
// - For genuinely complex boundary data (frequency-domain complex amplitude), this symmetry does NOT hold.
// - Use `BEM_FMM_REALFIELD_M_CONJ` to override at runtime (0/1).
// - Default is controlled by `BEM_FMM_REALFIELD_M_CONJ_DEFAULT` (compile-time, 0/1).
#ifndef BEM_FMM_REALFIELD_M_CONJ_DEFAULT
#define BEM_FMM_REALFIELD_M_CONJ_DEFAULT 1
#endif

inline bool bem_fmm_use_realfield_m_conj() {
  static const bool enabled = [] {
    if (const char* env = std::getenv("BEM_FMM_REALFIELD_M_CONJ")) {
      return std::atoi(env) != 0;
    }
    return BEM_FMM_REALFIELD_M_CONJ_DEFAULT != 0;
  }();
  return enabled;
}
} // namespace

namespace {
// SimpleM2L term ordering
//
// In SimpleM2L, we accumulate many terms per (j,k) destination coefficient.
// The ordering can impact:
// - numerical error (small-to-large sum is slightly more stable)
// - memory locality (grouping by source bucket can reduce cache misses)
//
// Runtime control:
//   BEM_SIMPLEM2L_TERM_ORDER=abs   (default) sort by |AAAY| ascending (more stable)
//   BEM_SIMPLEM2L_TERM_ORDER=none  keep insertion order (often better locality)
//   BEM_SIMPLEM2L_TERM_ORDER=src   stable-sort by source pointer (group sources)
enum class SimpleM2LTermOrder {
  Abs,
  None,
  Src
};

inline SimpleM2LTermOrder bem_simplem2l_term_order() {
  static const SimpleM2LTermOrder order = [] {
    const char* env = std::getenv("BEM_SIMPLEM2L_TERM_ORDER");
    if (env == nullptr || *env == '\0')
      return SimpleM2LTermOrder::Abs;
    std::string s(env);
    for (char& ch : s)
      ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (s == "none")
      return SimpleM2LTermOrder::None;
    if (s == "src" || s == "source")
      return SimpleM2LTermOrder::Src;
    return SimpleM2LTermOrder::Abs;
  }();
  return order;
}
} // namespace

// Project specific dependencies (provide SphericalCoordinates, Fourier2DConvolution, constants, helpers)
#include "lib_Fourier.hpp"             // Fourier2DConvolution
#include "lib_multipole_expansion.hpp" // make_zero_MM, make_nm_set, computeSizeM2M, AAA_* tables, i_absk_A_FMM

// PlaneWaveM2L: include plane-wave expansion library (must be before class definition)
#if defined(USE_PlaneWaveM2L)
#include "lib_plane_wave_m2l.hpp" // Full plane wave with 6-direction decomposition and realfield optimization
#endif

// Metal M2L: GPU-accelerated M2L transformation (optional)
#if defined(USE_METAL_M2L)
#include "metal_m2l_wrapper.hpp"
#endif

// SIMD nearfield: NEON-optimized direct integration (ARM64 only)
#if defined(__aarch64__)
#include "setDirectIntegration_SIMD.hpp"
#include "setDirectIntegration_SIMD_double.hpp"
#endif

// Metal nearfield: GPU-accelerated nearfield direct integration (optional)
#if defined(USE_METAL_NEARFIELD)
#include "metal_nearfield_wrapper.hpp"
#endif

// Forward declaration of Buckets to allow pointer members without full definition
template <typename T, int N>
struct Buckets;

// Moments class extracted from Buckets<T,N>
// T: underlying object type used by Buckets (only appears indirectly via Buckets<T,N>* caches)
// N: truncation order of multipole/local expansions

//@ -------------------------------------------------------------------------- */
//@              モーメント:　ある観測点で，ソース点を近似するための係数配列             */
//@ -------------------------------------------------------------------------- */
template <typename T, int N>
struct Moments {
  /*
  N = 0 [ 0,  ,  ,  ,  ] n=0   total 1 = (N+1)^2             center : 0, in 1D 0
  N = 1 [-1, 0, 1,  ,  ] n=1   total 1 + 3 = 4 = (N+1)^2     center : 1, in 1D 2
  N = 2 [-2,-1, 0, 1, 2] n=2=N total 4 + 5 = 9 = (N+1)^2     center : 2, in 1D 6
        <-- 2*N+1=5 -->
  */
  std::array<double, 3> X;
  using type_MM = std::array<std::array<cmplx, 2>, (N + 1) * (N + 1)>;
  type_MM MM_ = make_zero_MM<N>();
  //! Mは，sum {q_i*\rho_i^n*Y_n^-m}
  //! Greengard & Rokhlin(1997)の式(3.19)
  constexpr static std::array<std::array<int, 2>, (N + 1) * (N + 1)> nm_set = make_nm_set<N>();

  Moments() : X{0.0, 0.0, 0.0} {}
  Moments(const std::array<double, 3>& XIN) : X(XIN) {}

  const cmplx zero = {0., 0.};

  std::size_t index(int n, int m) const { return n * (n + 1) + m; }

  /* -------------------------------------------------------------------------- */

  struct SparseOp {
    int j;
    int k;
    int n;
    int m;
    std::uint16_t dst;
    std::uint16_t src;
  };

  static_assert((N + 1) * (N + 1) <= 65536, "N too large for uint16_t coefficient indices");

  static const std::vector<SparseOp>& m2m_ops() {
    static const std::vector<SparseOp> ops = []() {
      std::vector<SparseOp> v;
      v.reserve(static_cast<std::size_t>(computeSizeM2M(static_cast<std::size_t>(N))));
      for (int j = N; j >= 0; --j) {
        for (int k = -j; k <= j; ++k) {
          auto& AAA_M2M_FMM_j_k = AAA_M2M_FMM[j][k + N_AAA_M2M_FMM];
          for (int n = j; n >= 0; --n) {
            for (int m = -n; m <= n; ++m) {
              const auto AAA = AAA_M2M_FMM_j_k[n][m + N_AAA_M2M_FMM];
              if (AAA.real() == 0.0 && AAA.imag() == 0.0)
                continue;
              const std::size_t dst = static_cast<std::size_t>(j * (j + 1) + k);
              const std::size_t src = static_cast<std::size_t>((j - n) * ((j - n) + 1) + (k - m));
              v.push_back({j, k, n, m, static_cast<std::uint16_t>(dst), static_cast<std::uint16_t>(src)});
            }
          }
        }
      }
      return v;
    }();
    return ops;
  }

  static const std::vector<SparseOp>& m2l_ops() {
    static const std::vector<SparseOp> ops = []() {
      std::vector<SparseOp> v;
      v.reserve(static_cast<std::size_t>((N + 1) * (N + 1) * (N + 1) * (N + 1)));
      for (int j = N; j >= 0; --j) {
        for (int k = -j; k <= j; ++k) {
          auto& AAA_M2L_FMM_j_k = AAA_M2L_FMM[j][k + N_AAA_M2L_FMM];
          for (int n = N; n >= 0; --n) {
            for (int m = -n; m <= n; ++m) {
              const auto AAA = AAA_M2L_FMM_j_k[n][m + N_AAA_M2L_FMM];
              if (AAA.real() == 0.0 && AAA.imag() == 0.0)
                continue;
              const std::size_t dst = static_cast<std::size_t>(j * (j + 1) + k);
              const std::size_t src = static_cast<std::size_t>(n * (n + 1) + m);
              v.push_back({j, k, n, m, static_cast<std::uint16_t>(dst), static_cast<std::uint16_t>(src)});
            }
          }
        }
      }
      return v;
    }();
    return ops;
  }

  static const std::vector<SparseOp>& l2l_ops() {
    static const std::vector<SparseOp> ops = []() {
      std::vector<SparseOp> v;
      v.reserve(static_cast<std::size_t>((N + 1) * (N + 1) * (N + 1) * (N + 1)));
      for (int j = N; j >= 0; --j) {
        for (int k = -j; k <= j; ++k) {
          auto& AAA_L2L_FMM_j_k = AAA_L2L_FMM[j][k + N_AAA_L2L_FMM];
          for (int n = N; n >= j; --n) {
            for (int m = -n; m <= n; ++m) {
              const auto AAA = AAA_L2L_FMM_j_k[n][m + N_AAA_L2L_FMM];
              if (AAA.real() == 0.0 && AAA.imag() == 0.0)
                continue;
              const std::size_t dst = static_cast<std::size_t>(j * (j + 1) + k);
              const std::size_t src = static_cast<std::size_t>(n * (n + 1) + m);
              v.push_back({j, k, n, m, static_cast<std::uint16_t>(dst), static_cast<std::uint16_t>(src)});
            }
          }
        }
      }
      return v;
    }();
    return ops;
  }

  void initialize() { this->MM_ = make_zero_MM<N>(); }

  void initialize(const std::array<double, 3>& XIN) {
    this->X = XIN;
    this->MM_ = make_zero_MM<N>();
  }

  // Real-field only optimization:
  // For time-domain runs the boundary values/densities are real, and the spherical-harmonic coefficients
  // in this codebase satisfy:
  //   MM_(n,-m) = conj(MM_(n,m))   (m>=1)
  // This lets us compute only m>=0 coefficients and reconstruct m<0 cheaply.
  //
  // NOTE: Do NOT enable this for genuinely complex boundary data (frequency-domain complex amplitude),
  // because the symmetry generally does not hold.
  void enforce_conjugate_m_symmetry_if_enabled() {
    if (!bem_fmm_use_realfield_m_conj())
      return;
    for (int n = 0; n <= N; ++n) {
      for (int m = 1; m <= n; ++m) {
        auto& pos = this->MM_[this->index(n, m)];
        auto& neg = this->MM_[this->index(n, -m)];
        neg[0] = std::conj(pos[0]);
        neg[1] = std::conj(pos[1]);
      }
    }
  }

  /* -------------------------------------------------------------------------- */

  using TupleType = std::tuple<std::array<cmplx, 2>, std::array<cmplx, 2>*>;
  using ArrayType = std::array<TupleType, (N + 1) * (N + 1)>;
  std::vector<std::tuple<Tdd*, ArrayType>> vector_weightedSourceDensities_RRn_MM_;

  //^OK
  template <typename TYPE>
  void increment_M(const std::vector<TYPE>& elements) {
    this->vector_weightedSourceDensities_RRn_MM_.clear();
    // 要素数 × 内部ソース数 でリザーブ
    std::size_t total_sources = 0;
    for (const auto& elem : elements) {
      if (!elem->p2m_sources.empty())
        total_sources += elem->p2m_sources.size();
      else
        total_sources += 1; // 後方互換: p2m_sources 未使用の場合
    }
    this->vector_weightedSourceDensities_RRn_MM_.reserve(total_sources);
    for (const auto& elem : elements) {
      if (!elem->p2m_sources.empty()) {
        // 多点P2M: 各内部ソースについて球面調和関数を計算
        for (auto& isrc : elem->p2m_sources) {
          auto [w_phi, w_phin] = isrc.weighted_source_densities;
          SphericalCoordinates kernel(isrc.X - this->X);
          ArrayType RRn_MM_;
          std::array<cmplx, 2> RRn;
          for (std::size_t ind = nm_set.size(); ind-- > 0;) {
            auto& [n, m] = nm_set[ind];
            RRn = kernel.p2mFunction(n, m, isrc.normal);
            std::get<0>(MM_[ind]) += w_phin * std::get<0>(RRn);
            std::get<1>(MM_[ind]) += w_phi * std::get<1>(RRn);
            RRn_MM_[ind] = {RRn, &MM_[ind]};
          }
          this->vector_weightedSourceDensities_RRn_MM_.emplace_back(
              &isrc.weighted_source_densities, RRn_MM_);
        }
      } else {
        // 後方互換: 旧式の単一ソース (weighted_source_densities 直接参照)
        auto [w_phi, w_phin] = elem->weighted_source_densities;
        SphericalCoordinates kernel(elem->X - this->X);
        ArrayType RRn_MM_;
        std::array<cmplx, 2> RRn;
        for (std::size_t ind = nm_set.size(); ind-- > 0;) {
          auto& [n, m] = nm_set[ind];
          RRn = kernel.p2mFunction(n, m, elem->normal);
          std::get<0>(MM_[ind]) += w_phin * std::get<0>(RRn);
          std::get<1>(MM_[ind]) += w_phi * std::get<1>(RRn);
          RRn_MM_[ind] = {RRn, &MM_[ind]};
        }
        this->vector_weightedSourceDensities_RRn_MM_.emplace_back(
            &elem->weighted_source_densities, RRn_MM_);
      }
    }
    this->enforce_conjugate_m_symmetry_if_enabled();
  }

  //^OK
  void increment_M_reuse() {
    std::array<double, 2> w_phi_w_phin;
    for (const auto& [weighted_source_densities, RRn_MM_] : this->vector_weightedSourceDensities_RRn_MM_) {
      w_phi_w_phin = *weighted_source_densities;
      for (auto& [RRn, mm_] : RRn_MM_) {
        std::get<0>(*mm_) += std::get<1>(w_phi_w_phin) * std::get<0>(RRn);
        std::get<1>(*mm_) += std::get<0>(w_phi_w_phin) * std::get<1>(RRn);
      }
    }
    this->enforce_conjugate_m_symmetry_if_enabled();
  }

  //^OK
  std::array<double, 2> L2P(const std::array<double, 3>& a) const {
    SphericalCoordinates P(a - this->X);
    std::array<cmplx, 2> ret = {0, 0};
    cmplx Ynm_rhon1;
    for (std::size_t ind = nm_set.size(); ind-- > 0;) {
      auto& [n, m] = nm_set[ind];
      Ynm_rhon1 = P.l2pFunction(n, m);
      std::get<0>(ret) += Ynm_rhon1 * std::get<0>(this->MM_[ind]);
      std::get<1>(ret) += Ynm_rhon1 * std::get<1>(this->MM_[ind]);
    }
    //   return {std::get<0>(ret).real(), std::get<1>(ret).real()};
    return {std::get<0>(ret).real(), std::get<1>(ret).real()};
  }

  //$ ----------------------------------- M2M ---------------------------------- */

  /*DOC_EXTRACT M

  ### 定数の読み込み

  `AAA_M2M_FMM`，`AAA_M2L_FMM`，`AAA_L2L_FMM`は．`定数.nb`あらかじめ計算しておいたものを読み込む．

  */

  // AppleClang/libc++ struggles with constant-evaluating the huge precomputed tables (AAA_*),
  // and with constexpr std::vector in this context. Fall back to runtime initialization.
  std::vector<std::tuple<std::unique_ptr<SphericalCoordinates>, Moments*>> m2m_cache;

  //! これは計算毎に変わる
  template <typename TYPE>
  void set_m2m(const std::vector<TYPE>& children) {
    this->m2m_cache.resize(children.size()); // これでOK
    int i = 0;
    for (const auto& b : children) {
      //  auto sph = std::make_unique<SphericalCoordinates>(b->X - this->X);
      if (std::get<0>(this->m2m_cache[i]) == nullptr) {
        auto sph = std::make_unique<SphericalCoordinates>((b->X - this->X));
        sph->precompute_sph_rho(N);
        this->m2m_cache[i] = {std::move(sph), &b->MomentsMultipoleExpansion};
      } else {
        std::get<0>(this->m2m_cache[i])->initialize((b->X - this->X));
        std::get<0>(this->m2m_cache[i])->precompute_sph_rho(N);
        std::get<1>(this->m2m_cache[i]) = &b->MomentsMultipoleExpansion;
      }
      i++;
    }
  }

  void m2m() {
    const bool realfield_m_conj = bem_fmm_use_realfield_m_conj();
    cmplx R;
    for (const auto& [sph, moments_class] : this->m2m_cache) {
      for (const auto& op : m2m_ops()) {
        if (realfield_m_conj && op.k < 0)
          continue; // (real-field) MM_(j,-k) is reconstructed by conjugation
        // R = AAA * sph->sph_harmonics_rho(n, -m);
        R = sph->m2mFunction(op.j, op.k, op.n, op.m);
        auto& dst = this->MM_[op.dst];
        const auto& src = moments_class->MM_[op.src];
        std::get<0>(dst) += R * std::get<0>(src);
        std::get<1>(dst) += R * std::get<1>(src);
      }
    }
    if (realfield_m_conj)
      this->enforce_conjugate_m_symmetry_if_enabled();
  }

  //! 4重和のインデックスとそれがさす自身Moment係数へのポインターを，あらかじめ計算しておく．
  //% DFT畳み込み用
  //% PlaneWaveM2L: Plane-wave (exponential) expansion based M2L (diagonal translation)

#if !defined(SimpleM2L) && !defined(FourierM2L_double) && !defined(FourierM2L_DoubleDouble) && !defined(FourierM2L_BlockDecomposition) && !defined(USE_PlaneWaveM2L)
#define SimpleM2L
#endif

#if defined(SimpleM2L)
  using TREATED = double;
#elif defined(FourierM2L_double)
  using TREATED = double;
#elif defined(FourierM2L_DoubleDouble)
  using TREATED = DoubleDouble;
#elif defined(FourierM2L_BlockDecomposition)
  using TREATED = double;
//% --------------------------------
#define N_block 4
#elif defined(USE_PlaneWaveM2L)
  using TREATED = double;
  // Full plane-wave expansion with quadrature tables (Greengard & Rokhlin 1997)
  // Uses 6-direction decomposition and realfield optimization
  using PlaneWaveM2LConverter = PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>;
#endif

  // using TREATED = std::float128_t;
  using INPUT = std::complex<double>; //! これは固定とする．なぜならstd::float128_tならこれで問題なく収束したから
  using INPUT2D = std::array<std::array<INPUT, 2 * N + 1>, N + 1>;
  using fourier2DConvolutionINPUT = Fourier2DConvolution<INPUT, double, N + 1, 2 * N + 1, N + 1, 2 * N + 1>;
  using fourier2DConvolution = Fourier2DConvolution<INPUT, TREATED, N + 1, 2 * N + 1, N + 1, 2 * N + 1>;

#if defined(FourierM2L_BlockDecomposition)
  std::array<std::vector<std::tuple<std::shared_ptr<fourier2DConvolution>, Buckets<T, N>*>>, N_block> Yqp_Onm0_Onm1_cache_blocked;
#endif
  std::vector<std::tuple<std::shared_ptr<fourier2DConvolutionINPUT>, Buckets<T, N>*>> Yqp_Onm0_Onm1_cache;
  fourier2DConvolutionINPUT Onm0_DFT_M2L, Onm1_DFT_M2L;
  using cd = std::complex<double>;
  struct M2LTerm {
    cd AAAY;
    const std::array<cd, 2>* src_MM;
    const void* src_bucket; // source bucket pointer (for Metal M2L setup)
    int src_coeff_idx;      // index within source bucket's MM_ array
  };
  struct M2LRow {
    std::array<cd, 2>* dst_MM;
    std::size_t offset;
    std::size_t len;
  };
  std::vector<M2LRow> m2l_rows;
  std::vector<M2LTerm> m2l_terms;

  /* --------------------------------------------------------------------------- */
  // \label{Moments::set_m2l}

  template <typename TYPE>
  void set_m2l(std::vector<TYPE> children) {
    const bool realfield_m_conj = bem_fmm_use_realfield_m_conj();
    std::array<std::vector<M2LTerm>, (N + 1) * (N + 1)> terms_per_row{};

    this->Yqp_Onm0_Onm1_cache.assign(children.size(), {nullptr, nullptr});

    INPUT2D Ypq{};
    // Yqp_Onm0_Onm1_cache_blocked
#if defined(FourierM2L_BlockDecomposition)
    for (auto& v : this->Yqp_Onm0_Onm1_cache_blocked)
      v.clear();
    for (auto& v : this->Yqp_Onm0_Onm1_cache_blocked)
      v.reserve(children.size() / N_block + 1);
#elif defined(USE_PlaneWaveM2L)
    // Full plane wave M2L with quadrature tables (no initialization needed)
    this->pw_m2l_cache.clear();
    this->pw_m2l_cache.reserve(children.size());
#endif
    for (auto& v : terms_per_row)
      v.clear();
    // 削減率を出そう．
    int total_operations = 0;
    int total_reduced_operations = 0;
    const double threshold_AAAY = 0.; // 係数の閾値,スレッショルド
    std::sort(children.begin(), children.end(), [X = this->X](const auto& lhs, const auto& rhs) { return Norm(lhs->X - X) > Norm(rhs->X - X); });
    //
    int i = 0;
    for (auto& b : children) {
      auto sph = SphericalCoordinates(b->MomentsMultipoleExpansion.X - this->X);
      sph.precompute_sph_div_rhon1(2 * N);
#if defined(FourierM2L_double) || defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition)
      sph.precompute_sph_Ypq(N);
#endif
      //! 4重和に対応
      for (const auto& op : m2l_ops()) {
        if (realfield_m_conj && op.k < 0)
          continue; // (real-field) MM_(j,-k) is reconstructed by conjugation
        auto AAAY = sph.m2lFunction(op.j, op.k, op.n, op.m);
        if (threshold_AAAY < std::abs(AAAY)) {
          terms_per_row[op.dst].push_back({AAAY, &b->MomentsMultipoleExpansion.MM_[op.src], b, op.src});
          total_reduced_operations++;
        }
        total_operations++;
      }
      //$ ===================================== */
#if defined(FourierM2L_double) || defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition)
      for (int p = 0; p <= N; ++p)
        for (int q = -p; q <= p; ++q)
          Ypq[p][q + N] = sph.Ypq(p, -q);
#endif
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
      if (std::get<0>(this->Yqp_Onm0_Onm1_cache[i]) == nullptr)
        this->Yqp_Onm0_Onm1_cache[i] = {std::make_shared<fourier2DConvolutionINPUT>(Ypq, true, false), b};
      else {
        std::get<0>(this->Yqp_Onm0_Onm1_cache[i])->initialize(Ypq, true, false);
        std::get<1>(this->Yqp_Onm0_Onm1_cache[i]) = b;
      }
#elif defined(FourierM2L_BlockDecomposition)
      std::array<INPUT2D, N_block> Yqp_blocked{};
      for (int p = 0; p <= N; ++p)
        for (int q = -p; q <= p; ++q)
          Yqp_blocked[getBlockIndex(p)][p][q + N] = Ypq[p][q + N];
      for (int i = 0; i < N_block; ++i)
        this->Yqp_Onm0_Onm1_cache_blocked[i].emplace_back(std::make_shared<fourier2DConvolutionINPUT>(Yqp_blocked[i], true, false), b);
#elif defined(USE_PlaneWaveM2L)
      // Store translation vector for plane-wave M2L
      // Note: In M2L, this->X is target local expansion center, b->X is source multipole center
      std::array<double, 3> translation = {
          this->X[0] - b->MomentsMultipoleExpansion.X[0],
          this->X[1] - b->MomentsMultipoleExpansion.X[1],
          this->X[2] - b->MomentsMultipoleExpansion.X[2]};
      this->pw_m2l_cache.push_back({b, translation});
#endif
      //$ ===================================== */
      i++;
    }

    std::size_t total_terms = 0;
    for (const auto& v : terms_per_row)
      total_terms += v.size();
    this->m2l_rows.clear();
    this->m2l_rows.reserve((N + 1) * (N + 1));
    this->m2l_terms.clear();
    this->m2l_terms.reserve(total_terms);
    std::size_t offset = 0;
    for (std::size_t local_index = 0; local_index < terms_per_row.size(); ++local_index) {
      auto& v = terms_per_row[local_index];
      if (v.empty())
        continue;
      // Optional reordering (see bem_simplem2l_term_order()).
      switch (bem_simplem2l_term_order()) {
      case SimpleM2LTermOrder::None:
        break;
      case SimpleM2LTermOrder::Src:
        std::stable_sort(v.begin(), v.end(), [](const M2LTerm& lhs, const M2LTerm& rhs) { return std::less<const void*>{}(lhs.src_MM, rhs.src_MM); });
        break;
      case SimpleM2LTermOrder::Abs:
      default:
        std::sort(v.begin(), v.end(), [](const M2LTerm& lhs, const M2LTerm& rhs) { return std::abs(lhs.AAAY) < std::abs(rhs.AAAY); });
        break;
      }
      this->m2l_rows.push_back({&this->MM_[local_index], offset, v.size()});
      this->m2l_terms.insert(this->m2l_terms.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
      offset += v.size();
    }
  }

  /* --------------------------------------------------------------------------- */
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
  Fourier2DConvolution<INPUT, double, N + 1, 2 * N + 1, N + 1, 2 * N + 1> accum_Onm0Yqp, accum_Onm1Yqp;
#elif defined(FourierM2L_BlockDecomposition)
  Fourier2DConvolution<INPUT, double, N + 1, 2 * N + 1, N + 1, 2 * N + 1> accum_Onm0Yqp, accum_Onm1Yqp;
  std::array<fourier2DConvolution, N_block> Onm0_DFT_M2L_Blocked, Onm1_DFT_M2L_Blocked;
#elif defined(USE_PlaneWaveM2L)
  // Cache of source boxes and their translation vectors for plane-wave M2L
  struct PlaneWaveM2LEntry {
    Buckets<T, N>* source_bucket;
    std::array<double, 3> translation; // target.X - source.X
  };
  std::vector<PlaneWaveM2LEntry> pw_m2l_cache;
#endif

  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */

#if defined(FourierM2L_double) || defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition)
  void update_m2l_cache_realfield() {
    // Real-field optimization:
    // When bem_fmm_use_realfield_m_conj() is enabled, MM_(n,-m) == conj(MM_(n,m)).
    // Here A(n,m) is the same for ±m (depends only on |m|), and is purely real for even |m|
    // and purely imaginary for odd |m|. Therefore:
    //   Onm(n,-m) = (-1)^m * conj(Onm(n,m))
    //
    // This lets us build the full (n,m) triangular region from m>=0 only.
    INPUT2D Onm0{}, Onm1{};
    std::complex<TREATED> A;
    for (int n = 0; n <= N; ++n) {
      // m = 0
      {
        const int m = 0;
        A = static_cast<std::complex<TREATED>>(i_absk_A_FMM[n][m + N_i_absk_A_FMM]) * static_cast<TREATED>(1 - ((n & 1) << 1));
        const auto mm = this->MM_[this->index(n, m)];
        Onm0[n][m + N] = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<0>(mm)));
        Onm1[n][m + N] = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<1>(mm)));
      }
      for (int m = 1; m <= n; ++m) {
        A = static_cast<std::complex<TREATED>>(i_absk_A_FMM[n][m + N_i_absk_A_FMM]) * static_cast<TREATED>(1 - ((n & 1) << 1));
        const auto mm_pos = this->MM_[this->index(n, m)];
        const auto v0 = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<0>(mm_pos)));
        const auto v1 = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<1>(mm_pos)));
        Onm0[n][m + N] = v0;
        Onm1[n][m + N] = v1;

        const auto v0c = std::conj(v0);
        const auto v1c = std::conj(v1);
        if (m & 1) {
          Onm0[n][-m + N] = -v0c;
          Onm1[n][-m + N] = -v1c;
        } else {
          Onm0[n][-m + N] = v0c;
          Onm1[n][-m + N] = v1c;
        }
      }
    }

#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
    this->Onm0_DFT_M2L.clear();
    this->Onm1_DFT_M2L.clear();
    this->Onm0_DFT_M2L.add(Onm0);
    this->Onm1_DFT_M2L.add(Onm1);
#elif defined(FourierM2L_BlockDecomposition)
    /* ------------------------ distribute into 4 blocks ------------------------ */
    std::array<INPUT2D, N_block> Onm0Blocked{}, Onm1Blocked{};
    for (int n = 0; n <= N; ++n) {
      int index = getBlockIndex(n);
      for (int m = -n; m <= n; ++m) {
        Onm0Blocked[index][n][m + N] = Onm0[n][m + N];
        Onm1Blocked[index][n][m + N] = Onm1[n][m + N];
      }
    }
    for (int i = 0; i < N_block; ++i) {
      this->Onm0_DFT_M2L_Blocked[i].clear();
      this->Onm1_DFT_M2L_Blocked[i].clear();
      this->Onm0_DFT_M2L_Blocked[i].add(Onm0Blocked[i]);
      this->Onm1_DFT_M2L_Blocked[i].add(Onm1Blocked[i]);
    }
    /* --------------------------------------------------------------------------- */
#endif
  }
#endif

  void update_m2l_cache() {
#if defined(FourierM2L_double) || defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition)
    if (bem_fmm_use_realfield_m_conj()) {
      this->update_m2l_cache_realfield();
      return;
    }
    // Must be zero-initialized: outside the triangular (n,|m|<=n) region coefficients are implicitly 0.
    // Leaving them uninitialized contaminates the DFT and breaks accuracy (especially in block mode).
    INPUT2D Onm0{}, Onm1{};
    std::complex<TREATED> A;
    for (int n = 0; n <= N; ++n)
      for (int m = -n; m <= n; ++m) {
        A = static_cast<std::complex<TREATED>>(i_absk_A_FMM[n][m + N_i_absk_A_FMM]) * static_cast<TREATED>(1 - ((n & 1) << 1));
        auto mm = this->MM_[this->index(n, m)];
        Onm0[n][m + N] = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<0>(mm)));
        Onm1[n][m + N] = static_cast<INPUT>(A * static_cast<std::complex<TREATED>>(std::get<1>(mm)));
      }
#endif
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
    this->Onm0_DFT_M2L.clear();
    this->Onm1_DFT_M2L.clear();
    this->Onm0_DFT_M2L.add(Onm0);
    this->Onm1_DFT_M2L.add(Onm1);
#elif defined(FourierM2L_BlockDecomposition)
    /* ------------------------ distribute into 4 blocks ------------------------ */
    std::array<INPUT2D, N_block> Onm0Blocked{}, Onm1Blocked{};
    for (int n = 0; n <= N; ++n) {
      int index = getBlockIndex(n);
      for (int m = -n; m <= n; ++m) {
        Onm0Blocked[index][n][m + N] = Onm0[n][m + N];
        Onm1Blocked[index][n][m + N] = Onm1[n][m + N];
      }
    }
    for (int i = 0; i < N_block; ++i) {
      this->Onm0_DFT_M2L_Blocked[i].clear();
      this->Onm1_DFT_M2L_Blocked[i].clear();
      this->Onm0_DFT_M2L_Blocked[i].add(Onm0Blocked[i]);
      this->Onm1_DFT_M2L_Blocked[i].add(Onm1Blocked[i]);
    }
    /* --------------------------------------------------------------------------- */
#endif
  }

  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */
  //% -------------------------------------------------------------------------- */

  void m2l() {
#ifdef SimpleM2L
    for (const auto& row : this->m2l_rows) {
      double sum0_re = 0.0, sum0_im = 0.0;
      double sum1_re = 0.0, sum1_im = 0.0;
      const M2LTerm* terms = this->m2l_terms.data() + row.offset;
#if defined(_OPENMP)
#pragma omp simd reduction(+ : sum0_re, sum0_im, sum1_re, sum1_im)
#endif
      for (std::size_t i = 0; i < row.len; ++i) {
        const auto& t = terms[i];
        const cd& a = t.AAAY;
        const cd& b0 = (*t.src_MM)[0];
        const cd& b1 = (*t.src_MM)[1];

        const double ar = a.real();
        const double ai = a.imag();

        const double b0r = b0.real();
        const double b0i = b0.imag();
        const double b1r = b1.real();
        const double b1i = b1.imag();

        sum0_re += ar * b0r - ai * b0i;
        sum0_im += ar * b0i + ai * b0r;
        sum1_re += ar * b1r - ai * b1i;
        sum1_im += ar * b1i + ai * b1r;
      }
      (*row.dst_MM)[0] += cd(sum0_re, sum0_im);
      (*row.dst_MM)[1] += cd(sum1_re, sum1_im);
    }
    // (real-field) fill m<0 from m>0
    this->enforce_conjugate_m_symmetry_if_enabled();

#elif defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition) || defined(FourierM2L_double)
    //$ ===================================== */
    //  clear this->MM_
    this->MM_.fill({});

    // \label{DFT2D_Onm0_Onm1_Yqp}
    accum_Onm0Yqp.clear();
    accum_Onm1Yqp.clear();

    /* -------------------------------------------------------------------------- */
#if defined(FourierM2L_DoubleDouble) || defined(FourierM2L_double)
    for (const auto& [Ypq_ptr, b] : this->Yqp_Onm0_Onm1_cache) {
      accum_Onm0Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm0_DFT_M2L.MatrixDFT);
      accum_Onm1Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm1_DFT_M2L.MatrixDFT);
    }

    if constexpr (std::same_as<INPUT, std::complex<DoubleDouble>>) {
      accum_Onm0Yqp.copyDoubleDoubleToDouble();
      accum_Onm1Yqp.copyDoubleDoubleToDouble();
    }

    accum_Onm0Yqp.convolve();
    accum_Onm1Yqp.convolve();
    /* ------------------------------- block版あくまでテスト ----------------------- */
#elif defined(FourierM2L_BlockDecomposition)
    // Block decomposition must preserve the full bilinear product:
    //   DFT(Y) * DFT(M) = (sum_i DFT(Y_i)) * (sum_j DFT(M_j)) = sum_{i,j} DFT(Y_i) * DFT(M_j)
    // Dropping pairs (e.g. by restricting to i+j==k with k in [0,N_block)) loses information and degrades accuracy.
    for (int i = 0; i < N_block; ++i)
      for (const auto& [Ypq_ptr, b] : this->Yqp_Onm0_Onm1_cache_blocked[i])
        for (int j = 0; j < N_block; ++j) {
          accum_Onm0Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm0_DFT_M2L_Blocked[j].MatrixDFT);
          accum_Onm1Yqp.Hadamard_product_add(Ypq_ptr->MatrixDFT, b->MomentsMultipoleExpansion.Onm1_DFT_M2L_Blocked[j].MatrixDFT);
        }

    accum_Onm0Yqp.convolve();
    accum_Onm1Yqp.convolve();
#endif
    /* -------------------------------------------------------------------------- */
    std::complex<double> A; //! ここはdoubleでも良いようだ．

    for (int j = N; j >= 0; --j) {
      for (int k = -j; k <= j; ++k) {
        A = i_absk_A_FMM[j][k + N_i_absk_A_FMM];
        auto& MM = this->MM_[this->index(j, k)];
        std::get<0>(MM) += static_cast<cd>(A * accum_Onm0Yqp.convolution_T[k + 2 * N][N - j]);
        std::get<1>(MM) += static_cast<cd>(A * accum_Onm1Yqp.convolution_T[k + 2 * N][N - j]);
      }
    }

    //$ ===================================== */
#elif defined(USE_PlaneWaveM2L)
    //$ ===================================== */
    // Plane-wave (exponential expansion) based M2L
    // This achieves O(p^2) translation instead of O(p^4)
    //
    // Algorithm:
    // 1. For each source in interaction list:
    //    a. Convert source multipole to exponential expansion (C_MX): O(p^3)
    //    b. Translate exponential expansion diagonally: O(p^2) -- the key advantage!
    //    c. Accumulate into target exponential expansion
    // 2. Convert accumulated exponential to local expansion (C_XL): O(p^3)
    //$ ===================================== */
    const bool realfield_m_conj = bem_fmm_use_realfield_m_conj();

    // Clear this->MM_ before accumulating M2L contributions
    this->MM_.fill({});

    // Process each source in the M2L interaction list
    for (const auto& entry : this->pw_m2l_cache) {
      auto& src_moments = entry.source_bucket->MomentsMultipoleExpansion;

      // Full plane wave M2L with 6-direction decomposition
      // Automatically handles rotation and uses optimal quadrature
      Tddd translation_vec = {entry.translation[0], entry.translation[1], entry.translation[2]};
      // Type conversion: std::array<std::array<cmplx,2>,N> -> std::array<std::tuple<cmplx,cmplx>,N>
      // These have identical memory layout, so reinterpret_cast is safe
      using PWE_MM_Type = typename PlaneWaveM2LConverter::MM_Type;
      using PWE_L_Type = typename PlaneWaveM2LConverter::L_Type;
      PlaneWaveM2LConverter::translateMultipoleToLocal(
          reinterpret_cast<const PWE_MM_Type&>(src_moments.MM_),
          reinterpret_cast<PWE_L_Type&>(this->MM_),
          translation_vec, realfield_m_conj);
    }

    // (real-field) fill m<0 from m>0
    this->enforce_conjugate_m_symmetry_if_enabled();
    //$ ===================================== */
#endif
  }

  //^ ----------------------------------- L2L ---------------------------------- */

  std::tuple<std::unique_ptr<SphericalCoordinates>, Moments*> l2l_cache;

  //! これは計算毎に変わる
  void set_l2l(Moments* MomentsLocalExpansion) {
    //   auto sph = std::make_unique<SphericalCoordinates>(MomentsLocalExpansion->X - this->X);
    auto sph = std::make_unique<SphericalCoordinates>((MomentsLocalExpansion->X - this->X));
    sph->precompute_sph_rho(N);
    this->l2l_cache = {std::move(sph), MomentsLocalExpansion};
  }

  void l2l() {
    const bool realfield_m_conj = bem_fmm_use_realfield_m_conj();
    cmplx R;
    for (const auto& op : l2l_ops()) {
      if (realfield_m_conj && op.k < 0)
        continue; // (real-field) MM_(j,-k) is reconstructed by conjugation
      // R = AAA * std::get<0>(l2l_cache)->sph_harmonics_rho(n - j, m - k);
      R = std::get<0>(l2l_cache)->l2lFunction(op.j, op.k, op.n, op.m);
      auto& dst = this->MM_[op.dst];
      const auto& src = std::get<1>(l2l_cache)->MM_[op.src];
      std::get<0>(dst) += R * std::get<0>(src);
      std::get<1>(dst) += R * std::get<1>(src);
    }
    if (realfield_m_conj)
      this->enforce_conjugate_m_symmetry_if_enabled();
  }
};

//! ========================================================================== */
//! ========================================================================== */
//! ========================================================================== */

#include "lib_spatial_partitioning.hpp"

//! ========================================================================== */
//! ========================================================================== */
//! ========================================================================== */

//% -------------------------------------------------------------------------- */
//%        Functions Related to Multipole Expansion Utilizing Buckets          */
//% -------------------------------------------------------------------------- */
// #define _DEBUG_FMM_
inline void setM2M(auto& B_poles) {
  TimeWatch tw, tw_all;
  for (int level = B_poles.max_level - 1; level >= 0; level--) {
    B_poles.forEachAtLevel({level}, [&](auto* B) { B->MomentsMultipoleExpansion.set_m2m(B->getAllChildren()); });
#if defined(_DEBUG_FMM_)
    std::cout << yellow << "setM2M" << ", level=" << level << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
#endif
  }
  std::cout << Yellow << "setM2M" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
}

inline void M2M(auto& B_poles) {
  for (int level = B_poles.max_level - 1; level >= 0; level--)
    B_poles.forEachAtLevel({level}, [&](auto* B) { B->MomentsMultipoleExpansion.m2m(); });
};

/* -------------------------------------------------------------------------- */

inline void setL2L(auto& B_poles) {
  for (auto& buckets_from_top_level : B_poles.level_buckets)
    for (auto& parent : buckets_from_top_level)
      parent->traverseChildren([&](auto* child) { child->MomentsLocalExpansion.set_l2l(&parent->MomentsLocalExpansion); });
}

inline void L2L(auto& B_poles) {
  for (auto& buckets_from_top_level : B_poles.level_buckets)
    for (auto& parent : buckets_from_top_level)
      parent->traverseChildren([&](auto* child) { child->MomentsLocalExpansion.l2l(); });
}

/* -------------------------------------------------------------------------- */

// Box-based MAC (旧実装 - 一時的に有効化)
// scaledBoundsが返すのは，
// xmin <-- L -->xmax　として，
// xmin <-- L/2 -- xcenter -- L/2 -->xmax
// xmin <-- scale*L/2 -- xcenter -- scale*L/2 -->xmax
// xmin <-- scale*L -->xmax と結果的になる．

// inline const double scale = 3.;
// inline bool isFar(auto *A, auto *B) { return !InsideQ(B->X, A->scaledBounds(scale)); }
// inline bool isNear(auto *A, auto *B) { return InsideQ(B->X, A->scaledBounds(scale)); }

// Distance-based MAC (新実装 - 一時的に無効化):
// セル中心間距離 > mac_distance_threshold * dL なら Far (M2L対象)
// θ = √3·L / R において、R = 3.1L のとき θ ≈ 0.558 < 1 で収束保証

extern double g_mac_theta;             // settings.json の "mac_theta" から設定（デフォルト 0.25）
extern int g_p2m_quadrature_points;    // settings.json の "p2m_quadrature_points" から設定（デフォルト 6）
inline bool isFar(auto* A, auto* B) {
  // 内部ソースのはみ出し分をバケツの有効サイズに加算
  double effective_dL = std::max(A->dL + 2.0 * A->max_internal_offset,
                                 B->dL + 2.0 * B->max_internal_offset);
  return g_mac_theta > std::sqrt(3.0) * effective_dL / Norm(B->X - A->X);
}
inline bool isNear(auto* A, auto* B) { return !isFar(A, B); }

// inline bool isFar(const auto *A, const auto *B);  // 不要な前方宣言

/*
M2Lの対象の条件
1. 同じレベル
2. 十分離れているisFar==true, また直接近傍ではないisNear==false
3. 親同士はisNear==true
*/

template <typename T>
void buildNearAndM2LLists(T A, T B) {
  // isNear(A,B) && isFar(A->c,B->c)となるA->cとB->cを探し，A->c->buckets_for_M2LにB->cを追加する．
  if (isNear(A, B)) {
    if (A != B)
      A->buckets_near.emplace_back(B);
    A->traverseChildren([&](auto* A_c) { B->traverseChildren([&](auto* B_c) {
        if (A_c != B_c && isFar(A_c, B_c))
          A_c->buckets_for_M2L.emplace_back(B_c);
        else
          buildNearAndM2LLists(A_c, B_c); }); });
  }
}

// 極が少ないバケツは，M2Lの時に省略する．

inline void setBucketsForM2L(auto& B_poles) {

  B_poles.traverseTree([&](auto* A) {
    A->buckets_for_M2L.clear();
    A->buckets_near.clear();
    A->buckets_for_L2M.clear(); });

  B_poles.traverseChildren([&](auto* A) { B_poles.traverseChildren([&](auto* B) { buildNearAndM2LLists(A, B); }); });

  for (auto& buckets_from_top_level : B_poles.level_buckets) {
    for (auto& A : buckets_from_top_level)
      for (auto& B : A->buckets_for_M2L)
        B->buckets_for_L2M.emplace_back(A);
  }
}

inline void setM2L(auto& B_poles) {
  TimeWatch tw, tw_all;
  // std::cout << "各レベルの各セルのM2Lの相手を保存する" << std::endl;
  setBucketsForM2L(B_poles);
  // A -> M2L -> B
  int level = 0;
  for (auto& buckets_at_a_level : B_poles.level_buckets) {
#pragma omp parallel for
    for (auto& A : buckets_at_a_level)
      A->MomentsLocalExpansion.set_m2l(A->buckets_for_L2M);
    std::cout << yellow << "setM2L" << ", level=" << level << ", buckets_at_a_level.size()=" << buckets_at_a_level.size() << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
    level++;
  }
  std::cout << Yellow << "setM2L" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
};

/* -------------------------------------------------------------------------- */
/* Metal M2L GPU acceleration (optional)                                       */
/* -------------------------------------------------------------------------- */
#if defined(USE_METAL_M2L)

// Global Metal M2L wrapper (lazy initialization)
inline std::unique_ptr<MetalM2LWrapper> g_metal_m2l_wrapper;

// Runtime enable/disable flag for Metal M2L
// This can be toggled during GMRES to switch from GPU (float+Kahan) to CPU (double) M2L
inline bool g_metal_m2l_active = true;

// Check if Metal M2L is available and enabled
inline bool isMetalM2LAvailable() {
  return metal_m2l_is_available() != 0;
}

// Check if Metal M2L is currently in use (wrapper initialized AND active flag set)
inline bool isMetalM2LInUse() {
  return g_metal_m2l_wrapper && g_metal_m2l_active;
}

// Global Metal M2L settings (read from settings.json)
inline bool g_metal_m2l_threadgroup = false; // true: use threadgroup parallelization
inline bool g_metal_m2l_sort_terms = false;  // Sort terms for improved memory locality

// Initialize Metal M2L from settings (called from BEM solver)
// Parameters:
//   enabled: use_metal_m2l from settings.json
//   threadgroup: metal_m2l_threadgroup from settings.json (true/false)
//   sort_terms: metal_m2l_sort_terms from settings.json
inline void initMetalM2L(bool enabled, bool threadgroup = false, bool sort_terms = false) {
  if (g_metal_m2l_wrapper)
    return; // Already initialized

  if (!enabled) {
    std::cout << "[Metal M2L] Disabled (set use_metal_m2l=true in settings.json to enable)" << std::endl;
    return;
  }

  if (!isMetalM2LAvailable()) {
    std::cerr << "[Metal M2L] Not available on this system" << std::endl;
    return;
  }

  // Store settings for later use
  g_metal_m2l_threadgroup = threadgroup;
  g_metal_m2l_sort_terms = sort_terms;

  try {
    g_metal_m2l_wrapper = std::make_unique<MetalM2LWrapper>(threadgroup);
    std::cout << "[Metal M2L] Initialized [float+Kahan mode, TG=" << threadgroup << "]" << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "[Metal M2L] Initialization failed: " << e.what() << std::endl;
    g_metal_m2l_wrapper.reset();
  }
}

// Setup Metal M2L with bucket data (call after setM2L)
template <typename BPoles>
inline void setupMetalM2L(BPoles& B_poles) {
  if (!g_metal_m2l_wrapper)
    return;

  TimeWatch tw;

  // Collect all buckets from all levels
  using BucketType = std::decay_t<decltype(*B_poles.level_buckets[0][0])>;
  std::vector<BucketType*> all_buckets;

  for (auto& buckets_at_level : B_poles.level_buckets) {
    for (auto& bucket : buckets_at_level) {
      all_buckets.push_back(bucket);
    }
  }

  try {
    g_metal_m2l_wrapper->setup(all_buckets, g_metal_m2l_sort_terms);
    std::cout << Yellow << "setupMetalM2L" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "[Metal M2L] Setup failed: " << e.what() << std::endl;
    g_metal_m2l_wrapper.reset();
  }
}

// Execute Metal M2L (call instead of CPU M2L)
template <typename BPoles>
inline void metalM2L(BPoles& B_poles) {
  if (!g_metal_m2l_wrapper)
    return;

  // Collect all buckets
  using BucketType = std::decay_t<decltype(*B_poles.level_buckets[0][0])>;
  std::vector<BucketType*> all_buckets;
  for (auto& buckets_at_level : B_poles.level_buckets) {
    for (auto& bucket : buckets_at_level) {
      all_buckets.push_back(bucket);
    }
  }

  // Copy source MM values to GPU
  g_metal_m2l_wrapper->copySourceMM(all_buckets);

  // Execute GPU compute
  g_metal_m2l_wrapper->compute();

  // Write results back to bucket arrays
  g_metal_m2l_wrapper->writeResultsBack();

  // Handle realfield conjugate symmetry if enabled
  if (bem_fmm_use_realfield_m_conj()) {
    for (auto& buckets_at_level : B_poles.level_buckets) {
      for (auto& bucket : buckets_at_level) {
        bucket->MomentsLocalExpansion.enforce_conjugate_m_symmetry_if_enabled();
      }
    }
  }
}

#endif // USE_METAL_M2L

/* -------------------------------------------------------------------------- */

inline void M2L(auto& B_poles) {
#if defined(USE_METAL_M2L)
  // Use Metal GPU acceleration if available and active
  // g_metal_m2l_active can be set to false during GMRES to switch to CPU (double) M2L
  if (g_metal_m2l_wrapper && g_metal_m2l_active) {
    metalM2L(B_poles);
    return;
  }
#endif

  // CPU fallback
#if defined(_OPENMP)
  // BEM_OMP_M2L_SCHEDULE: "dynamic,1" (default), "static", "guided,8", etc.
  static const auto m2l_schedule = parse_omp_schedule(std::getenv("BEM_OMP_M2L_SCHEDULE"), omp_sched_dynamic, 1);
  omp_set_schedule(m2l_schedule.kind, m2l_schedule.chunk);
#endif
  // Note:
  // - Putting `#pragma omp parallel` outside reduces the overhead of repeatedly launching thread teams
  //   for each level.
  // - Keeping the work grouped per level tends to have better cache locality than fully flattening.
#if defined(_OPENMP)
#pragma omp parallel
  {
    // update_m2l_cache() is only needed for Fourier-based M2L (DFT preparation).
    // For SimpleM2L it is a no-op; skip the full-tree traversal to avoid overhead.
    // Bench (DeepCwind, 3 steps avg; M2L per mat-vec; static=1.00, lower is better):
    // +-------------+--------+
    // | schedule    | ratio  |
    // +-------------+--------+
    // | static      | 1.00   |
    // | dynamic,1   | 0.83   |
    // | dynamic,8   | 0.94   |
    // | guided,8    | 0.95   |
    // +-------------+--------+
#if defined(FourierM2L_double) || defined(FourierM2L_DoubleDouble) || defined(FourierM2L_BlockDecomposition)
    for (auto& buckets_at_a_level : B_poles.level_buckets) {
#pragma omp for schedule(runtime)
      for (auto& A : buckets_at_a_level)
        A->MomentsMultipoleExpansion.update_m2l_cache();
    }
#pragma omp barrier
#endif

    for (auto& buckets_at_a_level : B_poles.level_buckets) {
#pragma omp for schedule(runtime)
      for (auto& A : buckets_at_a_level)
        A->MomentsLocalExpansion.m2l();
    }
  }
#else
  for (auto& buckets_at_a_level : B_poles.level_buckets)
    for (auto& A : buckets_at_a_level)
      A->MomentsLocalExpansion.m2l();
#endif
};

/* -------------------------------------------------------------------------- */

/*
source4FMMは，
Tddd X;
Tdd weights;
Tddd normal;
std::function<Tdd()> getValues;
を持つ．これらを使うと，直接数値面積分ができる．
*/

/*

BIE：alpha_phi = w_Gn_phin - w_Gn_phiの，w_Gn_phin，w_Gn_phiをまず計算する．

ここで，w_Gn_phi，w_Gn_phiはそれぞれ，

$$
\begin{align*}
\int\int_{\Gamma} G \nabla \phi \cdot {\mathbfit n} dS\\
\int\int_{\Gamma} \phi \nabla G \cdot {\mathbfit n} dS
\end{align*}
$$

w_Gn_phiのincrementにおいて，符号が負になっている．
これは，勘違いして修正してしまいそうだが，w_Gn_phi自体が負になっているためで，
BIEの符号と混同しないように注意する．

BIEは次の形である．

alpha_phi = w_Gn_phin - w_Gn_phi

多重極展開（see L2P）も，w_Gn_phin，w_Gn_phiの近似を計算しており，直接積分しているものと対応している．

*/

/* -------------------------------------------------------------------------- */

// Metal M2L settings structure (read from settings.json)
struct MetalM2LSettings {
  bool enabled = false;
  bool threadgroup = false; // true: use threadgroup parallelization
  bool sort_terms = false;
};

inline void initializeFMM(auto& B_poles, const auto& targets, bool reuse_static_tree_cache,
                          MetalM2LSettings metal_m2l_settings = {},
                          const std::string& nearfield_mode = "scalar") {

  TimeWatch log_stage_watch;
  auto log_stage = [&](const std::string& label) { if (!bem_verbose()) return; std::cout << Magenta << " [FMM:init] " << Cyan << label << Green << " elapsed=" << log_stage_watch() << colorReset << std::endl; };

// M2Lの方法を表示
#if defined(SimpleM2L)
  std::cout << "SimpleM2L" << std::endl;
#elif defined(FourierM2L_double)
  std::cout << "FourierM2L_double" << std::endl;
#elif defined(FourierM2L_DoubleDouble)
  std::cout << "FourierM2L_DoubleDouble" << std::endl;
#elif defined(FourierM2L_BlockDecomposition)
  std::cout << "FourierM2L_BlockDecomposition" << std::endl;
#elif defined(USE_PlaneWaveM2L)
  std::cout << "PlaneWaveM2L (exponential expansion, diagonal translation)" << std::endl;
#endif

  // FMMパラメータを出力
  constexpr int N = std::decay_t<decltype(B_poles)>::expansion_order;
  {
    std::stringstream ss;
    ss << "[FMM:init] Expansion order N = " << N << " (coefficients: " << (N + 1) * (N + 1) << ")";
    log_stage(ss.str());
  }
  {
    std::stringstream ss;
    ss << "[FMM:init] MAC theta = " << g_mac_theta
       << " (far if R > " << std::sqrt(3.0) / g_mac_theta << " * dL)";
    log_stage(ss.str());
  }
  {
    std::stringstream ss;
    ss << "[FMM:init] P2M quadrature points = " << g_p2m_quadrature_points;
    log_stage(ss.str());
  }

  // ==========================================================================
  // Steps 0-2.5: ソース位置依存 → 毎回実行（ソースが移動する可能性があるため）
  // Steps 3-4:   バケット位置依存 → ツリー構造不変なら再利用可能
  // ==========================================================================

  // ------------------------------------------------------------------
  // Step 0: 事前準備
  //  - data1D_vector に入っている全 source の値を最新化
  //  - ここでの source->update() は，各 source が保持する (phi, phin) 等の動的量を再計算し，多重極係数生成に必要な weighted_source_densities を一貫させる目的。
  // ------------------------------------------------------------------
  {
    std::stringstream ss;
    ss << "\t[FMM:init] Step0: update poles start, #poles=" << B_poles.data1D_vector.size();
    log_stage(ss.str());
  }
  for (auto source : B_poles.data1D_vector)
    source->updateDensity();
  log_stage("Step0: update poles done");
  // ------------------------------------------------------------------
  // Step 1: 全ノードのモーメントをゼロ初期化
  //  - Multipole(Local) Moments の係数配列 MM_ をクリア
  //  - FMM の各フェーズ (M2M, M2L, L2L) を再構築する前に必須
  // ------------------------------------------------------------------
  {
    size_t total_nodes_before_reset = 0;
    for (auto& lv : B_poles.level_buckets)
      total_nodes_before_reset += lv.size();
    std::cout << "\t[FMM:init] Step1: reset moments for total nodes=" << total_nodes_before_reset << ", levels=" << B_poles.level_buckets.size() << std::endl;
  }
  B_poles.forEachAll([&](auto* b) {
    b->MomentsMultipoleExpansion.initialize(); //! MとM_の初期化
    b->MomentsLocalExpansion.initialize();     //! MとM_の初期化
  });
  log_stage("Step1: reset moments done");
  // ------------------------------------------------------------------
  // Step 2: 葉ノード(最深バケツ)で P2M (increment_M) を実行
  //  - 各葉が保持するソース点集合 data1D_vector から多重極係数を構築
  //  - increment_M は再利用キャッシュ(vector_weightedSourceDensities_RRn_MM_) を再生成
  //  - ソース位置 (isrc.X) が変わるとキャッシュ (RRn) が古くなるため，毎回再構築が必要
  //  - 並列 (par_unseq/OMP) で高速化
  // ------------------------------------------------------------------
  {
    size_t n_leaves = B_poles.deepest_level_buckets.size();
    size_t sources_in_leaves = 0;
    for (auto* leaf : B_poles.deepest_level_buckets)
      sources_in_leaves += leaf->data1D_vector.size();
    std::cout << "\t[FMM:init] Step2: P2M begin, #leaves=" << n_leaves << ", total leaf sources=" << sources_in_leaves << std::endl;
  }
  B_poles.forEachAtDeepestParallel([&](auto* b) {
    b->MomentsMultipoleExpansion.increment_M(b->data1D_vector); // !vector_weightedSourceDensities_RRn_MMの初期化
  });
  log_stage("Step2: P2M (increment_M) done");
  // ------------------------------------------------------------------
  // Step 2.5: 各バケツの max_internal_offset を計算（MAC拡張用）
  //  - 各リーフバケツ内の全要素の max_source_offset の最大値を集約
  //  - 親バケツへはボトムアップで伝播（子の max_internal_offset の最大値）
  // ------------------------------------------------------------------
  // リーフバケツ: 要素の max_source_offset から集約
  B_poles.forEachAtDeepest([&](auto* b) {
    b->max_internal_offset = 0.0;
    for (const auto& elem : b->data1D_vector)
      b->max_internal_offset = std::max(b->max_internal_offset, elem->max_source_offset);
  });
  // 親バケツへボトムアップ伝播（リーフでないバケツのみ）
  for (int level = B_poles.max_level - 1; level >= 0; level--)
    B_poles.forEachAtLevel({level}, [&](auto* b) {
      if (!b->hasChildren()) return; // リーフは Step2.5 前半で計算済み
      b->max_internal_offset = 0.0;
      b->traverseChildren([&](auto* child) {
        b->max_internal_offset = std::max(b->max_internal_offset, child->max_internal_offset);
      });
    });
  {
    double global_max_offset = 0.0;
    double leaf_dL = 0.0;
    for (auto* leaf : B_poles.deepest_level_buckets) {
      global_max_offset = std::max(global_max_offset, leaf->max_internal_offset);
      leaf_dL = leaf->dL;
    }
    std::stringstream ss;
    ss << "Step2.5: max_internal_offset = " << global_max_offset
       << ", leaf dL = " << leaf_dL
       << ", offset/dL = " << (leaf_dL > 0 ? global_max_offset / leaf_dL : 0.0);
    log_stage(ss.str());
  }

  // ==========================================================================
  // Steps 3-4: バケット間の相対位置に依存するオペレータ
  //  - M2M/M2L/L2L オペレータはバケット中心間の相対位置 (child->X - parent->X) で決まる
  //  - ツリー構造（バケット配置）が不変なら再利用可能
  // ==========================================================================
  if (reuse_static_tree_cache) {
    // M2M/L2L はバケツ中心間の相対位置のみに依存 → ツリー構造不変なら再利用可能
    log_stage("[FMM:init] reuse static M2M/L2L operators; skip M2M/L2L");
    // M2L ペアリストは max_internal_offset に依存 → 毎回現在値で再構築
    setM2L(B_poles);
    log_stage("Step3: setM2L rebuilt with current max_internal_offset");
#if defined(USE_METAL_M2L)
    initMetalM2L(metal_m2l_settings.enabled, metal_m2l_settings.threadgroup, metal_m2l_settings.sort_terms);
    if (g_metal_m2l_wrapper) {
      setupMetalM2L(B_poles);
      log_stage("Step3.5: setupMetalM2L rebuilt");
    }
#endif

    {
      std::cout << Magenta << " [FMM:stat] level | buckets |    dL     | mean_M2L | mean_near | max_offset/dL" << colorReset << std::endl;
      int level = 0;
      for (auto& buckets_at_a_level : B_poles.level_buckets) {
        if (buckets_at_a_level.empty()) { level++; continue; }
        double sum_m2l = 0., sum_near = 0., max_offset_ratio = 0.;
        double level_dL = buckets_at_a_level[0]->dL;
        for (auto& A : buckets_at_a_level) {
          sum_m2l += A->buckets_for_M2L.size();
          sum_near += A->buckets_near.size();
          if (level_dL > 0)
            max_offset_ratio = std::max(max_offset_ratio, A->max_internal_offset / level_dL);
        }
        size_t n = buckets_at_a_level.size();
        std::cout << " [FMM:stat] "
                  << std::setw(5) << level
                  << " | " << std::setw(7) << n
                  << " | " << std::setw(9) << std::scientific << std::setprecision(2) << level_dL
                  << " | " << std::setw(8) << std::fixed << std::setprecision(1) << sum_m2l / n
                  << " | " << std::setw(9) << std::fixed << std::setprecision(1) << sum_near / n
                  << " | " << std::setw(13) << std::fixed << std::setprecision(3) << max_offset_ratio
                  << std::endl;
        level++;
      }
    }
  } else {
    // ------------------------------------------------------------------
    // Step 3: ツリー構造が変わった → 全て再構築
    //   M2M 準備 (setM2M) → 親方向への集約で必要な幾何キャッシュ生成
    //   M2L 準備 (setM2L) → 相互作用(遠方)リストを構築し M2L キャッシュ生成
    // ------------------------------------------------------------------
    setM2M(B_poles);
    log_stage("Step3: setM2M done");
    setM2L(B_poles);
    log_stage("Step3: setM2L done");
#if defined(USE_METAL_M2L)
    // ------------------------------------------------------------------
    // Step 3.5: Metal M2L GPU acceleration setup (optional)
    //  - Initialize Metal M2L wrapper if use_metal_m2l=true in settings.json
    //  - Setup GPU buffers with M2L data structure
    // ------------------------------------------------------------------
    initMetalM2L(metal_m2l_settings.enabled, metal_m2l_settings.threadgroup, metal_m2l_settings.sort_terms);
    if (g_metal_m2l_wrapper) {
      setupMetalM2L(B_poles);
      log_stage("Step3.5: setupMetalM2L done");
    }
#endif

    {
      std::cout << Magenta << " [FMM:stat] level | buckets |    dL     | mean_M2L | mean_near | max_offset/dL" << colorReset << std::endl;
      int level = 0;
      for (auto& buckets_at_a_level : B_poles.level_buckets) {
        if (buckets_at_a_level.empty()) { level++; continue; }
        double sum_m2l = 0., sum_near = 0., max_offset_ratio = 0.;
        double level_dL = buckets_at_a_level[0]->dL;
        for (auto& A : buckets_at_a_level) {
          sum_m2l += A->buckets_for_M2L.size();
          sum_near += A->buckets_near.size();
          if (level_dL > 0)
            max_offset_ratio = std::max(max_offset_ratio, A->max_internal_offset / level_dL);
        }
        size_t n = buckets_at_a_level.size();
        std::cout << " [FMM:stat] "
                  << std::setw(5) << level
                  << " | " << std::setw(7) << n
                  << " | " << std::setw(9) << std::scientific << std::setprecision(2) << level_dL
                  << " | " << std::setw(8) << std::fixed << std::setprecision(1) << sum_m2l / n
                  << " | " << std::setw(9) << std::fixed << std::setprecision(1) << sum_near / n
                  << " | " << std::setw(13) << std::fixed << std::setprecision(3) << max_offset_ratio
                  << std::endl;
        level++;
      }
    }

    // ------------------------------------------------------------------
    // Step 4: L2L 準備 (setL2L)
    //  - 親 Local Moment を子へ伝搬するための幾何キャッシュ (spherical harmonics)を設定
    // ------------------------------------------------------------------
    setL2L(B_poles);
    log_stage("Step4: setL2L done");
  }

  // ------------------------------------------------------------------
  // Step 5: L2P 用ターゲット前処理
  //  - 各ターゲット点 t に対し，必要な (Y_nm ρ^{n+1}) * M_local のペアを前計算
  //  - integrateFarField() が O(N^2) 反復で足し合わせるだけになる
  // ------------------------------------------------------------------
  TimeWatch tw;
  // std::cout << "各レベルの各セルのL2Pの相手を保存する (#targets=" << targets.size() << ")" << std::endl;

  _Pragma("omp parallel for") for (auto& t : targets)
      t->setL2P(B_poles); // !L2Pの初期化
  log_stage("Step5: setL2P done");

  // ------------------------------------------------------------------
  // Step 6: 直接積分 (near field) 用前処理
  //  - 近傍バケツ (buckets_near) とそのソースを列挙し，ターゲットごとにソース寄与をキャッシュ (near_indices/near_weights_*)
  //  - GMRES の反復で再構築不要としコスト削減
  //  - nearfield_mode: "scalar" (default), "simd" (NEON 4-target), "metal" (GPU, future)
  // ------------------------------------------------------------------
  {
    std::stringstream ss_begin;
    ss_begin << "Step6: " << "setDirectIntegration (" << nearfield_mode << ") begin. "
             << "Total target size = " << targets.size();
    log_stage(ss_begin.str());
  }
#if defined(USE_METAL_NEARFIELD)
  if (nearfield_mode == "metal") {
    setDirectIntegrationMetal_linear(B_poles, targets);
  } else
#endif
#if defined(__aarch64__)
  if (nearfield_mode == "simd") {
    setDirectIntegrationSIMD_linear(B_poles, targets);
  } else if (nearfield_mode == "simd_double") {
    setDirectIntegrationSIMD_double_linear(B_poles, targets, g_p2m_quadrature_points);
  } else if (nearfield_mode == "cell_scalar") {
    setDirectIntegrationCellScalar_linear(B_poles, targets, g_p2m_quadrature_points);
  } else if (nearfield_mode == "flat_scalar") {
    setDirectIntegrationFlatScalar_linear(B_poles, targets, g_p2m_quadrature_points);
  } else
#endif
  {
    _Pragma("omp parallel for") for (auto& t : targets) {
      t->setDirectIntegration(B_poles);
    }
  }
  double cell_count = 0.0;
  for (auto& t : targets)
    cell_count += t->near_cell_count;
  std::stringstream ss;
  ss << "Step6: " << "setDirectIntegration (" << nearfield_mode << ")" << " done. "
     << "Total target size = " << targets.size()
     << ". Mean near cells = " << Red << (cell_count / targets.size()) << colorReset;
  log_stage(ss.str());

  // BEM行列のフィンガープリント: max_level間での一致を検証
  {
    std::size_t total_nnz = 0;
    double sum_wg = 0.0, sum_wgn = 0.0;
    double max_wg = 0.0, max_wgn = 0.0;
    for (const auto& t : targets) {
      total_nnz += t->near_indices.size();
      for (std::size_t i = 0; i < t->near_weights_phi.size(); ++i) {
        sum_wg += t->near_weights_phi[i];
        sum_wgn += t->near_weights_phin[i];
        max_wg = std::max(max_wg, std::abs(t->near_weights_phi[i]));
        max_wgn = std::max(max_wgn, std::abs(t->near_weights_phin[i]));
      }
    }
    std::cout << " [FMM:matrix] total_nnz=" << total_nnz
              << " mean_nnz/target=" << std::fixed << std::setprecision(1) << (double)total_nnz / targets.size()
              << " sum_wG=" << std::scientific << std::setprecision(10) << sum_wg
              << " sum_wGn=" << sum_wgn
              << " max_wG=" << max_wg
              << " max_wGn=" << max_wgn
              << std::defaultfloat << std::endl;
  }
};

inline void initializeFMM(auto& B_poles, const auto& targets) { initializeFMM(B_poles, targets, false); }

/* -------------------------------------------------------------------------- */

inline void updateFMM(auto& B_poles, std::array<double, 6>& elapsed_time) {
  TimeWatch tw;
  for (auto source : B_poles.data1D_vector)
    source->updateDensity();

  elapsed_time[0] = tw()[0];

  // !この操作は，Momentsのvector_weightedSourceDensities_RRn_MM_を削除しない．
  // !vector_weightedSourceDensities_RRn_MM_の保存するポインタ変数を初期化するだけで，配列内のポインタは生きたまま残る
  // !increment_Mを使うと，vector_weightedSourceDensities_RRn_MM_は初期化され再計算される
  B_poles.forEachAllParallel([&](auto* b) {
    b->MomentsMultipoleExpansion.initialize();
    b->MomentsLocalExpansion.initialize(); });
  elapsed_time[1] = tw()[0];

  B_poles.forEachAtDeepestParallel([&](auto* b) { b->MomentsMultipoleExpansion.increment_M_reuse(); });

  elapsed_time[2] = tw()[0];
  M2M(B_poles);

  elapsed_time[3] = tw()[0];
  M2L(B_poles);

  elapsed_time[4] = tw()[0];
  L2L(B_poles);

  elapsed_time[5] = tw()[0];
};

inline void updateFMM(auto& B_poles) {
  TimeWatch tw;
  for (auto source : B_poles.data1D_vector)
    source->updateDensity();

  std::cout << Yellow << "update poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

  // !この操作は，Momentsのvector_weightedSourceDensities_RRn_MM_を削除しない．
  // !vector_weightedSourceDensities_RRn_MM_の保存するポインタ変数を初期化するだけで，配列内のポインタは生きたまま残る
  // !increment_Mを使うと，vector_weightedSourceDensities_RRn_MM_は初期化され再計算される
  B_poles.forEachAllParallel([&](auto* b) {
    b->MomentsMultipoleExpansion.initialize();
    b->MomentsLocalExpansion.initialize(); });

  std::cout << Yellow << "reset Moments" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

  B_poles.forEachAtDeepestParallel([&](auto* b) { b->MomentsMultipoleExpansion.increment_M_reuse(); });
  std::cout << Yellow << "increment_M_reuse" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

  M2M(B_poles);
  std::cout << Magenta << "M2M" << ", Elapsed time : " << tw() << colorReset << std::endl;

  M2L(B_poles);
  std::cout << Magenta << "M2L" << ", Elapsed time : " << tw() << colorReset << std::endl;

  L2L(B_poles);
  std::cout << Magenta << "L2L" << ", Elapsed time : " << tw() << colorReset << std::endl;

  std::cout << Yellow << "updateFMM" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
};

/* -------------------------------------------------------------------------- */

inline void write2Dcsv(const std::string& filename, const std::vector<std::vector<double>>& data) {
  std::ofstream ofs(filename);
  for (const auto& row : data) {
    for (size_t j = 0; j < row.size(); ++j) {
      ofs << std::setprecision(15) << row[j];
      if (j < row.size() - 1)
        ofs << ",";
    }
    ofs << std::endl;
  }
  std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
  ofs.close();
};

inline void write2Dcsv(const std::string& filename, const std::vector<std::vector<std::complex<double>>>& data) {
  std::ofstream ofs(filename);
  for (const auto& row : data) {
    for (size_t j = 0; j < row.size(); ++j) {
      ofs << std::setprecision(15) << row[j].real();
      if (j < row.size() - 1)
        ofs << ",";
    }
    ofs << std::endl;
  }
  std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
  ofs.close();
};

static inline std::vector<std::vector<double>> Re(const auto& data) {
  std::vector<std::vector<double>> result(data.size(), std::vector<double>(data[0].size()));
  for (size_t i = 0; i < data.size(); ++i)
    for (size_t j = 0; j < data[0].size(); ++j)
      result[i][j] = data[i][j].real();
  return result;
}
inline std::vector<std::vector<double>> Im(const auto& data) {
  std::vector<std::vector<double>> result(data.size(), std::vector<double>(data[0].size()));
  for (size_t i = 0; i < data.size(); ++i)
    for (size_t j = 0; j < data[0].size(); ++j)
      result[i][j] = data[i][j].imag();
  return result;
}
