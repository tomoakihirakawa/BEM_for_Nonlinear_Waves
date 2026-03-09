#ifndef lib_plane_wave_expansion_H
#define lib_plane_wave_expansion_H

/*
 * Plane-wave (Exponential) Expansion for Fast Multipole Method
 * Based on Greengard & Rokhlin (1997) "A New Version of the Fast Multipole Method"
 *
 * SIMPLIFIED VERSION: Uses the same m2lFunction as SimpleM2L to ensure correctness.
 * This version directly computes the M2L using the existing SphericalCoordinates
 * infrastructure, which is known to work correctly.
 *
 * NOTE: A true O(p²) exponential expansion implementation for Laplace equation
 * requires the Sommerfeld identity and careful treatment of the 6-direction
 * decomposition. The formulas are significantly different from Helmholtz case.
 *
 * For now, this implementation provides the same interface as the future
 * optimized version, but uses O(p⁴) direct M2L internally.
 */

#include <array>
#include <cmath>
#include <complex>
#include <mutex>
#include <numbers>
#include <tuple>
#include <vector>

#include "lib_multipole_expansion.hpp"

namespace PlaneWave {

using cmplx = std::complex<double>;
constexpr double PI = std::numbers::pi;

/*
 * Direction indices for 6-direction decomposition (for future optimization)
 */
enum Direction { Up = 0, Down = 1, North = 2, South = 3, East = 4, West = 5 };

/*
 * PlaneWaveExpansion class - stores accumulated local expansion during M2L
 *
 * In this simplified version, we bypass the full plane-wave machinery and
 * directly accumulate the local expansion coefficients.
 */
template <int N, int N_k = 16, int N_alpha = 2 * N + 2>
struct PlaneWaveExpansion {
    static constexpr int num_k = N_k;
    static constexpr int num_alpha = N_alpha;
    static constexpr int num_points = num_k * num_alpha;

    bool initialized = false;

    // Direct storage for accumulated local expansion during M2L
    static constexpr int MM_SIZE = (N + 1) * (N + 1);
    std::array<std::tuple<cmplx, cmplx>, MM_SIZE> accumulated_L;

    PlaneWaveExpansion() {
        clear();
    }

    void clear() {
        initialized = false;
        for (auto& elem : accumulated_L) {
            std::get<0>(elem) = cmplx(0.0);
            std::get<1>(elem) = cmplx(0.0);
        }
    }
};

/*
 * PlaneWaveConverter - placeholder for conversion matrices (not used in simplified version)
 */
template <int N, int N_k = 16, int N_alpha = 2 * N + 2>
class PlaneWaveConverter {
public:
    using PWE = PlaneWaveExpansion<N, N_k, N_alpha>;
    static constexpr int MM_SIZE = (N + 1) * (N + 1);

    bool initialized = false;

    static int linear_nm(int n, int m) {
        return n * n + n + m;
    }

    PlaneWaveConverter() {}

    void initialize() {
        if (initialized) return;
        initialized = true;
    }
};

/*
 * M2L operation descriptor for the quadruple sum
 */
struct M2LOp {
    int j, k, n, m;  // indices
    int src;         // source index = n*n + n + m
    int dst;         // destination index = j*j + j + k
};

/*
 * Generate M2L operations for order N
 * Thread-safe initialization using lambda pattern (C++11 magic statics)
 */
template <int N>
inline const std::vector<M2LOp>& m2l_ops_planewave() {
    static const std::vector<M2LOp> ops = []() {
        std::vector<M2LOp> result;
        for (int j = 0; j <= N; ++j) {
            for (int k = -j; k <= j; ++k) {
                for (int n = 0; n <= N; ++n) {
                    for (int m = -n; m <= n; ++m) {
                        int p = j + n;
                        int q = m - k;
                        if (p <= 2 * N && std::abs(q) <= p) {
                            int src = n * n + n + m;
                            int dst = j * j + j + k;
                            result.push_back({j, k, n, m, src, dst});
                        }
                    }
                }
            }
        }
        return result;
    }();
    return ops;
}

/*
 * PlaneWaveM2LOperator: Main interface for plane-wave based M2L
 *
 * SIMPLIFIED VERSION: Uses the existing SphericalCoordinates::m2lFunction
 * to compute M2L translations. This ensures correctness by using the same
 * formula as SimpleM2L.
 */
template <int N, int N_k = 16, int N_alpha = 2 * N + 2>
class PlaneWaveM2LOperator {
public:
    using PWE = PlaneWaveExpansion<N, N_k, N_alpha>;
    using Converter = PlaneWaveConverter<N, N_k, N_alpha>;
    static constexpr int MM_SIZE = (N + 1) * (N + 1);

    // Thread-safe converter initialization using call_once
    static Converter& getConverter() {
        static Converter converter;
        static std::once_flag init_flag;
        std::call_once(init_flag, []() {
            converter.initialize();
        });
        return converter;
    }

    static void initializeConverter() {
        getConverter(); // Triggers thread-safe initialization
    }

    /*
     * Convert multipole expansion to intermediate representation
     * In this simplified version, this is a no-op - we just initialize the target.
     */
    template <typename MM_Type>
    static void convertMultipoleToExponential(
        const MM_Type& MM,
        PWE& exp,
        bool realfield_m_conj = false
    ) {
        initializeConverter();
        exp.clear();
        exp.initialized = true;
    }

    /*
     * Translate and accumulate M2L contribution from source multipole
     * to target local expansion.
     *
     * Uses the existing SphericalCoordinates::m2lFunction for correctness.
     * Note: For real-field, the symmetry L_{j,-k} = (-1)^k * conj(L_{j,k})
     * is enforced by the calling code (enforce_conjugate_m_symmetry_if_enabled).
     */
    template <typename MM_Type>
    static void translateExponentialWithSource(
        const MM_Type& source_MM,
        PWE& target_exp,
        const std::array<double, 3>& translation,
        bool realfield_m_conj = false
    ) {
        initializeConverter();

        // Create SphericalCoordinates for the translation vector
        // and precompute spherical harmonics
        Tddd trans_vec = {translation[0], translation[1], translation[2]};
        SphericalCoordinates sph(trans_vec);
        sph.precompute_sph_div_rhon1(2 * N);

        // Compute M2L using the same formula as SimpleM2L
        // Note: m2l_ops includes all (j,k,n,m) combinations
        // For real-field, we skip k < 0 (will be reconstructed by symmetry)
        // but we process all m values (source MM has all m values)
        for (const auto& op : m2l_ops_planewave<N>()) {
            if (realfield_m_conj && op.k < 0)
                continue;

            // Use the existing m2lFunction which includes AAA_M2L_FMM coefficients
            cmplx AAAY = sph.m2lFunction(op.j, op.k, op.n, op.m);

            // Get source multipole coefficient
            const cmplx& M0 = std::get<0>(source_MM[op.src]);
            const cmplx& M1 = std::get<1>(source_MM[op.src]);

            // Accumulate into target local expansion
            std::get<0>(target_exp.accumulated_L[op.dst]) += AAAY * M0;
            std::get<1>(target_exp.accumulated_L[op.dst]) += AAAY * M1;
        }
    }

    /*
     * Wrapper for translateExponential (no-op in simplified version)
     */
    static void translateExponential(
        const PWE& source_exp,
        PWE& target_exp,
        const std::array<double, 3>& translation,
        bool realfield_m_conj = false
    ) {
        // No-op: actual translation is done in translateExponentialWithSource
    }

    /*
     * Convert accumulated expansion to local expansion
     * Simply copies the accumulated values.
     */
    template <typename L_Type>
    static void convertExponentialToLocal(
        const PWE& exp,
        L_Type& L,
        bool realfield_m_conj = false
    ) {
        initializeConverter();

        // Copy accumulated local expansion to output
        for (int j = 0; j <= N; ++j) {
            for (int k = -j; k <= j; ++k) {
                int jk_idx = j * j + j + k;
                std::get<0>(L[jk_idx]) += std::get<0>(exp.accumulated_L[jk_idx]);
                std::get<1>(L[jk_idx]) += std::get<1>(exp.accumulated_L[jk_idx]);
            }
        }
    }
};

// Note: Static converter is now managed by getConverter() using std::call_once for thread safety

} // namespace PlaneWave

#endif // lib_plane_wave_expansion_H
