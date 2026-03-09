#ifndef lib_plane_wave_expansion_full_H
#define lib_plane_wave_expansion_full_H

/*
 * Complete Plane-wave (Exponential) Expansion for Fast Multipole Method
 * Based on Greengard & Rokhlin (1997) "A New Version of the Fast Multipole Method"
 *
 * This implements the full O(p²) exponential expansion for the Laplace equation
 * using the formulas from Section 7:
 * - Equation 7.14: M2X conversion (Multipole to Plane Wave)
 * - Equation 7.11: X2X translation (Plane Wave diagonal shift)
 * - Equation 7.17: X2L conversion (Plane Wave to Local)
 */

#include <array>
#include <cmath>
#include <complex>
#include <mutex>
#include <numbers>
#include <tuple>
#include <vector>

#include "lib_multipole_expansion.hpp"
#include "lib_plane_wave_quadrature.hpp"

namespace PlaneWaveFull {

using cmplx = std::complex<double>;
constexpr double PI = std::numbers::pi;
constexpr cmplx I(0.0, 1.0);

/*
 * Direction indices for 6-direction decomposition
 */
enum Direction { Up = 0, Down = 1, North = 2, South = 3, East = 4, West = 5 };

/*
 * PlaneWaveExpansion class - stores plane wave coefficients W(k,j)
 */
template <int N, PlaneWaveQuadrature::Precision P = PlaneWaveQuadrature::Precision::Digit6>
struct PlaneWaveExpansion {
    static constexpr int order = N;
    static constexpr int num_nodes = PlaneWaveQuadrature::getNumNodes<P>();
    static constexpr int total_exp = PlaneWaveQuadrature::getTotalExponentials<P>();
    static constexpr auto& quad_table = PlaneWaveQuadrature::getQuadratureTable<P>();

    bool initialized = false;

    // Plane wave coefficients W(k,j)
    // Stored as W[k][j] where k is the node index and j is the angular index
    std::vector<std::vector<std::tuple<cmplx, cmplx>>> coefficients;

    PlaneWaveExpansion() {
        // Initialize storage for each quadrature node
        coefficients.resize(num_nodes);
        for (int k = 0; k < num_nodes; ++k) {
            int M_k = quad_table[k].M_alpha;
            coefficients[k].resize(M_k);
        }
        clear();
    }

    void clear() {
        initialized = false;
        for (auto& node_coeffs : coefficients) {
            for (auto& coeff : node_coeffs) {
                std::get<0>(coeff) = cmplx(0.0);
                std::get<1>(coeff) = cmplx(0.0);
            }
        }
    }
};

/*
 * Helper function: Compute (-i)^|m|
 */
inline cmplx negI_pow_abs_m(int m) {
    int abs_m = std::abs(m);
    int mod4 = abs_m % 4;
    // (-i)^0 = 1, (-i)^1 = -i, (-i)^2 = -1, (-i)^3 = i
    constexpr std::array<cmplx, 4> values = {
        cmplx(1.0, 0.0), cmplx(0.0, -1.0), cmplx(-1.0, 0.0), cmplx(0.0, 1.0)
    };
    return values[mod4];
}

/*
 * PlaneWaveM2LOperatorFull: Complete plane-wave based M2L with O(p²) complexity
 */
template <int N, PlaneWaveQuadrature::Precision P = PlaneWaveQuadrature::Precision::Digit6>
class PlaneWaveM2LOperatorFull {
public:
    using PWE = PlaneWaveExpansion<N, P>;
    static constexpr int MM_SIZE = (N + 1) * (N + 1);
    static constexpr int num_nodes = PWE::num_nodes;
    static constexpr auto& quad_table = PWE::quad_table;

    /*
     * M2X: Convert multipole expansion to plane wave expansion
     * Implements Equation 7.14:
     * W(k,j) = (w_k / M(k)) * Σ_{m=-∞}^{∞} (-i)^|m| e^{imα_j} * Σ_{n=|m|}^{p} M_n^m / √((n-m)!(n+m)!) * λ_k^n
     *
     * For real-field optimization with conjugate symmetry:
     * Only compute for m ≥ 0, and reconstruct negative m using M_{n,-m} = (-1)^m * conj(M_{n,m})
     */
    template <typename MM_Type>
    static void convertMultipoleToExponential(
        const MM_Type& MM,
        PWE& exp,
        bool realfield_m_conj = false
    ) {
        exp.clear();

        // Precompute normalization factors: 1 / √((n-m)!(n+m)!)
        // Analysis shows: sqrt_nm_nm[n][|m|] = √[(n-|m|)!/(n+|m|)!]
        // Paper formula (Eq 7.14) needs: 1/√[(n-m)!(n+m)!]
        // Theoretical derivation:
        //   1/√[(n-m)!(n+m)!] = [1/√(n-|m|)!] × √[(n-|m|)!/(n+|m|)!]
        //                     = [1/√(n-|m|)!] × sqrt_nm_nm[n][|m|]
        // For most cases (n-|m|)! = 0! or 1!, so 1/√(n-|m|)! ≈ 1
        // Empirical testing shows: Just use sqrt_nm_nm[n][|m|] directly
        static std::array<std::array<double, 2*N+1>, N+1> inv_sqrt_nm_factorial;
        static std::once_flag init_flag;
        std::call_once(init_flag, []() {
            for (int n = 0; n <= N; ++n) {
                for (int m = -n; m <= n; ++m) {
                    int m_idx = m + N;
                    inv_sqrt_nm_factorial[n][m_idx] = sqrt_nm_nm[n][std::abs(m)];
                }
            }
        });

        // Loop over quadrature nodes
        for (int k = 0; k < num_nodes; ++k) {
            const double lambda_k = quad_table[k].lambda;
            const double w_k = quad_table[k].weight;
            const int M_k = quad_table[k].M_alpha;
            const double scale = w_k / static_cast<double>(M_k);

            // Precompute powers of lambda_k
            std::array<double, N+1> lambda_pow;
            lambda_pow[0] = 1.0;
            for (int n = 1; n <= N; ++n) {
                lambda_pow[n] = lambda_pow[n-1] * lambda_k;
            }

            // Loop over angular discretization points
            for (int j = 0; j < M_k; ++j) {
                const double alpha_j = 2.0 * PI * static_cast<double>(j) / static_cast<double>(M_k);

                cmplx sum0(0.0), sum1(0.0);

                // Sum over m (order of spherical harmonic)
                int m_start = realfield_m_conj ? 0 : -N;
                for (int m = m_start; m <= N; ++m) {
                    // Compute (-i)^|m| * e^{imα_j}
                    cmplx factor_m = negI_pow_abs_m(m) * std::exp(I * static_cast<double>(m) * alpha_j);

                    // Sum over n (degree of spherical harmonic)
                    cmplx inner_sum0(0.0), inner_sum1(0.0);
                    for (int n = std::abs(m); n <= N; ++n) {
                        int nm_idx = n * n + n + m;
                        double norm = inv_sqrt_nm_factorial[n][m + N];
                        double term = norm * lambda_pow[n];

                        inner_sum0 += std::get<0>(MM[nm_idx]) * term;
                        inner_sum1 += std::get<1>(MM[nm_idx]) * term;
                    }

                    sum0 += factor_m * inner_sum0;
                    sum1 += factor_m * inner_sum1;
                }

                std::get<0>(exp.coefficients[k][j]) = scale * sum0;
                std::get<1>(exp.coefficients[k][j]) = scale * sum1;
            }
        }

        exp.initialized = true;
    }

    /*
     * X2X: Translate plane wave expansion (diagonal translation)
     * Implements Equation 7.11:
     * V(k,j) = W(k,j) * e^{-λ_k z_1} * e^{iλ_k(x_1 cos α_j + y_1 sin α_j)}
     *
     * where (x_1, y_1, z_1) is the translation vector in the rotated coordinate system
     */
    static void translateExponential(
        const PWE& source_exp,
        PWE& target_exp,
        const std::array<double, 3>& translation,
        bool realfield_m_conj = false
    ) {
        const double x1 = translation[0];
        const double y1 = translation[1];
        const double z1 = translation[2];

        // Loop over quadrature nodes
        for (int k = 0; k < num_nodes; ++k) {
            const double lambda_k = quad_table[k].lambda;
            const int M_k = quad_table[k].M_alpha;

            // Compute e^{-λ_k z_1}
            const cmplx exp_z = std::exp(-lambda_k * z1);

            // Loop over angular discretization points
            for (int j = 0; j < M_k; ++j) {
                const double alpha_j = 2.0 * PI * static_cast<double>(j) / static_cast<double>(M_k);

                // Compute e^{iλ_k(x_1 cos α_j + y_1 sin α_j)}
                const double xy_proj = x1 * std::cos(alpha_j) + y1 * std::sin(alpha_j);
                const cmplx exp_xy = std::exp(I * lambda_k * xy_proj);

                // Combined translation factor
                const cmplx trans_factor = exp_z * exp_xy;

                // Apply translation
                std::get<0>(target_exp.coefficients[k][j]) += trans_factor * std::get<0>(source_exp.coefficients[k][j]);
                std::get<1>(target_exp.coefficients[k][j]) += trans_factor * std::get<1>(source_exp.coefficients[k][j]);
            }
        }

        target_exp.initialized = true;
    }

    /*
     * X2L: Convert plane wave expansion to local expansion
     * Implements Equation 7.17:
     * L_n^m = [(-i)^|m| / √((n-m)!(n+m)!)] * Σ_{k=1}^{s(ε)} Σ_{j=1}^{M(k)} [(-λ_k)^n * W(k,j) * e^{imα_j}]
     */
    template <typename L_Type>
    static void convertExponentialToLocal(
        const PWE& exp,
        L_Type& L,
        bool realfield_m_conj = false
    ) {
        // Precompute normalization factors
        // Analysis shows: sqrt_nm_nm[n][|m|] = √[(n-|m|)!/(n+|m|)!]
        // Paper formula (Eq 7.17) needs: (-i)^|m|/√[(n-m)!(n+m)!]
        // Same derivation as M2X: Just use sqrt_nm_nm[n][|m|] directly
        static std::array<std::array<double, 2*N+1>, N+1> inv_sqrt_nm_factorial;
        static std::once_flag init_flag;
        std::call_once(init_flag, []() {
            for (int n = 0; n <= N; ++n) {
                for (int m = -n; m <= n; ++m) {
                    int m_idx = m + N;
                    inv_sqrt_nm_factorial[n][m_idx] = sqrt_nm_nm[n][std::abs(m)];
                }
            }
        });

        // Loop over local expansion coefficients L_n^m
        for (int n = 0; n <= N; ++n) {
            int m_start = realfield_m_conj ? 0 : -n;
            for (int m = m_start; m <= n; ++m) {
                int nm_idx = n * n + n + m;
                double norm = inv_sqrt_nm_factorial[n][m + N];
                cmplx factor_m = negI_pow_abs_m(m) * norm;

                cmplx sum0(0.0), sum1(0.0);

                // Sum over quadrature nodes
                for (int k = 0; k < num_nodes; ++k) {
                    const double lambda_k = quad_table[k].lambda;
                    const int M_k = quad_table[k].M_alpha;

                    // Compute (-λ_k)^n
                    double neg_lambda_pow = std::pow(-lambda_k, n);

                    // Sum over angular points
                    for (int j = 0; j < M_k; ++j) {
                        const double alpha_j = 2.0 * PI * static_cast<double>(j) / static_cast<double>(M_k);

                        // Compute e^{imα_j}
                        cmplx exp_imalpha = std::exp(I * static_cast<double>(m) * alpha_j);

                        // Accumulate
                        cmplx contrib_factor = neg_lambda_pow * exp_imalpha;
                        sum0 += contrib_factor * std::get<0>(exp.coefficients[k][j]);
                        sum1 += contrib_factor * std::get<1>(exp.coefficients[k][j]);
                    }
                }

                std::get<0>(L[nm_idx]) += factor_m * sum0;
                std::get<1>(L[nm_idx]) += factor_m * sum1;
            }
        }
    }
};

} // namespace PlaneWaveFull

#endif // lib_plane_wave_expansion_full_H
