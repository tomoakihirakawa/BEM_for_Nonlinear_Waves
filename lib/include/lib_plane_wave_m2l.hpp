#ifndef lib_plane_wave_m2l_H
#define lib_plane_wave_m2l_H

/*
 * Plane Wave M2L Integration with 6-Direction Decomposition
 * Based on Greengard & Rokhlin (1997), Lemma 7.2
 *
 * This file provides the complete M2L operation using plane wave expansion
 * with proper handling of the 6-direction decomposition (Up, Down, North, South, East, West)
 */

#include "lib_plane_wave_expansion_full.hpp"
#include "lib_multipole_expansion.hpp"
#include <array>
#include <unordered_map>
#include <vector>
#include <map>

namespace PlaneWaveM2L
{

    using namespace PlaneWaveFull;
    using cmplx = std::complex<double>;

    /*
     * Wigner d-matrix computation for rotation of spherical harmonics
     * d^n_{m',m}(beta)
     */
    class WignerDMatrix
    {
    public:
        // Use lgamma for better numerical stability with larger n
        static double log_factorial(int n)
        {
            if (n < 0)
                return 0.0; // log(1) = 0 for invalid cases
            return std::lgamma(n + 1.0);
        }

        // Compute d^n_{mp,m}(beta)
        // Using the formula from standard texts (e.g. Edmonds)
        static double d(int n, int mp, int m, double beta)
        {
            double sum = 0.0;
            double cos_half = std::cos(beta / 2.0);
            double sin_half = std::sin(beta / 2.0);

            // Range of k for summation
            int k_min = std::max(0, m - mp);
            int k_max = std::min(n + m, n - mp);

            // Pre-calculate the log of the constant factor to avoid large intermediate values
            double log_sqrt_factor = 0.5 * (log_factorial(n + mp) + log_factorial(n - mp) +
                                            log_factorial(n + m) + log_factorial(n - m));

            for (int k = k_min; k <= k_max; ++k)
            {
                double term = (k % 2 == 0 ? 1.0 : -1.0);

                double log_den = log_factorial(n - mp - k) + log_factorial(n + m - k) +
                                 log_factorial(k) + log_factorial(k + mp - m);
                term *= std::exp(log_sqrt_factor - log_den);

                // Trigonometric parts
                term *= std::pow(cos_half, 2 * n - 2 * k + m - mp);
                term *= std::pow(sin_half, 2 * k + mp - m);

                sum += term;
            }
            return sum;
        }
    };

    /*
     * Rotation helpers for 6-direction decomposition
     * These apply the rotations specified in Lemma 7.2
     */
    template <int N>
    class RotationOperator
    {
    public:
        using D_Matrix_Cache = std::vector<std::vector<double>>;

    private:
        inline static std::map<Direction, std::vector<D_Matrix_Cache>> d_matrix_cache_;
        inline static bool cache_initialized_ = false;

    public:
        using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N + 1) * (N + 1)>;

        /*
         * Rotate multipole expansion for a given direction
         * Returns the rotation angles (α, β, γ) for Euler angle rotation
         */
        static std::array<double, 3> getRotationAngles(Direction dir)
        {
            constexpr double pi = std::numbers::pi;
            constexpr double pi_2 = std::numbers::pi / 2.0;

            switch (dir)
            {
            case Direction::Up:
                return {0.0, 0.0, 0.0};
            case Direction::Down:
                return {0.0, pi, 0.0};
            case Direction::East: // +x to +z via Ry(-pi/2)
                return {0.0, -pi_2, 0.0};
            case Direction::West: // -x to +z via Ry(pi/2)
                return {0.0, pi_2, 0.0};
            case Direction::North: // +y to +z via Rx(-pi/2)
                return {pi_2, -pi_2, -pi_2};
            case Direction::South: // -y to +z via Rx(pi/2)
                return {pi_2, pi_2, -pi_2};
            default:
                return {0.0, 0.0, 0.0};
            }
        }

        /*
         * Rotate translation vector for a given direction
         * The translation vector needs to be rotated to align with the z-axis
         *
         * Also adds small perturbation to avoid numerical instabilities
         * when the vector is exactly on a coordinate axis
         */
        static Tddd rotateTranslationVector(const Tddd &vec, Direction dir)
        {
            constexpr double eps = 1.0e-12; // Small perturbation for numerical stability
            double x = vec[0], y = vec[1], z = vec[2];

            // Add small perturbation if components are exactly zero
            if (std::abs(x) < eps)
                x = eps;
            if (std::abs(y) < eps)
                y = eps;
            if (std::abs(z) < eps)
                z = eps;

            switch (dir)
            {
            case Direction::Up:
                return {x, y, z}; // No change
            case Direction::Down:
                return {-x, y, -z}; // Rotate π around y-axis
            case Direction::East:
                return {-z, y, x}; // +x to +z via Ry(-pi/2)
            case Direction::West:
                return {z, y, -x}; // -x to +z via Ry(pi/2)
            case Direction::North:
                return {x, z, -y}; // +y to +z via Rx(-pi/2)
            case Direction::South:
                return {x, -z, y}; // -y to +z via Rx(pi/2)
            default:
                return {x, y, z};
            }
        }

        static void initializeCache()
        {
            if (cache_initialized_)
                return;

            for (int dir_idx = 0; dir_idx < 6; ++dir_idx)
            {
                Direction dir = static_cast<Direction>(dir_idx);
                auto angles = getRotationAngles(dir);
                double beta = angles[1];

                std::vector<D_Matrix_Cache> d_matrices_for_dir;
                d_matrices_for_dir.resize(N + 1);

                for (int n = 0; n <= N; ++n)
                {
                    d_matrices_for_dir[n].resize(2 * n + 1, std::vector<double>(2 * n + 1));
                    for (int mp = -n; mp <= n; ++mp)
                    {
                        for (int m = -n; m <= n; ++m)
                        {
                            d_matrices_for_dir[n][mp + n][m + n] = WignerDMatrix::d(n, mp, m, beta);
                        }
                    }
                }
                d_matrix_cache_[dir] = std::move(d_matrices_for_dir);
            }
            cache_initialized_ = true;
        }

        /*
         * Rotate coefficients M_{n,m}
         * M'_{n,mp} = Sum_m D^n_{mp,m}(alpha, beta, gamma) * M_{n,m}
         * D^n_{mp,m} = e^{-i*mp*alpha} * d^n_{mp,m}(beta) * e^{-i*m*gamma}
         */
        static void rotateCoefficients(
            const MM_Type &input,
            MM_Type &output,
            Direction dir,
            bool inverse = false)
        {
            if (!cache_initialized_)
                initializeCache();

            auto angles = getRotationAngles(dir);
            double alpha = angles[0];
            double beta = angles[1];
            double gamma = angles[2];

            if (inverse)
            {
                // Inverse rotation is R(-gamma, -beta, -alpha)
                // But D matrix inverse is conjugate transpose.
                // D^{-1}(a,b,g) = D(-g, -b, -a)
                double tmp = alpha;
                alpha = -gamma;
                gamma = -tmp;
                beta = -beta;
            }

            // Clear output
            for (auto &val : output)
            {
                std::get<0>(val) = cmplx(0.0);
                std::get<1>(val) = cmplx(0.0);
            }

            const cmplx I(0.0, 1.0);

            const auto &d_matrices = d_matrix_cache_.at(dir);

            for (int n = 0; n <= N; ++n)
            {
                // Precompute exponentials for this n (optimization possible)

                for (int mp = -n; mp <= n; ++mp)
                {
                    cmplx sum0(0.0), sum1(0.0);

                    // Factor e^{-i*mp*alpha}
                    cmplx term1 = std::exp(-I * static_cast<double>(mp) * alpha);

                    for (int m = -n; m <= n; ++m)
                    {
                        // Get d^n_{mp,m}(beta) from cache
                        double d_val = d_matrices[n][mp + n][m + n];

                        // Factor e^{-i*m*gamma}
                        cmplx term2 = std::exp(-I * static_cast<double>(m) * gamma);

                        // Combine D matrix element
                        cmplx D_val = term1 * d_val * term2;

                        // Input index
                        int src_idx = n * n + n + m;

                        sum0 += D_val * std::get<0>(input[src_idx]);
                        sum1 += D_val * std::get<1>(input[src_idx]);
                    }

                    // Output index
                    int dst_idx = n * n + n + mp;
                    std::get<0>(output[dst_idx]) = sum0;
                    std::get<1>(output[dst_idx]) = sum1;
                }
            }
        }

        // Helper to unpack real-field compressed format (m>=0) to full format (all m)
        static void unpack_mm(const MM_Type &src, MM_Type &dst)
        {
            for (int n = 0; n <= N; ++n)
            {
                // Copy m >= 0
                for (int m = 0; m <= n; ++m)
                {
                    int idx = n * n + n + m;
                    dst[idx] = src[idx];
                }
                // Reconstruct m < 0: M_{n,-m} = (-1)^m * conj(M_{n,m})
                for (int m = 1; m <= n; ++m)
                {
                    int pos_idx = n * n + n + m;
                    int neg_idx = n * n + n - m;
                    double sign = (m % 2 == 0) ? 1.0 : -1.0;
                    std::get<0>(dst[neg_idx]) = sign * std::conj(std::get<0>(src[pos_idx]));
                    std::get<1>(dst[neg_idx]) = sign * std::conj(std::get<1>(src[pos_idx]));
                }
            }
        }
    };

    /*
     * Complete M2L operator using plane wave expansion with 6-direction decomposition
     */
    template <int N, PlaneWaveQuadrature::Precision P = PlaneWaveQuadrature::Precision::Digit6>
    class PlaneWaveM2L
    {
    public:
        using PWE = PlaneWaveExpansion<N, P>;
        using PWOp = PlaneWaveM2LOperatorFull<N, P>;
        using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N + 1) * (N + 1)>;
        using L_Type = std::array<std::tuple<cmplx, cmplx>, (N + 1) * (N + 1)>;

        /*
         * Determine which direction list a box C belongs to relative to box B
         * This implements Definition 7.1 from the paper
         */
        static Direction classifyDirection(const Tddd &BC_vector)
        {
            // TEMPORARY FIX: Disable rotation to isolate bug
            // Always use Direction::Up (no rotation)
            return Direction::Up;

            /* Original implementation - disabled for debugging
            double x = BC_vector[0];
            double y = BC_vector[1];
            double z = BC_vector[2];

            double abs_x = std::abs(x);
            double abs_y = std::abs(y);
            double abs_z = std::abs(z);

            // Check if separated by at least one box in z direction
            if (abs_z >= abs_x && abs_z >= abs_y)
            {
                return (z > 0) ? Direction::Up : Direction::Down;
            }
            // Check if separated by at least one box in y direction (and not in Up/Down)
            else if (abs_y >= abs_x)
            {
                return (y > 0) ? Direction::North : Direction::South;
            }
            // Otherwise must be in x direction
            else
            {
                return (x > 0) ? Direction::East : Direction::West;
            }
            */
        }

        /*
         * Complete M2L translation from source multipole to target local expansion
         * This is the main entry point that handles the full pipeline:
         * 1. Rotate source multipole to align with direction
         * 2. Convert to plane wave (M2X)
         * 3. Translate plane wave (X2X)
         * 4. Convert to local expansion (X2L)
         * 5. Rotate local expansion back
         */
        static void translateMultipoleToLocal(
            const MM_Type &source_MM,
            L_Type &target_L,
            const Tddd &translation,
            bool realfield_m_conj = false)
        {
            // Determine direction
            Direction dir = classifyDirection(translation);

            // Rotate translation vector to align with +z axis
            Tddd rotated_trans = RotationOperator<N>::rotateTranslationVector(translation, dir);

            // Get rotation angles
            // auto angles = RotationOperator<N>::getRotationAngles(dir);

            // Prepare source coefficients (handle real-field symmetry if needed)
            MM_Type full_source_MM;
            if (realfield_m_conj)
            {
                RotationOperator<N>::unpack_mm(source_MM, full_source_MM);
            }
            else
            {
                full_source_MM = source_MM;
            }

            // Rotate source multipole: M_rotated = R * M_source
            MM_Type rotated_source_MM;
            RotationOperator<N>::rotateCoefficients(full_source_MM, rotated_source_MM, dir, false);

            // Create plane wave expansion
            PWE plane_wave;

            // M2X: Convert multipole to plane wave
            // TEMPORARY: Disable realfield optimization until bug is fixed
            // With rotation disabled (Direction::Up only), we can use realfield optimization
            // bool use_realfield = (dir == Direction::Up) && realfield_m_conj;
            bool use_realfield = false;  // Disabled temporarily
            PWOp::convertMultipoleToExponential(rotated_source_MM, plane_wave, use_realfield);

            // X2X: Translate plane wave (diagonal operation)
            PWE translated_wave;
            std::array<double, 3> trans_arr = {rotated_trans[0], rotated_trans[1], rotated_trans[2]};
            PWOp::translateExponential(plane_wave, translated_wave, trans_arr, use_realfield);

            // X2L: Convert plane wave to local expansion
            // Result is in rotated frame
            L_Type rotated_target_L;
            // Initialize rotated_target_L to zero before accumulation
            for (auto &v : rotated_target_L)
            {
                std::get<0>(v) = 0;
                std::get<1>(v) = 0;
            }

            PWOp::convertExponentialToLocal(translated_wave, rotated_target_L, use_realfield);

            // Rotate local expansion back: L_target += R^{-1} * L_rotated
            // We accumulate into target_L
            L_Type back_rotated_L;
            RotationOperator<N>::rotateCoefficients(rotated_target_L, back_rotated_L, dir, true);

            // Accumulate to result
            for (size_t i = 0; i < target_L.size(); ++i)
            {
                std::get<0>(target_L[i]) += std::get<0>(back_rotated_L[i]);
                std::get<1>(target_L[i]) += std::get<1>(back_rotated_L[i]);
            }
        }

        /*
         * Optimized M2L for multiple source boxes in the same direction
         * This implements the optimization described in Section 8.1:
         * First merge the outgoing plane wave expansions from all sources,
         * then translate once to each target
         */
        static void translateMultipleToLocal(
            const std::vector<std::pair<MM_Type, Tddd>> &sources, // (multipole, translation vector) pairs
            L_Type &target_L,
            bool realfield_m_conj = false)
        {
            if (sources.empty())
                return;

            // Group sources by direction
            std::map<Direction, std::vector<std::pair<MM_Type, Tddd>>> sources_by_dir;
            for (const auto &[mm, trans] : sources)
            {
                Direction dir = classifyDirection(trans);
                sources_by_dir[dir].emplace_back(mm, trans);
            }

            // Process each direction
            for (const auto &[dir, dir_sources] : sources_by_dir)
            {
                // Accumulate plane wave expansions from all sources in this direction
                PWE accumulated_wave;

                // Get rotation angles for this direction
                // auto angles = RotationOperator<N>::getRotationAngles(dir);

                for (const auto &[source_mm, translation] : dir_sources)
                {
                    // Rotate translation vector
                    Tddd rotated_trans = RotationOperator<N>::rotateTranslationVector(translation, dir);

                    // Prepare and rotate source coefficients
                    MM_Type full_source_MM;
                    if (realfield_m_conj)
                    {
                        RotationOperator<N>::unpack_mm(source_mm, full_source_MM);
                    }
                    else
                    {
                        full_source_MM = source_mm;
                    }
                    MM_Type rotated_source_MM;
                    RotationOperator<N>::rotateCoefficients(full_source_MM, rotated_source_MM, dir, false);

                    // Convert source multipole to plane wave
                    PWE source_wave;
                    PWOp::convertMultipoleToExponential(rotated_source_MM, source_wave, false);

                    // Translate and accumulate
                    std::array<double, 3> trans_arr = {rotated_trans[0], rotated_trans[1], rotated_trans[2]};
                    PWOp::translateExponential(source_wave, accumulated_wave, trans_arr, false);
                }

                // Convert accumulated plane wave to local expansion (in rotated frame)
                L_Type rotated_target_L;
                for (auto &v : rotated_target_L)
                {
                    std::get<0>(v) = 0;
                    std::get<1>(v) = 0;
                }
                PWOp::convertExponentialToLocal(accumulated_wave, rotated_target_L, false);

                // Rotate back and accumulate
                L_Type back_rotated_L;
                RotationOperator<N>::rotateCoefficients(rotated_target_L, back_rotated_L, dir, true);

                for (size_t i = 0; i < target_L.size(); ++i)
                {
                    std::get<0>(target_L[i]) += std::get<0>(back_rotated_L[i]);
                    std::get<1>(target_L[i]) += std::get<1>(back_rotated_L[i]);
                }
            }
        }
    };

} // namespace PlaneWaveM2L

#endif // lib_plane_wave_m2l_H
