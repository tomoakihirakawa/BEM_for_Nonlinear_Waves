/*
 * Test program for complete plane wave expansion implementation
 * Validates correctness by comparing against the simplified (known-correct) implementation
 * and measures performance improvement
 */

#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include "../../include/lib_plane_wave_expansion.hpp"
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;
constexpr int N = 6;  // Expansion order

/*
 * Create synthetic multipole expansion for testing
 * Uses simple random coefficients that are physically plausible
 */
template <int Order>
void createSyntheticMultipoleExpansion(
    std::array<std::tuple<cmplx, cmplx>, (Order+1)*(Order+1)>& MM,
    unsigned int seed = 42
) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    // Zero out
    for (auto& m : MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }

    // Create synthetic coefficients
    // Monopole (n=0) should dominate, higher orders smaller
    for (int n = 0; n <= Order; ++n) {
        // Decay factor: higher orders get smaller coefficients
        double decay = 1.0 / (1.0 + static_cast<double>(n));

        for (int m = -n; m <= n; ++m) {
            int idx = n * n + n + m;
            double real_part = dis(gen) * decay;
            double imag_part = dis(gen) * decay;
            std::get<0>(MM[idx]) = cmplx(real_part, imag_part);

            // Second component (for gradient)
            real_part = dis(gen) * decay * 0.5;
            imag_part = dis(gen) * decay * 0.5;
            std::get<1>(MM[idx]) = cmplx(real_part, imag_part);
        }
    }
}

/*
 * Test M2L: Compare simplified vs full plane wave expansion
 */
void testM2L() {
    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    std::cout << "=== Testing M2L Translation ===" << std::endl;
    std::cout << "Expansion order N = " << N << std::endl;

    // Create source multipole expansion with synthetic coefficients
    MM_Type source_MM;
    createSyntheticMultipoleExpansion<N>(source_MM);

    // Test translation vector (box in interaction list)
    Tddd translation = {2.5, 0.5, 3.0};  // Clearly in Uplist

    std::cout << "\nTranslation vector: (" << translation[0] << ", "
              << translation[1] << ", " << translation[2] << ")" << std::endl;

    // Method 1: Simplified (reference)
    std::cout << "\n--- Method 1: Simplified O(p^4) ---" << std::endl;
    L_Type target_L_simple;
    for (auto& l : target_L_simple) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    auto start_simple = std::chrono::high_resolution_clock::now();

    PlaneWave::PlaneWaveExpansion<N> simple_exp;
    PlaneWave::PlaneWaveM2LOperator<N>::convertMultipoleToExponential(source_MM, simple_exp, false);

    std::array<double, 3> trans_arr = {translation[0], translation[1], translation[2]};
    PlaneWave::PlaneWaveM2LOperator<N>::translateExponentialWithSource(
        source_MM, simple_exp, trans_arr, false);

    PlaneWave::PlaneWaveM2LOperator<N>::convertExponentialToLocal(
        simple_exp, target_L_simple, false);

    auto end_simple = std::chrono::high_resolution_clock::now();
    auto duration_simple = std::chrono::duration_cast<std::chrono::microseconds>(end_simple - start_simple);

    // Debug: Check simplified result magnitude
    double simple_mag = 0.0;
    for (const auto& l : target_L_simple) {
        simple_mag += std::abs(std::get<0>(l)) + std::abs(std::get<1>(l));
    }
    std::cout << "Simplified L magnitude: " << simple_mag << std::endl;

    std::cout << "Time: " << duration_simple.count() << " μs" << std::endl;

    // Method 2: Full plane wave expansion
    std::cout << "\n--- Method 2: Full O(p^2) Plane Wave ---" << std::endl;
    L_Type target_L_full;
    for (auto& l : target_L_full) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    // Debug: Check source MM magnitude
    double source_mag = 0.0;
    for (const auto& m : source_MM) {
        source_mag += std::abs(std::get<0>(m)) + std::abs(std::get<1>(m));
    }
    std::cout << "Source MM magnitude: " << source_mag << std::endl;

    auto start_full = std::chrono::high_resolution_clock::now();

    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        source_MM, target_L_full, translation, false);

    auto end_full = std::chrono::high_resolution_clock::now();
    auto duration_full = std::chrono::duration_cast<std::chrono::microseconds>(end_full - start_full);

    // Debug: Check result magnitude
    double result_mag = 0.0;
    for (const auto& l : target_L_full) {
        result_mag += std::abs(std::get<0>(l)) + std::abs(std::get<1>(l));
    }
    std::cout << "Result L magnitude: " << result_mag << std::endl;

    std::cout << "Time: " << duration_full.count() << " μs" << std::endl;

    // Compare results
    std::cout << "\n--- Comparison ---" << std::endl;
    double max_error = 0.0;
    double total_magnitude = 0.0;

    for (int j = 0; j <= N; ++j) {
        for (int k = -j; k <= j; ++k) {
            int idx = j * j + j + k;
            cmplx L_simple_0 = std::get<0>(target_L_simple[idx]);
            cmplx L_simple_1 = std::get<1>(target_L_simple[idx]);
            cmplx L_full_0 = std::get<0>(target_L_full[idx]);
            cmplx L_full_1 = std::get<1>(target_L_full[idx]);

            double error_0 = std::abs(L_simple_0 - L_full_0);
            double error_1 = std::abs(L_simple_1 - L_full_1);
            double mag_0 = std::abs(L_simple_0);
            double mag_1 = std::abs(L_simple_1);

            max_error = std::max(max_error, std::max(error_0, error_1));
            total_magnitude += mag_0 + mag_1;
        }
    }

    double relative_error = max_error / (total_magnitude / ((N+1)*(N+1) * 2));

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "Maximum absolute error: " << max_error << std::endl;
    std::cout << "Relative error: " << relative_error << std::endl;

    if (relative_error < 1e-4) {
        std::cout << "✓ PASS: Errors within tolerance" << std::endl;
    } else {
        std::cout << "✗ FAIL: Errors too large" << std::endl;
    }

    std::cout << std::fixed << std::setprecision(2);
    double speedup = static_cast<double>(duration_simple.count()) / static_cast<double>(duration_full.count());
    std::cout << "Speedup: " << speedup << "x" << std::endl;
}

/*
 * Benchmark different precision levels
 */
void benchmarkPrecisionLevels() {
    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    std::cout << "\n=== Benchmark: Precision Levels ===" << std::endl;

    MM_Type source_MM;
    createSyntheticMultipoleExpansion<N>(source_MM);
    Tddd translation = {2.5, 0.5, 3.0};

    // 3-digit
    {
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        auto start = std::chrono::high_resolution_clock::now();
        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit3>::translateMultipoleToLocal(
            source_MM, target_L, translation, false);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "3-digit (109 exp):  " << duration.count() << " μs" << std::endl;
    }

    // 6-digit
    {
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        auto start = std::chrono::high_resolution_clock::now();
        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, translation, false);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "6-digit (558 exp):  " << duration.count() << " μs" << std::endl;
    }

    // 10-digit
    {
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        auto start = std::chrono::high_resolution_clock::now();
        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit10>::translateMultipoleToLocal(
            source_MM, target_L, translation, false);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "10-digit (1751 exp): " << duration.count() << " μs" << std::endl;
    }
}

int main() {
    std::cout << "╔═══════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Complete Plane Wave Expansion Validation Test       ║" << std::endl;
    std::cout << "║  Based on Greengard & Rokhlin (1997)                  ║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════════════════╝" << std::endl;
    std::cout << std::endl;

    try {
        testM2L();
        benchmarkPrecisionLevels();

        std::cout << "\n╔═══════════════════════════════════════════════════════╗" << std::endl;
        std::cout << "║  All tests completed successfully!                   ║" << std::endl;
        std::cout << "╚═══════════════════════════════════════════════════════╝" << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
