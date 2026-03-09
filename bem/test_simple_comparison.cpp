/*
 * Simple comparison test: Compare full plane wave vs simplified version
 * Using a general position to avoid numerical issues
 */

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include "../../include/lib_plane_wave_expansion.hpp"
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;
constexpr int N = 4;  // Lower order for debugging

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== Simple Comparison Test ===" << std::endl;

    // Create synthetic multipole expansion
    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    MM_Type source_MM;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-0.1, 0.1);

    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(dis(gen), dis(gen));
        std::get<1>(m) = cmplx(dis(gen) * 0.1, dis(gen) * 0.1);
    }

    // Make monopole dominant
    std::get<0>(source_MM[0]) = cmplx(1.0, 0.0);

    double source_mag = 0.0;
    for (const auto& m : source_MM) {
        source_mag += std::abs(std::get<0>(m)) + std::abs(std::get<1>(m));
    }
    std::cout << "Source MM magnitude: " << source_mag << std::endl;

    // Translation vector: well-separated, not on any axis
    Tddd translation = {2.3, 1.7, 2.9};  // Generic position
    std::cout << "Translation: (" << translation[0] << ", "
              << translation[1] << ", " << translation[2] << ")" << std::endl;

    // Method 1: Simplified (reference)
    std::cout << "\n--- Simplified M2L ---" << std::endl;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    L_Type target_L_simple;
    for (auto& l : target_L_simple) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    PlaneWave::PlaneWaveExpansion<N> exp_simple;
    PlaneWave::PlaneWaveM2LOperator<N>::convertMultipoleToExponential(
        source_MM, exp_simple, false);

    std::array<double, 3> trans_arr = {translation[0], translation[1], translation[2]};
    PlaneWave::PlaneWaveM2LOperator<N>::translateExponentialWithSource(
        source_MM, exp_simple, trans_arr, false);

    PlaneWave::PlaneWaveM2LOperator<N>::convertExponentialToLocal(
        exp_simple, target_L_simple, false);

    double simple_mag = 0.0;
    for (const auto& l : target_L_simple) {
        simple_mag += std::abs(std::get<0>(l)) + std::abs(std::get<1>(l));
    }
    std::cout << "Result L magnitude: " << simple_mag << std::endl;

    // Print first few coefficients
    std::cout << "L_0^0 = " << std::get<0>(target_L_simple[0]) << std::endl;
    if (N >= 1) {
        std::cout << "L_1^0 = " << std::get<0>(target_L_simple[2]) << std::endl;
    }

    // Method 2: Full plane wave
    std::cout << "\n--- Full Plane Wave M2L ---" << std::endl;
    L_Type target_L_full;
    for (auto& l : target_L_full) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit3>::translateMultipoleToLocal(
        source_MM, target_L_full, translation, false);

    double full_mag = 0.0;
    for (const auto& l : target_L_full) {
        full_mag += std::abs(std::get<0>(l)) + std::abs(std::get<1>(l));
    }
    std::cout << "Result L magnitude: " << full_mag << std::endl;

    // Print first few coefficients
    std::cout << "L_0^0 = " << std::get<0>(target_L_full[0]) << std::endl;
    if (N >= 1) {
        std::cout << "L_1^0 = " << std::get<0>(target_L_full[2]) << std::endl;
    }

    // Compare
    std::cout << "\n--- Comparison ---" << std::endl;
    double max_error = 0.0;
    for (int j = 0; j <= N; ++j) {
        for (int k = -j; k <= j; ++k) {
            int idx = j * j + j + k;
            cmplx L_s0 = std::get<0>(target_L_simple[idx]);
            cmplx L_f0 = std::get<0>(target_L_full[idx]);
            double err = std::abs(L_s0 - L_f0);
            max_error = std::max(max_error, err);

            if (err > 1e-6 && std::abs(L_s0) > 1e-10) {
                std::cout << "L_" << j << "^" << k << ": simple=" << L_s0
                          << ", full=" << L_f0 << ", error=" << err << std::endl;
            }
        }
    }

    std::cout << "\nMax absolute error: " << max_error << std::endl;
    std::cout << "Relative error: " << (max_error / (simple_mag / ((N+1)*(N+1)))) << std::endl;

    if (!std::isnan(simple_mag) && !std::isnan(full_mag)) {
        std::cout << "✓ Both methods produced finite results" << std::endl;
    } else {
        std::cout << "✗ NaN detected" << std::endl;
    }

    return 0;
}
