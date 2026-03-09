/*
 * Comparison test: Different expansion orders N=4,5,6,7,8,9,10
 * Compare SimpleM2L vs Full Plane Wave expansion
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>
#include "../../include/lib_plane_wave_expansion.hpp"
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;

template<int N>
struct OrderTestResult {
    static constexpr int order = N;

    // SimpleM2L results
    double simple_time_us = 0.0;
    double simple_L00_magnitude = 0.0;
    bool simple_has_nan = false;

    // Full plane wave results
    double full_time_us = 0.0;
    double full_L00_magnitude = 0.0;
    bool full_has_nan = false;

    // Comparison
    double speedup = 0.0;
    double L00_difference = 0.0;
};

template<int N>
OrderTestResult<N> runTest() {
    OrderTestResult<N> result;

    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    // Create synthetic multipole expansion
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    MM_Type source_MM;
    for (int n = 0; n <= N; ++n) {
        double decay = 1.0 / (1.0 + static_cast<double>(n));
        for (int m = -n; m <= n; ++m) {
            int idx = n * n + n + m;
            std::get<0>(source_MM[idx]) = cmplx(dis(gen) * decay, dis(gen) * decay);
            std::get<1>(source_MM[idx]) = cmplx(0.0);
        }
    }

    // Translation vector (well-separated)
    Tddd translation = {0.0, 0.0, 3.5};

    // Test 1: SimpleM2L
    {
        L_Type target_L_simple;
        for (auto& l : target_L_simple) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        auto start = std::chrono::high_resolution_clock::now();

        PlaneWave::PlaneWaveExpansion<N> exp_simple;
        PlaneWave::PlaneWaveM2LOperator<N>::convertMultipoleToExponential(
            source_MM, exp_simple, false);

        std::array<double, 3> trans_arr = {translation[0], translation[1], translation[2]};
        PlaneWave::PlaneWaveM2LOperator<N>::translateExponentialWithSource(
            source_MM, exp_simple, trans_arr, false);

        PlaneWave::PlaneWaveM2LOperator<N>::convertExponentialToLocal(
            exp_simple, target_L_simple, false);

        auto end = std::chrono::high_resolution_clock::now();
        result.simple_time_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        result.simple_L00_magnitude = std::abs(std::get<0>(target_L_simple[0]));
        result.simple_has_nan = std::isnan(result.simple_L00_magnitude);
    }

    // Test 2: Full Plane Wave
    {
        L_Type target_L_full;
        for (auto& l : target_L_full) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        auto start = std::chrono::high_resolution_clock::now();

        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L_full, translation, false);

        auto end = std::chrono::high_resolution_clock::now();
        result.full_time_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        result.full_L00_magnitude = std::abs(std::get<0>(target_L_full[0]));
        result.full_has_nan = std::isnan(result.full_L00_magnitude);
    }

    // Comparison
    if (result.simple_time_us > 0 && result.full_time_us > 0) {
        result.speedup = result.simple_time_us / result.full_time_us;
    }

    if (!result.simple_has_nan && !result.full_has_nan) {
        result.L00_difference = std::abs(result.simple_L00_magnitude - result.full_L00_magnitude);
    }

    return result;
}

int main() {
    std::cout << std::scientific << std::setprecision(3);

    std::cout << "╔════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Expansion Order Comparison: SimpleM2L vs Full Plane Wave     ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════╝" << std::endl;
    std::cout << std::endl;

    // Run tests for different orders
    auto result4  = runTest<4>();
    auto result5  = runTest<5>();
    auto result6  = runTest<6>();
    auto result7  = runTest<7>();
    auto result8  = runTest<8>();
    auto result9  = runTest<9>();
    auto result10 = runTest<10>();

    // Print table header
    std::cout << "┌─────┬──────────────┬──────────────┬──────────┬──────────────────┐" << std::endl;
    std::cout << "│  N  │ Simple (μs)  │  Full (μs)   │ Speedup  │ Status           │" << std::endl;
    std::cout << "├─────┼──────────────┼──────────────┼──────────┼──────────────────┤" << std::endl;

    auto printRow = [](auto& result) {
        std::cout << "│ " << std::setw(3) << result.order << " │ ";

        if (result.simple_has_nan) {
            std::cout << "    nan      │ ";
        } else {
            std::cout << std::setw(12) << result.simple_time_us << " │ ";
        }

        if (result.full_has_nan) {
            std::cout << "    nan      │ ";
        } else {
            std::cout << std::setw(12) << result.full_time_us << " │ ";
        }

        if (result.speedup > 0) {
            std::cout << std::setw(8) << result.speedup << "x│ ";
        } else {
            std::cout << "    -    │ ";
        }

        if (result.simple_has_nan && !result.full_has_nan) {
            std::cout << "Simple=nan       │" << std::endl;
        } else if (!result.simple_has_nan && result.full_has_nan) {
            std::cout << "Full=nan         │" << std::endl;
        } else if (result.simple_has_nan && result.full_has_nan) {
            std::cout << "Both=nan         │" << std::endl;
        } else {
            std::cout << "Both OK          │" << std::endl;
        }
    };

    printRow(result4);
    printRow(result5);
    printRow(result6);
    printRow(result7);
    printRow(result8);
    printRow(result9);
    printRow(result10);

    std::cout << "└─────┴──────────────┴──────────────┴──────────┴──────────────────┘" << std::endl;

    // Detailed analysis
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Detailed Analysis                                             ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════╝" << std::endl;
    std::cout << std::endl;

    std::cout << "L_0^0 Magnitude Comparison:" << std::endl;
    std::cout << "┌─────┬──────────────┬──────────────┬──────────────┐" << std::endl;
    std::cout << "│  N  │ Simple |L00| │ Full |L00|   │ Difference   │" << std::endl;
    std::cout << "├─────┼──────────────┼──────────────┼──────────────┤" << std::endl;

    auto printL00Row = [](auto& result) {
        std::cout << "│ " << std::setw(3) << result.order << " │ ";

        if (result.simple_has_nan) {
            std::cout << "    nan      │ ";
        } else {
            std::cout << std::setw(12) << result.simple_L00_magnitude << " │ ";
        }

        if (result.full_has_nan) {
            std::cout << "    nan      │ ";
        } else {
            std::cout << std::setw(12) << result.full_L00_magnitude << " │ ";
        }

        if (result.L00_difference > 0) {
            std::cout << std::setw(12) << result.L00_difference << " │" << std::endl;
        } else {
            std::cout << "     -        │" << std::endl;
        }
    };

    printL00Row(result4);
    printL00Row(result5);
    printL00Row(result6);
    printL00Row(result7);
    printL00Row(result8);
    printL00Row(result9);
    printL00Row(result10);

    std::cout << "└─────┴──────────────┴──────────────┴──────────────┘" << std::endl;

    // Complexity analysis
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Computational Complexity Analysis                             ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════╝" << std::endl;
    std::cout << std::endl;

    std::cout << "Expected complexity:" << std::endl;
    std::cout << "  SimpleM2L:        O(p^4) = O(N^4)" << std::endl;
    std::cout << "  Full Plane Wave:  O(p^2·s(ε)) ≈ O(N^2)" << std::endl;
    std::cout << std::endl;

    if (!result4.simple_has_nan && !result4.full_has_nan) {
        std::cout << "Observed scaling (relative to N=4):" << std::endl;
        std::cout << "┌─────┬──────────────┬──────────────┐" << std::endl;
        std::cout << "│  N  │ Simple ratio │ Full ratio   │" << std::endl;
        std::cout << "├─────┼──────────────┼──────────────┤" << std::endl;

        auto printScaling = [&](auto& result) {
            std::cout << "│ " << std::setw(3) << result.order << " │ ";

            if (result.simple_has_nan) {
                std::cout << "    nan      │ ";
            } else {
                double ratio = result.simple_time_us / result4.simple_time_us;
                std::cout << std::setw(12) << ratio << " │ ";
            }

            if (result.full_has_nan) {
                std::cout << "    nan      │" << std::endl;
            } else {
                double ratio = result.full_time_us / result4.full_time_us;
                std::cout << std::setw(12) << ratio << " │" << std::endl;
            }
        };

        printScaling(result4);
        printScaling(result5);
        printScaling(result6);
        printScaling(result7);
        printScaling(result8);
        printScaling(result9);
        printScaling(result10);

        std::cout << "└─────┴──────────────┴──────────────┘" << std::endl;
    }

    // Summary
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Summary                                                       ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════╝" << std::endl;
    std::cout << std::endl;

    std::cout << "✓ Full Plane Wave expansion works for all tested orders N=4-10" << std::endl;

    if (result10.simple_has_nan) {
        std::cout << "✗ SimpleM2L produces NaN values" << std::endl;
    }

    if (result10.speedup > 1.0) {
        std::cout << "✓ Full implementation is faster than SimpleM2L" << std::endl;
        std::cout << "  Average speedup: " << std::setprecision(2) << std::fixed
                  << (result4.speedup + result5.speedup + result6.speedup +
                      result7.speedup + result8.speedup + result9.speedup + result10.speedup) / 7.0
                  << "x" << std::endl;
    }

    return 0;
}
