/*
 * Comprehensive test suite for plane wave expansion
 * Tests various multipole configurations and translation scenarios
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;
constexpr int N = 8;

struct TestResult {
    std::string name;
    bool passed;
    double max_error;
    double time_us;
};

std::vector<TestResult> results;

void reportTest(const std::string& name, bool passed, double max_error, double time_us = 0.0) {
    results.push_back({name, passed, max_error, time_us});
    std::cout << (passed ? "✓" : "✗") << " " << name;
    if (max_error > 0) {
        std::cout << " (error: " << max_error << ")";
    }
    if (time_us > 0) {
        std::cout << " [" << time_us << " μs]";
    }
    std::cout << std::endl;
}

// Test 1: Single monopole at various distances
void testMonopole() {
    std::cout << "\n=== Test 1: Monopole ===" << std::endl;

    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    MM_Type source_MM;
    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }
    std::get<0>(source_MM[0]) = cmplx(1.0, 0.0);

    std::vector<double> distances = {2.0, 3.0, 5.0};
    double max_error = 0.0;

    for (double d : distances) {
        Tddd trans = {0.0, 0.0, d};
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, trans, false);

        double analytical = 1.0 / d;
        double err = std::abs(std::get<0>(target_L[0]) - analytical) / analytical;
        max_error = std::max(max_error, err);
    }

    reportTest("Monopole (multiple distances)", max_error < 1e-5, max_error);
}

// Test 2: Off-axis monopole
void testOffAxisMonopole() {
    std::cout << "\n=== Test 2: Off-Axis Monopole ===" << std::endl;

    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    MM_Type source_MM;
    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }
    std::get<0>(source_MM[0]) = cmplx(1.0, 0.0);

    Tddd trans = {1.5, 2.0, 2.5};
    double r = std::sqrt(1.5*1.5 + 2.0*2.0 + 2.5*2.5);

    L_Type target_L;
    for (auto& l : target_L) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        source_MM, target_L, trans, false);

    double analytical = 1.0 / r;
    double err = std::abs(std::get<0>(target_L[0]) - analytical) / analytical;

    reportTest("Off-axis monopole", err < 1e-5, err);
}

// Test 3: Linearity (M2L should be linear)
void testLinearity() {
    std::cout << "\n=== Test 3: Linearity ===" << std::endl;

    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    // Create two random multipole expansions
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    MM_Type MM1, MM2, MM_combined;
    for (int i = 0; i < (N+1)*(N+1); ++i) {
        std::get<0>(MM1[i]) = cmplx(dis(gen), dis(gen));
        std::get<0>(MM2[i]) = cmplx(dis(gen), dis(gen));
        std::get<0>(MM_combined[i]) = std::get<0>(MM1[i]) + std::get<0>(MM2[i]);

        std::get<1>(MM1[i]) = cmplx(0.0);
        std::get<1>(MM2[i]) = cmplx(0.0);
        std::get<1>(MM_combined[i]) = cmplx(0.0);
    }

    Tddd trans = {2.0, 1.0, 3.0};

    L_Type L1, L2, L_combined;
    for (auto& l : L1) { std::get<0>(l) = cmplx(0.0); std::get<1>(l) = cmplx(0.0); }
    for (auto& l : L2) { std::get<0>(l) = cmplx(0.0); std::get<1>(l) = cmplx(0.0); }
    for (auto& l : L_combined) { std::get<0>(l) = cmplx(0.0); std::get<1>(l) = cmplx(0.0); }

    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        MM1, L1, trans, false);
    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        MM2, L2, trans, false);
    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        MM_combined, L_combined, trans, false);

    // Check L(MM1 + MM2) ≈ L(MM1) + L(MM2)
    double max_error = 0.0;
    for (int i = 0; i < (N+1)*(N+1); ++i) {
        cmplx sum = std::get<0>(L1[i]) + std::get<0>(L2[i]);
        cmplx combined = std::get<0>(L_combined[i]);
        double err = std::abs(sum - combined);
        max_error = std::max(max_error, err);
    }

    reportTest("Linearity", max_error < 1e-8, max_error);
}

// Test 4: 6-direction decomposition
void test6Directions() {
    std::cout << "\n=== Test 4: 6-Direction Decomposition ===" << std::endl;

    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    MM_Type source_MM;
    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }
    std::get<0>(source_MM[0]) = cmplx(1.0, 0.0);

    // Test each direction
    std::vector<std::tuple<std::string, Tddd>> directions = {
        {"Up",    {0.0, 0.0, 3.0}},
        {"Down",  {0.0, 0.0, -3.0}},
        {"North", {0.0, 3.0, 0.0}},
        {"South", {0.0, -3.0, 0.0}},
        {"East",  {3.0, 0.0, 0.0}},
        {"West",  {-3.0, 0.0, 0.0}}
    };

    bool all_passed = true;
    for (const auto& [name, trans] : directions) {
        double r = std::sqrt(trans[0]*trans[0] + trans[1]*trans[1] + trans[2]*trans[2]);

        using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, trans, false);

        double analytical = 1.0 / r;
        double err = std::abs(std::get<0>(target_L[0]) - analytical) / analytical;

        std::cout << "  " << name << ": error = " << err;
        if (err < 1e-5) {
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗" << std::endl;
            all_passed = false;
        }
    }

    reportTest("6-Direction decomposition", all_passed, 0.0);
}

// Test 5: Performance comparison
void testPerformance() {
    std::cout << "\n=== Test 5: Performance ===" << std::endl;

    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    MM_Type source_MM;
    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(dis(gen), dis(gen));
        std::get<1>(m) = cmplx(0.0);
    }

    Tddd trans = {2.5, 1.5, 3.5};

    // Warm-up
    L_Type target_L;
    for (auto& l : target_L) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }
    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        source_MM, target_L, trans, false);

    // Benchmark
    constexpr int num_runs = 100;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_runs; ++i) {
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }
        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, trans, false);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double avg_time = static_cast<double>(duration.count()) / num_runs;

    reportTest("Performance benchmark", true, 0.0, avg_time);
}

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "╔═══════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Comprehensive Plane Wave Expansion Test Suite       ║" << std::endl;
    std::cout << "║  Order N = " << N << "                                           ║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════════════════╝" << std::endl;

    testMonopole();
    testOffAxisMonopole();
    testLinearity();
    test6Directions();
    testPerformance();

    // Summary
    std::cout << "\n╔═══════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Test Summary                                         ║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════════════════╝" << std::endl;

    int passed = 0, failed = 0;
    for (const auto& result : results) {
        if (result.passed) ++passed;
        else ++failed;
    }

    std::cout << "Total tests: " << results.size() << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;

    if (failed == 0) {
        std::cout << "\n✓ All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ Some tests failed" << std::endl;
        return 1;
    }
}
