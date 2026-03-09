/*
 * Simple test with a single monopole charge to validate formulas
 * Analytical solution: Φ(r) = q/|r-r0|
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include "../../include/lib_multipole_expansion.hpp"
#include "../../include/lib_plane_wave_expansion.hpp"
#include "../../include/lib_plane_wave_expansion_full.hpp"

using cmplx = std::complex<double>;
constexpr int N = 4;  // Lower order for easier debugging

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== Single Monopole Test ===" << std::endl;
    std::cout << "Expansion order N = " << N << std::endl;

    // Single monopole at origin with charge q = 1.0
    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    MM_Type source_MM;

    // Initialize to zero
    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }

    // For a monopole at origin, only M_0^0 is non-zero
    // M_0^0 = q * Y_0^0(0,0) * r^0 = q * Y_0^0
    // Y_0^0 = 1/(2*sqrt(π)) according to normalization
    double Y00 = 1.0 / (2.0 * std::sqrt(std::numbers::pi));
    std::get<0>(source_MM[0]) = cmplx(1.0 * Y00, 0.0);  // q=1, only real part

    std::cout << "\nMonopole M_0^0 = " << std::get<0>(source_MM[0]) << std::endl;

    // Evaluation point (target box center)
    Tddd target_pos = {2.0, 0.0, 0.0};  // Distance = 2.0 from origin
    double distance = std::sqrt(target_pos[0]*target_pos[0] +
                                target_pos[1]*target_pos[1] +
                                target_pos[2]*target_pos[2]);

    std::cout << "Target position: (" << target_pos[0] << ", "
              << target_pos[1] << ", " << target_pos[2] << ")" << std::endl;
    std::cout << "Distance: " << distance << std::endl;

    // Analytical solution for monopole: Φ = q/r = 1/2 = 0.5
    double phi_analytical = 1.0 / distance;
    std::cout << "\nAnalytical potential: " << phi_analytical << std::endl;

    // Test 1: Direct evaluation using M2L
    std::cout << "\n--- Direct M2L (Simplified) ---" << std::endl;

    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    L_Type target_L_direct;
    for (auto& l : target_L_direct) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    // Use simplified M2L
    PlaneWave::PlaneWaveExpansion<N> exp_simple;
    PlaneWave::PlaneWaveM2LOperator<N>::convertMultipoleToExponential(
        source_MM, exp_simple, false);

    std::array<double, 3> trans_arr = {target_pos[0], target_pos[1], target_pos[2]};
    PlaneWave::PlaneWaveM2LOperator<N>::translateExponentialWithSource(
        source_MM, exp_simple, trans_arr, false);

    PlaneWave::PlaneWaveM2LOperator<N>::convertExponentialToLocal(
        exp_simple, target_L_direct, false);

    // Evaluate potential at target center (r=0 in local coordinates)
    // For local expansion: Φ(r) = Σ L_n^m * r^n * Y_n^m
    // At r=0, only L_0^0 contributes
    cmplx L00 = std::get<0>(target_L_direct[0]);
    std::cout << "L_0^0 = " << L00 << std::endl;

    // Φ(0) = L_0^0 * Y_0^0 = L_0^0 * 1/(2√π)
    double phi_m2l = std::real(L00) * Y00;
    std::cout << "M2L potential: " << phi_m2l << std::endl;
    std::cout << "Error: " << std::abs(phi_m2l - phi_analytical) << std::endl;
    std::cout << "Relative error: " << std::abs((phi_m2l - phi_analytical) / phi_analytical) << std::endl;

    // Print all local expansion coefficients to see the pattern
    std::cout << "\nAll local expansion coefficients:" << std::endl;
    for (int n = 0; n <= std::min(2, N); ++n) {
        for (int m = -n; m <= n; ++m) {
            int idx = n * n + n + m;
            cmplx L_nm = std::get<0>(target_L_direct[idx]);
            if (std::abs(L_nm) > 1e-10) {
                std::cout << "L_" << n << "^" << m << " = " << L_nm
                          << " (mag: " << std::abs(L_nm) << ")" << std::endl;
            }
        }
    }

    return 0;
}
