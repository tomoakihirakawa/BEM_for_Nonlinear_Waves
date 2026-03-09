/*
 * Test single monopole with full plane wave expansion
 * Compare against analytical solution
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;
constexpr int N = 8;  // Higher order for better accuracy

int main() {
    std::cout << std::scientific << std::setprecision(8);
    std::cout << "=== Single Monopole Test (Full Plane Wave) ===" << std::endl;
    std::cout << "Expansion order N = " << N << std::endl << std::endl;

    // Single monopole: M_0^0 = 1, all others = 0
    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    MM_Type source_MM;
    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }

    double monopole_strength = 1.0;
    std::get<0>(source_MM[0]) = cmplx(monopole_strength, 0.0);  // M_0^0 = 1

    std::cout << "Monopole M_0^0 = " << std::get<0>(source_MM[0]) << std::endl;

    // Test multiple translation distances
    std::vector<double> distances = {2.0, 3.0, 5.0, 10.0};

    for (double dist : distances) {
        std::cout << "\n--- Distance r = " << dist << " ---" << std::endl;

        // Translation vector along z-axis (numerically stable)
        Tddd translation = {0.0, 0.0, dist};

        // Analytical solution: φ = M_0^0 / r
        // Local expansion at origin: L_0^0 = M_0^0 / r
        double analytical_L00 = monopole_strength / dist;
        std::cout << "Analytical L_0^0 = " << analytical_L00 << std::endl;

        // Compute using full plane wave expansion
        using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, translation, false);

        cmplx computed_L00 = std::get<0>(target_L[0]);
        std::cout << "Computed L_0^0  = " << computed_L00 << std::endl;

        double abs_error = std::abs(computed_L00 - analytical_L00);
        double rel_error = abs_error / std::abs(analytical_L00);

        std::cout << "Absolute error  = " << abs_error << std::endl;
        std::cout << "Relative error  = " << rel_error << std::endl;

        // Check other coefficients (should be zero for monopole)
        double max_spurious = 0.0;
        for (int n = 1; n <= N; ++n) {
            for (int m = -n; m <= n; ++m) {
                int idx = n * n + n + m;
                double mag = std::abs(std::get<0>(target_L[idx]));
                max_spurious = std::max(max_spurious, mag);
            }
        }
        std::cout << "Max spurious coefficient: " << max_spurious << std::endl;

        if (rel_error < 1e-6) {
            std::cout << "✓ PASS" << std::endl;
        } else if (rel_error < 1e-3) {
            std::cout << "⚠ MARGINAL" << std::endl;
        } else {
            std::cout << "✗ FAIL" << std::endl;
        }
    }

    // Test off-axis translation
    std::cout << "\n--- Off-axis translation ---" << std::endl;
    Tddd translation_offaxis = {1.5, 2.0, 2.5};
    double r_offaxis = std::sqrt(1.5*1.5 + 2.0*2.0 + 2.5*2.5);
    std::cout << "Translation: (" << translation_offaxis[0] << ", "
              << translation_offaxis[1] << ", " << translation_offaxis[2] << ")" << std::endl;
    std::cout << "Distance r = " << r_offaxis << std::endl;

    double analytical_L00_offaxis = monopole_strength / r_offaxis;
    std::cout << "Analytical L_0^0 = " << analytical_L00_offaxis << std::endl;

    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    L_Type target_L_offaxis;
    for (auto& l : target_L_offaxis) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        source_MM, target_L_offaxis, translation_offaxis, false);

    cmplx computed_L00_offaxis = std::get<0>(target_L_offaxis[0]);
    std::cout << "Computed L_0^0  = " << computed_L00_offaxis << std::endl;

    double abs_error_offaxis = std::abs(computed_L00_offaxis - analytical_L00_offaxis);
    double rel_error_offaxis = abs_error_offaxis / std::abs(analytical_L00_offaxis);

    std::cout << "Absolute error  = " << abs_error_offaxis << std::endl;
    std::cout << "Relative error  = " << rel_error_offaxis << std::endl;

    if (rel_error_offaxis < 1e-3) {
        std::cout << "✓ PASS" << std::endl;
    } else {
        std::cout << "✗ FAIL" << std::endl;
    }

    return 0;
}
