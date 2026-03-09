/*
 * Test dipole moment with full plane wave expansion
 * Dipole has M_1^0 (z-direction) or M_1^±1 (x,y-directions)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;
constexpr int N = 8;
constexpr double PI = std::numbers::pi;

int main() {
    std::cout << std::scientific << std::setprecision(8);
    std::cout << "=== Dipole Moment Test (Full Plane Wave) ===" << std::endl;
    std::cout << "Expansion order N = " << N << std::endl << std::endl;

    // Test z-direction dipole: M_1^0 = 1
    {
        std::cout << "--- Z-direction Dipole (M_1^0 = 1) ---" << std::endl;

        using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        MM_Type source_MM;
        for (auto& m : source_MM) {
            std::get<0>(m) = cmplx(0.0);
            std::get<1>(m) = cmplx(0.0);
        }

        // M_1^0 corresponds to index: n=1, m=0 -> idx = 1*1 + 1 + 0 = 2
        double dipole_strength = 1.0;
        std::get<0>(source_MM[2]) = cmplx(dipole_strength, 0.0);  // M_1^0 = 1

        std::cout << "Source: M_1^0 = " << std::get<0>(source_MM[2]) << std::endl;

        // Translation along z-axis
        Tddd translation = {0.0, 0.0, 3.0};
        double r = translation[2];

        std::cout << "Translation: (0, 0, " << r << ")" << std::endl;

        // For dipole along z-axis:
        // φ(r) = M_1^0 * cos(θ) / r²
        // At z-axis (θ=0): φ = M_1^0 / r²
        // Local expansion: L_0^0 = M_1^0 / r², L_1^0 = -M_1^0 / r³

        double analytical_L00 = dipole_strength * 1.0 / (r * r);  // cos(0) = 1
        double analytical_L10 = -dipole_strength / (r * r * r);

        std::cout << "\nAnalytical:" << std::endl;
        std::cout << "  L_0^0 = " << analytical_L00 << std::endl;
        std::cout << "  L_1^0 = " << analytical_L10 << std::endl;

        // Compute using full plane wave expansion
        using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, translation, false);

        cmplx computed_L00 = std::get<0>(target_L[0]);     // idx = 0*0 + 0 + 0 = 0
        cmplx computed_L10 = std::get<0>(target_L[2]);     // idx = 1*1 + 1 + 0 = 2

        std::cout << "\nComputed:" << std::endl;
        std::cout << "  L_0^0 = " << computed_L00 << std::endl;
        std::cout << "  L_1^0 = " << computed_L10 << std::endl;

        double err_L00 = std::abs(computed_L00 - analytical_L00) / std::abs(analytical_L00);
        double err_L10 = std::abs(computed_L10 - analytical_L10) / std::abs(analytical_L10);

        std::cout << "\nRelative errors:" << std::endl;
        std::cout << "  L_0^0: " << err_L00 << std::endl;
        std::cout << "  L_1^0: " << err_L10 << std::endl;

        // Check other significant coefficients
        std::cout << "\nOther coefficients:" << std::endl;
        for (int n = 0; n <= std::min(3, N); ++n) {
            for (int m = -n; m <= n; ++m) {
                int idx = n * n + n + m;
                cmplx val = std::get<0>(target_L[idx]);
                if (std::abs(val) > 1e-10) {
                    std::cout << "  L_" << n << "^" << m << " = " << val << std::endl;
                }
            }
        }

        if (err_L00 < 1e-4 && err_L10 < 1e-4) {
            std::cout << "\n✓ PASS" << std::endl;
        } else {
            std::cout << "\n✗ FAIL" << std::endl;
        }
    }

    // Test off-axis dipole
    {
        std::cout << "\n\n--- Off-axis Z-Dipole ---" << std::endl;

        using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        MM_Type source_MM;
        for (auto& m : source_MM) {
            std::get<0>(m) = cmplx(0.0);
            std::get<1>(m) = cmplx(0.0);
        }

        double dipole_strength = 1.0;
        std::get<0>(source_MM[2]) = cmplx(dipole_strength, 0.0);  // M_1^0 = 1

        // Translation off z-axis
        Tddd translation = {1.0, 0.0, 2.0};
        double x = translation[0], y = translation[1], z = translation[2];
        double r = std::sqrt(x*x + y*y + z*z);
        double theta = std::acos(z / r);
        double cos_theta = z / r;

        std::cout << "Translation: (" << x << ", " << y << ", " << z << ")" << std::endl;
        std::cout << "r = " << r << ", θ = " << theta << " rad" << std::endl;

        // Dipole potential: φ = M_1^0 * cos(θ) / r²
        double analytical_potential = dipole_strength * cos_theta / (r * r);
        std::cout << "Analytical potential at target: " << analytical_potential << std::endl;

        // Compute using full plane wave expansion
        using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, translation, false);

        // The L_0^0 coefficient gives the leading order potential
        cmplx computed_L00 = std::get<0>(target_L[0]);
        std::cout << "Computed L_0^0: " << computed_L00 << std::endl;

        double err = std::abs(computed_L00 - analytical_potential) / std::abs(analytical_potential);
        std::cout << "Relative error in L_0^0: " << err << std::endl;

        if (err < 1e-3) {
            std::cout << "✓ PASS" << std::endl;
        } else {
            std::cout << "✗ FAIL" << std::endl;
        }
    }

    // Test combined monopole + dipole
    {
        std::cout << "\n\n--- Combined Monopole + Dipole ---" << std::endl;

        using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        MM_Type source_MM;
        for (auto& m : source_MM) {
            std::get<0>(m) = cmplx(0.0);
            std::get<1>(m) = cmplx(0.0);
        }

        double monopole = 2.0;
        double dipole = 0.5;
        std::get<0>(source_MM[0]) = cmplx(monopole, 0.0);  // M_0^0
        std::get<0>(source_MM[2]) = cmplx(dipole, 0.0);    // M_1^0

        Tddd translation = {0.0, 0.0, 4.0};
        double r = translation[2];

        std::cout << "Source: M_0^0 = " << monopole << ", M_1^0 = " << dipole << std::endl;
        std::cout << "Translation: (0, 0, " << r << ")" << std::endl;

        // At z-axis: φ = M_0^0/r + M_1^0/r²
        double analytical = monopole / r + dipole / (r * r);
        std::cout << "Analytical potential: " << analytical << std::endl;

        using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
        L_Type target_L;
        for (auto& l : target_L) {
            std::get<0>(l) = cmplx(0.0);
            std::get<1>(l) = cmplx(0.0);
        }

        PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
            source_MM, target_L, translation, false);

        cmplx computed_L00 = std::get<0>(target_L[0]);
        std::cout << "Computed L_0^0: " << computed_L00 << std::endl;

        double err = std::abs(computed_L00 - analytical) / std::abs(analytical);
        std::cout << "Relative error: " << err << std::endl;

        if (err < 1e-6) {
            std::cout << "✓ PASS" << std::endl;
        } else {
            std::cout << "✗ FAIL" << std::endl;
        }
    }

    return 0;
}
