/*
 * Direct test of spherical harmonics computation
 */

#include <iostream>
#include <iomanip>
#include "../../include/lib_multipole_expansion.hpp"

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== Spherical Harmonics Test ===" << std::endl << std::endl;

    // Test case 1: Z-axis
    {
        Tddd trans_vec = {0.0, 0.0, 2.0};
        SphericalCoordinates sph(trans_vec);
        std::cout << "Z-axis translation (0, 0, 2):" << std::endl;
        std::cout << "  rho=" << sph.rho << ", theta=" << sph.theta << ", phi=" << sph.phi << std::endl;

        // Direct computation without precomputation
        auto Y00 = sph.sph_harmonics_(0, 0);
        auto Y10 = sph.sph_harmonics_(1, 0);
        auto Y11 = sph.sph_harmonics_(1, 1);
        auto Y1m1 = sph.sph_harmonics_(1, -1);

        std::cout << "  Y(0,0)  = " << Y00 << std::endl;
        std::cout << "  Y(1,0)  = " << Y10 << std::endl;
        std::cout << "  Y(1,1)  = " << Y11 << (std::isnan(std::abs(Y11)) ? " <- NaN" : "") << std::endl;
        std::cout << "  Y(1,-1) = " << Y1m1 << (std::isnan(std::abs(Y1m1)) ? " <- NaN" : "") << std::endl;

        // With precomputation
        sph.precompute_sph_div_rhon1(2);
        auto Y11_div = sph.sph_harmonics_div_rhon1(1, 1);
        std::cout << "  Y(1,1)/rho^2 (precomputed) = " << Y11_div
                  << (std::isnan(std::abs(Y11_div)) ? " <- NaN" : "") << std::endl;
        std::cout << std::endl;
    }

    // Test case 2: Off-axis
    {
        Tddd trans_vec = {1.0, 1.0, 2.0};
        SphericalCoordinates sph(trans_vec);
        std::cout << "Off-axis translation (1, 1, 2):" << std::endl;
        std::cout << "  rho=" << sph.rho << ", theta=" << sph.theta << ", phi=" << sph.phi << std::endl;

        auto Y00 = sph.sph_harmonics_(0, 0);
        auto Y10 = sph.sph_harmonics_(1, 0);
        auto Y11 = sph.sph_harmonics_(1, 1);
        auto Y1m1 = sph.sph_harmonics_(1, -1);

        std::cout << "  Y(0,0)  = " << Y00 << std::endl;
        std::cout << "  Y(1,0)  = " << Y10 << std::endl;
        std::cout << "  Y(1,1)  = " << Y11 << (std::isnan(std::abs(Y11)) ? " <- NaN" : "") << std::endl;
        std::cout << "  Y(1,-1) = " << Y1m1 << (std::isnan(std::abs(Y1m1)) ? " <- NaN" : "") << std::endl;

        // With precomputation
        sph.precompute_sph_div_rhon1(2);
        auto Y11_div = sph.sph_harmonics_div_rhon1(1, 1);
        std::cout << "  Y(1,1)/rho^2 (precomputed) = " << Y11_div
                  << (std::isnan(std::abs(Y11_div)) ? " <- NaN" : "") << std::endl;
        std::cout << std::endl;
    }

    // Test case 3: Check AAA_M2L_FMM coefficient
    {
        Tddd trans_vec = {0.0, 0.0, 2.0};
        SphericalCoordinates sph(trans_vec);
        sph.precompute_sph_div_rhon1(2);

        std::cout << "M2L Function Components:" << std::endl;
        int j = 1, k = 0, n = 0, m = 1;
        int p = j + n;  // 1
        int q = m - k;  // 1

        std::cout << "  Testing (j=" << j << ", k=" << k << ", n=" << n << ", m=" << m << ")" << std::endl;
        std::cout << "  → p = j+n = " << p << ", q = m-k = " << q << std::endl;

        auto AAA_coeff = AAA_M2L_FMM[j][k + N_AAA_M2L_FMM][n][m + N_AAA_M2L_FMM];
        std::cout << "  AAA_M2L_FMM[" << j << "][" << k << "][" << n << "][" << m << "] = "
                  << AAA_coeff << std::endl;

        auto Y_pq = sph.sph_harmonics_div_rhon1(p, q);
        std::cout << "  Y(" << p << "," << q << ")/rho^" << (p+1) << " = " << Y_pq
                  << (std::isnan(std::abs(Y_pq)) ? " <- NaN" : "") << std::endl;

        auto result = AAA_coeff * Y_pq;
        std::cout << "  Result = AAA * Y = " << result
                  << (std::isnan(std::abs(result)) ? " <- NaN" : "") << std::endl;
        std::cout << std::endl;
    }

    return 0;
}
