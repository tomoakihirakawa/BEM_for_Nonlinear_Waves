/*
 * Test SphericalCoordinates initialization and m2lFunction
 */

#include <iostream>
#include <iomanip>
#include "../../include/lib_multipole_expansion.hpp"

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== SphericalCoordinates Test ===" << std::endl;

    // Test translation vector (use z-axis for numerical stability)
    Tddd trans_vec = {0.0, 0.0, 2.0};
    std::cout << "Translation vector: (" << trans_vec[0] << ", "
              << trans_vec[1] << ", " << trans_vec[2] << ")" << std::endl;

    SphericalCoordinates sph(trans_vec);
    std::cout << "rho (radius): " << sph.rho << std::endl;
    std::cout << "theta: " << sph.theta << std::endl;
    std::cout << "phi: " << sph.phi << std::endl;

    // Precompute spherical harmonics
    int max_order = 8;
    std::cout << "\nPrecomputing spherical harmonics up to order " << max_order << std::endl;
    sph.precompute_sph_div_rhon1(max_order);

    // Test m2lFunction for simple cases
    std::cout << "\nTesting m2lFunction:" << std::endl;

    // For j=0, k=0, n=0, m=0:
    // This should compute Y_(0+0)^(0-0) / r^(0+0+1) = Y_0^0 / r
    auto m2l_0000 = sph.m2lFunction(0, 0, 0, 0);
    std::cout << "m2lFunction(0,0,0,0) = " << m2l_0000 << std::endl;
    std::cout << "Expected ~= Y_0^0/r = " << (1.0/(2.0*std::sqrt(std::numbers::pi))/2.0) << std::endl;

    // Test direct spherical harmonic evaluation
    std::cout << "\nDirect spherical harmonic tests:" << std::endl;
    auto Y_1_1 = sph.sph_harmonics_div_rhon1(1, 1);
    std::cout << "Y_1^1 / r^2 = " << Y_1_1 << std::endl;

    auto Y_1_0 = sph.sph_harmonics_div_rhon1(1, 0);
    std::cout << "Y_1^0 / r^2 = " << Y_1_0 << std::endl;

    // Test some other values
    std::cout << "\nm2lFunction tests:" << std::endl;
    auto m2l_0000_again = sph.m2lFunction(0, 0, 0, 0);
    std::cout << "m2lFunction(0,0,0,0) = " << m2l_0000_again << std::endl;

    auto m2l_1000 = sph.m2lFunction(1, 0, 0, 0);
    std::cout << "m2lFunction(1,0,0,0) = " << m2l_1000 << std::endl;

    // This one causes nan
    auto m2l_1001 = sph.m2lFunction(1, 0, 0, 1);
    std::cout << "m2lFunction(1,0,0,1) = " << m2l_1001 << " (causes nan)" << std::endl;

    // Check if AAA_M2L_FMM coefficients are accessible
    #ifdef GreenGardAndRokhlin1997
    std::cout << "\nGreenGardAndRokhlin1997 is defined" << std::endl;
    std::cout << "AAA_M2L_FMM[0][N_AAA_M2L_FMM][0][N_AAA_M2L_FMM] = "
              << AAA_M2L_FMM[0][N_AAA_M2L_FMM][0][N_AAA_M2L_FMM] << std::endl;
    #else
    std::cout << "\nGreenGardAndRokhlin1997 is NOT defined" << std::endl;
    #endif

    return 0;
}
