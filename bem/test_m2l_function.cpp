/*
 * Simple test of SphericalCoordinates::m2lFunction
 */

#include <iostream>
#include <iomanip>
#include "../../include/lib_multipole_expansion.hpp"

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== m2lFunction Test ===" << std::endl;

    // Translation vector
    Tddd trans_vec = {0.0, 0.0, 2.0};
    std::cout << "Translation: (" << trans_vec[0] << ", "
              << trans_vec[1] << ", " << trans_vec[2] << ")" << std::endl;

    SphericalCoordinates sph(trans_vec);
    std::cout << "rho: " << sph.rho << ", theta: " << sph.theta
              << ", phi: " << sph.phi << std::endl << std::endl;

    // Test m2lFunction for various (j,k,n,m)
    std::cout << "Testing m2lFunction(j, k, n, m):" << std::endl;

    auto test = [&](int j, int k, int n, int m) {
        auto result = sph.m2lFunction(j, k, n, m);
        std::cout << "  (" << j << "," << k << "," << n << "," << m << "): "
                  << result;
        if (std::isnan(std::abs(result))) {
            std::cout << " <- NaN!";
        }
        std::cout << std::endl;
    };

    test(0, 0, 0, 0);
    test(1, 0, 0, 0);
    test(0, 0, 1, 0);
    test(1, 0, 1, 0);
    test(1, 0, 0, 1);  // This one caused nan before
    test(1, 1, 0, 0);
    test(1, -1, 0, 0);

    // Check if macro is defined
    #ifdef GreenGardAndRokhlin1997
    std::cout << "\n✓ GreenGardAndRokhlin1997 is defined" << std::endl;
    #else
    std::cout << "\n✗ GreenGardAndRokhlin1997 is NOT defined" << std::endl;
    #endif

    #ifdef Nishimura2002
    std::cout << "✓ Nishimura2002 is defined" << std::endl;
    #else
    std::cout << "✗ Nishimura2002 is NOT defined" << std::endl;
    #endif

    return 0;
}
