/*
 * Test to verify that axis-aligned translations cause NaN
 * while off-axis translations work correctly
 */

#include <iostream>
#include <iomanip>
#include "../../include/lib_multipole_expansion.hpp"

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== Axis Alignment Test ===" << std::endl << std::endl;

    // Test various translation vectors
    std::vector<std::pair<std::string, Tddd>> test_cases = {
        {"Z-axis (0, 0, 2.0)", {0.0, 0.0, 2.0}},
        {"Off-axis (0.1, 0.1, 2.0)", {0.1, 0.1, 2.0}},
        {"Off-axis (1.0, 1.0, 2.0)", {1.0, 1.0, 2.0}},
        {"X-axis (2.0, 0, 0)", {2.0, 0.0, 0.0}},
        {"Y-axis (0, 2.0, 0)", {0.0, 2.0, 0.0}},
        {"Diagonal (1.0, 1.0, 1.0)", {1.0, 1.0, 1.0}}
    };

    for (const auto& [name, trans_vec] : test_cases) {
        std::cout << "Testing: " << name << std::endl;
        std::cout << "  Vector: (" << trans_vec[0] << ", "
                  << trans_vec[1] << ", " << trans_vec[2] << ")" << std::endl;

        SphericalCoordinates sph(trans_vec);
        std::cout << "  rho: " << sph.rho << ", theta: " << sph.theta
                  << ", phi: " << sph.phi << std::endl;

        // Test critical cases
        bool has_nan = false;

        // Case 1: (j=1, k=0, n=0, m=1) - requires Y(1,1)
        auto result1 = sph.m2lFunction(1, 0, 0, 1);
        if (std::isnan(std::abs(result1))) {
            std::cout << "  ✗ m2lFunction(1,0,0,1) = NaN" << std::endl;
            has_nan = true;
        } else {
            std::cout << "  ✓ m2lFunction(1,0,0,1) = " << result1 << std::endl;
        }

        // Case 2: (j=1, k=1, n=0, m=0) - requires Y(1,-1)
        auto result2 = sph.m2lFunction(1, 1, 0, 0);
        if (std::isnan(std::abs(result2))) {
            std::cout << "  ✗ m2lFunction(1,1,0,0) = NaN" << std::endl;
            has_nan = true;
        } else {
            std::cout << "  ✓ m2lFunction(1,1,0,0) = " << result2 << std::endl;
        }

        // Case 3: (j=0, k=0, n=0, m=0) - always works (monopole)
        auto result3 = sph.m2lFunction(0, 0, 0, 0);
        if (std::isnan(std::abs(result3))) {
            std::cout << "  ✗ m2lFunction(0,0,0,0) = NaN" << std::endl;
            has_nan = true;
        } else {
            std::cout << "  ✓ m2lFunction(0,0,0,0) = " << result3 << std::endl;
        }

        if (has_nan) {
            std::cout << "  → FAILED: Contains NaN" << std::endl;
        } else {
            std::cout << "  → PASSED: All finite" << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "╔════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Summary                                                       ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════╝" << std::endl;
    std::cout << std::endl;
    std::cout << "Expected behavior:" << std::endl;
    std::cout << "  ✗ Z-axis, X-axis, Y-axis: Should produce NaN for m≠0 terms" << std::endl;
    std::cout << "  ✓ Off-axis, Diagonal: Should work correctly" << std::endl;
    std::cout << std::endl;
    std::cout << "This demonstrates that the bug is in SphericalCoordinates," << std::endl;
    std::cout << "specifically when theta=0 or theta=π/2 (coordinate axes)." << std::endl;

    return 0;
}
