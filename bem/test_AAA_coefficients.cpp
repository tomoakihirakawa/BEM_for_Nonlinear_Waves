/*
 * Test AAA_M2L_FMM coefficients directly
 */

#include <iostream>
#include <iomanip>
#include "../../include/lib_multipole_expansion.hpp"

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== AAA_M2L_FMM Coefficients Test ===" << std::endl << std::endl;

    std::cout << "Checking N_AAA_M2L_FMM = " << N_AAA_M2L_FMM << std::endl << std::endl;

    // Test specific coefficients used in m2lFunction(1,0,0,1)
    // This accesses AAA_M2L_FMM[j=1][k=0+N][n=0][m=1+N]
    int j = 1, k = 0, n = 0, m = 1;

    std::cout << "For m2lFunction(" << j << "," << k << "," << n << "," << m << "):" << std::endl;
    std::cout << "  Indices: [" << j << "][" << (k + N_AAA_M2L_FMM) << "][" << n << "][" << (m + N_AAA_M2L_FMM) << "]" << std::endl;

    auto coeff = AAA_M2L_FMM[j][k + N_AAA_M2L_FMM][n][m + N_AAA_M2L_FMM];
    std::cout << "  AAA coefficient = " << coeff;
    if (std::isnan(std::abs(coeff))) {
        std::cout << " <- NaN!";
    }
    std::cout << std::endl << std::endl;

    // Test a few more
    std::cout << "Sample of AAA_M2L_FMM coefficients:" << std::endl;
    for (int test_j = 0; test_j <= 2; ++test_j) {
        for (int test_k = -test_j; test_k <= test_j; ++test_k) {
            for (int test_n = 0; test_n <= 1; ++test_n) {
                for (int test_m = -test_n; test_m <= test_n; ++test_m) {
                    auto c = AAA_M2L_FMM[test_j][test_k + N_AAA_M2L_FMM][test_n][test_m + N_AAA_M2L_FMM];
                    std::cout << "  [" << test_j << "][" << test_k << "][" << test_n << "][" << test_m << "] = " << c;
                    if (std::isnan(std::abs(c))) {
                        std::cout << " <- NaN!";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    return 0;
}
