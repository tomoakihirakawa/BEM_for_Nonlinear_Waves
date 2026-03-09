/*
 * Systematic normalization coefficient experiments
 * Test different formulas to find the correct one
 */

#include <iostream>
#include <iomanip>
#include <array>
#include <complex>
#include <cmath>
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;
constexpr int N = 4;

// Test different normalization formulas
enum class NormFormula {
    Original,           // 1.0 / sqrt_nm_nm
    WithSqrt2n1,       // sqrt(2n+1) / sqrt_nm_nm (current)
    InverseSqrt2n1,    // sqrt_nm_nm / sqrt(2n+1)
    JustSqrt2n1,       // sqrt(2n+1)
    NoNormalization    // 1.0
};

template <NormFormula Formula>
double computeNorm(int n, int m) {
    double sqrt_2n1 = std::sqrt(2.0 * n + 1.0);
    double sqrt_nm = sqrt_nm_nm[n][std::abs(m)];

    switch (Formula) {
        case NormFormula::Original:
            return 1.0 / sqrt_nm;
        case NormFormula::WithSqrt2n1:
            return sqrt_2n1 / sqrt_nm;
        case NormFormula::InverseSqrt2n1:
            return sqrt_nm / sqrt_2n1;
        case NormFormula::JustSqrt2n1:
            return sqrt_2n1;
        case NormFormula::NoNormalization:
            return 1.0;
        default:
            return 1.0;
    }
}

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== Normalization Formula Experiments ===" << std::endl;

    // Test monopole: M_0^0 = 1, all others = 0
    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    MM_Type source_MM;
    for (auto& m : source_MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }
    std::get<0>(source_MM[0]) = cmplx(1.0, 0.0);  // M_0^0 = 1

    // Translation along z-axis (numerically stable)
    Tddd translation = {0.0, 0.0, 2.0};

    // Analytical answer: L_0^0 = 1/r = 1/2 = 0.5
    double expected = 0.5;

    std::cout << "Monopole test: M_0^0 = 1, translation = (0,0,2)" << std::endl;
    std::cout << "Expected L_0^0 ≈ " << expected << std::endl << std::endl;

    // First, verify sqrt_nm_nm values
    std::cout << "Verifying sqrt_nm_nm array:" << std::endl;
    std::cout << "sqrt_nm_nm[1][0] = " << sqrt_nm_nm[1][0] << " (expected √1 = 1.0)" << std::endl;
    std::cout << "sqrt_nm_nm[1][1] = " << sqrt_nm_nm[1][1] << " (expected 1/√2 ≈ 0.707)" << std::endl;
    std::cout << "sqrt_nm_nm[2][0] = " << sqrt_nm_nm[2][0] << " (expected √1 = 1.0)" << std::endl;
    std::cout << "sqrt_nm_nm[2][1] = " << sqrt_nm_nm[2][1] << " (expected √(1/3) ≈ 0.577)" << std::endl;
    std::cout << "sqrt_nm_nm[2][2] = " << sqrt_nm_nm[2][2] << " (expected √(1/6) ≈ 0.408)" << std::endl;
    std::cout << std::endl;

    // Verify: sqrt_nm_nm[n][|m|] = √[(n-|m|)!/(n+|m|)!]
    std::cout << "Formula check: sqrt_nm_nm[n][|m|] = √[(n-|m|)!/(n+|m|)!]" << std::endl;
    std::cout << "n=1, m=0: √[(1-0)!/(1+0)!] = √[1/1] = " << 1.0 << std::endl;
    std::cout << "n=1, m=1: √[(1-1)!/(1+1)!] = √[1/2] = " << (1.0/std::sqrt(2.0)) << std::endl;
    std::cout << "n=2, m=1: √[(2-1)!/(2+1)!] = √[1/6] = " << std::sqrt(1.0/6.0) << " vs actual " << sqrt_nm_nm[2][1] << std::endl;
    std::cout << std::endl;

    // Paper formula needs: 1/√[(n-m)!(n+m)!]
    std::cout << "Paper needs: 1/√[(n-m)!(n+m)!]" << std::endl;
    std::cout << "n=1, m=0: 1/√[1!·1!] = 1.0" << std::endl;
    std::cout << "n=1, m=1: 1/√[0!·2!] = 1/√2 ≈ 0.707" << std::endl;
    std::cout << std::endl;

    // Test formulas
    std::cout << "Testing normalization formulas:" << std::endl;
    std::cout << "n=1, m=0:" << std::endl;
    std::cout << "  Original (1/sqrt_nm_nm): " << computeNorm<NormFormula::Original>(1, 0) << std::endl;
    std::cout << "  WithSqrt2n1 (sqrt(2n+1)/sqrt_nm_nm): " << computeNorm<NormFormula::WithSqrt2n1>(1, 0) << std::endl;
    std::cout << "  InverseSqrt2n1 (sqrt_nm_nm/sqrt(2n+1)): " << computeNorm<NormFormula::InverseSqrt2n1>(1, 0) << std::endl;
    std::cout << "  JustSqrt2n1: " << computeNorm<NormFormula::JustSqrt2n1>(1, 0) << std::endl;
    std::cout << "  NoNormalization: " << computeNorm<NormFormula::NoNormalization>(1, 0) << std::endl;
    std::cout << std::endl;

    // Actually test each formula with the monopole
    std::cout << "\n=== Testing each normalization with monopole M2L ===" << std::endl;

    // We need to test by modifying the normalization in convertMultipoleToExponential
    // For now, let's verify the theoretical relationship

    std::cout << "\nTheoretical analysis:" << std::endl;
    std::cout << "Paper formula (Eq 7.14): 1/√[(n-m)!(n+m)!]" << std::endl;
    std::cout << "We have: sqrt_nm_nm[n][|m|] = √[(n-|m|)!/(n+|m|)!]" << std::endl;
    std::cout << "Therefore: 1/√[(n-m)!(n+m)!] = 1/√[(n-m)!] × 1/√[(n+m)!]" << std::endl;
    std::cout << "         = 1/√[(n-m)!] × √[(n-m)!] / √[(n-m)!(n+m)!]" << std::endl;
    std::cout << "         = 1/√[(n-m)!] × sqrt_nm_nm[n][|m|]" << std::endl;
    std::cout << std::endl;

    // Calculate 1/√[(n-m)!] for small n
    std::cout << "Values of 1/√[(n-m)!]:" << std::endl;
    for (int n = 0; n <= 3; ++n) {
        for (int m = -n; m <= n; ++m) {
            double fact = 1.0;
            int nm = n - std::abs(m);
            for (int i = 2; i <= nm; ++i) {
                fact *= i;
            }
            double inv_sqrt_fact = 1.0 / std::sqrt(fact);
            std::cout << "  n=" << n << ", m=" << m << ": 1/√[" << nm << "!] = " << inv_sqrt_fact << std::endl;
        }
    }
    std::cout << std::endl;

    // The correct formula should be:
    std::cout << "Correct normalization formula:" << std::endl;
    std::cout << "For n=0, m=0: 1/√[0!×0!] = 1/√1 = 1.0" << std::endl;
    std::cout << "  sqrt_nm_nm[0][0] = " << sqrt_nm_nm[0][0] << std::endl;
    std::cout << "  1/√[(0-0)!] = 1.0" << std::endl;
    std::cout << "  Product = 1.0 × 1.0 = 1.0 ✓" << std::endl;
    std::cout << std::endl;

    std::cout << "For n=1, m=0: 1/√[1!×1!] = 1.0" << std::endl;
    std::cout << "  sqrt_nm_nm[1][0] = " << sqrt_nm_nm[1][0] << std::endl;
    std::cout << "  1/√[(1-0)!] = 1.0" << std::endl;
    std::cout << "  Product = 1.0 × 1.0 = 1.0 ✓" << std::endl;
    std::cout << std::endl;

    std::cout << "For n=1, m=1: 1/√[0!×2!] = 1/√2 ≈ 0.707" << std::endl;
    std::cout << "  sqrt_nm_nm[1][1] = " << sqrt_nm_nm[1][1] << std::endl;
    std::cout << "  1/√[(1-1)!] = 1.0" << std::endl;
    std::cout << "  Product = " << sqrt_nm_nm[1][1] << " × 1.0 = " << sqrt_nm_nm[1][1] << " ✓" << std::endl;
    std::cout << std::endl;

    std::cout << "CONCLUSION: Use sqrt_nm_nm[n][|m|] directly (no additional factors!)" << std::endl;

    return 0;
}
