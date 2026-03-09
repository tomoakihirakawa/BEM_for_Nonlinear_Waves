#ifndef lib_plane_wave_quadrature_H
#define lib_plane_wave_quadrature_H

/*
 * Quadrature Tables for Plane Wave Expansion
 * Based on Greengard & Rokhlin (1997), Section 12, Tables 5-7
 *
 * These tables provide discretization for the outer integral in equation (7.2)
 * using generalized Gaussian quadrature rules designed for FMM interaction lists.
 */

#include <array>
#include <vector>

namespace PlaneWaveQuadrature {

/*
 * Quadrature node structure
 */
struct QuadNode {
    double lambda;  // Node value λ_k
    double weight;  // Weight w_k
    int M_alpha;    // Number of angular discretization points M(k)
};

/*
 * Table 5: 3-digit accuracy (10^-3)
 * 9 nodes, 109 total exponentials
 */
constexpr std::array<QuadNode, 9> Table_3digit = {{
    {0.09927399673971, 0.24776441819008, 4},
    {0.47725674637049, 0.49188566500464, 7},
    {1.05533661382183, 0.65378749137677, 11},
    {1.76759343354008, 0.76433038408784, 15},
    {2.57342629351471, 0.84376180565628, 20},
    {3.44824339201583, 0.90445883985098, 20},
    {4.37680983554726, 0.95378613136833, 24},
    {5.34895757205460, 0.99670261613218, 7},
    {6.35765785313375, 1.10429422730252, 1}
}};

/*
 * Table 6: 6-digit accuracy (10^-6)
 * 18 nodes, 558 total exponentials
 */
constexpr std::array<QuadNode, 18> Table_6digit = {{
    {0.05278852766117, 0.13438265914335, 5},
    {0.26949859838931, 0.29457752727395, 8},
    {0.63220353174689, 0.42607819361148, 12},
    {1.11307564277608, 0.53189220776549, 16},
    {1.68939496140213, 0.61787306245538, 20},
    {2.34376200469530, 0.68863156078905, 25},
    {3.06269982907806, 0.74749099381426, 29},
    {3.83562941265296, 0.79699192718599, 34},
    {4.65424734321562, 0.83917454386997, 38},
    {5.51209386593581, 0.87570092283745, 43},
    {6.40421268377278, 0.90792943590067, 47},
    {7.32688001906175, 0.93698393742461, 51},
    {8.27740099258238, 0.96382546688788, 56},
    {9.25397180602489, 0.98932985769673, 59},
    {10.25560272374640, 1.01438284597917, 59},
    {11.28208829787774, 1.04003654374165, 51},
    {12.33406790967692, 1.06815489269567, 4},
    {13.41492024017240, 1.10907580975537, 1}
}};

/*
 * Table 7: 10-digit accuracy (10^-10)
 * 30 nodes, 1751 total exponentials
 */
constexpr std::array<QuadNode, 30> Table_10digit = {{
    {0.03239542384523, 0.08289159611006, 7},
    {0.16861844033714, 0.18838810673274, 10},
    {0.40611377169029, 0.28485143005306, 14},
    {0.73466473057596, 0.37041553715895, 18},
    {1.14340561998398, 0.44539043894975, 22},
    {1.62232408412252, 0.51100452150290, 26},
    {2.16276138867422, 0.56865283856139, 30},
    {2.75739199003682, 0.61958013174010, 35},
    {3.40002470112078, 0.66481004321965, 39},
    {4.08539104793552, 0.70517204769960, 43},
    {4.80897515497095, 0.74134967169016, 48},
    {5.56688915983444, 0.77392103530415, 53},
    {6.35578243654166, 0.80338600122756, 57},
    {7.17277232990713, 0.83018277269650, 62},
    {8.01538803542112, 0.85469824839953, 66},
    {8.88152313049502, 0.87727539085565, 71},
    {9.76939480982937, 0.89821948245755, 76},
    {10.67750922034750, 0.91780416582368, 80},
    {11.60463289992789, 0.93627766216629, 85},
    {12.54977061299652, 0.95386940504388, 89},
    {13.51215012257297, 0.97079739700556, 94},
    {14.49121482655196, 0.98727684670885, 97},
    {15.48662587630224, 1.00353112433459, 103},
    {16.49827659770404, 1.01980697905712, 107},
    {17.52632405530625, 1.03639774457222, 110},
    {18.57124579700721, 1.05368191266322, 112},
    {19.63393428118300, 1.07219343903929, 108},
    {20.71585163675095, 1.09278318162014, 84},
    {21.81939113866225, 1.11737373706779, 4},
    {22.95080495008893, 1.15786184931141, 1}
}};

/*
 * Precision level enum
 */
enum class Precision {
    Digit3,   // 10^-3 accuracy, 109 exponentials
    Digit6,   // 10^-6 accuracy, 558 exponentials
    Digit10   // 10^-10 accuracy, 1751 exponentials
};

/*
 * Get quadrature table for desired precision
 */
template <Precision P>
inline constexpr auto& getQuadratureTable() {
    if constexpr (P == Precision::Digit3) {
        return Table_3digit;
    } else if constexpr (P == Precision::Digit6) {
        return Table_6digit;
    } else {
        return Table_10digit;
    }
}

/*
 * Compute total number of exponentials for a given precision
 */
template <Precision P>
inline constexpr int getTotalExponentials() {
    constexpr auto& table = getQuadratureTable<P>();
    int total = 0;
    for (const auto& node : table) {
        total += node.M_alpha;
    }
    return total;
}

/*
 * Get number of nodes for a given precision
 */
template <Precision P>
inline constexpr int getNumNodes() {
    return getQuadratureTable<P>().size();
}

} // namespace PlaneWaveQuadrature

#endif // lib_plane_wave_quadrature_H
