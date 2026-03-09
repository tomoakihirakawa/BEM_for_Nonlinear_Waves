#include "rootFinding.hpp"

double norm_f(const T4d &Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    //ゼロにしたい関数
    auto [Ax, Ay, Az] = A;     // measured
    auto [Mx, My, Mz] = M;     // measured
    auto [gx0, gy0, gz0] = G0; //
    auto [mx0, my0, mz0] = M0; //
    auto [a, b, c, d] = Q;     // 計算しているクォータニオン
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    return (pow(Ax - a2 * gx0 - b2 * gx0 + c2 * gx0 + d2 * gx0 - 2 * b * c * gy0 - 2 * a * d * gy0 + 2 * a * c * gz0 - 2 * b * d * gz0, 2) + pow(Az - 2 * a * c * gx0 - 2 * b * d * gx0 + 2 * a * b * gy0 - 2 * c * d * gy0 - a2 * gz0 + b2 * gz0 + c2 * gz0 - d2 * gz0, 2) + pow(Ay + 2 * a * d * gx0 - a2 * gy0 + b2 * gy0 - c2 * gy0 + d2 * gy0 - 2 * c * d * gz0 - 2 * b * (c * gx0 + a * gz0), 2) + pow(Mx - a2 * mx0 - b2 * mx0 + c2 * mx0 + d2 * mx0 - 2 * b * c * my0 - 2 * a * d * my0 + 2 * a * c * mz0 - 2 * b * d * mz0, 2) + pow(-2 * a * c * mx0 - 2 * b * d * mx0 + 2 * a * b * my0 - 2 * c * d * my0 + Mz - a2 * mz0 + b2 * mz0 + c2 * mz0 - d2 * mz0, 2) + pow(2 * a * d * mx0 + My - a2 * my0 + b2 * my0 - c2 * my0 + d2 * my0 - 2 * c * d * mz0 - 2 * b * (c * mx0 + a * mz0), 2)) / 4.;
};

T4d DDq_norm_f(const T4d &Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0, const T3Tddd &Tmag)
{
    //ゼロにしたい関数のqに関する微分
    auto [Tmag0, Tmag1, Tmag2] = Tmag;
    auto [Tmag00, Tmag01, Tmag02] = Tmag0;
    auto [Tmag10, Tmag11, Tmag12] = Tmag1;
    auto [Tmag20, Tmag21, Tmag22] = Tmag2;
    auto [Ax, Ay, Az] = A;     // measured
    auto [Mx, My, Mz] = M;     // measured
    auto [a, b, c, d] = Q;     // 計算しているクォータニオン
    auto [gx0, gy0, gz0] = G0; //
    auto [mx0, my0, mz0] = M0; //
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    double a2b2c2d2 = a2 + b2 + c2 + d2;
    double normGM2 = gx0 * gx0 + gy0 * gy0 + gz0 * gz0 + mx0 * mx0 + my0 * my0 + mz0 * mz0;
    // return {Az * (-(c * gx0) + b * gy0 - a * gz0) + Ay * (d * gx0 - a * gy0 - b * gz0) + Ax * (-(a * gx0) - d * gy0 + c * gz0) + Mz * (-(c * mx0) + b * my0 - a * mz0) + My * (d * mx0 - a * my0 - b * mz0) + Mx * (-(a * mx0) - d * my0 + c * mz0) + a * a2b2c2d2 * normGM2, Ay * (-(c * gx0) + b * gy0 - a * gz0) + Az * (-(d * gx0) + a * gy0 + b * gz0) + Ax * (-(b * gx0) - c * gy0 - d * gz0) + My * (-(c * mx0) + b * my0 - a * mz0) + Mz * (-(d * mx0) + a * my0 + b * mz0) + Mx * (-(b * mx0) - c * my0 - d * mz0) + a2b2c2d2 * b * normGM2, Ax * (c * gx0 - b * gy0 + a * gz0) + Az * (-(a * gx0) - d * gy0 + c * gz0) + Ay * (-(b * gx0) - c * gy0 - d * gz0) + Mx * (c * mx0 - b * my0 + a * mz0) + Mz * (-(a * mx0) - d * my0 + c * mz0) + My * (-(b * mx0) - c * my0 - d * mz0) + a2b2c2d2 * c * normGM2, Ax * (d * gx0 - a * gy0 - b * gz0) + Ay * (a * gx0 + d * gy0 - c * gz0) + Az * (-(b * gx0) - c * gy0 - d * gz0) + Mx * (d * mx0 - a * my0 - b * mz0) + My * (a * mx0 + d * my0 - c * mz0) + Mz * (-(b * mx0) - c * my0 - d * mz0) + a2b2c2d2 * d * normGM2};

    return {Az * (-(c * gx0) + b * gy0 - a * gz0) + Ay * (d * gx0 - a * gy0 - b * gz0) + Ax * (-(a * gx0) - d * gy0 + c * gz0) + Mz * (-(c * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) + b * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - a * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + My * (d * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - a * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - b * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + Mx * (-(a * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) - d * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) + c * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + a * a2b2c2d2 * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))), Ay * (-(c * gx0) + b * gy0 - a * gz0) + Az * (-(d * gx0) + a * gy0 + b * gz0) + Ax * (-(b * gx0) - c * gy0 - d * gz0) + My * (-(c * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) + b * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - a * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + Mz * (-(d * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) + a * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) + b * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + Mx * (-(b * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) - c * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - d * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + a2b2c2d2 * b * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))), Ax * (c * gx0 - b * gy0 + a * gz0) + Az * (-(a * gx0) - d * gy0 + c * gz0) + Ay * (-(b * gx0) - c * gy0 - d * gz0) + Mx * (c * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - b * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) + a * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + Mz * (-(a * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) - d * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) + c * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + My * (-(b * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) - c * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - d * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + a2b2c2d2 * c * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))), Ax * (d * gx0 - a * gy0 - b * gz0) + Ay * (a * gx0 + d * gy0 - c * gz0) + Az * (-(b * gx0) - c * gy0 - d * gz0) + Mx * (d * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - a * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - b * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + My * (a * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) + d * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - c * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + Mz * (-(b * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02)) - c * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - d * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22)) + a2b2c2d2 * d * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2)))};
};

T4T4d D2D2q_norm_f(const T4d &Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0, const T3Tddd &Tmag)
{
    //ゼロにしたい関数のqに関する微分
    auto [Tmag0, Tmag1, Tmag2] = Tmag;
    auto [Tmag00, Tmag01, Tmag02] = Tmag0;
    auto [Tmag10, Tmag11, Tmag12] = Tmag1;
    auto [Tmag20, Tmag21, Tmag22] = Tmag2;
    //ゼロにしたい関数の微分の微分
    auto [Ax, Ay, Az] = A;     // measured
    auto [Mx, My, Mz] = M;     // measured
    auto [a, b, c, d] = Q;     // 計算しているクォータニオン
    auto [gx0, gy0, gz0] = G0; //
    auto [mx0, my0, mz0] = M0; //
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    double a2b2c2d2 = a2 + b2 + c2 + d2;
    double normGM2 = gx0 * gx0 + gy0 * gy0 + gz0 * gz0 + mx0 * mx0 + my0 * my0 + mz0 * mz0;
    return {{-(Ax * gx0) - Ay * gy0 - Az * gz0 - Mx * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - My * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) -
                 Mz * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) + (3 * a2 + b2 + c2 + d2) * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             Az * gy0 - Ay * gz0 + Mz * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - My * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) +
                 2 * a * b * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             -(Az * gx0) + Ax * gz0 - Mz * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) + Mx * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) +
                 2 * a * c * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             Ay * gx0 - Ax * gy0 + My * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - Mx * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) +
                 2 * a * d * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2)))},
            {Az * gy0 - Ay * gz0 + Mz * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) -
                 My * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) +
                 2 * a * b * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             -(Ax * gx0) + Ay * gy0 + Az * gz0 - Mx * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) + My * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) + Mz * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) + (a2 + 3 * b2 + c2 + d2) * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))), -(Ay * gx0) - Ax * gy0 - My * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - Mx * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) + 2 * b * c * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))), -(Az * gx0) - Ax * gz0 - Mz * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - Mx * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) + 2 * b * d * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2)))},
            {-(Az * gx0) + Ax * gz0 - Mz * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) +
                 Mx * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) +
                 2 * a * c * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             -(Ay * gx0) - Ax * gy0 - My * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - Mx * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) +
                 2 * b * c * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             Ax * gx0 - Ay * gy0 + Az * gz0 + Mx * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - My * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) + Mz * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) + (a2 + b2 + 3 * c2 + d2) * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))), -(Az * gy0) - Ay * gz0 - Mz * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - My * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) + 2 * c * d * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2)))},
            {Ay * gx0 - Ax * gy0 + My * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) -
                 Mx * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) +
                 2 * a * d * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             -(Az * gx0) - Ax * gz0 - Mz * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) - Mx * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) +
                 2 * b * d * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             -(Az * gy0) - Ay * gz0 - Mz * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - My * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) +
                 2 * c * d * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2))),
             Ax * gx0 + Ay * gy0 - Az * gz0 + Mx * (mx0 * Tmag00 + my0 * Tmag01 + mz0 * Tmag02) + My * (mx0 * Tmag10 + my0 * Tmag11 + mz0 * Tmag12) - Mz * (mx0 * Tmag20 + my0 * Tmag21 + mz0 * Tmag22) + (a2 + b2 + c2 + 3 * d2) * (2 * my0 * mz0 * (Tmag01 * Tmag02 + Tmag11 * Tmag12 + Tmag21 * Tmag22) + 2 * mx0 * (my0 * (Tmag00 * Tmag01 + Tmag10 * Tmag11 + Tmag20 * Tmag21) + mz0 * (Tmag00 * Tmag02 + Tmag10 * Tmag12 + Tmag20 * Tmag22)) + pow(gx0, 2) + pow(gy0, 2) + pow(gz0, 2) + pow(mx0, 2) * (pow(Tmag00, 2) + pow(Tmag10, 2) + pow(Tmag20, 2)) + pow(my0, 2) * (pow(Tmag01, 2) + pow(Tmag11, 2) + pow(Tmag21, 2)) + pow(mz0, 2) * (pow(Tmag02, 2) + pow(Tmag12, 2) + pow(Tmag22, 2)))}};
};

T7d DDq_norm_f_extra(const T4d &Q, const Tddd &ABODY, const T4d &approx_Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    //ゼロにしたい関数のqに関する微分
    auto [Ax, Ay, Az] = A;                 // measured
    auto [Mx, My, Mz] = M;                 // measured
    auto [a, b, c, d] = Q;                 // 計算しているクォータニオン
    auto [Aa, Ab, Ac, Ad] = approx_Q;      // 計算しているクォータニオン
    auto [gx0, gy0, gz0] = G0;             //
    auto [mx0, my0, mz0] = M0;             //
    auto [ABODYx, ABODYy, ABODYz] = ABODY; // 計算しているクォータニオン
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    double a3 = a2 * a, b3 = b2 * b, c3 = c2 * c, d3 = d2 * d;
    double a2b2c2d2 = a2 + b2 + c2 + d2;
    double normGM2 = gx0 * gx0 + gy0 * gy0 + gz0 * gz0 + mx0 * mx0 + my0 * my0 + mz0 * mz0;
    return {ABODYz * c * gx0 - ABODYy * d * gx0 - ABODYz * b * gy0 + ABODYx * d * gy0 + ABODYy * b * gz0 - ABODYx * c * gz0 + Az * (-(c * gx0) + b * gy0 - a * gz0) + Ay * (d * gx0 - a * gy0 - b * gz0) + Ax * (-(a * gx0) - d * gy0 + c * gz0) + d * mx0 * My - d * Mx * my0 - c * mx0 * Mz + b * my0 * Mz + c * Mx * mz0 - b * My * mz0 + a3 * normGM2 + a * (ABODYx * gx0 + ABODYy * gy0 + ABODYz * gz0 - Mx * mx0 - My * my0 - Mz * mz0 + (b2 + c2 + d2) * normGM2), ABODYz * d * gx0 + a2 * b * pow(gx0, 2) + b3 * pow(gx0, 2) + b * c2 * pow(gx0, 2) + b * d2 * pow(gx0, 2) - a * ABODYz * gy0 + a2 * b * pow(gy0, 2) + b3 * pow(gy0, 2) + b * c2 * pow(gy0, 2) + b * d2 * pow(gy0, 2) - ABODYz * b * gz0 + a2 * b * pow(gz0, 2) + b3 * pow(gz0, 2) + b * c2 * pow(gz0, 2) + b * d2 * pow(gz0, 2) + Ay * (-(c * gx0) + b * gy0 - a * gz0) + ABODYy * (c * gx0 - b * gy0 + a * gz0) + Az * (-(d * gx0) + a * gy0 + b * gz0) + Ax * (-(b * gx0) - c * gy0 - d * gz0) + ABODYx * (b * gx0 + c * gy0 + d * gz0) - b * Mx * mx0 + a2 * b * pow(mx0, 2) + b3 * pow(mx0, 2) + b * c2 * pow(mx0, 2) + b * d2 * pow(mx0, 2) - c * mx0 * My - c * Mx * my0 + b * My * my0 + a2 * b * pow(my0, 2) + b3 * pow(my0, 2) + b * c2 * pow(my0, 2) + b * d2 * pow(my0, 2) - d * mx0 * Mz + a * my0 * Mz - (d * Mx + a * My - b * Mz) * mz0 + a2b2c2d2 * b * pow(mz0, 2), -(ABODYx * c * gx0) + b2 * c * pow(gx0, 2) + c2 * pow(gx0, 2) + c * d2 * pow(gx0, 2) + ABODYx * b * gy0 + ABODYz * d * gy0 + b2 * c * pow(gy0, 2) + c2 * pow(gy0, 2) + c * d2 * pow(gy0, 2) - ABODYz * c * gz0 + b2 * c * pow(gz0, 2) + c2 * pow(gz0, 2) + c * d2 * pow(gz0, 2) + Ax * (c * gx0 - b * gy0 + a * gz0) + Az * (-(a * gx0) - d * gy0 + c * gz0) + Ay * (-(b * gx0) - c * gy0 - d * gz0) + ABODYy * (b * gx0 + c * gy0 + d * gz0) + c * Mx * mx0 + b2 * c * pow(mx0, 2) + c2 * pow(mx0, 2) + c * d2 * pow(mx0, 2) - b * mx0 * My - b * Mx * my0 - c * My * my0 + b2 * c * pow(my0, 2) + c2 * pow(my0, 2) + c * d2 * pow(my0, 2) - d * my0 * Mz - d * My * mz0 + c * Mz * mz0 + c * (b2 + c2 + d2) * pow(mz0, 2) + a * (ABODYz * gx0 - ABODYx * gz0 - mx0 * Mz + Mx * mz0) + a2 * c * normGM2, -(ABODYx * d * gx0) + b2 * d * pow(gx0, 2) + c2 * d * pow(gx0, 2) + d3 * pow(gx0, 2) - ABODYy * d * gy0 + b2 * d * pow(gy0, 2) + c2 * d * pow(gy0, 2) + d3 * pow(gy0, 2) + ABODYx * b * gz0 + ABODYy * c * gz0 + b2 * d * pow(gz0, 2) + c2 * d * pow(gz0, 2) + d3 * pow(gz0, 2) + Ax * (d * gx0 - a * gy0 - b * gz0) + Ay * (a * gx0 + d * gy0 - c * gz0) + Az * (-(b * gx0) - c * gy0 - d * gz0) + ABODYz * (b * gx0 + c * gy0 + d * gz0) + d * Mx * mx0 + b2 * d * pow(mx0, 2) + c2 * d * pow(mx0, 2) + d3 * pow(mx0, 2) + d * My * my0 + b2 * d * pow(my0, 2) + c2 * d * pow(my0, 2) + d3 * pow(my0, 2) + a * (-(ABODYy * gx0) + ABODYx * gy0 + mx0 * My - Mx * my0) - b * mx0 * Mz - c * my0 * Mz - (b * Mx + c * My + d * Mz) * mz0 + d * (b2 + c2 + d2) * pow(mz0, 2) + a2 * d * normGM2, -Ax + (2 * ABODYx + a2 * gx0 + pow(Aa, 2) * gx0 + pow(Ab, 2) * gx0 - pow(Ac, 2) * gx0 - pow(Ad, 2) * gx0 + b2 * gx0 - c2 * gx0 - d2 * gx0 + 2 * Ab * Ac * gy0 + 2 * Aa * Ad * gy0 + 2 * b * c * gy0 + 2 * a * d * gy0 + 2 * (-(Aa * Ac) + Ab * Ad - a * c + b * d) * gz0) / 2., ABODYy - Ay + (Ab * Ac - Aa * Ad + b * c - a * d) * gx0 + ((a2 + pow(Aa, 2) - pow(Ab, 2) + pow(Ac, 2) - pow(Ad, 2) - b2 + c2 - d2) * gy0) / 2. + (Aa * Ab + Ac * Ad + a * b + c * d) * gz0, ABODYz - Az + Aa * Ac * gx0 + Ab * Ad * gx0 + a * c * gx0 + b * d * gx0 - Aa * Ab * gy0 + Ac * Ad * gy0 - a * b * gy0 + c * d * gy0 + ((a2 + pow(Aa, 2) - pow(Ab, 2) - pow(Ac, 2) + pow(Ad, 2) - b2 - c2 + d2) * gz0) / 2.};
};

T7T7d D2D2q_norm_f_extra(const T4d &Q, const Tddd &ABODY, const T4d &approx_Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    //ゼロにしたい関数の微分の微分
    auto [Ax, Ay, Az] = A;                 // measured
    auto [Mx, My, Mz] = M;                 // measured
    auto [a, b, c, d] = Q;                 // 計算しているクォータニオン
    auto [Aa, Ab, Ac, Ad] = approx_Q;      // 計算しているクォータニオン
    auto [gx0, gy0, gz0] = G0;             //
    auto [mx0, my0, mz0] = M0;             //
    auto [ABODYx, ABODYy, ABODYz] = ABODY; // 計算しているクォータニオン
    double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
    double a3 = a2 * a, b3 = b2 * b, c3 = c2 * c, d3 = d2 * d;
    double a2b2c2d2 = a2 + b2 + c2 + d2;
    double normGM2 = gx0 * gx0 + gy0 * gy0 + gz0 * gz0 + mx0 * mx0 + my0 * my0 + mz0 * mz0;
    return {{ABODYx * gx0 - Ax * gx0 + ABODYy * gy0 - Ay * gy0 + ABODYz * gz0 - Az * gz0 - Mx * mx0 - My * my0 - Mz * mz0 + (3 * a2 + b2 + c2 + d2) * normGM2,
             -(ABODYz * gy0) + Az * gy0 + ABODYy * gz0 - Ay * gz0 + my0 * Mz - My * mz0 + 2 * a * b * normGM2,
             ABODYz * gx0 - Az * gx0 - ABODYx * gz0 + Ax * gz0 - mx0 * Mz + Mx * mz0 + 2 * a * c * normGM2,
             -(ABODYy * gx0) + Ay * gx0 + ABODYx * gy0 - Ax * gy0 + mx0 * My - Mx * my0 + 2 * a * d * normGM2,
             a * gx0 + d * gy0 - c * gz0,
             -(d * gx0) + a * gy0 + b * gz0,
             c * gx0 - b * gy0 + a * gz0},
            {-(ABODYz * gy0) + Az * gy0 + ABODYy * gz0 - Ay * gz0 + my0 * Mz - My * mz0 + 2 * a * b * normGM2,
             ABODYx * gx0 - Ax * gx0 - ABODYy * gy0 + Ay * gy0 - ABODYz * gz0 + Az * gz0 - Mx * mx0 + My * my0 + Mz * mz0 + (a2 + 3 * b2 + c2 + d2) * normGM2,
             ABODYy * gx0 - Ay * gx0 + ABODYx * gy0 - Ax * gy0 - mx0 * My - Mx * my0 + 2 * b * c * normGM2,
             ABODYz * gx0 - Az * gx0 + ABODYx * gz0 - Ax * gz0 - mx0 * Mz - Mx * mz0 + 2 * b * d * normGM2,
             b * gx0 + c * gy0 + d * gz0,
             c * gx0 - b * gy0 + a * gz0,
             d * gx0 - a * gy0 - b * gz0},
            {ABODYz * gx0 - Az * gx0 - ABODYx * gz0 + Ax * gz0 - mx0 * Mz + Mx * mz0 + 2 * a * c * normGM2,
             ABODYy * gx0 - Ay * gx0 + ABODYx * gy0 - Ax * gy0 - mx0 * My - Mx * my0 + 2 * b * c * normGM2,
             -(ABODYx * gx0) + Ax * gx0 + ABODYy * gy0 - Ay * gy0 - ABODYz * gz0 + Az * gz0 + Mx * mx0 - My * my0 + Mz * mz0 + (a2 + b2 + 3 * c2 + d2) * normGM2,
             ABODYz * gy0 - Az * gy0 + ABODYy * gz0 - Ay * gz0 - my0 * Mz - My * mz0 + 2 * c * d * normGM2,
             -(c * gx0) + b * gy0 - a * gz0,
             b * gx0 + c * gy0 + d * gz0,
             a * gx0 + d * gy0 - c * gz0},
            {-(ABODYy * gx0) + Ay * gx0 + ABODYx * gy0 - Ax * gy0 + mx0 * My - Mx * my0 + 2 * a * d * normGM2,
             ABODYz * gx0 - Az * gx0 + ABODYx * gz0 - Ax * gz0 - mx0 * Mz - Mx * mz0 + 2 * b * d * normGM2,
             ABODYz * gy0 - Az * gy0 + ABODYy * gz0 - Ay * gz0 - my0 * Mz - My * mz0 + 2 * c * d * normGM2,
             -(ABODYx * gx0) + Ax * gx0 - ABODYy * gy0 + Ay * gy0 + ABODYz * gz0 - Az * gz0 + Mx * mx0 + My * my0 - Mz * mz0 + (a2 + b2 + c2 + 3 * d2) * normGM2,
             -(d * gx0) + a * gy0 + b * gz0,
             -(a * gx0) - d * gy0 + c * gz0,
             b * gx0 + c * gy0 + d * gz0},
            {a * gx0 + d * gy0 - c * gz0, b * gx0 + c * gy0 + d * gz0, -(c * gx0) + b * gy0 - a * gz0, -(d * gx0) + a * gy0 + b * gz0,
             1,
             0,
             0},
            {-(d * gx0) + a * gy0 + b * gz0, c * gx0 - b * gy0 + a * gz0, b * gx0 + c * gy0 + d * gz0, -(a * gx0) - d * gy0 + c * gz0,
             0,
             1,
             0},
            {c * gx0 - b * gy0 + a * gz0,
             d * gx0 - a * gy0 - b * gz0,
             a * gx0 + d * gy0 - c * gz0,
             b * gx0 + c * gy0 + d * gz0,
             0,
             0,
             1}};
};

T7d DDq_norm_f_extra(const T7d &Q_ABODY, const T4d &approx_Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    auto [q0, q1, q2, q3, ABODYx, ABODYy, ABODYz] = Q_ABODY;
    return DDq_norm_f_extra({q0, q1, q2, q3}, {ABODYx, ABODYy, ABODYz}, approx_Q, A, M, G0, M0);
};
T7T7d D2D2q_norm_f_extra(const T7d &Q_ABODY, const T4d &approx_Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
{
    auto [q0, q1, q2, q3, ABODYx, ABODYy, ABODYz] = Q_ABODY;
    return D2D2q_norm_f_extra({q0, q1, q2, q3}, {ABODYx, ABODYy, ABODYz}, approx_Q, A, M, G0, M0);
};
struct Fusion
{
    std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd> time_Q0;
    std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd> time_Q1;
    std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd> time_Q2;
    std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd> time_Q3;
    Tddd G0;
    Tddd M0;
    /* ---------------------- 磁場の構成パラメタ --------------------- */
    Tdd minmax_Mx;
    Tdd minmax_My;
    Tdd minmax_Mz;
    bool activateOffsetM, activateScaleM;
    Tddd offsetM, scaleM;
    T3Tddd magTransMat;
    /* ------------------------------------------------------ */
    NewtonRaphson<T4d> NR;
    NewtonRaphson<T7d> NR_T7d;
    Fusion(const Tddd &G0_IN, const Tddd &M0_IN)
        : G0(G0_IN), M0(M0_IN),
          NR({1., 0., 0., 0.}),
          NR_T7d({1., 0., 0., 0., 0., 0., 0.}),
          minmax_Mx({1E+50, -1E+50}),
          minmax_My({1E+50, -1E+50}),
          minmax_Mz({1E+50, -1E+50}),
          offsetM({0., 0., 0.}),
          activateOffsetM(false),
          activateScaleM(false),
          magTransMat({1, 0, 0}, {0, 1, 0}, {0, 0, 1}),
          time_Q0({0., Quaternion(), {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}),
          time_Q1({1., Quaternion(), {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}),
          time_Q2({2., Quaternion(), {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}),
          time_Q3({3., Quaternion(), {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}){};
    // gxS, gx, u2, u3, u4, u5はセンサーが実際に計測した値 重力と磁場
    void setG(const Tddd &G0_IN) { this->G0 = G0_IN; };
    void setM(const Tddd &M0_IN) { this->M0 = M0_IN; };

    void setScaleM(const Tddd &scaleM_IN)
    {
        this->scaleM = scaleM_IN;
        this->activateScaleM = true;
    };
    void setOffsetM(const Tddd &offsetM_IN)
    {
        this->offsetM = offsetM_IN;
        this->activateOffsetM = true;
    };
    Tddd calculateOffsetM() const
    {
        return {(std::get<0>(minmax_Mx) + std::get<1>(minmax_Mx)) / 2.,
                (std::get<0>(minmax_My) + std::get<1>(minmax_My)) / 2.,
                (std::get<0>(minmax_Mz) + std::get<1>(minmax_Mz)) / 2.};
    };

    T6d original_eq(const T4d &Q, const Tddd &A, const Tddd &M, const Tddd &G0, const Tddd &M0)
    {
        auto [Ax, Ay, Az] = A;     // measured
        auto [Mx, My, Mz] = M;     // measured
        auto [gx0, gy0, gz0] = G0; //
        auto [mx0, my0, mz0] = M0; //
        auto [a, b, c, d] = Q;     // 計算しているクォータニオン
        double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
        return {pow(Ax - a2 * gx0 - b2 * gx0 + (c2 + d2) * gx0 + a * (-2 * d * gy0 + 2 * c * gz0) - 2 * b * (c * gy0 + d * gz0), 2),
                pow(Ay + 2 * a * d * gx0 - a2 * gy0 + b2 * gy0 - c2 * gy0 + d2 * gy0 - 2 * c * d * gz0 - 2 * b * (c * gx0 + a * gz0), 2),
                pow(Az - 2 * (a * c * gx0 + b * d * gx0 - a * b * gy0 + c * d * gy0) + (-a2 + b2 + c2 - d2) * gz0, 2),
                pow(Mx - a2 * mx0 - b2 * mx0 + (c2 + d2) * mx0 + a * (-2 * d * my0 + 2 * c * mz0) - 2 * b * (c * my0 + d * mz0), 2),
                pow(2 * a * d * mx0 + My - a2 * my0 + b2 * my0 - c2 * my0 + d2 * my0 - 2 * c * d * mz0 - 2 * b * (c * mx0 + a * mz0), 2),
                pow(-2 * (a * c * mx0 + b * d * mx0 - a * b * my0 + c * d * my0) + Mz + (-a2 + b2 + c2 - d2) * mz0, 2)};
    };

    void setMagTransMat(const T3Tddd &transMat)
    {
        magTransMat = transMat;
    };

    Quaternion updateStandard(const Tddd &A, const Tddd &M_IN, const Tddd &W, const double t, const double alpha)
    {
        /*
            alphaは，加速度と磁気から予測した姿勢の重み
            (1-alpha)は，ジャイロと１ステップ前の姿勢から予測された現在の姿勢の重み
        */
        Tddd M = M_IN; //現在の磁気ベクトル
        Tddd M0_ = M0; //初期の磁気ベクトル
        if (this->activateOffsetM)
        {
            // 磁気センサーの補正
            M -= this->offsetM;
            M0_ -= this->offsetM;
        }

        if (true)
        {
            // 磁気センサーの補正
            M = Normalize(M);
            M0_ = Normalize(M0_);
        }

        auto dt = std::abs(t - std::get<0>(time_Q0));
        auto W0 = std::get<4>(time_Q0);
        auto Q0 = std::get<1>(time_Q0);
        auto Q1 = std::get<1>(time_Q1);
        auto Q2 = std::get<1>(time_Q2);
        auto Q3 = std::get<1>(time_Q3);

        // １つ過去の結果より現在のクォータニオンを予測
        Quaternion approx_Q(Q0 + Q0.d_dt((W + W0) / 2. * dt));
        // Quaternion approx_Q(Q0);

        NR.X = Normalize(Q0() + Q1() + Q2() + Q3());
        for (auto i = 0; i < 40; ++i)
        {
            NR.update(DDq_norm_f(NR.X, A, M, G0, M0_, this->magTransMat),
                      D2D2q_norm_f(NR.X, A, M, G0, M0_, this->magTransMat)); /*initial X is updated*/
            NR.X = Normalize(NR.X);
            if (Norm(NR.dX) < 1E-6)
                break;
        }

        Quaternion ans(alpha * NR.X + (1. - alpha) * approx_Q());

        time_Q3 = time_Q2;
        time_Q2 = time_Q1;
        time_Q1 = time_Q0;
        std::get<0>(time_Q0) = t;
        std::get<1>(time_Q0).set(ans);
        std::get<2>(time_Q0) = A;
        std::get<3>(time_Q0) = M;
        std::get<4>(time_Q0) = W;
        std::get<5>(time_Q0) = A - ans.Rs(G0);
        return ans;
    };
    Quaternion solveForQuaternion(const Tddd &A, const Tddd &M_IN, const Tddd &W, const double t = 0.)
    {
        /* ------------------------------------------------------ */
        // 磁場のオフセットの計算用．この結果を使うには，calculateOffsetMを実行し，setOffsetMする必要がある．
        auto [Mx, My, Mz] = M_IN;
        if (Mx < std::get<0>(minmax_Mx))
            std::get<0>(minmax_Mx) = Mx;
        else if (std::get<1>(minmax_Mx) < Mx)
            std::get<1>(minmax_Mx) = Mx;

        if (My < std::get<0>(minmax_My))
            std::get<0>(minmax_My) = My;
        else if (std::get<1>(minmax_My) < My)
            std::get<1>(minmax_My) = My;

        if (Mz < std::get<0>(minmax_Mz))
            std::get<0>(minmax_Mz) = Mz;
        else if (std::get<1>(minmax_Mz) < Mz)
            std::get<1>(minmax_Mz) = Mz;
        /* ------------------------------------------------------ */
        Tddd M = M_IN;
        Tddd M0_ = M0;
        // オフセットが有効になった場合のみノーマライズする
        if (this->activateOffsetM)
        {
            M -= this->offsetM;
            M0_ -= this->offsetM;
            // M = Normalize(M);
            // M0_ = Normalize(M0_);
        }
        // if (this->activateScaleM)
        // {
        //     M *= this->scaleM;
        //     M0_ *= this->scaleM;
        // }
        /* ------------------------------------------------------ */
        auto dt = std::abs(t - std::get<0>(time_Q0));
        // Magnetについては，Normalizeした方がいいかもしれない．
        //  NRは前回の結果を引き継いでいる
        auto W2 = std::get<4>(time_Q2);
        auto W1 = std::get<4>(time_Q1);
        auto W0 = std::get<4>(time_Q0);
        auto Q0 = std::get<1>(time_Q0)();
        auto Q1 = std::get<1>(time_Q1)();
        auto Q2 = std::get<1>(time_Q2)();
        auto approx_Q = (Q0 + Q1) / 2. + std::get<1>(time_Q1).d_dt((W + W0 + W1 + W2) / 4. * dt)();

        auto Ab0 = std::get<5>(time_Q0);
        auto Ab1 = std::get<5>(time_Q1);

        auto A_ = A;
        NR.X = Q0;
        for (auto i = 0; i < 50; ++i)
        {
            NR.update(DDq_norm_f(NR.X, A_, M, G0, M0_, this->magTransMat),
                      D2D2q_norm_f(NR.X, A_, M, G0, M0_, this->magTransMat)); /*initial X is updated*/
            NR.X = Normalize(NR.X);
            if (Norm(NR.dX) < 1E-9)
                break;
        }

        Quaternion ans;
        if (Norm(W - W0) > 5. /*角加速が大きすぎる*/)
            ans.set(Quaternion(0.9 * NR.X + 0.1 * Q0));
        else
            ans.set(Quaternion(0.9 * NR.X + 0.1 * approx_Q));

        time_Q3 = time_Q2;
        time_Q2 = time_Q1;
        time_Q1 = time_Q0;
        std::get<0>(time_Q0) = t;
        std::get<1>(time_Q0).set(ans);
        std::get<2>(time_Q0) = A;
        std::get<3>(time_Q0) = M;
        std::get<4>(time_Q0) = W;
        std::get<5>(time_Q0) = A - ans.Rs(G0);
        return std::get<1>(time_Q0);
    };

    // Quaternion solveForQuaternion(const Tddd &A, const Tddd &M_IN, const Tddd &W, const double t = 0.)
    // {
    //     /* ------------------------------------------------------ */
    //     // 磁場のオフセットの計算用．この結果を使うには，calculateOffsetMを実行し，setOffsetMする必要がある．
    //     auto [Mx, My, Mz] = M_IN;
    //     if (Mx < std::get<0>(minmax_Mx))
    //         std::get<0>(minmax_Mx) = Mx;
    //     else if (std::get<1>(minmax_Mx) < Mx)
    //         std::get<1>(minmax_Mx) = Mx;

    //     if (My < std::get<0>(minmax_My))
    //         std::get<0>(minmax_My) = My;
    //     else if (std::get<1>(minmax_My) < My)
    //         std::get<1>(minmax_My) = My;

    //     if (Mz < std::get<0>(minmax_Mz))
    //         std::get<0>(minmax_Mz) = Mz;
    //     else if (std::get<1>(minmax_Mz) < Mz)
    //         std::get<1>(minmax_Mz) = Mz;
    //     /* ------------------------------------------------------ */
    //     Tddd M = M_IN;
    //     Tddd M0_ = M0;
    //     // オフセットが有効になった場合のみノーマライズする
    //     if (this->activateOffsetM)
    //     {
    //         M -= this->offsetM;
    //         M0_ -= this->offsetM;
    //         M = Normalize(M);
    //         M0_ = Normalize(M0_);
    //     }
    //     // if (this->activateScaleM)
    //     // {
    //     //     M *= this->scaleM;
    //     //     M0_ *= this->scaleM;
    //     // }
    //     /* ------------------------------------------------------ */
    //     auto dt = std::abs(t - std::get<0>(time_Q0));
    //     // Magnetについては，Normalizeした方がいいかもしれない．
    //     //  NRは前回の結果を引き継いでいる
    //     auto W2 = std::get<4>(time_Q2);
    //     auto W1 = std::get<4>(time_Q1);
    //     auto W0 = std::get<4>(time_Q0);
    //     auto Q0 = std::get<1>(time_Q0)();
    //     auto Q1 = std::get<1>(time_Q1)();
    //     auto Q2 = std::get<1>(time_Q2)();
    //     auto approx_Q = (Q0 + Q1) / 2. + std::get<1>(time_Q1).d_dt((W + W0 + W1 + W2) / 4. * dt)();

    //     auto Ab0 = std::get<5>(time_Q0);
    //     auto Ab1 = std::get<5>(time_Q1);

    //     auto A_ = A;
    //     NR.X = Q0;
    //     for (auto i = 0; i < 50; ++i)
    //     {
    //         NR.update(DDq_norm_f(NR.X, A_, M, G0, M0_, this->magTransMat),
    //                   D2D2q_norm_f(NR.X, A_, M, G0, M0_, this->magTransMat)); /*initial X is updated*/
    //         NR.X = Normalize(NR.X);
    //         if (Norm(NR.dX) < 1E-9)
    //             break;
    //     }

    //     Quaternion ans;
    //     if (Norm(W - W0) > 5. /*角加速が大きすぎる*/)
    //         ans.set(Quaternion(0.9 * NR.X + 0.1 * Q0));
    //     else
    //         ans.set(Quaternion(0.9 * NR.X + 0.1 * approx_Q));

    //     time_Q3 = time_Q2;
    //     time_Q2 = time_Q1;
    //     time_Q1 = time_Q0;
    //     std::get<0>(time_Q0) = t;
    //     std::get<1>(time_Q0).set(ans);
    //     std::get<2>(time_Q0) = A;
    //     std::get<3>(time_Q0) = M;
    //     std::get<4>(time_Q0) = W;
    //     std::get<5>(time_Q0) = A - ans.Rs(G0);
    //     return std::get<1>(time_Q0);
    // };

    std::vector<std::tuple<double, Quaternion, Tddd, Tddd, Tddd, Tddd>> history() const
    {
        return {time_Q0, time_Q1, time_Q2, time_Q3};
    };

    std::vector<std::tuple<double, Quaternion>> historyQ() const
    {
        return {{std::get<0>(time_Q0), std::get<1>(time_Q0)},
                {std::get<0>(time_Q1), std::get<1>(time_Q1)},
                {std::get<0>(time_Q2), std::get<1>(time_Q2)},
                {std::get<0>(time_Q3), std::get<1>(time_Q3)}};
    };
    std::vector<std::tuple<double, T4d>> historyQ_tuple() const
    {
        return {{std::get<0>(time_Q0), std::get<1>(time_Q0)()},
                {std::get<0>(time_Q1), std::get<1>(time_Q1)()},
                {std::get<0>(time_Q2), std::get<1>(time_Q2)()},
                {std::get<0>(time_Q3), std::get<1>(time_Q3)()}};
    };

    std::vector<std::tuple<double, Tddd>> historyA() const
    {
        return {{std::get<0>(time_Q0), std::get<2>(time_Q0)},
                {std::get<0>(time_Q1), std::get<2>(time_Q1)},
                {std::get<0>(time_Q2), std::get<2>(time_Q2)},
                {std::get<0>(time_Q3), std::get<2>(time_Q3)}};
    };

    std::vector<std::tuple<double, Tddd>> historyM() const
    {
        return {{std::get<0>(time_Q0), std::get<3>(time_Q0)},
                {std::get<0>(time_Q1), std::get<3>(time_Q1)},
                {std::get<0>(time_Q2), std::get<3>(time_Q2)},
                {std::get<0>(time_Q3), std::get<3>(time_Q3)}};
    };

    std::vector<std::tuple<double, Tddd>> historyW() const
    {
        return {{std::get<0>(time_Q0), std::get<4>(time_Q0)},
                {std::get<0>(time_Q1), std::get<4>(time_Q1)},
                {std::get<0>(time_Q2), std::get<4>(time_Q2)},
                {std::get<0>(time_Q3), std::get<4>(time_Q3)}};
    };

    std::vector<std::tuple<double, Tddd>> historyAbody() const
    {
        return {{std::get<0>(time_Q0), std::get<5>(time_Q0)},
                {std::get<0>(time_Q1), std::get<5>(time_Q1)},
                {std::get<0>(time_Q2), std::get<5>(time_Q2)},
                {std::get<0>(time_Q3), std::get<5>(time_Q3)}};
    };

    std::tuple<Quaternion, Tddd, Tddd, Tddd, Tddd> interp(const double t) const
    {
        return {Quaternion(InterpolationLagrange(historyQ_tuple())(t)),
                InterpolationLagrange(historyA())(t),
                InterpolationLagrange(historyM())(t),
                InterpolationLagrange(historyW())(t),
                InterpolationLagrange(historyAbody())(t)};
    };

    Quaternion interpQ(const double t) const
    {
        return Quaternion(InterpolationLagrange(historyQ_tuple())(t));
    };

    T4d interpQtuple(const double t) const
    {
        return InterpolationLagrange(historyQ_tuple())(t);
    };

    Tddd interpA(const double t) const
    {
        return InterpolationLagrange(historyA())(t);
    };

    Tddd interpM(const double t) const
    {
        return InterpolationLagrange(historyM())(t);
    };

    Tddd interpW(const double t) const
    {
        return InterpolationLagrange(historyW())(t);
    };

    Tddd interpAbody(const double t) const
    {
        return InterpolationLagrange(historyAbody())(t);
    };

    Tddd interpYPR(const double t) const
    {
        return Quaternion(InterpolationLagrange(historyQ_tuple())(t)).YPR();
    };

    Quaternion solveForQuaternionModified(const Tddd &A, const Tddd &M_IN, const Tddd &W, const double t = 0.)
    {
        /* ------------------------------------------------------ */
        // 磁場のオフセットの計算用．この結果を使うには，calculateOffsetMを実行し，setOffsetMする必要がある．
        auto [Mx, My, Mz] = M_IN;
        if (Mx < std::get<0>(minmax_Mx))
            std::get<0>(minmax_Mx) = Mx;
        else if (std::get<1>(minmax_Mx) < Mx)
            std::get<1>(minmax_Mx) = Mx;

        if (My < std::get<0>(minmax_My))
            std::get<0>(minmax_My) = My;
        else if (std::get<1>(minmax_My) < My)
            std::get<1>(minmax_My) = My;

        if (Mz < std::get<0>(minmax_Mz))
            std::get<0>(minmax_Mz) = Mz;
        else if (std::get<1>(minmax_Mz) < Mz)
            std::get<1>(minmax_Mz) = Mz;
        /* ------------------------------------------------------ */
        Tddd M = M_IN;
        Tddd M0_ = M0;
        // オフセットが有効になった場合のみノーマライズする
        if (this->activateOffsetM)
        {
            M -= this->offsetM;
            M0_ -= this->offsetM;
        }
        // Magnetについては，Normalizeした方がいいかもしれない．
        //  NRは前回の結果を引き継いでいる
        auto dt = std::abs(t - std::get<0>(time_Q0));
        auto ApproxQ = (std::get<0>(this->interp(t)))();
        if (1E-10 < dt && dt < 0.06)
        {
            // Magnetについては，Normalizeした方がいいかもしれない．
            //  NRは前回の結果を引き継いでいる
            auto W_ = (W + std::get<4>(time_Q0)) / 2.;
            auto Q = Quaternion(std::get<1>(time_Q0) + std::get<1>(time_Q0).d_dt(W_ * dt));
            ApproxQ = Q();
            // NR.X = Q();
        };

        auto init = NR_T7d.X;
        int i = 0;
        while (i++ < 20)
        {
            NR_T7d.update(DDq_norm_f_extra(NR_T7d.X, ApproxQ, A, M, G0, M0_),
                          D2D2q_norm_f_extra(NR_T7d.X, ApproxQ, A, M, G0, M0_)); /*initial X is updated*/
            NR_T7d.X = Normalize(NR_T7d.X);
            if (Norm(NR_T7d.dX) < 1E-4)
                break;
        }
        if (i <= 20)
        {
            time_Q3 = time_Q2;
            time_Q2 = time_Q1;
            time_Q1 = time_Q0;
            std::get<0>(time_Q0) = t;
            std::get<1>(time_Q0).set(T4d{std::get<0>(NR_T7d.X), std::get<1>(NR_T7d.X), std::get<2>(NR_T7d.X), std::get<3>(NR_T7d.X)});
            std::get<2>(time_Q0) = A;
            std::get<3>(time_Q0) = M;
            std::get<4>(time_Q0) = W;
            std::get<5>(time_Q0) = {std::get<4>(NR_T7d.X), std::get<5>(NR_T7d.X), std::get<6>(NR_T7d.X)};
        }
        else
        {
            NR_T7d.X = init;
        }
        return std::get<1>(time_Q0);
    };
};