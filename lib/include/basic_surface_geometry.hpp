#ifndef basic_surface_geometry_H
#define basic_surface_geometry_H

#include "basic_vectors.hpp"
#include "basic_linear_systems.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

// グリッド探索の初期値から、曲面上の真の最近点パラメタを Newton 法で精緻化
// 勾配 g = {(X-P)·∂X/∂t0, (X-P)·∂X/∂t1} = 0 を解く
// 2x2 Hessian をクラメルで直接解く（スレッドセーフ、LAPACK 不要）
template <typename XFunc, typename DXFunc>
inline Tdd refineNearestParam(const Tddd& P, Tdd param,
                              XFunc X_func, DXFunc DX_func,
                              int max_iter = 5, double tol = 1e-10) {
  auto [t0, t1] = param;
  for (int iter = 0; iter < max_iter; ++iter) {
    Tddd X = X_func(t0, t1);
    Tddd diff = X - P;

    Tddd dXdt0 = DX_func(t0, t1, 1, 0);
    Tddd dXdt1 = DX_func(t0, t1, 0, 1);

    double g0 = Dot(diff, dXdt0);
    double g1 = Dot(diff, dXdt1);

    if (std::abs(g0) < tol && std::abs(g1) < tol)
      break;

    Tddd d2Xdt0dt0 = DX_func(t0, t1, 2, 0);
    Tddd d2Xdt1dt1 = DX_func(t0, t1, 0, 2);
    Tddd d2Xdt0dt1 = DX_func(t0, t1, 1, 1);

    double H00 = Dot(dXdt0, dXdt0) + Dot(diff, d2Xdt0dt0);
    double H11 = Dot(dXdt1, dXdt1) + Dot(diff, d2Xdt1dt1);
    double H01 = Dot(dXdt0, dXdt1) + Dot(diff, d2Xdt0dt1);

    double det = H00 * H11 - H01 * H01;
    if (std::abs(det) < 1e-20)
      break;
    double dt0 = -(H11 * g0 - H01 * g1) / det;
    double dt1 = -(H00 * g1 - H01 * g0) / det;

    t0 += dt0;
    t1 += dt1;

    // 参照三角形内にクランプ: t0 >= 0, t1 >= 0, t0+t1 <= 1
    t0 = std::max(0.0, t0);
    t1 = std::max(0.0, t1);
    if (t0 + t1 > 1.0) {
      double s = t0 + t1;
      t0 /= s;
      t1 /= s;
    }
  }
  return {t0, t1};
}


// ============================================================================
// 曲面幾何に関する基本関数群
//
// ここには Network や networkPoint などに依存しない純粋関数だけを置く。
// 入力は生の座標列 Tddd を前提とし，メッシュ固有のラッパは別の場所で持つ。
//
// このファイル内の約束:
//   - cotangentWeight は (cot α + cot β) / 2 を返す
//     つまり half-sum convention を採用する
//   - 2次フィットから曲率を出すときは，
//       S = I^{-1} II
//     を用いる
//     ここで I は第一基本形式 [[E,F],[F,G]]，
//     II は第二基本形式 [[L,M],[M,N]] である
//   - したがって，Monge patch
//       w(u,v) = au^2 + buv + cv^2 + du + ev
//     における傾き d,e も補正に入る
// ============================================================================

namespace surface_geometry {

// ---------------------------------------------------------------------------
// 局所座標系
// ---------------------------------------------------------------------------

struct LocalFrame {
  Tddd u_axis = {1, 0, 0};
  Tddd v_axis = {0, 1, 0};
  Tddd n_axis = {0, 0, 1};
  bool valid = false;
};

// ---------------------------------------------------------------------------
// 2次曲面フィットの結果
// ---------------------------------------------------------------------------

struct QuadricFitResult {
  double a = 0, b = 0, c = 0, d = 0, e = 0; // w = au^2 + buv + cv^2 + du + ev
  LocalFrame frame;                          // フィット時に使った局所座標系
  bool valid = false;
};

// ---------------------------------------------------------------------------
// 主曲率・主方向
// ---------------------------------------------------------------------------

struct PrincipalCurvatureResult {
  double k1 = 0, k2 = 0; // 主曲率．|k1| >= |k2| となるように並べる
  Tddd PD1 = {0, 0, 0};  // k1 に対応する主方向（3D，接平面内）
  Tddd PD2 = {0, 0, 0};  // k2 に対応する主方向（3D，接平面内）
  double kmax = 0;       // max(|k1|, |k2|)
  bool valid = false;
};

// ---------------------------------------------------------------------------
// 局所 2次曲面フィット
// ---------------------------------------------------------------------------
//
// 評価点 origin のまわりで，法線方向を高さとする局所グラフ
//
//   w(u,v) = au^2 + buv + cv^2 + du + ev
//
// を最小二乗で当てはめる。
//
// origin:
//   曲率を評価したい点．局所座標の原点になる
// u_axis, v_axis, n_axis:
//   局所直交座標系の軸（正規化済み）
// neighbor_positions:
//   近傍点群．未知数が 5 個なので最低 5 点必要

inline QuadricFitResult fitQuadricLocal(const Tddd& origin,
                                        const Tddd& u_axis, const Tddd& v_axis, const Tddd& n_axis,
                                        const std::vector<Tddd>& neighbor_positions) {
  QuadricFitResult result;
  const int N = static_cast<int>(neighbor_positions.size());
  if (N < 5)
    return result;

  result.frame = {u_axis, v_axis, n_axis, true};

  // 正規方程式 A^T A x = A^T b に直接累積する（5×5 固定、ヒープ割当無し、スレッド安全）
  // A は各行が [u², uv, v², u, v]、未知数は 5 個。行列 A を明示的に保持する必要はなく、
  // 1 点処理するたびに対応する row を stack 上の array に作って A^T A / A^T b に加算すればよい。
  // Cholesky が下三角しか参照しないので、対称行列は下三角 (k ≤ j) だけ埋める。
  std::array<std::array<double, 5>, 5> ATA = {};
  std::array<double, 5> ATb = {};
  for (int i = 0; i < N; ++i) {
    const Tddd dX = neighbor_positions[i] - origin;
    const double u = Dot(dX, u_axis);
    const double v = Dot(dX, v_axis);
    const double w = Dot(dX, n_axis);
    const std::array<double, 5> row = {u * u, u * v, v * v, u, v};
    for (int j = 0; j < 5; ++j) {
      ATb[j] += row[j] * w;
      for (int k = 0; k <= j; ++k)
        ATA[j][k] += row[j] * row[k];
    }
  }

  // Cholesky 分解 (5×5)
  std::array<std::array<double, 5>, 5> L = {};
  bool cholesky_ok = true;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = ATA[i][j];
      for (int k = 0; k < j; ++k)
        s -= L[i][k] * L[j][k];
      if (i == j) {
        if (s <= 1e-20) { cholesky_ok = false; break; }
        L[i][j] = std::sqrt(s);
      } else {
        L[i][j] = s / L[j][j];
      }
    }
    if (!cholesky_ok) break;
  }
  if (!cholesky_ok)
    return result;

  // 前進代入 L y = ATb
  std::array<double, 5> y = {};
  for (int i = 0; i < 5; ++i) {
    double s = ATb[i];
    for (int k = 0; k < i; ++k)
      s -= L[i][k] * y[k];
    y[i] = s / L[i][i];
  }
  // 後退代入 L^T x = y
  std::array<double, 5> coeffs = {};
  for (int i = 4; i >= 0; --i) {
    double s = y[i];
    for (int k = i + 1; k < 5; ++k)
      s -= L[k][i] * coeffs[k];
    coeffs[i] = s / L[i][i];
  }

  result.a = coeffs[0];
  result.b = coeffs[1];
  result.c = coeffs[2];
  result.d = coeffs[3];
  result.e = coeffs[4];
  result.valid = true;
  return result;
}

// ---------------------------------------------------------------------------
// fitQuadricLocal_unroll: 完全手動展開バージョン
// ---------------------------------------------------------------------------
//
// fitQuadricLocal と数学的に等価だが、5 行列の全ループを展開し、配列アクセスを
// 全て スカラ変数 (register) に置き換えた版。N についての外側ループだけが残る。
//
// 効果の源:
//   - A^T A 下三角の 15 要素、ATb の 5 要素、Cholesky L の 15 要素、forward / backward
//     の y, x をすべて独立スカラにしているので、コンパイラがレジスタに載せやすい
//   - fma を明示 (row 累積部) し、FPU スループット向上
//   - 配列索引の読み替え・境界判定が消えて命令数が減る
//   - 元コードと同じ "Cholesky 失敗時は invalid を返す" セマンティクスを維持
//
// 数値的に同一結果を期待するが、fma / 非 fma の違いによる丸め差が出うる。

inline QuadricFitResult fitQuadricLocal_unroll(const Tddd& origin,
                                                const Tddd& u_axis, const Tddd& v_axis, const Tddd& n_axis,
                                                const std::vector<Tddd>& neighbor_positions) {
  QuadricFitResult result;
  const int N = static_cast<int>(neighbor_positions.size());
  if (N < 5)
    return result;
  result.frame = {u_axis, v_axis, n_axis, true};

  // A^T A 下三角 (15 要素) + A^T b (5 要素) をスカラで累積
  double A00 = 0.0;
  double A10 = 0.0, A11 = 0.0;
  double A20 = 0.0, A21 = 0.0, A22 = 0.0;
  double A30 = 0.0, A31 = 0.0, A32 = 0.0, A33 = 0.0;
  double A40 = 0.0, A41 = 0.0, A42 = 0.0, A43 = 0.0, A44 = 0.0;
  double b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;

  for (int i = 0; i < N; ++i) {
    const Tddd dX = neighbor_positions[i] - origin;
    const double u = Dot(dX, u_axis);
    const double v = Dot(dX, v_axis);
    const double w = Dot(dX, n_axis);
    // row = [u², uv, v², u, v]
    const double r0 = u * u;
    const double r1 = u * v;
    const double r2 = v * v;
    const double r3 = u;
    const double r4 = v;

    // A^T A 下三角の 15 要素を累積 (すべて fma)
    A00 = std::fma(r0, r0, A00);
    A10 = std::fma(r1, r0, A10); A11 = std::fma(r1, r1, A11);
    A20 = std::fma(r2, r0, A20); A21 = std::fma(r2, r1, A21); A22 = std::fma(r2, r2, A22);
    A30 = std::fma(r3, r0, A30); A31 = std::fma(r3, r1, A31); A32 = std::fma(r3, r2, A32); A33 = std::fma(r3, r3, A33);
    A40 = std::fma(r4, r0, A40); A41 = std::fma(r4, r1, A41); A42 = std::fma(r4, r2, A42); A43 = std::fma(r4, r3, A43); A44 = std::fma(r4, r4, A44);

    // A^T b
    b0 = std::fma(r0, w, b0);
    b1 = std::fma(r1, w, b1);
    b2 = std::fma(r2, w, b2);
    b3 = std::fma(r3, w, b3);
    b4 = std::fma(r4, w, b4);
  }

  // Cholesky 分解 A = L L^T (下三角)
  //   L[i][i] = sqrt(A[i][i] - Σ_{k<i} L[i][k]²)
  //   L[i][j] = (A[i][j] - Σ_{k<j} L[i][k]·L[j][k]) / L[j][j]   (j<i)
  // 対角が非正 or 極小なら fit 失敗として早期 return (元コードと同一セマンティクス)
  constexpr double kMinPivot = 1e-20;
  if (A00 <= kMinPivot) return result;
  const double L00 = std::sqrt(A00);

  const double L10 = A10 / L00;
  const double d11 = A11 - L10 * L10;
  if (d11 <= kMinPivot) return result;
  const double L11 = std::sqrt(d11);

  const double L20 = A20 / L00;
  const double L21 = (A21 - L20 * L10) / L11;
  const double d22 = A22 - L20 * L20 - L21 * L21;
  if (d22 <= kMinPivot) return result;
  const double L22 = std::sqrt(d22);

  const double L30 = A30 / L00;
  const double L31 = (A31 - L30 * L10) / L11;
  const double L32 = (A32 - L30 * L20 - L31 * L21) / L22;
  const double d33 = A33 - L30 * L30 - L31 * L31 - L32 * L32;
  if (d33 <= kMinPivot) return result;
  const double L33 = std::sqrt(d33);

  const double L40 = A40 / L00;
  const double L41 = (A41 - L40 * L10) / L11;
  const double L42 = (A42 - L40 * L20 - L41 * L21) / L22;
  const double L43 = (A43 - L40 * L30 - L41 * L31 - L42 * L32) / L33;
  const double d44 = A44 - L40 * L40 - L41 * L41 - L42 * L42 - L43 * L43;
  if (d44 <= kMinPivot) return result;
  const double L44 = std::sqrt(d44);

  // 前進代入: L y = ATb
  const double y0 = b0 / L00;
  const double y1 = (b1 - L10 * y0) / L11;
  const double y2 = (b2 - L20 * y0 - L21 * y1) / L22;
  const double y3 = (b3 - L30 * y0 - L31 * y1 - L32 * y2) / L33;
  const double y4 = (b4 - L40 * y0 - L41 * y1 - L42 * y2 - L43 * y3) / L44;

  // 後退代入: L^T x = y
  const double x4 = y4 / L44;
  const double x3 = (y3 - L43 * x4) / L33;
  const double x2 = (y2 - L42 * x4 - L32 * x3) / L22;
  const double x1 = (y1 - L41 * x4 - L31 * x3 - L21 * x2) / L11;
  const double x0 = (y0 - L40 * x4 - L30 * x3 - L20 * x2 - L10 * x1) / L00;

  result.a = x0;
  result.b = x1;
  result.c = x2;
  result.d = x3;
  result.e = x4;
  result.valid = true;
  return result;
}


// ---------------------------------------------------------------------------
// 2次曲面フィット結果から主曲率・主方向を計算
// ---------------------------------------------------------------------------
//
// 単純に Hessian [[2a,b],[b,2c]] の固有値を使うのではなく，
//
//   S = I^{-1} II
//
// を使う。
//
// I  : 第一基本形式
// II : 第二基本形式
// S  : shape operator
//
// こうしておくと，Monge patch の傾き d,e が 0 でない場合でも，
// その傾きを含めて主曲率を補正できる。

inline PrincipalCurvatureResult principalCurvaturesFromQuadric(const QuadricFitResult& qf) {
  PrincipalCurvatureResult result;
  if (!qf.valid)
    return result;

  double a = qf.a, b = qf.b, c = qf.c, d = qf.d, e = qf.e;
  const Tddd& u_axis = qf.frame.u_axis;
  const Tddd& v_axis = qf.frame.v_axis;
  const Tddd& n_w = qf.frame.n_axis;

  // 第一基本形式 I = [[E,F],[F,G]]
  // w_u(0,0)=d, w_v(0,0)=e より
  double E = 1.0 + d * d;
  double F = d * e;
  double G = 1.0 + e * e;
  double det_I = E * G - F * F;
  if (!(det_I > 1e-20))
    return result;

  // 第二基本形式 II = [[L,M],[M,N]]
  // 単位法線の z 成分による補正が 1/sqrt(1+d^2+e^2) に相当する
  double denom = std::sqrt(1.0 + d * d + e * e);
  double L = 2.0 * a / denom;
  double M = b / denom;
  double N = 2.0 * c / denom;

  // shape operator S = I^{-1} II （一般に非対称）
  // 2x2 行列を明示的に展開している
  double s00 = (L * G - M * F) / det_I;
  double s01 = (G * M - F * N) / det_I;
  double s10 = (M * E - L * F) / det_I;
  double s11 = (N * E - M * F) / det_I;

  // shape operator の固有値 = 主曲率
  // det(S) = s00*s11 - s01*s10 = (LN - M^2) / (EG - F^2)
  double trace = s00 + s11;
  double disc = std::sqrt(std::max(trace * trace - 4.0 * (s00 * s11 - s01 * s10), 0.0));
  double lambda1 = 0.5 * (trace + disc);
  double lambda2 = 0.5 * (trace - disc);

  // 固有ベクトルを (u,v) 平面で作ってから 3D へ戻す
  // (S - λI)v = 0 の第2行: s10*v_u + (s11-λ)*v_v = 0 → v = (λ-s11, s10)
  Tddd PD1_3d, PD2_3d;
  {
    double ev_u = lambda1 - s11;
    double ev_v = s10;
    double ev_len = std::sqrt(ev_u * ev_u + ev_v * ev_v);
    if (ev_len > 1e-12) {
      PD1_3d = (ev_u / ev_len) * u_axis + (ev_v / ev_len) * v_axis;
    } else {
      PD1_3d = u_axis; // 等方的で方向が決めにくいので任意に u 軸を返す
    }
    PD2_3d = CrossDouble(n_w, PD1_3d);
  }

  // |k1| >= |k2| となるように並べ替える
  if (std::abs(lambda1) >= std::abs(lambda2)) {
    result.k1 = lambda1;
    result.k2 = lambda2;
    result.PD1 = PD1_3d;
    result.PD2 = PD2_3d;
  } else {
    result.k1 = lambda2;
    result.k2 = lambda1;
    result.PD1 = PD2_3d;
    result.PD2 = PD1_3d;
  }
  result.kmax = std::max(std::abs(result.k1), std::abs(result.k2));
  result.valid = true;
  return result;
}

// ---------------------------------------------------------------------------
// 点群から直接，主曲率・主方向を計算する簡易関数
// ---------------------------------------------------------------------------

inline PrincipalCurvatureResult computePrincipalCurvatures(const Tddd& origin,
                                                           const Tddd& u_axis, const Tddd& v_axis, const Tddd& n_axis,
                                                           const std::vector<Tddd>& neighbor_positions) {
  // 完全手動展開版 (15 スカラの Cholesky + fma) を使用。
  // 元の fitQuadricLocal と数学的に等価だが、A^T A / Cholesky / 前進後退代入の
  // 全ループを展開してレジスタ割当と fma 発行を最大化している。
  // remesh scenarios ループの内側で頂点ごとに呼ばれるため、この単位での高速化が効く。
  auto qf = fitQuadricLocal_unroll(origin, u_axis, v_axis, n_axis, neighbor_positions);
  return principalCurvaturesFromQuadric(qf);
}

// ---------------------------------------------------------------------------
// 辺の角度カバー量 θ（曲面忠実度の指標）
// ---------------------------------------------------------------------------
//
// 頂点 p における辺 e の角度カバー量 [rad]:
//   θ^(p)(e) = ||S(e_t)|| = sqrt(k1²(e·d1)² + k2²(e·d2)²)
//
// 形状作用素 S を辺ベクトルに作用させた結果のノルム。
// 辺の長さ × その方向の曲率 = その辺がカバーする角度（ラジアン）。

inline double edgeCurvatureAngle(const Tddd& edge,
                                  double k1, double k2,
                                  const Tddd& PD1, const Tddd& PD2) {
  double c1 = Dot(edge, PD1);
  double c2 = Dot(edge, PD2);
  return std::sqrt(k1 * k1 * c1 * c1 + k2 * k2 * c2 * c2);
}

// ---------------------------------------------------------------------------
// 辺の θ = 両端点の θ の max（安全側）
// ---------------------------------------------------------------------------
//
// 各頂点で独立に θ を計算し、大きい方を採用する。
// PD1 の符号合わせや平均は不要（2乗するので符号が消える）。

inline double edgeThetaMax(const Tddd& p0_pos, const Tddd& p1_pos,
                           const PrincipalCurvatureResult& p0_curv,
                           const PrincipalCurvatureResult& p1_curv) {
  Tddd edge = p1_pos - p0_pos;
  double theta0 = edgeCurvatureAngle(edge, p0_curv.k1, p0_curv.k2,
                                      p0_curv.PD1, p0_curv.PD2);
  double theta1 = edgeCurvatureAngle(edge, p1_curv.k1, p1_curv.k2,
                                      p1_curv.PD1, p1_curv.PD2);
  return std::max(theta0, theta1);
}

// ---------------------------------------------------------------------------
// cotan weight
// ---------------------------------------------------------------------------
//
// 辺 (a,b) を共有する 2 つの三角形
//
//   (a,b,p), (a,b,q)
//
// に対し，辺 (a,b) の向かい角を α, β とすると，
//
//   (cot α + cot β) / 2
//
// を返す。
//
// ここでは half-sum convention を採用しているので，
// Laplace-Beltrami を組む側で係数 1/2 を二重に掛けないこと。

inline double cotangentWeight(const Tddd& a, const Tddd& b, const Tddd& p, const Tddd& q) {
  Tddd pa = a - p, pb = b - p;
  Tddd qa = a - q, qb = b - q;
  double cross_p = Norm(Cross(pa, pb));
  double cross_q = Norm(Cross(qa, qb));
  double cot_alpha = (cross_p > 1e-20) ? Dot(pa, pb) / cross_p : 0.0;
  double cot_beta = (cross_q > 1e-20) ? Dot(qa, qb) / cross_q : 0.0;
  return 0.5 * (cot_alpha + cot_beta);
}

// ---------------------------------------------------------------------------
// 頂点の mixed area への 1 三角形からの寄与
// ---------------------------------------------------------------------------
//
// Meyer et al. 2003 の Voronoi-capped mixed area に沿った実装。
// 三角形が鈍角かどうかで寄与式を切り替える。

inline double mixedAreaContribution(const Tddd& p, const Tddd& a, const Tddd& b) {
  Tddd pa = a - p, pb = b - p, ab = b - a;
  double dot_p = Dot(pa, pb);
  double dot_a = Dot(-pa, ab);
  double dot_b = Dot(-pb, -ab);
  double area = 0.5 * Norm(Cross(pa, pb));

  if (dot_p < 0.0) {
    return area * 0.5; // 頂点 p が鈍角
  } else if (dot_a < 0.0 || dot_b < 0.0) {
    return area * 0.25; // a または b が鈍角
  } else {
    double cross_a = Norm(Cross(-pa, ab));
    double cross_b = Norm(Cross(-pb, -ab));
    double cot_a = (cross_a > 1e-20) ? dot_a / cross_a : 0.0;
    double cot_b = (cross_b > 1e-20) ? dot_b / cross_b : 0.0;
    return (cot_a * Dot(pb, pb) + cot_b * Dot(pa, pa)) / 8.0;
  }
}

// ---------------------------------------------------------------------------
// 角度欠損
// ---------------------------------------------------------------------------
//
// 内点:
//   2π - Σ angle
// 境界点:
//   π - Σ angle
//
// これを mixed area で割ると，離散ガウス曲率の近似になる。

inline double angularDeficit(const std::vector<double>& angles_at_vertex, bool is_boundary_vertex = false) {
  double sum = 0.0;
  for (double a : angles_at_vertex)
    sum += a;
  return (is_boundary_vertex ? M_PI : 2.0 * M_PI) - sum;
}

// ---------------------------------------------------------------------------
// edge collapse 後の法線変化を予測
// ---------------------------------------------------------------------------
//
// 1 つの頂点を targetX へ動かしたときの新旧法線の内積を返す。
// 値の意味:
//
//   1.0 に近い : ほとんど変化しない
//   0 付近      : 大きく回転する
//   負          : 面反転の可能性が高い
//
// 面積がほぼ 0 になる場合は -1.0 を返す。

inline double predictNormalChangeAfterCollapse(const Tddd& fp0, const Tddd& fp1, const Tddd& fp2,
                                               const Tddd& old_normal,
                                               const Tddd& targetX,
                                               int moving_vertex_index) {
  Tddd x0 = (moving_vertex_index == 0) ? targetX : fp0;
  Tddd x1 = (moving_vertex_index == 1) ? targetX : fp1;
  Tddd x2 = (moving_vertex_index == 2) ? targetX : fp2;
  Tddd new_normal = Cross(x1 - x0, x2 - x0);
  double new_area = Norm(new_normal);
  if (new_area < 1e-20)
    return -1.0;
  new_normal = new_normal / new_area;
  return Dot(new_normal, old_normal);
}

} // namespace surface_geometry

#endif
