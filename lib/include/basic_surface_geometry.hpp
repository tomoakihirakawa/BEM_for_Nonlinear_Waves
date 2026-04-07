#ifndef basic_surface_geometry_H
#define basic_surface_geometry_H

#include "basic_vectors.hpp"
#include "basic_linear_systems.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

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

  // 近傍点を局所座標系へ写し，最小二乗系 A x = b を作る
  std::vector<std::vector<double>> A(N, std::vector<double>(5, 0.0));
  std::vector<double> w_vals(N, 0.0);

  for (int i = 0; i < N; ++i) {
    Tddd dX = neighbor_positions[i] - origin;
    double u = Dot(dX, u_axis);
    double v = Dot(dX, v_axis);
    double w = Dot(dX, n_axis);
    A[i] = {u * u, u * v, v * v, u, v};
    w_vals[i] = w;
  }

  // SVD で解く
  std::vector<double> coeffs(5, 0.0);
  try {
    lapack_svd_solve(A, coeffs, w_vals);
  } catch (...) {
    return result;
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
    PD2_3d = Cross(n_w, PD1_3d);
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
  auto qf = fitQuadricLocal(origin, u_axis, v_axis, n_axis, neighbor_positions);
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
