#pragma once

#include "basic_linear_systems.hpp"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <filesystem>
#include <string>

extern std::string solver_type;
extern std::string coupling_type;
extern std::string preconditioner_type;
extern std::string ilu_neighborhood_type;
extern int ilu_kring_num;
extern double milu_omega;
extern double ilut_drop_tol;
extern int ilut_max_entries_per_row;
extern double ilut_pivot_min;
extern int schwarz_core_k;
extern int schwarz_overlap_k;
extern int schwarz_max_core_size;
extern int schwarz_max_block_size;
extern double schwarz_pivot_min;
extern double schwarz_diag_shift;
extern double solver_tol;
extern int solver_max_iter;
extern int solver_restart;
#if defined(USE_METAL_M2L)
extern bool g_metal_m2l_active; // true: use Metal GPU M2L (float), false: use CPU M2L (double)
// Metal M2L settings (read from settings.json)
extern bool use_metal_m2l;         // Enable Metal GPU acceleration
extern bool metal_m2l_threadgroup; // true: use threadgroup parallelization
extern bool metal_m2l_sort_terms;  // Sort terms for improved memory locality
#endif
extern std::vector<double> coupling_params;
extern bool use_pseudo_quadratic_element;
extern bool use_true_quadratic_element;
extern std::string nearfield_mode;  // "scalar" (default), "simd", "metal"
extern int g_p2m_quadrature_points; // P2M Dunavant quadrature points (1, 3, 6, 7)
extern double g_mac_theta;          // MAC criterion parameter (0.25 default)
#include "BEM_BoundaryValues.hpp"
#include "Network.hpp"
#include "dunavant_rules.hpp"
#include "lib_multipole_expansion.hpp"
#include "searcher.hpp"

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

## 境界値問題

### 基礎方程式

非粘性渦なし流れを仮定し，ラプラス方程式を満たす速度ポテンシャル$\phi(t,\bf{x})$によって流れ場$\bf{u}(t,\bf{x})=\nabla\phi(t,\bf{x})$を表す．水面，壁面，浮体表面における境界条件は，

$$
\begin{align}
\nabla\cdot\nabla \phi& = 0&&\text{in}&&{\bf x} \in \Omega(t),\\
\frac{\partial\phi}{\partial t} +\frac{1}{2}\nabla\phi\cdot\nabla\phi + g z &=0 &&\text{on}&&{\bf x} \in \Gamma^{\rm D}(t),\\
\phi_n + {{\bf u}_b}\cdot{{\bf n}_b} &=0&&\text{on}&&{\bf x}\in \Gamma^{\rm N}(t),
\end{align}
$$

ここで，
${\bf x} ={(x,y,z)}$は空間座標，${\bf u}_b$は物体の流速，
${\bf n}_b$は物体の外向き単位法線ベクトル，
$\nabla=(\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z})$
である．
また，$\phi_n$は境界面上での外向き法線方向の流速を表し，
境界面上の外向き単位法線ベクトル$\bf n$を使えば$\phi_n ={\nabla\phi}\cdot {\bf n}$で表される．

### 境界積分方程式（BIE）

**グリーンの定理**

任意の$\phi$，$G$に対して次が成り立つ（**グリーンの定理**）．

$$
\iiint_\Omega \left(G({\bf x},{\bf a})\nabla^2 \phi({\bf x}) - \phi({\bf x})\nabla^2 G({\bf x},{\bf a})\right)dV
= \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
$$


$\phi$がラプラス方程式$\nabla^2\phi=0$を満たし，$G=1/\|{\bf x}-{\bf a}\|$とすると，
グリーンの定理から$\phi$と$\phi_n$の関係式，BIEが得られる．

$$
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
$$

ここで，${\bf a}$は境界面上の位置ベクトルであり，この原点${\bf a}$を固定し${\bf x}$について面積分される．
$G$は任意のスカラー関数で$G=1/\|{\bf x}-{\bf a}\|$とすることで，グリーンの定理の体積積分が消え，BIEの左辺のように，
原点での立体角$\alpha\left( {\bf{a}} \right)$とポテンシャル$\phi( {\bf{a}})$の積だけが残る．

<img src="schematic_BIE.png" width="400px">

この式は，流体内部では，$\alpha ({\bf{a}})$は$1$とできる．
この式は，$\bf{a}$におけるポテンシャル$\phi ({\bf{a}})$が，右辺の１重層ポテンシャルと２重層ポテンシャルの和で表されることを示している．
$G=1/\|{\bf x}-{\bf a}\|$がラプラス方程式の基本解であり，$\phi$は境界におけるポテンシャルの分布である．

*/

struct calculateFluidInteraction {
  const Network* PasObj;
  std::unordered_set<networkFace*> actingFaces;
  Tddd simplified_drag, simplified_drag_torque;
  double area;
  calculateFluidInteraction(const auto& faces /*waterfaces*/, const Network* PasObjIN) : PasObj(PasObjIN) {
    // PasObjと接したfaceの頂点にpressureが設定されている前提
    int count = 0;
    for (const auto& f : faces)
      if (f->Neumann) {
        // A face is acting on the body if all its points are either in contact with the body or have penetrated it.
        if (std::ranges::all_of(f->getPoints(), [&](const auto& p) {
              auto effectiveFaces = getEffectiveContactFaces(p);
              return std::ranges::any_of(effectiveFaces, [&](const auto& F) { return F->getNetwork() == PasObj; });
            })) {
          auto [p0, p1, p2] = f->getPoints();
          auto result = this->actingFaces.emplace(f);
          if (result.second)
            count++;
        }
      }

    // calculate area
    std::array<double, 3> P012;
    std::array<std::array<double, 3>, 3> X012;
    area = 0.;
    for (const auto& f : this->actingFaces) {
      auto [p0, p1, p2] = f->getPoints();
      P012 = {p0->pressure, p1->pressure, p2->pressure};
      X012 = {p0->X, p1->X, p2->X};
      auto intpX = interpolationTriangleLinear0101(X012);
      for (const auto& [x0, x1, w0w1] : __GWGW10__Tuple)
        area += intpX.J(x0, x1) * w0w1;
    }
    std::cout << "接触している面の数:" << count << " 表面積:" << area << std::endl;
  };

  // \label{BEM:surfaceIntegralOfPressure}
  std::array<Tddd, 2> surfaceIntegralOfPressure() {
    Tddd force = {0., 0., 0.}, torque = {0., 0., 0.};
    std::array<double, 3> P012;
    std::array<std::array<double, 3>, 3> X012;
    for (const auto& f : this->actingFaces) {
      auto [p0, p1, p2] = f->getPoints();
      P012 = {p0->pressure, p1->pressure, p2->pressure};
      X012 = {p0->X, p1->X, p2->X}; // auto [pre0, pre1, pre2] = P012;
      // auto [X0, X1, X2] = X012;
      // this->force += (p0 + p1 + p2) / 3. * Cross(X1 - X0, X2 - X0) / 2.;
      auto intpP = interpolationTriangleLinear0101(P012);
      auto intpX = interpolationTriangleLinear0101(X012);
      auto n = TriangleNormal(X012);
      Tddd F_tmp = {0., 0., 0.}, T_tmp = {0., 0., 0.};
      for (const auto& [x0, x1, w0w1] : __GWGW10__Tuple) {
        auto f = intpP(x0, x1) * intpX.J(x0, x1) * w0w1 * n;
        F_tmp += f;
        T_tmp += Cross(intpX(x0, x1) - this->PasObj->COM, f);
      }
      force += F_tmp;
      torque += T_tmp;
    }
    return {force, torque};
  };

  std::array<Tddd, 2> surfaceIntegralOfVerySimplifiedDrag() {
    this->simplified_drag.fill(0.);
    this->simplified_drag_torque.fill(0.);
    std::array<std::array<double, 3>, 3> X012;
    for (const auto& f : this->actingFaces) {
      auto [p0, p1, p2] = f->getPoints();
      std::array<double, 3> P012 = {p0->pressure, p1->pressure, p2->pressure};
      auto X0 = p0->X;
      auto X1 = p1->X;
      auto X2 = p2->X;
      X012 = {X0, X1, X2};
      const Tddd relative_U0 = p0->u_potential_BEM - PasObj->velocityRigidBody(X0);
      const Tddd relative_U1 = p1->u_potential_BEM - PasObj->velocityRigidBody(X1);
      const Tddd relative_U2 = p2->u_potential_BEM - PasObj->velocityRigidBody(X2);
      // this->force += (p0 + p1 + p2) / 3. * Cross(X1 - X0, X2 - X0) / 2.;
      auto intpRelativeVelocity = interpolationTriangleLinear0101(T3Tddd{relative_U0, relative_U1, relative_U2});
      auto intpX = interpolationTriangleLinear0101(X012);
      const double nu = _WATER_NU_10deg_; // m2 /s
      Tddd drag_f;
      for (const auto& [x0, x1, w0w1] : __GWGW10__Tuple) {
        drag_f = nu * intpRelativeVelocity(x0, x1) * intpX.J(x0, x1) * w0w1;
        this->simplified_drag += drag_f;
        this->simplified_drag_torque += Cross(intpX(x0, x1) - this->PasObj->COM, drag_f);
      }
    }
    return {this->simplified_drag, this->simplified_drag_torque};
  };
};

// b@ -------------------------------------------------------------------------- */
// b@                                   BEM_BVP                                  */
// b@ -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

### BIEの離散化

$$
\alpha ({\bf a})\phi({\bf a})
= \iint_\Gamma {\left({
\frac{1}{\|{\bf x}-{\bf a}\|}
\nabla \phi ({\bf{x}}) + \phi ({\bf{x}})
\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}}
\right) \cdot {\bf{n}}({\bf{x}})dS}
$$

面は面上の節点を使って補間され，面積分はこの補間された面上に沿って行われる．
面の法線ベクトル${\bf n}=\frac{\frac{{\partial {\bf x}}}{{\partial \xi_0}}\times\frac{{\partial {\bf x}}}{{\partial \xi_1}}}{\left\|\frac{{\partial {\bf x}}}{{\partial \xi_0}}\times\frac{{\partial {\bf x}}}{{\partial \xi_1}}\right\|}$を代入し，BIEをGauss-Legendre積分で離散化すると，

$$
\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} {\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}_{k _\vartriangle}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}\left\|\frac{{\partial{{\bf x}_{k _\vartriangle}}}}{{\partial{\xi_0}}} \times \frac{{\partial{\bf{x}}_{k _\vartriangle}}}{{\partial{\xi_1}}}\right\|} \right)} }=
$$

$$
\alpha_{i_\circ}(\phi)_{i_\circ}-\sum\limits_{k_\vartriangle}\sum\limits_{{\xi_1},{w_1}} \sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{\bf x}_{k \vartriangle}}({\pmb{\xi}})-{{\bf x}_{i_\circ} }}{{{{\| {{{\bf x}_{k_\vartriangle}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}_{k_\vartriangle}}}{{\partial
{\xi_0}}}\times\frac{{\partial {\bf{x}}_{k_\vartriangle}}}{{\partial {\xi_1}}}\right)}\right)}
$$

離散化では，$\phi_{i_\circ}$と${\phi_n}_{i_\circ}$の係数を知りたいので，
$\phi_{k_\vartriangle}({\pmb{\xi}})$と${\phi_n}_{k_\vartriangle}({\pmb{\xi}})$と書くのではなく，
$\phi_{i_\circ}$と${\phi_n}_{i_\circ}$が見えるように$\phi_{k_\vartriangle}({\pmb{\xi}})$と${\phi_n}_{k_\vartriangle}({\pmb{\xi}})$の補間を書いている．

ここで，$\phi_{k_\vartriangle,j}$における$k_\vartriangle$は三角形要素の番号，$j$は三角形要素の頂点番号．
$N_j$は三角形要素の形状関数，${\pmb{\xi}}$は三角形要素の内部座標，$w_0,w_1$はGauss-Legendre積分の重み，$\alpha_{i_\circ}$は原点$i_\circ$における立体角，$\phi$はポテンシャル，$\phi_n$は法線方向のポテンシャル，$\bf{x}$は空間座標，${\bf x}_{i_\circ}$は原点の空間座標である．

* $\phi_{k_\vartriangle}$は補間で作った関数
* $\phi_{k_\vartriangle,j}$は補間を構成する節点$j$での値
* $\phi_{i_\circ}$はより直接的にある節点$i_\circ$での値

NOTE: この段階ではまだ，1.数値積分のパラメタと，2.形状関数のパラメタと元々の面都の対応関係は，指定していない．例えば，やり方によっては$\xi_1$のパラメタは，$\xi_0$に依存するかもしれない．

補間に使うパラメタを${\bf \xi}=(\xi_0, \xi_1)$として，よく使われる３節点を使う線形補間を使うことにする．
元の面に対応する，線形補間面は，パラメタ上では${\xi_0 + \xi_1 = 1}$を満たす範囲なので，
積分範囲は例えば$0\leq \xi_0 \leq 1, 0\leq \xi_1 \leq 1-\xi_0$となる．
しかし，数値積分につかう変数と重みの組み合わせは，コンパイルタイムに決めておき計算を効率化したいので，
この点で，変化する積分範囲は数値積分との相性が悪い．

### 線形三角要素

#### 線形三角要素

<img src="./img/schematic_linear_triangle_element.png" width="400px">

形状関数${\pmb N}_j({\pmb \xi}),{\pmb \xi}=(\xi_0,\xi_1)$は，$\xi_0,\xi_1$が$0$から$1$動くことで，範囲で三角要素全体を動くように定義している．

$$
{\pmb N}({\pmb \xi}) = (N_0({\pmb \xi}),N_1({\pmb \xi}),N_2({\pmb \xi})) = (\xi_0, - \xi_1 (\xi_0 - 1), (\xi_0-1)(\xi_1-1))
$$

####  線形三角要素のヤコビアン

線形三角要素のヤコビアンは，$\|\frac{\partial {\bf{x}}}{\partial {\xi_0}} \times \frac{\partial {\bf{x}}}{\partial {\xi_1}}\|$である．

```Mathematica
shape[t0_, t1_] := With[{t2 = 1 - t0 - t1, t0m1 = t0 - 1, t1m1 = t1 - 1}, {t0, -t1*t0m1, t0m1*t1m1}];
D0shape[t0_, t1_] = (D[shape[T0, t1], T0] /. T0 -> t0);
D1shape[t0_, t1_] = (D[shape[t0, T1], T1] /. T1 -> t1);
{a, b, c} = {{x0, y0, z0}, {x1, y1, z1}, {x2, y2, z2}}
FullSimplify[Cross[Dot[D[shape[T0, t1], T0], {a, b, c}], Dot[D[shape[t0, T1], T1], {a, b, c}]]]
FullSimplify[Cross[Dot[D[shape[T0, t1], T0], {a, b, c}], Dot[D[shape[t0, T1], T1], {a, b, c}]]/Cross[b - a, c - a]]
```

上の結果は，$1-\xi_0$となる．つまり，線形補間の場合，ヤコビアン内の外積は次のように，節点位置を使ってシンプルに計算できる．

$$
\frac{\partial {{\bf x}_{{k _\vartriangle}}}}{\partial {\xi_0}} \times \frac{\partial {{\bf x}_{{k _\vartriangle}}}}{\partial {\xi_1}} = (1-\xi_0) (({{\bf x} _{{k _\vartriangle}_1}}-{{\bf x} _{{k _\vartriangle}_0}})\times({{\bf x} _{{k _\vartriangle}_2}}-{{\bf x} _{{k _\vartriangle} _0}}))
= 2(1-\xi_0)A_{k_\vartriangle}{\bf n}_{k_\vartriangle}
$$

これを使えば，BIEは次のように簡単になる．

$$
\sum\limits_{k_\vartriangle}{2A_{k_\vartriangle}}
\sum\limits_{{\xi_1},{w_1}}
{\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left( {\sum\limits_{j=0}^2 {{{\left( {{\phi_n}} \right)}_{k_\vartriangle,j }}{N_{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}
(1-\xi_0)
} \right)} }=
$$

$$
\alpha_{i_\circ}(\phi)_{i_\circ}
-\sum\limits_{k_\vartriangle}{2A_{k_\vartriangle}{\bf n}_{k_\vartriangle}}\cdot
\sum\limits_{{\xi_1},{w_1}}
\sum\limits_{{\xi_0},{w_0}} {\left( {{w_0}{w_1}\left({\sum\limits_{j =0}^2{{{\left( \phi  \right)}_{k_\vartriangle,j }}{N_{j}}\left( \pmb{\xi } \right)} } \right)\frac{{{\bf x}_{k _\vartriangle}}({\pmb{\xi}})-{{\bf x}_{i_\circ} }}{{{{\| {{{\bf x}_{k_\vartriangle}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}}\|}^3}}}} (1-\xi_0)\right)}
$$

NOTE: ちなみに，$\frac{1-\xi_0}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x}_{i_\circ}}} \|}}$の分子に$1-\xi_0$があることで，
関数の特異的な変化を抑えることができる．プログラム上ではこの性質が利用できるように，この分数をまとめて計算している．

*/

std::array<double, 3> weight(double t0, double t1) {
  auto shape = [](double t0, double t1) -> std::array<double, 3> { return {t0, t1, 1 - t0 - t1}; };

  static const std::array<std::array<double, 2>, 3> vertex = {{{std::cos(M_PI / 2), std::sin(M_PI / 2)}, {std::cos(7 * M_PI / 6), std::sin(7 * M_PI / 6)}, {std::cos(11 * M_PI / 6), std::sin(11 * M_PI / 6)}}};
  auto s = shape(t0, t1);
  auto s_half_half = shape(2. / 3., 2. / 3.);
  auto s_minus_half_half = shape(-1. / 3., 2. / 3.);
  auto s_half_minus_half = shape(2. / 3., -1. / 3.);

  std::array<double, 4> distances;

  distances[0] = Norm(Dot(s, vertex));
  distances[1] = Norm(Dot(s, vertex) - Dot(s_half_half, vertex));
  distances[2] = Norm(Dot(s, vertex) - Dot(s_minus_half_half, vertex));
  distances[3] = Norm(Dot(s, vertex) - Dot(s_half_minus_half, vertex));

  auto weight_func = [&](const double x) {
    const double h = 1 + 1E-10;
    return std::pow(w_Linear(x, h), 1);
    // return std::pow(2. * x, -1.);
    // return std::pow(std::cosh(M_PI * x), -1.);
  };

  // double w0 = weight_func(distances[0]);
  double w1 = weight_func(distances[1]);
  double w2 = weight_func(distances[2]);
  double w3 = weight_func(distances[3]);
  double total = w1 + w2 + w3;
  // w0 /= total;
  w1 /= total;
  w2 /= total;
  w3 /= total;

  return {w1, w2, w3};
}

struct BEM_BVP {
  const bool Neumann = false;
  const bool Dirichlet = true;
  bool use_rigid_mode = true;
  std::vector<Network*> WATERS;
  lapack_lu* lu = nullptr;
  std::unique_ptr<ILU0_CRS> ilu_preconditioner = nullptr;
  std::unique_ptr<ILUT_CRS> ilut_preconditioner = nullptr;
  std::unique_ptr<SchwarzAdditive_CRS> schwarz_preconditioner = nullptr;
  V_d ilu_diag_inv_preconditioner;

  // Cached neighborhood for ILU sparsity (built once per BVP / topology state)
  std::vector<std::vector<networkFace*>> ilu_faces_to_integrate_by_point;
  std::vector<std::vector<networkFace*>> ilu_faces_to_integrate_by_midpoint; // true_quadratic midpoint k-ring
  int ilu_k_ring_num_cached = 0;
  std::size_t ilu_topology_hash_cached = 0;
  // Cache the built preconditioner across timesteps when geometry is unchanged.
  std::size_t ilu_geometry_hash_cached = 0;
  std::string preconditioner_type_built;
  std::string ilu_neighborhood_type_built;
  int ilu_k_ring_num_built = -1;
  int ilu_matrix_size_built = -1;
  double ilut_drop_tol_built = -1.0;
  int ilut_max_entries_per_row_built = -1;
  double ilut_pivot_min_built = -1.0;
  int schwarz_core_k_built = -1;
  int schwarz_overlap_k_built = -1;
  int schwarz_max_core_size_built = -1;
  int schwarz_max_block_size_built = -1;
  double schwarz_pivot_min_built = -1.0;
  double schwarz_diag_shift_built = -1.0;

  // Coordinate scaling for FMM numerical stability
  double coordinate_scale_factor_ = 1.0;
  bool use_coordinate_scaling_ = true;

  // Diagnostics for the most recent GMRES solve
  int last_gmres_total_iter = 0;
  double last_gmres_residual_norm = std::numeric_limits<double>::infinity();
  int last_converged_k = 10; // default value
  bool last_gmres_converged = false;
  struct GMRESResidualHotspot {
    int row = -1;
    int face_index = -1;
    int point_index = -1;
    int midpoint_index = -1;
    bool is_midpoint = false;
    bool dirichlet = false;
    bool neumann = false;
    bool corner = false;
    bool is_multiple_node = false;
    int contact_faces_count = 0;
    double residual = 0.0;
    double residual_abs = 0.0;
    double rhs = 0.0;
    double solution = 0.0;
    double diag_coeff = 0.0;
    Tddd position = {0., 0., 0.};
    Tddd endpoint0 = {0., 0., 0.};
    Tddd endpoint1 = {0., 0., 0.};
    std::string adjacent_face_flags;
    std::string dof_summary;
    std::string contact_summary;
  };
  std::vector<GMRESResidualHotspot> last_gmres_hotspots;
  // Diagnostics for preconditioner / sparsity (per most recent solve)
  double last_ilu_build_time = 0.0;
  double last_ilu_apply_time_sum = 0.0;
  double last_gmres_iter_time_sum = 0.0;
  double last_A_sparse_nnz = 0.0;
  double last_A_sparse_avg_nnz = 0.0;

  using T_PBF = std::tuple<netP*, netF*>;
  using mapTPBF_Tdd = std::map<T_PBF, Tdd>;
  using mapTPBF_mapTPBF_Tdd = std::map<T_PBF /*タプル*/, mapTPBF_Tdd>;
  using map_P_Vd = std::map<netP*, V_d>;

  //@ 各バケツでのモーメントを次数別に保存する．(ユニーク) p->{k,m,Yn,Y}ベクトル
  using uo_P_uoTiiTdd = std::unordered_map<networkPoint*, std::unordered_map<Tii /*k,m*/, Tdd /*YYn*/>>;
  using V_uo_P_uoTiiTdd = std::vector<uo_P_uoTiiTdd>;
  using VV_uo_P_uoTiiTdd = std::vector<V_uo_P_uoTiiTdd>;
  using VVV_uo_P_uoTiiTdd = std::vector<VV_uo_P_uoTiiTdd>;
  VV_d mat_ukn, mat_kn;
  V_d knowns, b_RHS;
  V_d diag_coeffs;
  std::vector<std::vector<std::array<double, 2>>> IGIGn;

  // GMRES / FMM related members
  std::vector<networkPoint*> points;
  // True quadratic midpoint lines (networkLine inherits target4FMM, used directly as FMM targets)
  std::vector<networkLine*> midpoint_lines;
  std::vector<std::array<double*, 2>> copy_map, copy_map_t;

  // Optional profiling hooks for FMM mat-vec breakdown (enabled by solveSystemGMRES).
  double* profile_fmm_update_time_sum = nullptr;
  double* profile_fmm_near_time_sum = nullptr;
  double* profile_fmm_far_time_sum = nullptr;
  std::size_t* profile_fmm_matvec_calls = nullptr;
  std::array<double, 6>* profile_fmm_update_step_time_sum = nullptr;

  std::filesystem::path output_directory; // settings.json の output_directory

  BEM_BVP(std::vector<Network*> WATERS) : WATERS(WATERS) {};
  ~BEM_BVP() {
    if (this->lu)
      delete this->lu;
  };

  void logLastGmresHotspots(std::ostream& os = std::cerr, int max_entries = -1) const {
    if (last_gmres_hotspots.empty())
      return;
    const int count = (max_entries < 0) ? static_cast<int>(last_gmres_hotspots.size()) : std::min(max_entries, static_cast<int>(last_gmres_hotspots.size()));
    os << Yellow << "[GMRES:debug] top residual rows=" << count << colorReset << std::endl;
    for (int i = 0; i < count; ++i) {
      const auto& h = last_gmres_hotspots[i];
      os << Yellow << "  #" << i + 1 << colorReset
         << " row=" << Green << h.row << colorReset
         << " |r|=" << Red << h.residual_abs << colorReset
         << " r=" << Red << h.residual << colorReset
         << " rhs=" << Cyan << h.rhs << colorReset
         << " x=" << Cyan << h.solution << colorReset
         << " diag=" << Cyan << h.diag_coeff << colorReset
         << " type=" << (h.is_midpoint ? "midpoint" : "vertex")
         << " face=" << h.face_index;
      if (h.is_midpoint)
        os << " midpoint_index=" << h.midpoint_index << " mid=" << h.position << " p0=" << h.endpoint0 << " p1=" << h.endpoint1;
      else
        os << " point=" << h.point_index << " pos=" << h.position;
      os << " flags={D:" << h.dirichlet
         << ",N:" << h.neumann
         << ",C:" << h.corner
         << ",M:" << h.is_multiple_node
         << "}"
         << " contact_faces=" << h.contact_faces_count;
      if (!h.adjacent_face_flags.empty())
        os << " adjacent_faces=" << h.adjacent_face_flags;
      os << std::endl;
      if (!h.dof_summary.empty())
        os << Yellow << "     dofs=" << colorReset << h.dof_summary << std::endl;
      if (!h.contact_summary.empty())
        os << Yellow << "     contact=" << colorReset << h.contact_summary << std::endl;
    }
  }

// Solver-specific implementations are split into dedicated headers.
#include "BEM_solveBVP_GMRES_FMM.hpp"
#include "BEM_solveBVP_LU.hpp"

  /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

  ### 高速多重極展開との関係

  GMRES法は，$A\cdot x$の計算を何度も行い，その線形和で解を近似するので，$A$をプログラム中で保持せずとも，$A\cdot x$を計算することができれば解を求めることができる．
  高速多重極展開は，この$A\cdot x$を高速に計算するための手法である．$A\cdot {\bf x}={\bf b}$のある行において，具体的な計算を考えてみる．

  $$
  \begin{align*}
  A\cdot x &= b \\
  \sum\limits_{j=0}^{N-1} A_{i_\circ,j}x_j &= b_{i_\circ} \\
  \end{align*}t
  $$
  　　
  $\sum\limits_{j=0}^{N-1} A_{{i_\circ},j}x_j = b_{i_\circ}$は，節点${i_\circ}$を原点節点としてBIEを離散化したものである．

  $A_{i,j}({\bf a}_i)$は，${\bf a}_i$に依存しており，${\bf a}_i$が変わると$A_{i,j}({\bf a}_i)$も変わる．
  しかし，これをソース点と観測点の関数の積と和の形に変形することできる．
  また，展開中心をソース点付近にとれば，ある変数が小さい場合限っては，その展開は早く収束する．
  ある変数とは具体的には，展開中心からソース点までの距離/展開中心から観測点までの距離である．

  */

  /* -------------------------------------------------------------------------- */

  void isSolutionFinite(const std::vector<double>& solution) const {
    for (const auto& val : solution)
      if (!isFinite(val)) {
        std::stringstream ss;
        ss << val;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
  };

  void isSolutionFinite(const auto& water) const {
    for (const auto& p : water.getBoundaryPoints()) {
      if (!isFinite(p->phiphin)) {
        std::stringstream ss;
        ss << p->phiphin;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
      for (const auto& [f, phi] : p->phiOnFace)
        if (!isFinite(phi)) {
          std::stringstream ss;
          for (const auto& [f, phi] : p->phiOnFace)
            ss << phi << ", ";
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
        }
      for (const auto& [f, phin] : p->phinOnFace)
        if (!isFinite(phin)) {
          std::stringstream ss;
          for (const auto& [f, phin] : p->phiOnFace)
            ss << phin << ", ";
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
        }
    }
  };

  // b! -------------------------------------------------------------------------- */
  // b!                                    solve        　　　　　　　　　　         */
  // b! -------------------------------------------------------------------------- */

  V_d ans;
  // V_d tmp_b_RHS;

  int matrix_size = 0;

  /*展開次数*/
  Buckets<std::shared_ptr<source4FMM<target4FMM>>, 10 /*展開項数*/> B_poles;

  std::array<double, 3> solve() {
    TimeWatch watch, watch_from_start;
    last_ilu_build_time = 0.0;
    last_ilu_apply_time_sum = 0.0;
    last_gmres_iter_time_sum = 0.0;
    last_A_sparse_nnz = 0.0;
    last_A_sparse_avg_nnz = 0.0;

    double time_setPhiPhinOnFace = watch()[0];

    this->matrix_size = setNodeFaceIndices(WATERS);
    // Re-populate midpoint per-face maps now that f2Index is (re-)populated.
    // On first call, setPhiPhinOnFace was called before setNodeFaceIndices, so f2Index was
    // empty and only bootstrap values ({nullptr → phi_mid/phin_mid}) were stored.
    // This call ensures per-face maps match the current f2Index (including CORNER 2-DOF keys).
    setPhiPhinOnFace(WATERS);
    std::cout << Red << "   unknown size : " << this->matrix_size << colorReset << std::endl;

    {
      int n_linear = 0, n_pseudo = 0, n_true = 0, n_total = 0;
      for (const auto& water : WATERS)
        for (const auto& f : water->getBoundaryFaces()) {
          if (f->isLinearElement)
            n_linear++;
          else if (f->isPseudoQuadraticElement)
            n_pseudo++;
          else if (f->isTrueQuadraticElement)
            n_true++;
          n_total++;
        }
      if (n_total > 0)
        std::cout << "   element types : linear=" << n_linear << " (" << 100. * n_linear / n_total << "%)"
                  << "  pseudo_quad=" << n_pseudo << " (" << 100. * n_pseudo / n_total << "%)"
                  << "  true_quad=" << n_true << " (" << 100. * n_true / n_total << "%)"
                  << "  total=" << n_total << std::endl;
    }

    double unknownsizse = this->matrix_size;

    // If mesh topology / indexing changes across timesteps, an old ILU is invalid.
    if (ilu_preconditioner && ilu_preconditioner->n != static_cast<int>(this->matrix_size)) {
      ilu_preconditioner.reset();
    }
    if (ilut_preconditioner && ilut_preconditioner->n != static_cast<int>(this->matrix_size)) {
      ilut_preconditioner.reset();
    }
    if (schwarz_preconditioner && schwarz_preconditioner->n != static_cast<int>(this->matrix_size)) {
      schwarz_preconditioner.reset();
    }
    if (!ilu_diag_inv_preconditioner.empty() && ilu_diag_inv_preconditioner.size() != this->matrix_size) {
      ilu_diag_inv_preconditioner.clear();
    }

    double time_setup = 0.;
    double time_solve = 0.0;

    if (solver_type == "LU") {
      solveLU(watch, time_setup, time_solve);
    } else if (solver_type == "GMRES") {
      solveGMRES(watch, time_setup, time_solve);
    }

    std::cout << Green << "Elapsed time for solving BIE: " << Red << watch_from_start() << colorReset << " s\n";

    return {time_setup, time_solve, unknownsizse};
  };

  // b! ------------------------------------------------------------------------------ */
  // b!                             solve phi_t and phi_n_t                            */
  // b! ------------------------------------------------------------------------------ */

  /*DOC_EXTRACT 0_4_0_1_FLOATING_BODY_SIMULATION

  ## 浮体動揺解析

  BEM-MELで浮体動揺解析ができるようにするのは簡単ではない．
  浮体に掛かる圧力の計算に必要な$\phi_t$が簡単には求まらないためである．
  これに関しては，\cite{Wu2003}や\cite{Ma2009a}が参考になる．

  ### 浮体の運動方程式

  <img src="schematic_float.png" width="400px" />

  浮体の重心の運動方程式：

  $$
  m \frac{d {\mathbfit U}_{\rm c}}{d t} = {\mathbfit F}_{\text {ext }}+{\mathbfit F}_{\text {hydro }}, \quad
  {\mathbfit I} \frac{d {\mathbfit \Omega}_{\rm c}}{d t} = {\mathbfit T}_{\text {ext }}+{\mathbfit T}_{\text {hydro }}
  $$

  NOTE:これらの変数は固定されたグローバル座標系上での変数である．大抵の場合$m$は座標系によらないが，${\mathbfit I}$は，浮体の基本姿勢において定義されたものであり，（浮体の座標系においては変化しないがグローバル座標系においては）浮体の姿勢によって変化する．

  ${\mathbfit U}_{\rm c}$は浮体の移動速度．
  ${\mathbfit F}_{\text {ext }}$は重力などの外力，${\mathbfit F}_{\text {hydro }}$は水の力，${\mathbfit T}_{\text {ext }}$は外力によるトルク，${\mathbfit T}_{\text {hydro }}$は水の力によるトルク．
  浮体が流体から受ける力${\mathbfit F}_{\text {hydro }}$は，浮体表面の圧力$p$を積分することで得られ，
  また圧力$p$は速度ポテンシャル$\phi$を用いて，以下のように書ける．

  \ref{BEM:surfaceIntegralOfPressure}{圧力積分}と
  \ref{BEM:surfaceIntegralOfTorque}{トルクの積分}：

  $$
  {\mathbfit F}_{\text {hydro }}=\iint_{\Gamma_{\rm float}} p{\mathbfit n}  d S, \quad
  {\mathbfit T}_{\text {hydro }}=\iint_{\Gamma_{\rm float}} ({\bf x}-{\bf x}_{\rm c})\times (p{\mathbfit n})  d S, \quad
  p= p({\bf x}) =-\rho\left(\frac{\partial \phi}{\partial t}+\frac{1}{2} \|\nabla \phi\|^{2}+g z\right)
  $$

  $\frac{\partial \phi}{\partial t}$を$\phi_t$と書くことにする．この$\phi_t$は陽には求められない．
  そこで，$\phi$と似た方法，BIEを使った方法で$\phi_t$を求める．$\phi$と$\phi_n$の間に成り立つ境界積分方程式と全く同じ式が，$\phi_t$と$\phi_{nt}$の間にも成り立つ：

  $$
  \alpha ({\bf{a}})\phi_t ({\bf{a}}) = \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi_t ({\bf{x}}) - \phi_t ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
  \quad\text{on}\quad{\bf x} \in \Gamma(t).
  $$

  */

  /*DOC_EXTRACT 0_4_1_FLOATING_BODY_SIMULATION

  ### 加速度の計算の難しさ

  これは，浮体表面の圧力の計算の困難，もっと言えば$\phi_t$の計算の困難に起因する． \cite{Ma2009}によると，$\phi_t$の計算方法として以下の４つの方法が提案されている．

  1. 間接的法 (indirect method) ：補助関数を使う方法
  2. モード分解法 (mode-decomposition method)
  3. Dalena and Tanizawa's method
  4. Caoの反復法 (iterative method) \cite{Cao1994}
  5. Maの方法 \cite{Ma2009}

  #### 間接的法，モード分解法

  間接的法，モード分解法は，$\phi$に関するBVPと似ているが異なる新たなBVPを解く必要がある．
  この新たなBIEの境界条件は違うが，係数行列は同じ（私は違うと思うのだが）らしい．
  ただ悩ましいので，LU分解のような直接法なら逆行列を保持するので余計な計算が発生しないが，直接法はそもそも遅い．
  そこで反復法を使いたいが，反復法は逆行列を保持しないので，毎回係数行列を計算する必要がある，というジレンマがある．

  #### Dalena and Tanizawa's method

  Dalena and Tanizawa's methodは，$\phi$に関するBVPと全く違うBVPを解く必要があり，係数行列も違うので，新たに行列を構成する必要がある．
  \cite{Feng2017}によると，この方法は境界面の局所的な曲率を用いる必要があるため，３次元解析に利用することがとても難しい．

  #### 反復法

  Caoの反復法は，新たなBVPを解く必要がないので，上の問題はないらしい\cite{Ma2009}．
  （これは間違いで，直接法を使った場合はそうだが，反復法（GMRESのような）を使うなら，このCaoの反復法の内部で反復法（GMRESなど）をする必要があり時間がかかり，
  初めの２つに優っているとは言えない．同じ程度の時間がかかる．）

  #### Maの反復法

  | Method |  |
  |:---:|:---:|
  | Indirect method | 浮体１つに対して６つ，新しいBIEを立てる．新たに解く必要があり遅い |
  | Mode-decomposition method | 浮体１つに対して7つ，新しいBIEを立てる．新たに解く必要があり遅い |
  | Dalena and Tanizawa's method | ? |
  | Cao's iterative method | 直接法で解くなら同じBIE係数行列を使えるので，速い．反復法なら，反復法の内部で反復法をするので遅い．|
  | Ma's iterative method | 直接法で解くなら同じBIE係数行列を使えるので，速い．反復法なら，反復法の内部で反復法をするので遅い．|

  ### $\phi_t$と$\phi_{nt}$に関するBIEの解き方（と$\phi_{nt}$の与え方）

  $\phi_t$と$\phi_{nt}$に関するBIEを解くためには，ディリクレ境界には$\phi_t$を，ノイマン境界には$\phi_{nt}$を与える．

  #### ディリクレ節点の$\phi_{nt}$の与え方(水面：圧力が既知，$\phi$が既知)

  このディリクレ境界では，圧力が与えられていないので，このBiEにおいては，ノイマン境界条件を与える．
  ただし，壁が完全に固定されている場合，$\phi_{nt}$は0とする．

  #### ディリクレ節点の$\phi_{t}$の与え方($\phi$を与える造波装置：圧力が未知，$\phi$が既知)

  ディリクレ境界では$\phi_t$は，圧力が大気圧と決まっているので，ベルヌーイの圧力方程式から$\phi_t$を求めることができる．

  #### ノイマン節点での$\phi_{nt}$の与え方

  境界面が静止しているかどうかに関わらず，流体と物体との境界では，境界法線方向速度が一致する．
  浮体重心${\bf x}_c$から境界面上の点$\bf x$までの位置ベクトルを$\mathbfit r = {\bf x} - {\bf x}_c$とする．
  表面上のある点の移動速度$\frac{d\mathbfit r}{dt}$と流体粒子の流速$\nabla \phi$の間には，次の境界条件が成り立つ．

  $$
  {\bf n}\cdot\frac{d\mathbfit r}{dt} =  {\bf n} \cdot \nabla \phi,\quad \frac{d\mathbfit r}{dt} = \mathbfit U_{\rm c} + {\mathbfit \Omega}_{\rm c} \times \mathbfit r
  $$

  物体上のある点ではこれが常に成り立つ．

  これを微分することで，$\phi_{nt}$を$\phi$と加速度$\frac{d{\mathbfit U}_{\rm c}}{dt}$と角加速度$\frac{d{\mathbfit \Omega}_{\rm c}}{dt}$を使って表すことができる．
  \cite{Wu1998}

  $$
  \begin{aligned}
  &\rightarrow& 0& =\frac{d}{dt}\left({\bf n}\cdot \left(\frac{d\mathbfit r}{dt}-\nabla \phi\right)\right) \\
  &\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\mathbfit r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \frac{d}{dt}\left(\frac{d\mathbfit r}{dt}-\nabla \phi\right)\\
  &\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\mathbfit r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\mathbfit r}{dt^2}-\left(\frac{\partial}{\partial t}+\frac{d{\mathbfit r}}{dt}\cdot\nabla\right)\nabla \phi\right)\\
  &\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\mathbfit r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\mathbfit r}{dt^2}- {\nabla \phi_t - \left(\frac{d\mathbfit r}{dt} \cdot \nabla\right)\nabla \phi}\right)\\
  &\rightarrow& \phi_{nt}& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\mathbfit r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\mathbfit r}{dt^2} - \frac{d\mathbfit r}{dt} \cdot (\nabla\otimes\nabla \phi) \right)
  \end{aligned}
  $$

  ここの$\frac{d{\bf n}}{dt}$と$\frac{d^2\mathbfit r}{dt^2}$は，${\mathbfit U}_{\rm c}$と$\mathbfit \Omega_{\rm c}$を用いて，

  $$
  \frac{d^2\mathbfit r}{dt^2}
  = \frac{d}{dt}\left({\mathbfit U}_{\rm c} + \mathbfit \Omega_{\rm c} \times \mathbfit r\right)
  = \frac{d{\mathbfit U}_{\rm c}}{dt} + \frac{d{\mathbfit \Omega_{\rm c}}}{dt} \times \mathbfit r + \mathbfit \Omega_{\rm c} \times \frac{d\mathbfit r}{dt}
  ,\quad \frac{d{\bf n}}{dt} = {\mathbfit \Omega}_{\rm c}\times{\bf n}
  $$

  $\frac{d \mathbfit r}{dt}$は\ref{velocityRigidBody}{`velocityRigidBody`}
  $\frac{d^2 \mathbfit r}{dt^2}$は\ref{accelRigidBody}{`accelRigidBody`}で計算する．

  \ref{BEM:phint_Neumann}{`phin_Neuamnn`}で$\phi_{nt}$を計算する．これは\ref{BEM:setPhiPhin_t}{`setPhiPhin_t`}で使っている．

  $\frac{d^2\mathbfit r}{dt^2}$を上の式に代入し，$\phi_{nt}$を求め，
  次にBIEから$\phi_t$を求め，次に圧力$p$を求める．
  そして，浮体の重さと慣性モーメントを考慮して圧力から求めた$\frac{d^2\mathbfit r}{dt^2}$は，
  入力した$\frac{d^2\mathbfit r}{dt^2}$と一致しなければならない．

  現状を整理すると，この浮体動揺解析において，知りたい未知変数は，浮体の加速度と角加速度だけ．
  しかし，浮体の没水面上にある節点での圧力$p$が得られないと，${\mathbfit F}_{\text {hydro }}$が得られず，運動方程式から浮体加速度が計算できない．
  圧力を計算するためには，$\phi_t$が必要で，$\phi_t$は簡単には得られない，という状況．

  物体の加速度は， 節点における$\{\phi_{nt0},\phi_{nt1},\phi_{nt2},..\} = \Phi_{nt}$が分かれば求まるが，
  逆に$\phi_{nt}$は$\frac{d\mathbfit U_{\rm c}}{dt}$と$\frac{d {\mathbfit \Omega} _{\rm c}}{d t}$が分かれば求まる．また，物体の角加速度に関しても同様である．

  $$
  m \frac{d\mathbfit U_{\rm c}}{dt} = {\mathbfit F} _{\text {ext }}+ F_{\text {hydro}}\left(\Phi_{nt}\left(\frac{d\mathbfit U_{\rm c}}{dt},\frac{d {\mathbfit \Omega} _{\rm c}}{d t}\right)\right),\quad
  {\mathbfit I} \frac{d {\mathbfit \Omega} _{\rm c}}{d t} = {\mathbfit T} _{\text {ext }}+{\mathbfit T} _{\text {hydro }}\left(\Phi_{nt}\left(\frac{d\mathbfit U_{\rm c}}{dt},\frac{d {\mathbfit \Omega} _{\rm c}}{d t}\right)\right)
  $$

  これを満たすように，$\Phi_{nt}$を求める．これは次のように書き換えて，根探し問題として解く．
  このプログラムでは，\ref{quasi_newton:broyden}{Broyden法}を使って，根探している．

  $$
  {\mathbfit 0} = m \frac{d\mathbfit U_{\rm c}}{dt} - {\mathbfit F} _{\text {ext }} - F_{\text {hydro}}\left(\Phi_{nt}\left(\frac{d\mathbfit U_{\rm c}}{dt},\frac{d {\mathbfit \Omega} _{\rm c}}{d t}\right)\right),\quad
  {\mathbfit 0} = {\mathbfit I} \frac{d {\mathbfit \Omega} _{\rm c}}{d t} - {\mathbfit T} _{\text {ext }} - {\mathbfit T} _{\text {hydro }}\left(\Phi_{nt}\left(\frac{d\mathbfit U_{\rm c}}{dt},\frac{d {\mathbfit \Omega} _{\rm c}}{d t} \right)\right)
  $$

  この式を，${\mathbfit Q}\left(\dfrac{d {\mathbfit U} _{\rm c}}{d t}, \dfrac{d {\mathbfit \Omega} _{\rm c}}{d t}\right)=(0,0,0,0,0,0)$
  として，これを満たすような$\dfrac{d {\mathbfit U} _{\rm c}}{d t}$と$\dfrac{d {\mathbfit \Omega} _{\rm c}}{d t}$を求める．
  $\phi_{nt}$はこれを満たした$\dfrac{d {\mathbfit U} _{\rm c}}{d t}$と$\dfrac{d {\mathbfit \Omega} _{\rm c}}{d t}$を用いて求める．

  $\phi_{nt}$は，\ref{BEM:setphint}{ここ}で与えている．

  この方法は，基本的には\cite{Cao1994}と同じ方法である．

  */

  /*DOC_EXTRACT 0_4_2_FLOATING_BODY_SIMULATION

  ### 補助関数を使った方法

  浮体動揺解析で問題となったのは，圧力の計算に使う$\phi_t\,{\rm on}\,🚢$が簡単には求まらないことであったが，
  $\iint_{\Gamma_{🚢}} \phi_t{\bf n}dS$と$\iint_{\Gamma_{🚢}}\phi_{t}({\bf x}-{\bf x}_c)\times{\bf n}dS$がわかればある場所の圧力はわからないが，
  🚢にかかる力は計算できるのでそれでも問題ない．

  体積積分がゼロとなるように，領域内でラプラス方程式を満たすような$\varphi$，
  そして$\Gamma _{🚢}$上ではこちらが望む$\varphi_n$となり，また$\Gamma \rm other$上では$\varphi=0$となる
  そんな$\varphi$をBIEを使って計算する．この$\varphi$を使うと次の式が成り立つ．
  （NOTE：境界上の全ての節点上で$\varphi$と$\varphi_n$が求まったとする）

  $$
  \begin{align*}
  0 &= \iint _\Gamma {\left( {\varphi\nabla {\phi_t} ({\bf{x}}) - {\phi_t} ({\bf{x}})\nabla \varphi} \right) \cdot {\bf{n}}({\bf{x}})dS}\\
  \rightarrow 0 &= \iint _{\Gamma _{🚢}+\Gamma _{🌊}+\Gamma _{\rm wall}} \varphi {\phi_{nt}} dS - \iint _{\Gamma _{🚢}+\Gamma _{🌊}+\Gamma _{\rm wall}} {\phi_t} \varphi_n dS\\
  \rightarrow 0 &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} \varphi {\phi_{nt}} dS - \iint _{\Gamma _{🚢}+\Gamma _{🌊}} {\phi_t} \varphi_n dS\\
  \rightarrow \iint _{\Gamma _{🚢}} {\phi_t} \varphi_n dS &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} \varphi {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} \varphi_n dS\\
  \rightarrow \iint_{\Gamma_{🚢}} \phi_t
  \begin{bmatrix}
  {\mathbfit n} \\
  ({\mathbfit x} - {\mathbfit x}_c) \times {\mathbfit n}
  \end{bmatrix} dS
  &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} {{\mathbfit\varphi}_{1-6}} {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} {{\mathbfit\varphi}_n}_{1-6} dS\\
  \end{align*}
  $$

  つまり，$\varphi_n$を適当に選べば，左辺は知りたかった積分となり，右辺の積分で計算できることになる．

  もし浮体がもう一つあると

  $$
  \begin{align*}
  \iint_{\Gamma_{🚢}} \phi_t
  \begin{bmatrix}
  {\mathbfit n} \\
  ({\mathbfit x} - {\mathbfit x}_c) \times {\mathbfit n}
  \end{bmatrix} dS
  & = \iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {{\mathbfit\varphi}_{1-6}} {\phi_{nt}} dS - \iint _{\Gamma_{🚤}+\Gamma _{🌊}} {\phi_t} {{\mathbfit\varphi}_n}_{1-6} dS\\
  \rightarrow \iint_{\Gamma_{🚢}} \phi_t
  \begin{bmatrix}
  {\mathbfit n} \\
  ({\mathbfit x} - {\mathbfit x}_c) \times {\mathbfit n}
  \end{bmatrix} dS
  & = \iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {{\mathbfit\varphi}_{1-6}} {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} {{\mathbfit\varphi}_n}_{1-6} dS
  \end{align*}
  $$

  同じように

  $$
  \begin{align*}
  \iint_{\Gamma_{🚤}} \phi_t
  \begin{bmatrix}
  {\mathbfit n} \\
  {(\mathbfit x} - {\mathbfit x}_c) \times {\mathbfit n}
  \end{bmatrix} dS
  & = \iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {{\mathbfit\varphi}_{7-12}} {\phi_{nt}} dS - \iint _{\Gamma _{🌊}} {\phi_t} {{\mathbfit\varphi}_n}_{7-12} dS
  \end{align*}
  $$

  $\iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {{\mathbfit\varphi}_{1-6}} {\phi_{nt}} dS$や
  $\iint _{\Gamma _{🚢}+\Gamma_{🚤}+\Gamma _{\rm wall}} {{\mathbfit\varphi}_{7-12}} {\phi_{nt}} dS$
  は加速度行列とある既知変数から成る行列の積で表される．こうして，運動方程式の${\mathbfit F}_{\text {hydro }}$と${\mathbfit T}_{\text {hydro }}$を加速度によって表すことができ，
  運動方程式は加速度だけに関する連立方程式となる．

  この方法は，\cite{Wu1996}，\cite{Kashiwagi2000}，\cite{Wu2003}で使用されている．
  この方法は，複数の浮体を考えていないが，\cite{Feng2017}はこれを基にして２浮体の場合でも動揺解析を行っている．

  */

  /* ------------------------------------------------------ */

  V_d initializeAcceleration(const std::vector<Network*>& rigidbodies) {
    V_d ACCELS_init;
    for (const auto& net : rigidbodies) {
      // if (net->interp_accel.size() > 3) {
      //    std::cout << Red << "interp_accel" << colorReset << std::endl;
      //    std::ranges::for_each(net->interp_accel(net->RK_Q.get_t()), [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
      // } else
      //   std::ranges::for_each(net->acceleration, [&](const auto &a_w) { ACCELS_init.emplace_back(a_w); });
      for (const auto& a_w : net->acceleration)
        ACCELS_init.emplace_back(a_w);
    }
    return ACCELS_init;
  }

  void insertAcceleration(const std::vector<Network*>& rigidbodies, const V_d& BM_X) {
    int i = 0;
    for (const auto& net : rigidbodies) {
      if (net->isFloatingBody) {
        double start_time = 0;
        if (net->inputJSON.at("velocity").size() > 1)
          start_time = std::stod(net->inputJSON.at("velocity")[1]);
        if (simulation_time < start_time) {
          i += net->acceleration.size();
        } else {
          net->acceleration[0] = BM_X[i++];
          net->acceleration[1] = BM_X[i++];
          net->acceleration[2] = BM_X[i++];
          net->acceleration[3] = BM_X[i++];
          net->acceleration[4] = BM_X[i++];
          net->acceleration[5] = BM_X[i++];
        }
      } else {
        i += net->acceleration.size();
      }
    }
  }

  /* -------------------------------------------------------------------------- */
  /* -------------------------------------------------------------------------- */
  /* -------------------------------------------------------------------------- */
  V_d Func(V_d ACCELS_IN, const std::vector<Network*> WATERS, const std::vector<Network*>& rigidbodies, V_d& ACCELS_OUT, const char* tag = nullptr) {
    TimeWatch watch;
    auto ACCELS = ACCELS_IN;

    //@ -------------------------------------------------------------------------- */
    //@   Coordinate Scaling for phi_t (same as phi solve)                          */
    //@ -------------------------------------------------------------------------- */
    use_coordinate_scaling_ = true;
    auto obj = WATERS[0];
    obj->setGeometricPropertiesForce();
    coordinate_scale_factor_ = obj->computeCharacteristicLength();
    if (coordinate_scale_factor_ > 1e-10) {
      std::cout << Magenta << "[BEM:phi_t] " << Cyan << "Coordinate scaling enabled" << Green << " L=" << coordinate_scale_factor_ << colorReset << std::endl;
      for (auto& net : WATERS)
        net->applyScaling(coordinate_scale_factor_);
    } else {
      use_coordinate_scaling_ = false;
      coordinate_scale_factor_ = 1.0;
    }

    //
    //* --------------------------------------------------- */
    //*                  加速度 --> phiphin_t                */
    //* --------------------------------------------------- */
    setPhiPhin_t(WATERS);

    // True quadratic: swap phi_mid/phin_mid with phi_t_mid/phin_t_mid
    // so that the solver transparently reads/writes phi_t values via phi_mid/phin_mid.
    if (use_true_quadratic_element) {
      for (auto* water : WATERS) {
        for (auto* l : water->getBoundaryLines()) {
          std::swap(l->phiphin[0], l->phiphin_t[0]);
          std::swap(l->phiphin[1], l->phiphin_t[1]);
          std::swap(l->phiOnFace, l->phitOnFace);
          std::swap(l->phinOnFace, l->phintOnFace);
        }
      }
    }

    knowns = getVectorFromBoundary(WATERS, this->matrix_size, [](const networkPoint* p, networkFace* f) { return p->phitOnFace.at(f); }, [](const networkPoint* p, networkFace* f) { return p->phintOnFace.at(f); });
    ans.resize(knowns.size(), 0);

    // Scale Neumann BC (phint) by L for coordinate scaling: ∂φ_t/∂n' = L × ∂φ_t/∂n
    if (use_coordinate_scaling_ && coordinate_scale_factor_ > 1e-10) {
      for (auto& p : this->points) {
        for (auto& [f, i] : p->f2Index) {
          if (isNeumannID_BEM(p, f)) {
            knowns[i] *= coordinate_scale_factor_;
          }
        }
      }
      // Midpoint Neumann BC scaling is handled by copyToFMM() inside solveSystemGMRES.
      // No pre-scaling of phinOnFace needed here.
    }

    if (solver_type == "LU") {
      static std::size_t phi_t_call_id = 0;
      const std::size_t call_id = ++phi_t_call_id;
      std::cout << Magenta << "  [BVP:phi_t] " << Cyan << "call=" << Yellow << call_id << Cyan;
      if (tag && tag[0] != '\0')
        std::cout << "  tag=" << Green << tag << Cyan;
      std::cout << "  solver=" << Yellow << "LU" << colorReset << std::endl;
      std::cout << "Solving for phi_t using LU..." << std::flush;
      _Pragma("omp parallel for") for (size_t i = 0; i < mat_kn.size(); ++i) ans[i] = Dot(mat_kn[i], knowns);
      if (this->lu != nullptr)
        this->lu->solve(ans /*解*/);

    } else if (solver_type == "GMRES") {
      static std::size_t phi_t_call_id = 0;
      const std::size_t call_id = ++phi_t_call_id;
      std::cout << Magenta << "  [BVP:phi_t] " << Cyan << "call=" << Yellow << call_id << Cyan;
      if (tag && tag[0] != '\0')
        std::cout << "  tag=" << Green << tag << Cyan;
      std::cout << "  solver=" << Yellow << "GMRES" << colorReset << std::endl;
      std::cout << "Solving for phi_t using GMRES..." << std::flush;
      //   前回の答えを初期値に使う
      std::vector<double> x0 = getVectorFromBoundary(WATERS, knowns.size(), [](const networkPoint* p, networkFace* f) { return p->phintOnFace.at(f); }, [](const networkPoint* p, networkFace* f) { return p->phitOnFace.at(f); });
      ans = solveSystemGMRES(true, x0);
    }

    // std::cout << "solved" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

    //@ -------------------------------------------------------------------------- */
    //@   Unscale computed ∂φ_t/∂n for coordinate scaling                           */
    //@ -------------------------------------------------------------------------- */
    // The computed ∂φ_t/∂n' (in scaled coords) needs to be DIVIDED by L
    if (use_coordinate_scaling_ && coordinate_scale_factor_ > 1e-10) {
      std::cout << Magenta << "[BEM:phi_t] " << Cyan << "Unscaling computed phint by L=" << coordinate_scale_factor_ << colorReset << std::endl;
      for (auto& p : this->points) {
        for (auto& [f, i] : p->f2Index) {
          if (isDirichletID_BEM(p, f)) {
            ans[i] /= coordinate_scale_factor_;
          }
        }
      }
      // Also unscale midpoint Dirichlet unknowns (phin_t is the unknown for Dirichlet DOFs)
      if (use_true_quadratic_element) {
        for (auto* water : WATERS) {
          for (auto* l : water->getBoundaryLines()) {
            for (const auto& [f, idx] : l->f2Index) {
              if (idx >= 0 && idx < static_cast<int>(ans.size()) && isDirichletID_BEM(l, f))
                ans[idx] /= coordinate_scale_factor_;
            }
          }
        }
      }
      // Midpoint Neumann phin_t unscaling not needed — phinOnFace was never pre-scaled
      // (scaling is handled only in copyToFMM on the FMM copy).
    }

    storePhiPhin_t(WATERS, ans);

    // True quadratic: swap back phi_mid/phin_mid to restore phi values.
    // After storePhiPhin_t, the per-face maps contain solved phi_t values.
    // After swap-back, phitOnFace/phintOnFace will have the solved phi_t values,
    // and phiOnFace/phinOnFace will be restored to original phi values.
    if (use_true_quadratic_element) {
      for (auto* water : WATERS) {
        for (auto* l : water->getBoundaryLines()) {
          std::swap(l->phiphin[0], l->phiphin_t[0]);
          std::swap(l->phiphin[1], l->phiphin_t[1]);
          std::swap(l->phiOnFace, l->phitOnFace);
          std::swap(l->phinOnFace, l->phintOnFace);
        }
      }
    }

    copyToFMM(true);

    // std::cout << Green << "storePhiPhin_t" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

    //@ -------------------------------------------------------------------------- */
    //@   Remove coordinate scaling before pressure calculation                     */
    //@ -------------------------------------------------------------------------- */
    if (use_coordinate_scaling_) {
      std::cout << Magenta << "[BEM:phi_t] " << Cyan << "Removing coordinate scaling" << colorReset << std::endl;
      for (auto& net : WATERS)
        net->removeScaling();
      use_coordinate_scaling_ = false;
      coordinate_scale_factor_ = 1.0;
    }

    //* --------------------------------------------------- */
    //*                 phiphin_t --> 圧力                   */
    //* --------------------------------------------------- */

    for (const auto water : WATERS)
      for (const auto& p : water->getBoundaryPoints())
        for (const auto& [f, i] : p->f2Index) {
          if (isDirichletID_BEM(p, f))
            p->pressure = p->pressure_BEM = 0;
          else
            p->pressure = p->pressure_BEM = -_WATER_DENSITY_ * (std::get<0>(p->phiphin_t) + 0.5 * Dot(p->u_total, p->u_total) + _GRAVITY_ * p->height());
        }

    // std::cout << Green << "pressure" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

    //* --------------------------------------------------- */
    //*              圧力 ---> 力 --> 加速度                  */
    //* --------------------------------------------------- */
    /*DOC_EXTRACT 0_4_0_2_FLOATING_BODY_SIMULATION

    実際の実験では，浮体のある基本的な姿勢における主慣性モーメントが与えられる．${\mathbfit I}$を主慣性モーメントテンソルとする．

    $$
    {\mathbfit I} = \begin{pmatrix}
    I_x & 0 & 0 \\
    0 & I_y & 0 \\
    0 & 0 & I_z
    \end{pmatrix}
    $$

    global座標における浮体の慣性モーメントテンソルを求めるには，次のように考えればいい．

    $$
    \begin{aligned}
    {\mathbfit I}\frac{d{\bf \Omega}_{\rm L}}{dt} &= {\bf T}_{\rm L}\\
    {\mathbfit I}{\rm R}_{g2l} \frac{d{\bf \Omega}_{\rm G}}{dt}& = {\rm R}_{g2l}{\bf T}_{\rm G}\\
    {\rm R}_{g2l}^{-1}{\mathbfit I}{\rm R}_{g2l} \frac{d{\bf \Omega}_{\rm G}}{dt}& = {\bf T}_{\rm G}\\
    \end{aligned}
    $$

    このことから，global座標における慣性モーメントテンソルは，次のようになる．

    $$
    {\mathbfit I}_{\rm G} = {\rm R}_{g2l}^{-1}{\mathbfit I}{\rm R}_{g2l}
    $$

    $R^{-1}$が$R^{\top}$であることと，$I^{-1}$が対角成分のみの行列であることを利用すれば，次のように書ける．

    $$
    \begin{align*}
    {\left({{I}_{G}}\right)}_{il}&={\left({{R}^{-{1}}}\right)}_{ij}{I}_{jk}{R}_{kl}={R}_{ji}{I}_{jj}{R}_{jl}={R}_{1i}{R}_{1l}{I}_{x}+{R}_{2i}{R}_{2l}{I}_{y}+{R}_{3i}{R}_{3l}{I}_{z}\\
    {\left({{I}_{G}^{-{1}}}\right)}_{il}&={\left({{R}^{-{1}}}\right)}_{ij}{\left({{I}^{-{1}}}\right)}_{jk}{R}_{kl}={R}_{ji}{\left({{I}^{-{1}}}\right)}_{jj}{R}_{jl}=\frac{{R}_{1i}{R}_{1l}}{{I}_{x}}+\frac{{R}_{2i}{R}_{2l}}{{I}_{y}}+\frac{{R}_{3i}{R}_{3l}}{{I}_{z}}
    \end{align*}
    $$

    この運動方程式から，求めたいのは$\frac{d{\bf \Omega}_{\rm G}}{dt}$である．これはとても簡単で，次のように求めることができる．

    $$
    \frac{d{\bf \Omega}_{\rm G}}{dt} = {\rm R}_{g2l}^{-1}{\mathbfit I}^{-1}{\rm R}_{g2l} {\bf T}_{\rm G}
    $$

    */
    int i = 0;
    std::vector<networkFace*> all_faces;
    for (auto& water : WATERS) {
      auto surfaces = water->getBoundaryFaces();
      all_faces.insert(all_faces.end(), surfaces.begin(), surfaces.end());
    }
    for (const auto& body : rigidbodies)
      if (body->isFloatingBody) {
        //$ ------------------------------ 係留索から受ける力とトルク ----------------------------- */
        //$ フェアリードの節点が隣の線要素から受けている張力ベクトル --> 浮体が受ける力とトルク
        std::array<double, 3> F_mooring = {0., 0., 0.}, T_mooring = {0., 0., 0.};
        //! simulateはアップデートの際に行なっておく．
        for (auto& mooring_line : body->mooringLines) {
          F_mooring += mooring_line->lastPoint->getForce();
          T_mooring += Cross(mooring_line->lastPoint->X - body->COM, mooring_line->lastPoint->getForce());
        }

        //@ ------------------------------ 浮体が流体力とトルク ----------------------------- */

        auto F_ext = _GRAVITY3_ * body->mass; // 外力（重力）
        auto tmp = calculateFluidInteraction(all_faces, body);
        auto [F_hydro, T_hydro] = tmp.surfaceIntegralOfPressure();
        std::cout << Red << "F_hydro = " << F_hydro << ", T_hydro = " << T_hydro << colorReset << std::endl;
        auto F = F_ext + F_hydro;
        auto T_GLOBAL = T_hydro;
        F += F_mooring;
        T_GLOBAL += T_mooring;

        for (const auto& [key, vec_string] : body->inputJSON.map_S_S) {
          if (key.contains("spring")) {
            // std::cout << "Applying Spring Force defined by [" << key << "] with params: " << vec_string << std::endl;
            //% ---------------------- 重心の並進移動によって伸びる線形バネによる係留 --------------- */
            // std::cout << "simple spring mooring" << std::endl;
            auto X_k = stod(vec_string);
            std::array<double, 3> init_fairleader_position = {X_k[0], X_k[1], X_k[2]};
            std::array<double, 3> current_fairleader_position = body->rigidTransformation(init_fairleader_position);
            std::array<double, 3> anchor = {X_k[3], X_k[4], X_k[5]};
            std::array<double, 3> dir = Normalize(anchor - current_fairleader_position);
            double L0 = Norm(init_fairleader_position - anchor);   // デフォルトの自然長は，初期の長さとする
            double L = Norm(current_fairleader_position - anchor); // 現在の長さ
            double kx, ky, kz;
            if (X_k.size() == 10) {
              kx = X_k[6];
              ky = X_k[7];
              kz = X_k[8];
              L0 = X_k[9]; // 自然長が与えられている場合は，L0を上書きする．
            } else if (X_k.size() == 9) {
              kx = X_k[6];
              ky = X_k[7];
              kz = X_k[8];
            } else if (X_k.size() == 8) {
              kz = ky = kx = X_k[6];
              L0 = X_k[7];
            } else if (X_k.size() == 7)
              kz = ky = kx = X_k[6];
            std::array<double, 3> f = Tddd{kx, ky, kz} * (L - L0) * dir;
            F += f;
            // std::cout << "Spring tension: " << f << ", Norm: " << Norm(f) << ", L: " << L << ", L0: " << L0 << ", diff: " << (L - L0) << std::endl;
          } else if (key.contains("linear_cable")) {
            // std::cout << "Applying Linear Cable Force defined by [" << key << "] with params: " << vec_string << std::endl;
            auto X_k = stod(vec_string);
            std::array<double, 3> init_fairleader_position = {X_k[0], X_k[1], X_k[2]};
            std::array<double, 3> current_fairleader_position = body->rigidTransformation(init_fairleader_position);
            std::array<double, 3> anchor = {X_k[3], X_k[4], X_k[5]};
            double kx, ky, kz;
            if (X_k.size() >= 9) {
              kx = X_k[6];
              ky = X_k[7];
              kz = X_k[8];
            } else if (X_k.size() == 7)
              kz = ky = kx = X_k[6];
            double diff = Norm(current_fairleader_position - anchor) - Norm(init_fairleader_position - anchor);
            Tddd direction = Normalize(anchor - current_fairleader_position);
            std::array<double, 3> f = Tddd{kx, ky, kz} * (anchor - current_fairleader_position);
            T_GLOBAL += Cross(current_fairleader_position - body->COM, f);
            F += f;
          } else if (key == "damping") {
            const auto c = stod(vec_string);
            std::array<double, 3> c_xyz = {c[0], c[1], c[2]};
            std::array<double, 3> c_abc = {c[3], c[4], c[5]};
            double start_t = 0, end_t = 1E10;
            if (c.size() == 7)
              start_t = c[6];
            else if (c.size() == 8) {
              start_t = c[6];
              end_t = c[7];
            }
            if (start_t <= simulation_time && simulation_time <= end_t) {
              F -= c_xyz * body->velocityTranslational();
              T_GLOBAL -= c_abc * body->velocityRotational();
            }
          }
        }

        //^ ------------------------------------------------------ */

        double gamma = 0.5;
        {
          const auto [mx, my, mz, Ix, Iy, Iz, IG, inv_IG] = body->getInertiaGC(gamma);
          F += std::array<double, 3>{mx * ACCELS_IN[i], my * ACCELS_IN[i + 1], mz * ACCELS_IN[i + 2]};
          T_GLOBAL += Dot(IG, std::array<double, 3>{ACCELS_IN[i + 3], ACCELS_IN[i + 4], ACCELS_IN[i + 5]});
        }

        const auto [mx, my, mz, Ix, Iy, Iz, IG, inv_IG] = body->getInertiaGC(1. + gamma);

        double a0 = F[0] / mx;
        double a1 = F[1] / my;
        double a2 = F[2] / mz;

        Tddd w = body->velocityRotational();
        Tddd Iw = Dot(IG, w);
        T_GLOBAL -= Cross(w, Iw);

        Tddd A_rotaton = Dot(inv_IG, T_GLOBAL);
        auto [a3, a4, a5] = A_rotaton;
        // std::ranges::for_each(T6d{a0, a1, a2, a3, a4, a5}, [&](const auto &a_w) { ACCELS[i++] = a_w; }); // 複数浮体がある場合があるので．
        // 閾値を大きく設定（1E20に対してはfalse、1E12程度の実船に対してはtrueになるように）
        const double threshold = 1E17;

        ACCELS[i] = std::abs(mx) < threshold ? a0 : 0;
        ACCELS[i + 1] = std::abs(my) < threshold ? a1 : 0;
        ACCELS[i + 2] = std::abs(mz) < threshold ? a2 : 0;
        ACCELS[i + 3] = std::abs(Ix) < threshold ? a3 : 0;
        ACCELS[i + 4] = std::abs(Iy) < threshold ? a4 : 0;
        ACCELS[i + 5] = std::abs(Iz) < threshold ? a5 : 0;
      } else
        i += 6;
    // std::cout << Green << "other" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
    ACCELS_OUT = ACCELS;
    std::cout << "acceleration: " << ACCELS_OUT << std::endl;
    return ACCELS - ACCELS_IN;
  };

  /* -------------------------------------------------------------------------- */

  V_d getVectorFromBoundary(const std::vector<Network*>& waters, int size, std::function<double(const networkPoint*, networkFace*)> getDirichletVal, std::function<double(const networkPoint*, networkFace*)> getNeumannVal) {
    V_d vec(size);
    for (const auto water : waters) {
      const auto& pts = water->getBoundaryPoints();
#pragma omp parallel for
      for (size_t k = 0; k < pts.size(); ++k) {
        const auto& p = pts[k];
        for (const auto& [f, i] : p->f2Index) {
          if (isDirichletID_BEM(p, f))
            vec[i] = getDirichletVal(p, f);
          else if (isNeumannID_BEM(p, f))
            vec[i] = getNeumannVal(p, f);
        }
      }
      // Include midpoint DOFs for true quadratic elements (per-face maps)
      if (use_true_quadratic_element) {
        for (auto* l : water->getBoundaryLines()) {
          for (const auto& [f, idx] : l->f2Index) {
            if (idx >= 0 && idx < size) {
              if (isDirichletID_BEM(l, f))
                vec[idx] = l->phinOnFace.at(f);
              else if (isNeumannID_BEM(l, f))
                vec[idx] = l->phiOnFace.at(f);
            }
          }
        }
      }
    }
    return vec;
  }

  //@ --------------------------------------------------- */
  //@        加速度 --> phiphin_t --> 圧力 --> 加速度     */
  //@ --------------------------------------------------- */
  std::vector<double> solveForPhiPhin_t(const std::vector<Network*>& rigidbodies) {
    std::vector<double> convergence;
    for (auto& water : WATERS)
      water->setGeometricPropertiesForce();

    // updateGeometricProperties();

    auto ACCELS_init = initializeAcceleration(rigidbodies);

    if (ACCELS_init.empty()) {
      return convergence;
    }

    int count = 0;
    double alpha = 1.;

    V_d ACCELS_OUT = ACCELS_init;

    auto normalize_flag = [](std::string s) {
      while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())))
        s.erase(s.begin());
      while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))
        s.pop_back();
      for (auto& ch : s)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
      return s;
    };
    const std::string coupling = normalize_flag(coupling_type);
    if (coupling.empty() || coupling == "none" || coupling == "off") {
      insertAcceleration(rigidbodies, ACCELS_init);
      auto func = Func(ACCELS_init, WATERS, rigidbodies, ACCELS_OUT, "uncoupled");
      convergence.push_back(Norm(func));
      return convergence;
    }

    auto solve_coupling = [&](auto& solver, const std::string& name, auto print_extra_info) {
      V_d current_accels = ACCELS_init;
      insertAcceleration(rigidbodies, current_accels);

      int seems_converged = 0;
      auto solver_tol_copy = solver_tol;
      for (auto j = 0; j < 50; ++j) {

        solver_tol = coupling_tol * 0.5; // 内部解法の許容誤差を厳しくする

        const std::string tag = name + " j=" + std::to_string(j) + " (residual)";
        auto func = Func(current_accels, WATERS, rigidbodies, ACCELS_OUT, tag.c_str());
        convergence.push_back(Norm(func));

        V_d next_accels = solver.compute_next(current_accels, func);
        V_d dX = next_accels - current_accels;

        std::cout << "[" << name << "] j = " << j;
        print_extra_info();
        std::cout << ", Norm(func) = " << Norm(func) << ", " << Red << "Norm(dX) = " << Norm(dX) << colorReset << std::endl;

        current_accels = next_accels;
        insertAcceleration(rigidbodies, current_accels);

        if (!isFinite(Norm(current_accels))) {
          current_accels = ACCELS_init;
          solver.reset();
          insertAcceleration(rigidbodies, current_accels);
          std::cout << "Re-initialize " << name << std::endl;
        }

        if (Norm(func) < coupling_tol && j > 1) {
          seems_converged++;
          if (seems_converged == 2)
            break;
        }
      }
      solver_tol = solver_tol_copy;
    };

    if (coupling == "anderson" || coupling == "aa") {
      int history_size = 10;
      if (!coupling_params.empty())
        history_size = static_cast<int>(coupling_params[0]);

      AndersonAcceleration<V_d> AA(history_size);
      solve_coupling(AA, "AA", []() {});
    } else if (coupling == "aitken") {
      double relaxation = 0.5;
      if (!coupling_params.empty())
        relaxation = coupling_params[0];

      AitkenAcceleration<V_d> Aitken(relaxation);
      solve_coupling(Aitken, "Aitken", [&]() { std::cout << ", omega = " << Aitken.omega; });
    } else if (coupling == "newton") {
      // Newton-Raphson method
      // Since the fluid force is linear w.r.t acceleration (added mass), the Jacobian is constant.
      // This is equivalent to calculating the added mass matrix and solving the linear system directly.
      auto solve_newton = [&](const std::string& name) {
        auto solver_tol_copy = solver_tol;
        solver_tol = coupling_tol * 0.1; // Tighten tolerance for Jacobian calculation

        // Cache Jacobian once per time step to avoid repeating FD GMRES solves.
        static bool jacobian_cached = false;
        static double jacobian_time = 0.0;
        static std::size_t jacobian_size = 0;
        static std::vector<V_d> jacobian;

        V_d current_accels = ACCELS_init;
        insertAcceleration(rigidbodies, current_accels);

        for (int iter = 0; iter < 20; ++iter) {
          // 1. Evaluate residual F(x)
          const std::string tag = name + " iter=" + std::to_string(iter) + " (residual)";
          auto func = Func(current_accels, WATERS, rigidbodies, ACCELS_OUT, tag.c_str());
          double res = Norm(func);
          convergence.push_back(res);
          std::cout << "[" << name << "] iter=" << iter << ", Norm(func) = " << res << std::endl;

          if (res < coupling_tol) {
            std::cout << "[" << name << "] Converged." << std::endl;
            break;
          }

          int n = current_accels.size();
          const bool need_jacobian = (!jacobian_cached) || (jacobian_time != simulation_time) || (jacobian_size != static_cast<std::size_t>(n));
          if (need_jacobian) {
            jacobian.assign(n, V_d(n));
            double epsilon = 1e-6;
            std::cout << Cyan << "  [" << name << "] Jacobian via FD: n=" << n << " eps=" << epsilon << " -> " << n << " extra GMRES solves" << colorReset << std::endl;

            // 2. Compute Jacobian J = dF/dx via finite differences (once per time step)
            for (int j = 0; j < n; ++j) {
              V_d perturbed = current_accels;
              perturbed[j] += epsilon;
              insertAcceleration(rigidbodies, perturbed);
              V_d dummy;
              const std::string tag_p = name + " iter=" + std::to_string(iter) + " (FD j=" + std::to_string(j) + "/" + std::to_string(n - 1) + ")";
              auto func_p = Func(perturbed, WATERS, rigidbodies, dummy, tag_p.c_str());
              for (int i = 0; i < n; ++i)
                jacobian[i][j] = (func_p[i] - func[i]) / epsilon;
            }
            jacobian_cached = true;
            jacobian_time = simulation_time;
            jacobian_size = static_cast<std::size_t>(n);
          }

          // 3. Solve J * dx = -F using a linear solver from basic_linear_systems.hpp
          V_d dx(n);
          V_d rhs = -func;
          lapack_lu lu_solver(jacobian);
          lu_solver.solve(rhs, dx);

          for (int i = 0; i < n; ++i)
            current_accels[i] += dx[i];
          insertAcceleration(rigidbodies, current_accels);

          std::cout << "[" << name << "] step, Norm(dX) = " << Norm(dx) << std::endl;
        }
        solver_tol = solver_tol_copy;
      };
      solve_newton("Newton");
    } else {
      struct BroydenWrapper {
        BroydenMethod<V_d> BM;
        V_d F_prev;
        bool first = true;
        double alpha_init = 0.1;
        double alpha_update = 1.0;

        BroydenWrapper(const V_d& init_val) {
          if (!coupling_params.empty()) {
            alpha_init = coupling_params[0];
            if (coupling_params.size() > 1)
              alpha_update = coupling_params[1];
          }
          BM.initialize(init_val, V_d(init_val.size(), 0.0));
        }

        void reset() {
          first = true;
          // Re-initialize to reset J and Inv_J to Identity
          BM.initialize(BM.X, V_d(BM.X.size(), 0.0));
        }

        V_d compute_next(const V_d& X, const V_d& F) {
          if (first) {
            first = false;
            F_prev = F;
            BM.X = X;
            // Initial step: X_{k+1} = X_k - alpha * F_k (assuming J=I)
            BM.dX = V_d(X.size());
            for (size_t i = 0; i < X.size(); ++i)
              BM.dX[i] = -alpha_init * F[i];
            BM.X += BM.dX;
            return BM.X;
          } else {
            // BM.dX already holds the step (X - X_prev) from the previous iteration/return
            // BM.updateGoodBroyden updates Inv_J using (dX, dF) and then computes new dX, X
            BM.updateGoodBroyden(F, F_prev, alpha_update);
            F_prev = F;
            return BM.X;
          }
        }
      };

      BroydenWrapper BM_Wrapper(ACCELS_init);
      solve_coupling(BM_Wrapper, "Broyden", []() {});
    }
    return convergence;
  };
};
