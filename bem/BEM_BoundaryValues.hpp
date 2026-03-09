#ifndef BEM_utilities_H
#define BEM_utilities_H

#include "BEM_WaveGenerator.hpp"
#include "Hadzic2005.hpp"
#include "Network.hpp"

using V_i = std::vector<int>;
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_Netp = std::vector<Network*>;
using V_netFp = std::vector<networkFace*>;
using VV_netFp = std::vector<V_netFp>;

/*

# BEMの計算における物理量の導出

具体的には、以下の3つの主要な目的を担っています。

1. 物体表面の速度・加速度の計算:
velocity_of_Body, accel_of_Body といった関数が、剛体・柔体を問わず、流体が接触している物体表面の正確な速度や加速度を計算します。これはNeumann境界条件 phi_n や phi_nt を決定する上で不可欠です。

2. Neumann境界条件の決定:
contactNormalVelocity, accelNeumann, contactPureVelocity などの関数群が、上記で計算された物体表面の速度・加速度を元に、流体節点における法線方向の境界条件を決定します。多重節点や複数の物体との接触も考慮する、複雑で重要なロジックです。

3. BVP（境界値問題）求解後の物理量導出:
gradPhi (流速の計算)、HessianOfPhi (ヘッセ行列の計算)、phint_Neumann (ポテンシャルの時間微分の法線成分) といった関数が、BVPを解いて得られたポテンシャル phi を使って、シミュレーションの次のステップに必要な物理量を計算します。
このように、ファイルの内容は**「BEMの境界条件と物理量の計算」**というテーマに沿って一貫性が保たれています。

*/

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_1_2_BOUNDARY_CONDITIONS

### `getContactFaces()`や`getNearestContactFace()`の利用

#### `contact_angle`と`isInContact()`

\insert{networkPoint::contact_angle}

#### `addContactFaces()`

\insert{networkPoint::addContactFaces()}

#### 呼び出し方法

\insert{networkPoint::getContactFaces()}

これらは，`contactNormalVelocity()`や`accelNeumann()`で利用される．

### `contactNormalVelocity()`と`accelNeumann()`

接触している物体が，剛体でない場合，
`velocity_of_Body`は，物体の節点（ `networkPoint` ）の速度（加速度）を元にして速度（加速度）を計算する．
そのため，`networkPoint::velocity`や`networkPoint::accel`を設定しておく必要がある．

`contactNormalVelocity(p, const adjacent_f)`や`accelNeumann(p, const adjacent_f)`
を使う時は，必ず`adjacent_f`が`p`に**隣接面するノイマン面**であることを確認する．

*/

/* -------------------------------------------------------------------------- */
// 接触面、または貫入している物体の最近傍面を取得するヘルパー関数
// networkPoint / networkLine 両対応（ContactDetectable 経由）
template <typename Node>
inline std::vector<networkFace*> getEffectiveContactFaces(const Node* p) {
  auto faces = p->getContactFaces();
  if (!faces.empty())
    return faces;

  // 接触面がないが、物体内部に貫入している場合は、その物体の最近傍面を返す
  if (p->penetratedBody) {
    auto [near_f, _] = p->penetratedBody->Nearest(p->getPosition());
    if (near_f)
      faces.push_back(near_f);
  }
  return faces;
}

/* -------------------------------------------------------------------------- */
// 接触面、または貫入している物体の最近傍面と最近傍点を取得するヘルパー関数
// networkPoint / networkLine 両対応（ContactDetectable 経由）
template <typename Node>
inline std::tuple<networkFace*, Tddd, double> getEffectiveNearestContactFace(const Node* p) {
  // 1. 既存の接触面リストから最近傍のものを取得
  auto [f, X, dist] = p->getNearestContactFace();
  if (f) {
    return {f, X, dist};
  }

  // 2. 接触面がないが、物体内部に貫入している場合は、その物体の最近傍面を返す
  if (p->penetratedBody) {
    auto [near_f, near_X] = p->penetratedBody->Nearest(p->getPosition());
    if (near_f) {
      return {near_f, near_X, Norm(p->getPosition() - near_X)};
    }
  }

  // 3. 接触も貫入もしていない
  return {nullptr, {}, 1E+20};
}

/* -------------------------------------------------------------------------- */

//$ --------------------------------------------------------------- */

template <typename Func>
inline Tddd interpolateOnFace(const networkFace* const face, const Tddd& X_contact, Func getter) {
  auto [p0, p1, p2] = face->getPoints();
  Tddd v0 = p1->X - p0->X, v1 = p2->X - p0->X, v2 = X_contact - p0->X;
  double d00 = Dot(v0, v0), d01 = Dot(v0, v1), d11 = Dot(v1, v1), d20 = Dot(v2, v0), d21 = Dot(v2, v1);
  double denom = d00 * d11 - d01 * d01;
  if (std::abs(denom) < 1e-20)
    return (getter(p0) + getter(p1) + getter(p2)) / 3.0;
  double v = (d11 * d20 - d01 * d21) / denom, w = (d00 * d21 - d01 * d20) / denom, u = 1.0 - v - w;
  return u * getter(p0) + v * getter(p1) + w * getter(p2);
}

inline Tddd velocity_of_Body(const networkFace* const face, const Tddd& X_contact) {
  if (face->getNetwork()->isRigidBody)
    return face->getNetwork()->velocityRigidBody(X_contact);
  if (face->getNetwork()->isSoftBody)
    return interpolateOnFace(face, X_contact, [](const networkPoint* p) { return p->velocityTranslational(); });
  return {0., 0., 0.};
};

inline Tddd accel_of_Body(const networkFace* const face, const Tddd& X_contact) {
  if (face->getNetwork()->isRigidBody)
    return face->getNetwork()->accelRigidBody(X_contact);
  if (face->getNetwork()->isSoftBody)
    return interpolateOnFace(face, X_contact, [](const networkPoint* p) { return p->accelTranslational(); });
  return {0., 0., 0.};
};

inline Tddd propertyNeumann(const networkPoint* const p, std::function<Tddd(const networkFace*, const Tddd&)> propertyFunc) {
  std::vector<Tddd> Directions;
  std::vector<double> Vsample;
  Tddd Vinit = {0., 0., 0.}, V;

  for (const auto& [f, contact_face_of_body_X] : p->getNearestContactFaces()) {
    auto [contact_face_of_body, X, _] = contact_face_of_body_X;
    if (contact_face_of_body) {
      V = propertyFunc(contact_face_of_body, X);
      Vsample.emplace_back(Dot(V, f->normal));
      Vinit += V;
      Directions.emplace_back(f->normal);
    }
  }
  if (!Vsample.empty()) {
    Vinit /= Vsample.size();
    auto ret = optimalVector(Vsample, Directions, Vinit);
    if (isFinite(ret))
      return ret;
  }
  return Vinit;
};

inline Tddd contactNormalVelocity(const networkPoint* const p) { return propertyNeumann(p, velocity_of_Body); };
inline Tddd accelNeumann(const networkPoint* const p) { return propertyNeumann(p, accel_of_Body); };

inline Tddd propertyNeumann(const networkPoint* const p, const networkFace* const adjacent_f, std::function<Tddd(const networkFace*, const Tddd&)> propertyFunc) {
  if (!adjacent_f) {
    return propertyNeumann(p, propertyFunc);
  }
  auto [contact_face_of_body, X_contact, dist] = p->getNearestContactFace_(adjacent_f);
  return contact_face_of_body ? propertyFunc(contact_face_of_body, X_contact) : Tddd{0., 0., 0.};
};

inline Tddd contactNormalVelocity(const networkPoint* const p, const networkFace* const adjacent_f) { return propertyNeumann(p, adjacent_f, velocity_of_Body); };
inline Tddd accelNeumann(const networkPoint* const p, const networkFace* const adjacent_f) { return propertyNeumann(p, adjacent_f, accel_of_Body); };

// --- networkLine overloads (midpoint uses ContactDetectable base class) ---

inline Tddd propertyNeumann(const networkLine* const l, std::function<Tddd(const networkFace*, const Tddd&)> propertyFunc) {
  std::vector<Tddd> Directions;
  std::vector<double> Vsample;
  Tddd Vinit = {0., 0., 0.}, V;
  for (const auto& [f, contact_face_of_body_X] : l->getNearestContactFaces()) {
    auto [contact_face_of_body, X, _] = contact_face_of_body_X;
    if (contact_face_of_body) {
      V = propertyFunc(contact_face_of_body, X);
      Vsample.emplace_back(Dot(V, f->normal));
      Vinit += V;
      Directions.emplace_back(f->normal);
    }
  }
  if (!Vsample.empty()) {
    Vinit /= Vsample.size();
    auto ret = optimalVector(Vsample, Directions, Vinit);
    if (isFinite(ret))
      return ret;
  }
  return Vinit;
};

inline Tddd propertyNeumann(const networkLine* const l, const networkFace* const adjacent_f, std::function<Tddd(const networkFace*, const Tddd&)> propertyFunc) {
  if (!adjacent_f)
    return propertyNeumann(l, propertyFunc);
  auto [contact_face_of_body, X_contact, dist] = l->getNearestContactFace_(adjacent_f);
  return contact_face_of_body ? propertyFunc(contact_face_of_body, X_contact) : Tddd{0., 0., 0.};
};

inline Tddd contactNormalVelocity(const networkLine* const l) { return propertyNeumann(l, velocity_of_Body); };
inline Tddd contactNormalVelocity(const networkLine* const l, const networkFace* const adjacent_f) { return propertyNeumann(l, adjacent_f, velocity_of_Body); };

inline Tddd getNormalNeumann_mid(const networkLine* const l) {
  Tddd normal = {0., 0., 0.};
  for (auto* f : l->getBoundaryFaces())
    if (f && f->Neumann)
      normal += f->area * f->normal;
  return Normalize(normal);
};

template <typename BodyFunc>
inline Tddd contactPureProperty(const networkPoint* const p, BodyFunc bodyFunc) {
  V_d Vsample;
  std::vector<Tddd> Directions;
  Tddd Vinit = {0., 0., 0.}, u, ex = {1., 0., 0.}, ey = {0., 1., 0.}, ez = {0., 0., 1.};
  auto nearestFaces = p->getNearestContactFaces();

  if (!nearestFaces.empty()) {
    for (const auto& [f, contact_face_of_body_X] : nearestFaces) {
      auto [contact_face_of_body, X, _] = contact_face_of_body_X;
      if (contact_face_of_body) {
        u = bodyFunc(contact_face_of_body, X);
        Vsample.insert(Vsample.end(), {Dot(u, ex), Dot(u, ey), Dot(u, ez)});
        Directions.insert(Directions.end(), {ex, ey, ez});
      }
    }
  } else {
    // 変更箇所: ヘルパーを利用
    auto [near_f, near_X, _] = getEffectiveNearestContactFace(p);
    if (near_f) {
      u = bodyFunc(near_f, near_X);
      Vsample.insert(Vsample.end(), {Dot(u, ex), Dot(u, ey), Dot(u, ez)});
      Directions.insert(Directions.end(), {ex, ey, ez});
    }
  }
  if (!Vsample.empty())
    return optimalVector(Vsample, Directions, Vinit);
  else
    return Vinit;
}

inline Tddd contactPureVelocity(const networkPoint* const p) { return contactPureProperty(p, velocity_of_Body); }
inline Tddd contactPureAccel(const networkPoint* const p) { return contactPureProperty(p, accel_of_Body); }

template <typename BodyFunc, typename PureFunc>
inline Tddd contactPureProperty(const networkPoint* const p, const networkFace* const adjacent_f, BodyFunc bodyFunc, PureFunc pureFunc) {
  if (!adjacent_f)
    return pureFunc(p); // fallback
  auto [contact_face_of_body, X_contact, _] = p->getNearestContactFace_(adjacent_f);
  if (contact_face_of_body)
    return bodyFunc(contact_face_of_body, X_contact);
  if (adjacent_f->penetratedBody) {
    auto [near_f, near_X] = adjacent_f->penetratedBody->Nearest(getPosition(p));
    if (near_f)
      return bodyFunc(near_f, near_X);
  }
  return pureFunc(p); // fallback
}

inline Tddd contactPureVelocity(const networkPoint* const p, const networkFace* const adjacent_f) {
  return contactPureProperty(p, adjacent_f, velocity_of_Body, [](const networkPoint* p) { return contactPureVelocity(p); });
}
inline Tddd contactPureAccel(const networkPoint* const p, const networkFace* const adjacent_f) {
  return contactPureProperty(p, adjacent_f, accel_of_Body, [](const networkPoint* p) { return contactPureAccel(p); });
}

//$ --------------------------------------------------------------- */

using map_P_d = std::map<netP*, double>;
using map_P_Vd = std::map<netP*, V_d>;
using map_P_VVd = std::map<netP*, VV_d>;
using map_F_P_Vd = std::map<netF*, map_P_Vd>;
using map_P_P_Vd = std::map<netP*, map_P_Vd>;
using pair_PB = std::pair<netP*, bool>;
using map_pairPB_Tdd = std::unordered_map<std::tuple<netP*, bool, netF*>, Tdd>;
using map_pairPB_pairPB_Tdd = std::unordered_map<std::tuple<netP*, bool, netF*>, map_pairPB_Tdd>;
using map_P_P_Tdd = std::map<netP*, std::map<netP*, Tdd>>;
using map_P_F_P_Vd = std::map<netP*, map_F_P_Vd>;
using VV_SorIorMap = std::vector<std::vector<std::variant<std::string, int, map_P_Vd>>>;

inline V_netFp takeFaces(const V_Netp& nets) {
  V_netFp ret({});
  for (const auto& n : nets)
    ret.insert(ret.end(), n->getBoundaryFaces().begin(), n->getBoundaryFaces().end());
  return DeleteDuplicates(ret);
};

inline Tddd grad_phi_tangential(const networkFace* const f) {
  auto [p0, p1, p2] = f->getPoints();
  Tddd PHI = {std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)};
  return gradP1(T3Tddd{{p0->X, p1->X, p2->X}}, PHI);
};

// T3Tddd grad_LinearElement(const T3Tddd &F012, const T3Tddd &X012, const T3Tddd &F_n) {
//    //! 三角要素の節点の情報変数F0,F1,F2から，三角要素上でのgrad(F)を計算する．
//    return {grad_LinearElement(std::get<0>(F012), X012, std::get<0>(F_n)),
//            grad_LinearElement(std::get<1>(F012), X012, std::get<1>(F_n)),
//            grad_LinearElement(std::get<2>(F012), X012, std::get<2>(F_n))};
// };

inline T3Tddd OrthogonalBasis(const Tddd& n_IN) {
  auto n = Normalize(n_IN);
  Tddd s0 = Chop(Tddd{1, 0, 0}, n);
  if (Norm(s0) < 1E-3)
    s0 = Chop(Tddd{0, 1, 0}, n);
  s0 = Normalize(s0);
  Tddd s1 = Normalize(Cross(n, s0));
  return {n, s0, s1};
};

// Unified getPhi/getPhin: works for both networkPoint* and networkLine*
template <typename Node>
inline double getPhi(const Node* n, const networkFace* f) {
  auto iter = n->phiOnFace.find(const_cast<networkFace*>(f));
  if (iter != n->phiOnFace.end())
    return iter->second;
  if (n->phiOnFace.find(nullptr) != n->phiOnFace.end())
    return n->phiOnFace.at(nullptr);
  return std::get<0>(n->phiphin);
}

template <typename Node>
inline double getPhin(const Node* n, const networkFace* f) {
  auto iter = n->phinOnFace.find(const_cast<networkFace*>(f));
  if (iter != n->phinOnFace.end())
    return iter->second;
  if (n->phinOnFace.find(nullptr) != n->phinOnFace.end())
    return n->phinOnFace.at(nullptr);
  return std::get<1>(n->phiphin);
}

inline Tddd gradPhiQuadElement(const networkPoint* p, networkFace* f) {
  try {
    //* p will be set as node 4

    // DodecaPoints dodecapoint(f, p, [](const networkLine *line) -> bool { return !line->CORNER; });

    int index = -1;
    for (auto i = 0; i < 3; ++i)
      if (f->Points[i] == p) {
        index = i;
        break;
      }
    if (index == -1)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error");

    auto ToPhi = [&](const networkPoint* p) -> double { return std::get<0>(p->phiphin); };
    auto ToX = [&](const networkPoint* p) -> Tddd { return p->X; };

    const double phi_t0 = f->dodecaPoints[index]->D_interpolate<1, 0>(1., 0., ToPhi); //! at 4
    const double phi_t1 = f->dodecaPoints[index]->D_interpolate<0, 1>(1., 0., ToPhi); //! at 4
    const double phi_n = getPhin(p, f);

    const Tddd dX_t0 = f->dodecaPoints[index]->D_interpolate<1, 0>(1., 0., ToX); //! at 4
    const Tddd dX_t1 = f->dodecaPoints[index]->D_interpolate<0, 1>(1., 0., ToX); //! at 4
    const auto Nxyz = Normalize(Cross(dX_t0, dX_t1));

    Tddd grad_phi;
    lapack_svd_solve(T3Tddd{dX_t0, dX_t1, Nxyz}, grad_phi, Tddd{phi_t0, phi_t1, phi_n});
    return grad_phi;

    /*check!

    `Dot[{{Sx, Sy, Sz}, {T0, T1, T2}, {N0, N1, N2}}, DPHI]`この計算は，`{Sx, Sy, Sz}`などの成分を取り出す計算．

    ```Mathematica
    DPHI = {phix, phiy, phiz}
    Dot[{{Sx, Sy, Sz}, {T0, T1, T2}, {N0, N1, N2}}, DPHI]
    Dot[{Sx, Sy, Sz}, DPHI]
    ```

    */
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in gradPhiQuadElement");
  }
};

/* -------------------------------------------------------------------------- */

// Compute gradient of phi at vertex p on face f using 6-node quadratic (TriShape<6>) interpolation.
// Parametric coordinates of a vertex on a face
// vertex 0 -> (1,0), vertex 1 -> (0,1), vertex 2 -> (0,0)
inline std::pair<double, double> getParametricCoords(const networkPoint* p, const networkFace* f) {
  constexpr double t0[3] = {1., 0., 0.};
  constexpr double t1[3] = {0., 1., 0.};
  for (int i = 0; i < 3; ++i)
    if (f->Points[i] == p)
      return {t0[i], t1[i]};
  return {-1., -1.};
}

// Parametric coordinates of an edge midpoint on a face
// l0 (p0-p1) -> (0.5,0.5), l1 (p1-p2) -> (0,0.5), l2 (p2-p0) -> (0.5,0)
inline std::pair<double, double> getParametricCoords(const networkLine* l, const networkFace* f) {
  constexpr double t0[3] = {0.5, 0.0, 0.5};
  constexpr double t1[3] = {0.5, 0.5, 0.0};
  for (int i = 0; i < 3; ++i)
    if (f->Lines[i] == l)
      return {t0[i], t1[i]};
  return {-1., -1.};
}

// Compute grad(phi) at a node (vertex or edge midpoint) on face f
// using 6-node quadratic shape functions D_TriShape<6>.
// Solves: [dX/dt0; dX/dt1; n] . grad_phi = [dphi/dt0; dphi/dt1; phin]
template <typename Node>
inline Tddd gradPhiTrueQuadElement(const Node* node, const networkFace* f) {
  auto [t0, t1] = getParametricCoords(node, f);
  if (t0 < 0.)
    return {0., 0., 0.};

  const auto [p0, l0, p1, l1, p2, l2] = f->PLPLPL;
  const std::array<double, 6> phi6 = {getPhi(p0, f), getPhi(p1, f), getPhi(p2, f), getPhi(l0, f), getPhi(l1, f), getPhi(l2, f)};
  const auto dN_dt0 = D_TriShape<6, 1, 0>(t0, t1);
  const auto dN_dt1 = D_TriShape<6, 0, 1>(t0, t1);

  double dphi_dt0 = 0., dphi_dt1 = 0.;
  for (int j = 0; j < 6; ++j) {
    dphi_dt0 += dN_dt0[j] * phi6[j];
    dphi_dt1 += dN_dt1[j] * phi6[j];
  }

  // Jacobian vectors (6-node quadratic geometry)
  const std::array<Tddd, 6> X6 = {p0->X, p1->X, p2->X, l0->X_mid, l1->X_mid, l2->X_mid};
  const Tddd dX_dt0 = Dot(dN_dt0, X6);
  const Tddd dX_dt1 = Dot(dN_dt1, X6);
  const Tddd n_local = Normalize(Cross(dX_dt0, dX_dt1));

  Tddd grad_phi;
  lapack_svd_solve(T3Tddd{dX_dt0, dX_dt1, n_local}, grad_phi, Tddd{dphi_dt0, dphi_dt1, getPhin(node, f)});
  return grad_phi;
}

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_3_0_INITIAL_VALUE_PROBLEM

## 初期値問題

節点の位置と速度ポテンシャル$\phi$に関する初期値問題を解いて行くことが，シミュレーションである．
言い換えると，節点位置$\frac{d\bf x}{dt}$と速度ポテンシャル$\frac{d\phi}{dt}$を少しずつ$\Delta t$ずつ時間積分することが，シミュレーションである．
ちなみに，$\frac{d\bf x}{dt}$や$\frac{d\phi}{dt}$を計算するには，境界値問題を解く必要がある．

ある時刻において，境界値問題が解けたら，$\frac{d\bf x}{dt}$と$\frac{d\phi}{dt}$はどのように計算できるだろうか．

### 流速$\frac{d\bf x}{dt}$の計算

ある三角要素上の接線流速$\nabla \phi_{\parallel}$は，線形三角要素補間を使って次のように計算する．

$$
\nabla \phi_{\parallel} = \frac{\bf n}{2A} \times (({\bf x}_2 - {\bf x}_1) \phi_0 +({\bf x}_0 - {\bf x}_2) \phi_1 + ({\bf x}_1 - {\bf x}_0) \phi_2)
\\= \frac{\bf n}{2A} \times (({\bf x}_0,{\bf x}_1,{\bf x}_2)\cdot(\phi_1-\phi_2,\phi_2-\phi_0,\phi_0-\phi_1))
$$

三角要素上の流速$\nabla \phi$は，次のように計算する．

$$
\nabla \phi = \frac{(\phi_n)_0+(\phi_n)_1+(\phi_n)_2}{3} {\bf n} + \nabla \phi_{\parallel}
$$

### $\frac{d\phi}{dt}$の計算

ある流体粒子に乗ってみたときの，速度ポテンシャルの時間変化$\frac{D \phi}{D t}$は，次のように計算できる．

$$
\frac{D \phi}{D t} = \frac{\partial \phi}{\partial t} + \nabla \phi \cdot \nabla \phi
$$

<details style="background-color: rgba(144, 238, 144, 0.2);">
<summary>
NOTE: オイラー的記述
</summary>

$\phi=\phi(t,{\bf x})$のように書き表し，位置と空間を独立させ分けて考える方法を，オイラー的記述という．こう書くと，$\frac{d \phi}{d t}$は，$\frac{\partial \phi}{\partial t}$であり，これは，速度ポテンシャルの純粋な時間変化ではない．純粋な，ある流体粒子の速度ポテンシャルの時間変化を表すためには，位置が時間によって変わると考え，つまり$\phi=\phi(t,{\bf x}(t))$と一時的に考えなおし，そして，時間微分する．そうすると$\frac{d\phi}{dt} = \frac{\partial \phi}{\partial t} + \frac{d\bf x}{dt}\cdot \nabla
\phi$となる．

</details>

ここの$\frac{\partial \phi}{\partial t}$の計算は簡単ではない．そこで，ベルヌーイの式（大気圧と接する水面におけるベルヌーイの式は圧力を含まず簡単）を使って，$\frac{\partial \phi}{\partial t}$を消去する．

*/

inline Tddd gradPhi(const networkFace* const f) {
  double phi_n = 0;
  for (const auto& p : f->getPoints())
    phi_n += getPhin(p, f);
  auto [p0, p1, p2] = f->getPoints();
  Tddd PHI = {std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)};
  return phi_n / 3. * f->normal + gradP1(T3Tddd{p0->X, p1->X, p2->X}, PHI);
}

//! use simple gradient method but Newton method
inline Tddd gradPhi(const networkPoint* const p, std::array<double, 3>& convergence_info) {

  Tddd u;
  const Tddd ex = {1., 0., 0.}, ey = {0., 1., 0.}, ez = {0., 0., 1.};
  const auto faces = p->getBoundaryFaces();
  auto s = faces.size();
  if (s == 0)
    return {0., 0., 0.};
  thread_local V_Tddd Directions;
  Directions.clear();
  thread_local V_d W, Vsample;
  W.clear();
  Vsample.clear();
  for (const auto& f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    Tddd PHI = {std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin)};
    Tddd u;
    if (f->isPseudoQuadraticElement)
      u = gradPhiQuadElement(p, f);
    else if (f->isTrueQuadraticElement)
      u = gradPhiTrueQuadElement(p, f);
    else
      u = gradP1(T3Tddd{{p0->X, p1->X, p2->X}}, PHI) + getPhin(p, f) * f->normal;

    Vsample.emplace_back(Dot(u, ex));
    Directions.emplace_back(ex);
    W.emplace_back(f->area);
    Vsample.emplace_back(Dot(u, ey));
    Directions.emplace_back(ey);
    W.emplace_back(f->area);
    Vsample.emplace_back(Dot(u, ez));
    Directions.emplace_back(ez);
    W.emplace_back(f->area);
  }

  double meanW = std::accumulate(W.begin(), W.end(), 0.) / W.size();

  /*
  このようなCORNERが必要なのは，水面の喫水線において，接線方向に流れが大きくなり，構造物にめり込んでしまうことが生じるため．
  構造物にはめり込まないような，phinが与えられているはずだが，最小値問題において，めり込む方が最小になるのだろう．
  以下を加えれば，めり込まない流れの方が最小となりやすくはなるだろう，
  接線流速の精度が良くないとい，ということもできるだろう．

  ->
  2025/01/19
  下はいらない．メッシュ解像度を上げることで，不安定はおさっまった．
  下をつけると，運動しないが，phinを与えたい境界条件のphinが０になる．
  */

  if (p->CORNER) {
    for (const auto& f : faces) {
      if (!f || !f->Neumann)
        continue;
      const Tddd n = f->normal;
      if (!isFinite(n) || !(Norm(n) > 0.))
        continue;
      Vsample.emplace_back(getPhin(p, f));
      Directions.emplace_back(n);
      W.emplace_back(10. * meanW);
    }
  }

  return optimalVector(Vsample, Directions, Tddd{0., 0., 0.}, W, convergence_info);
};

inline Tddd gradPhi(const networkPoint* const p) {
  std::array<double, 3> convergence_info;
  return gradPhi(p, convergence_info);
};

// Edge-midpoint version: weighted average of per-face gradients at the midpoint
// (analogous to the vertex gradPhi above).
inline Tddd gradPhi(const networkLine* const l) {
  const auto faces = l->getBoundaryFaces();
  if (faces.empty())
    return {0., 0., 0.};

  const Tddd ex = {1., 0., 0.}, ey = {0., 1., 0.}, ez = {0., 0., 1.};
  thread_local V_Tddd Directions;
  Directions.clear();
  thread_local V_d W, Vsample;
  W.clear();
  Vsample.clear();

  for (const auto& f : faces) {
    Tddd u;
    if (f->isTrueQuadraticElement)
      u = gradPhiTrueQuadElement(l, f);
    else {
      // Fallback: linear average of endpoint gradients
      auto [pA, pB] = l->getPoints();
      u = 0.5 * (gradPhi(pA) + gradPhi(pB));
    }
    Vsample.emplace_back(Dot(u, ex));
    Directions.emplace_back(ex);
    W.emplace_back(f->area);
    Vsample.emplace_back(Dot(u, ey));
    Directions.emplace_back(ey);
    W.emplace_back(f->area);
    Vsample.emplace_back(Dot(u, ez));
    Directions.emplace_back(ez);
    W.emplace_back(f->area);
  }

  double meanW = std::accumulate(W.begin(), W.end(), 0.) / W.size();
  if (l->CORNER) {
    for (const auto& f : faces) {
      if (!f || !f->Neumann)
        continue;
      const Tddd n = f->normal;
      if (!isFinite(n) || !(Norm(n) > 0.))
        continue;
      Vsample.emplace_back(getPhin(l, f));
      Directions.emplace_back(n);
      W.emplace_back(10. * meanW);
    }
  }

  std::array<double, 3> convergence_info;
  return optimalVector(Vsample, Directions, Tddd{0., 0., 0.}, W, convergence_info);
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

#### $\phi$のヘッセ行列の計算

$$
\nabla\otimes{\bf u} = \nabla \otimes \nabla \phi =
\begin{bmatrix} \phi_{xx} & \phi_{xy} & \phi_{xz} \\
　　　　　　　　　　\phi_{yx} & \phi_{yy} & \phi_{yz} \\
　　　　　　　　　　\phi_{zx} & \phi_{zy} & \phi_{zz}
\end{bmatrix}
$$

ヘッセ行列の計算には，要素における変数の勾配の接線成分を計算する\ref{BEM:HessianOfPhi}{`HessianOfPhi`}を用いる．
節点における変数を$v$とすると，$\nabla v-{\bf n}({\bf n}\cdot\nabla v)$が計算できる．
要素の法線方向${\bf n}$が$x$軸方向${(1,0,0)}$である場合，$\nabla v - (\frac{\partial}{\partial x},0,0)v$なので，
$(0,\frac{\partial v}{\partial y},\frac{\partial v}{\partial z})$が得られる．
ただし，これは位置座標の基底を変えた後で使用する．

*/

/* -------------------------------------------------------------------------- */
// \label{BEM:HessianOfPhi}
inline T3Tddd HessianOfPhi(auto F, const T3Tddd& basis) {
  //! 位置座標変換
  auto [P0, P1, P2] = F->getPoints();
  auto X012_for_s0s1s2 = T3Tddd{Dot(basis, P0->X), Dot(basis, P1->X), Dot(basis, P2->X)};
  //! 速度ポテンシャルの勾配の座標変換
  auto [g0_s0, g0_s1, g0_s2] = Dot(basis, P0->u_potential_BEM);
  auto [g1_s0, g1_s1, g1_s2] = Dot(basis, P1->u_potential_BEM);
  auto [g2_s0, g2_s1, g2_s2] = Dot(basis, P2->u_potential_BEM);
  // auto [g_s0s0, g_s0s1, g_s0s2] = gradTangential_LinearElement(Tddd{getPhin(P0, F), getPhin(P1, F), getPhin(P2, F)}, X012);
  auto [g_s0s0, g_s0s1, g_s0s2] = gradP1(X012_for_s0s1s2, Tddd{g0_s0, g1_s0, g2_s0}); // 速度勾配テンソル　\nabla \otimes u or \nabla u or d/dx_j u_i
  auto [g_s1s0, g_s1s1, g_s1s2] = gradP1(X012_for_s0s1s2, Tddd{g0_s1, g1_s1, g2_s1}); // 速度勾配テンソル
  auto [g_s2s0, g_s2s1, g_s2s2] = gradP1(X012_for_s0s1s2, Tddd{g0_s2, g1_s2, g2_s2}); // 速度勾配テンソル
  // return T3Tddd{{{-g_s1s1 - g_s2s2, g_s0s1 /*wont be used*/, g_s0s2 /*wont be used*/},
  //                {g_s1s0, g_s1s1 /*wont be used*/, g_s1s2 /*wont be used*/},
  //                {g_s2s0, g_s2s1 /*wont be used*/, g_s2s2 /*wont be used*/}}};
  return T3Tddd{{{-g_s1s1 - g_s2s2, (g_s0s1 + g_s1s0) * 0.5 /*wont be used*/, (g_s0s2 + g_s2s0) * 0.5 /*wont be used*/}, {(g_s0s1 + g_s1s0) * 0.5, g_s1s1 /*wont be used*/, (g_s2s1 + g_s1s2) * 0.5 /*wont be used*/}, {(g_s0s2 + g_s2s0) * 0.5, (g_s2s1 + g_s1s2) * 0.5 /*wont be used*/, g_s2s2 /*wont be used*/}}};
  // return T3Tddd{{{-g_s1s1 - g_s2s2, g_s0s1 /*wont be used*/, g_s0s2 /*wont be used*/},
  //                {g_s0s1, g_s1s1 /*wont be used*/, g_s1s2 /*wont be used*/},
  //                {g_s0s2, g_s2s1 /*wont be used*/, g_s2s2 /*wont be used*/}}};
};

/*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

### $\phi_{nt}$の計算で必要となる${\bf n}\cdot \left({\frac{d\mathbfit r}{dt}  \cdot \nabla\otimes\nabla \phi}\right)$について．

$\nabla$を，$(x,y,z)$の座標系ではなく，
面の法線方向${\bf n}$を$x$の代わりにとり，
面に水平な方向を$t_0,t_1$とする座標系で考えることにして，$\nabla^*$と書くことにする．
${\bf n}\cdot \left({\frac{d\mathbfit r}{dt}  \cdot \nabla\otimes\nabla \phi}\right)$では，${\bf n}$方向成分だけをとる操作をしているので，
新しい座標系でも同じようにすれば，結果は変わらない．

$$
{\bf n}\cdot \left({\frac{d\mathbfit r}{dt}  \cdot \nabla\otimes\nabla \phi}\right) =  {(1,0,0)}\cdot\left({\frac{d{\mathbfit r}^\ast}{dt} \cdot \nabla^* \otimes\nabla^* \phi}\right).
\quad \nabla^* \otimes\nabla^* \phi =
\begin{bmatrix}
\phi_{nn} & \phi_{nt_0} & \phi_{nt_1} \\
\phi_{t_0n} & \phi_{t_0t_0} & \phi_{t_0t_1} \\
\phi_{t_1n} & \phi_{t_1t_0} & \phi_{t_1t_1}
\end{bmatrix}
$$

最後に第１成分だけが残るので，

$$
{(1,0,0)}\cdot\left({\frac{d{\mathbfit r}^\ast}{dt}  \cdot \nabla^* \otimes\nabla^* \phi}\right) = \frac{d{\mathbfit r}^\ast}{dt} \cdot (\phi_{nn}, \phi_{t_0n}, \phi_{t_1n})
$$

$\phi_{nn}$は，直接計算できないが，ラプラス方程式から$\phi_{nn}=- \phi_{t_0t_0}- \phi_{t_1t_1}$となるので，水平方向の勾配の計算から求められる．

*/

// \label{BEM:phint_Neumann}

inline double phint_Neumann(const networkPoint* const p, networkFace* F) {
  try {
    //$ faceがNeumannである条件は，faceの持つpointがすべて，外部の面と接触している場合である．
    //$ なので，{p,f}は，かならずp->getNearestContactFace(F)を持つ．
    auto structure_f = p->getNearestContactFace(F);
    Network* net_ref = nullptr;
    if (!structure_f)
      net_ref = F->penetratedBody;
    else
      net_ref = structure_f->getNetwork();

    if (!net_ref)
      return 0.;

    Tddd Omega = (net_ref)->velocityRotational();
    Tddd n = F->normal;
    auto dndt = Cross(Omega, n);
    auto drdt = contactPureVelocity(p, F);
    auto dr2dt2 = contactPureAccel(p, F);
    auto basis = OrthogonalBasis(n);
    Tddd tmp = {Dot(basis[0], drdt), Dot(basis[1], drdt), Dot(basis[2], drdt)};
    auto phint = Dot(n, dr2dt2) + Dot(dndt, drdt - p->u_total) - Dot(tmp, HessianOfPhi(F, basis))[0];
    return phint;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 0.;
  }
};

// double phint_Neumann(const networkPoint *const p, networkFace *F) {
//    //$ faceがNeumannである条件は，faceの持つpointがすべて，外部の面と接触している場合である．
//    //$ なので，{p,f}は，かならずp->getNearestContactFace(F)を持つ．
//    auto f = p->getNearestContactFace(F);
//    if (f) {
//       Tddd Omega = (f->getNetwork())->velocityRotational();
//       auto grad_phi = gradPhi(F);
//       // auto grad_phi = gradPhi(p, F);
//       auto U_body = contactPureVelocity(p, F);
//       auto dndt = Cross(Omega, F->normal);
//       auto ret = Dot(dndt, U_body - grad_phi);
//       ret += Dot(F->normal, accelNeumann(p, F));
//       auto basis = OrthogonalBasis(F->normal);
//       // ret -= Dot(Dot(basis, F->normal) /*=(1,0,0)*/, Dot(Dot(basis, U_body), HessianOfPhi(F, basis)));
//       ret -= Dot(Dot(basis, U_body), HessianOfPhi(F, basis))[0];
//       // ret -= Dot(Dot(basis, F->normal) /*=(1,0,0)*/, Dot(Dot(basis, grad_phi), HessianOfPhi(F, basis)));
//       return ret;
//    } else
//       throw std::runtime_error("p is not in contact with F");
// };

inline double phint_Neumann(const networkPoint* const p) {
  // V_d Phin, W;
  // std::vector<Tddd> Direcctions;
  // double total = 0;
  // Tddd normal = {0., 0., 0.};
  // Phin.reserve(10);
  // W.reserve(10);
  // Direcctions.reserve(10);
  // for (const auto &f : p->getBoundaryFaces())
  //    if (f->Neumann) {
  //       Phin.emplace_back(phint_Neumann(p, f));
  //       Direcctions.emplace_back(f->normal);
  //       W.emplace_back(f->area);
  //       normal += f->normal * f->area;
  //       total += f->area;
  //    }

  // return Dot(optimalVector(Phin, Direcctions, {0., 0., 0.}, W), normal / total);

  double total = 0, phin = 0;
  for (const auto& f : p->getBoundaryFaces())
    if (f->Neumann) {
      phin += phint_Neumann(p, f) * f->area;
      total += f->area;
    }
  return total ? (phin / total) : 0.;
};

// --- networkLine overloads of phint_Neumann ---

inline double phint_Neumann(const networkLine* const l, networkFace* F) {
  try {
    auto structure_f = l->getNearestContactFace(F);
    Network* net_ref = nullptr;
    if (!structure_f)
      net_ref = F->penetratedBody;
    else
      net_ref = structure_f->getNetwork();

    if (!net_ref)
      return 0.;

    Tddd Omega = net_ref->velocityRotational();
    Tddd n = F->normal;
    auto dndt = Cross(Omega, n);
    auto drdt = contactNormalVelocity(l, F);
    auto dr2dt2 = propertyNeumann(l, F, accel_of_Body);
    auto basis = OrthogonalBasis(n);
    Tddd tmp = {Dot(basis[0], drdt), Dot(basis[1], drdt), Dot(basis[2], drdt)};
    return Dot(n, dr2dt2) + Dot(dndt, drdt - l->u_total) - Dot(tmp, HessianOfPhi(F, basis))[0];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 0.;
  }
};

inline double phint_Neumann(const networkLine* const l) {
  double total = 0, phint_val = 0;
  for (auto* f : l->getBoundaryFaces())
    if (f && f->Neumann) {
      phint_val += phint_Neumann(l, f) * f->area;
      total += f->area;
    }
  return total ? (phint_val / total) : 0.;
};

#endif
