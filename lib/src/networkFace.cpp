// #pragma once
#include "Network.hpp"
#include "pch.hpp"

bool networkFace::replace(netP* const oldP, netP* const newP) {
  bool found = false;
  if (std::get<0>(this->Points) == oldP) {
    std::get<0>(this->Points) = newP;
    found = true;
  } else if (std::get<1>(this->Points) == oldP) {
    std::get<1>(this->Points) = newP;
    found = true;
  } else if (std::get<2>(this->Points) == oldP) {
    std::get<2>(this->Points) = newP;
    found = true;
  }
  if (found) {
    if (std::get<0>(this->PLPLPL) == oldP)
      std::get<0>(this->PLPLPL) = newP;
    else if (std::get<2>(this->PLPLPL) == oldP)
      std::get<2>(this->PLPLPL) = newP;
    else if (std::get<4>(this->PLPLPL) == oldP)
      std::get<4>(this->PLPLPL) = newP;

    // 面の全頂点を invalidate（法線・面積が変わるため）
    auto [p0, p1, p2] = this->Points;
    if (p0) p0->geom_curvature.valid = false;
    if (p1) p1->geom_curvature.valid = false;
    if (p2) p2->geom_curvature.valid = false;
    if (oldP) oldP->geom_curvature.valid = false;
  }
  return found;
};

namespace {

bool isPointFaceNeumannLocal(const networkPoint* p, const networkFace* f) {
  if (!p || !f)
    return false;
  const auto* d = p->findContactState(f);
  if (d && d->detached_by_pressure)
    return false;
  return d && d->nearestContactFace() != nullptr;
}

bool isFaceNeumannFromPointFaceStatesLocal(const networkFace* f) {
  if (!f)
    return false;
  const auto pts = f->getPoints();
  return std::ranges::all_of(pts, [&](const auto* p) {
    return isPointFaceNeumannLocal(p, f);
  });
}

} // namespace

void networkFace::setDodecaPoints() {
  try {
    // setIntegrationInfo用に統一した条件（useOppositeFace）を使用
    auto condition = [](const networkLine *line) -> bool { return useOppositeFace(line, M_PI / 3); };
    std::get<0>(this->dodecaPoints) = std::make_shared<DodecaPoints>(this, std::get<0>(this->Points), condition);
    std::get<1>(this->dodecaPoints) = std::make_shared<DodecaPoints>(this, std::get<1>(this->Points), condition);
    std::get<2>(this->dodecaPoints) = std::make_shared<DodecaPoints>(this, std::get<2>(this->Points), condition);
  } catch (...) {
    // OMP 並列内で throw は禁止。エラーを無視。
  }
}

/*
   原点を変えて積分する際に，積分のどこが，原点に依存し，どこが依存しないかを把握しておく．
   原点に依存しない部分は，事前に計算しておくことで，計算量を削減できる．
   また，積分で面をトラバースする際に，係数行列のどの行に（どの点に）重みをかけるかを把握しておく．
   積分では，原点（行）・面・点（列）が重要である
   `map_Point_BEM_IGIGn_info_init`とは，
*/

void networkFace::setIntegrationInfo() {
  // set linear integration info
  this->map_Point_BEM_IGIGn_info_init.clear();
  this->map_Point_LinearIntegrationInfo_vector.clear();
  this->map_Point_LinearIntegrationInfo_vector.resize(3);
  this->map_Point_PseudoQuadraticIntegrationInfo_vector.clear();
  this->map_Point_PseudoQuadraticIntegrationInfo_vector.resize(3);

  // 最適化: 既存のdodecaPointsを再利用（なければ作成）
  if (!std::get<0>(this->dodecaPoints))
    this->setDodecaPoints();

  //@ -------------------------------------------------------------------------- */
  auto addIntegrationInfo = [&](networkPoint *const p, const DodecaPoints &dodecapoint) {
    //% -------------------------------------------------------------------------- */
    std::vector<BEM_IGIGn_info_type> temp;
    temp.reserve(24); // 6点 × 4 PseudoQuadPatch
    for (const auto &[pt, f] : dodecapoint.quadpoint.points_faces)
      temp.push_back({pt, f, 0., 0.});
    for (const auto &[pt, f] : dodecapoint.quadpoint_l0.points_faces)
      temp.push_back({pt, f, 0., 0.});
    for (const auto &[pt, f] : dodecapoint.quadpoint_l1.points_faces)
      temp.push_back({pt, f, 0., 0.});
    for (const auto &[pt, f] : dodecapoint.quadpoint_l2.points_faces)
      temp.push_back({pt, f, 0., 0.});
    this->map_Point_BEM_IGIGn_info_init[p] = std::move(temp);
    //% -------------------------------------------------------------------------- */
    auto X012 = ToX(this->getPoints(p));
    // 線形要素用のcrossは定数なので事前計算
    const auto cross_linear = Cross(X012[1] - X012[0], X012[2] - X012[0]);
    const auto norm_cross_linear = Norm(cross_linear);

    auto add = [&](const int i, const auto &GWGW) {
      std::vector<linear_triangle_integration_info> info_linears;
      std::vector<pseudo_quadratic_triangle_integration_info> info_quadratics;
      info_linears.reserve(GWGW.size());
      info_quadratics.reserve(GWGW.size());

      for (const auto &[t0, t1, ww] : GWGW) {
        auto N012_geometry = ModTriShape<3>(t0, t1);
        auto X = Dot(N012_geometry, X012);
        const double weight = ww * (1. - t0);
        info_linears.emplace_back(linear_triangle_integration_info{Tdd{t0, t1}, weight, N012_geometry, X, cross_linear, norm_cross_linear});
        //$ ------------------------------------ */
        auto [xi0, xi1, xi2] = N012_geometry;
        // 最適化: X_N6_cross()で一度に計算
        auto [X_quad, Nc_N0_N1_N2, cross_quad] = dodecapoint.X_N6_cross(xi0, xi1);
        info_quadratics.emplace_back(pseudo_quadratic_triangle_integration_info{Tdd{t0, t1}, weight, Nc_N0_N1_N2, X_quad, cross_quad, Norm(cross_quad)});
      }
      this->map_Point_LinearIntegrationInfo_vector[i][p] = std::move(info_linears);
      this->map_Point_PseudoQuadraticIntegrationInfo_vector[i][p] = std::move(info_quadratics);
    };
    add(0, __array_GW1xGW1__);
    add(1, __array_GW5xGW5__);
    add(2, __array_GW10xGW10__);
  };
  //@ -------------------------------------------------------------------------- */

  auto [p0, p1, p2] = this->getPoints();
  addIntegrationInfo(p0, *std::get<0>(this->dodecaPoints));
  addIntegrationInfo(p1, *std::get<1>(this->dodecaPoints));
  addIntegrationInfo(p2, *std::get<2>(this->dodecaPoints));
};

Tddd networkFace::normalVelocityRigidBody(const Tddd &X) const { return this->normal * Dot(this->normal, this->network->velocityRigidBody(X)); };

// コンストラクタ
networkFace::networkFace(Network *network_IN, networkPoint *p0, networkLine *l0, networkPoint *p1, networkLine *l1, networkPoint *p2, networkLine *l2) : Triangle(p0->X, p1->X, p2->X), network(network_IN) {
  try {
    this->network->add(this);
    this->setPoints(p0, l0, p1, l1, p2, l2);
    l0->add(this);
    l1->add(this);
    l2->add(this);
    p0->add(this);
    p1->add(this);
    p2->add(this);
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};

// コンストラクタ
networkFace::networkFace(Network *network_IN, networkPoint *p0, networkPoint *p1, networkPoint *p2) : Triangle(p0->X, p1->X, p2->X), network(network_IN) {
  try {
    this->network->add(this);
    auto l0 = link(p0, p1, network_IN);
    auto l1 = link(p1, p2, network_IN);
    auto l2 = link(p2, p0, network_IN);
    this->setPoints(p0, l0, p1, l1, p2, l2);
    l0->add(this);
    l1->add(this);
    l2->add(this);
    p0->add(this);
    p1->add(this);
    p2->add(this);
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};
// コンストラクタ
networkFace::networkFace(const netFp f) : Triangle(extractXtuple(f)), network(f->network), Lines(f->Lines), Points(f->Points), PLPLPL(f->PLPLPL) {
  this->network->add(this);
  this->penetratedBody = f->penetratedBody;
  this->normal = f->normal;
  this->angles = f->angles;
  this->area = f->area;
};
// b% -------------------------------------------------------------------------- */
// b% particlizeDetailsは普通のparticlizeに詳しい情報を加えて返す．
// b% 深さ毎に，面の頂点をシフトしてm線形補間に利用する．2021/11/17
// b% -------------------------------------------------------------------------- */

std::unordered_set<networkPoint *> networkFace::particlize(const double dx, const V_d &depth_list) {
  // depth_list: 法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など

  std::unordered_set<networkPoint *> ret;
  // double alpha;
  // T3Tddd X0X1X2;
  int count = 0;
  networkPoint *p0, *p1, *p2;
  for (const auto &d : depth_list /*double実際の長さ*/) {
    auto [p0_, p1_, p2_] = this->getPoints();
    if (count % 3 == 1) {
      p0 = p2_;
      p1 = p0_;
      p2 = p1_;
    } else if (count % 3 == 2) {
      p0 = p1_;
      p1 = p2_;
      p2 = p0_;
    } else {
      p0 = p0_;
      p1 = p1_;
      p2 = p2_;
    }
    T3Tddd X0X1X2 = {p0->getXtuple() + d / Dot(p0->getNormalTuple(), this->normal) * p0->getNormalTuple(), p1->getXtuple() + d / Dot(p1->getNormalTuple(), this->normal) * p1->getNormalTuple(), p2->getXtuple() + d / Dot(p2->getNormalTuple(), this->normal) * p2->getNormalTuple()};
    for (const auto &[xyz, t0t1] : triangleIntoPoints(X0X1X2, dx)) {
      auto p = new networkPoint(this->getNetwork(), xyz);
      p->particlize_info = {this, {p0, p1, p2}, t0t1, d, dx};
      ret.emplace(p);
      this->addParametricPoints(p);
    }
    count++;
  }
  return ret;
};
