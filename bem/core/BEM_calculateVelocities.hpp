#ifndef BEM_calculateVelocities_H
#define BEM_calculateVelocities_H

#include "BEM_BoundaryValues.hpp"
#include "BEM_midpoint_debug.hpp"
#include "Network.hpp"
#include "minMaxOfFunctions.hpp"
#include <filesystem>

//$ -------------------------------------------------------------------------- */
//$                         calculateVecToSurface                              */
//$ -------------------------------------------------------------------------- */

inline void writeNodeRelocationQuadraticDebugVTU(
    const std::filesystem::path& filename,
    const V_netFp& faces,
    const std::function<Tddd(const networkPoint*)>& point_position_of,
    const std::unordered_map<networkPoint*, Tddd>& point_correction,
    const std::unordered_map<networkPoint*, Tddd>& point_nearest,
    const std::unordered_map<networkLine*, Tddd>& line_correction,
    const std::unordered_map<networkLine*, Tddd>& line_nearest,
    const double stage_id) {
  FILE* fp = fopen(filename.string().c_str(), "wb");
  if (!fp)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename.string() + " can not be opened");

  struct DebugPoint {
    Tddd position;
    Tddd correction;
    Tddd nearest;
  };

  std::vector<DebugPoint> all_points;
  all_points.reserve(faces.size() * 6);

  for (const auto* f : faces) {
    auto [p0, l0, p1, l1, p2, l2] = f->PLPLPL;

    auto push_vertex = [&](networkPoint* p) {
      const auto pos = point_position_of(p);
      const auto it_corr = point_correction.find(p);
      const auto it_near = point_nearest.find(p);
      all_points.push_back({pos, (it_corr != point_correction.end()) ? it_corr->second : Tddd{0., 0., 0.}, (it_near != point_nearest.end()) ? it_near->second : Tddd{0., 0., 0.}});
    };
    auto push_midpoint = [&](networkLine* l) {
      auto [pa, pb] = l->getPoints();
      const auto pos = 0.5 * (point_position_of(pa) + point_position_of(pb));
      const auto it_corr = line_correction.find(l);
      const auto it_near = line_nearest.find(l);
      all_points.push_back({pos, (it_corr != line_correction.end()) ? it_corr->second : Tddd{0., 0., 0.}, (it_near != line_nearest.end()) ? it_near->second : Tddd{0., 0., 0.}});
    };

    push_vertex(p0);
    push_vertex(p1);
    push_vertex(p2);
    push_midpoint(l0);
    push_midpoint(l1);
    push_midpoint(l2);
  }

  fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", (int)faces.size(), (int)all_points.size());

  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for (const auto& p : all_points)
    fprintf(fp, "%s %s %s ", NumtoString(std::get<0>(p.position)).c_str(), NumtoString(std::get<1>(p.position)).c_str(), NumtoString(std::get<2>(p.position)).c_str());
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<PointData>\n");

  fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='stage' format='ascii'>\n");
  for (std::size_t i = 0; i < all_points.size(); ++i)
    fprintf(fp, "%s ", NumtoString(stage_id).c_str());
  fprintf(fp, "\n</DataArray>\n");

  auto write_vec = [&](const char* name, const auto getter) {
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='%s' format='ascii'>\n", name);
    for (const auto& p : all_points) {
      const auto v = getter(p);
      fprintf(fp, "%s %s %s ", NumtoString(std::get<0>(v)).c_str(), NumtoString(std::get<1>(v)).c_str(), NumtoString(std::get<2>(v)).c_str());
    }
    fprintf(fp, "\n</DataArray>\n");
  };
  write_vec("correction", [](const auto& p) { return p.correction; });
  write_vec("vectorToNextSurface", [](const auto& p) { return p.nearest; });

  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for (int i = 0; i < (int)faces.size(); ++i)
    fprintf(fp, "%d %d %d %d %d %d ", 6 * i + 0, 6 * i + 1, 6 * i + 2, 6 * i + 3, 6 * i + 4, 6 * i + 5);
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for (int i = 0; i < (int)faces.size(); ++i)
    fprintf(fp, "%d ", 6 * (i + 1));
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for (std::size_t i = 0; i < faces.size(); ++i)
    fprintf(fp, "22 ");
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

inline void logNodeRelocationMidpointStats(
    const int steps,
    const std::unordered_map<networkLine*, Tddd>& midpoint_line_nearest_raw,
    const std::unordered_map<networkLine*, Tddd>& midpoint_line_correction,
    const std::vector<std::tuple<networkLine*, networkPoint*, Tddd, networkPoint*, Tddd>>& linear_midpoint_endpoint_corrections,
    const std::vector<Tddd>& midpoint_stage_delta) {
  std::size_t line_nearest_nonzero = 0, line_corr_nonzero = 0, penaltyA_nonzero = 0, penaltyB_nonzero = 0, point_delta_nonzero = 0;
  double line_nearest_sum = 0., line_nearest_max = 0.;
  double line_corr_sum = 0., line_corr_max = 0.;
  double penaltyA_sum = 0., penaltyA_max = 0.;
  double penaltyB_sum = 0., penaltyB_max = 0.;
  double point_delta_sum = 0., point_delta_max = 0.;

  for (const auto& [l, nearest] : midpoint_line_nearest_raw) {
    const double n = Norm(nearest);
    if (n > 1e-12) {
      ++line_nearest_nonzero;
      line_nearest_sum += n;
      line_nearest_max = std::max(line_nearest_max, n);
    }
  }
  for (const auto& [l, corr] : midpoint_line_correction) {
    const double n = Norm(corr);
    if (n > 1e-12) {
      ++line_corr_nonzero;
      line_corr_sum += n;
      line_corr_max = std::max(line_corr_max, n);
    }
  }
  for (const auto& [l, pA, linear_midpoint_penaltyA, pB, linear_midpoint_penaltyB] : linear_midpoint_endpoint_corrections) {
    const double nA = Norm(linear_midpoint_penaltyA);
    const double nB = Norm(linear_midpoint_penaltyB);
    if (nA > 1e-12) {
      ++penaltyA_nonzero;
      penaltyA_sum += nA;
      penaltyA_max = std::max(penaltyA_max, nA);
    }
    if (nB > 1e-12) {
      ++penaltyB_nonzero;
      penaltyB_sum += nB;
      penaltyB_max = std::max(penaltyB_max, nB);
    }
  }
  for (const auto& delta : midpoint_stage_delta) {
    const double n = Norm(delta);
    if (n > 1e-12) {
      ++point_delta_nonzero;
      point_delta_sum += n;
      point_delta_max = std::max(point_delta_max, n);
    }
  }

  auto mean_or_zero = [](const double sum, const std::size_t count) {
    return count ? sum / static_cast<double>(count) : 0.0;
  };

  std::cout << "[midpoint debug] iter=" << steps
            << " lines=" << linear_midpoint_endpoint_corrections.size()
            << " nearest_nonzero=" << line_nearest_nonzero
            << " nearest_mean=" << mean_or_zero(line_nearest_sum, line_nearest_nonzero)
            << " nearest_max=" << line_nearest_max
            << " corr_nonzero=" << line_corr_nonzero
            << " corr_mean=" << mean_or_zero(line_corr_sum, line_corr_nonzero)
            << " corr_max=" << line_corr_max
            << " penaltyA_nonzero=" << penaltyA_nonzero
            << " penaltyA_mean=" << mean_or_zero(penaltyA_sum, penaltyA_nonzero)
            << " penaltyA_max=" << penaltyA_max
            << " penaltyB_nonzero=" << penaltyB_nonzero
            << " penaltyB_mean=" << mean_or_zero(penaltyB_sum, penaltyB_nonzero)
            << " penaltyB_max=" << penaltyB_max
            << " point_nonzero=" << point_delta_nonzero
            << " point_mean=" << mean_or_zero(point_delta_sum, point_delta_nonzero)
            << " point_max=" << point_delta_max
            << std::endl;
}

inline void setNodeRelocationStageBuffer(
    const std::vector<networkPoint*>& points,
    const std::vector<Tddd>& point_stage_delta,
    const std::unordered_map<networkLine*, Tddd>&) {
  for (std::size_t i = 0; i < points.size(); ++i)
    points[i]->vecToSurface_BUFFER = (i < point_stage_delta.size()) ? point_stage_delta[i] : Tddd{0., 0., 0.};
}

inline Tddd RK_without_Ubuff(const networkPoint* p) {
  return p->RK_X.getX(p->u_reloc);
};
inline T3Tddd RK_without_Ubuff(const networkPoint* p0, const networkPoint* p1, const networkPoint* p2) { return {RK_without_Ubuff(p0), RK_without_Ubuff(p1), RK_without_Ubuff(p2)}; };
inline T3Tddd RK_without_Ubuff(const T_PPP& p012) { return RK_without_Ubuff(std::get<0>(p012), std::get<1>(p012), std::get<2>(p012)); };
inline T3Tddd RK_without_Ubuff(const networkFace* f) { return RK_without_Ubuff(f->getPoints()); };
inline T2Tddd RK_without_Ubuff(const networkPoint* p0, const networkPoint* p1) { return {RK_without_Ubuff(p0), RK_without_Ubuff(p1)}; };
inline Tddd RK_without_Ubuff_Normal(const networkFace* f) { return TriangleNormal(RK_without_Ubuff(f->getPoints())); };
inline Tddd RK_without_Ubuff_Normal(const networkPoint* p) {
  Tddd normal = {0., 0., 0.};
  double a = 0, total = 0;
  for (const auto& f : p->getBoundaryFaces()) {
    a = TriangleArea(RK_without_Ubuff(f));
    FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
    total += a;
  }
  if (!(total > 0) || !std::isfinite(total) || !isFinite(normal))
    return {0., 0., 0.};
  return Normalize(normal / total);
};
inline Tddd getNextNormalDirichlet_BEM(const networkPoint* p) {
  Tddd normal = {0., 0., 0.};
  double a = 0, total = 0;
  for (const auto& f : p->getFacesDirichlet()) {
    a = TriangleArea(RK_without_Ubuff(f));
    // normal += a * RK_without_Ubuff_Normal(f);
    FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
    total += a;
  }
  if (!(total > 0) || !std::isfinite(total) || !isFinite(normal))
    return {0., 0., 0.};
  return Normalize(normal / total);
};
inline Tddd getNextNormalNeumann_BEM(const networkPoint* p) {
  Tddd normal = {0., 0., 0.};
  double a = 0, total = 0;
  for (const auto& f : p->getFacesNeumann()) {
    a = TriangleArea(RK_without_Ubuff(f));
    // normal += a * RK_without_Ubuff_Normal(f);
    FusedMultiplyIncrement(a, RK_without_Ubuff_Normal(f), normal);
    total += a;
  }
  if (!(total > 0) || !std::isfinite(total) || !isFinite(normal))
    return {0., 0., 0.};
  return Normalize(normal / total);
};

// Tddd RK_with_Ubuff(const networkPoint *p) { return p->RK_X.getX(p->u_reloc + p->vecToSurface / p->RK_X.getdt()); };
inline Tddd RK_with_Ubuff(const networkPoint* p) {
  // return p->RK_X.toReachAtNextTimeQ(p->RK_X.getX(p->u_reloc) + p->vecToSurface);
  const auto U = p->RK_X.getVectorToReachAtNextTimeQ(RK_without_Ubuff(p) + p->vecToSurface);
  return p->RK_X.getX(U);
};
inline Tddd RK_with_Ubuff(const networkPoint* p, const Tddd& vecToSurface) {
  // return p->RK_X.getX(p->u_reloc + vecToSurface / p->RK_X.getdt());
  const auto U = p->RK_X.getVectorToReachAtNextTimeQ(RK_without_Ubuff(p) + vecToSurface);
  return p->RK_X.getX(U);
};
inline T3Tddd RK_with_Ubuff(const networkPoint* p0, const networkPoint* p1, const networkPoint* p2) { return {RK_with_Ubuff(p0), RK_with_Ubuff(p1), RK_with_Ubuff(p2)}; };
inline T3Tddd RK_with_Ubuff(const T_PPP& p012) { return {RK_with_Ubuff(p012[0]), RK_with_Ubuff(p012[1]), RK_with_Ubuff(p012[2])}; };
inline T3Tddd RK_with_Ubuff(const networkFace* f) { return RK_with_Ubuff(f->getPoints()); };
inline Tddd RK_with_Ubuff_Normal(const networkFace* f) { return TriangleNormal(RK_with_Ubuff(f->getPoints())); };
inline double RK_with_Ubuff_Area(const networkFace* f) { return TriangleArea(RK_with_Ubuff(f->getPoints())); };
inline Tddd RK_with_Ubuff_Normal(const networkPoint* p) {
  Tddd normal = {0., 0., 0.};
  for (const auto& f : p->getBoundaryFaces())
    normal += RK_with_Ubuff_Area(f) * RK_with_Ubuff_Normal(f);
  return Normalize(normal);
};

inline void add_vecToSurface_BUFFER_to_vecToSurface(const auto& p, const double a = 1.) {
  if (isFinite(p->vecToSurface_BUFFER))
    p->vecToSurface += p->vecToSurface_BUFFER * a;
};

/* ---------- midpoint RK helpers ---------- */
inline Tddd RK_without_Ubuff(const networkLine* l) {
  if (l->RK_X.steps == 0)
    return l->X_mid;
  return l->RK_X.getX(l->u_reloc);
}
inline Tddd RK_with_Ubuff(const networkLine* l) {
  if (l->RK_X.steps == 0)
    return l->X_mid + l->vecToSurface;
  const auto U = l->RK_X.getVectorToReachAtNextTimeQ(RK_without_Ubuff(l) + l->vecToSurface);
  return l->RK_X.getX(U);
}

// mooringで利用
inline std::array<double, 3> nextPositionOnBody(Network* net, networkPoint* p) {
  if (net->isSoftBody || net->inputJSON.find("relative_velocity")) {
    auto v = p->velocityTranslational();

    // 軸ごとの固定を適用（SoftBodyでもRigidBodyと同じ期待挙動にする）
    for (auto i = 0; i < 3; ++i)
      if (net->isFixed[i] == true)
        v[i] = 0;

    // 3並進が全部固定ならゼロ
    if (net->isFixed.size() == 1)
      if (net->isFixed[0] == true)
        v.fill(0.0);

    return p->RK_X.getX(v);
  } else if (net->isRigidBody) {
    auto velocity = net->velocityTranslational();
    auto rotation = net->velocityRotational();
    for (auto i = 0; i < 3; ++i)
      if (net->isFixed[i])
        velocity[i] = 0;
    for (auto i = 3; i < 6; ++i)
      if (net->isFixed[i])
        rotation[i - 3] = 0;
    auto COM_next = net->RK_COM.getX(velocity);
    auto Q_next = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->RK_Q.getX())));
    return rigidTransformation(net->ICOM, COM_next, Q_next.Rv(), p->initialX);
  }
  return getPosition(p);
}

/* -------------------------------------------------------------------------- */

inline std::vector<std::tuple<T3Tddd, T3Tddd>> nextBodyVertices(const auto& Fs) {
  std::vector<std::tuple<T3Tddd, T3Tddd>> ret(Fs.size());
  int i = -1;
  for (auto& f : Fs) {
    i++;
    auto net = f->getNetwork();
    auto [p0, p1, p2] = f->getPoints();

    if (net->isSoftBody || net->inputJSON.find("relative_velocity")) {
      auto v0 = p0->velocityTranslational();
      auto v1 = p1->velocityTranslational();
      auto v2 = p2->velocityTranslational();
      if (std::ranges::all_of(net->isFixed, [](bool isfixed) { return isfixed; })) {
        v0.fill(0.);
        v1.fill(0.);
        v2.fill(0.);
      }
      auto X0 = p0->RK_X.getX(v0);
      auto X1 = p1->RK_X.getX(v1);
      auto X2 = p2->RK_X.getX(v2);
      ret[i] = {{p0->X, p1->X, p2->X}, {X0, X1, X2}};
    } else if (net->isRigidBody) {
      auto velocity = net->velocityTranslational();
      auto rotation = net->velocityRotational();

      if (std::ranges::all_of(velocity, [](double v) { return v == 0; }) && std::ranges::all_of(rotation, [](double v) { return v == 0; })) {
        ret[i] = {{p0->X, p1->X, p2->X}, {p0->X, p1->X, p2->X}};
        continue;
      }

      for (auto k = 0; k < 3; ++k) {
        if (net->isFixed[k])
          velocity[k] = 0;
        if (net->isFixed[k + 3])
          rotation[k] = 0;
      }

      auto next_COM = net->RK_COM.getX(velocity);
      auto next_Q = Quaternion(net->RK_Q.getX(AngularVelocityTodQdt(rotation, net->RK_Q.getX())));
      auto X_next = [&](const auto& p) { return rigidTransformation(net->ICOM, next_COM, next_Q.Rv(), p->initialX); };

      ret[i] = {{p0->X, p1->X, p2->X}, {X_next(p0), X_next(p1), X_next(p2)}};
    } else {
      ret[i] = {{p0->X, p1->X, p2->X}, {p0->X, p1->X, p2->X}};
    }
  }
  return ret;
};

/* -------------------------------------------------------------------------- */

// やはりこの辺りなのだろう

// \label{BEM:vectorToNextSurface}
// 参照三角形上の均等格子点（TriShape<6>/TriShape<3> の入力用、コンパイル時計算）
inline constexpr auto t0t1 = UniformPointsOnTriangle_00_10_01<6>();       // 28点
inline constexpr auto t0t1_high = UniformPointsOnTriangle_00_10_01<15>(); // 136点

/* -------------------------------------------------------------------------- */
/*     Dirichlet面上の最近点探索                                               */
/* -------------------------------------------------------------------------- */

// グリッド探索の初期値から、曲面上の真の最近点パラメタをNewton法で精緻化
// 勾配 g = {(X-P)·∂X/∂t0, (X-P)·∂X/∂t1} = 0 を解く
// factor 2 は g と H 両方にかかるため省略
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

// pseudo-quadratic補間によるDirichlet面上の最近点
// dodecaPoints[0] の補間で曲面を構成し、均一格子点でサンプルして最近点を探す
// 均一分布なので頂点インデックスは0固定で面全体をカバーできる
// out_param: 非nullなら最近点の(t0, t1)パラメタを書き出す
inline Tddd NearestOnDirichletFace_pseudo(const Tddd& X_shifted, const networkFace* f, const auto& t0t1_param, Tdd* out_param, auto getPos) {
  auto points = f->dodecaPoints[0]->interpolate(t0t1_param, [&](networkPoint* q) -> Tddd { return getPos(q); });
  double nearest_dist = 1E+20, dist = 0.;
  Tddd X_nearest = {1E+20, 1E+20, 1E+20};
  std::size_t best_idx = 0;
  for (std::size_t idx = 0; idx < points.size(); ++idx) {
    if ((dist = Norm(points[idx] - X_shifted)) < nearest_dist) {
      nearest_dist = dist;
      X_nearest = points[idx];
      best_idx = idx;
    }
  }
  Tdd best_param = t0t1_param[best_idx];
  // Newton精緻化
  auto& dode = f->dodecaPoints[0];
  best_param = refineNearestParam(X_shifted, best_param, [&](double t0, double t1) -> Tddd { return dode->interpolate(t0, t1, getPos); }, [&](double t0, double t1, int i, int j) -> Tddd {
        if (i == 1 && j == 0) return dode->template D_interpolate<1, 0>(t0, t1, getPos);
        if (i == 0 && j == 1) return dode->template D_interpolate<0, 1>(t0, t1, getPos);
        if (i == 2 && j == 0) return dode->template D_interpolate<2, 0>(t0, t1, getPos);
        if (i == 0 && j == 2) return dode->template D_interpolate<0, 2>(t0, t1, getPos);
        return dode->template D_interpolate<1, 1>(t0, t1, getPos); });
  X_nearest = dode->interpolate(best_param[0], best_param[1], getPos);
  if (out_param)
    *out_param = best_param;
  return X_nearest;
}

// true-quadratic補間によるDirichlet面上の最近点
// TriShape<6>で6節点（頂点3+midpoint3）の二次補間した曲面上でX_shiftedに最も近い点を返す
// out_param: 非nullなら最近点の(t0, t1)パラメタを書き出す
inline Tddd NearestOnDirichletFace_true_quad(const Tddd& X_shifted, const networkFace* f, const auto& t0t1_param, Tdd* out_param, auto getPos) {
  auto [p0, l0, p1, l1, p2, l2] = f->PLPLPL;
  T6Tddd X123456 = {getPos(p0), getPos(p1), getPos(p2), getPos(l0), getPos(l1), getPos(l2)};
  double nearest_dist = 1E+20, dist = 0.;
  Tddd X_nearest = {1E+20, 1E+20, 1E+20}, X_interp;
  Tdd best_param = {0., 0.};
  for (const Tdd& param : t0t1_param) {
    X_interp = Dot(TriShape<6>(param), X123456);
    if ((dist = Norm(X_interp - X_shifted)) < nearest_dist) {
      nearest_dist = dist;
      X_nearest = X_interp;
      best_param = param;
    }
  }
  // Newton精緻化
  best_param = refineNearestParam(X_shifted, best_param, [&](double t0, double t1) -> Tddd { return Dot(TriShape<6>(t0, t1), X123456); }, [&](double t0, double t1, int i, int j) -> Tddd {
        if (i == 1 && j == 0) return Dot(D_TriShape<6, 1, 0>(t0, t1), X123456);
        if (i == 0 && j == 1) return Dot(D_TriShape<6, 0, 1>(t0, t1), X123456);
        if (i == 2 && j == 0) return Dot(D_TriShape<6, 2, 0>(t0, t1), X123456);
        if (i == 0 && j == 2) return Dot(D_TriShape<6, 0, 2>(t0, t1), X123456);
        return Dot(D_TriShape<6, 1, 1>(t0, t1), X123456); });
  X_nearest = Dot(TriShape<6>(best_param), X123456);
  if (out_param)
    *out_param = best_param;
  return X_nearest;
}

// Dirichlet面上の最近点探索（モード自動選択）
// out_param: 非nullなら最近点の(t0, t1)パラメタを書き出す
inline Tddd NearestOnDirichletFace(const Tddd& X_shifted, const networkFace* f, Tdd* out_param, auto getPos) {
  switch (node_relocation_surface) {
  case NodeRelocationSurface::linear: {
    auto [q0, q1, q2] = f->getPoints();
    auto [wa, wb, X_near, normal] = Nearest_(X_shifted, T3Tddd{getPos(q0), getPos(q1), getPos(q2)});
    if (out_param)
      *out_param = {wa, wb};
    return X_near;
  }
  case NodeRelocationSurface::true_quadratic:
    return NearestOnDirichletFace_true_quad(X_shifted, f, t0t1, out_param, getPos);
  default:
    return NearestOnDirichletFace_pseudo(X_shifted, f, t0t1, out_param, getPos);
  }
}

// Default overload: uses RK_without_Ubuff for node positions (for use during RK4)
inline Tddd NearestOnDirichletFace(const Tddd& X_shifted, const networkFace* f, Tdd* out_param = nullptr) {
  return NearestOnDirichletFace(X_shifted, f, out_param, [](const auto* q) -> Tddd { return RK_without_Ubuff(q); });
}

/* -------------------------------------------------------------------------- */

// networkPoint* / networkLine* 共通のスナッピング関数
// X_shifted: 補正済み予測位置（点: RK_with_Ubuff(p), 辺: 頂点平均）
// entity: networkPoint* or networkLine* （Dirichlet/Neumann/CORNER/contact_range等の共通インターフェースを持つ）
template <typename Entity>
inline Tddd vectorToNextSurface(Entity entity, Tddd X_shifted) {

  try {
    const bool has_dirichlet_state = hasAnyDirichletBoundaryState(entity);
    const bool has_neumann_state = hasAnyNeumannBoundaryState(entity);
    const bool is_mixed_state = has_dirichlet_state && has_neumann_state;

    if (!has_dirichlet_state && !has_neumann_state)
      return {0., 0., 0.};

    Tddd X_not_shifted = RK_without_Ubuff(entity);
    Tddd X, vecToDirichlet = {1E+20, 1E+20, 1E+20};

    struct DirectionInfo {
      double distance;
      Tddd direction;
    };
    std::vector<DirectionInfo> direction_infos;

    auto faces = entity->getBoundaryFaces();
    if (faces.empty())
      return {0., 0., 0.};

    networkFace* closest_face = nullptr;
    Tdd best_face_param = {0., 0.};

    //! CORNERは２段階処理する必要がある．まずDirichlet面に張り付く．その後，Neumann面に張り付く．
    //! Dirichlet境界面への接近ベクトルを計算
    if (has_dirichlet_state) {
      //! Neumann面には後で境界条件が付与されるので問題ない
      for (const auto& f : faces) {
        if (getNodeFaceBoundaryType(entity, f) == NodeFaceBoundaryType::Dirichlet) {
          Tdd face_param;
          X = NearestOnDirichletFace(X_shifted, f, &face_param);

          if (Norm(vecToDirichlet) >= Norm(X - X_shifted)) {
            vecToDirichlet = X - X_shifted;
            closest_face = f;
            best_face_param = face_param;
          }
        }
      }
      if (closest_face != nullptr) {
        entity->relocation_face = closest_face;
        entity->relocation_param = best_face_param;
        if (!has_neumann_state)
          return vecToDirichlet;
        X_shifted += vecToDirichlet; // for CORNER, shift the position to the closest Dirichlet surface
      } else {
        if (!has_neumann_state)
          return {0., 0., 0.};
        vecToDirichlet = {0., 0., 0.};
      }
    }

    auto add_vector = [&](const Tddd& V, Tddd n) {
      if (!isFinite(V) || !isFinite(n) || !(Norm(n) > 0))
        return;
      n = Normalize(n);
      double new_dist = Dot(V, n);
      // delete if the new vector is closer in the same direction
      std::erase_if(direction_infos, [&](const auto& info) { return Dot(n, info.direction) > std::cos(M_PI / 180. * 20.) && std::abs(info.distance) > std::abs(Dot(V, n)); });
      // Add the new vector if no existing better vector is in the same direction
      if (std::ranges::none_of(direction_infos, [&](const auto& info) { return Dot(n, info.direction) > std::cos(M_PI / 180. * 20.); })) {
        direction_infos.push_back({new_dist, n});
      }
    };

    // ! Neumann境界面への接近ベクトルを計算
    std::vector<T3Tddd> next_triangles;
    int isInContact_pass_count = 0;
    if (has_neumann_state) {

      const double short_range = 0.01; //@ short_range * radius is the detection range of the structure face

      bool has_contact = std::ranges::any_of(entity->dofs, [](const auto& kv) { return !kv.second.contact_opponent_faces.empty(); });
      if (!has_contact && entity->penetratedBody != nullptr) {
        auto [f, X_nearest] = entity->penetratedBody->Nearest(X_shifted, [&](const networkPoint* p) { return RK_without_Ubuff(p); });
        if (f != nullptr) {
          add_vector(X_nearest - X_shifted, Normalize(X_nearest - X_shifted));
        }
      }

      auto Vec_CurrentX012_NextX012 = nextBodyVertices(bfs(getEffectiveContactFaces(entity), 2));
      entity->debug_body_vertices_count = static_cast<int>(Vec_CurrentX012_NextX012.size());

      if (!Vec_CurrentX012_NextX012.empty()) {
        for (const auto& f : faces)
          if (getNodeFaceBoundaryType(entity, f) == NodeFaceBoundaryType::Neumann) {
            for (const auto& [struct_vertex, next_struct_vertex] : Vec_CurrentX012_NextX012) {
              if (!isInContact(X_not_shifted, f->normal, struct_vertex, entity->contact_range))
                continue;
              ++isInContact_pass_count;
              next_triangles.emplace_back(next_struct_vertex);
              auto [t0, t1, X_nearest, n] = Nearest_(X_shifted, next_struct_vertex);
              Tddd To_nearest = X_nearest - X_shifted;
              double distance = Norm(To_nearest);
              if ((isFlat(n, To_nearest, contactAcceptanceAngle(distance, entity->contact_range)) || isFlat(n, -To_nearest, contactAcceptanceAngle(distance, entity->contact_range))) && (entity->contact_range >= distance))
                add_vector(To_nearest, n);
              else if ((short_range * entity->contact_range >= distance))
                add_vector(To_nearest, n);
            }
          }
      }
    }
    std::vector<double> distances;
    std::vector<Tddd> directions;
    for (const auto& info : direction_infos) {
      distances.push_back(info.distance);
      directions.push_back(info.direction);
    }

    entity->debug_direction_info_count = static_cast<int>(direction_infos.size());
    entity->debug_contact_faces_count = static_cast<int>(std::ranges::count_if(entity->dofs, [](const auto& kv) { return !kv.second.contact_opponent_faces.empty(); }));
    entity->debug_isInContact_pass_count = isInContact_pass_count;

    if (distances.empty())
      return (is_mixed_state ? vecToDirichlet : Tddd{0., 0., 0.});

    auto V_opt = optimalVector(distances, directions, {0., 0., 0.});
    if (!next_triangles.empty()) {
      auto X_target = X_shifted + V_opt;
      auto X_projected = Nearest(X_target, next_triangles);
      V_opt = X_projected - X_shifted;
    }
    // relocation_face / relocation_param は calculateVecToSurface 後に
    // setRelocationParam で最終位置から再計算されるため、ここでは更新しない。
    return V_opt + (is_mixed_state ? vecToDirichlet : Tddd{0., 0., 0.});
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  }
}

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_3_1_INITIAL_VALUE_PROBLEM

### Arbitrary Lagrangian–Eulerian Methods (ALE)

水面を移動する物体が存在する場合，格子が潰れてしまうため，何らかの方法で格子を綺麗に保たなければ，長時間の計算は不可能である．
その方法の一つが，節点の位置を流速$\nabla \phi$ではなく任意のベクトル${\bf v}$で移動させる，ALEである．
${\bf v}$で移動する節点位置を${\mathbfit\chi}(t)$と置くと，
$\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\nabla\phi\cdot\nabla\phi$の代わりに，
$\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\frac{d{\mathbfit\chi}}{dt} \cdot\nabla\phi,\frac{d{\mathbfit\chi}}{dt} = \bf v$
を使って$\phi$を時間発展させると次時刻の節点位置${\mathbfit\chi}(t+\delta t)$での$\phi$が得られる．

ディリクレ節点（水面）：

求めた流速から，次の時刻の境界面$\Omega(t+\Delta t)$を見積もり，その面上で節点を移動させ歪さを解消する．
修正ベクトルは，$\Delta t$で割り，求めた流速$\nabla \phi$に足し合わせて，節点を時間発展させる．

ノイマン節点：

ノイマン節点も修正流速を加え時間発展させる．
ただし，ノイマン節点の修正流速に対しては，節点が水槽の角から離れないように，工夫を施している．

\ref{BEM:calculateVecToSurface}{`calculateVecToSurface`}で$\Omega(t+\Delta t)$上へのベクトルを計算する．

1. まず，\ref{BEM:vectorTangentialShift}{`vectorTangentialShift`}で接線方向にシフトし，
2. \ref{BEM:vectorToNextSurface}{`vectorToNextSurface`}で近くの$\Omega(t+\Delta t)$上へのベクトルを計算する．

#### 使っているALEの手法

`CircumradiusToInradius`をもとに，節点に隣接する面の歪みを評価し，その歪みに応じて節点を修正する．
面の歪みを修正するための節点の移動方向は，節点の向かい側にある辺の中点から辺が作る正三角形の高さの方向としている．
歪みを小さくするための方向を歪みの`CircumradiusToInradius`の微分から求めることもできるだろうが，プログラムが複雑になるためこの方法を採用している．

* 重みの最大値と最小値を設定している
* よりディリクレ面の歪みを緩和するように重みを大きくしている
* より喫水線の歪みを緩和するように重みを大きくしている

*/

// \label{BEM:calculateVecToSurface}
Tddd DistorsionMeasureWeightedSmoothingVector_modified(const networkPoint* p, const std::array<double, 3>& current_pX, std::function<Tddd(const networkPoint*)> position) {
  const int max_sum_depth = 20;
  auto faces = p->CORNER ? p->getFacesDirichlet() : p->getBoundaryFaces();
  thread_local std::vector<double> weights;
  weights.clear();
  thread_local std::vector<Tddd> positions;
  positions.clear();
  Tddd X0, X1, X2, X3, Xmid, vertical, X_ideal, To_ideal;
  double W, height, normalized_discrepancy;

  // tetraがある場合
  std::vector<networkTetra*> tetras = p->getTetras();

  auto isSurface = [](const networkTetra* tet) -> bool {
    for (const auto& f : tet->Faces)
      if (f->Tetras[0] == nullptr || f->Tetras[1] == nullptr)
        return true;
    return false;
  };

  if (!tetras.empty()) {
    for (const auto& tet : tetras) {
      auto [p0, p1, p2, p3] = tet->getPoints(p);
      X0 = current_pX;
      X1 = position(p1);
      X2 = position(p2);
      X3 = position(p3);
      W = CircumradiusToInradius(X0, X1, X2, X3) - 3;
      //! 　面それぞれに対する，修正方向は正三角形の高さの方向としている
      auto Xmid = (X1 + X2 + X3) / 3.;
      auto cross = Cross(X2 - X1, X3 - X1);
      const auto cross_norm = Norm(cross);
      if (!(cross_norm > 0) || !std::isfinite(cross_norm))
        continue;
      auto vertical = cross / cross_norm;
      if (Dot(vertical, X0 - Xmid) < 0)
        vertical = -vertical;
      auto area = cross_norm * 0.5;
      auto approx_edge_length = (Norm(X1 - X2) + Norm(X2 - X3) + Norm(X3 - X1)) / 3.;
      auto height = approx_edge_length * std::sqrt(6.) / 3.; // 正四面体の高さ
      if (!(height > 0) || !std::isfinite(height))
        continue;
      if (Dot(vertical, current_pX - Xmid) < 0)
        vertical = -vertical;
      auto To_ideal = height * vertical + (Xmid - current_pX);
      auto normalized_discrepancy = Norm(To_ideal) / height;
      double a = 1.;
      weights.push_back(a * (normalized_discrepancy + W));
      positions.emplace_back(To_ideal);
    }
  }
  //   else
  //     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Tetrahedra are required for DistorsionMeasureWeightedSmoothingVector_modified.");

  {
    for (const auto& f : faces) {
      auto [p0, p1, p2] = f->getPoints(p);
      X0 = current_pX;
      X1 = position(p1);
      X2 = position(p2);
      //! 重みの計算
      W = CircumradiusToInradius(X0, X1, X2) - 2.; // CircumradiusToInradius(X0, X1, X2)は最小で2
      //! 　面それぞれに対する，修正方向は正三角形の高さの方向としている
      Xmid = (X2 + X1) * 0.5;
      const auto e = X2 - X1;
      const auto e_norm = Norm(e);
      if (!(e_norm > 0) || !std::isfinite(e_norm))
        continue;
      auto v_raw = Chop(X0 - Xmid, e);
      const auto v_norm = Norm(v_raw);
      if (!(v_norm > 0) || !std::isfinite(v_norm))
        continue;
      vertical = v_raw / v_norm;
      height = e_norm * std::sqrt(3.) * 0.5;
      if (!(height > 0) || !std::isfinite(height))
        continue;
      X_ideal = height * vertical + Xmid;
      To_ideal = X_ideal - current_pX;
      normalized_discrepancy = Norm(To_ideal) / height;
      //! 重みと位置を保存
      double a = 5.;
      //   a *= std::pow(1.1, std::ranges::count_if(f->getPoints(), [](const auto &p) { return p->isMultipleNode; })); // 重みは，面に含まれるマルチプルノードの数に応じて増加する
      weights.push_back(a * (normalized_discrepancy + W));
      positions.emplace_back(To_ideal);
    }
  }

  if (weights.empty())
    return {0., 0., 0.};

  double sum = Sum(weights);
  if (sum < 1E-12)
    weights /= static_cast<double>(weights.size());
  else
    weights /= sum;
  Tddd V = {0., 0., 0.};
  for (int i = 0; i < weights.size(); ++i)
    V += weights[i] * positions[i];
  return V;
};

/* -------------------------------------------------------------------------- */

void calculateVecToSurface(const Network& net,
                           const int loop,
                           const double coefIN,
                           const std::filesystem::path* debug_output_directory = nullptr) {
  auto points = ToVector(net.getPoints());
  const std::size_t n_points = points.size();
  const bool enable_stage_debug = loop > 0 && debug_output_directory && !debug_output_directory->empty();

  //! 初期化
  for (const auto& p : points) {
    p->vecToSurface.fill(0.);
    p->vecToSurface_BUFFER.fill(0.);
    p->vecToSurface_BUFFER_BUFFER.fill(0.);
  }

  for (const auto& l : net.getBoundaryLines()) {
    l->vecToSurface.fill(0.);
  }

  // 最適化: RK_with_Ubuff結果をキャッシュ
  std::vector<Tddd> cached_positions(n_points);
  // Shift limit should be based on the *un-corrected* (without vecToSurface) geometry.
  // Using RK_with_Ubuff() here creates a positive feedback loop:
  // once a point drifts, local element sizes appear larger -> larger allowed shifts -> runaway.
  std::vector<Tddd> base_positions(n_points);

  _Pragma("omp parallel for") for (std::size_t i = 0; i < n_points; ++i) base_positions[i] = RK_without_Ubuff(points[i]);

  auto update_cache = [&]() { _Pragma("omp parallel for") for (std::size_t i = 0; i < n_points; ++i) cached_positions[i] = RK_with_Ubuff(points[i]); };

  // ポイントからインデックスへのマップ（隣接点の位置取得用）
  std::unordered_map<const networkPoint*, std::size_t> point_to_index;
  point_to_index.reserve(n_points);
  for (std::size_t i = 0; i < n_points; ++i)
    point_to_index[points[i]] = i;

  auto get_cached_position = [&](const networkPoint* p) -> Tddd {
    auto it = point_to_index.find(p);
    if (it != point_to_index.end())
      return cached_positions[it->second];
    return RK_with_Ubuff(p); // フォールバック（他のネットワークの点など）
  };

  auto get_base_position = [&](const networkPoint* p) -> Tddd {
    auto it = point_to_index.find(p);
    if (it != point_to_index.end())
      return base_positions[it->second];
    return RK_without_Ubuff(p); // フォールバック（他のネットワークの点など）
  };

  auto lines = net.getBoundaryLines();

  PVDWriter node_relocation_debug_pvd(((enable_stage_debug ? *debug_output_directory : std::filesystem::path{}) / (net.getName() + "_node_relocation_debug.pvd")).string());
  auto write_node_relocation_stage = [&](const std::string& stage_name,
                                         const int stage_id,
                                         const int iter,
                                         const double stage_time,
                                         const auto& position_of,
                                         const std::vector<Tddd>& stage_delta,
                                         const std::vector<Tddd>& stage_nearest,
                                         const std::unordered_map<networkLine*, Tddd>& line_correction) {
    // stage ごとの補正量を実際の BUFFER に載せる。
    // debug 出力を切っていても，現在反復の補正量はオブジェクト側に残す。
    setNodeRelocationStageBuffer(points, stage_delta, line_correction);
    if (!enable_stage_debug)
      return;
    std::unordered_map<networkPoint*, Tddd> point_correction;
    std::unordered_map<networkPoint*, Tddd> point_nearest;
    point_correction.reserve(n_points);
    point_nearest.reserve(n_points);
    for (std::size_t i = 0; i < n_points; ++i) {
      point_correction[points[i]] = stage_delta[i];
      point_nearest[points[i]] = stage_nearest[i];
    }
    std::unordered_map<networkLine*, Tddd> line_nearest;
    line_nearest.reserve(lines.size());
    for (const auto* l : lines) {
      auto* mutable_line = const_cast<networkLine*>(l);
      auto [pA, pB] = l->getPoints();
      Tddd nearest = vectorToNextSurface(mutable_line, 0.5 * (position_of(pA) + position_of(pB)));
      if (!isFinite(nearest))
        nearest.fill(0.);
      line_nearest[mutable_line] = nearest;
    }
    const std::string filename = net.getName() + "_node_relocation_" + std::to_string(iter) + "_" + stage_name + ".vtu";
    writeNodeRelocationQuadraticDebugVTU(
        *debug_output_directory / filename,
        net.getBoundaryFaces(),
        std::function<Tddd(const networkPoint*)>(position_of),
        point_correction,
        point_nearest,
        line_correction,
        line_nearest,
        static_cast<double>(stage_id));
    node_relocation_debug_pvd.push(filename, stage_time);
  };

  auto compute_point_stage_vectorToNextSurface = [&](const auto& position_of, std::vector<Tddd>& point_nearest) {
    if (!enable_stage_debug)
      return;
    _Pragma("omp parallel for") for (std::size_t i = 0; i < n_points; ++i) {
      Tddd nearest = vectorToNextSurface(points[i], position_of(points[i]));
      if (!isFinite(nearest))
        nearest.fill(0.);
      point_nearest[i] = nearest;
    }
  };

  /* -------------------------------------------------------------------------- */

  auto shiftable_distance = [&](const networkPoint* p, std::size_t idx) -> double {
    double ret = 1E+20;
    auto tet = p->getTetras();
    if (!tet.empty()) {
      for (const auto& t : tet) {
        auto [p0, p1, p2, p3] = t->Points;
        // double dist = Inradius(get_cached_position(p0), get_cached_position(p1), get_cached_position(p2), get_cached_position(p3));
        double dist = (Norm(get_cached_position(p0) - get_cached_position(p1)) + Norm(get_cached_position(p1) - get_cached_position(p2)) + Norm(get_cached_position(p2) - get_cached_position(p0)) + Norm(get_cached_position(p0) - get_cached_position(p3)) + Norm(get_cached_position(p1) - get_cached_position(p3)) + Norm(get_cached_position(p2) - get_cached_position(p3))) / 6.0;
        ret = std::min(ret, dist);
      }
    } else {
      for (const auto& f : p->getBoundaryFaces()) {
        auto [p0, p1, p2] = f->getPoints();
        // double dist = Inradius(get_cached_position(p0), get_cached_position(p1), get_cached_position(p2));
        double dist = (Norm(get_cached_position(p0) - get_cached_position(p1)) + Norm(get_cached_position(p1) - get_cached_position(p2)) + Norm(get_cached_position(p2) - get_cached_position(p0))) / 3.0;
        ret = std::min(ret, dist);
      }
    }
    const double a = 0.3; // 安全率
    return a * ret;
  };

  /* ----------------------- 要素を整えるためのベクトル -------------------------- */

  double coef = coefIN;
  constexpr double convergence_tol = 1e-10; // 全体収束判定の閾値

  auto noThroughCondition = [&](const networkPoint* p, Tddd vec) {
    if (p->Neumann || p->CORNER)
      vec = Chop(vec, RK_with_Ubuff_Normal(p));
    return vec;
  };

  // Anderson acceleration用のフラット化されたベクトル
  // X = [x0, y0, z0, x1, y1, z1, ...]
  const std::size_t vec_size = n_points * 3;
  AndersonAcceleration<std::vector<double>> anderson(5); // history_size = 5
  const bool enable_anderson = []() {
    // true quadratic node relocation is more sensitive to overshoot in accelerated fixed-point updates.
    // Keep Anderson off by default there; allow explicit opt-in.
    if (node_relocation_surface == NodeRelocationSurface::true_quadratic) {
      if (const char* env = std::getenv("BEM_ALE_TRUEQ_ANDERSON"))
        return std::string(env) == "1";
      return false;
    }
    return true;
  }();

  // vecToSurfaceをフラットベクトルに変換
  auto flattenVecToSurface_into = [&](std::vector<double>& X) {
    for (std::size_t i = 0; i < n_points; ++i) {
      X[3 * i + 0] = points[i]->vecToSurface[0];
      X[3 * i + 1] = points[i]->vecToSurface[1];
      X[3 * i + 2] = points[i]->vecToSurface[2];
    }
  };

  // フラットベクトルからvecToSurfaceに復元
  auto unflatten = [&](const std::vector<double>& X) {
    for (std::size_t i = 0; i < n_points; ++i) {
      points[i]->vecToSurface[0] = X[3 * i + 0];
      points[i]->vecToSurface[1] = X[3 * i + 1];
      points[i]->vecToSurface[2] = X[3 * i + 2];
    }
  };

  std::vector<double> X_curr(vec_size), X_next(vec_size), F(vec_size);

  std::vector<std::tuple<networkLine*, networkPoint*, Tddd, networkPoint*, Tddd>> linear_midpoint_endpoint_corrections; // (line, pA, penaltyA, pB, penaltyB)
  linear_midpoint_endpoint_corrections.reserve(lines.size());
  for (const auto* l : lines) {
    if (l->Neumann || l->CORNER)
      linear_midpoint_endpoint_corrections.emplace_back(const_cast<networkLine*>(l), l->Point_A, std::array<double, 3>{0., 0., 0.}, l->Point_B, std::array<double, 3>{0., 0., 0.});
  }

  for (auto steps = 0; steps < loop; ++steps) {
    coef = coefIN;
    std::vector<Tddd> quality_stage_delta(n_points, Tddd{0., 0., 0.});
    std::vector<Tddd> quality_stage_nearest(n_points, Tddd{0., 0., 0.});
    std::vector<Tddd> midpoint_stage_delta(n_points, Tddd{0., 0., 0.});
    std::vector<Tddd> midpoint_stage_nearest(n_points, Tddd{0., 0., 0.});
    std::vector<Tddd> cling_stage_delta(n_points, Tddd{0., 0., 0.});
    std::vector<Tddd> cling_stage_nearest(n_points, Tddd{0., 0., 0.});

    // 現在のvecToSurfaceを保存
    flattenVecToSurface_into(X_curr);

    /* -------------------------------------------------------------------------- */

    // キャッシュを更新
    update_cache();

    _Pragma("omp parallel for") for (std::size_t i = 0; i < n_points; ++i) {
      const auto& p = points[i];
      auto current_pX = cached_positions[i];
      Tddd V = DistorsionMeasureWeightedSmoothingVector_modified(p, current_pX, get_cached_position);
      V = noThroughCondition(p, V);
      double v_norm = Norm(V);
      if (v_norm > 1E-12) {
        double shift_limit = coef * shiftable_distance(p, i);
        p->vecToSurface_BUFFER = std::min(v_norm, shift_limit) * Normalize(V);
      } else
        p->vecToSurface_BUFFER.fill(0.);
      quality_stage_delta[i] = p->vecToSurface_BUFFER;
      points[i]->vecToSurface += p->vecToSurface_BUFFER;
    }
    update_cache();
    compute_point_stage_vectorToNextSurface([&](const networkPoint* p) { return RK_with_Ubuff(p); }, quality_stage_nearest);
    write_node_relocation_stage("quality", 0, steps, 3.0 * steps + 0.0, [&](const networkPoint* p) { return RK_without_Ubuff(p) + p->vecToSurface; }, quality_stage_delta, quality_stage_nearest, {});

    /* -------------------------------------------------------------------------- */

    // キャッシュを更新（clungSurface計算用）
    update_cache();

    _Pragma("omp parallel for") for (std::size_t i = 0; i < n_points; ++i) {
      const auto& p = points[i];
      p->clungSurface = vectorToNextSurface(p, RK_with_Ubuff(p));
      if (!isFinite(p->clungSurface))
        p->clungSurface.fill(0.);
    }

    for (std::size_t i = 0; i < n_points; ++i) {
      const auto& p = points[i];
      double clung_norm = Norm(p->clungSurface);
      if (clung_norm > 1E-12) {
        double shift_limit = shiftable_distance(p, i);
        cling_stage_delta[i] = std::min(clung_norm, shift_limit) * Normalize(p->clungSurface);
        p->vecToSurface += cling_stage_delta[i];
      }
      p->clungSurface.fill(0.);
    }
    update_cache();
    compute_point_stage_vectorToNextSurface([&](const networkPoint* p) { return RK_with_Ubuff(p); }, cling_stage_nearest);
    write_node_relocation_stage("cling", 2, steps, 3.0 * steps + 2.0, [&](const networkPoint* p) { return RK_without_Ubuff(p) + p->vecToSurface; }, cling_stage_delta, cling_stage_nearest, {});

    /* -------------------------------------------------------------------------- */

    // Anderson acceleration: 残差 F = G(X) - X を計算
    flattenVecToSurface_into(X_next);
    for (std::size_t i = 0; i < vec_size; ++i)
      F[i] = X_next[i] - X_curr[i];

    // Anderson accelerationで次の推定値を計算
    if (enable_anderson) {
      std::vector<double> X_accelerated = anderson.compute_next(X_curr, F);
      bool finite_ok = true;
      for (const auto& v : X_accelerated) {
        if (!std::isfinite(v)) {
          finite_ok = false;
          break;
        }
      }
      if (finite_ok)
        unflatten(X_accelerated);
      else
        unflatten(X_next); // fallback to non-accelerated update
    } else {
      unflatten(X_next);
    }

    /* -------------------------------------------------------------------------- */
  }

  if (enable_stage_debug) {
    node_relocation_debug_pvd.output();
  }

  /* -------------------------------------------------------------------------- */
  /*  midpoint relocation: 頂点の品質改善後、midpointを面上にスナッピング         */
  /*  linear 辺も含む: 端点平均が凸壁を貫入するのを防止                          */
  /* -------------------------------------------------------------------------- */
  for (auto* l : net.getBoundaryLines()) {
    l->vecToSurface.fill(0.);
    l->clungSurface.fill(0.);
  }

  // midpointの面スナッピング: 頂点の補正に追従 + Neumann面制約
  int corner_mid_count = 0, dirichlet_mid_count = 0;
  for (auto* l : net.getBoundaryLines()) {
    auto [pA, pB] = l->getPoints();
    auto X_shidted = 0.5 * (RK_with_Ubuff(pA) + RK_with_Ubuff(pB));

    l->clungSurface = vectorToNextSurface(l, X_shidted);
    if (!isFinite(l->clungSurface))
      l->clungSurface.fill(0.);
    l->vecToSurface = l->clungSurface + (X_shidted - RK_without_Ubuff(l));

    // Diagnostic: measure how far midpoint is shifted from linear center
    double edge_len = Norm(RK_with_Ubuff(pA) - RK_with_Ubuff(pB));
    Tddd final_pos = RK_without_Ubuff(l) + l->vecToSurface;
    double shift = Norm(final_pos - X_shidted);
    double relative_shift = (edge_len > 1e-20) ? shift / edge_len : 0.;
    if (l->CORNER) {
      ++corner_mid_count;
    } else if (l->Dirichlet) {
      ++dirichlet_mid_count;
    }
  }
};

// ! -------------------------------------------------------------------------- */
// !                               calculateVelocities                          */
// ! -------------------------------------------------------------------------- */

using FuncU = std::function<Tddd(const networkPoint*)>;

inline void set_u_potential_BEM(const Network& net) {
  _Pragma("omp parallel for") for (const auto& p : ToVector(net.getPoints())) {
    p->u_potential_BEM = gradPhi(p);
  }

  // Compute edge midpoint potential velocity.
  // For true quadratic elements, use the proper quadratic gradient at the midpoint
  // instead of the linear average of endpoint velocities.
  for (auto* l : net.getBoundaryLines()) {
    auto [pA, pB] = l->getPoints();
    l->u_potential_BEM = 0.5 * (pA->u_potential_BEM + pB->u_potential_BEM); // default: linear
    if (l->hasActiveBieDof())
      l->u_potential_BEM = gradPhi(l);
  }
}

inline void set_u_potential_BEM(Network* net) { set_u_potential_BEM(*net); }

inline void set_u_potential_BEM(const std::vector<Network*>& nets) {
  for (auto* net : nets)
    set_u_potential_BEM(*net);
}

inline double DphiDt_at_midpoint(const networkLine* l) {
  double aphiat = -0.5 * Dot(l->u_potential_BEM, l->u_potential_BEM) - _GRAVITY_ * l->X_mid[2];
  return aphiat + Dot(l->u_reloc, l->u_potential_BEM);
}

/* -------------------------------------------------------------------------- */

inline void set_u_total(const Network& net) {

  _Pragma("omp parallel for") for (const auto& p : ToVector(net.getPoints())) { p->u_total = p->u_potential_BEM + p->u_omega_VPM; }

  for (auto* l : net.getBoundaryLines()) {
    auto [pA, pB] = l->getPoints();
    l->u_total = 0.5 * (pA->u_total + pB->u_total); // default: linear
    if (l->hasActiveBieDof())
      l->u_total = l->u_potential_BEM + l->u_omega_VPM;
  }
}

inline void set_u_total(Network* net) { set_u_total(*net); }

inline void set_u_total(const std::vector<Network*>& nets) {
  for (auto* net : nets)
    set_u_total(*net);
}

/* -------------------------------------------------------------------------- */

inline void setNodeVelocity(const Network& net,
                            const int loop = 0,
                            const double coef = 0.,
                            const std::filesystem::path* debug_output_directory = nullptr) {
  for (const auto& p : ToVector(net.getPoints())) {
    p->u_reloc = p->u_total;
    // CAUTIION : 以下のようにNeuamnnだけことなる時間発展の方法を使うのはよくないようだ．
    // 基本的にラグランジュ的に時間発展させ，修正も同じ状態に対して行うのが安全だ．
    // CORNERで見られた，内部との若干のズレもこのあたりが原因の可能性がある．
    // if (p->Neumann)
    //   p->u_reloc = velocity_of_Body(std::get<0>(getEffectiveNearestContactFace(p)), getPosition(p));
    p->vecToSurface.fill(0.);
    p->vecToSurface_BUFFER.fill(0.);
    p->vecToSurface_BUFFER_BUFFER.fill(0.);
  }

  for (const auto& l : net.getBoundaryLines()) {
    // もし，lineがtrue_quadratic_elementを持っていなるなら，意味がある．linear elementの場合は，そもそもu_relocは使わない．
    // true quadraticの場合，ラグランジュ的時間発展で予測される曲面は，２次補間で予測される曲面でなくては，２次補間を利用する意味がない．
    l->u_reloc = l->u_total;
    // if (l->Neumann) {
    //   auto [f, X, dist] = getEffectiveNearestContactFace(l);
    //   if (f)
    //     l->u_reloc = velocity_of_Body(f, getPosition(l));
    // }
    l->vecToSurface.fill(0.);
  }

  for (const auto& p : net.getPoints()) {
    if (!isFinite(p->u_reloc)) {
      std::cout << "p->X = " << p->X << std::endl;
      std::cout << "p->u_potential_BEM = " << p->u_potential_BEM << ", gradPhi(p) = " << gradPhi(p) << std::endl;
      std::cout << "p->u_reloc = " << p->u_reloc << std::endl;
      std::cout << "p->vecToSurface = " << p->vecToSurface << std::endl;
      std::cout << "p->Dirichlet = " << p->Dirichlet << std::endl;
      std::cout << "p->Neumann = " << p->Neumann << std::endl;
      std::cout << "p->CORNER = " << p->CORNER << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    }
  }

  if (loop <= 0 || !(coef > 0.))
    return;

  calculateVecToSurface(net, loop, coef, debug_output_directory);

  // Compute relocation velocity/target from vecToSurface.
  // ALE:
  //   use u_reloc during RK4 push so position and phi are advanced with the
  //   same modified velocity.
  // interpolation:
  //   keep RK4 evolution on u_total, and use X_reloc only after RK4 to move
  //   nodes and re-interpolate phi.
  // This block does not modify the RK push choice; it only updates u_reloc and
  // X_reloc from the current relocation target.
  auto setRelocFromVecToSurface = [](auto* entity) {
    const Tddd X_target = entity->vecToSurface + RK_without_Ubuff(entity);
    const Tddd u_target = entity->RK_X.getVectorToReachAtNextTimeQ(X_target);
    if (isFinite(u_target)) {
      entity->u_reloc = u_target;
      entity->X_reloc = entity->RK_X.getX(u_target);
    } else {
      entity->u_reloc = isFinite(entity->u_total) ? entity->u_total : Tddd{0., 0., 0.};
      entity->X_reloc = RK_without_Ubuff(entity);
    }
  };

  for (const auto& p : net.getPoints())
    setRelocFromVecToSurface(p);

  for (auto* l : net.getBoundaryLines())
    if (l->hasActiveBieDof())
      setRelocFromVecToSurface(l);

  dumpDebugMidpointLineState(&net, "post-ale-relocation", -1, -1);
}

inline void setNodeVelocity(const std::vector<Network*>& nets,
                            const int loop = 0,
                            const double coef = 0.,
                            const std::filesystem::path* debug_output_directory = nullptr) {
  for (auto* net : nets)
    setNodeVelocity(*net, loop, coef, debug_output_directory);
}

/* ========================================================================== */

/*DOC_EXTRACT 0_7_OTHERS

### エネルギー保存則（計算精度のチェックに利用できる）

流体全体の運動エネルギーは，ラプラス���程式と発散定理を使うと，次のように境界面に沿った積分で表される．

$$
E_K =\frac{\rho}{2} \iint_\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
$$

また，流体の位置エネルギーは，次のように表される．

$$
E_P = \frac{\rho}{2} \iint_\Gamma (0,0,g(z - z_0)^2) \cdot {\bf n} d\Gamma
$$


???   NOTE: なぜか？

      テンソルを使って考えてみると

      $$
      \begin{align*}
      \nabla \cdot (\phi\nabla\phi) &= \frac{\partial\phi}{\partial x_i} \frac{\partial\phi}{\partial x_i} + \phi \frac{\partial^2\phi}{\partial x_i \partial x_i}\\
      &= \nabla \phi \cdot \nabla \phi + \phi \nabla^2 \phi\\
      &= \nabla \phi \cdot \nabla \phi
      \end{align*}
      $$

      よって，

      $$
      \iiint_\Omega \nabla\phi\cdot\nabla\phi d\Omega = \iiint_\Omega \nabla \cdot (\phi\nabla\phi) d\Omega = \iint_\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
      $$

      ---

      $$
      E_P = \rho g \iiint_\Omega (z - z_0) d\Omega
      = \rho g \iiint_\Omega \frac{1}{2} \nabla \cdot (0,0,(z - z_0)^2) d\Omega
      = \rho g \iint_\Gamma \frac{1}{2} (0,0,(z - z_0)^2) \cdot {\bf n} d\Gamma
      = \frac{1}{2}\rho g \iint_\Gamma (z - z_0)^2 n_z d\Gamma
      $$

      ---

*/

inline double KinematicEnergy(const auto& faces) {
  double EK = 0;
  for (const auto& f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    EK += Dot((p0->phiphin[0] + p1->phiphin[0] + p2->phiphin[0]) / 3. * gradPhi(f), f->normal) * f->area;
  }
  return _WATER_DENSITY_ * EK / 2.;
};

inline double PotentialEnergy(const auto& faces) {
  double EP = 0;
  for (const auto& f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    auto intpX = interpolationTriangleLinear0101(T3Tddd{p0->X, p1->X, p2->X});
    for (const auto& [x0, x1, w0w1] : __GWGW10__Tuple)
      EP += std::pow(intpX(x0, x1)[2], 2) * f->normal[2] * w0w1 * intpX.J(x0, x1);
  }
  return _WATER_DENSITY_ * EP / 2.;
};

inline double TotalEnergy(const auto& faces) { return KinematicEnergy(faces) + PotentialEnergy(faces); };

/*DOC_EXTRACT 0_7_OTHERS

### 内部流速の計算方法（使わなくてもいい）

[Fochesato2005](https://onlinelibrary.wiley.com/doi/10.1002/fld.838)にあるように，
流体内部の流速$\nabla \phi$は，BIEを微分して求めることができる．

$$
u({\bf a}) = \nabla\phi({\bf a}) = \int_{\partial \Omega} \frac{\partial Q}{\partial n} ({\bf x})Q({\bf x}, {\bf a}) - \phi({\bf x}) \frac{\partial Q}{\partial n} ({\bf x}, {\bf a}) d\Gamma
$$

$$
Q({\bf x},{\bf a}) = \frac{{\bf r}}{4\pi r^3}, \quad \frac{\partial Q}{\partial n} ({\bf x},{\bf a}) = \frac{1}{4\pi r^3} (3 \mathbf{n} - (\mathbf{r} \cdot \mathbf{n}) \frac{\mathbf{r}}{r^2})
$$

*/

#endif
