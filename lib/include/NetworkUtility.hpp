#ifndef NetworkUtility_H
#define NetworkUtility_H

#include "Network.hpp"

namespace network_utility_detail {

inline bool isLineFaceDetachedByPressureLocal(const networkLine* l, const networkFace* f) {
  if (!l || !f)
    return false;
  const auto* d = l->findContactState(f);
  return d && d->detached_by_pressure;
}

inline bool isLineFaceNeumannLocal(const networkLine* l, const networkFace* f) {
  if (!l || !f || !f->BoundaryQ())
    return false;
  if (isLineFaceDetachedByPressureLocal(l, f))
    return false;
  const auto* d = l->findContactState(f);
  return d && d->nearestContactFace() != nullptr;
}

inline bool lineIsDirichletLocal(const networkLine* l) {
  if (!l)
    return false;
  const auto bfs = l->getBoundaryFaces();
  if (bfs.empty())
    return false;
  return std::ranges::all_of(bfs, [&](const auto* f) { return !isLineFaceNeumannLocal(l, f); });
}

}  // namespace network_utility_detail

/* -------------------------------------------------------------------------- */

inline void writeVertices(std::ofstream &ofs, std::map<networkPoint *, int> &P_i, const auto &points) {
  int i = 0;
  for (const auto &p : points) {
    P_i[p] = ++i;
    auto [X0, X1, X2] = ToX(p);
    ofs << "v " << X0 << " " << X1 << " " << X2 << std::endl;
  }
  ofs << std::endl;
}

inline void writeFaces(std::ofstream &ofs, const std::map<networkPoint *, int> &P_i, const auto &faces) {
  for (const auto &f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    ofs << "f " << P_i.at(p0) << " " << P_i.at(p1) << " " << P_i.at(p2) << std::endl;
  }
}

inline void createOBJ(std::ofstream &ofs, Network &net) {
  std::map<networkPoint *, int> pointIndices;
  writeVertices(ofs, pointIndices, net.getPoints());
  writeFaces(ofs, pointIndices, net.getBoundaryFaces());
}

//% ------------------------------------------------------ */
//%                         反射点の計算                     */
//% ------------------------------------------------------ */

// bool isInConvexPolygon(const networkPoint *const p) {
//    try {
//       auto ps = p->getNeighborsSort();  // ソートできない場合エラーとなる
//       if (ps.size() < 2)
//          return false;
//       if (ps.size() == 3)
//          return true;
//       if (ps.size() != p->getNeighbors().size())
//          return false;  // 面が重複する場合を覗くために設置

//       return isConvexPolygon(extX(ps), p->getNormalTuple());
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorReset << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };

// bool isStraight(const V_d &v0, const V_d &v1, const double angle = 1E-2 /*閾値*/) {
//    //        ^
//    //      v1|angle (positive)
//    // ---v0-->----
//    //      v1|angle (negative)
//    //        v
//    if (std::abs(MyVectorAngle(v0, v1)) < angle)
//       return true;
//    return false;
// };

// bool isStraight(const netPp p0, const netPp p1, const netPp p2, const double angle = 1E-2 /*閾値*/) {
//    return (isStraight(ToVector(ToX(p0) - ToX(p1)), ToVector(ToX(p1) - ToX(p2))) < angle);
// };

inline bool isFlat(const netPp p, double minangle = M_PI / 180.) {
  if (p->getLines().empty())
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p->getLines().empty()");
  auto faces = p->getBoundaryFaces();
  for (auto i = 0; i < faces.size(); i++)
    for (auto j = i; j < faces.size(); j++)
      if (!isFlat(ToX(faces[i]), ToX(faces[j]), minangle))
        return false;
  return true;
};

inline bool isFlat(const netLp line, double minangle = M_PI / 180.) {
  auto fs = line->getBoundaryFaces();
  if (fs.size() != 2)
    return false;
  else
    return isFlat(ToX(fs[0]), ToX(fs[1]), minangle);
};

/* -------------------------------------------------------------------------- */

inline Tddd exact_along_surface(const networkPoint *const p, Tddd VECTOR) {
  for (const auto &f : p->getBoundaryFaces())
    VECTOR = Chop(VECTOR, f->normal);
  return VECTOR;
  // auto faces = ToVector(p->getFaces());
  // std::vector<std::tuple<Tddd, Tddd>> normals;
  // normals.reserve(2 * faces.size());
  // for (auto i = 0; i < faces.size(); i++)
  //    for (auto j = i; j < faces.size(); j++)
  //       if (!isFlat(faces[i]->normal, faces[j]->normal, 1E-2 * M_PI / 180.))
  //          normals.push_back({faces[i]->normal, faces[j]->normal});
  // for (const auto &[N0, N1] : normals)
  //    VECTOR = Projection(VECTOR, Cross(N0, N1));
  // return VECTOR;
};

inline bool isFlat(const auto p, const double lim_rad = 1E-2 * M_PI / 180.) {
  auto faces = p->getBoundaryFaces();
  for (auto i = 0; i < faces.size(); i++)
    for (auto j = i; j < faces.size(); j++)
      if (!isFlat(faces[i]->normal, faces[j]->normal, lim_rad))
        return false;
  return true;
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT smoothing_vector

### 格子の平滑化

境界面を格子を綺麗な状態に整えるために，境界面の点を移動させる．
移動方法は，様々考えられる．

| 移動方法 | 説明 |
|:-------|:---------|
| Laplacian smoothing                     | 境界面の点をその点の近傍の点の重心に移動させる．          |
| Area weighted smoothing                 | 隣接面の面積の重みを掛けて，面の中心に移動させる．        |
| Distorsion measure weighted smoothing   | 隣接面の歪みに関する係数を重みとして掛けて，面の中心に移動させる．|

共通点は，許されない移動を防止する移動前に，移動後の形状をチェックすることである．
`canFlip`と同様，\ref{isValidTriangle}{`isValidTriangle`}を用いる．

徐々に移動させる場合，誤差の蓄積，条件の変化を把握するのが難しい．
大きな変化は防げても，小さな変化には対応できない場合が考えられる．

*/
// Assuming Tddd is a type that can be used in a constexpr context

inline std::array<std::tuple<double, double, Tddd>, 25> t0t1ModTriShape25(std::array<double, 2> t0_range, std::array<double, 2> t1_range) {
  auto subdiv_for_t1 = Subdivide<4>(t1_range[0], t1_range[1]);
  std::array<std::tuple<double, double, Tddd>, 25> result{};
  size_t index = 0;
  for (const auto &t0 : Subdivide<4>(t0_range[0], t0_range[1])) {
    for (const auto &t1 : subdiv_for_t1) {
      if (index < 25) {
        result[index++] = std::make_tuple(t0, t1, ModTriShape<3>(t0, t1));
      }
    }
  }
  return result;
}

inline std::array<std::tuple<double, double, Tddd>, 400> t0t1ModTriShape400(std::array<double, 2> t0_range, std::array<double, 2> t1_range) {
  auto subdiv_for_t1 = Subdivide<19>(t1_range[0], t1_range[1]);
  std::array<std::tuple<double, double, Tddd>, 400> result{};
  size_t index = 0;
  for (const auto &t0 : Subdivide<19>(t0_range[0], t0_range[1])) {
    for (const auto &t1 : subdiv_for_t1) {
      if (index < 40) {
        result[index++] = std::make_tuple(t0, t1, ModTriShape<3>(t0, t1));
      }
    }
  }
  return result;
}

// AppleClang/libc++ cannot constant-evaluate some of the helper routines used here.
#if defined(_LIBCPP_VERSION)
inline auto precalculateXtrys400() {
  const auto subdiv_for_t0 = Subdivide<199>(0.8, 1.);
  const auto subdiv_for_t1 = Subdivide<199>(0., 1.);
  std::array<std::tuple<double, double, Tddd>, 400> result{};

  size_t index = 0;
  for (const auto &t0 : subdiv_for_t0) {
    for (const auto &t1 : subdiv_for_t1) {
      if (index < 400) {
        result[index++] = std::make_tuple(t0, t1, ModTriShape<3>(t0, t1));
      }
    }
  }

  return result;
}

inline const auto evaluation_points_for_findOptimalPositionByCriteria400 = precalculateXtrys400();
#else
inline constexpr std::array<std::tuple<double, double, Tddd>, 400> precalculateXtrys400() {
  constexpr auto subdiv_for_t0 = Subdivide<199>(0.8, 1.);
  constexpr auto subdiv_for_t1 = Subdivide<199>(0., 1.);
  std::array<std::tuple<double, double, Tddd>, 400> result{};

  size_t index = 0;
  for (const auto &t0 : subdiv_for_t0) {
    for (const auto &t1 : subdiv_for_t1) {
      if (index < 400) {
        result[index++] = std::make_tuple(t0, t1, ModTriShape<3>(t0, t1));
      }
    }
  }

  return result;
}

// Variable Initialization
inline constexpr auto evaluation_points_for_findOptimalPositionByCriteria400 = precalculateXtrys400();
#endif

inline Tddd findOptimalPositionByCriteriaOn(const networkPoint *p, const std::function<double(const Tddd &)> &criteria, std::function<Tddd(const networkPoint *)> get_triangle_base_X /*optimum position will be found on this triangle*/) {

  constexpr double acceptable_change_angle = 1E-2 * M_PI / 180.;
  const auto all_face = p->getFaces();
  std::array<double, 3> where_the_min_score = get_triangle_base_X(p);
  double min_score = criteria(where_the_min_score), s;
  if (!isFinite(min_score))
    min_score = 1E+10;
  bool found = false;

  auto check_if_next_faces_valid = [&](const Tddd &replace_X) {
    return std::ranges::all_of(all_face, [&](auto f) {
      auto [p0, p1, p2] = f->getPoints(p);
      Tddd base_X0 = get_triangle_base_X(p0), base_X1 = get_triangle_base_X(p1), base_X2 = get_triangle_base_X(p2);
      if (Dot(TriangleNormal(base_X0, base_X1, base_X2), TriangleNormal(replace_X, base_X1, base_X2)) < 0. /* opposite direction */)
        return false;
      // if (isFlat(TriangleNormal(base_X0, base_X1, base_X2), TriangleNormal(replace_X, base_X1, base_X2), acceptable_change_angle))
      if (p0 != p && p1 != p && p2 != p)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p0 != p && p1 != p && p2 != p");
      return true;
    });
  };

  double min_t0, min_t1;
  networkFace *min_face = all_face[0];
  auto check_and_update_min_score = [&](const Tddd &X0_candidate, const double t0, const double t1, networkFace *f) {
    if (check_if_next_faces_valid(X0_candidate)) {
      s = criteria(X0_candidate);
      if (isFinite(s) && min_score > s) {
        min_score = s;
        where_the_min_score = X0_candidate;
        min_t0 = t0;
        min_t1 = t1;
        found = true;
        min_face = f;
      }
    }
  };

  int itr = 5;
  T3Tddd base_X012;

  auto faces_will_be_updated = p->getFaces();
  for (const auto &f : faces_will_be_updated) {
    std::array<double, 2> t0_minmax = {0.9, 1.};
    std::array<double, 2> t1_minmax = {0., 1.};
    auto subdiv_for_t1 = Subdivide<5>(t1_minmax[0], t1_minmax[1]);
    auto subdiv_for_t0 = Subdivide<5>(t0_minmax[0], t0_minmax[1]);
    for (int i = 0; i < itr; ++i) {
      auto [p0, p1, p2] = f->getPoints(p);
      base_X012 = {get_triangle_base_X(p0), get_triangle_base_X(p1), get_triangle_base_X(p2)};
      for (const auto &t0 : subdiv_for_t0)
        for (const auto &t1 : subdiv_for_t1)
          check_and_update_min_score(Dot(ModTriShape<3>(t0, t1), base_X012), t0, t1, f);
      // faces_will_be_updated = {min_face};
      auto del_t = std::pow(0.2, i + 2);
      t0_minmax = {std::clamp(min_t0 - del_t, 0.0, 1.0), std::clamp(min_t0 + del_t, 0.0, 1.0)};
      t1_minmax = {std::clamp(min_t1 - del_t, 0.0, 1.0), std::clamp(min_t1 + del_t, 0.0, 1.0)};
    }
  }

  return where_the_min_score;
};

inline double calculateScoreIfPointMoved(const networkPoint *p, const Tddd &X_candidate, std::function<Tddd(const networkPoint *)> get_position_to_be_updated) {
  double total = 0.;
  for (const auto &f : p->getFaces()) {
    auto [_, p1, p2] = f->getPoints(p);
    Tddd X1 = get_position_to_be_updated(p1), X2 = get_position_to_be_updated(p2);
    total += CircumradiusToInradius(X_candidate, X1, X2);
  }
  return total * total;
};

inline void SmoothingPreserveShape(networkPoint *p) {
  auto calculateScoreIfPointMoved = [](const networkPoint *p, const Tddd &X_candidate, std::function<Tddd(const networkPoint *)> get_position_to_be_updated) {
    const auto p_faces = p->getFaces();
    double n = p_faces.size(), area;
    const double acceptable_change_angle = 1E-2 * M_PI / 180.;
    double total = 0;
    for (const auto &f : p_faces) {
      auto [p0, p1, p2] = f->getPoints(p);
      Tddd X1 = get_position_to_be_updated(p1);
      Tddd X2 = get_position_to_be_updated(p2);
      total += CircumradiusToInradius(X_candidate, X1, X2) * std::pow(Distance(X_candidate, X1) + Distance(X_candidate, X2) + Distance(X1, X2), 2) / TriangleArea(X_candidate, X1, X2);
      total += differenceFromEquilateralTriangle(X_candidate, X1, X2);
      // total += std::pow(Circumradius(X_candidate, X1, X2), 2) * M_PI / TriangleArea(X_candidate, X1, X2);
      if (!isFlat(X_candidate, p1->X, p2->X, p0->X, p1->X, p2->X, acceptable_change_angle))
        return 1E+30;
    }
    return total;
  };
  auto get_base_X = [](const networkPoint *p) { return p->X; };
  auto criteria = [&](const Tddd &X) { return calculateScoreIfPointMoved(p, X, get_base_X); };
  p->X_temporary = 0.5 * (findOptimalPositionByCriteriaOn(p, criteria, get_base_X) - p->X) + p->X;
};

// void SmoothingPreserveShape(netPp p, const std::function<Tddd(const netPp, const std::array<double, 3> &)> &SmoothingVector) {
//    try {
//       if (!isEdgePoint(p)) { /*端の点はsmoothingしない*/
//          auto ps = p->getNeighbors();
//          // 無視できる角度
//          const double acceptable_change_angle = 1E-6 * M_PI / 180.;
//          const auto faces = ToVector(p->getFaces());
//          double Atot = 0.;
//          for (const auto &f : faces)
//             Atot += f->area;
//          const double length_scale = std::sqrt(Atot / faces.size());
//          Tddd X_candidate, X1, X2, cross_before, cross_after, result_V;

//          double min_angle_before = 1E+10;
//          double min_angle_after = 1E+10;
//          auto check_if_next_faces_valid = [&](const Tddd &replace_X) {
//             for (const auto &f : faces) {
//                auto [p0, p1, p2] = f->getPoints();
//                if (p0 != p && p1 != p && p2 != p)
//                   throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p0 != p && p1 != p && p2 != p");

//                X0 = p0->X;
//                X1 = p1->X;
//                X2 = p2->X;

//                auto a = Min(TriangleAngles(X0, X1, X2));
//                if (min_angle_before < a)
//                   min_angle_before = a;

//                cross_before = Cross(X1 - X0, X2 - X0);

//                if (p0 == p)
//                   X0 = replace_X;
//                else if (p1 == p)
//                   X1 = replace_X;
//                else if (p2 == p)
//                   X2 = replace_X;

//                cross_after = Cross(X1 - X0, X2 - X0);

//                a = Min(TriangleAngles(X0, X1, X2));
//                if (min_angle_after < a)
//                   min_angle_after = a;

//                //! どうやらisValidTriangleはいらないようだ．
//                // if (!isValidTriangle(T3Tddd{X0, X1, X2}, 1E-2 * M_PI / 180.))
//                //    return false;
//                if (Dot(cross_before, cross_after) <= 0.)
//                   return false;

//                // if (!isFlat(cross_before, cross_after, acceptable_change_angle))
//                //    return false;
//                if (!isFlat(p0->X, p1->X, p2->X, X0, X1, X2, acceptable_change_angle))
//                   return false;
//             }

//             if (min_angle_before > min_angle_after)
//                return false;

//             return true;
//          };

//          if (ps.size() > 2) {  // 2点の場合は2点の中点に動いてしまうので，実行しない
//             Tddd X_next = p->X, V_try, X_try;
//             const double a_default = 0.01;
//             double a = a_default;
//             for (auto iter = 0; iter < 10; ++iter) {
//                X_try = X_next + a * exact_along_surface(p, SmoothingVector(p, X_next));
//                if (check_if_next_faces_valid(X_try)) {
//                   p->setXSingle(X_try);
//                   X_next = X_try;
//                   a = a_default;
//                } else
//                   a *= 0.5;
//             }
//          }
//       }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorReset << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };

/* ------------------------------------------------------ */

inline Tddd DistorsionMeasureWeightedSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V = {0., 0., 0.};
  if (!isEdgePoint(p)) {
    Tddd X = {0., 0., 0.}, Xmid, vertical;
    double Wtot = 0, W, height;
    for (const auto &f : p->getFaces()) {
      auto [X0, X1, X2] = ToX(f->getPoints(p));
      X0 = current_pX;
      // W = std::log10(CircumradiusToInradius(X0, X1, X2));
      W = std::pow(CircumradiusToInradius(X0, X1, X2), 2);
      W = std::clamp(W, 4., 1000.);
      /*
      CircumradiusToInradius(X0, X1, X2)は最小で2なので，log2は1以上の値をとる．
      */
      // W = std::log2(CircumradiusToInradius(X0, X1, X2));
      // W = std::clamp(W, 1., 2.);
      Wtot += W;
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3.) * 0.5;
      FusedMultiplyIncrement(W, height * vertical + Xmid, V);
      // X = height * vertical + Xmid;
      // V += W * X;
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

inline Tddd DistorsionMeasureWeightedSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX, std::function<Tddd(const networkPoint *)> position) {
  Tddd V = {0., 0., 0.};
  if (!isEdgePoint(p)) {
    Tddd X = {0., 0., 0.}, Xmid, vertical;
    double Wtot = 0, W, height;
    for (const auto &f : p->getFaces()) {
      auto [_, p1, p2] = f->getPoints(p);
      auto X0 = current_pX;
      auto X1 = position(p1);
      auto X2 = position(p2);
      // W = std::log10(CircumradiusToInradius(X0, X1, X2));
      W = std::pow(CircumradiusToInradius(X0, X1, X2), 2);
      W = std::clamp(W, 4., 20.);
      /*
      CircumradiusToInradius(X0, X1, X2)は最小で2なので，log2は1以上の値をとる．
      */
      // W = std::log2(CircumradiusToInradius(X0, X1, X2));
      // W = std::clamp(W, 1., 2.);
      Wtot += W;
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3.) * 0.5;
      FusedMultiplyIncrement(W, height * vertical + Xmid, V);
      // X = height * vertical + Xmid;
      // V += W * X;
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

// void DistorsionMeasureWeightedSmoothingPreserveShape(const auto &ps, const int iteration = 1) {
//    for (auto i = 0; i < iteration; ++i)
//       for (const auto &p : ps) SmoothingPreserveShape(p, DistorsionMeasureWeightedSmoothingVector);
// };

/* ------------------------------------------------------ */

inline Tddd AreaWeightedSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V = {0., 0., 0.};
  if (!isEdgePoint(p)) {
    double Wtot = 0, W;
    for (const auto &f : p->getFaces()) {
      auto [X0, X1, X2] = ToX(f->getPoints(p));
      X0 = current_pX;
      Wtot += (W = TriangleArea(X0, X1, X2));
      FusedMultiplyIncrement(W, Incenter(X0, X1, X2), V);
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

// void AreaWeightedSmoothingPreserveShape(const auto &ps, const int iteration = 1) {
//    for (auto i = 0; i < iteration; ++i)
//       for (const auto &p : ps)
//          SmoothingPreserveShape(p, AreaWeightedSmoothingVector);
// };
/* ------------------------------------------------------ */

inline Tddd ArithmeticWeightedIncenterSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V{};
  double Wtot = 0;
  if (!isEdgePoint(p)) {
    for (const auto &q : p->getNeighbors()) {
      V += q->X;
      Wtot += 1.;
    }
    for (const auto &face : p->getFaces()) {
      auto [X0, X1, X2] = ToX(face->getPoints(p));
      X0 = current_pX;
      V += Incenter(X0, X1, X2);
      Wtot += 1.;
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

// void ArithmeticWeightedIncenterSmoothingPreserveShape(const auto &ps, const int iteration = 1) {
//    for (auto i = 0; i < iteration; ++i)
//       for (const auto &p : ps) SmoothingPreserveShape(p, ArithmeticWeightedIncenterSmoothingVector);
// };

/* ------------------------------------------------------ */

inline Tddd ArithmeticWeightedSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX, const double power = 1.) {
  Tddd V = {0., 0, 0.};
  double Wtot = 0, norm;
  if (!isEdgePoint(p)) {
    for (const auto &q : p->getNeighbors()) {
      V += std::pow(Norm(q->X - current_pX), power - 1.) * (q->X - current_pX);

      // this is equivalent to the following
      // norm = Norm(q->X - current_pX);
      // V += std::pow(norm, power) * (q->X - current_pX) / norm;

      Wtot += 1.;
    }
    return V / Wtot;
  } else
    return V;
};

// void LaplacianSmoothingPreserveShape(const auto &ps, const int iteration = 1) {
//    for (auto i = 0; i < iteration; ++i)
//       for (const auto &p : ps) SmoothingPreserveShape(p, ArithmeticWeightedSmoothingVector);
// };

inline Tddd ArithmeticWeightedSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX, std::function<Tddd(const networkPoint *)> position, const double power = 1.) {
  Tddd V = {0., 0, 0.}, qX;
  double Wtot = 0, norm;
  if (!isEdgePoint(p)) {
    for (const auto &q : p->getNeighbors()) {
      qX = position(q);
      V += std::pow(Norm(qX - current_pX), power - 1.) * (qX - current_pX);
      Wtot += 1.;
    }
    return V / Wtot;
  } else
    return V;
};

/* -------------------------------------------------------------------------- */

inline Tddd NeighborAverageSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX, const double power = 1.) {
  Tddd V = {0., 0, 0.};
  double Wtot = 0, norm;
  if (!isEdgePoint(p)) {
    for (const auto &q : p->getNeighbors()) {
      V += std::pow(Norm(q->X - current_pX), power - 1.) * (q->X - current_pX);
      Wtot += 1.;
    }
    return V / Wtot;
  } else
    return V;
};

inline Tddd TriangleCentroidAverageSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V = {0., 0, 0.};
  double Wtot = 0, W;
  if (!isEdgePoint(p)) {
    for (const auto &face : p->getFaces()) {
      auto [X0, X1, X2] = ToX(face->getPoints(p));
      X0 = current_pX;
      Wtot += 1.;
      V += W * Centroid(X0, X1, X2);
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

inline Tddd TriangleAreaAverageSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V = {0., 0, 0.};
  double Wtot = 0, W;
  if (!isEdgePoint(p)) {
    for (const auto &face : p->getFaces()) {
      auto [X0, X1, X2] = ToX(face->getPoints(p));
      X0 = current_pX;
      Wtot += (W = TriangleArea(X0, X1, X2));
      V += W * Centroid(X0, X1, X2);
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

inline Tddd IncenterAverageSmoothingVector(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V = {0., 0, 0.};
  double Wtot = 0;
  if (!isEdgePoint(p)) {
    for (const auto &face : p->getFaces()) {
      auto [X0, X1, X2] = ToX(face->getPoints(p));
      X0 = current_pX;
      Wtot += 1.;
      V += Incenter(X0, X1, X2);
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

inline Tddd EquilateralVertexAveragingVector(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V = {0., 0., 0.};
  if (!isEdgePoint(p)) {
    Tddd X = {0., 0., 0.}, Xmid, vertical;
    double Wtot = 0, W, height;
    for (const auto &f : p->getFaces()) {
      auto [X0, X1, X2] = ToX(f->getPoints(p));
      X0 = current_pX;
      Wtot += 1.;
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3) * 0.5;
      V += height * vertical + Xmid;
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

inline Tddd EquilateralVertexAveragingVector(const networkPoint *p, const std::array<double, 3> &current_pX, std::function<double(const networkFace *)> weight) {
  Tddd V = {0., 0., 0.};
  if (!isEdgePoint(p)) {
    Tddd X = {0., 0., 0.}, Xmid, vertical;
    double Wtot = 0, W, height;
    for (const auto &f : p->getFaces()) {
      auto [X0, X1, X2] = ToX(f->getPoints(p));
      X0 = current_pX;
      Wtot += weight(f);
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3) * 0.5;
      V += Wtot * (height * vertical + Xmid);
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

inline Tddd EquilateralVertexAveragingVector2(const networkPoint *p, const std::array<double, 3> &current_pX) {
  Tddd V = {0., 0., 0.};
  if (!isEdgePoint(p)) {
    Tddd X = {0., 0., 0.}, Xmid, vertical;
    double Wtot = 0, W, height;
    for (const auto &f : p->getFaces()) {
      auto [X0, X1, X2] = ToX(f->getPoints(p));
      X0 = current_pX;
      // W = std::log10(std::clamp(CircumradiusToInradius(X0, X1, X2), 1.0, 1E+10) - 1);
      W = log10_CircumradiusToInradius(X0, X1, X2);
      Wtot += W;
      Xmid = (X2 + X1) * 0.5;
      vertical = Normalize(Chop(X0 - Xmid, X2 - X1));
      height = Norm(X2 - X1) * std::sqrt(3) * 0.5;
      V += W * (height * vertical + Xmid);
    }
    return V / Wtot - current_pX;
  } else
    return V;
};

/* ------------------------------------------------------ */

inline void flipIf(networkPoint *p, const Tdd &limit_Dirichlet, const Tdd &limit_Neumann, bool force = false, int s_mean = 0) {
  try {
    auto [target_of_max_normal_diffD, acceptable_normal_change_by_flipD] = limit_Dirichlet;
    auto [target_of_max_normal_diffN, acceptable_normal_change_by_flipN] = limit_Neumann;
    // 2022/04/13こっちにBEMのmainから持ってきた
    // std::cout << "flipIf" << std::endl;
    for (auto &f : p->getFaces())
      f->setGeometricProperties();

    bool isfound = false;
    int count = 0;
    auto V = p->getLines();
    // ランダムにソート
    for (const auto &l : V) {
      auto [p0, p1] = l->getPoints();
      if (!l->CORNER) {
        if (force) {
          if (network_utility_detail::lineIsDirichletLocal(l)) {
            isfound = l->flipIfTopologicallyBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, s_mean == 0 ? 5 : s_mean);
            if (isfound)
              count++;
          } else {
            isfound = l->flipIfTopologicallyBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, s_mean == 0 ? 4 : s_mean);
            if (isfound)
              count++;
          }
        } else {
          if (network_utility_detail::lineIsDirichletLocal(l)) {
            //! 最小の変の数を３としている．もしこれを増やすと，柔軟に対応でいなくなる．特に角．
            isfound = l->flipIfBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, 5);
            if (isfound)
              count++;
          } else {
            isfound = l->flipIfBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, 4);
            if (isfound)
              count++;
          }
        }
      }
    }
    for (auto &f : p->getFaces())
      f->setGeometricProperties();

  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};

inline void flipIf(networkPoint *p, const Tdd &limit, bool force = false, int iteration = 0) { flipIf(p, limit, limit, force, iteration); };

inline void flipIf(Network &water, const Tdd &limit_Dirichlet, const Tdd &limit_Neumann, bool force = false, int iteration = 0) {
  try {
    auto [target_of_max_normal_diffD, acceptable_normal_change_by_flipD] = limit_Dirichlet;
    auto [target_of_max_normal_diffN, acceptable_normal_change_by_flipN] = limit_Neumann;
    // 2022/04/13こっちにBEMのmainから持ってきた
    std::cout << "flipIf" << std::endl;
    water.setGeometricPropertiesForce();
    double mean_length = Mean(extLength(water.getLines()));
    bool isfound = false, ismerged = false;
    int count = 0;
    auto V = ToVector(water.getLines());
    // ランダムにソート
    for (const auto &l : RandomSample(V)) {
      auto [p0, p1] = l->getPoints();
      if (!l->CORNER) {
        if (force && (iteration == 0 || count < iteration)) {
          if (network_utility_detail::lineIsDirichletLocal(l)) {
            isfound = l->flipIfTopologicallyBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD);
            if (isfound)
              count++;
          } else {
            isfound = l->flipIfTopologicallyBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN);
            if (isfound)
              count++;
          }
        } else {
          if (network_utility_detail::lineIsDirichletLocal(l)) {
            //! 最小の変の数を３としている．もしこれを増やすと，柔軟に対応でいなくなる．特に角．
            isfound = l->flipIfBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, 5);
            if (isfound)
              count++;
          } else {
            isfound = l->flipIfBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, 4);
            if (isfound)
              count++;
          }
        }
      }
    }
    if (count > 0) {
      std::cout << Green << "  " << count << " edges flipped." << colorReset << std::endl;
      water.setGeometricPropertiesForce();
    } else
      std::cout << Green << "  No edges flipped." << colorReset << std::endl;
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};

inline void flipIf(Network &water, const Tdd &limit, bool force = false, int iteration = 0) { flipIf(water, limit, limit, force, iteration); };

inline void flipIf(Network &water, double limit_angle = M_PI / 180., bool force = false, int iteration = 0) { flipIf(water, limit_angle, force, iteration); };

// thresholdの定義を「最短辺 / 他の2辺の平均長」に変更した改善版
inline void mergeIf(Network &water, double threshold = 0.1, int repetitions = 10) {
  // 無効な閾値の場合は何もしない
  if (threshold <= 0)
    return;

  std::cout << "mergeIf started with aspect ratio threshold < " << threshold << std::endl;

  for (int pass = 0; pass < repetitions; ++pass) {
    bool merged_in_this_pass = false;

    while (true) {
      bool merged_in_this_iteration = false;
      networkLine *line_to_merge = nullptr;

      auto faces = water.getBoundaryFaces();
      if (faces.empty())
        break;

      for (auto &f : faces) {
        auto lines = f->getLines();
        if (lines.size() != 3)
          continue;

        double lengths[3] = {lines[0]->length(), lines[1]->length(), lines[2]->length()};

        // 最短の辺の長さとインデックスを見つける
        double min_length = lengths[0];
        int min_index = 0;
        if (lengths[1] < min_length) {
          min_length = lengths[1];
          min_index = 1;
        }
        if (lengths[2] < min_length) {
          min_length = lengths[2];
          min_index = 2;
        }

        // 他の2辺の平均長を計算
        double others_mean_length = 0;
        for (int i = 0; i < 3; ++i)
          if (i != min_index)
            others_mean_length += lengths[i];

        others_mean_length /= 2.0;

        // 条件A: この三角形の形状が悪いか？ (新しいアスペクト比の定義)
        if (others_mean_length > 0 && ((min_length / others_mean_length < (1. - threshold)) || (min_length / others_mean_length > (1. + threshold)))) {
          // 【修正】変数名を min_index に修正
          networkLine *candidate_line = lines[min_index];

          // 条件B: この候補の辺は、マージしても安全か？
          if (candidate_line && !candidate_line->CORNER && candidate_line->isMergeable()) {
            line_to_merge = candidate_line;
            break;
          }
        }
      }

      if (line_to_merge) {
        std::cout << Red << "  Merging line (length: " << line_to_merge->length() << ")" << colorReset << std::endl;
        line_to_merge->Collapse();
        delete line_to_merge; // マージ後にラインを削除
        merged_in_this_iteration = true;
        merged_in_this_pass = true;
      }

      if (!merged_in_this_iteration) {
        break;
      }
    }

    if (!merged_in_this_pass) {
      std::cout << Green << "  No more lines to merge in this pass." << colorReset << std::endl;
      break;
    }
  }

  water.setGeometricPropertiesForce();
  std::cout << "mergeIf finished." << std::endl;
};

/* ------------------------------------------------------ */

// void LaplacianSmoothingIfFlat(netPp p) {
//    try {
//       if (!isEdgePoint(p) && isFlat(p) && isInConvexPolygon(p)) { /*端の点はsmoothingしない*/
//          auto ps = p->getNeighbors();
//          if (ps.size() > 2)  // 2点の場合は2点の中点に動いてしまうので，実行しない
//             p->setX(Mean(extractX(ps)));
//       }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorReset << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };
//////////////////////////////////////////////////////////
// void LaplacianSmoothingIfFlat(V_netPp ps /*copy*/) {
//    for (const auto &p : ps)
//       LaplacianSmoothingIfFlat(p);
// };
///////////////////////////////////////////////////////////
// void LaplacianSmoothingIfFlat_ExceptIntX(V_netPp ps /*copy*/) {
//    try {
//       std::shuffle(std::begin(ps), std::end(ps), std::default_random_engine());
//       for (const auto &p : ps)
//          if (!isEdgePoint(p) && isFlat(p) && isInConvexPolygon(p)) {  // 端の点はsmoothingしない
//             auto countintx = 0;
//             for (const auto &l : p->getLines())
//                if (l->isIntxn())
//                   countintx++;
//             if (!(countintx > 1)) {
//                auto tmp = Mean(extractX(p->getNeighbors()));
//                p->setX(tmp);
//             }
//          }
//    } catch (std::exception &e) {
//       std::cerr << e.what() << colorReset << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//    };
// };
///////////////////////////////////////////////////////////
inline void LaplacianSmoothingIfOnStraightLine(V_netPp ps) {
  try {
    std::shuffle(std::begin(ps), std::end(ps), std::default_random_engine());
    for (auto p : ps)
      if (!isEdgePoint(p) && !isFlat(p)) // 端の点はsmoothingしない
      {
        std::map<netPp, V_d> mapP_vec;
        for (auto q : p->getNeighbors())
          mapP_vec[q] = ToVector(q->getXtuple() - ToX(p)) / Norm(q->getXtuple() - ToX(p));

        netPp p0, p1;
        bool isStraight = false;
        auto Ps = TakeFirst(mapP_vec);
        double minarea = 1E-2;
        for (auto i = 0; i < Ps.size() - 1; i++)
          for (auto j = i + 1; j < Ps.size(); j++) {
            // (Ps[i]) <--------- (p) <-------- (Ps[j]) : angle = 0
            double abs_angle = std::abs(VectorAngle(Ps[i]->getXtuple() - ToX(p), ToX(p) - Ps[j]->getXtuple()));
            if (abs_angle < minarea) {
              p0 = Ps[i];
              p1 = Ps[j];
              isStraight = true;
              minarea = abs_angle;
            }
          }
        if (isStraight)
          p->setX(ToVector(ToX(p0) + ToX(p1)) / 2.);
      }
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};

inline double Distance(const Tddd &X, const networkPoint *const p) { return Norm(ToX(p) - X); };
inline double Distance(const networkPoint *const p, const Tddd &X) { return Norm(ToX(p) - X); };
inline double Distance(const networkPoint *const p, const networkPoint *const q) { return Norm(p->X - q->X); };
inline double Distance(const networkLine *const l, const networkPoint *const p) { return Norm(ToX(p) - l->getXtuple()); };
inline double Distance(const networkPoint *const p, const networkLine *const l) { return Norm(ToX(p) - l->getXtuple()); };
inline V_d Distance(const networkPoint *const p, const std::vector<networkPoint *> &ps) {
  V_d ret(ps.size());
  int i = 0;
  for (const auto &q : ps)
    ret[i++] = Distance(q, p);
  return ret;
};
inline V_d Distance(const std::vector<networkPoint *> &ps, const networkPoint *p) { return Distance(p, ps); };
/*
sortは，
comp(a,b)->trueなら
comp(b,a)->false
でないといけない
*/
inline void sortByLength(V_netLp &lines) {
  std::sort(lines.begin(), lines.end(), [](const auto &v, const auto &w) {
    if (v->length() - w->length() <= 1E-10)
      return v < w;
    else
      return v->length() < w->length();
  });
};

inline void reverseSortByLength(V_netLp &lines) {
  std::sort(lines.begin(), lines.end(), [](const auto &w, const auto &v) {
    if (v->length() - w->length() <= 1E-10)
      return v < w;
    else
      return v->length() < w->length();
  });
};

// void reverseSortByLength(const std::unordered_set<networkLine *> &uo_lines) {
//    V_netLp lines(uo_lines.begin(), uo_lines.end());
//    std::sort(lines.begin(), lines.end(),
//              [](const auto &w, const auto &v) {
//                 if (v->length() - w->length() <= 1E-10)
//                    return v < w;
//                 else
//                    return v->length() < w->length();
//              });
// };

inline void sortByDistance(V_netPp &points, const netPp origin) {
  std::sort(points.begin(), points.end(), [&origin](const netPp &v, const netPp &w) {
    if (Norm(ToX(v) - ToX(w)) <= 1E-13)
      return v < w;
    else
      return Norm(ToX(v) - ToX(origin)) < Norm(ToX(w) - ToX(origin));
  });
};

inline void sortByDistance(V_netPp &points, const Tddd &origin) {
  std::sort(points.begin(), points.end(), [&](const netPp &v, const netPp &w) {
    if (Norm(ToX(v) - ToX(w)) <= 1E-13)
      return v < w;
    else
      return Norm(ToX(v) - origin) < Norm(ToX(w) - origin);
  });
};

// void Merge(const Network *net, const std::function<bool(const networkLine *)> &condition) {
//    bool found = true;
//    int count = 0;
//    while (found && count++ < 5) {
//       found = false;
//       for (const auto &l : net->Lines)
//          if (!found && condition(l) && l->Collapse()) {
//             found = true;
//             break;
//          }
//    }
// };
inline void Merge(const std::unordered_set<networkLine *> &uo_lines, const std::function<bool(const networkLine *)> &func) {
  for (const auto &l : uo_lines)
    if (func(l)) {
      l->Collapse();
      return;
    }
  // std::unordered_set<std::shared_ptr<networkLine>> shared_pointers;
  // std::transform(uo_lines.begin(), uo_lines.end(), std::inserter(shared_pointers, shared_pointers.begin()),
  //                [](networkLine *raw_ptr) { return std::make_shared<networkLine>(*raw_ptr); });
  // std::vector<std::shared_ptr<networkLine>> shared_pointers_copy(shared_pointers.begin(), shared_pointers.end());
  // for (const auto &l : shared_pointers_copy) {
  //    if (func(l.get())) {
  //       l->Collapse();
  //       shared_pointers.erase(l);
  //    }
  // }
};

inline void Divide(std::unordered_set<networkLine *> uo_lines, const std::function<bool(networkLine *)> &func) {
  for (const auto &l : uo_lines)
    if (func(l))
      l->Split();
};

inline void Divide(const std::unordered_set<networkLine *> &uo_lines, const double lim_len) {
  std::cout << Red << "|";
  for (const auto &l : uo_lines)
    if (l->length() >= lim_len) {
      std::cout << Green << "|";
      l->Split();
    }
  std::cout << Blue << "|" << colorReset;
};

//////////////////////////////////////////////////////////////
//////////////////////////// 外部関数 /////////////////////////
//////////////////////////////////////////////////////////////
inline std::vector<netPp> getPointsWithoutLines(const std::vector<netPp> &obj) {
  std::vector<netPp> ret;
  for (const auto &o : obj)
    if (o->getLines().empty())
      ret.emplace_back(o);
  return ret;
};

inline std::vector<netPp> getPointsWithoutFaces(const std::vector<netPp> &obj) {
  std::vector<netPp> ret;
  V_netLp ls;
  bool facefound;
  for (const auto &o : obj) {
    ls = o->getLines();
    if (ls.empty())
      ret.emplace_back(o);
    else {
      facefound = false;
      for (const auto &l : ls)
        if (!l->getFaces().empty()) {
          facefound = true;
          break;
        }
      if (!facefound)
        ret.emplace_back(o);
    }
  }
  return ret;
};

inline V_d areas(const V_netFp &fs) {
  V_d ret(fs.size());
  for (auto i = 0; i < fs.size(); i++)
    ret[i] = fs[i]->area;
  return ret;
};
/* ------------------------------------------------------ */
inline std::vector<std::string> getNames(const V_netFp &fs) {
  std::vector<std::string> ret;
  for (const auto &f : fs)
    ret.emplace_back(f->getNetwork()->getName());
  return ret;
};

inline std::vector<std::string> getNames(const V_netPp &ps) {
  std::vector<std::string> ret;
  for (const auto &p : ps)
    ret.emplace_back(p->getNetwork()->getName());
  return ret;
};

inline std::vector<std::string> getNames(const V_Netp &nets) {
  std::vector<std::string> ret;
  for (const auto &n : nets)
    ret.emplace_back(n->getName());
  return ret;
};

inline V_netLp takeIntxn(const V_netLp &ls) {
  V_netLp ret({});
  for (const auto &l : ls)
    if (l->isIntxn())
      ret.emplace_back(l);
  return ret;
};

template <typename T> std::vector<T *> takeNetwork(const std::vector<T *> &ps, const std::vector<Network *> net) {
  std::vector<T *> ret = {};
  for (const auto &p : ps)
    if (MemberQ(net, p->getNetwork()))
      ret.emplace_back(p);
  return ret;
};

////////////////////////////////////////////////////////////////
////////////////////////// dislay //////////////////////////////
////////////////////////////////////////////////////////////////

template <typename T> void displayNames(const std::vector<T *> &ps) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::map<Network *, int> storages;
  for (const auto &p : ps)
    if (storages.find(p->getNetwork()) != storages.end())
      storages[p->getNetwork()] += 1;
    else
      storages[p->getNetwork()] = 1;
  for (auto [n, i] : storages)
    std::cout << std::setw(8) << n->getName() << "  " << std::setw(8) << i << std::endl;
  if (storages.empty())
    std::cout << Red << "not stored !?" << colorReset << std::endl;
};

template <typename T> void displayStorages(const std::vector<T *> &ps) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::map<Network *, int> storages;
  for (const auto &p : ps)
    if (storages.find(p->getStorage()) != storages.end())
      storages[p->getStorage()] += 1;
    else
      storages[p->getStorage()] = 1;
  for (auto [n, i] : storages)
    std::cout << std::setw(8) << n->getName() << "  " << std::setw(8) << i << std::endl;
  if (storages.empty())
    std::cout << Red << "not stored !?" << colorReset << std::endl;
};

inline void display(const std::unordered_set<networkLine *> &ls) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::map<size_t, int> data;
  int intx = 0;
  for (const auto &l : ls) {
    auto s = l->getFaces().size();
    if (l->isIntxn())
      intx++;
    if (data.find(s) != data.end())
      data[s] += 1;
    else
      data[s] = 1;
  }
  for (auto [n, i] : data)
    std::cout << std::setw(8) << "Faces size = " << n << "  " << std::setw(8) << i << std::endl;

  std::cout << std::setw(8) << "intx = " << intx << "  " << std::setw(8) << intx << std::endl;
};

inline void display(Network *net) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  V_i num(15, 0), pointsLines(15, 0), linesFaces(15, 0);
  ;
  for (const auto &p : net->getPoints())
    for (auto i = 0; i < 15; i++)
      if (p->getLines().size() == (size_t)i)
        pointsLines[i]++;
  // for (const auto &f : net->getFaces())
  // {
  // 	auto f_Points = f->getPoints();
  // 	for (auto i = 0; i < 15; i++)
  // 	{
  // 		if (f_Points.size() == (size_t)i)
  // 			facesPoints[i]++;
  // 		// if (f->getLines().size() == (size_t)i)
  // 		// 	facesLines[i]++;
  // 	}
  // }
  for (const auto &l : net->getLines()) {
    auto fs = l->getFaces();
    for (auto i = 0; i < 15; i++)
      if (fs.size() == (size_t)i)
        linesFaces[i]++;
  }
  std::cout << "Network name : " << net->getName() << std::endl;
  std::cout << "--------------------Names---------------------------" << std::endl;
  std::cout << Blue << "Points : " << colorReset;
  // displayNames(net->getPoints());
  std::cout << Magenta << " Faces : " << colorReset;
  // displayNames(net->getFaces());
  std::cout << "------------------Storages--------------------------" << std::endl;
  std::cout << Blue << "Points : " << colorReset;
  // displayStorages(net->getPoints());
  std::cout << Magenta << " Faces : " << colorReset;
  // displayStorages(net->getFaces());
  std::cout << "--------------------Lines--------------------------" << std::endl;
  display(net->getLines());
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << Magenta << "Points.size() : " << net->getPoints().size() << std::endl;
  std::cout << Magenta << " Faces.size() : " << net->getFaces().size() << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  int i = 0;
  for (auto &n : num) {
    n = i++;
  };
  std::cout << Red << "                  ";
  for (const auto &n : num) {
    std::cout << std::setw(6) << n;
  };
  std::cout << colorReset << std::endl;

  std::cout << magenta << "Lines of points : ";
  for (const auto &n : pointsLines) {
    std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
  }
  std::cout << colorReset << std::endl;

  // std::cout << magenta << "Points of faces : ";
  // for (const auto &n : facesPoints)
  // {
  // 	std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
  // };
  // std::cout << colorReset << std::endl;

  // std::cout << magenta << " Lines of faces : ";
  // for (const auto &n : facesLines)
  // {
  // 	std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
  // };
  std::cout << colorReset << std::endl;

  std::cout << magenta << " Faces of Lines : ";
  for (const auto &n : linesFaces) {
    std::cout << (n == 0 ? magenta : Magenta) << std::setw(6) << n;
  };
  std::cout << colorReset << std::endl;

  V_i connection(3);
  // net->setLinesStatus(true);
  for (const auto &p : net->getPoints())
    for (const auto &line : p->getLines()) {
      if ((line->getFaces()).size() == 0)
        connection[0]++;
      else if ((line->getFaces()).size() == 1)
        connection[1]++;
      else if ((line->getFaces()).size() == 2)
        connection[2]++;
    }
  std::cout << magenta << " Connection of faces : " << connection;
  if (connection[0] == 0 && connection[1] == 0 && connection[2] != 0)
    std::cout << Blue << " <- Network is a closed surface";
  std::cout << colorReset << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
};

//* ------------------------------------------------------ */
//*                     粒子法のためのユーティリティ            */
//* ------------------------------------------------------ */
inline std::vector<Tddd> fullparticlize(const Tddd &X0, const Tddd &X1, const Tddd &X2, const double dx) {
  /*
  |-o--o--o--o-| n = 4
  x,yはパラメタ
  */
  interpolationTriangleLinearTuple intp({X0, X1, X2});
  auto y_list = [&intp, &dx](const double x) {
    if (1. - x < 1E-10)
      return V_d{0.};
    int n = std::round(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
    if (n > 0)
      return Subdivide(0., 1. - x, n);
    else
      return V_d{};
  };
  /* ------------------------------------------------------ */
  Tddd v = intp(0., 1.) - intp(0., 0.);
  Tddd u = intp(1., 0.) - intp(0., 0.);
  double height = Norm(u - Dot(Normalize(v), u) * Normalize(v));
  int n = std::round(height / dx);
  /* ------------------------------------------------------ */
  std::vector<Tddd> ret = {};
  if (n > 0) {
    for (const auto &x : Subdivide(0., 1., n))
      for (const auto &y : y_list(x))
        ret.emplace_back(intp(x, y));
  }
  return ret;
};
inline std::vector<Tddd> fullparticlize(networkFace const *f, const double dx) {
  // auto ps = f->getPoints();
  // return fullparticlize(ps[0]->getXtuple(), ps[1]->getXtuple(), ps[2]->getXtuple(), dx);

  auto [p0, p1, p2] = f->getPoints();
  return fullparticlize(ToX(p0), ToX(p1), ToX(p2), dx);
};

/*
        % particlizeは座標とパラメタを返す
        % 面をdx間隔で粒子化し，各点の情報をベクトルで返す
*/

// b% -------------- networkFaceを使ったparticlize ------------- */

inline std::vector<std::tuple<Tddd, Tdd>> particlize(const networkFace *const f, const double dx) { return particlize(f->getXVertices(), dx); };
inline std::vector<std::tuple<Tddd, Tdd>> particlize(const V_netFp &faces, const double dx) {
  std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
  for (const auto &f : faces) {
    tmp = particlize(f, dx);
    ret.insert(ret.end(), tmp.begin(), tmp.end());
  }
  return ret;
};
inline std::vector<std::tuple<Tddd, Tdd>> particlize(const std::unordered_set<networkFace *> &faces, const double dx) {
  std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
  for (const auto &f : faces) {
    tmp = particlize(f, dx);
    ret.insert(ret.end(), tmp.begin(), tmp.end());
  }
  return ret;
};

/* ------------------------ 深さあり ------------------------ */

inline std::vector<std::tuple<Tddd, Tdd>> particlize(const networkFace *const f, const double dx, const V_d &depth_list /*法線方向にdx*depthだけ動かす{-1,0,1,2,3,..}など*/) {
  std::vector<std::tuple<Tddd, Tdd>> ret, tmp;
  double alpha;
  Tddd X0, X1, X2;
  for (const auto &d : depth_list /*double実際の長さ*/) {
    auto [p0, p1, p2] = f->getPoints();
    X0 = ToX(p0) + d / Dot(p0->getNormalTuple(), f->normal) * p0->getNormalTuple();
    X1 = ToX(p1) + d / Dot(p1->getNormalTuple(), f->normal) * p1->getNormalTuple();
    X2 = ToX(p2) + d / Dot(p2->getNormalTuple(), f->normal) * p2->getNormalTuple();
    tmp = particlize({X0, X1, X2}, dx);
    ret.insert(ret.end(), tmp.begin(), tmp.end());
  }
  return ret;
};

inline Tddd particlize(const networkFace *f, const Tdd &t0t1, const double d) {
  auto [p0, p1, p2] = f->getPoints();
  return Dot(Tddd{std::get<0>(t0t1), std::get<1>(t0t1), 1. - std::get<0>(t0t1) - std::get<1>(t0t1)}, T3Tddd{ToX(p0) + d / Dot(p0->getNormalTuple(), f->normal) * p0->getNormalTuple(), ToX(p1) + d / Dot(p1->getNormalTuple(), f->normal) * p1->getNormalTuple(), ToX(p2) + d / Dot(p2->getNormalTuple(), f->normal) * p2->getNormalTuple()});
};

// b% ------------------------------------------------------ */
// Tddd oppositeX(const std::tuple<networkFace * /*補間に使った三角形の頂点*/,
//                                 T_PPP /*補間に使った三角形の頂点*/,
//                                 Tdd /*パラメタt0,t1*/,
//                                 double /*深さ方向距離*/,
//                                 double /*粒子間隔*/> &particlize_info) {
//    // 粒子化された点の面にたいする反対側の位置を返す．
//    // 単純にp+2*Dot(f->p0 - p,n)*nでは，角の面において重なりが生じてしまう．
//    // パラメタを使って計算すれば重なりは生じない．
//    auto [f, p0p1p2, t0t1, d, dx] = particlize_info;
//    if (f) {
//       auto [p0, p1, p2] = p0p1p2;  // f->getPoints();
//       return Dot(Tddd{std::get<0>(t0t1), std::get<1>(t0t1), 1. - std::get<0>(t0t1) - std::get<1>(t0t1)},
//                  T3Tddd{ToX(p0) + d / Dot(p0->getNormalTuple(), f->normal) * p0->getNormalTuple(),
//                         ToX(p1) + d / Dot(p1->getNormalTuple(), f->normal) * p1->getNormalTuple(),
//                         ToX(p2) + d / Dot(p2->getNormalTuple(), f->normal) * p2->getNormalTuple()});
//    } else {
//       std::stringstream ss;
//       // ss << "particlize_info = " << f << ", " << t0 << ", " << t1 << ", " << d << ", " << dx << std::endl;
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
//    };
// };

// Tddd oppositeX(const networkPoint *p) {
//    return oppositeX(p->particlize_info);
// };
/* ------------------------------------------------------ */

inline std::vector<networkPoint *> operator+(std::vector<networkPoint *> A /*copy and return*/, const std::vector<networkPoint *> &B) {
  A.insert(A.end(), B.begin(), B.end());
  return A;
};
inline std::vector<networkPoint *> &operator+=(std::vector<networkPoint *> &v, const std::vector<networkPoint *> &w) { return (v = v + w); };
#endif
