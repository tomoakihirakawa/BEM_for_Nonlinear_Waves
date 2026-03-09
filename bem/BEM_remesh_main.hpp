#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BEM_collision.hpp"

struct SubsurfaceAltitudeCheckResult {
  std::size_t checked_lines = 0;
  std::size_t angle_lines = 0;
  double worst_altitude_rel = 1E+100;
  double worst_angle_deg = 0.;
  networkLine* worst_line = nullptr;
  networkFace* worst_face = nullptr;
  int worst_subface_index = -1;
};

struct TinyFaceCheckResult {
  std::size_t checked_faces = 0;
  double worst_area_ratio = 1E+100;
  double worst_area = 0.0;
  double worst_local_mean_area = 0.0;
  networkFace* worst_face = nullptr;
  networkLine* worst_shortest_line = nullptr;
  int worst_subface_index = -1;
};

using FaceTriangle = std::array<Tddd, 3>;

struct FaceTriangleWithBase {
  FaceTriangle tri;
  int subface_index = -1;
};

struct FaceAltitudeDetail {
  double altitude_rel = 1E+100;
  int subface_index = -1;
};

inline std::vector<FaceTriangle> qualityTriangles(networkFace* face) {
  std::vector<FaceTriangle> tris;
  if (!face)
    return tris;

  auto [p0, p1, p2] = face->getPoints();
  if (!p0 || !p1 || !p2)
    return tris;

  if (!face->isTrueQuadraticElement) {
    tris.push_back({p0->X, p1->X, p2->X});
    return tris;
  }

  auto [l0, l1, l2] = face->getLines();
  if (!l0 || !l1 || !l2) {
    tris.push_back({p0->X, p1->X, p2->X});
    return tris;
  }

  tris.push_back({p0->X, l0->X_mid, l2->X_mid});
  tris.push_back({l0->X_mid, p1->X, l1->X_mid});
  tris.push_back({l2->X_mid, l1->X_mid, p2->X});
  tris.push_back({l0->X_mid, l1->X_mid, l2->X_mid});
  return tris;
}

inline std::vector<FaceTriangleWithBase> qualityTrianglesAlongSharedLine(networkFace* face, networkLine* shared_line) {
  std::vector<FaceTriangleWithBase> tris;
  if (!face || !shared_line)
    return tris;

  auto [p0, p1, p2] = face->getPoints();
  if (!p0 || !p1 || !p2)
    return tris;

  if (!face->isTrueQuadraticElement)
    return tris;

  auto [l0, l1, l2] = face->getLines();
  if (!l0 || !l1 || !l2) {
    tris.push_back({{p0->X, p1->X, p2->X}, -1});
    return tris;
  }

  if (shared_line == l0) {
    tris.push_back({{p0->X, l0->X_mid, l2->X_mid}, 0});
    tris.push_back({{l0->X_mid, p1->X, l1->X_mid}, 1});
  } else if (shared_line == l1) {
    tris.push_back({{p1->X, l1->X_mid, l0->X_mid}, 0});
    tris.push_back({{l1->X_mid, p2->X, l2->X_mid}, 1});
  } else if (shared_line == l2) {
    tris.push_back({{p2->X, l2->X_mid, l1->X_mid}, 0});
    tris.push_back({{l2->X_mid, p0->X, l0->X_mid}, 1});
  } else {
    tris.push_back({{p0->X, p1->X, p2->X}, -1});
  }
  return tris;
}

inline double triangleArea(const FaceTriangle& tri) {
  return TriangleArea(tri[0], tri[1], tri[2]);
}

inline double triangleMinimumInteriorAngleDeg(const FaceTriangle& tri) {
  const double l01 = Norm(tri[0] - tri[1]);
  const double l12 = Norm(tri[1] - tri[2]);
  const double l20 = Norm(tri[2] - tri[0]);
  if (l01 <= 1e-20 || l12 <= 1e-20 || l20 <= 1e-20)
    return 0.0;
  auto clamp_cos = [](double x) { return std::max(-1.0, std::min(1.0, x)); };
  const double ang0 = std::acos(clamp_cos(Dot((tri[1] - tri[0]) / l01, (tri[2] - tri[0]) / l20))) * 180.0 / M_PI;
  const double ang1 = std::acos(clamp_cos(Dot((tri[0] - tri[1]) / l01, (tri[2] - tri[1]) / l12))) * 180.0 / M_PI;
  const double ang2 = 180.0 - ang0 - ang1;
  return std::min({ang0, ang1, ang2});
}

inline double triangleAspectRatio(const FaceTriangle& tri) {
  const double l01 = Norm(tri[0] - tri[1]);
  const double l12 = Norm(tri[1] - tri[2]);
  const double l20 = Norm(tri[2] - tri[0]);
  const double max_edge = std::max({l01, l12, l20});
  if (!(max_edge > 1e-20))
    return 1E+100;
  const double area = triangleArea(tri);
  const double altitude = (area > 0.0) ? 2.0 * area / max_edge : 0.0;
  return (altitude > 1e-20) ? max_edge / altitude : 1E+100;
}

inline double triangleAltitudeRelativeToBase(const FaceTriangle& tri) {
  const double base = Norm(tri[1] - tri[0]);
  if (!(base > 0.0) || !std::isfinite(base))
    return 0.0;
  const double area = triangleArea(tri);
  if (!(area > 0.0) || !std::isfinite(area))
    return 0.0;
  return (2.0 * area / base) / base;
}

inline double boundaryFaceArea(networkFace* face) {
  if (!face)
    return 0.0;
  auto [p0, p1, p2] = face->getPoints();
  if (!p0 || !p1 || !p2)
    return 0.0;
  return TriangleArea(p0->X, p1->X, p2->X);
}

inline double localMeanFaceArea(networkFace* face) {
  if (!face)
    return 0.0;
  std::unordered_set<networkFace*> adjacent_faces;
  for (auto* p : face->getPoints()) {
    if (!p)
      continue;
    for (auto* f : p->getBoundaryFaces())
      if (f && f != face)
        adjacent_faces.insert(f);
  }
  if (adjacent_faces.empty())
    return 0.0;

  double sum = 0.0;
  std::size_t count = 0;
  for (auto* f : adjacent_faces) {
    const double area = boundaryFaceArea(f);
    if (!(area > 0.0) || !std::isfinite(area))
      continue;
    sum += area;
    ++count;
  }
  return (count > 0) ? sum / static_cast<double>(count) : 0.0;
}

inline double localMeanFaceAreaForPoints(const std::vector<networkPoint*>& points,
                                         const std::unordered_set<networkFace*>& exclude = {}) {
  std::unordered_set<networkFace*> adjacent_faces;
  for (auto* p : points) {
    if (!p)
      continue;
    for (auto* f : p->getBoundaryFaces())
      if (f && !exclude.contains(f))
        adjacent_faces.insert(f);
  }
  if (adjacent_faces.empty())
    return 0.0;

  double sum = 0.0;
  std::size_t count = 0;
  for (auto* f : adjacent_faces) {
    const double area = boundaryFaceArea(f);
    if (!(area > 0.0) || !std::isfinite(area))
      continue;
    sum += area;
    ++count;
  }
  return (count > 0) ? sum / static_cast<double>(count) : 0.0;
}

inline double localLineLength(const networkLine* l) {
  if (!l)
    return 0.0;
  std::unordered_set<networkPoint*> adjacent_points;
  for (const auto& f : l->getBoundaryFaces()) {
    auto points = f->getPoints();
    adjacent_points.insert(points.begin(), points.end());
  }
  std::unordered_set<networkLine*> adjacent_lines;
  for (const auto& p : adjacent_points)
    for (const auto& L : p->getBoundaryLines())
      if (L != l)
        adjacent_lines.insert(L);
  if (adjacent_lines.empty())
    return 0.0;
  double ret = 0.0;
  for (const auto& L : adjacent_lines)
    ret += L->length();
  return ret / adjacent_lines.size();
}

inline double minimumInteriorAngleDeg(networkFace* face) {
  if (!face)
    return 1E+100;
  const auto tris = qualityTriangles(face);
  if (tris.empty())
    return 1E+100;
  double worst = 1E+100;
  for (const auto& tri : tris)
    worst = std::min(worst, triangleMinimumInteriorAngleDeg(tri));
  return worst;
}

inline double faceAspectRatio(networkFace* face) {
  if (!face)
    return 1E+100;
  const auto tris = qualityTriangles(face);
  if (tris.empty())
    return 1E+100;
  double worst = 0.0;
  for (const auto& tri : tris)
    worst = std::max(worst, triangleAspectRatio(tri));
  return (worst > 0.0) ? worst : 1E+100;
}

inline bool isPostTypeCollapseTargetLine(const networkLine* l) {
  return l && !l->CORNER;
}

struct CornerConnectedNeumannCollapseResult {
  bool changed = false;
  int collapsed = 0;
};

inline CornerConnectedNeumannCollapseResult collapseCornerConnectedNeumannLinesAfterBoundaryTypes(Network& water,
                                                                                                  const int time_step) {
  CornerConnectedNeumannCollapseResult result;

  auto log_candidate = [&](const char* event, networkLine* l, double local_mean_len, double alt0, double alt1,
                           double area_ratio0, double area_ratio1, double min_angle_deg0, double min_angle_deg1) {
    if (!l)
      return;
    const auto faces = l->getBoundaryFaces();
    if (faces.size() != 2 || !faces[0] || !faces[1])
      return;
    const auto pts = l->getPoints();
    std::cout << Magenta << "[corner_neumann_post_type] event=" << event
              << " time_step=" << time_step
              << " len=" << l->length()
              << " local_mean_len=" << local_mean_len
              << " alt_ratio0=" << ((local_mean_len > 0.0) ? alt0 / local_mean_len : 1E+100)
              << " alt_ratio1=" << ((local_mean_len > 0.0) ? alt1 / local_mean_len : 1E+100)
              << " area_ratio0=" << area_ratio0
              << " area_ratio1=" << area_ratio1
              << " min_angle_deg0=" << min_angle_deg0
              << " min_angle_deg1=" << min_angle_deg1
              << " line_flags={D:" << l->Dirichlet << ",N:" << l->Neumann << ",C:" << l->CORNER << "}"
              << " endpoint_corner={" << (pts[0] ? pts[0]->CORNER : false) << "," << (pts[1] ? pts[1]->CORNER : false) << "}"
              << " x0=" << (pts[0] ? pts[0]->X : Tddd{0., 0., 0.})
              << " x1=" << (pts[1] ? pts[1]->X : Tddd{0., 0., 0.})
              << colorReset << std::endl;
  };

  for (int iter = 0; iter < 8; ++iter) {
    struct Candidate {
      networkLine* line;
      double len;
      double local_mean_len;
      double alt0;
      double alt1;
      double area_ratio0;
      double area_ratio1;
      double min_angle_deg0;
      double min_angle_deg1;
    };

    std::vector<Candidate> candidates;
    for (auto* l : water.getBoundaryLines()) {
      if (!isPostTypeCollapseTargetLine(l))
        continue;
      const auto faces = l->getBoundaryFaces();
      if (faces.size() != 2 || !faces[0] || !faces[1])
        continue;
      const double local_mean_len = localLineLength(l);
      if (!(local_mean_len > 0.0) || !std::isfinite(local_mean_len))
        continue;
      const double area0 = boundaryFaceArea(faces[0]);
      const double area1 = boundaryFaceArea(faces[1]);
      const double mean_area0 = localMeanFaceArea(faces[0]);
      const double mean_area1 = localMeanFaceArea(faces[1]);
      const double max_edge0 = std::max({Norm(faces[0]->getPoints()[0]->X - faces[0]->getPoints()[1]->X),
                                         Norm(faces[0]->getPoints()[1]->X - faces[0]->getPoints()[2]->X),
                                         Norm(faces[0]->getPoints()[2]->X - faces[0]->getPoints()[0]->X)});
      const double max_edge1 = std::max({Norm(faces[1]->getPoints()[0]->X - faces[1]->getPoints()[1]->X),
                                         Norm(faces[1]->getPoints()[1]->X - faces[1]->getPoints()[2]->X),
                                         Norm(faces[1]->getPoints()[2]->X - faces[1]->getPoints()[0]->X)});
      const double alt0 = (max_edge0 > 1e-20) ? 2.0 * area0 / max_edge0 : 0.0;
      const double alt1 = (max_edge1 > 1e-20) ? 2.0 * area1 / max_edge1 : 0.0;
      const double area_ratio0 = (mean_area0 > 0.0) ? area0 / mean_area0 : 1E+100;
      const double area_ratio1 = (mean_area1 > 0.0) ? area1 / mean_area1 : 1E+100;
      const double min_angle_deg0 = minimumInteriorAngleDeg(faces[0]);
      const double min_angle_deg1 = minimumInteriorAngleDeg(faces[1]);
      const double alt_ratio_min = std::min(alt0, alt1) / local_mean_len;
      if (alt_ratio_min >= 0.2)
        continue;
      if (std::min(area_ratio0, area_ratio1) >= 0.2)
        continue;
      if (std::min(min_angle_deg0, min_angle_deg1) >= 20.0)
        continue;
      candidates.push_back({l, l->length(), local_mean_len, alt0, alt1, area_ratio0, area_ratio1, min_angle_deg0, min_angle_deg1});
    }

    if (candidates.empty())
      break;

    std::ranges::sort(candidates, [](const Candidate& a, const Candidate& b) { return a.len < b.len; });

    bool collapsed_this_iter = false;
    for (const auto& c : candidates) {
      log_candidate("attempt", c.line, c.local_mean_len, c.alt0, c.alt1, c.area_ratio0, c.area_ratio1, c.min_angle_deg0, c.min_angle_deg1);
      if (c.line->Collapse()) {
        result.changed = true;
        ++result.collapsed;
        collapsed_this_iter = true;
        std::cout << Magenta << "[corner_neumann_post_type] event=success time_step=" << time_step
                  << " collapsed_len=" << c.len << colorReset << std::endl;
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        break;
      }
      log_candidate("failed", c.line, c.local_mean_len, c.alt0, c.alt1, c.area_ratio0, c.area_ratio1, c.min_angle_deg0, c.min_angle_deg1);
    }

    if (!collapsed_this_iter)
      break;
  }

  return result;
}

inline Tddd collapseTargetPointOnLine(const networkLine* line) {
  if (!line)
    return {0., 0., 0.};
  auto [pA, pB] = line->getPoints();
  if (!pA || !pB)
    return {0., 0., 0.};
  if (pA->CORNER && !pB->CORNER)
    return pA->X;
  if (pB->CORNER && !pA->CORNER)
    return pB->X;
  return 0.5 * (pA->X + pB->X);
}

inline double signedPatchVolumeAroundLineAfterCollapse(const networkLine* line, const bool after_collapse) {
  if (!line)
    return 0.0;
  auto [pA, pB] = line->getPoints();
  if (!pA || !pB)
    return 0.0;

  const Tddd target = collapseTargetPointOnLine(line);
  const Tddd x_ref = target;
  std::unordered_set<networkFace*> patch_faces;
  for (auto* p : {pA, pB}) {
    for (auto* f : p->getBoundaryFaces())
      if (f)
        patch_faces.insert(f);
  }

  double volume = 0.0;
  for (auto* face : patch_faces) {
    if (!face)
      continue;
    if (after_collapse) {
      const auto lines = face->getLines();
      if (std::ranges::find(lines, const_cast<networkLine*>(line)) != lines.end())
        continue;
    }

    auto pts = face->getPoints();
    std::array<Tddd, 3> xyz = {pts[0]->X, pts[1]->X, pts[2]->X};
    if (after_collapse) {
      for (int i = 0; i < 3; ++i)
        if (pts[i] == pA || pts[i] == pB)
          xyz[i] = target;
      if (Norm(Cross(xyz[1] - xyz[0], xyz[2] - xyz[0])) < 1E-20)
        continue;
    }
    volume += Dot(xyz[0] - x_ref, Cross(xyz[1] - x_ref, xyz[2] - x_ref)) / 6.0;
  }
  return volume;
}

inline double worstTinyFaceAreaRatioOnLine(networkLine* line) {
  if (!line)
    return 1E+100;
  double worst_ratio = 1E+100;
  for (auto* face : line->getBoundaryFaces()) {
    if (!face)
      continue;
    const double area = boundaryFaceArea(face);
    const double mean_area = localMeanFaceArea(face);
    if (!(area > 0.0) || !(mean_area > 0.0))
      continue;
    worst_ratio = std::min(worst_ratio, area / mean_area);
  }
  return worst_ratio;
}

inline double predictedWorstTinyFaceAreaRatioAfterFlip(networkLine* line) {
  if (!line)
    return 1E+100;
  auto faces = line->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1])
    return 1E+100;

  auto* fA = faces[0];
  auto* fB = faces[1];
  auto [p0, p1, p2] = fA->getPoints(line);
  auto [q0, q1, q2] = fB->getPoints(line);
  if (!p0 || !p1 || !p2 || !q0 || !q1 || !q2)
    return 1E+100;

  std::unordered_set<networkFace*> exclude = {fA, fB};
  const double areaA = TriangleArea(q2->X, p2->X, p1->X);
  const double areaB = TriangleArea(p2->X, q2->X, p0->X);
  if (!(areaA > 0.0) || !(areaB > 0.0))
    return 1E+100;

  const double meanA = localMeanFaceAreaForPoints({q2, p2, p1}, exclude);
  const double meanB = localMeanFaceAreaForPoints({p2, q2, p0}, exclude);
  if (!(meanA > 0.0) || !(meanB > 0.0))
    return 1E+100;

  return std::min(areaA / meanA, areaB / meanB);
}

inline bool flipTinyFaceIfImproves(networkLine* line, const double min_face_area_ratio_rel) {
  if (!line)
    return false;
  if (!line->canFlip(20.0 * M_PI / 180.0))
    return false;

  const double current_ratio = worstTinyFaceAreaRatioOnLine(line);
  const double predicted_ratio = predictedWorstTinyFaceAreaRatioAfterFlip(line);
  if (!std::isfinite(current_ratio) || !std::isfinite(predicted_ratio))
    return false;
  if (!(predicted_ratio > current_ratio + 1E-12))
    return false;
  if (!(current_ratio < min_face_area_ratio_rel || predicted_ratio > current_ratio * 1.05))
    return false;

  return line->Flip(false);
}

inline TinyFaceCheckResult checkTinyFacesRelativeToLocalMean(Network& water) {
  TinyFaceCheckResult out;
  for (auto* face : water.getBoundaryFaces()) {
    if (!face)
      continue;
    const double mean_area = localMeanFaceArea(face);
    if (!(mean_area > 0.0) || !std::isfinite(mean_area))
      continue;
    ++out.checked_faces;

    const double subface_mean_area = face->isTrueQuadraticElement ? mean_area * 0.25 : mean_area;
    if (!(subface_mean_area > 0.0) || !std::isfinite(subface_mean_area))
      continue;

    const auto tris = qualityTriangles(face);
    for (std::size_t i = 0; i < tris.size(); ++i) {
      const double area = triangleArea(tris[i]);
      if (!(area > 0.0) || !std::isfinite(area))
        continue;
      const double ratio = area / subface_mean_area;
      if (ratio >= out.worst_area_ratio)
        continue;

      out.worst_area_ratio = ratio;
      out.worst_area = area;
      out.worst_local_mean_area = subface_mean_area;
      out.worst_face = face;
      out.worst_shortest_line = nullptr;
      out.worst_subface_index = face->isTrueQuadraticElement ? static_cast<int>(i) : -1;

      double shortest_len = 1E+100;
      for (auto* l : face->getLines()) {
        if (!l)
          continue;
        const double len = l->length();
        if (!(len > 0.0) || !std::isfinite(len))
          continue;
        if (len < shortest_len) {
          shortest_len = len;
          out.worst_shortest_line = l;
        }
      }
    }
  }
  return out;
}

inline void throwIfTinyFaceRelativeToLocalMean(
    Network& water,
    const int time_step,
    const std::optional<int> rk_step,
    const double min_face_area_ratio_rel = 0.05) {
  const auto result = checkTinyFacesRelativeToLocalMean(water);
  if (!result.worst_face)
    return;
  if (!(result.worst_area_ratio < min_face_area_ratio_rel))
    return;

  const std::string rk_suffix = rk_step ? (" RK_step " + std::to_string(*rk_step)) : "";
  throw step_failure("tiny face area ratio " + std::to_string(result.worst_area_ratio) +
                     " < " + std::to_string(min_face_area_ratio_rel) +
                     " on " + water.getName() +
                     " at time_step " + std::to_string(time_step) + rk_suffix +
                     " (face_index=" + std::to_string(result.worst_face->index) +
                     ", subface_index=" + std::to_string(result.worst_subface_index) +
                     ", face_area=" + std::to_string(result.worst_area) +
                     ", local_mean_area=" + std::to_string(result.worst_local_mean_area) + ")");
}

inline void monitorTinyFaceRelativeToLocalMean(
    Network& water,
    const int time_step,
    const std::optional<int> rk_step,
    const double warn_face_area_ratio_rel = 0.1) {
  const auto result = checkTinyFacesRelativeToLocalMean(water);
  if (!result.worst_face)
    return;
  if (!(result.worst_area_ratio < warn_face_area_ratio_rel))
    return;

  std::cout << Yellow << "[tiny_face] time_step " << time_step;
  if (rk_step)
    std::cout << " RK_step " << *rk_step;
  std::cout << ": face=" << result.worst_face->index
            << " subface=" << result.worst_subface_index
            << " area_ratio=" << result.worst_area_ratio
            << " area=" << result.worst_area
            << " local_mean_area=" << result.worst_local_mean_area
            << colorReset << std::endl;
}

inline bool collapseFaceByIndexIfPossible(Network& water, const int face_index) {
  networkFace* target_face = nullptr;
  for (auto* face : water.getBoundaryFaces()) {
    if (face && face->index == face_index) {
      target_face = face;
      break;
    }
  }
  if (!target_face)
    return false;

  auto try_collapse = [&](const bool allow_corner_lines) {
    std::vector<networkLine*> candidates;
    for (auto* l : target_face->getLines()) {
      if (!l)
        continue;
      if (!allow_corner_lines && l->CORNER)
        continue;
      candidates.emplace_back(l);
    }
    std::ranges::sort(candidates, [](const auto* a, const auto* b) { return a->length() < b->length(); });
    for (auto* l : candidates)
      if (l && l->Collapse())
        return true;
    return false;
  };

  if (try_collapse(false))
    return true;
  return try_collapse(true);
}

inline FaceAltitudeDetail faceAltitudeRelativeToSharedEdgeDetail(networkFace* face, networkLine* shared_line) {
  FaceAltitudeDetail out;
  if (!face || !shared_line)
    return out;

  if (!face->isTrueQuadraticElement) {
    auto [p0, p1] = shared_line->getPoints();
    if (!p0 || !p1)
      return out;
    const double base = Norm(p1->X - p0->X);
    if (!(base > 0.) || !std::isfinite(base))
      return out;
    auto [q0, q1, q2] = face->getPoints();
    const double area = TriangleArea(q0->X, q1->X, q2->X);
    out.altitude_rel = ((area > 0.) && std::isfinite(area)) ? (2. * area / base) / base : 0.;
    return out;
  }

  const auto tris = qualityTrianglesAlongSharedLine(face, shared_line);
  for (const auto& tri : tris) {
    const double rel = triangleAltitudeRelativeToBase(tri.tri);
    if (rel < out.altitude_rel) {
      out.altitude_rel = rel;
      out.subface_index = tri.subface_index;
    }
  }
  return out;
}

inline double faceAltitudeRelativeToSharedEdge(networkFace* face, networkLine* shared_line) {
  return faceAltitudeRelativeToSharedEdgeDetail(face, shared_line).altitude_rel;
}

inline SubsurfaceAltitudeCheckResult checkSubsurfaceFaceAltitude(
    Network& water,
    const SimulationSettings::RemeshingSettings::SubsurfaceAltitudeRejectSettings& settings) {
  SubsurfaceAltitudeCheckResult out;
  if (!settings.enabled)
    return out;

  for (auto* line : water.getBoundaryLines()) {
    if (!line)
      continue;
    auto faces = line->getBoundaryFaces();
    if (faces.size() != 2 || !faces[0] || !faces[1])
      continue;
    ++out.checked_lines;

    const double dot = std::clamp(Dot(faces[0]->normal, faces[1]->normal), -1., 1.);
    const double theta_deg = 180. - std::acos(dot) * 180. / M_PI;
    if (!(settings.min_edge_angle_deg < theta_deg && theta_deg < settings.max_edge_angle_deg))
      continue;
    ++out.angle_lines;

    const auto detail0 = faceAltitudeRelativeToSharedEdgeDetail(faces[0], line);
    const auto detail1 = faceAltitudeRelativeToSharedEdgeDetail(faces[1], line);
    const double line_worst = std::min(detail0.altitude_rel, detail1.altitude_rel);
    if (line_worst < out.worst_altitude_rel) {
      out.worst_altitude_rel = line_worst;
      out.worst_angle_deg = theta_deg;
      out.worst_line = line;
      if (detail0.altitude_rel <= detail1.altitude_rel) {
        out.worst_face = faces[0];
        out.worst_subface_index = detail0.subface_index;
      } else {
        out.worst_face = faces[1];
        out.worst_subface_index = detail1.subface_index;
      }
    }
  }
  return out;
}

inline void throwIfSubsurfaceFaceAltitudeTooSmall(
    Network& water,
    const int time_step,
    const std::optional<int> rk_step,
    const SimulationSettings::RemeshingSettings::SubsurfaceAltitudeRejectSettings& settings) {
  if (!settings.enabled)
    return;
  const auto result = checkSubsurfaceFaceAltitude(water, settings);
  if (result.worst_line == nullptr)
    return;
  if (!(result.worst_altitude_rel < settings.min_face_altitude_rel))
    return;

  auto [p0, p1] = result.worst_line->getPoints();
  const std::string rk_suffix = rk_step ? (" RK_step " + std::to_string(*rk_step)) : "";
  throw step_failure("subsurface face altitude ratio " + std::to_string(result.worst_altitude_rel) +
                     " < " + std::to_string(settings.min_face_altitude_rel) +
                     " on " + water.getName() +
                     " at time_step " + std::to_string(time_step) + rk_suffix +
                     " (edge_angle_deg=" + std::to_string(result.worst_angle_deg) +
                     ", point_indices=" + std::to_string(p0 ? p0->index : -1) + "," + std::to_string(p1 ? p1->index : -1) +
                     ", face_index=" + std::to_string(result.worst_face ? result.worst_face->index : -1) +
                     ", subface_index=" + std::to_string(result.worst_subface_index) + ")");
}

inline bool flipIfOnce(Network& water, const Tdd& limit_Dirichlet, const Tdd& limit_Neumann, bool force = false, int iteration = 0,
                       const std::unordered_set<networkLine*>* skip_lines = nullptr) {
  try {
    auto [target_of_max_normal_diffD, acceptable_normal_change_by_flipD] = limit_Dirichlet;
    auto [target_of_max_normal_diffN, acceptable_normal_change_by_flipN] = limit_Neumann;
    std::cout << "flipIf" << std::endl;
    water.setGeometricPropertiesForce();
    int count = 0;
    auto V = ToVector(water.getLines());
    for (const auto& l : RandomSample(V)) {
      if (!l->CORNER) {
        if (skip_lines && skip_lines->count(l))
          continue;
        if (force && (iteration == 0 || count < iteration)) {
          if (l->Dirichlet) {
            if (l->flipIfTopologicallyBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD))
              count++;
          } else {
            if (l->flipIfTopologicallyBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN))
              count++;
          }
        } else {
          if (l->Dirichlet) {
            if (l->flipIfBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, 5))
              count++;
          } else {
            if (l->flipIfBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, 4))
              count++;
          }
        }
      }
    }
    if (count > 0) {
      std::cout << Green << "  " << count << " edges flipped." << colorReset << std::endl;
      water.setGeometricPropertiesForce();
      water.checkConnectivity();
      return true;
    }
    std::cout << Green << "  No edges flipped." << colorReset << std::endl;
    return false;
  } catch (std::exception& e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
}

inline bool flipIfBatched(Network& water, const Tdd& limit_Dirichlet, const Tdd& limit_Neumann, const char* tag = nullptr,
                          const std::unordered_set<networkLine*>* skip_lines = nullptr) {
  auto [target_of_max_normal_diffD, acceptable_normal_change_by_flipD] = limit_Dirichlet;
  auto [target_of_max_normal_diffN, acceptable_normal_change_by_flipN] = limit_Neumann;
  std::cout << "flipIf (batched)" << (tag ? std::string(" ") + tag : "") << std::endl;
  water.setGeometricPropertiesForce();
  auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };
  auto collect_non_adjacent = [&](const std::vector<networkLine*>& candidates) {
    std::vector<networkLine*> batch;
    std::unordered_set<networkPoint*> used_points;
    batch.reserve(candidates.size());
    used_points.reserve(candidates.size() * 2);
    for (auto* l : candidates) {
      if (!l)
        continue;
      auto [p0, p1] = l->getPoints();
      if (!p0 || !p1)
        continue;
      if (used_points.contains(p0) || used_points.contains(p1))
        continue;
      used_points.insert(p0);
      used_points.insert(p1);
      batch.emplace_back(l);
    }
    return batch;
  };
  auto gather_neighbor_lines = [&](const networkLine* l, std::unordered_set<networkLine*>& out) {
    if (!l)
      return;
    auto [p0, p1] = l->getPoints();
    auto add_from_point = [&](networkPoint* p) {
      if (!p)
        return;
      for (auto* nl : p->getBoundaryLines())
        if (nl)
          out.insert(nl);
    };
    add_from_point(p0);
    add_from_point(p1);
    for (auto* f : l->getBoundaryFaces())
      if (f)
        for (auto* nl : f->getLines())
          if (nl)
            out.insert(nl);
  };

  int total_flipped = 0;
  std::unordered_set<networkLine*> dirty;
  for (auto* l : water.getBoundaryLines())
    dirty.insert(l);
  for (auto iter = 0; iter < 20; iter++) {
    std::vector<networkLine*> candidates;
    candidates.reserve(dirty.size());
    for (auto* l : dirty) {
      if (!line_alive(l) || l->CORNER)
        continue;
      candidates.emplace_back(l);
    }
    if (candidates.empty())
      break;
    auto shuffled = RandomSample(candidates);
    auto batch = collect_non_adjacent(shuffled);
    if (batch.empty())
      break;
    int flipped_in_batch = 0;
    std::unordered_set<networkLine*> touched;
    for (auto* l : batch) {
      if (!line_alive(l) || l->CORNER)
        continue;
      if (skip_lines && skip_lines->count(l))
        continue;
      bool flipped = false;
      if (l->Dirichlet) {
        if (l->flipIfBetter(target_of_max_normal_diffD, acceptable_normal_change_by_flipD, 5))
          flipped = true;
      } else {
        if (l->flipIfBetter(target_of_max_normal_diffN, acceptable_normal_change_by_flipN, 4))
          flipped = true;
      }
      if (flipped) {
        flipped_in_batch++;
        total_flipped++;
        gather_neighbor_lines(l, touched);
      }
    }
    if (flipped_in_batch == 0)
      break;
    std::cout << Green << "  [batched] iter=" << iter << " batch_size=" << batch.size() << " flipped=" << flipped_in_batch << colorReset << std::endl;
    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    dirty = std::move(touched);
    if (dirty.empty())
      break;
  }
  if (total_flipped > 0) {
    std::cout << Green << "  " << total_flipped << " edges flipped (batched)." << colorReset << std::endl;
    return true;
  }
  std::cout << Green << "  No edges flipped (batched)." << colorReset << std::endl;
  return false;
}

inline void remesh_for_main_loop(Network& water, const int time_step, const double min_edge_length, const bool tetrahedralize, const bool surface_flip,
                                 const SimulationSettings::RemeshingSettings::CollisionSettings& collision_settings = {},
                                 const bool surface_split = true, const bool surface_collapse = true) {
  const double rad = M_PI / 180.0;
  const double cos_3rad = std::cos(3.0 * rad);
  const double cos_rad = std::cos(rad);
  const double global_mean_len = Mean(extLength(water.getLines()));
  const double limit_len = (min_edge_length > 0.0) ? min_edge_length : global_mean_len * 0.1;
  constexpr double min_local_face_area_ratio = 0.05;

  /* --------------------------------- 四面体の削除 --------------------------------- */
  std::cout << "孤立した四面体の削除を開始" << std::endl;
  water.DeleteIsolatedTetras();
  /* -------------------------------------------------------------------------- */
  std::cout << "内部要素の削除を開始" << std::endl;
  water.DeleteInteriorTetras();

  // Surface collision detection and resolution
  // protected_lines: lines near unresolved collision zones — do not split/collapse/flip these.
  std::unordered_set<networkLine*> protected_lines;
  bool collision_ok = detectAndResolveCollisions(water, time_step, collision_settings, protected_lines);
  if (!collision_ok) {
    // TetGen repair failed. Check mesh usability via fold ratio AND collision zone severity.
    constexpr double max_fold_ratio = 0.15;
    constexpr double max_protected_ratio = 0.5;
    auto folded = detectFoldedFaces(water, collision_settings.normal_reversal_cos, global_mean_len);
    size_t n_boundary_faces = water.getBoundaryFaces().size();
    size_t n_boundary_lines = water.getBoundaryLines().size();
    double fold_ratio = (n_boundary_faces > 0) ? static_cast<double>(folded.size()) / n_boundary_faces : 0.0;
    double protected_ratio = (n_boundary_lines > 0) ? static_cast<double>(protected_lines.size()) / n_boundary_lines : 0.0;
    if (fold_ratio > max_fold_ratio) {
      throw step_failure("collision unresolved + fold ratio " + std::to_string(fold_ratio) + " > " + std::to_string(max_fold_ratio) + " at time_step " + std::to_string(time_step));
    }
    if (protected_ratio > max_protected_ratio) {
      throw step_failure("collision unresolved + protected ratio " + std::to_string(protected_ratio) + " > " + std::to_string(max_protected_ratio) + " (" + std::to_string(protected_lines.size()) + " / " + std::to_string(n_boundary_lines) + " lines)" + " at time_step " + std::to_string(time_step));
    }
    std::cout << Yellow << "[remesh] time_step " << time_step
              << ": collision unresolved (fold_ratio=" << fold_ratio
              << ", protected_ratio=" << protected_ratio << "), continuing"
              << colorReset << std::endl;
  }

  auto build_protection_halo = [&](const std::unordered_set<networkLine*>& seeds) {
    std::unordered_set<networkLine*> halo = seeds;
    std::vector<networkLine*> queue(seeds.begin(), seeds.end());
    for (auto* l : queue) {
      if (!l)
        continue;
      auto [p0, p1] = l->getPoints();
      auto add_from_point = [&](networkPoint* p) {
        if (!p)
          return;
        for (auto* nl : p->getBoundaryLines())
          if (nl)
            halo.insert(nl);
      };
      add_from_point(p0);
      add_from_point(p1);
      for (auto* f : l->getBoundaryFaces())
        if (f)
          for (auto* nl : f->getLines())
            if (nl)
              halo.insert(nl);
    }
    return halo;
  };
  const auto protected_halo_lines = build_protection_halo(protected_lines);

  water.setGeometricPropertiesForce();
  water.checkConnectivity();

  auto is_near_protected = [&](networkLine* l) {
    if (!l)
      return false;
    return protected_halo_lines.count(l) > 0;
  };

  const std::size_t boundary_line_count = water.getBoundaryLines().size();
  // Log protected region size but do NOT skip remesh.
  // During wave breaking the wave approaches itself, causing non-adjacent collision
  // to protect most of the mesh. Skipping all remesh leads to mesh quality collapse.
  if (boundary_line_count > 0 && protected_lines.size() > 0) {
    std::cout << Yellow << "[remesh] time_step " << time_step
              << ": protected region is " << protected_lines.size() << " / " << boundary_line_count
              << " boundary lines (" << protected_halo_lines.size() << " incl. 1-ring halo)"
              << colorReset << std::endl;
  }
  const bool heavy_collision_protection = false; // never skip remesh entirely

  // Pre-remesh topology check: if lines have topology errors, skip split/collapse/flip.
  {
    int n_pre_topo_err = 0;
    for (auto* l : water.getBoundaryLines())
      if (!l->checkTopology())
        n_pre_topo_err++;
    if (n_pre_topo_err > 0) {
      std::cerr << Red << "[remesh] time_step " << time_step
                << ": " << n_pre_topo_err << " pre-existing topology errors detected"
                << colorReset << std::endl;
      throw step_failure("pre-existing topology errors (" + std::to_string(n_pre_topo_err) + " lines) at time_step " + std::to_string(time_step));
    }
  }

  const int iter_divide_collapse = 3;
  for (auto i = 0; i < iter_divide_collapse; i++) {

    // 適当な周辺の線分長さの平均を計算
    // Bench (BEM_BENCH_STEPS=3, remesh_total avg, baseline=1.00):
    // - collapse local_mean_len cache: 1.14 (slower)
    // - local_line_length vector+sort unique: 1.02 (slower)
    // - skip post-split/post-collapse setGeometricProperties/checkConnectivity when no change: 0.96 (faster)
    struct CollapseImpactMetrics {
      double length_ratio = 1E+100;
      double max_face_area_ratio = 1E+100;
      double volume_change_ratio = 1E+100;
      bool valid = false;
    };
    auto collapse_impact_metrics = [&](const networkLine* l) {
      CollapseImpactMetrics out;
      if (!l)
        return out;
      const auto faces = l->getBoundaryFaces();
      if (faces.size() != 2 || !faces[0] || !faces[1])
        return out;

      const double local_mean_len = localLineLength(l);
      const double len = l->length();
      if (!(local_mean_len > 0.0) || !(len > 0.0) || !std::isfinite(local_mean_len) || !std::isfinite(len))
        return out;

      const double area0 = boundaryFaceArea(faces[0]);
      const double area1 = boundaryFaceArea(faces[1]);
      const double mean_area0 = localMeanFaceArea(faces[0]);
      const double mean_area1 = localMeanFaceArea(faces[1]);
      if (!(area0 > 0.0) || !(area1 > 0.0) || !(mean_area0 > 0.0) || !(mean_area1 > 0.0))
        return out;
      if (!std::isfinite(area0) || !std::isfinite(area1) || !std::isfinite(mean_area0) || !std::isfinite(mean_area1))
        return out;

      const double ratio0 = area0 / mean_area0;
      const double ratio1 = area1 / mean_area1;
      out.length_ratio = len / local_mean_len;
      out.max_face_area_ratio = std::max(ratio0, ratio1);
      const double volume_before = signedPatchVolumeAroundLineAfterCollapse(l, false);
      const double volume_after = signedPatchVolumeAroundLineAfterCollapse(l, true);
      out.volume_change_ratio = std::abs(volume_after - volume_before) / std::pow(local_mean_len, 3);
      out.valid = std::isfinite(out.length_ratio) && std::isfinite(out.max_face_area_ratio) && std::isfinite(out.volume_change_ratio);
      return out;
    };
    constexpr double small_impact_length_ratio = 0.25;
    constexpr double small_impact_volume_change_ratio = 0.02;
    auto is_small_impact_collapse_candidate = [&](const networkLine* l) {
      const auto metrics = collapse_impact_metrics(l);
      if (!metrics.valid)
        return false;
      return metrics.length_ratio < small_impact_length_ratio &&
             metrics.volume_change_ratio < small_impact_volume_change_ratio;
    };
    /* --------------------------------- split (batched) --------------------------------- */
    auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };
    auto collect_non_adjacent = [&](const std::vector<networkLine*>& candidates) {
      std::vector<networkLine*> batch;
      std::unordered_set<networkPoint*> used_points;
      batch.reserve(candidates.size());
      used_points.reserve(candidates.size() * 2);
      for (auto* l : candidates) {
        if (!l)
          continue;
        auto [p0, p1] = l->getPoints();
        if (!p0 || !p1)
          continue;
        if (used_points.contains(p0) || used_points.contains(p1))
          continue;
        used_points.insert(p0);
        used_points.insert(p1);
        batch.emplace_back(l);
      }
      return batch;
    };
    auto gather_neighbor_lines = [&](const networkLine* l, std::unordered_set<networkLine*>& out) {
      if (!l)
        return;
      auto [p0, p1] = l->getPoints();
      auto add_from_point = [&](networkPoint* p) {
        if (!p)
          return;
        for (auto* nl : p->getBoundaryLines())
          if (nl)
            out.insert(nl);
      };
      add_from_point(p0);
      add_from_point(p1);
      for (auto* f : l->getBoundaryFaces())
        if (f)
          for (auto* nl : f->getLines())
            if (nl)
              out.insert(nl);
    };
    auto log_corner_connected_neumann_line = [&](const char* phase, const networkLine* l, const char* event = nullptr) {
      if (!isPostTypeCollapseTargetLine(l))
        return;
      auto [p0, p1] = l->getPoints();
      const auto faces = l->getBoundaryFaces();
      const int n_faces = static_cast<int>(faces.size());
      const double len = l->length();
      const double local_mean_len = (n_faces == 2 && faces[0] && faces[1]) ? localLineLength(l) : -1.0;
      double alt0 = -1.0, alt1 = -1.0, alt_threshold = -1.0;
      double alt_ratio0 = -1.0, alt_ratio1 = -1.0;
      double aspect_ratio0 = -1.0, aspect_ratio1 = -1.0;
      double min_angle_deg0 = -1.0, min_angle_deg1 = -1.0;
      double area_ratio0 = -1.0, area_ratio1 = -1.0, normal_dot = -2.0;
      int common_points = -1;
      int p0_lines = p0 ? static_cast<int>(p0->getLines().size()) : -1;
      int p1_lines = p1 ? static_cast<int>(p1->getLines().size()) : -1;
      int opp0_lines = -1, opp1_lines = -1;
      bool face0_dirichlet = false, face0_neumann = false;
      bool face1_dirichlet = false, face1_neumann = false;
      if (n_faces == 2 && faces[0] && faces[1] && p0 && p1) {
        auto* f0 = faces[0];
        auto* f1 = faces[1];
        face0_dirichlet = f0->Dirichlet;
        face0_neumann = f0->Neumann;
        face1_dirichlet = f1->Dirichlet;
        face1_neumann = f1->Neumann;
        auto [a, this0, b, l1, p2, l2] = f0->getPointsAndLines(const_cast<networkLine*>(l));
        auto [q0, this1, q1, e1, q2, e2] = f1->getPointsAndLines(const_cast<networkLine*>(l));
        if (this0 == l && this1 == l && a == q1 && b == q0) {
          const double area0 = boundaryFaceArea(f0);
          const double area1 = boundaryFaceArea(f1);
          const double max_edge0 = std::max({Norm(f0->getPoints()[0]->X - f0->getPoints()[1]->X),
                                             Norm(f0->getPoints()[1]->X - f0->getPoints()[2]->X),
                                             Norm(f0->getPoints()[2]->X - f0->getPoints()[0]->X)});
          const double max_edge1 = std::max({Norm(f1->getPoints()[0]->X - f1->getPoints()[1]->X),
                                             Norm(f1->getPoints()[1]->X - f1->getPoints()[2]->X),
                                             Norm(f1->getPoints()[2]->X - f1->getPoints()[0]->X)});
          alt0 = (max_edge0 > 1e-20) ? 2.0 * area0 / max_edge0 : 0.0;
          alt1 = (max_edge1 > 1e-20) ? 2.0 * area1 / max_edge1 : 0.0;
          alt_threshold = (local_mean_len > 0.0) ? 0.1 * local_mean_len : -1.0;
          alt_ratio0 = (local_mean_len > 0.0) ? alt0 / local_mean_len : -1.0;
          alt_ratio1 = (local_mean_len > 0.0) ? alt1 / local_mean_len : -1.0;
          aspect_ratio0 = faceAspectRatio(f0);
          aspect_ratio1 = faceAspectRatio(f1);
          min_angle_deg0 = minimumInteriorAngleDeg(f0);
          min_angle_deg1 = minimumInteriorAngleDeg(f1);
          const double mean_area0 = localMeanFaceArea(f0);
          const double mean_area1 = localMeanFaceArea(f1);
          area_ratio0 = (mean_area0 > 0.0) ? area0 / mean_area0 : -1.0;
          area_ratio1 = (mean_area1 > 0.0) ? area1 / mean_area1 : -1.0;
          normal_dot = Dot(f0->normal, f1->normal);
          common_points = static_cast<int>(Intersection(p0->getNeighborPointsOnSurfaces(),
                                                        p1->getNeighborPointsOnSurfaces()).size());
          opp0_lines = p2 ? static_cast<int>(p2->getLines().size()) : -1;
          opp1_lines = q2 ? static_cast<int>(q2->getLines().size()) : -1;
        }
      }
      std::cout << Magenta << "[corner_neumann_debug] " << phase;
      if (event)
        std::cout << " event=" << event;
      std::cout << " p0=" << (p0 ? p0->index : -1)
                << " p1=" << (p1 ? p1->index : -1)
                << " len=" << len
                << " local_mean_len=" << local_mean_len
                << " faces=" << n_faces
                << " alt0=" << alt0
                << " alt1=" << alt1
                << " alt_threshold=" << alt_threshold
                << " alt_ratio0=" << alt_ratio0
                << " alt_ratio1=" << alt_ratio1
                << " aspect_ratio0=" << aspect_ratio0
                << " aspect_ratio1=" << aspect_ratio1
                << " min_angle_deg0=" << min_angle_deg0
                << " min_angle_deg1=" << min_angle_deg1
                << " area_ratio0=" << area_ratio0
                << " area_ratio1=" << area_ratio1
                << " normal_dot=" << normal_dot
                << " common_points=" << common_points
                << " line_flags={D:" << l->Dirichlet << ",N:" << l->Neumann << ",C:" << l->CORNER << "}"
                << " face0_flags={D:" << face0_dirichlet << ",N:" << face0_neumann << "}"
                << " face1_flags={D:" << face1_dirichlet << ",N:" << face1_neumann << "}"
                << " endpoint0_flags={D:" << (p0 ? p0->Dirichlet : false)
                << ",N:" << (p0 ? p0->Neumann : false)
                << ",C:" << (p0 ? p0->CORNER : false) << "}"
                << " endpoint1_flags={D:" << (p1 ? p1->Dirichlet : false)
                << ",N:" << (p1 ? p1->Neumann : false)
                << ",C:" << (p1 ? p1->CORNER : false) << "}"
                << " endpoint_lines={" << p0_lines << "," << p1_lines << "}"
                << " opposite_lines={" << opp0_lines << "," << opp1_lines << "}"
                << " x0=" << (p0 ? p0->X : Tddd{0., 0., 0.})
                << " x1=" << (p1 ? p1->X : Tddd{0., 0., 0.})
                << colorReset << std::endl;
    };
    std::vector<networkLine*> corner_connected_neumann_shortlist;
    for (auto* l : water.getBoundaryLines())
      if (isPostTypeCollapseTargetLine(l))
        corner_connected_neumann_shortlist.emplace_back(l);
    std::ranges::sort(corner_connected_neumann_shortlist,
                      [](const auto* a, const auto* b) { return a->length() < b->length(); });
    if (corner_connected_neumann_shortlist.size() > 5)
      corner_connected_neumann_shortlist.resize(5);
    for (auto* l : corner_connected_neumann_shortlist)
      log_corner_connected_neumann_line("shortlist", l);
    auto should_split = [&](networkLine* l, double local_mean_len) {
      if (!l)
        return false;
      if (is_near_protected(l))
        return false;
      if (local_mean_len <= 0.)
        return false;
      auto len = l->length();
      if (len <= 2.0 * limit_len || len <= 1.4 * local_mean_len)
        return false;
      auto surfaces = l->getBoundaryFaces();
      if (surfaces.size() != 2)
        return false;

      // Split連鎖を防ぐ: 隣接線との長さ比が大きすぎる場合のみsplit
      // 最短の隣接線の2倍以上でなければsplitしない
      double min_neighbor_len = 1e20;
      for (const auto& f : surfaces) {
        for (const auto& nl : f->getLines()) {
          if (nl != l && nl->length() > 0.)
            min_neighbor_len = std::min(min_neighbor_len, static_cast<double>(nl->length()));
        }
      }
      if (min_neighbor_len < 1e19 && len < 2.0 * min_neighbor_len)
        return false;

      return Dot(surfaces[0]->normal, surfaces[1]->normal) < cos_3rad;
    };

    if (surface_split) {
      int total_splits = 0;
      const int max_splits_per_step = std::clamp(static_cast<int>(water.getBoundaryLines().size()) / 50, 10, 100);
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);
      for (auto iter = 0; iter < 20 && total_splits < max_splits_per_step; iter++) {
        std::vector<networkLine*> candidates;
        candidates.reserve(dirty.size());
        std::unordered_map<const networkLine*, double> local_mean_cache;
        local_mean_cache.reserve(dirty.size());
        for (auto* l : dirty) {
          if (!line_alive(l))
            continue;
          auto it = local_mean_cache.find(l);
          auto local_mean_len = (it != local_mean_cache.end()) ? it->second : localLineLength(l);
          if (it == local_mean_cache.end())
            local_mean_cache.emplace(l, local_mean_len);
          if (should_split(l, local_mean_len))
            candidates.emplace_back(l);
        }
        if (candidates.empty())
          break;
        auto batch = collect_non_adjacent(candidates);
        if (batch.empty())
          break;
        bool divided_any = false;
        bool split_topo_error = false;
        std::unordered_set<networkLine*> touched;
        for (auto* l : batch) {
          if (!line_alive(l))
            continue;
          auto len = l->length();
          auto local_mean_len = localLineLength(l);
          if (!should_split(l, local_mean_len))
            continue;
          // Check topology of adjacent lines BEFORE split
          bool pre_ok = true;
          for (auto* f : l->getBoundaryFaces())
            for (auto* nl : f->getLines())
              if (nl && !nl->checkTopology()) {
                pre_ok = false;
                break;
              }
          if (!pre_ok) {
            std::cerr << Yellow << "[remesh] time_step " << time_step
                      << ": skipping split (neighbor topology already bad)" << colorReset << std::endl;
            continue;
          }
          auto [sp0, sp1] = l->getPoints();
          bool near_collision = is_near_protected(l);
          std::cout << Red << "time_step " << time_step << ": splitting line"
                    << " len=" << len << " local_mean=" << local_mean_len
                    << " p0=" << sp0->X << " p1=" << sp1->X
                    << (near_collision ? " [NEAR_COLLISION_ZONE]" : "")
                    << colorReset << std::endl;
          l->Split();
          divided_any = true;
          ++total_splits;
          gather_neighbor_lines(l, touched);
          if (total_splits >= max_splits_per_step)
            break;
        }
        if (!divided_any)
          break;
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        // Post-batch topology check — stop split iterations if broken
        for (const auto& l : water.getLines()) {
          if (!l->checkTopology()) {
            split_topo_error = true;
            break;
          }
        }
        if (split_topo_error) {
          std::cerr << Red << "[remesh] time_step " << time_step
                    << ": topology error after split batch, stopping further splits" << colorReset << std::endl;
          break;
        }
        dirty = std::move(touched);
        if (dirty.empty())
          break;
      }
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    {
      bool post_split_ok = true;
      for (const auto& l : water.getLines()) {
        if (!l->checkTopology()) {
          post_split_ok = false;
          break;
        }
      }
      if (!post_split_ok) {
        throw step_failure("topology error after division at time_step " + std::to_string(time_step));
      }
    }

    /* ---------------------------------------------------------------------------- */

    if (surface_flip && !heavy_collision_protection) {
      flipIfBatched(water, {10 * rad /*target n diff*/, 10 * rad /*change n diff*/}, {3 * rad, 3 * rad}, "pre-collapse", &protected_halo_lines); // Stricter for Neumann
    }

    /* ---------------------------------- collapse --------------------------------- */

    if (surface_collapse)
      for (auto iter = 0; iter < 20; iter++) {
        bool found_small_line = false;

        // Face-quality flip: before changing DOFs with collapse, try topology-only
        // flips around locally tiny faces and keep only flips that improve the
        // local area ratio.
        if (surface_flip) {
          for (auto tiny_flip_iter = 0; tiny_flip_iter < 10; ++tiny_flip_iter) {
            std::vector<networkLine*> tiny_flip_candidates;
            std::unordered_set<networkLine*> seen;
            for (auto* face : water.getBoundaryFaces()) {
              if (!face)
                continue;
              const double area = boundaryFaceArea(face);
              if (!(area > 0.0) || !std::isfinite(area))
                continue;
              const double mean_area = localMeanFaceArea(face);
              if (!(mean_area > 0.0) || !std::isfinite(mean_area))
                continue;
              const double ratio = area / mean_area;
              if (!(ratio < min_local_face_area_ratio))
                continue;

              for (auto* l : face->getLines()) {
                if (!l)
                  continue;
                if (!line_alive(l) || l->CORNER || is_near_protected(l))
                  continue;
                if (seen.insert(l).second)
                  tiny_flip_candidates.emplace_back(l);
              }
            }
            if (tiny_flip_candidates.empty())
              break;

            auto batch = collect_non_adjacent(RandomSample(tiny_flip_candidates));
            if (batch.empty())
              break;

            bool changed_tiny_faces = false;
            std::unordered_set<networkLine*> touched;
            for (auto* l : batch) {
              if (!line_alive(l) || l->CORNER || is_near_protected(l))
                continue;
              const double before_ratio = worstTinyFaceAreaRatioOnLine(l);
              const double predicted_ratio = predictedWorstTinyFaceAreaRatioAfterFlip(l);
              if (flipTinyFaceIfImproves(l, min_local_face_area_ratio)) {
                found_small_line = true;
                changed_tiny_faces = true;
                gather_neighbor_lines(l, touched);
                std::cout << "time_step " << time_step
                          << ": line flipped due to tiny local-area face. ratio_before="
                          << before_ratio << " ratio_after_pred=" << predicted_ratio << std::endl;
              }
            }
            if (!changed_tiny_faces)
              break;
            water.setGeometricPropertiesForce();
            water.checkConnectivity();
          }
        }

        // Face-quality collapse: collapse lines whose geometric impact is small.
        // Primary target is the shortest edge of tiny faces, but also allow
        // directly collapsing very short lines whose both adjacent faces are
        // small relative to their local neighborhoods.
        for (auto tiny_iter = 0; tiny_iter < 10; ++tiny_iter) {
          std::vector<networkLine*> tiny_candidates;
          std::unordered_set<networkLine*> seen;
          for (auto* l : water.getBoundaryLines()) {
            if (!l)
              continue;
            if (!line_alive(l) || l->CORNER || is_near_protected(l))
              continue;
            if (!is_small_impact_collapse_candidate(l))
              continue;
            if (seen.insert(l).second)
              tiny_candidates.emplace_back(l);
          }
          for (auto* face : water.getBoundaryFaces()) {
            if (!face)
              continue;
            const double area = boundaryFaceArea(face);
            if (!(area > 0.0) || !std::isfinite(area))
              continue;
            const double mean_area = localMeanFaceArea(face);
            if (!(mean_area > 0.0) || !std::isfinite(mean_area))
              continue;
            const double ratio = area / mean_area;
            if (!(ratio < min_local_face_area_ratio))
              continue;

            networkLine* shortest = nullptr;
            double shortest_len = 1E+100;
            for (auto* l : face->getLines()) {
              if (!l)
                continue;
              if (!line_alive(l) || l->CORNER || is_near_protected(l))
                continue;
              const double len = l->length();
              if (!(len > 0.0) || !std::isfinite(len))
                continue;
              if (len < shortest_len) {
                shortest_len = len;
                shortest = l;
              }
            }
            if (shortest && seen.insert(shortest).second)
              tiny_candidates.emplace_back(shortest);
          }
          if (tiny_candidates.empty())
            break;

          auto batch = collect_non_adjacent(tiny_candidates);
          if (batch.empty())
            break;

          bool changed_tiny_faces = false;
          for (auto* l : batch) {
            if (!line_alive(l) || l->CORNER || is_near_protected(l))
              continue;
            const auto attached_faces = l->getBoundaryFaces();
            double worst_ratio = 1E+100;
            for (auto* face : attached_faces) {
              const double area = boundaryFaceArea(face);
              const double mean_area = localMeanFaceArea(face);
              if (!(area > 0.0) || !(mean_area > 0.0))
                continue;
              worst_ratio = std::min(worst_ratio, area / mean_area);
            }
            const auto metrics = collapse_impact_metrics(l);
            if (l->Collapse()) {
              found_small_line = true;
              changed_tiny_faces = true;
              std::cout << "time_step " << time_step
                        << ": line merged due to tiny local-area face. worst_ratio="
                        << worst_ratio;
              if (metrics.valid)
                std::cout << " len_ratio=" << metrics.length_ratio
                          << " max_face_ratio=" << metrics.max_face_area_ratio
                          << " volume_change_ratio=" << metrics.volume_change_ratio;
              std::cout << std::endl;
            }
          }
          if (!changed_tiny_faces)
            break;
          water.setGeometricPropertiesForce();
          water.checkConnectivity();
        }

        std::unordered_set<networkLine*> dirty;
        for (auto* l : water.getBoundaryLines())
          dirty.insert(l);
        while (!dirty.empty()) {
          std::vector<networkLine*> candidates;
          candidates.reserve(dirty.size());
          for (auto* l : dirty) {
            if (!line_alive(l))
              continue;
            candidates.emplace_back(l);
          }
          if (candidates.empty())
            break;
          auto batch = collect_non_adjacent(candidates);
          if (batch.empty())
            break;
          bool changed_in_batch = false;
          std::unordered_set<networkLine*> touched;
          for (auto* l : batch) {
            if (!line_alive(l))
              continue;

            // Skip CORNER lines (similar to flip operations)
            if (l->CORNER)
              continue;
            // Skip lines near unresolved collision zones
            if (is_near_protected(l))
              continue;

            auto surfaces = l->getBoundaryFaces();
            if (surfaces.size() != 2)
              continue;
            auto f0 = surfaces[0];
            auto f1 = surfaces[1];
            if (!f0 || !f1)
              continue;
            auto len = l->length();

            auto [a, b, c] = f0->getPoints(l);
            auto [_, __, d] = f1->getPoints(l);

            auto local_mean_len = localLineLength(l);

            // Aspect-ratio collapse: if either adjacent triangle has very low height
            // (altitude < 10% of local mean edge length), collapse regardless of Neumann.
            // This prevents mesh crushing near wave-making paddles and structures.
            {
              double area0 = TriangleArea(f0->getPoints()[0]->X, f0->getPoints()[1]->X, f0->getPoints()[2]->X);
              double area1 = TriangleArea(f1->getPoints()[0]->X, f1->getPoints()[1]->X, f1->getPoints()[2]->X);
              // altitude = 2 * area / base (longest edge of the triangle)
              double max_edge0 = std::max({Norm(f0->getPoints()[0]->X - f0->getPoints()[1]->X),
                                           Norm(f0->getPoints()[1]->X - f0->getPoints()[2]->X),
                                           Norm(f0->getPoints()[2]->X - f0->getPoints()[0]->X)});
              double max_edge1 = std::max({Norm(f1->getPoints()[0]->X - f1->getPoints()[1]->X),
                                           Norm(f1->getPoints()[1]->X - f1->getPoints()[2]->X),
                                           Norm(f1->getPoints()[2]->X - f1->getPoints()[0]->X)});
              double alt0 = (max_edge0 > 1e-20) ? 2.0 * area0 / max_edge0 : 0.0;
              double alt1 = (max_edge1 > 1e-20) ? 2.0 * area1 / max_edge1 : 0.0;
              double alt_threshold = local_mean_len * 0.1;
              double mean_area0 = localMeanFaceArea(f0);
              double mean_area1 = localMeanFaceArea(f1);
              double area_ratio0 = (mean_area0 > 0.0) ? area0 / mean_area0 : 1E+100;
              double area_ratio1 = (mean_area1 > 0.0) ? area1 / mean_area1 : 1E+100;
              double alt_ratio0 = (local_mean_len > 0.0) ? alt0 / local_mean_len : 1E+100;
              double alt_ratio1 = (local_mean_len > 0.0) ? alt1 / local_mean_len : 1E+100;
              double min_angle_deg0 = minimumInteriorAngleDeg(f0);
              double min_angle_deg1 = minimumInteriorAngleDeg(f1);
              if (!l->CORNER &&
                  std::min(alt_ratio0, alt_ratio1) < 0.2 &&
                  std::min(area_ratio0, area_ratio1) < 0.2 &&
                  std::min(min_angle_deg0, min_angle_deg1) < 20.0) {
                gather_neighbor_lines(l, touched);
                if (std::ranges::find(corner_connected_neumann_shortlist, l) != corner_connected_neumann_shortlist.end())
                  log_corner_connected_neumann_line("collapse", l, "corner_neumann_area_alt_attempt");
                if (l->Collapse()) {
                  changed_in_batch = true;
                  if (std::ranges::find(corner_connected_neumann_shortlist, l) != corner_connected_neumann_shortlist.end())
                    std::cout << Magenta << "[corner_neumann_debug] collapse event=corner_neumann_area_alt_success" << colorReset << std::endl;
                  continue;
                }
                if (std::ranges::find(corner_connected_neumann_shortlist, l) != corner_connected_neumann_shortlist.end())
                  std::cout << Magenta << "[corner_neumann_debug] collapse event=corner_neumann_area_alt_failed"
                            << " alt_ratio_min=" << std::min(alt_ratio0, alt_ratio1)
                            << " area_ratio_min=" << std::min(area_ratio0, area_ratio1)
                            << " min_angle_deg_min=" << std::min(min_angle_deg0, min_angle_deg1)
                            << colorReset << std::endl;
              }
              if (alt0 < alt_threshold || alt1 < alt_threshold) {
                gather_neighbor_lines(l, touched);
                if (std::ranges::find(corner_connected_neumann_shortlist, l) != corner_connected_neumann_shortlist.end())
                  log_corner_connected_neumann_line("collapse", l, "low_altitude_attempt");
                if (l->Collapse()) {
                  found_small_line = true;
                  changed_in_batch = true;
                  std::cout << "time_step " << time_step << ": line merged due to low altitude."
                            << " alt0=" << alt0 << " alt1=" << alt1
                            << " threshold=" << alt_threshold << std::endl;
                  continue;
                }
                if (std::ranges::find(corner_connected_neumann_shortlist, l) != corner_connected_neumann_shortlist.end())
                  log_corner_connected_neumann_line("collapse", l, "low_altitude_failed");
              }
            }

            // Neumann lines should only be collapsed if extremely small to prevent large normal changes
            if (l->Neumann && len > limit_len * 0.5)
              continue;

            // 線の長さが平均の0.4倍以下ならマージ (dtの低下を防ぐため閾値を上げる)
            if ((len < local_mean_len * 0.4 || len < limit_len)) {
              gather_neighbor_lines(l, touched);
              l->Collapse();
              found_small_line = true;
              changed_in_batch = true;
              std::cout << "time_step " << time_step << ": line merged due to small length. length = " << len << ", local mean length = " << local_mean_len << std::endl;
              continue;
            }
            //
            double volume = TetrahedronVolume(a->X, b->X, c->X, d->X);
            double ref_volume = std::pow(local_mean_len, 3) * 0.5;
            (void)volume;
            (void)ref_volume;

            auto [p0, p1, p2] = f0->getPoints(l);
            auto ci0 = CircumradiusToInradius(p0->X, p1->X, p2->X);
            auto [q0, q1, q2] = f1->getPoints(l);
            auto ci1 = CircumradiusToInradius(q0->X, q1->X, q2->X);

            // For Neumann lines, check that normals are nearly opposite before allowing collapse
            const double cos_5rad = std::cos(5.0 * rad);
            if (l->Neumann) {
              // Only allow collapse if normals are truly opposite (angle > 175 degrees)
              if (Dot(f0->normal, -f1->normal) < cos_5rad)
                continue;
            }

            // 接する面の法線が反対方向でかつ近ければマージ
            if (l->length() < 0.1 * (Norm(p2->X - p1->X) + Norm(p2->X - p0->X) + Norm(q2->X - q1->X) + Norm(q2->X - q0->X)) / 4.)
              if (l->Collapse()) {
                found_small_line = true;
                changed_in_batch = true;
                std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
                continue;
              }

            if (ci0 > 50. || ci1 > 50.)
              if (l->Collapse()) {
                found_small_line = true;
                changed_in_batch = true;
                std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
                continue;
              }
            if (Norm(p2->X - q2->X) < std::sqrt(TriangleArea(p0->X, p1->X, p2->X) + TriangleArea(q0->X, q1->X, q2->X)) / 10.) {
              l->Flip(true); // できればでいい
              l->Collapse();
              found_small_line = true;
              changed_in_batch = true;
              std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
              continue;
            }
            if (Dot(f0->normal, -f1->normal) > cos_rad)
              if (l->Collapse()) {
                found_small_line = true;
                changed_in_batch = true;
                std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
                continue;
              }
          }
          if (changed_in_batch) {
            water.setGeometricPropertiesForce();
            water.checkConnectivity();
          }
          if (!changed_in_batch)
            break;
          dirty = std::move(touched);
        }
        if (!found_small_line)
          break;
      }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();

    {
      bool post_merge_ok = true;
      for (const auto& l : water.getLines()) {
        if (!l->checkTopology()) {
          post_merge_ok = false;
          break;
        }
      }
      if (!post_merge_ok) {
        throw step_failure("topology error after merge at time_step " + std::to_string(time_step));
      }
    }
    // end if (surface_collapse) — the 'if' on line 452 controls the outer for-loop only;
    // topology check and setGeometricProperties are inside the for-loop body so they are also gated.

    /* --------------------------- エッジフリップによる表面メッシュの改善-------------------------- */

    if (surface_flip && !heavy_collision_protection) {
      flipIfBatched(water, {20 * rad /*target n diff*/, 20 * rad /*change n diff*/}, {5 * rad, 5 * rad}, "post-collapse", &protected_halo_lines); // Stricter for Neumann
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();

    {
      bool post_flip_ok = true;
      for (const auto& l : water.getLines()) {
        if (!l->checkTopology()) {
          post_flip_ok = false;
          break;
        }
      }
      if (!post_flip_ok) {
        throw step_failure("topology error after flip at time_step " + std::to_string(time_step));
      }
    }
  }

  // Post-remesh fold check: detect folds introduced by split/collapse/flip
  {
    constexpr double max_fold_ratio = 0.15;
    auto folded = detectFoldedFaces(water, collision_settings.normal_reversal_cos, global_mean_len);
    size_t n_boundary_faces = water.getBoundaryFaces().size();
    double fold_ratio = (n_boundary_faces > 0) ? static_cast<double>(folded.size()) / n_boundary_faces : 0.0;
    if (!folded.empty()) {
      std::cout << Yellow << "[remesh] time_step " << time_step
                << ": post-remesh fold check: " << folded.size() << " / " << n_boundary_faces
                << " folded faces (ratio=" << fold_ratio << ")"
                << colorReset << std::endl;
    }
    if (fold_ratio > max_fold_ratio) {
      throw step_failure("post-remesh fold ratio " + std::to_string(fold_ratio) + " > " + std::to_string(max_fold_ratio) + " at time_step " + std::to_string(time_step));
    }
  }
  {
    const auto tiny = checkTinyFacesRelativeToLocalMean(water);
    if (tiny.worst_face && tiny.worst_area_ratio < min_local_face_area_ratio) {
      std::cout << Yellow << "[remesh] time_step " << time_step
                << ": tiny face check: face=" << tiny.worst_face->index
                << " area_ratio=" << tiny.worst_area_ratio
                << " area=" << tiny.worst_area
                << " local_mean_area=" << tiny.worst_local_mean_area
                << colorReset << std::endl;
      throw step_failure("post-remesh tiny face area ratio " + std::to_string(tiny.worst_area_ratio) + " < " + std::to_string(min_local_face_area_ratio) +
                         " at time_step " + std::to_string(time_step));
    }
  }

#ifdef USE_TETGEN
  if (tetrahedralize)
    water.tetrahedralize();
#endif
  water.setGeometricPropertiesForce();

  water.improveTetrahedraDelaunay();
}
