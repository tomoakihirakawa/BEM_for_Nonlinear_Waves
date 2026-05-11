#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

#include "BEM_BoundaryValues.hpp"
#include "BEM_mesh_adjustment.hpp"
#include "BEM_setBoundaryTypes.hpp"
#include "Network.hpp"
#include "basic.hpp"
#include "basic_geometry.hpp"
#include "searcher.hpp"

namespace BEMMeshPipeline {

struct WaterlineQuery {
  networkLine* l = nullptr;
  Tddd X_linear = {0., 0., 0.};
  double move_limit_factor = 0.3;
  int max_iter = 4;
  double tol_relative = 1e-6;
  double max_body_gap_factor = 5.0;
  double contact_range_factor = 2.0;
};

enum class WaterlineStatus {
  ok,
  ok_via_fallback,
  skipped_not_corner,
  no_dirichlet_face,
  no_body_candidate,
  no_convergence,
  body_gap_too_large,
  no_move_possible,
  invalid_projection
};

struct WaterlineResult {
  WaterlineStatus status = WaterlineStatus::skipped_not_corner;
  Tddd target_X = {0., 0., 0.};
  Tddd target_X_clamped = {0., 0., 0.};
  Tddd delta_clamped = {0., 0., 0.};
  double free_surface_gap = 0.0;
  double body_gap = 0.0;
  double body_gap_threshold = 0.0;
  double shift_limit = 0.0;
  double move_ratio = 0.0;
  int iterations_used = 0;
  bool used_fallback = false;
};

inline const char* toString(WaterlineStatus status) {
  switch (status) {
  case WaterlineStatus::ok:
    return "ok";
  case WaterlineStatus::ok_via_fallback:
    return "ok_via_fallback";
  case WaterlineStatus::skipped_not_corner:
    return "skipped_not_corner";
  case WaterlineStatus::no_dirichlet_face:
    return "no_dirichlet_face";
  case WaterlineStatus::no_body_candidate:
    return "no_body_candidate";
  case WaterlineStatus::no_convergence:
    return "no_convergence";
  case WaterlineStatus::body_gap_too_large:
    return "body_gap_too_large";
  case WaterlineStatus::no_move_possible:
    return "no_move_possible";
  case WaterlineStatus::invalid_projection:
    return "invalid_projection";
  }
  return "unknown";
}

inline T3Tddd faceTriangle(const networkFace* f) {
  auto [p0, p1, p2] = f->getPoints();
  return {p0->X, p1->X, p2->X};
}

inline Tddd nearestOnFace(const networkFace* f, const Tddd& X) {
  return Nearest(X, faceTriangle(f));
}

inline double distanceToFace(const networkFace* f, const Tddd& X) {
  return Norm(nearestOnFace(f, X) - X);
}

inline std::vector<networkFace*> uniqueFaces(const std::vector<networkFace*>& faces) {
  std::vector<networkFace*> out;
  std::unordered_set<networkFace*> seen;
  out.reserve(faces.size());
  for (auto* f : faces) {
    if (!f || seen.contains(f))
      continue;
    seen.insert(f);
    out.push_back(f);
  }
  return out;
}

inline std::vector<networkFace*> faceSetToVector(const std::unordered_set<networkFace*>& faces) {
  std::vector<networkFace*> out;
  out.reserve(faces.size());
  for (auto* f : faces) {
    if (f)
      out.push_back(f);
  }
  return out;
}

inline std::vector<networkFace*> unionFaces(std::vector<networkFace*> a, const std::vector<networkFace*>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return uniqueFaces(a);
}

inline std::pair<networkFace*, double> nearestFaceByDistance(const std::vector<networkFace*>& faces, const Tddd& X) {
  networkFace* best = nullptr;
  double best_dist = std::numeric_limits<double>::infinity();
  for (auto* f : faces) {
    if (!f)
      continue;
    const double d = distanceToFace(f, X);
    if (std::isfinite(d) && d < best_dist) {
      best_dist = d;
      best = f;
    }
  }
  return {best, best_dist};
}

inline Tddd nearestOnFaces(const std::vector<networkFace*>& faces, const Tddd& X, double* distance = nullptr) {
  auto [best, best_dist] = nearestFaceByDistance(faces, X);
  if (!best) {
    if (distance)
      *distance = std::numeric_limits<double>::infinity();
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  }
  if (distance)
    *distance = best_dist;
  return nearestOnFace(best, X);
}

inline double localMeanEdgeLength(const networkLine* l) {
  if (!l)
    return 1.0;
  double sum = 0.0;
  int count = 0;
  for (const auto* f : l->getBoundaryFaces()) {
    if (!f)
      continue;
    auto [p0, p1, p2] = f->getPoints();
    const std::array<double, 3> lengths = {
        Norm(p0->X - p1->X),
        Norm(p1->X - p2->X),
        Norm(p2->X - p0->X)};
    for (const double len : lengths) {
      if (std::isfinite(len) && len > 0.0) {
        sum += len;
        ++count;
      }
    }
  }
  if (count > 0)
    return sum / static_cast<double>(count);
  auto [pA, pB] = l->getPoints();
  const double len = (pA && pB) ? Norm(pA->X - pB->X) : 0.0;
  return (std::isfinite(len) && len > 0.0) ? len : 1.0;
}

inline std::vector<networkFace*> dirichletFacesOfLine(const networkLine* l) {
  std::vector<networkFace*> faces;
  if (!l)
    return faces;
  for (auto* f : l->getBoundaryFaces()) {
    if (f && getNodeFaceBoundaryType(l, f) == NodeFaceBoundaryType::Dirichlet)
      faces.push_back(f);
  }
  return uniqueFaces(faces);
}

inline std::vector<networkFace*> endpointContactFaces(const networkLine* l) {
  if (!l)
    return {};
  auto [pA, pB] = l->getPoints();
  std::vector<networkFace*> faces;
  if (pA)
    faces = unionFaces(faces, getEffectiveContactFaces(pA));
  if (pB)
    faces = unionFaces(faces, getEffectiveContactFaces(pB));
  return uniqueFaces(faces);
}

inline std::vector<networkFace*> bfsExpandedFaces(const std::vector<networkFace*>& faces, short unsigned int depth) {
  if (faces.empty())
    return {};
  return faceSetToVector(bfs(faces, depth));
}

inline std::vector<networkFace*> selectBodyCandidates(const networkLine* l,
                                                      const Tddd& X_linear,
                                                      double body_gap_threshold,
                                                      bool* used_fallback) {
  if (used_fallback)
    *used_fallback = false;

  const std::vector<networkFace*> line_contacts = uniqueFaces(getEffectiveContactFaces(l));
  auto [line_best, line_dist] = nearestFaceByDistance(line_contacts, X_linear);
  if (line_best && line_dist <= body_gap_threshold)
    return line_contacts;

  const std::vector<networkFace*> endpoint_contacts = endpointContactFaces(l);
  auto [endpoint_best, endpoint_dist] = nearestFaceByDistance(endpoint_contacts, X_linear);
  if (endpoint_best && endpoint_dist <= body_gap_threshold) {
    if (used_fallback)
      *used_fallback = true;
    return endpoint_contacts;
  }

  std::unordered_set<networkFace*> seed;
  seed.insert(line_contacts.begin(), line_contacts.end());
  seed.insert(endpoint_contacts.begin(), endpoint_contacts.end());
  const std::vector<networkFace*> expanded(seed.begin(), seed.end());
  const std::vector<networkFace*> bfs_expanded = bfsExpandedFaces(expanded, 2);
  auto [expanded_best, expanded_dist] = nearestFaceByDistance(bfs_expanded, X_linear);
  if (expanded_best && expanded_dist <= body_gap_threshold) {
    if (used_fallback)
      *used_fallback = true;
    return bfs_expanded;
  }

  return {};
}

inline WaterlineResult queryWaterlineGeometry(const WaterlineQuery& q) {
  WaterlineResult result;
  if (!q.l || !q.l->BCInterface)
    return result;

  const auto dirichlet_faces = dirichletFacesOfLine(q.l);
  if (dirichlet_faces.empty()) {
    result.status = WaterlineStatus::no_dirichlet_face;
    return result;
  }

  const double eps = 1e-12;
  result.shift_limit = currentShiftLimit(q.l, q.move_limit_factor);
  if (!std::isfinite(result.shift_limit) || result.shift_limit <= eps) {
    result.status = WaterlineStatus::no_move_possible;
    return result;
  }

  const double contact_range = q.l->contact_range;
  if (std::isfinite(contact_range) && contact_range > 0.0)
    result.body_gap_threshold = std::min(q.max_body_gap_factor * result.shift_limit,
                                         q.contact_range_factor * contact_range);
  else
    result.body_gap_threshold = q.max_body_gap_factor * result.shift_limit;

  bool used_fallback = false;
  const auto body_candidates = selectBodyCandidates(q.l, q.X_linear, result.body_gap_threshold, &used_fallback);
  if (body_candidates.empty()) {
    result.status = WaterlineStatus::no_body_candidate;
    return result;
  }
  result.used_fallback = used_fallback;

  WaterlineStatus status = WaterlineStatus::ok;
  Tddd X = q.X_linear;
  const double tol = std::max(q.tol_relative * localMeanEdgeLength(q.l), eps);
  const int max_iter = std::max(1, q.max_iter);
  bool converged = false;

  for (int iter = 0; iter < max_iter; ++iter) {
    double free_gap = 0.0;
    const Tddd X_D = nearestOnFaces(dirichlet_faces, X, &free_gap);
    double body_gap = 0.0;
    const Tddd X_B = nearestOnFaces(body_candidates, X_D, &body_gap);
    if (!isFinite(X_D) || !isFinite(X_B) || !std::isfinite(free_gap) || !std::isfinite(body_gap)) {
      result.status = WaterlineStatus::invalid_projection;
      return result;
    }
    result.iterations_used = iter + 1;
    result.free_surface_gap = free_gap;
    result.body_gap = body_gap;
    if (Norm(X_B - X) < tol) {
      X = X_B;
      converged = true;
      break;
    }
    X = X_B;
  }

  if (!converged)
    status = WaterlineStatus::no_convergence;

  result.target_X = X;
  auto [final_body_face, final_body_dist] = nearestFaceByDistance(body_candidates, X);
  if (!final_body_face || !std::isfinite(final_body_dist)) {
    result.status = WaterlineStatus::invalid_projection;
    return result;
  }
  result.body_gap = final_body_dist;
  double free_gap = 0.0;
  (void)nearestOnFaces(dirichlet_faces, X, &free_gap);
  result.free_surface_gap = free_gap;
  if (!std::isfinite(result.body_gap) || result.body_gap > result.body_gap_threshold) {
    result.status = WaterlineStatus::body_gap_too_large;
    return result;
  }

  const Tddd delta = X - q.X_linear;
  const double delta_norm = Norm(delta);
  if (std::isfinite(delta_norm) && delta_norm > result.shift_limit)
    result.delta_clamped = result.shift_limit * Normalize(delta);
  else
    result.delta_clamped = delta;
  result.target_X_clamped = q.X_linear + result.delta_clamped;
  result.move_ratio = Norm(result.delta_clamped) / result.shift_limit;
  if (!isFinite(result.target_X_clamped) || !isFinite(result.delta_clamped) || !std::isfinite(result.move_ratio)) {
    result.status = WaterlineStatus::invalid_projection;
    return result;
  }

  result.status = (status == WaterlineStatus::no_convergence)
                      ? WaterlineStatus::no_convergence
                      : (used_fallback ? WaterlineStatus::ok_via_fallback : WaterlineStatus::ok);
  return result;
}

} // namespace BEMMeshPipeline
