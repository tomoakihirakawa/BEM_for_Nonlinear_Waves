#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "BEM_mesh_adjustment.hpp"
#include "BEM_reference_state.hpp"

namespace BEMMeshPipeline {

enum class GeometryProjectionMode {
  LinearFace,
  // Uses the 6-node geometry stored in ReferenceState.  For linear elements
  // this currently evaluates the saved edge-midpoint geometry directly; a true
  // local WLS patch fit can replace this mode later without changing callers.
  ReferenceQuadraticFace,
  // Compatibility alias for the original plan wording.  It is intentionally
  // not used as the default name until a real patch-fit implementation exists.
  QuadraticFitFromLinearVertices = ReferenceQuadraticFace,
  TrueQuadraticFace
};

enum class GeometryTargetStatus {
  ok,
  ok_via_fallback,
  no_convergence,
  not_applicable,
  no_reference_surface,
  no_body_surface,
  body_gap_too_large,
  no_move_possible,
  invalid_projection
};

enum class GeometryFallbackLevel {
  none,
  linear_face,
  endpoint_contacts,
  bfs_contacts,
  global_body
};

struct GeometryProjectorOptions {
  GeometryProjectionMode mode = GeometryProjectionMode::ReferenceQuadraticFace;
  double move_limit_factor = 0.3;
  int max_iter = 4;
  double tol_relative = 1e-6;
  double max_body_gap_factor = 5.0;
  double contact_range_factor = 2.0;
  double trust_radius_factor = 5.0;
  int body_bfs_fallback_depth = 0;
  bool enable_global_body_fallback = false;
  double global_body_fallback_range_factor = 1.0;
  double feature_angle_rad = 60.0 * M_PI / 180.0;
  bool anisotropic_body_nearest_for_bcinterface = false;
  double anisotropic_body_nearest_tangent_axis_factor = 0.75;
  double anisotropic_body_nearest_normal_axis_factor = 1.0;
  bool bcinterface_local_dirichlet_first = true;
};

struct GeometryProjectorQuery {
  networkPoint* p = nullptr;
  networkLine* l = nullptr;
  const Network* water = nullptr;
  Tddd X_in = {0., 0., 0.};
};

struct GeometryTarget {
  GeometryTargetStatus status = GeometryTargetStatus::not_applicable;
  GeometryFallbackLevel fallback = GeometryFallbackLevel::none;
  Tddd target_X = {0., 0., 0.};
  Tddd target_X_clamped = {0., 0., 0.};
  Tddd delta_clamped = {0., 0., 0.};
  const PreModificationSnapshot::FaceData* reference_face = nullptr;
  Tdd param_uv = {0., 0.};
  networkFace* body_face = nullptr;
  double dirichlet_gap = 0.0;
  double body_gap = 0.0;
  double body_gap_threshold = 0.0;
  double shift_limit = 0.0;
  double move_ratio = 0.0;
  int iterations = 0;
  bool used_quadratic_fit = false;
};

struct TargetThetaResult {
  bool valid = false;
  double theta = 0.0;
  bool used_reference = false;
  bool used_body = false;
};

enum class ProviderFieldKind {
  Phi,
  PhiT
};

struct ProviderFieldSample {
  bool ok = false;
  ProviderFieldKind kind = ProviderFieldKind::Phi;
  double value = 0.0;
};

namespace GeometryProjectorDetail {

inline double localLength(const networkPoint* p) {
  if (!p)
    return 1.0;
  double sum = 0.0;
  int count = 0;
  for (auto* l : p->getBoundaryLines()) {
    if (!l)
      continue;
    auto [p0, p1] = l->getPoints();
    const double len = (p0 && p1) ? Norm(p0->X - p1->X) : 0.0;
    if (std::isfinite(len) && len > 0.0) {
      sum += len;
      ++count;
    }
  }
  return count > 0 ? sum / static_cast<double>(count) : 1.0;
}

inline double localLength(const networkLine* l) {
  if (!l)
    return 1.0;
  auto [p0, p1] = l->getPoints();
  const double len = (p0 && p1) ? Norm(p0->X - p1->X) : 0.0;
  return (std::isfinite(len) && len > 0.0) ? len : 1.0;
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
  for (auto* f : faces)
    if (f)
      out.push_back(f);
  return out;
}

template <class Entity>
inline bool isBCInterfaceEntity(const Entity* entity) {
  if constexpr (requires { entity->BCInterface; })
    return entity && entity->BCInterface;
  else if constexpr (requires { entity->BCInterface; })
    return entity && entity->BCInterface;
  else
    return false;
}

template <class Entity>
inline bool hasSharpGeometry(const Entity* entity) {
  if constexpr (requires { entity->SharpQ(); })
    return entity && entity->SharpQ();
  else
    return false;
}

inline std::vector<networkFace*> endpointContactFaces(const networkLine* l) {
  if (!l)
    return {};
  auto [p0, p1] = l->getPoints();
  std::vector<networkFace*> faces;
  if (p0) {
    const auto a = getEffectiveContactFaces(p0);
    faces.insert(faces.end(), a.begin(), a.end());
  }
  if (p1) {
    const auto b = getEffectiveContactFaces(p1);
    faces.insert(faces.end(), b.begin(), b.end());
  }
  return uniqueFaces(faces);
}

template <class Entity>
inline bool useAnisotropicBodyNearest(const Entity* entity,
                                      const GeometryProjectorOptions* options) {
  return options && options->anisotropic_body_nearest_for_bcinterface &&
         isBCInterfaceEntity(entity) &&
         entity &&
         entity->contact_range > 0.0 &&
         std::isfinite(entity->contact_range);
}

template <class Entity>
inline double anisotropicBodyNearestTangentAxisLength(const Entity* entity,
                                                      const GeometryProjectorOptions* options) {
  if (!entity || !options)
    return 0.0;
  return std::max(options->anisotropic_body_nearest_tangent_axis_factor * entity->contact_range,
                  1e-12);
}

template <class Entity>
inline double anisotropicBodyNearestNormalAxisLength(const Entity* entity,
                                                     const GeometryProjectorOptions* options) {
  if (!entity || !options)
    return 0.0;
  return std::max(options->anisotropic_body_nearest_normal_axis_factor * entity->contact_range,
                  1e-12);
}

template <class Entity>
inline Tddd bodyFaceNearestPoint(const Entity* entity,
                                 const GeometryProjectorOptions* options,
                                 const Tddd& X,
                                 const networkFace* f) {
  if (f && useAnisotropicBodyNearest(entity, options))
    return NearestAnisotropic(X, f,
                              anisotropicBodyNearestNormalAxisLength(entity, options),
                              anisotropicBodyNearestTangentAxisLength(entity, options));
  return Nearest(X, f);
}

template <class Entity>
inline double bodyFaceSelectionDistance(const Entity* entity,
                                        const GeometryProjectorOptions* options,
                                        const Tddd& X,
                                        const networkFace* f,
                                        const Tddd& nearest) {
  const double euclidean = Norm(nearest - X);
  if (!useAnisotropicBodyNearest(entity, options))
    return euclidean;
  return contactAnisotropicDistance(X, nearest, f ? f->normal : Tddd{0., 0., 0.},
                                    entity->contact_range,
                                    options->anisotropic_body_nearest_normal_axis_factor,
                                    options->anisotropic_body_nearest_tangent_axis_factor);
}

struct BodyFaceNearest {
  networkFace* face = nullptr;
  Tddd X = {0., 0., 0.};
  double distance = std::numeric_limits<double>::infinity();
  double metric = std::numeric_limits<double>::infinity();
};

template <class Entity>
inline BodyFaceNearest nearestFaceByDistance(const std::vector<networkFace*>& faces,
                                             const Tddd& X,
                                             const GeometryProjectorOptions* options,
                                             const Entity* entity) {
  BodyFaceNearest best;
  for (auto* f : faces) {
    if (!f)
      continue;
    const auto nearest = bodyFaceNearestPoint(entity, options, X, f);
    const double d = Norm(nearest - X);
    const double metric = bodyFaceSelectionDistance(entity, options, X, f, nearest);
    if (std::isfinite(d) && std::isfinite(metric) &&
        (metric < best.metric || (metric == best.metric && d < best.distance))) {
      best.face = f;
      best.X = nearest;
      best.distance = d;
      best.metric = metric;
    }
  }
  return best;
}

constexpr double kBodyNormalClusterAngleDeg = 15.0;

inline bool sameNormalCluster(const Tddd& a,
                              const Tddd& b,
                              double angle_deg = kBodyNormalClusterAngleDeg);

struct ReferenceProjection {
  bool ok = false;
  Tddd X = {0., 0., 0.};
  double distance = std::numeric_limits<double>::infinity();
  const PreModificationSnapshot::FaceData* face = nullptr;
  Tdd param = {0., 0.};
  bool used_quadratic_fit = false;
  GeometryFallbackLevel fallback = GeometryFallbackLevel::none;
};

struct LocalDirichletProjection {
  bool ok = false;
  Tddd X = {0., 0., 0.};
  double distance = std::numeric_limits<double>::infinity();
  networkFace* face = nullptr;
  Tdd param = {0., 0.};
};

template <class Entity>
inline LocalDirichletProjection projectToLocalDirichletSurface(Entity* entity,
                                                               const Tddd& X) {
  LocalDirichletProjection out;
  if (!entity || !isFinite(X))
    return out;

  for (auto* f : entity->getBoundaryFaces()) {
    if (!f || getNodeFaceBoundaryType(entity, f) != NodeFaceBoundaryType::Dirichlet)
      continue;
    Tdd param = {0., 0.};
    const Tddd nearest = NearestOnDirichletFace(
        X, f, &param,
        [](const auto* q) -> Tddd { return q ? q->getPosition() : Tddd{0., 0., 0.}; });
    const double distance = Norm(nearest - X);
    if (!isFinite(nearest) || !std::isfinite(distance))
      continue;
    if (distance < out.distance) {
      out.ok = true;
      out.X = nearest;
      out.distance = distance;
      out.face = f;
      out.param = param;
    }
  }
  return out;
}

struct ReferenceNormalConstraint {
  double distance = 0.0;
  Tddd direction = {0., 0., 0.};
  double weight = 1.0;
};

inline double normalConstraintDistanceWeight(double gap, double scale) {
  if (!std::isfinite(gap) || !std::isfinite(scale) || scale <= 0.0)
    return 1.0;
  const double h = std::max(0.5 * scale, 1e-12);
  const double r = std::max(0.0, gap) / h;
  return 1.0 / (1.0 + r * r);
}

inline Tddd snapshotProjectionNormal(const PreModificationSnapshot::FaceData& face,
                                     const Tdd& param,
                                     GeometryProjectionMode mode) {
  Tddd n = {0., 0., 0.};
  if (mode == GeometryProjectionMode::LinearFace) {
    const auto& v = face.geometry.vertices_X;
    n = Cross(v[1] - v[0], v[2] - v[0]);
  } else {
    const Tddd du = evalSnapshotQuadDeriv(face, param[0], param[1], 1, 0);
    const Tddd dv = evalSnapshotQuadDeriv(face, param[0], param[1], 0, 1);
    n = Cross(du, dv);
    if (!isFinite(n) || Norm(n) <= 1e-14) {
      const auto& v = face.geometry.vertices_X;
      n = Cross(v[1] - v[0], v[2] - v[0]);
    }
  }
  const double nn = Norm(n);
  return (isFinite(n) && nn > 1e-14) ? n / nn : Tddd{0., 0., 0.};
}

inline std::vector<const PreModificationSnapshot::FaceData*> snapshotCandidateFaces(const SnapshotReferenceState& reference,
                                                                                    const Tddd& X) {
  std::vector<const PreModificationSnapshot::FaceData*> candidates;
  if (reference.index && reference.index->ready) {
    const auto bucket_hits = reference.index->buckets.getData(X, reference.index->query_range);
    candidates.reserve(bucket_hits.size());
    for (const auto* face : bucket_hits)
      if (face)
        candidates.push_back(face);
  }
  if (candidates.empty()) {
    candidates.reserve(reference.snap.faces.size());
    for (const auto& face : reference.snap.faces)
      candidates.push_back(&face);
  }
  return candidates;
}

inline ReferenceProjection projectToSnapshotReference(const SnapshotReferenceState& reference,
                                                      const Tddd& X,
                                                      GeometryProjectionMode mode,
                                                      double trust_radius_factor) {
  ReferenceProjection out;
  if (reference.snap.faces.empty())
    return out;

  if (mode == GeometryProjectionMode::LinearFace) {
    auto projection = findNearestSnapshotFace(X, reference.snap, NodeRelocationSurface::linear, reference.index.get());
    if (!projection.face)
      return out;
    out.ok = true;
    out.X = projection.nearest;
    out.distance = projection.distance;
    out.face = projection.face;
    out.param = projection.param;
    return out;
  }

  if (mode == GeometryProjectionMode::TrueQuadraticFace) {
    auto projection = findNearestSnapshotFace(X, reference.snap, NodeRelocationSurface::true_quadratic, reference.index.get());
    if (!projection.face)
      return out;
    out.ok = true;
    out.X = projection.nearest;
    out.distance = projection.distance;
    out.face = projection.face;
    out.param = projection.param;
    return out;
  }

  auto linear = findNearestSnapshotFace(X, reference.snap, NodeRelocationSurface::linear, reference.index.get());
  if (!linear.face)
    return out;
  const auto* face = linear.face;
  const auto X6 = snapshotFaceNodes(*face);
  if (!std::all_of(X6.begin(), X6.end(), [](const Tddd& q) { return isFinite(q); })) {
    auto fallback = projectToSnapshotReference(reference, X, GeometryProjectionMode::LinearFace, trust_radius_factor);
    fallback.fallback = GeometryFallbackLevel::linear_face;
    return fallback;
  }

  const Tddd fitted = evalSnapshotQuad(*face, linear.param[0], linear.param[1]);
  const double dist = Norm(fitted - X);
  const double trust = trust_radius_factor * std::max(snapshotFaceMaxEdge(*face), 1e-12);
  if (!isFinite(fitted) || !std::isfinite(dist) || dist > trust) {
    auto fallback = projectToSnapshotReference(reference, X, GeometryProjectionMode::LinearFace, trust_radius_factor);
    fallback.fallback = GeometryFallbackLevel::linear_face;
    return fallback;
  }

  out.ok = true;
  out.X = fitted;
  out.distance = dist;
  out.face = face;
  out.param = linear.param;
  out.used_quadratic_fit = true;
  return out;
}

template <class Entity>
inline ReferenceProjection projectToSnapshotReferenceFeatureAware(const SnapshotReferenceState& reference,
                                                                  const Tddd& X,
                                                                  GeometryProjectionMode mode,
                                                                  double trust_radius_factor,
                                                                  const Entity* entity) {
  auto nearest = projectToSnapshotReference(reference, X, mode, trust_radius_factor);
  if (!nearest.ok || !entity)
    return nearest;

  const double identity_tol = std::max(localLength(entity) * 1e-10, 1e-12);
  if (!hasSharpGeometry(entity) || nearest.distance <= identity_tol)
    return nearest;

  std::vector<ReferenceNormalConstraint> constraints;
  const auto candidates = snapshotCandidateFaces(reference, X);
  constraints.reserve(candidates.size());
  for (const auto* face : candidates) {
    if (!face)
      continue;
    const NodeRelocationSurface scheme = (mode == GeometryProjectionMode::LinearFace)
                                             ? NodeRelocationSurface::linear
                                             : NodeRelocationSurface::true_quadratic;
    auto projection = projectToSnapshotFace(X, *face, scheme);
    if (!projection.face || !std::isfinite(projection.distance))
      continue;
    const double trust = trust_radius_factor * std::max(snapshotFaceMaxEdge(*face), 1e-12);
    if (projection.distance > trust)
      continue;
    const Tddd n = snapshotProjectionNormal(*face, projection.param, mode);
    if (!isFinite(n) || Norm(n) <= 0.0)
      continue;
    const double dist = Dot(projection.nearest - X, n);
    const double weight = normalConstraintDistanceWeight(projection.distance, trust);

    bool merged = false;
    for (auto& existing : constraints) {
      if (!sameNormalCluster(existing.direction, n))
        continue;
      if (weight > existing.weight ||
          (weight == existing.weight && std::abs(dist) < std::abs(existing.distance))) {
        existing.distance = dist;
        existing.direction = n;
        existing.weight = weight;
      }
      merged = true;
      break;
    }
    if (!merged)
      constraints.push_back({dist, n, weight});
  }

  if (constraints.empty())
    return nearest;

  std::vector<double> distances;
  std::vector<Tddd> directions;
  std::vector<double> weights;
  distances.reserve(constraints.size());
  directions.reserve(constraints.size());
  weights.reserve(constraints.size());
  for (const auto& c : constraints) {
    if (!std::isfinite(c.distance) || !isFinite(c.direction) || Norm(c.direction) <= 0.0 ||
        !std::isfinite(c.weight) || c.weight <= 0.0)
      continue;
    distances.push_back(c.distance);
    directions.push_back(Normalize(c.direction));
    weights.push_back(c.weight);
  }
  if (distances.empty())
    return nearest;

  const Tddd delta = optimalVectorRegularized(distances, directions, Tddd{0., 0., 0.}, weights);
  if (!isFinite(delta))
    return nearest;
  auto constrained = projectToSnapshotReference(reference, X + delta, mode, trust_radius_factor);
  if (!constrained.ok)
    return nearest;
  constrained.distance = Norm(constrained.X - X);
  return constrained;
}

inline ReferenceProjection projectToAnalyticReference(const InitialConditionReferenceState& reference,
                                                      const Network* water,
                                                      const Tddd& X) {
  ReferenceProjection out;
  if (!water || !water->ic_eta)
    return out;
  Tddd target = X;
  target[2] = water->ic_eta(X, reference.simulation_time);
  if (!isFinite(target))
    return out;
  out.ok = true;
  out.X = target;
  out.distance = Norm(target - X);
  return out;
}

inline ReferenceProjection projectToReference(const ReferenceState& reference,
                                              const Network* water,
                                              const Tddd& X,
                                              GeometryProjectionMode mode,
                                              double trust_radius_factor) {
  return std::visit([&](const auto& concrete_reference) -> ReferenceProjection {
    using ConcreteReference = std::decay_t<decltype(concrete_reference)>;
    if constexpr (std::is_same_v<ConcreteReference, SnapshotReferenceState>)
      return projectToSnapshotReference(concrete_reference, X, mode, trust_radius_factor);
    else
      return projectToAnalyticReference(concrete_reference, water, X);
  },
                    reference);
}

template <class Entity>
inline ReferenceProjection projectToReferenceForEntity(const ReferenceState& reference,
                                                       const Network* water,
                                                       const Tddd& X,
                                                       GeometryProjectionMode mode,
                                                       double trust_radius_factor,
                                                       const Entity* entity) {
  return std::visit([&](const auto& concrete_reference) -> ReferenceProjection {
    using ConcreteReference = std::decay_t<decltype(concrete_reference)>;
    if constexpr (std::is_same_v<ConcreteReference, SnapshotReferenceState>)
      return projectToSnapshotReferenceFeatureAware(concrete_reference, X, mode,
                                                    trust_radius_factor, entity);
    else
      return projectToAnalyticReference(concrete_reference, water, X);
  },
                    reference);
}

struct BodyProjection {
  bool ok = false;
  Tddd X = {0., 0., 0.};
  double distance = std::numeric_limits<double>::infinity();
  networkFace* face = nullptr;
  GeometryFallbackLevel fallback = GeometryFallbackLevel::none;
};

struct BodyNormalConstraint {
  double distance = 0.0;
  Tddd direction = {0., 0., 0.};
  T3Tddd triangle = {};
  networkFace* face = nullptr;
  GeometryFallbackLevel fallback = GeometryFallbackLevel::none;
  double weight = 1.0;
};

template <class Entity>
inline BodyProjection nearestOnBodyFaces(const std::vector<networkFace*>& faces,
                                         const Tddd& X,
                                         GeometryFallbackLevel fallback,
                                         const GeometryProjectorOptions* options,
                                         const Entity* entity) {
  BodyProjection out;
  auto nearest = nearestFaceByDistance(faces, X, options, entity);
  if (!nearest.face || !std::isfinite(nearest.distance))
    return out;
  out.ok = true;
  out.face = nearest.face;
  out.distance = nearest.distance;
  out.X = nearest.X;
  out.fallback = fallback;
  return out;
}

inline BodyProjection nearestOnBodyFaces(const std::vector<networkFace*>& faces,
                                         const Tddd& X,
                                         GeometryFallbackLevel fallback) {
  return nearestOnBodyFaces(faces, X, fallback,
                            static_cast<const GeometryProjectorOptions*>(nullptr),
                            static_cast<const ContactDetectable*>(nullptr));
}

inline std::vector<networkFace*> expandBodyFacesOneRing(const std::vector<networkFace*>& seeds,
                                                        int depth) {
  std::vector<networkFace*> current = uniqueFaces(seeds);
  if (depth <= 0 || current.empty())
    return current;

  std::unordered_set<networkFace*> seen(current.begin(), current.end());
  std::vector<networkFace*> frontier = current;
  for (int d = 0; d < depth; ++d) {
    std::vector<networkFace*> next;
    for (auto* f : frontier) {
      if (!f)
        continue;
      for (auto* l : f->getLines()) {
        if (!l)
          continue;
        for (auto* adj : l->getBoundaryFaces()) {
          if (!adj || seen.contains(adj))
            continue;
          seen.insert(adj);
          next.push_back(adj);
          current.push_back(adj);
        }
      }
    }
    frontier = std::move(next);
    if (frontier.empty())
      break;
  }
  return current;
}

inline std::vector<networkFace*> expandBodyFacesAcrossNonSharp(networkFace* seed,
                                                               int depth,
                                                               double feature_angle_rad) {
  std::vector<networkFace*> current;
  if (!seed)
    return current;
  current.push_back(seed);
  if (depth <= 0)
    return current;

  std::unordered_set<networkFace*> seen;
  seen.insert(seed);
  std::vector<networkFace*> frontier{seed};
  for (int d = 0; d < depth; ++d) {
    std::vector<networkFace*> next;
    for (auto* f : frontier) {
      if (!f)
        continue;
      for (auto* l : f->getLines()) {
        if (!l || l->SharpQ(feature_angle_rad))
          continue;
        for (auto* adj : l->getBoundaryFaces()) {
          if (!adj || seen.contains(adj))
            continue;
          seen.insert(adj);
          next.push_back(adj);
          current.push_back(adj);
        }
      }
    }
    frontier = std::move(next);
    if (frontier.empty())
      break;
  }
  return current;
}

inline surface_geometry::PrincipalCurvatureResult principalCurvatureFromParametricDerivatives(
    const Tddd& Xu,
    const Tddd& Xv,
    const Tddd& Xuu,
    const Tddd& Xuv,
    const Tddd& Xvv) {
  surface_geometry::PrincipalCurvatureResult result;
  if (!isFinite(Xu) || !isFinite(Xv) || !isFinite(Xuu) || !isFinite(Xuv) || !isFinite(Xvv))
    return result;

  Tddd n = Cross(Xu, Xv);
  const double n_norm = Norm(n);
  if (!(n_norm > 1e-14) || !std::isfinite(n_norm))
    return result;
  n /= n_norm;

  const double E = Dot(Xu, Xu);
  const double F = Dot(Xu, Xv);
  const double G = Dot(Xv, Xv);
  const double det_I = E * G - F * F;
  if (!(det_I > 1e-20) || !std::isfinite(det_I))
    return result;

  const double L = Dot(Xuu, n);
  const double M = Dot(Xuv, n);
  const double N = Dot(Xvv, n);

  const double s00 = (L * G - M * F) / det_I;
  const double s01 = (G * M - F * N) / det_I;
  const double s10 = (M * E - L * F) / det_I;
  const double s11 = (N * E - M * F) / det_I;
  if (!std::isfinite(s00) || !std::isfinite(s01) ||
      !std::isfinite(s10) || !std::isfinite(s11))
    return result;

  const double trace = s00 + s11;
  const double det_S = s00 * s11 - s01 * s10;
  const double disc = std::sqrt(std::max(trace * trace - 4.0 * det_S, 0.0));
  const double lambda_a = 0.5 * (trace + disc);
  const double lambda_b = 0.5 * (trace - disc);

  auto direction_for = [&](double lambda) -> Tddd {
    double ev_u = lambda - s11;
    double ev_v = s10;
    double ev_len = std::sqrt(ev_u * ev_u + ev_v * ev_v);
    if (!(ev_len > 1e-12)) {
      ev_u = s01;
      ev_v = lambda - s00;
      ev_len = std::sqrt(ev_u * ev_u + ev_v * ev_v);
    }
    Tddd dir = (ev_len > 1e-12) ? (ev_u * Xu + ev_v * Xv) / ev_len : Xu;
    dir -= Dot(dir, n) * n;
    const double dir_norm = Norm(dir);
    return (dir_norm > 1e-14 && isFinite(dir)) ? dir / dir_norm : Tddd{0., 0., 0.};
  };

  Tddd pd_a = direction_for(lambda_a);
  Tddd pd_b = direction_for(lambda_b);
  if (!isFinite(pd_a) || Norm(pd_a) <= 1e-14) {
    pd_a = Xu - Dot(Xu, n) * n;
    const double pd_norm = Norm(pd_a);
    if (!(pd_norm > 1e-14) || !isFinite(pd_a))
      return result;
    pd_a /= pd_norm;
  }
  if (!isFinite(pd_b) || Norm(pd_b) <= 1e-14) {
    pd_b = Cross(n, pd_a);
    const double pd_norm = Norm(pd_b);
    if (!(pd_norm > 1e-14) || !isFinite(pd_b))
      return result;
    pd_b /= pd_norm;
  }

  if (std::abs(lambda_a) >= std::abs(lambda_b)) {
    result.k1 = lambda_a;
    result.k2 = lambda_b;
    result.PD1 = pd_a;
    result.PD2 = pd_b;
  } else {
    result.k1 = lambda_b;
    result.k2 = lambda_a;
    result.PD1 = pd_b;
    result.PD2 = pd_a;
  }
  result.kmax = std::max(std::abs(result.k1), std::abs(result.k2));
  result.valid = isFinite(result.PD1) && isFinite(result.PD2) &&
                 std::isfinite(result.k1) && std::isfinite(result.k2);
  return result;
}

inline surface_geometry::PrincipalCurvatureResult referenceCurvatureAt(
    const PreModificationSnapshot::FaceData* face,
    const Tdd& param,
    GeometryProjectionMode mode) {
  if (!face || mode == GeometryProjectionMode::LinearFace)
    return {};
  const Tddd Xu = evalSnapshotQuadDeriv(*face, param[0], param[1], 1, 0);
  const Tddd Xv = evalSnapshotQuadDeriv(*face, param[0], param[1], 0, 1);
  const Tddd Xuu = evalSnapshotQuadDeriv(*face, param[0], param[1], 2, 0);
  const Tddd Xuv = evalSnapshotQuadDeriv(*face, param[0], param[1], 1, 1);
  const Tddd Xvv = evalSnapshotQuadDeriv(*face, param[0], param[1], 0, 2);
  return principalCurvatureFromParametricDerivatives(Xu, Xv, Xuu, Xuv, Xvv);
}

inline surface_geometry::PrincipalCurvatureResult bodyCurvatureAt(
    networkFace* body_face,
    const Tddd& target_X,
    double feature_angle_rad,
    int bfs_depth = 2) {
  surface_geometry::PrincipalCurvatureResult result;
  if (!body_face || !isFinite(target_X))
    return result;
  Tddd n_axis = body_face->normal;
  const double n_norm = Norm(n_axis);
  if (!(n_norm > 1e-14) || !isFinite(n_axis))
    return result;
  n_axis /= n_norm;

  auto [p0, p1, p2] = body_face->getPoints();
  std::array<networkPoint*, 3> face_points{p0, p1, p2};
  Tddd u_axis = {0., 0., 0.};
  for (int i = 0; i < 3; ++i) {
    auto* a = face_points[i];
    auto* b = face_points[(i + 1) % 3];
    if (!a || !b)
      continue;
    Tddd candidate = b->X - a->X;
    candidate -= Dot(candidate, n_axis) * n_axis;
    const double len = Norm(candidate);
    if (len > 1e-14 && isFinite(candidate)) {
      u_axis = candidate / len;
      break;
    }
  }
  if (!(Norm(u_axis) > 1e-14) || !isFinite(u_axis)) {
    const Tddd ref = std::abs(n_axis[0]) < 0.9 ? Tddd{1., 0., 0.} : Tddd{0., 1., 0.};
    u_axis = ref - Dot(ref, n_axis) * n_axis;
    const double len = Norm(u_axis);
    if (!(len > 1e-14) || !isFinite(u_axis))
      return result;
    u_axis /= len;
  }
  Tddd v_axis = Cross(n_axis, u_axis);
  const double v_norm = Norm(v_axis);
  if (!(v_norm > 1e-14) || !isFinite(v_axis))
    return result;
  v_axis /= v_norm;

  const auto patch_faces = expandBodyFacesAcrossNonSharp(body_face, bfs_depth, feature_angle_rad);
  std::vector<Tddd> points;
  std::unordered_set<networkPoint*> seen_points;
  points.reserve(3 * patch_faces.size());
  for (auto* f : patch_faces) {
    if (!f)
      continue;
    for (auto* p : f->getPoints()) {
      if (!p || seen_points.contains(p) || !isFinite(p->X))
        continue;
      seen_points.insert(p);
      points.push_back(p->X);
    }
  }
  if (points.size() < 5)
    return result;

  const auto fit = surface_geometry::fitQuadricLocal_unroll(target_X, u_axis, v_axis, n_axis, points);
  if (!fit.valid)
    return result;
  result = surface_geometry::principalCurvaturesFromQuadric(fit);
  if (!result.valid || !std::isfinite(result.kmax))
    return {};
  return result;
}

inline void accumulateThetaFromCurvature(TargetThetaResult& out,
                                         const Tddd& edge,
                                         const surface_geometry::PrincipalCurvatureResult& curv,
                                         bool reference_target) {
  if (!curv.valid || !isFinite(edge) || Norm(edge) <= 1e-14)
    return;
  const double theta = surface_geometry::edgeCurvatureAngle(edge, curv.k1, curv.k2, curv.PD1, curv.PD2);
  if (!std::isfinite(theta) || theta < 0.0)
    return;
  out.valid = true;
  out.theta = std::max(out.theta, theta);
  if (reference_target)
    out.used_reference = true;
  else
    out.used_body = true;
}

template <class Entity>
inline BodyProjection nearestOnBodyObjects(const std::vector<Network*>& body_objects,
                                           const Tddd& X,
                                           double threshold,
                                           const GeometryProjectorOptions& options,
                                           const Entity* entity) {
  if (!useAnisotropicBodyNearest(entity, &options)) {
    BodyProjection best;
    for (auto* body : body_objects) {
      if (!body)
        continue;
      auto [face, nearest] = body->Nearest(X);
      if (!face || !isFinite(nearest))
        continue;
      const double distance = Norm(nearest - X);
      if (std::isfinite(distance) && distance < best.distance) {
        best.ok = true;
        best.face = face;
        best.X = nearest;
        best.distance = distance;
        best.fallback = GeometryFallbackLevel::global_body;
      }
    }
    if (best.ok && std::isfinite(best.distance) && best.distance <= threshold)
      return best;
    return {};
  }

  std::vector<networkFace*> candidates;
  for (auto* body : body_objects) {
    if (!body)
      continue;
    const auto faces = body->getBoundaryFaces();
    candidates.insert(candidates.end(), faces.begin(), faces.end());
  }
  auto projection = nearestOnBodyFaces(candidates, X, GeometryFallbackLevel::global_body, &options, entity);
  if (projection.ok && std::isfinite(projection.distance) && projection.distance <= threshold)
    return projection;
  return {};
}

inline bool sameNormalCluster(const Tddd& a, const Tddd& b, double angle_deg) {
  if (!isFinite(a) || !isFinite(b) || Norm(a) <= 0.0 || Norm(b) <= 0.0)
    return false;
  return Dot(Normalize(a), Normalize(b)) > std::cos(angle_deg * M_PI / 180.0);
}

inline int normalClusterCount(const std::vector<BodyNormalConstraint>& constraints) {
  std::vector<Tddd> representatives;
  for (const auto& c : constraints) {
    if (!isFinite(c.direction) || Norm(c.direction) <= 0.0)
      continue;
    const Tddd n = Normalize(c.direction);
    bool merged = false;
    for (auto& rep : representatives) {
      if (sameNormalCluster(rep, n)) {
        rep = Normalize(rep + n);
        merged = true;
        break;
      }
    }
    if (!merged)
      representatives.push_back(n);
  }
  return static_cast<int>(representatives.size());
}

inline void addNormalConstraintFromNearest(std::vector<BodyNormalConstraint>& constraints,
                                           const Tddd& X,
                                           const T3Tddd& tri,
                                           const Tddd& X_nearest,
                                           const Tddd& n_raw,
                                           networkFace* body_face,
                                           GeometryFallbackLevel fallback,
                                           double weight = 1.0) {
  if (!body_face)
    return;
  if (!isFinite(X_nearest) || !isFinite(n_raw) || Norm(n_raw) <= 0.0)
    return;
  if (!std::isfinite(weight) || weight <= 0.0)
    return;
  const Tddd n = Normalize(n_raw);
  const double dist = Dot(X_nearest - X, n);

  // 法線方向ごとに最も近い拘束だけを残す。flat な領域ではほぼ同じ法線が
  // 一つの代表拘束へまとまり、sharp/contact-interface では異なる法線方向の
  // 拘束が複数残る。これにより、角や喫水線では単純な最近傍面への貼り付きではなく、
  // 複数面の法線拘束を同時に満たす repair を行える。
  for (auto& existing : constraints) {
    if (!sameNormalCluster(existing.direction, n))
      continue;
    if (weight > existing.weight ||
        (weight == existing.weight && std::abs(dist) < std::abs(existing.distance))) {
      existing.distance = dist;
      existing.direction = n;
      existing.triangle = tri;
      existing.face = body_face;
      existing.fallback = fallback;
      existing.weight = weight;
    }
    return;
  }

  constraints.push_back({dist, n, tri, body_face, fallback, weight});
}

inline void addNormalConstraint(std::vector<BodyNormalConstraint>& constraints,
                                const Tddd& X,
                                networkFace* body_face,
                                GeometryFallbackLevel fallback,
                                double weight = 1.0) {
  if (!body_face)
    return;
  const T3Tddd tri = ToX(body_face);
  auto [t0, t1, X_nearest, n_raw] = Nearest_(X, tri);
  (void)t0;
  (void)t1;
  addNormalConstraintFromNearest(constraints, X, tri, X_nearest, n_raw, body_face, fallback, weight);
}

template <class Entity>
inline std::vector<BodyNormalConstraint> directContactNormalConstraints(const Entity* entity,
                                                                        const Tddd& X,
                                                                        GeometryFallbackLevel fallback,
                                                                        const GeometryProjectorOptions& options) {
  std::vector<BodyNormalConstraint> constraints;
  if (!entity)
    return constraints;

  for (auto* fluid_face : entity->getBoundaryFaces()) {
    if (!fluid_face || getNodeFaceBoundaryType(entity, fluid_face) != NodeFaceBoundaryType::Neumann)
      continue;
    const auto* state = entity->findContactState(fluid_face);
    if (!state || state->contact_opponent_faces.empty())
      continue;
    auto projection = nearestOnBodyFaces(state->contact_opponent_faces, X, fallback, &options, entity);
    if (!projection.ok)
      continue;
    const double weight = normalConstraintDistanceWeight(projection.distance, entity->contact_range);
    if (useAnisotropicBodyNearest(entity, &options)) {
      const T3Tddd tri = ToX(projection.face);
      addNormalConstraintFromNearest(constraints, X, tri, projection.X,
                                     projection.face ? projection.face->normal : Tddd{0., 0., 0.},
                                     projection.face, fallback, weight);
    } else {
      addNormalConstraint(constraints, X, projection.face, fallback, weight);
    }
  }
  return constraints;
}

template <class Entity>
inline std::vector<BodyNormalConstraint> contactPatchNormalConstraints(const Entity* entity,
                                                                       const Tddd& X,
                                                                       const std::vector<networkFace*>& body_faces,
                                                                       GeometryFallbackLevel fallback,
                                                                       const GeometryProjectorOptions& options) {
  std::vector<BodyNormalConstraint> constraints;
  if (!entity || body_faces.empty())
    return constraints;

  constexpr double short_range = 0.01;
  const Tddd X_contact_query = entity->getPosition();
  for (auto* fluid_face : entity->getBoundaryFaces()) {
    if (!fluid_face || getNodeFaceBoundaryType(entity, fluid_face) != NodeFaceBoundaryType::Neumann)
      continue;
    for (auto* body_face : body_faces) {
      if (!body_face)
        continue;
      const T3Tddd tri = ToX(body_face);
      if (!isInContact(X_contact_query, fluid_face->normal, tri, entity->contact_range))
        continue;
      Tddd X_nearest = {0., 0., 0.};
      Tddd n_raw = {0., 0., 0.};
      if (useAnisotropicBodyNearest(entity, &options)) {
        X_nearest = bodyFaceNearestPoint(entity, &options, X, body_face);
        n_raw = body_face->normal;
      } else {
        auto [t0, t1, euclidean_nearest, euclidean_normal] = Nearest_(X, tri);
        (void)t0;
        (void)t1;
        X_nearest = euclidean_nearest;
        n_raw = euclidean_normal;
      }
      if (!isFinite(X_nearest) || !isFinite(n_raw) || Norm(n_raw) <= 0.0)
        continue;
      const Tddd to_nearest = X_nearest - X;
      const double dist = Norm(to_nearest);
      const double metric_dist = bodyFaceSelectionDistance(entity, &options, X, body_face, X_nearest);
      const Tddd n = Normalize(n_raw);
      if (((isFlat(n, to_nearest, contactAcceptanceAngle(metric_dist, entity->contact_range)) ||
            isFlat(n, -to_nearest, contactAcceptanceAngle(metric_dist, entity->contact_range))) &&
           entity->contact_range >= metric_dist) ||
          short_range * entity->contact_range >= metric_dist) {
        const double weight = normalConstraintDistanceWeight(dist, entity->contact_range);
        if (useAnisotropicBodyNearest(entity, &options))
          addNormalConstraintFromNearest(constraints, X, tri, X_nearest, n_raw, body_face, fallback, weight);
        else
          addNormalConstraint(constraints, X, body_face, fallback, weight);
      }
    }
  }
  return constraints;
}

template <class Entity>
inline bool usesSharpNormalConstraints(const Entity* entity,
                                       const std::vector<BodyNormalConstraint>& constraints) {
  if (!entity)
    return false;
  if (isBCInterfaceEntity(entity))
    return true;
  if (hasSharpGeometry(entity))
    return true;
  return normalClusterCount(constraints) > 1;
}

inline BodyProjection projectFromNormalConstraints(const Tddd& X,
                                                   const std::vector<BodyNormalConstraint>& constraints,
                                                   GeometryFallbackLevel fallback) {
  BodyProjection out;
  if (constraints.empty())
    return out;

  std::vector<double> distances;
  std::vector<Tddd> directions;
  std::vector<double> weights;
  std::vector<T3Tddd> triangles;
  distances.reserve(constraints.size());
  directions.reserve(constraints.size());
  weights.reserve(constraints.size());
  triangles.reserve(constraints.size());
  for (const auto& c : constraints) {
    if (!isFinite(c.direction) || Norm(c.direction) <= 0.0 || !std::isfinite(c.distance) ||
        !std::isfinite(c.weight) || c.weight <= 0.0)
      continue;
    distances.push_back(c.distance);
    directions.push_back(Normalize(c.direction));
    weights.push_back(c.weight);
    triangles.push_back(c.triangle);
  }
  if (distances.empty())
    return out;

  const Tddd delta = optimalVectorRegularized(distances, directions, Tddd{0., 0., 0.}, weights);
  if (!isFinite(delta))
    return out;
  Tddd target = X + delta;
  if (!triangles.empty())
    target = Nearest(target, triangles);
  if (!isFinite(target))
    return out;

  networkFace* best_face = nullptr;
  double best_residual = std::numeric_limits<double>::infinity();
  for (const auto& c : constraints) {
    if (!c.face)
      continue;
    const double residual = Norm(Nearest(target, c.face) - target);
    if (std::isfinite(residual) && residual < best_residual) {
      best_residual = residual;
      best_face = c.face;
    }
  }

  out.ok = true;
  out.X = target;
  out.distance = Norm(target - X);
  out.face = best_face;
  out.fallback = fallback;
  return out;
}

template <class Entity>
inline BodyProjection projectBodyFromCandidates(const Entity* entity,
                                                const Tddd& X,
                                                const std::vector<networkFace*>& body_faces,
                                                std::vector<BodyNormalConstraint> constraints,
                                                GeometryFallbackLevel fallback,
                                                double threshold,
                                                const GeometryProjectorOptions& options,
                                                bool allow_nearest_without_constraints = false) {
  if (!entity || body_faces.empty())
    return {};

  if (constraints.empty())
    constraints = contactPatchNormalConstraints(entity, X, body_faces, fallback, options);

  const bool sharp_normal_constraints = usesSharpNormalConstraints(entity, constraints);
  if (sharp_normal_constraints && !constraints.empty()) {
    auto projection = projectFromNormalConstraints(X, constraints, fallback);
    if (projection.ok && std::isfinite(projection.distance) && projection.distance <= threshold)
      return projection;
    return {};
  }
  if (sharp_normal_constraints && constraints.empty() && !allow_nearest_without_constraints)
    return {};

  auto projection = nearestOnBodyFaces(body_faces, X, fallback, &options, entity);
  if (projection.ok && projection.distance <= threshold)
    return projection;
  return {};
}

template <class Entity>
inline BodyProjection projectToBodyCandidates(const Entity* entity,
                                              const Tddd& X,
                                              double threshold,
                                              const std::vector<Network*>& body_objects,
                                              const GeometryProjectorOptions& options) {
  if (!entity)
    return {};

  const auto line_contacts = uniqueFaces(getEffectiveContactFaces(entity));
  auto direct_constraints = directContactNormalConstraints(entity, X, GeometryFallbackLevel::none, options);
  auto line_projection = projectBodyFromCandidates(entity, X, line_contacts, direct_constraints,
                                                   GeometryFallbackLevel::none, threshold, options);
  if (line_projection.ok)
    return line_projection;

  // Repair projection is contact-local first.  If direct contact fails, expand
  // only from known local contact seeds before considering a tightly limited
  // global nearest fallback.
  std::vector<networkFace*> endpoint_contacts;
  if constexpr (std::is_same_v<Entity, networkLine>) {
    endpoint_contacts = endpointContactFaces(entity);
    auto endpoint_projection = projectBodyFromCandidates(entity, X, endpoint_contacts, {},
                                                         GeometryFallbackLevel::endpoint_contacts, threshold, options,
                                                         /*allow_nearest_without_constraints=*/true);
    if (endpoint_projection.ok)
      return endpoint_projection;
  }

  std::vector<networkFace*> seed_contacts = line_contacts;
  seed_contacts.insert(seed_contacts.end(), endpoint_contacts.begin(), endpoint_contacts.end());
  seed_contacts = uniqueFaces(seed_contacts);
  if (options.body_bfs_fallback_depth > 0 && !seed_contacts.empty()) {
    auto bfs_contacts = expandBodyFacesOneRing(seed_contacts, options.body_bfs_fallback_depth);
    auto bfs_projection = projectBodyFromCandidates(entity, X, bfs_contacts, {},
                                                    GeometryFallbackLevel::bfs_contacts, threshold, options,
                                                    /*allow_nearest_without_constraints=*/true);
    if (bfs_projection.ok)
      return bfs_projection;
  }

  const double contact_range = entity ? entity->contact_range : 0.0;
  if (options.enable_global_body_fallback &&
      std::isfinite(contact_range) && contact_range > 0.0) {
    const double global_threshold = std::min(
        threshold,
        std::max(0.0, options.global_body_fallback_range_factor) * contact_range);
    auto global_projection = nearestOnBodyObjects(body_objects, X, global_threshold, options, entity);
    if (global_projection.ok)
      return global_projection;
  }

  return {};
}

template <class Entity>
inline double bodyGapThreshold(const Entity* entity,
                               double shift_limit,
                               double max_body_gap_factor,
                               double contact_range_factor) {
  const double by_shift = max_body_gap_factor * std::max(shift_limit, 1e-12);
  const double contact_range = entity ? entity->contact_range : 0.0;
  if (std::isfinite(contact_range) && contact_range > 0.0)
    return std::min(by_shift, contact_range_factor * contact_range);
  return by_shift;
}

inline void applyClamp(GeometryTarget& target, const Tddd& X_in) {
  const Tddd delta = target.target_X - X_in;
  const double delta_norm = Norm(delta);
  if (std::isfinite(delta_norm) && delta_norm > target.shift_limit && target.shift_limit > 0.0)
    target.delta_clamped = target.shift_limit * Normalize(delta);
  else
    target.delta_clamped = delta;
  target.target_X_clamped = X_in + target.delta_clamped;
  target.move_ratio = (target.shift_limit > 0.0) ? Norm(target.delta_clamped) / target.shift_limit : 0.0;
}

} // namespace GeometryProjectorDetail

struct BoundaryGeometryTargetFlags {
  bool need_reference = false;
  bool need_body = false;
};

/*
 * BoundaryGeometryTargetProvider は、境界 DOF が従うべき幾何 target を決める
 * source of truth である。
 *
 * Dirichlet DOF は reference water/free surface を target にする。
 * Neumann DOF は接触している構造物 body surface を target にする。
 * BCInterface DOF は水面と構造物の接合線なので、reference と body の両方を
 * target にする。
 *
 * repair、waterline repair、remesh trial の検証、HD check は、必ずこの
 * provider を通して target を問い合わせる。repair と HD check が別々の面を
 * 正しい target と見なすと、正しく repair された点が HD で reject されたり、
 * 逆に間違った面へ repair された点が HD を通過したりする。
 *
 * repair では、Dirichlet-only entity は reference surface へ射影し、Neumann-only
 * entity は body surface へ射影する。BCInterface entity は、古い node relocation
 * と同じく現在の局所 Dirichlet 面へ一度貼り付け、その位置から body surface へ
 * 射影する Dirichlet-first 候補を優先する。局所 Dirichlet 面が使えない場合は、
 * reference へ射影した位置を起点に body surface へ射影する反復に fallback する。
 * 実際に点を動かすときは、ここで返す target_X をそのまま使うのではなく、
 * shift_limit で制限された target_X_clamped / delta_clamped を使う。つまり、この
 * provider は「どの面へ寄せるか」だけでなく、「今回の repair で許容される移動量内
 * ではどこまで寄せるか」も一貫して返す。
 *
 * body surface への repair 方法は、接触面の法線構造で変わる。ほぼ同じ法線方向の
 * face だけに接している flat な接触では、代表的な最近傍 body face へ射影する。
 * 一方、SharpQ な entity、BCInterface entity、または複数の法線クラスタに接している
 * entity では、複数の独立した法線拘束を同時に満たす移動量を解く。角や喫水線では
 * 「単に一番近い面へ貼る」と別の面へ吸い込まれることがあるため、sharp/contact
 * interface では複数法線拘束を優先する。
 *
 * body target は保守的に探す。まず entity 自身の contact face を使い、line の
 * 場合は endpoint の contact face も見る。まだ見つからなければ既知の local
 * contact seed から BFS=1 だけ広げる。最後に、contact_range 内にある場合に限り
 * global nearest body face を fallback として使う。この fallback policy はここに
 * 閉じ込め、呼び出し側がより広い探索や別ルールを勝手に実装しないようにする。
 *
 * このクラスは topology 操作、trial score、patch replace、補間、BVP の境界値更新
 * には責任を持たない。このクラスの責務は一つだけで、「この境界 entity/sample は
 * どの幾何面と比較され、どの幾何面へ repair されるべきか」を答えることである。
 */
struct BoundaryGeometryTargetProvider {
  const ReferenceState& reference;
  const std::vector<Network*>& body_objects;
  GeometryProjectorOptions options;
  std::uint64_t provider_epoch = 0;

  GeometryTarget queryEntityTarget(networkPoint* p,
                                   const Network* water,
                                   const Tddd& X_in) const;
  GeometryTarget queryEntityTarget(networkLine* l,
                                   const Network* water,
                                   const Tddd& X_in) const;
  GeometryTarget queryLineSampleTarget(networkLine* l,
                                       const Network* water,
                                       double s,
                                       const Tddd& X_seed) const;
  bool ensureProviderTargetCurvature(networkPoint* p,
                                     const Network* water) const;
  Tddd curvatureSeedForLineSample(networkLine* l,
                                  const Network* water,
                                  double s,
                                  const Tddd& X_seed,
                                  bool* used_curvature_seed = nullptr) const;
  ProviderFieldSample sampleReferenceField(const GeometryTarget& target,
                                           ProviderFieldKind kind) const;
  TargetThetaResult queryTargetTheta(networkPoint* p,
                                     const Network* water,
                                     const Tddd& X_in,
                                     const Tddd& edge) const;
  TargetThetaResult queryTargetTheta(networkLine* l,
                                     const Network* water,
                                     const Tddd& X_in,
                                     const Tddd& edge) const;

  GeometryProjectorDetail::ReferenceProjection queryReferenceDistance(const Network* water,
                                                                      const Tddd& X) const {
    return GeometryProjectorDetail::projectToReference(reference, water, X,
                                                       options.mode,
                                                       options.trust_radius_factor);
  }

  template <class Entity>
  GeometryProjectorDetail::BodyProjection queryBodyDistanceForEntity(const Entity* entity,
                                                                     const Tddd& X) const {
    if (!entity)
      return {};
    const double shift_limit = currentShiftLimit(entity, options.move_limit_factor);
    if (!std::isfinite(shift_limit) || shift_limit <= 1e-12)
      return {};
    const double threshold = GeometryProjectorDetail::bodyGapThreshold(entity, shift_limit,
                                                                       options.max_body_gap_factor,
                                                                       options.contact_range_factor);
    return GeometryProjectorDetail::projectToBodyCandidates(entity, X, threshold,
                                                            body_objects, options);
  }

  template <class Entity>
  BoundaryGeometryTargetFlags targetFlags(const Entity* entity) const {
    if (!entity)
      return {};
    const bool is_interface = GeometryProjectorDetail::isBCInterfaceEntity(entity);
    return {is_interface || hasAnyDirichletBoundaryState(entity),
            is_interface || hasAnyNeumannBoundaryState(entity)};
  }

  BoundaryGeometryTargetFlags targetFlags(networkFace* f) const {
    BoundaryGeometryTargetFlags flags;
    if (!f)
      return flags;
    auto visit_entity = [&](auto* entity) {
      if (!entity)
        return;
      // Face quadrature samples must follow the entity-face boundary type, not
      // the summary BCInterface flag.  A Dirichlet face that touches a waterline
      // edge should not be forced to satisfy the body target over the whole face;
      // BCInterface point/line constraints are checked separately via
      // targetFlags(entity) and queryEntityTarget(...).
      const auto bc = getNodeFaceBoundaryType(entity, f);
      if (bc == NodeFaceBoundaryType::Dirichlet)
        flags.need_reference = true;
      else if (bc == NodeFaceBoundaryType::Neumann)
        flags.need_body = true;
    };
    auto [p0, p1, p2] = f->getPoints();
    auto [l0, l1, l2] = f->getLines();
    visit_entity(p0);
    visit_entity(p1);
    visit_entity(p2);
    visit_entity(l0);
    visit_entity(l1);
    visit_entity(l2);
    return flags;
  }
};

inline const char* toString(GeometryTargetStatus status) {
  switch (status) {
  case GeometryTargetStatus::ok:
    return "ok";
  case GeometryTargetStatus::ok_via_fallback:
    return "ok_via_fallback";
  case GeometryTargetStatus::no_convergence:
    return "no_convergence";
  case GeometryTargetStatus::not_applicable:
    return "not_applicable";
  case GeometryTargetStatus::no_reference_surface:
    return "no_reference_surface";
  case GeometryTargetStatus::no_body_surface:
    return "no_body_surface";
  case GeometryTargetStatus::body_gap_too_large:
    return "body_gap_too_large";
  case GeometryTargetStatus::no_move_possible:
    return "no_move_possible";
  case GeometryTargetStatus::invalid_projection:
    return "invalid_projection";
  }
  return "unknown";
}

inline const char* toString(GeometryFallbackLevel fallback) {
  switch (fallback) {
  case GeometryFallbackLevel::none:
    return "none";
  case GeometryFallbackLevel::linear_face:
    return "linear_face";
  case GeometryFallbackLevel::endpoint_contacts:
    return "endpoint_contacts";
  case GeometryFallbackLevel::bfs_contacts:
    return "bfs_contacts";
  case GeometryFallbackLevel::global_body:
    return "global_body";
  }
  return "unknown";
}

template <class Entity>
inline GeometryTarget queryGeometryTargetForEntity(const BoundaryGeometryTargetProvider& provider,
                                                   Entity* entity,
                                                   const Network* water,
                                                   const Tddd& X_in) {
  using namespace GeometryProjectorDetail;
  GeometryTarget target;
  if (!entity) {
    target.status = GeometryTargetStatus::not_applicable;
    return target;
  }

  const bool is_interface = isBCInterfaceEntity(entity);
  const bool has_dirichlet = hasAnyDirichletBoundaryState(entity) || is_interface;
  const bool has_neumann = hasAnyNeumannBoundaryState(entity) || is_interface;
  if (!has_dirichlet && !has_neumann) {
    target.status = GeometryTargetStatus::not_applicable;
    return target;
  }

  target.shift_limit = currentShiftLimit(entity, provider.options.move_limit_factor);
  if (!std::isfinite(target.shift_limit) || target.shift_limit <= 1e-12) {
    target.status = GeometryTargetStatus::no_move_possible;
    return target;
  }
  target.body_gap_threshold = bodyGapThreshold(entity, target.shift_limit,
                                               provider.options.max_body_gap_factor,
                                               provider.options.contact_range_factor);

  if (has_dirichlet && !has_neumann) {
    const auto projection = projectToReferenceForEntity(provider.reference, water, X_in,
                                                        provider.options.mode,
                                                        provider.options.trust_radius_factor,
                                                        entity);
    if (!projection.ok) {
      target.status = GeometryTargetStatus::no_reference_surface;
      return target;
    }
    target.status = projection.fallback == GeometryFallbackLevel::none
                        ? GeometryTargetStatus::ok
                        : GeometryTargetStatus::ok_via_fallback;
    target.target_X = projection.X;
    target.reference_face = projection.face;
    target.param_uv = projection.param;
    target.dirichlet_gap = projection.distance;
    target.fallback = projection.fallback;
    target.used_quadratic_fit = projection.used_quadratic_fit;
    applyClamp(target, X_in);
    if (!isFinite(target.target_X_clamped) || !std::isfinite(target.move_ratio))
      target.status = GeometryTargetStatus::invalid_projection;
    return target;
  }

  if (!has_dirichlet && has_neumann) {
    const auto body = projectToBodyCandidates(entity, X_in, target.body_gap_threshold,
                                              provider.body_objects, provider.options);
    if (!body.ok) {
      target.status = GeometryTargetStatus::no_body_surface;
      return target;
    }
    target.status = body.fallback == GeometryFallbackLevel::none
                        ? GeometryTargetStatus::ok
                        : GeometryTargetStatus::ok_via_fallback;
    target.target_X = body.X;
    target.body_face = body.face;
    target.body_gap = body.distance;
    target.fallback = body.fallback;
    applyClamp(target, X_in);
    if (!isFinite(target.target_X_clamped) || !std::isfinite(target.move_ratio))
      target.status = GeometryTargetStatus::invalid_projection;
    return target;
  }

  const double tol = std::max(provider.options.tol_relative * localLength(entity), 1e-12);
  const int max_iter = std::max(1, provider.options.max_iter);

  if (provider.options.bcinterface_local_dirichlet_first && is_interface) {
    Tddd X = X_in;
    double best_score = std::numeric_limits<double>::infinity();
    GeometryTarget best = target;
    bool found = false;
    bool converged = false;

    for (int iter = 0; iter < max_iter; ++iter) {
      const auto local_dirichlet = projectToLocalDirichletSurface(entity, X);
      if (!local_dirichlet.ok)
        break;
      const auto body = projectToBodyCandidates(entity, local_dirichlet.X,
                                                target.body_gap_threshold,
                                                provider.body_objects, provider.options);
      if (!body.ok)
        break;
      const auto local_dirichlet_after_body = projectToLocalDirichletSurface(entity, body.X);
      const double final_dirichlet_gap = local_dirichlet_after_body.ok
                                             ? local_dirichlet_after_body.distance
                                             : local_dirichlet.distance;
      const Tddd reference_query_X = local_dirichlet_after_body.ok
                                         ? local_dirichlet_after_body.X
                                         : local_dirichlet.X;
      const auto reference_meta = projectToReferenceForEntity(provider.reference, water,
                                                              reference_query_X,
                                                              provider.options.mode,
                                                              provider.options.trust_radius_factor,
                                                              entity);
      if (!reference_meta.ok || body.distance > target.body_gap_threshold)
        break;

      const double score = final_dirichlet_gap * final_dirichlet_gap +
                           body.distance * body.distance;
      if (std::isfinite(score) && score < best_score) {
        best_score = score;
        found = true;
        best = target;
        best.target_X = body.X;
        best.reference_face = reference_meta.face;
        best.param_uv = reference_meta.param;
        best.body_face = body.face;
        best.dirichlet_gap = final_dirichlet_gap;
        best.body_gap = body.distance;
        best.fallback = body.fallback != GeometryFallbackLevel::none ? body.fallback : reference_meta.fallback;
        best.used_quadratic_fit = reference_meta.used_quadratic_fit;
        best.iterations = iter + 1;
      }

      if (Norm(body.X - X) < tol && final_dirichlet_gap < tol) {
        converged = true;
        break;
      }
      X = body.X;
    }

    if (found) {
      best.body_gap_threshold = target.body_gap_threshold;
      best.shift_limit = target.shift_limit;
      best.status = !converged ? GeometryTargetStatus::no_convergence
                               : (best.fallback == GeometryFallbackLevel::none
                                      ? GeometryTargetStatus::ok
                                      : GeometryTargetStatus::ok_via_fallback);
      applyClamp(best, X_in);
      if (!isFinite(best.target_X_clamped) || !isFinite(best.delta_clamped) ||
          !std::isfinite(best.move_ratio))
        best.status = GeometryTargetStatus::invalid_projection;
      return best;
    }
  }

  Tddd X = X_in;
  double best_score = std::numeric_limits<double>::infinity();
  GeometryTarget best = target;
  bool converged = false;

  for (int iter = 0; iter < max_iter; ++iter) {
    const auto dirichlet = projectToReferenceForEntity(provider.reference, water, X,
                                                       provider.options.mode,
                                                       provider.options.trust_radius_factor,
                                                       entity);
    if (!dirichlet.ok) {
      target.status = GeometryTargetStatus::no_reference_surface;
      return target;
    }
    const auto body = projectToBodyCandidates(entity, dirichlet.X, target.body_gap_threshold,
                                              provider.body_objects, provider.options);
    if (!body.ok) {
      target.status = GeometryTargetStatus::no_body_surface;
      return target;
    }
    const double score = dirichlet.distance * dirichlet.distance + body.distance * body.distance;
    if (std::isfinite(score) && score < best_score) {
      best_score = score;
      best = target;
      best.target_X = body.X;
      best.reference_face = dirichlet.face;
      best.param_uv = dirichlet.param;
      best.body_face = body.face;
      best.dirichlet_gap = dirichlet.distance;
      best.body_gap = body.distance;
      best.fallback = body.fallback != GeometryFallbackLevel::none ? body.fallback : dirichlet.fallback;
      best.used_quadratic_fit = dirichlet.used_quadratic_fit;
      best.iterations = iter + 1;
    }
    if (Norm(body.X - X) < tol) {
      converged = true;
      break;
    }
    X = body.X;
  }

  if (!isFinite(best.target_X) || !std::isfinite(best.body_gap)) {
    target.status = GeometryTargetStatus::invalid_projection;
    return target;
  }
  if (best.body_gap > target.body_gap_threshold) {
    target.status = GeometryTargetStatus::body_gap_too_large;
    return target;
  }
  best.body_gap_threshold = target.body_gap_threshold;
  best.shift_limit = target.shift_limit;
  best.status = !converged ? GeometryTargetStatus::no_convergence
                           : (best.fallback == GeometryFallbackLevel::none
                                  ? GeometryTargetStatus::ok
                                  : GeometryTargetStatus::ok_via_fallback);
  applyClamp(best, X_in);
  if (!isFinite(best.target_X_clamped) || !isFinite(best.delta_clamped) || !std::isfinite(best.move_ratio))
    best.status = GeometryTargetStatus::invalid_projection;
  return best;
}

inline GeometryTarget queryGeometryTarget(const BoundaryGeometryTargetProvider& provider,
                                          const GeometryProjectorQuery& query) {
  if (query.p)
    return queryGeometryTargetForEntity(provider, query.p, query.water, query.X_in);
  if (query.l)
    return queryGeometryTargetForEntity(provider, query.l, query.water, query.X_in);
  return {};
}

inline bool targetStatusAllowsCurvature(GeometryTargetStatus status) {
  return status == GeometryTargetStatus::ok ||
         status == GeometryTargetStatus::ok_via_fallback ||
         status == GeometryTargetStatus::no_convergence;
}

inline bool providerTargetCurvatureCacheUsable(const ProviderTargetCurvatureCache& cache) {
  return cache.valid && cache.curvature.valid &&
         isFinite(cache.source_X) && isFinite(cache.target_X) &&
         isFinite(cache.normal) && Norm(cache.normal) > 1e-14;
}

inline bool directionalCurvatureFromProviderCache(const ProviderTargetCurvatureCache& cache,
                                                  const Tddd& chord,
                                                  double& k_dir) {
  if (!providerTargetCurvatureCacheUsable(cache) || !isFinite(chord))
    return false;
  const double chord_len = Norm(chord);
  if (!(chord_len > 1e-14) || !std::isfinite(chord_len))
    return false;
  const Tddd n = Normalize(cache.normal);
  Tddd tangent = chord - Dot(chord, n) * n;
  const double tangent_len = Norm(tangent);
  if (!(tangent_len > 1e-14) || !isFinite(tangent))
    return false;
  tangent /= tangent_len;
  const double c1 = Dot(tangent, cache.curvature.PD1);
  const double c2 = Dot(tangent, cache.curvature.PD2);
  k_dir = cache.curvature.k1 * c1 * c1 + cache.curvature.k2 * c2 * c2;
  return std::isfinite(k_dir);
}

inline bool lineSampleIsMidpoint(double s) {
  return !std::isfinite(s) || std::abs(s - 0.5) <= 1e-12;
}

template <class Entity>
inline TargetThetaResult queryTargetThetaForEntity(const BoundaryGeometryTargetProvider& provider,
                                                   Entity* entity,
                                                   const Network* water,
                                                   const Tddd& X_in,
                                                   const Tddd& edge) {
  using namespace GeometryProjectorDetail;
  TargetThetaResult out;
  if (!entity || !isFinite(edge) || Norm(edge) <= 1e-14)
    return out;

  const auto target = queryGeometryTargetForEntity(provider, entity, water, X_in);
  if (!targetStatusAllowsCurvature(target.status))
    return out;

  const auto flags = provider.targetFlags(entity);
  if (flags.need_reference && target.reference_face) {
    const auto curv = referenceCurvatureAt(target.reference_face, target.param_uv, provider.options.mode);
    accumulateThetaFromCurvature(out, edge, curv, /*reference_target=*/true);
  }
  if (flags.need_body && target.body_face) {
    const auto curv = bodyCurvatureAt(target.body_face, target.target_X,
                                      provider.options.feature_angle_rad, /*bfs_depth=*/2);
    accumulateThetaFromCurvature(out, edge, curv, /*reference_target=*/false);
  }
  return out;
}

inline GeometryTarget BoundaryGeometryTargetProvider::queryEntityTarget(networkPoint* p,
                                                                        const Network* water,
                                                                        const Tddd& X_in) const {
  GeometryProjectorQuery query;
  query.p = p;
  query.water = water;
  query.X_in = X_in;
  return queryGeometryTarget(*this, query);
}

inline GeometryTarget BoundaryGeometryTargetProvider::queryEntityTarget(networkLine* l,
                                                                        const Network* water,
                                                                        const Tddd& X_in) const {
  return queryLineSampleTarget(l, water, 0.5, X_in);
}

inline bool BoundaryGeometryTargetProvider::ensureProviderTargetCurvature(networkPoint* p,
                                                                          const Network* water) const {
  using namespace GeometryProjectorDetail;
  if (!p || !isFinite(p->X))
    return false;

  const auto flags = targetFlags(p);
  auto& cache = p->provider_target_curvature;
  const double source_tol = std::max(localLength(p) * 1e-10, 1e-12);
  if (provider_epoch != 0 &&
      providerTargetCurvatureCacheUsable(cache) &&
      cache.provider_epoch == provider_epoch &&
      cache.need_reference == flags.need_reference &&
      cache.need_body == flags.need_body &&
      Norm(cache.source_X - p->X) <= source_tol) {
    return true;
  }

  cache.invalidate();
  cache.provider_epoch = provider_epoch;
  cache.need_reference = flags.need_reference;
  cache.need_body = flags.need_body;
  cache.source_X = p->X;
  if (!flags.need_reference && !flags.need_body)
    return false;

  const auto target = queryGeometryTargetForEntity(*this, p, water, p->X);
  if (!targetStatusAllowsCurvature(target.status))
    return false;

  surface_geometry::PrincipalCurvatureResult curv;
  Tddd normal = {0., 0., 0.};
  bool used_reference = false;
  bool used_body = false;

  if (flags.need_reference && target.reference_face) {
    curv = referenceCurvatureAt(target.reference_face, target.param_uv, options.mode);
    normal = snapshotProjectionNormal(*target.reference_face, target.param_uv, options.mode);
    used_reference = curv.valid;
  }
  if (!curv.valid && flags.need_body && target.body_face) {
    curv = bodyCurvatureAt(target.body_face, target.target_X,
                           options.feature_angle_rad, /*bfs_depth=*/2);
    normal = target.body_face->normal;
    used_body = curv.valid;
  }

  if (!curv.valid || !isFinite(normal) || Norm(normal) <= 1e-14)
    return false;

  cache.valid = true;
  cache.provider_epoch = provider_epoch;
  cache.need_reference = flags.need_reference;
  cache.need_body = flags.need_body;
  cache.used_reference = used_reference;
  cache.used_body = used_body;
  cache.source_X = p->X;
  cache.target_X = target.target_X;
  cache.normal = Normalize(normal);
  cache.curvature = curv;
  return true;
}

inline Tddd BoundaryGeometryTargetProvider::curvatureSeedForLineSample(networkLine* l,
                                                                       const Network* water,
                                                                       double s,
                                                                       const Tddd& X_seed,
                                                                       bool* used_curvature_seed) const {
  if (used_curvature_seed)
    *used_curvature_seed = false;
  if (!l || !lineSampleIsMidpoint(s) || !isFinite(X_seed))
    return X_seed;

  auto [p0, p1] = l->getPoints();
  if (!p0 || !p1 || !isFinite(p0->X) || !isFinite(p1->X))
    return X_seed;
  if (!ensureProviderTargetCurvature(p0, water) ||
      !ensureProviderTargetCurvature(p1, water))
    return X_seed;

  const auto& c0 = p0->provider_target_curvature;
  const auto& c1 = p1->provider_target_curvature;
  if (!providerTargetCurvatureCacheUsable(c0) ||
      !providerTargetCurvatureCacheUsable(c1))
    return X_seed;

  // The first use is intentionally Dirichlet-only.  BCInterface/body targets
  // have additional contact constraints, so the final Provider projection
  // remains the only authority there.
  if (c0.need_body || c1.need_body ||
      !c0.used_reference || !c1.used_reference)
    return X_seed;

  Tddd chord = p1->X - p0->X;
  const double L = Norm(chord);
  if (!(L > 1e-14) || !std::isfinite(L))
    return X_seed;

  Tddd n0 = Normalize(c0.normal);
  Tddd n1 = Normalize(c1.normal);
  if (Dot(n0, n1) <= 0.25)
    return X_seed;

  double k0 = 0.0;
  double k1 = 0.0;
  if (!directionalCurvatureFromProviderCache(c0, chord, k0) ||
      !directionalCurvatureFromProviderCache(c1, chord, k1))
    return X_seed;

  const double theta = std::max(std::abs(k0), std::abs(k1)) * L;
  if (!std::isfinite(theta) || theta > 0.75)
    return X_seed;

  Tddd n_mid = n0 + n1;
  const double n_mid_norm = Norm(n_mid);
  if (!(n_mid_norm > 1e-14) || !isFinite(n_mid))
    return X_seed;
  n_mid /= n_mid_norm;

  double sagitta = 0.125 * 0.5 * (k0 + k1) * L * L;
  const double max_sagitta = 0.25 * L;
  if (!std::isfinite(sagitta))
    return X_seed;
  sagitta = std::clamp(sagitta, -max_sagitta, max_sagitta);
  const Tddd seeded = X_seed + sagitta * n_mid;
  if (!isFinite(seeded))
    return X_seed;
  if (used_curvature_seed)
    *used_curvature_seed = std::abs(sagitta) > 1e-14;
  return seeded;
}

inline ProviderFieldSample BoundaryGeometryTargetProvider::sampleReferenceField(const GeometryTarget& target,
                                                                                ProviderFieldKind kind) const {
  ProviderFieldSample sample;
  sample.kind = kind;
  if (!target.reference_face || !isFinite(target.target_X) ||
      !std::isfinite(target.param_uv[0]) || !std::isfinite(target.param_uv[1]))
    return sample;

  const auto& values = (kind == ProviderFieldKind::Phi)
                           ? target.reference_face->field.phi6
                           : target.reference_face->field.phi_t6;
  const auto N = TriShape<6>(target.param_uv[0], target.param_uv[1]);
  double value = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i)
    value += N[i] * values[i];
  if (!std::isfinite(value))
    return sample;
  sample.ok = true;
  sample.value = value;
  return sample;
}

inline GeometryTarget BoundaryGeometryTargetProvider::queryLineSampleTarget(networkLine* l,
                                                                            const Network* water,
                                                                            double s,
                                                                            const Tddd& X_seed) const {
  Tddd X_in = X_seed;
  if (!isFinite(X_in) && l) {
    auto [p0, p1] = l->getPoints();
    if (p0 && p1) {
      const double a = std::clamp(std::isfinite(s) ? s : 0.5, 0.0, 1.0);
      X_in = (1.0 - a) * p0->X + a * p1->X;
    }
  }
  if (l && isFinite(X_in)) {
    auto [p0, p1] = l->getPoints();
    if (p0 && p1 && isFinite(p0->X) && isFinite(p1->X)) {
      const double a = std::clamp(std::isfinite(s) ? s : 0.5, 0.0, 1.0);
      const Tddd x_linear = (1.0 - a) * p0->X + a * p1->X;
      const double tol = std::max(Norm(p1->X - p0->X) * 1e-10, 1e-12);
      if (isFinite(x_linear) && Norm(X_in - x_linear) <= tol)
        X_in = curvatureSeedForLineSample(l, water, a, X_in);
    }
  }
  GeometryProjectorQuery query;
  query.l = l;
  query.water = water;
  query.X_in = X_in;
  return queryGeometryTarget(*this, query);
}

inline TargetThetaResult BoundaryGeometryTargetProvider::queryTargetTheta(networkPoint* p,
                                                                          const Network* water,
                                                                          const Tddd& X_in,
                                                                          const Tddd& edge) const {
  return queryTargetThetaForEntity(*this, p, water, X_in, edge);
}

inline TargetThetaResult BoundaryGeometryTargetProvider::queryTargetTheta(networkLine* l,
                                                                          const Network* water,
                                                                          const Tddd& X_in,
                                                                          const Tddd& edge) const {
  return queryTargetThetaForEntity(*this, l, water, X_in, edge);
}

} // namespace BEMMeshPipeline
