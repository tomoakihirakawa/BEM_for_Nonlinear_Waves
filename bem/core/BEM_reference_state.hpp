#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <variant>

#include "BEM_snapshot.hpp"
#include "BEM_time_domain_types.hpp"
#include "basic_surface_geometry.hpp"
#include "lib_spatial_partitioning.hpp"

namespace BEMMeshPipeline {

struct SnapshotProjection {
  const PreModificationSnapshot::FaceData* face = nullptr;
  Tdd param = {0., 0.};
  Tddd nearest = {0., 0., 0.};
  double distance = std::numeric_limits<double>::infinity();
};

struct SnapshotFaceIndex {
  Buckets<const PreModificationSnapshot::FaceData*> buckets;
  BoundingVolumeHierarchy<const PreModificationSnapshot::FaceData*> linear_bvh;
  double query_range = 0.;
  bool ready = false;
  bool linear_bvh_ready = false;
};

inline T6Tddd snapshotFaceNodes(const PreModificationSnapshot::FaceData& face) {
  auto [p0, p1, p2] = face.geometry.vertices_X;
  return {p0, p1, p2, face.geometry.mid_X[0], face.geometry.mid_X[1], face.geometry.mid_X[2]};
}

inline Tddd snapshotFaceCentroid(const PreModificationSnapshot::FaceData& face) {
  auto [p0, p1, p2] = face.geometry.vertices_X;
  return (p0 + p1 + p2) / 3.0;
}

inline double snapshotFaceMaxEdge(const PreModificationSnapshot::FaceData& face) {
  auto [p0, p1, p2] = face.geometry.vertices_X;
  return std::max({Norm(p0 - p1), Norm(p1 - p2), Norm(p2 - p0)});
}

inline void buildSnapshotFaceIndex(const PreModificationSnapshot& snap, SnapshotFaceIndex& index) {
  index.buckets.clear();
  index.linear_bvh.clear();
  index.query_range = 0.;
  index.ready = false;
  index.linear_bvh_ready = false;
  if (snap.faces.empty())
    return;

  CoordinateBounds bounds;
  double max_edge = 0.;
  for (const auto& face : snap.faces) {
    bounds += CoordinateBounds(face.geometry.vertices_X);
    for (const auto& mid : face.geometry.mid_X)
      bounds += CoordinateBounds(mid);
    max_edge = std::max(max_edge, snapshotFaceMaxEdge(face));
  }

  const double scale = bounds.getScale();
  if (!(max_edge > 0.) || !std::isfinite(max_edge))
    max_edge = (scale > 0. && std::isfinite(scale)) ? scale : 1.0;
  const double dL = std::max(max_edge, scale / 128.0);
  index.query_range = 3.0 * max_edge;
  index.buckets.initialize(bounds.scaledBounds(1.05), dL);
  for (const auto& face : snap.faces) {
    const Tddd centroid = snapshotFaceCentroid(face);
    if (!index.buckets.add(centroid, &face))
      index.buckets.add_bypass_insideQ(centroid, &face);
  }
  index.ready = true;

  std::vector<const PreModificationSnapshot::FaceData*> faces;
  faces.reserve(snap.faces.size());
  for (const auto& face : snap.faces)
    faces.push_back(&face);
  index.linear_bvh.setVector(faces, [](const PreModificationSnapshot::FaceData* face) {
    return face ? CoordinateBounds(face->geometry.vertices_X) : CoordinateBounds();
  });
  index.linear_bvh_ready = !index.linear_bvh.empty();
}

inline Tddd evalSnapshotQuad(const PreModificationSnapshot::FaceData& face, double t0, double t1) {
  return Dot(TriShape<6>(t0, t1), snapshotFaceNodes(face));
}

inline Tddd evalSnapshotQuadDeriv(const PreModificationSnapshot::FaceData& face, double t0, double t1, int di, int dj) {
  const auto X6 = snapshotFaceNodes(face);
  if (di == 1 && dj == 0)
    return Dot(D_TriShape<6, 1, 0>(t0, t1), X6);
  if (di == 0 && dj == 1)
    return Dot(D_TriShape<6, 0, 1>(t0, t1), X6);
  if (di == 2 && dj == 0)
    return Dot(D_TriShape<6, 2, 0>(t0, t1), X6);
  if (di == 0 && dj == 2)
    return Dot(D_TriShape<6, 0, 2>(t0, t1), X6);
  return Dot(D_TriShape<6, 1, 1>(t0, t1), X6);
}

inline SnapshotProjection projectToSnapshotFace(const Tddd& X,
                                                const PreModificationSnapshot::FaceData& face,
                                                NodeRelocationSurface scheme) {
  auto [t0_linear, t1_linear, near_linear, normal] = Nearest_(X, face.geometry.vertices_X);
  (void)normal;

  if (scheme == NodeRelocationSurface::linear) {
    return {&face, {t0_linear, t1_linear}, near_linear, Norm(near_linear - X)};
  }

  Tdd param = refineNearestParam(
      X, Tdd{t0_linear, t1_linear},
      [&](double t0, double t1) { return evalSnapshotQuad(face, t0, t1); },
      [&](double t0, double t1, int di, int dj) { return evalSnapshotQuadDeriv(face, t0, t1, di, dj); });
  const Tddd near_quad = evalSnapshotQuad(face, param[0], param[1]);
  return {&face, param, near_quad, Norm(near_quad - X)};
}

inline SnapshotProjection findNearestSnapshotFace(const Tddd& X,
                                                  const PreModificationSnapshot& snap,
                                                  NodeRelocationSurface scheme,
                                                  const SnapshotFaceIndex* index = nullptr) {
  SnapshotProjection best;
  if (scheme == NodeRelocationSurface::linear && index && index->linear_bvh_ready) {
    const auto hit = index->linear_bvh.nearest(X, [&](const PreModificationSnapshot::FaceData* face) {
      return face ? projectToSnapshotFace(X, *face, NodeRelocationSurface::linear).nearest
                  : Tddd{1E+20, 1E+20, 1E+20};
    });
    if (hit.found && hit.item)
      return projectToSnapshotFace(X, *hit.item, NodeRelocationSurface::linear);
  }

  if (index && index->ready) {
    const auto candidates = index->buckets.getData(X, index->query_range);
    for (const auto* face : candidates) {
      if (!face)
        continue;
      auto candidate = projectToSnapshotFace(X, *face, scheme);
      if (candidate.distance < best.distance)
        best = candidate;
    }
    if (best.face && best.distance <= index->query_range)
      return best;
  }

  for (const auto& face : snap.faces) {
    auto candidate = projectToSnapshotFace(X, face, scheme);
    if (candidate.distance < best.distance)
      best = candidate;
  }
  return best;
}

// SnapshotReferenceState is the current concrete ReferenceState implementation for
// non-initial steps. It owns the reference geometry (vertices + X_mid) and the
// reference field values used by GeometryProjector/FieldMapper.
struct SnapshotReferenceState {
  PreModificationSnapshot snap;
  std::unique_ptr<SnapshotFaceIndex> index = std::make_unique<SnapshotFaceIndex>();

  SnapshotReferenceState() = default;
  SnapshotReferenceState(SnapshotReferenceState&&) noexcept = default;
  SnapshotReferenceState& operator=(SnapshotReferenceState&&) noexcept = default;
  SnapshotReferenceState(const SnapshotReferenceState&) = delete;
  SnapshotReferenceState& operator=(const SnapshotReferenceState&) = delete;
};

struct InitialConditionReferenceState {
  double simulation_time = 0.0;
};

using ReferenceState = std::variant<SnapshotReferenceState, InitialConditionReferenceState>;

inline bool hasAnyInitialConditionReference(const std::vector<Network*>& fluid_nets) {
  return std::ranges::any_of(fluid_nets, [](const Network* water) {
    return water && water->ic_eta && water->ic_phi;
  });
}

inline ReferenceState makeReferenceState(const std::vector<Network*>& fluid_nets,
                                         int time_step,
                                         int start_time_step,
                                         double simulation_time) {
  if (time_step == 0 && start_time_step == 0 && hasAnyInitialConditionReference(fluid_nets)) {
    std::cout << Green << "[mesh:reference] phase0 type=initial_condition time="
              << simulation_time << colorReset << std::endl;
    return InitialConditionReferenceState{simulation_time};
  }

  SnapshotReferenceState reference;
  reference.snap = takeSnapshot(fluid_nets);
  buildSnapshotFaceIndex(reference.snap, *reference.index);
  std::cout << Green << "[mesh:reference] phase0 type=snapshot faces="
            << reference.snap.faces.size() << colorReset << std::endl;
  return ReferenceState{std::move(reference)};
}

inline bool isSnapshotReference(const ReferenceState& reference) {
  return std::holds_alternative<SnapshotReferenceState>(reference);
}

inline bool isInitialConditionReference(const ReferenceState& reference) {
  return std::holds_alternative<InitialConditionReferenceState>(reference);
}

} // namespace BEMMeshPipeline
