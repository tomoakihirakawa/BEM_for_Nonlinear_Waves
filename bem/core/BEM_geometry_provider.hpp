#pragma once

#include <optional>

#include "BEM_mesh_adjustment.hpp"

namespace BEMMeshPipeline {

enum class GeometryStatus {
  ok,
  not_applicable,
  no_projection,
  invalid_query
};

struct GeometryQuery {
  networkPoint* p = nullptr;
  networkLine* l = nullptr;
  Tddd X = {0., 0., 0.};
  double move_limit_factor = 0.3;
};

struct GeometryResult {
  GeometryStatus status = GeometryStatus::invalid_query;
  // Raw vector returned by the current surface projection.
  Tddd vector_to_target = {0., 0., 0.};
  // Unclamped target, equal to X + vector_to_target.
  std::optional<Tddd> target_X;
  // Movement limited by currentShiftLimit(entity, move_limit_factor).
  std::optional<Tddd> delta_clamped;
  std::optional<Tddd> target_X_clamped;
  double gap = 0.0;
  std::optional<double> shift_limit;
};

inline GeometryResult evaluateCurrentSurfaceGeometry(const GeometryQuery& q) {
  auto fill = [&](auto* entity) {
    GeometryResult result;
    const Tddd vec = vectorToCurrentSurface(entity, q.X);
    if (!isFinite(vec))
      return GeometryResult{GeometryStatus::no_projection};
    result.status = GeometryStatus::ok;
    result.vector_to_target = vec;
    result.target_X = q.X + vec;
    result.gap = Norm(vec);

    const double limit = currentShiftLimit(entity, q.move_limit_factor);
    result.shift_limit = limit;
    if (result.gap > 1e-12 && limit > 0.0 && std::isfinite(limit)) {
      const Tddd delta = std::min(result.gap, limit) * Normalize(vec);
      result.delta_clamped = delta;
      result.target_X_clamped = q.X + delta;
    }
    return result;
  };

  if (q.p)
    return fill(q.p);
  if (q.l)
    return fill(q.l);
  return {GeometryStatus::invalid_query};
}

} // namespace BEMMeshPipeline
