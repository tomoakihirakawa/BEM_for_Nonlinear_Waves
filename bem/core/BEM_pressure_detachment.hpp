#pragma once

#include <algorithm>

#include "Network.hpp"

inline bool bodyIsFullyFixedForPressureDetachment(const Network* net) {
  if (!net)
    return false;
  return std::ranges::all_of(net->isFixed, [](bool is_fixed) { return is_fixed; });
}

inline bool pressureDetachmentEligible(const networkFace* f) {
  if (!f)
    return false;

  bool found_supporting_body = false;
  auto body_is_movable = [&](const Network* net) {
    if (!net)
      return false;
    found_supporting_body = true;
    return !bodyIsFullyFixedForPressureDetachment(net);
  };

  if (body_is_movable(f->penetratedBody))
    return true;

  auto [p0, p1, p2] = f->getPoints();
  for (const auto* p : {p0, p1, p2}) {
    if (!p)
      continue;
    if (body_is_movable(p->penetratedBody))
      return true;
    if (const auto* d = p->findContactState(f); body_is_movable(d && d->nearestContactFace() ? d->nearestContactFace()->getNetwork() : nullptr))
      return true;
  }

  // Rationale:
  // When a target wave is imposed by absorber/node-relocation, the free surface can
  // move numerically while the contacted wall itself remains fixed. In that case the
  // BCInterface should stay constrained to the wall geometry, and a transient negative
  // pressure must not flip the boundary to Dirichlet. Pressure-based peeling is
  // therefore allowed only when the supporting body is actually movable.
  return !found_supporting_body;
}

template <class Entity>
inline bool pressureDetachmentEligible(const Entity* entity, const networkFace* f) {
  if (!entity || !f)
    return false;

  bool found_supporting_body = false;
  auto body_is_movable = [&](const Network* net) {
    if (!net)
      return false;
    found_supporting_body = true;
    return !bodyIsFullyFixedForPressureDetachment(net);
  };

  if (body_is_movable(entity->penetratedBody))
    return true;

  if (const auto* d = entity->findContactState(f); body_is_movable(d && d->nearestContactFace() ? d->nearestContactFace()->getNetwork() : nullptr))
    return true;

  return !found_supporting_body;
}
