#pragma once

#include <array>
#include <limits>
#include <vector>

#include "BEM_setBoundaryTypes.hpp"
#include "Network.hpp"

namespace BEMMeshPipeline {

// ReferenceState stores the last known good geometry + field state used by
// geometry projection and field mapping.  X_mid is part of the reference
// geometry state, not a disposable visualization helper: for Lagrangian
// quadratic/true-quadratic geometry it has the same semantic weight as a
// vertex geometry DOF.
struct ReferenceFaceGeometry {
  T3Tddd vertices_X;
  std::array<Tddd, 3> mid_X;
};

struct ReferenceFaceField {
  std::array<double, 6> phi6;
  std::array<double, 6> phi_t6;
};

struct ReferenceFaceState {
  ReferenceFaceGeometry geometry;
  ReferenceFaceField field;
  NodeFaceBoundaryType bc_type = NodeFaceBoundaryType::Dirichlet;
};

struct PreModificationSnapshot {
  using FaceData = ReferenceFaceState;

  std::vector<FaceData> faces;
};

inline bool faceHasDirichletSource(const networkFace* f) {
  if (!f)
    return false;
  auto [p0, p1, p2] = f->getPoints();
  if (getNodeFaceBoundaryType(p0, f) == NodeFaceBoundaryType::Dirichlet)
    return true;
  if (getNodeFaceBoundaryType(p1, f) == NodeFaceBoundaryType::Dirichlet)
    return true;
  if (getNodeFaceBoundaryType(p2, f) == NodeFaceBoundaryType::Dirichlet)
    return true;
  auto [l0, l1, l2] = f->getLines();
  if (getNodeFaceBoundaryType(l0, f) == NodeFaceBoundaryType::Dirichlet)
    return true;
  if (getNodeFaceBoundaryType(l1, f) == NodeFaceBoundaryType::Dirichlet)
    return true;
  return getNodeFaceBoundaryType(l2, f) == NodeFaceBoundaryType::Dirichlet;
}

inline Tddd snapshotMidPosition(const networkLine* l) {
  if (l->midpoint_index >= 0)
    return l->X_mid;
  auto [pA, pB] = l->getPoints();
  return 0.5 * (pA->X + pB->X);
}

inline double snapshotMidPhi(const networkLine* l) {
  if (l->midpoint_index >= 0)
    return std::get<0>(l->phiphin);
  auto [pA, pB] = l->getPoints();
  return 0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
}

inline double snapshotMidPhiT(const networkLine* l) {
  if (l->midpoint_index >= 0)
    return std::get<0>(l->phiphin_t);
  auto [pA, pB] = l->getPoints();
  return 0.5 * (std::get<0>(pA->phiphin_t) + std::get<0>(pB->phiphin_t));
}

inline PreModificationSnapshot takeSnapshot(const std::vector<Network*>& fluid_nets) {
  PreModificationSnapshot snap;
  std::size_t reserve_count = 0;
  for (const auto* net : fluid_nets)
    reserve_count += net->getBoundaryFaces().size();
  snap.faces.reserve(reserve_count);

  for (const auto* net : fluid_nets) {
    for (auto* f : net->getBoundaryFaces()) {
      if (!faceHasDirichletSource(f))
        continue;
      auto [p0, p1, p2] = f->getPoints();
      auto [l0, l1, l2] = f->getLines();
      snap.faces.push_back({
          {{p0->X, p1->X, p2->X},
           {snapshotMidPosition(l0), snapshotMidPosition(l1), snapshotMidPosition(l2)}},
          {{std::get<0>(p0->phiphin), std::get<0>(p1->phiphin), std::get<0>(p2->phiphin),
            snapshotMidPhi(l0), snapshotMidPhi(l1), snapshotMidPhi(l2)},
           {std::get<0>(p0->phiphin_t), std::get<0>(p1->phiphin_t), std::get<0>(p2->phiphin_t),
            snapshotMidPhiT(l0), snapshotMidPhiT(l1), snapshotMidPhiT(l2)}},
          NodeFaceBoundaryType::Dirichlet});
    }
  }
  return snap;
}

} // namespace BEMMeshPipeline
