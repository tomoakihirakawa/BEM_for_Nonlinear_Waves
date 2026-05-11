#pragma once

#include "BEM_node_face_state.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_set>

namespace {

struct DebugMidpointTarget {
  bool enabled = false;
  Tddd x = {0., 0., 0.};
  double tol = 5e-2;
};

inline const DebugMidpointTarget& getDebugMidpointTarget() {
  static const DebugMidpointTarget target = [] {
    DebugMidpointTarget t;
    if (const char* env = std::getenv("BEM_DEBUG_MIDPOINT_LINE")) {
      std::stringstream ss(env);
      char comma0 = 0, comma1 = 0;
      if (ss >> t.x[0] >> comma0 >> t.x[1] >> comma1 >> t.x[2]) {
        if (comma0 == ',' && comma1 == ',')
          t.enabled = true;
      }
    }
    if (const char* env_tol = std::getenv("BEM_DEBUG_MIDPOINT_TOL")) {
      try {
        t.tol = std::stod(env_tol);
      } catch (...) {
      }
    }
    return t;
  }();
  return target;
}

inline std::string midpointDebugFaceKeySummary(networkFace* f) {
  if (!f)
    return "nullptr";
  auto [p0, p1, p2] = f->getPoints();
  std::array<int, 3> ids = {p0 ? p0->index : -1, p1 ? p1->index : -1, p2 ? p2->index : -1};
  std::sort(ids.begin(), ids.end());
  std::ostringstream oss;
  oss << "ptr=" << f
      << " idx=" << f->index
      << " pts={" << ids[0] << "," << ids[1] << "," << ids[2] << "}"
      << " centroid=" << f->centroid;
  return oss.str();
}

inline const networkLine* findDebugMidpointLine(const Network* water, double* nearest_dist = nullptr) {
  const auto& target = getDebugMidpointTarget();
  if (!target.enabled || water == nullptr)
    return nullptr;

  const networkLine* best = nullptr;
  double best_dist = std::numeric_limits<double>::max();
  for (const auto* l : water->getBoundaryLines()) {
    const double d = Norm(l->X_mid - target.x);
    if (d < best_dist) {
      best_dist = d;
      best = l;
    }
  }
  if (nearest_dist)
    *nearest_dist = best_dist;
  if (best_dist > target.tol)
    return nullptr;
  return best;
}

inline void dumpDebugMidpointLineState(const Network* water, const char* stage, int time_step, int rk_step) {
  const auto& target = getDebugMidpointTarget();
  if (!target.enabled || water == nullptr)
    return;

  double nearest_dist = std::numeric_limits<double>::max();
  const auto* l = findDebugMidpointLine(water, &nearest_dist);
  if (!l) {
    std::cout << Magenta << "[midpoint_debug] stage=" << stage
              << " water=" << water->getName()
              << " time_step=" << time_step
              << " rk_step=" << rk_step
              << " target=" << target.x
              << " no line within tol=" << target.tol
              << " nearest_dist=" << nearest_dist
              << colorReset << std::endl;
    return;
  }

  auto boundary_faces = l->getBoundaryFaces();
  std::unordered_set<networkFace*> boundary_face_set(boundary_faces.begin(), boundary_faces.end());
  std::stable_sort(boundary_faces.begin(), boundary_faces.end(), [](const auto* a, const auto* b) {
    return reinterpret_cast<std::uintptr_t>(a) < reinterpret_cast<std::uintptr_t>(b);
  });

  auto [pA, pB] = l->getPoints();
  const auto* d0 = l->findActiveBieDof(nullptr);
  double min_face_phin = std::numeric_limits<double>::max();
  double max_face_phin = -std::numeric_limits<double>::max();
  std::size_t face_phin_count = 0;
  for (const auto& [f, d] : l->dofs) {
    if (f == nullptr || d.index < 0)
      continue;
    min_face_phin = std::min(min_face_phin, d.phin);
    max_face_phin = std::max(max_face_phin, d.phin);
    ++face_phin_count;
  }
  const double face_phin_jump = (face_phin_count > 0) ? (max_face_phin - min_face_phin) : 0.;
  const double rep_minus_null = d0 ? (l->phiphin[1] - d0->phin) : 0.;
  const Tddd rk_base = (l->RK_X.steps == 0) ? l->X_mid : l->RK_X.getX(l->u_reloc);

  std::cout << Magenta << "[midpoint_debug] stage=" << stage
            << " water=" << water->getName()
            << " time_step=" << time_step
            << " rk_step=" << rk_step
            << " target=" << target.x
            << " tol=" << target.tol
            << " line_ptr=" << l
            << " endpts={" << (pA ? pA->index : -1) << "," << (pB ? pB->index : -1) << "}"
            << " X_mid=" << l->X_mid
            << " flags={D:" << l->Dirichlet << ",N:" << l->Neumann << ",C:" << l->BCInterface << ",M:" << l->isMultipleNode << "}"
            << " midpoint_index=" << l->midpoint_index
            << " face_count=" << boundary_faces.size()
            << " phiphin={" << l->phiphin[0] << "," << l->phiphin[1] << "}"
            << " nullptr_phin=" << (d0 ? d0->phin : 0.)
            << " rep-null=" << rep_minus_null
            << " face_phin_jump=" << face_phin_jump
            << " relocation_face={" << midpointDebugFaceKeySummary(l->relocation_face) << "}"
            << " relocation_param={" << l->relocation_param[0] << "," << l->relocation_param[1] << "}"
            << " rk_steps=" << l->RK_X.steps
            << " rk_base=" << rk_base
            << " vecToSurface=" << l->vecToSurface
            << " u_reloc=" << l->u_reloc
            << " dist=" << nearest_dist
            << colorReset << std::endl;

  std::cout << "  boundary_faces:" << std::endl;
  for (auto* f : boundary_faces) {
    const auto bc = getEdgeNodeFaceBoundaryType(l, f);
    const char* bc_name = (bc == NodeFaceBoundaryType::Dirichlet) ? "D" :
                          (bc == NodeFaceBoundaryType::Neumann) ? "N" : "U";
    std::cout << "    {" << midpointDebugFaceKeySummary(f) << "} bc=" << bc_name << std::endl;
  }

  std::vector<std::pair<networkFace*, const NodeFaceState*>> dof_entries;
  dof_entries.reserve(l->dofs.size());
  for (const auto& [f, d] : l->dofs)
    dof_entries.emplace_back(f, &d);
  std::stable_sort(dof_entries.begin(), dof_entries.end(), [](const auto& a, const auto& b) {
    return reinterpret_cast<std::uintptr_t>(a.first) < reinterpret_cast<std::uintptr_t>(b.first);
  });

  std::cout << "  dofs size=" << dof_entries.size() << std::endl;
  for (const auto& [f, d] : dof_entries) {
    const bool alive = (f == nullptr) || boundary_face_set.count(f);
    const auto bc = (f != nullptr) ? getEdgeNodeFaceBoundaryType(l, f) : NodeFaceBoundaryType::Undefined;
    const char* bc_name = (bc == NodeFaceBoundaryType::Dirichlet) ? "D" :
                          (bc == NodeFaceBoundaryType::Neumann) ? "N" :
                          (f == nullptr ? "-" : "U");
    std::cout << "    key={" << midpointDebugFaceKeySummary(f) << "} alive=" << alive
              << " bc=" << bc_name
              << " index=" << d->index
              << " phi=" << d->phi
              << " phin=" << d->phin
              << " phi_t=" << d->phi_t
              << " phin_t=" << d->phin_t
              << " contacts=" << d->contact_opponent_faces.size()
              << std::endl;
  }
}

} // namespace
