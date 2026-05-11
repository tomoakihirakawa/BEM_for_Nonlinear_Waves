#pragma once

#include <stdexcept>
#include <string>
#include <utility>

#include "BEM_field_provider.hpp"

namespace BEMMeshPipeline {

template <class Entity>
inline bool assignInterpolatedDirichletValue(Entity* entity,
                                             const Network* water,
                                             const ReferenceState& reference,
                                             NodeRelocationSurface scheme,
                                             std::size_t* updated_per_face_dofs = nullptr) {
  if (!hasAnyDirichletBoundaryState(entity))
    return false;
  FieldQuery query;
  query.water = water;
  query.X = entity->getPosition();
  query.bc_type = NodeFaceBoundaryType::Dirichlet;
  query.surface = scheme;

  auto result = evaluate(reference, query);
  if (result.status != FieldStatus::ok || !result.phi)
    return false;
  std::get<0>(entity->phiphin) = *result.phi;
  if (result.phi_t)
    std::get<0>(entity->phiphin_t) = *result.phi_t;
  for (auto& [face, dof] : entity->dofs) {
    if (!isDirichletBieDofKey(entity, face))
      continue;
    dof.phi = *result.phi;
    if (result.phi_t)
      dof.phi_t = *result.phi_t;
    if (updated_per_face_dofs)
      ++(*updated_per_face_dofs);
  }
  return true;
}

inline void validateFullPassInterpolationPolicy(const ReferenceState& reference,
                                                NodeRelocationSurface scheme,
                                                MeshPreparationPseudoQuadraticPolicy pseudo_policy) {
  if (scheme != NodeRelocationSurface::pseudo_quadratic || !isSnapshotReference(reference))
    return;
  if (pseudo_policy == MeshPreparationPseudoQuadraticPolicy::error) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                        "mesh preparation pipeline does not use pseudo_quadratic interpolation by default. "
                        "Set node_relocation surface to linear/true_quadratic, or explicitly set "
                        "mesh_preparation_pseudo_quadratic_policy=fallback_6node to use the 6-node fallback.");
  }
  static bool warned = false;
  if (!warned) {
    std::cout << Yellow
              << "[mesh:field] pseudo_quadratic fallback=6node"
              << colorReset << std::endl;
    warned = true;
  }
}

inline void fullPassInterpolation(const std::vector<Network*>& fluid_nets,
                                  const ReferenceState& reference,
                                  NodeRelocationSurface scheme,
                                  MeshPreparationPseudoQuadraticPolicy pseudo_policy = MeshPreparationPseudoQuadraticPolicy::error) {
  if (const auto* snapshot = std::get_if<SnapshotReferenceState>(&reference)) {
    if (snapshot->snap.faces.empty())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                          "mesh preparation pipeline snapshot has no Dirichlet reference faces");
  }
  validateFullPassInterpolationPolicy(reference, scheme, pseudo_policy);

  std::size_t assigned_points = 0;
  std::size_t assigned_lines = 0;
  std::size_t updated_per_face_dofs = 0;
  std::size_t missed = 0;

  for (auto* water : fluid_nets) {
    for (auto* p : water->getBoundaryPoints()) {
      if (!hasAnyDirichletBoundaryState(p))
        continue;
      if (assignInterpolatedDirichletValue(p, water, reference, scheme, &updated_per_face_dofs))
        ++assigned_points;
      else
        ++missed;
    }
    for (auto* l : water->getBoundaryLines()) {
      if (!hasAnyDirichletBoundaryState(l))
        continue;
      if (assignInterpolatedDirichletValue(l, water, reference, scheme, &updated_per_face_dofs))
        ++assigned_lines;
      else
        ++missed;
    }
  }

  if (missed > 0)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                        "fullPassInterpolation could not assign all Dirichlet entities");

  std::cout << Green << "[mesh:field] full_pass_interpolation"
            << " assigned_points=" << assigned_points
            << " assigned_midpoints=" << assigned_lines
            << " reference="
            << (isSnapshotReference(reference) ? std::to_string(std::get<SnapshotReferenceState>(reference).snap.faces.size()) + " snapshot faces" : "initial-condition reference")
            << " updated_per_face_dirichlet_dofs=" << updated_per_face_dofs
            << colorReset << std::endl;
}

inline void fullPassInterpolation(const std::vector<Network*>& fluid_nets,
                                  const PreModificationSnapshot& snap,
                                  NodeRelocationSurface scheme,
                                  MeshPreparationPseudoQuadraticPolicy pseudo_policy = MeshPreparationPseudoQuadraticPolicy::error) {
  SnapshotReferenceState reference;
  reference.snap = snap;
  buildSnapshotFaceIndex(reference.snap, *reference.index);
  fullPassInterpolation(fluid_nets, ReferenceState{std::move(reference)}, scheme, pseudo_policy);
}

} // namespace BEMMeshPipeline
