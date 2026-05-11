#pragma once

#include "BEM_full_interpolation.hpp"

namespace BEMMeshPipeline {
namespace FieldMapper {

// Thin wrapper by design: this names the pipeline role that maps a
// ReferenceState's field values onto the already-repaired geometry/topology.
// The current implementation delegates to fullPassInterpolation so default
// behavior stays unchanged while future mapper variants can be introduced here.
inline void mapFieldsFromReferenceState(
    const std::vector<Network*>& fluid_nets,
    const ReferenceState& reference,
    NodeRelocationSurface scheme,
    MeshPreparationPseudoQuadraticPolicy pseudo_policy = MeshPreparationPseudoQuadraticPolicy::error) {
  fullPassInterpolation(fluid_nets, reference, scheme, pseudo_policy);
}

inline void mapFieldsFromReferenceState(
    const std::vector<Network*>& fluid_nets,
    const PreModificationSnapshot& snap,
    NodeRelocationSurface scheme,
    MeshPreparationPseudoQuadraticPolicy pseudo_policy = MeshPreparationPseudoQuadraticPolicy::error) {
  fullPassInterpolation(fluid_nets, snap, scheme, pseudo_policy);
}

} // namespace FieldMapper
} // namespace BEMMeshPipeline
