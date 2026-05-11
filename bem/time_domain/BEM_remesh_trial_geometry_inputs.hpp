#pragma once

#include <cstdint>
#include <vector>

#include "../core/BEM_legacy_globals.hpp"
#include "../core/BEM_setBoundaryTypes.hpp"
#include "../core/BEM_geometry_projector.hpp"
#include "../core/BEM_inputfile_reader.hpp"

namespace BEMMeshPipeline {

// Inputs used by Trial remesh to build a BoundaryGeometryTargetProvider.
// The provider is then used on candidate patch copies to ask where a point,
// line midpoint, or sample should lie on the reference/body geometry.
struct RemeshTrialGeometryInputs {
  // Dirichlet/free-surface 側の参照形状。時間発展後の水面を
  // remesh trial 中に比較・投影するための source of truth。
  const ReferenceState* reference = nullptr;

  // Neumann/body 側の現在形状。BCInterface や body-side entity の
  // target を現在の剛体/接触面から決めるために使う。
  const std::vector<Network*>* body_objects = nullptr;

  // geometry_projector_v2 と移動制限など、target provider を作る設定。
  const SimulationSettings::TimeDomainSettings::MeshPreparationPipelineSettings* mesh_pipeline = nullptr;

  // free-surface/reference 側を linear 面で見るか quadratic 面で見るか。
  GeometryProjectionMode projection_mode = GeometryProjectionMode::ReferenceQuadraticFace;

  // provider 内の field/curvature cache を時間ステップごとに分ける世代番号。
  std::uint64_t provider_epoch = 0;

  bool enabled() const {
    return reference && body_objects && mesh_pipeline && mesh_pipeline->geometry_projector_v2;
  }
};

} // namespace BEMMeshPipeline
