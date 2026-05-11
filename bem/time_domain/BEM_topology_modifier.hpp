#pragma once

#include "BEM_remesh_main.hpp"
#include "BEM_remesh_trial_geometry_inputs.hpp"

namespace BEMMeshPipeline {
namespace TopologyModifier {

// Thin wrapper by design: split/collapse/flip changes topology, while
// GeometryRepair is responsible for projecting the resulting mesh back to the
// reference/body geometry.  Keeping this call boundary explicit makes the main
// pipeline order auditable without changing remesh behavior.
inline void applyLocalRemesh(Network& water,
                             int time_step,
                             const SimulationSettings::RemeshingSettings& settings,
                             bool skip_post_remesh_quality_rejects,
                             const std::string& patch_output_directory,
                             double simulation_time,
                             PVDWriter* candidate_patches_pvd,
                             PVDWriter* remeshed_patches_pvd,
                             PVDWriter* trigger_edges_pvd,
                             PVDWriter* edges_pvd,
                             const RemeshTrialGeometryInputs* trial_geometry_inputs = nullptr,
                             const std::string& phase_debug_tag = "",
                             PVDWriter* remeshed_water_unrepaired_pvd = nullptr,
                             int step_retry = 0,
                             RemeshWaterSnapshotWriter remeshed_water_snapshot_writer = RemeshWaterSnapshotWriter{}) {
  remesh_for_main_loop(water, time_step, settings, skip_post_remesh_quality_rejects,
                       patch_output_directory, simulation_time,
                       candidate_patches_pvd, remeshed_patches_pvd,
                       trigger_edges_pvd, edges_pvd, trial_geometry_inputs, phase_debug_tag,
                       remeshed_water_unrepaired_pvd, step_retry,
                       remeshed_water_snapshot_writer);
}

} // namespace TopologyModifier
} // namespace BEMMeshPipeline
