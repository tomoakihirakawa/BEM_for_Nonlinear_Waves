// BEM_remesh_main.cpp
//
// Dispatcher for remesh_for_main_loop().
//
// Trial/local_patch is the only active remeshing path.

#define BEM
#include "Network.hpp"
#include "BEM_remesh_main.hpp"
#include "BEM_remesh_common.hpp"
#include "BEM_inputfile_reader.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

namespace {

std::string to_lower_copy(const std::string& s) {
  std::string out = s;
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return out;
}

}  // namespace

void remesh_for_main_loop(Network& water, int time_step,
                          const SimulationSettings::RemeshingSettings& rs,
                          bool skip_post_remesh_quality_rejects,
                          const std::string& patch_output_directory,
                          double simulation_time,
                          PVDWriter* candidate_patches_pvd,
                          PVDWriter* remeshed_patches_pvd,
                          PVDWriter* trigger_edges_pvd,
                          PVDWriter* edges_pvd,
                          const BEMMeshPipeline::RemeshTrialGeometryInputs* trial_geometry_inputs,
                          const std::string& phase_debug_tag,
                          PVDWriter* remeshed_water_unrepaired_pvd,
                          int step_retry,
                          RemeshWaterSnapshotWriter remeshed_water_snapshot_writer) {
  const std::string m = to_lower_copy(rs.remesh_method);

  if (m.empty() || m == "trial" || m == "local_patch_trial_multi_objective") {
    return remesh_local_patch_trial_multi_objective(
        water, time_step, rs, skip_post_remesh_quality_rejects,
        patch_output_directory, simulation_time,
        candidate_patches_pvd, remeshed_patches_pvd,
        trigger_edges_pvd, edges_pvd, trial_geometry_inputs, phase_debug_tag,
        remeshed_water_unrepaired_pvd, step_retry, remeshed_water_snapshot_writer);
  }

  throw std::runtime_error(
      "unknown rs.remesh_method: '" + rs.remesh_method +
      "' — valid remesh_method is Trial.");
}
