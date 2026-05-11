// BEM_remesh_scenario_engine.hpp
//
// Per-edge multi-scenario comparison engine used by the Trial/local_patch
// remesh path. The caller delegates a target edge to scenario_engine::runTrials,
// which:
//
//   1. copies the 2-ring patch around the target edge
//   2. runs each scenario string (e.g. "P", "PFv", "PFv2", "PFvFd", "CFS") in parallel
//      on its own patch copy via runPatchOps
//   3. scores each post-op patch via patchQuality
//   4. returns the best-scoring valid patch (aggressive=true → any valid)
//
// The caller then commits the winning patch via water.replacePatch().
//
// Op codes (scenario strings):
//   'P'  = Split   (linear midpoint + linear phi fallback)
//   'C'  = Collapse (feature_aware_collapse_target + collapse_preserves_normals)
//   'S'  = Smooth  (settings remesh_smooth_mode + feature-aware + step cap)
//   'Sc' = Smooth with CircumradiusToInradius kernel
//   'Sl' = Smooth with AreaLaplacian kernel
//   'Ss', 'Scs', 'Sls'
//        = same smooth variants, but skip points touching sharp edges
//   'F'  = Patch-wide flip pass when valence improves AND the edge violates
//          Delaunay.
//   'Fv' = Patch-wide Flip-Valence  (valence improvement only)
//   'Fv2' = Patch-wide Flip-Valence v2 (SharpQ-sector valence only)
//   'Fd' = Patch-wide Flip-Delaunay (Delaunay violation only)
//   'Fs', 'Fvs', 'Fv2s', 'Fds'
//        = same flip variants, but skip sharp-edge flips.
//          All flip variants still obey boundary, normal, dihedral, and
//          midpoint-drift safety guards. They operate on all eligible non-rim
//          edges in the patch, not only the original trigger edge.
//
// Fold guard: uses a relative threshold computed from the source patch's worst
// face-adjacency dot product, so pre-existing sharp features are not confused
// with op-created folds.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "BEM_inputfile_reader.hpp"

class Network;
class networkLine;

namespace scenario_engine {

// Context gathers the non-RS parameters that all scenario helpers need.
// The caller (each base method's dispatcher) computes these once per step.
struct Context {
  const SimulationSettings::RemeshingSettings& rs;
  double free_surface_target_len;  // for fs-target log² term in patchQuality
  double theta_target;             // 2π / rs.theta_target_N
  double feature_angle;            // rs.feature_angle_deg * π / 180
};

struct TrialResult {
  std::string ops;
  double score = -1.0;
  bool valid = false;
  // 0=success, 1=no valid trial, 3=score not improved, 4=replace_failed
  int reject_code = 1;
  std::shared_ptr<Network> patch;         // post-op winner (pass to replacePatch)
  std::shared_ptr<Network> source_patch;  // pre-op reference (for debug output)
};

// Per-edge multi-scenario trial.
//   water       : mesh to copyLocalPatch from
//   target_line : edge whose 2-ring becomes the patch
//   scenarios   : list of op strings; empty → returns invalid with reject=1
//   aggressive  : true → accept any valid patch
//                 false → require score > score_before + rs.quality_score_improve_margin
//   ctx         : scenario-engine context
// Returns:
//   TrialResult with .patch ready for water.replacePatch() on success.
TrialResult runTrials(Network& water,
                       networkLine* target_line,
                       const std::vector<std::string>& scenarios,
                       bool aggressive,
                       const Context& ctx);

} // namespace scenario_engine
