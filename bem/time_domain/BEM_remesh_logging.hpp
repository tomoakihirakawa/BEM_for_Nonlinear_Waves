#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "BEM_remesh_debug_output.hpp"

namespace BEMMeshPipeline::RemeshLog {

inline double safeRatio(long long numerator, long long denominator) {
  return denominator > 0 ? static_cast<double>(numerator) / static_cast<double>(denominator) : 0.0;
}

inline const char* boundaryName(int bc) {
  switch (bc) {
    case 0: return "N";
    case 1: return "D";
    case 2: return "BCI";
    default: return "?";
  }
}

inline const char* rejectCodeName(int rc) {
  switch (rc) {
    case 0: return "success";
    case 1: return "no_valid";
    case 2: return "hd_exceeded";
    case 3: return "score_no_improve";
    case 4: return "replace_failed";
    case 5: return "repair_failed";
    case 6: return "patch_ops_failed";
    case 7: return "border_moved";
    case 8: return "subsurface_altitude";
    case 9: return "bcinterface_no_contact";
    case 10: return "exception";
    default: return "unknown";
  }
}

struct StageTiming {
  double split = 0.0;
  double collapse = 0.0;
  double smooth_ir = 0.0;
  double smooth_angle = 0.0;
  double curv = 0.0;
  double topo = 0.0;

  double total() const {
    return split + collapse + smooth_ir + smooth_angle + curv + topo;
  }
};

struct RunTrialTiming {
  long long calls = 0;
  double ref = 0.0;
  double copy = 0.0;
  double run = 0.0;
  double pick = 0.0;

  double total() const {
    return ref + copy + run + pick;
  }
};

struct HdLimitStats {
  long long checks = 0;
  long long early_exits = 0;
  long long samples = 0;
  long long nearest_calls = 0;
  long long target_missing = 0;
  double max_ref = 0.0;
  double max_body = 0.0;

  double earlyExitRatio() const {
    return safeRatio(early_exits, checks);
  }
};

struct RepairAwareTrialStats {
  bool enabled = false;
  long long source_calls = 0;
  long long repaired_points = 0;
  long long repaired_lines = 0;
  long long repair_failed = 0;
  long long altitude_rescue_attempts = 0;
  long long altitude_rescue_successes = 0;
  long long altitude_rescue_split_worst_line = 0;
  long long altitude_rescue_split_longest_face_edge = 0;
  long long subsurface_trigger_edges = 0;
  long long best_rejected_altitude_cases = 0;
  double best_rejected_altitude_before = 0.0;
  double best_rejected_altitude_after = 0.0;
  double best_rejected_altitude_gain = 0.0;
  double max_free_gap = 0.0;
  double max_body_gap = 0.0;
  double max_move_ratio = 0.0;
};

struct PhaseParts {
  double needs = 0.0;
  double run_trials = 0.0;
  double replace = 0.0;
  double postop = 0.0;
};

struct GeometryRepairSummary {
  std::string phase;
  int iterations = 0;
  double max_move = 0.0;
  int nodes_adjusted = 0;
  int lines_snapped = 0;
  int repaired_points = -1;
  int repaired_lines = -1;
  int quality_damaged = -1;
  int damaged_faces = -1;
  int penetration_rejected = -1;
};

inline void printPhaseParts(std::ostream& os, const char* label, const PhaseParts& p) {
  os << "[remesh:parts] " << label
     << " needs=" << p.needs << "s"
     << " runTrials=" << p.run_trials << "s"
     << " replace=" << p.replace << "s"
     << " postop=" << p.postop << "s"
     << std::endl;
}

inline void printLocalPatchProfile(int step,
                                   const StageTiming& stage,
                                   const RunTrialTiming& trial,
                                   const HdLimitStats& hd,
                                   const RepairAwareTrialStats& repair,
                                   const PhaseParts& split_parts,
                                   const PhaseParts& collapse_parts,
                                   std::ostream& os = std::cout) {
  os << "[remesh:profile] method=local_patch"
     << " step=" << step
     << " total=" << stage.total() << "s"
     << " split=" << stage.split << "s"
     << " collapse=" << stage.collapse << "s"
     << " smooth_ir=" << stage.smooth_ir << "s"
     << " smooth_angle=" << stage.smooth_angle << "s"
     << " curv=" << stage.curv << "s"
     << " topo=" << stage.topo << "s"
     << std::endl;

  os << "[remesh:trial] calls=" << trial.calls
     << " total=" << trial.total() << "s"
     << " ref=" << trial.ref << "s"
     << " copy=" << trial.copy << "s"
     << " run=" << trial.run << "s"
     << " pick=" << trial.pick << "s"
     << std::endl;

  os << "[remesh:hd] checks=" << hd.checks
     << " early=" << hd.early_exits
     << " ratio=" << hd.earlyExitRatio()
     << " samples=" << hd.samples
     << " nearest=" << hd.nearest_calls
     << " target_missing=" << hd.target_missing
     << " max_ref=" << hd.max_ref
     << " max_body=" << hd.max_body
     << std::endl;

  if (repair.enabled) {
    os << "[remesh:repair-trial]"
       << " source_calls=" << repair.source_calls
       << " repaired_points=" << repair.repaired_points
       << " repaired_lines=" << repair.repaired_lines
       << " repair_failed=" << repair.repair_failed
       << " altitude_rescue_attempts=" << repair.altitude_rescue_attempts
       << " altitude_rescue_successes=" << repair.altitude_rescue_successes
       << " split_worst_line=" << repair.altitude_rescue_split_worst_line
       << " split_longest_face_edge=" << repair.altitude_rescue_split_longest_face_edge
       << " subsurface_trigger_edges=" << repair.subsurface_trigger_edges
       << " best_rejected_altitude_cases=" << repair.best_rejected_altitude_cases
       << " best_before=" << repair.best_rejected_altitude_before
       << " best_after=" << repair.best_rejected_altitude_after
       << " best_gain=" << repair.best_rejected_altitude_gain
       << " max_free_gap=" << repair.max_free_gap
       << " max_body_gap=" << repair.max_body_gap
       << " max_move_ratio=" << repair.max_move_ratio
       << std::endl;
  }

  printPhaseParts(os, "split", split_parts);
  printPhaseParts(os, "collapse", collapse_parts);
}

inline void printRejectSummary(const std::vector<RemeshDebug::TriggerEdgeRecord>& records,
                               int reject_code_count,
                               std::ostream& os = std::cout) {
  constexpr int kBoundaryCount = 3;
  constexpr int kMaxRejectCodes = 16;
  std::array<int, kBoundaryCount> bc_count{0, 0, 0};
  std::array<int, kBoundaryCount> succ_by_bc{0, 0, 0};
  std::array<std::array<int, kMaxRejectCodes>, kBoundaryCount> reject_by_bc{};

  const int n_codes = std::min(reject_code_count, kMaxRejectCodes);
  for (const auto& r : records) {
    if (r.bc_type < 0 || r.bc_type >= kBoundaryCount)
      continue;
    ++bc_count[r.bc_type];
    if (r.success == 1) {
      ++succ_by_bc[r.bc_type];
    } else if (r.reject_code >= 0 && r.reject_code < n_codes) {
      ++reject_by_bc[r.bc_type][r.reject_code];
    }
  }

  os << "[remesh:reject] by_trigger_bc";
  for (int bc = 0; bc < kBoundaryCount; ++bc) {
    if (bc_count[bc] == 0)
      continue;
    os << " " << boundaryName(bc)
       << "{total=" << bc_count[bc]
       << ",success=" << succ_by_bc[bc];
    for (int rc = 1; rc < n_codes; ++rc)
      os << "," << rejectCodeName(rc) << "=" << reject_by_bc[bc][rc];
    os << "}";
  }
  os << std::endl;
}

inline void printPipelineGeometryRepair(const GeometryRepairSummary& s,
                                        std::ostream& os = std::cout) {
  os << "[mesh:geometry]"
     << " phase=" << s.phase
     << " iter=" << s.iterations
     << " max_move=" << s.max_move
     << " nodes=" << s.nodes_adjusted
     << " lines=" << s.lines_snapped;
  if (s.repaired_points >= 0)
    os << " repaired_points=" << s.repaired_points;
  if (s.repaired_lines >= 0)
    os << " repaired_lines=" << s.repaired_lines;
  if (s.quality_damaged >= 0)
    os << " quality_damaged=" << s.quality_damaged;
  if (s.damaged_faces >= 0)
    os << " damaged_faces=" << s.damaged_faces;
  if (s.penetration_rejected >= 0)
    os << " penetration_rejected=" << s.penetration_rejected;
  os << std::endl;
}

inline void printPipelineTopologyRun(const std::string& reason,
                                     std::ostream& os = std::cout) {
  os << "[mesh:topology] run reason=" << reason << std::endl;
}

inline void printPipelineTopologySummary(const std::string& phase,
                                         bool requested,
                                         int run_count,
                                         int remaining_damaged,
                                         std::ostream& os = std::cout) {
  os << "[mesh:topology]"
     << " phase=" << phase
     << " requested=" << (requested ? 1 : 0)
     << " run=" << run_count
     << " remaining_damaged=" << remaining_damaged
     << std::endl;
}

inline void printPipelineWaterlineTopologyRequest(int damaged_lines,
                                                  int damaged_faces,
                                                  std::ostream& os = std::cout) {
  os << "[mesh:waterline]"
     << " topology_repair_requested=1"
     << " damaged_lines=" << damaged_lines
     << " damaged_faces=" << damaged_faces
     << std::endl;
}

inline void printPipelineWaterlineRetrySummary(int remaining_damaged,
                                               int damaged_faces,
                                               std::ostream& os = std::cout) {
  os << "[mesh:waterline]"
     << " retry_remaining_damaged=" << remaining_damaged
     << " damaged_faces=" << damaged_faces
     << std::endl;
}

} // namespace BEMMeshPipeline::RemeshLog
