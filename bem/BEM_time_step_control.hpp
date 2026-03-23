#pragma once

// Time-step retry control, rejected face tracking, body state synchronization.
// Extracted from main_time_domain.cpp.

#include <cctype>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "Network.hpp"
#include "BEM_step_failure.hpp"
#include "BEM_remesh_main.hpp"

namespace {

// ============================================================================
//  Utility parsers
// ============================================================================

std::optional<int> parse_int_after(const std::string& text, const std::string& key) {
  const auto pos = text.find(key);
  if (pos == std::string::npos)
    return std::nullopt;
  std::size_t begin = pos + key.size();
  std::size_t end = begin;
  while (end < text.size() && std::isdigit(static_cast<unsigned char>(text[end])))
    ++end;
  if (end == begin)
    return std::nullopt;
  return std::stoi(text.substr(begin, end - begin));
}

std::optional<std::string> parse_water_name_from_reason(const std::string& reason) {
  const auto on_pos = reason.find(" on ");
  if (on_pos == std::string::npos)
    return std::nullopt;
  const auto at_pos = reason.find(" at time_step", on_pos + 4);
  if (at_pos == std::string::npos || at_pos <= on_pos + 4)
    return std::nullopt;
  return reason.substr(on_pos + 4, at_pos - (on_pos + 4));
}

// ============================================================================
//  Body state synchronization (RK → COM/Q/velocity)
// ============================================================================

void sync_body_states_from_rk(const std::vector<Network*>& RigidBodyObject,
                               const std::vector<Network*>& SoftBodyObject) {
  for (auto* net : RigidBodyObject) {
    if (net->RK_COM.steps > 0)
      net->COM = net->RK_COM.getX();
    if (net->RK_Q.steps > 0)
      net->Q = Normalize(net->RK_Q.getX());
    if (net->RK_Velocity.steps > 0)
      net->velocity = net->RK_Velocity.getX();
    net->setGeometricPropertiesForce();
  }

  for (auto* net : SoftBodyObject) {
    if (net->RK_COM.steps > 0)
      net->COM = net->RK_COM.getX();
    if (net->RK_Q.steps > 0)
      net->Q = Normalize(net->RK_Q.getX());
    if (net->RK_Velocity.steps > 0)
      net->velocity = net->RK_Velocity.getX();
    net->setGeometricPropertiesForce();
  }
}

// ============================================================================
//  Step retry state & control
// ============================================================================

enum class RetryAction { ContinueRetry,
                         BreakStep };

struct StepRetryState {
  static constexpr int max_step_retries = 5;
  static constexpr double dt_reduction_factor = 0.5;

  double dt_override = 0.0;
  bool degraded_mode = false;
  int current_rk_step = 0;
  std::unordered_map<std::string, std::unordered_map<int, int>> repeated_reject_faces;

  void remember_rejected_face(const step_failure& e) {
    const bool is_face_collapse_failure =
        e.reason.find("tiny face area ratio") != std::string::npos ||
        e.reason.find("subsurface face altitude ratio") != std::string::npos;
    if (!is_face_collapse_failure)
      return;
    const auto face_index = parse_int_after(e.reason, "face_index=");
    const auto water_name = parse_water_name_from_reason(e.reason);
    if (!face_index || !water_name)
      return;
    auto& count = repeated_reject_faces[*water_name][*face_index];
    ++count;
    std::cout << Yellow << "[step_reject] remember pathological face: water=" << *water_name
              << " face=" << *face_index << " count=" << count
              << colorReset << std::endl;
  }

  void collapse_repeatedly_rejected_faces(Network& water, const int step_retry) {
    if (step_retry <= 0)
      return;
    auto it = repeated_reject_faces.find(water.getName());
    if (it == repeated_reject_faces.end())
      return;
    bool changed = false;
    for (const auto& [face_index, count] : it->second) {
      if (count < 2)
        continue;
      if (collapseFaceByIndexIfPossible(water, face_index)) {
        changed = true;
        std::cout << Yellow << "[step_reject] collapsed repeatedly rejected face on retry: water="
                  << water.getName() << " face=" << face_index
                  << " count=" << count << " retry=" << step_retry
                  << colorReset << std::endl;
      }
    }
    if (changed) {
      water.setGeometricPropertiesForce();
      water.checkConnectivity();
    }
  }

  // handle_step_failure returns the action for the retry loop.
  // write_failure_snapshot_fn is a callback to dump VTU/checkpoint on max-retry.
  RetryAction handle_step_failure(const step_failure& e, const int step_retry,
                                  double dt, double last_successful_dt, int time_step,
                                  std::function<void(int, int)> write_failure_snapshot_fn) {
    std::cerr << Red << "[step_reject] time_step " << time_step << ": " << e.reason << colorReset << std::endl;
    remember_rejected_face(e);
    if (degraded_mode) {
      std::cerr << Yellow << "[step_reject] DEGRADED MODE also failed. Skipping time_step " << time_step
                << " and advancing." << colorReset << std::endl;
      return RetryAction::BreakStep;
    }
    if (step_retry >= max_step_retries) {
      std::cerr << Yellow << "[step_reject] writing failure snapshot before degraded continuation" << colorReset << std::endl;
      write_failure_snapshot_fn(current_rk_step, step_retry);
      std::cerr << Yellow << "[step_reject] max retries exceeded. Entering DEGRADED CONTINUATION MODE — "
                << "quality checks will be skipped for this step." << colorReset << std::endl;
      degraded_mode = true;
      dt_override = std::max(1E-13, dt_override > 0 ? dt_override * dt_reduction_factor : last_successful_dt * 0.01);
      return RetryAction::ContinueRetry;
    }
    if (dt_override <= 0)
      dt_override = std::max(1E-13, std::min(dt, last_successful_dt) * dt_reduction_factor);
    else
      dt_override = std::max(1E-13, dt_override * dt_reduction_factor);
    return RetryAction::ContinueRetry;
  }
};

} // namespace
