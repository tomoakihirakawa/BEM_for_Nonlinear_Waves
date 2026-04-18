#pragma once

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include "BEM_step_failure.hpp"
#include "Network.hpp"

void sync_body_states_from_rk(const std::vector<Network*>& RigidBodyObject,
                              const std::vector<Network*>& SoftBodyObject);

enum class RetryAction {
  ContinueRetry,
  BreakStep
};

struct StepRetryState {
  static constexpr int max_step_retries = 5;
  static constexpr double dt_reduction_factor = 0.5;

  double dt_override = 0.0;
  bool degraded_mode = false;
  int current_rk_step = 0;
  std::unordered_map<std::string, std::unordered_map<int, int>> repeated_reject_faces;

  void remember_rejected_face(const step_failure& e);
  void collapse_repeatedly_rejected_faces(Network& water, int step_retry);
  RetryAction handle_step_failure(const step_failure& e,
                                  int step_retry,
                                  double dt,
                                  double last_successful_dt,
                                  int time_step,
                                  std::function<void(int, int)> write_failure_snapshot_fn);
};
