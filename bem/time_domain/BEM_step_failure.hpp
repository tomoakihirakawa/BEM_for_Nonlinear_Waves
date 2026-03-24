#pragma once
#include <stdexcept>
#include <string>

// Exception thrown when a time step should be rejected and retried with smaller dt.
// This is NOT an unrecoverable error — it signals that the step produced unreliable results.
struct step_failure : std::runtime_error {
   explicit step_failure(const std::string& reason)
       : std::runtime_error("[step_failure] " + reason), reason(reason) {}
   std::string reason;
};
