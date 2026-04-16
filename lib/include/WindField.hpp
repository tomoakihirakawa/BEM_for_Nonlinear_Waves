#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <random>

#include "basic_arithmetic_array_operations.hpp"

/*
 * WindField.hpp — ambient fluid velocity generators for LumpedCable.
 *
 * Returns std::function<Tddd(const Tddd& X, double t)> that maps (position,
 * time) to a fluid velocity vector [m/s]. Used by LumpedCable::wind_field and
 * read by the 3 getDragForce call sites in LumpedCable.hpp.
 *
 * MVP scope (2026-04-15):
 *   - makeUniform:   time-/space-constant U
 *   - makeAR1:       time-varying Ornstein-Uhlenbeck turbulence, no spatial
 *                    variation (all nodes see the same U(t) — i.e. perfectly
 *                    correlated turbulence, which is a reasonable first-cut
 *                    for buffeting response of a single bridge cable).
 *
 * Out of scope (future): Davenport/von Karman spectrum synthesis, spatial
 *                        coherence Coh(f, |Δx|), height profiles.
 */

namespace WindField {

inline std::function<Tddd(const Tddd&, double)>
makeUniform(const Tddd& U_mean) {
   return [U_mean](const Tddd& /*X*/, double /*t*/) -> Tddd { return U_mean; };
}

/*
 * AR(1) / Ornstein-Uhlenbeck turbulence.
 *
 *   dx_i/dt = -x_i / T_L + sigma_u * sqrt(2 / T_L) * eta_i(t),   i = x, y, z
 *   U(t)    = U_mean + x(t)
 *
 * Discrete update with substep dt:
 *   x_{n+1} = x_n * exp(-dt / T_L) + sigma_u * sqrt(1 - exp(-2 dt / T_L)) * xi_n
 *
 * where xi_n ~ N(0, 1). sigma_u = TI * |U_mean| (turbulence intensity).
 *
 * The returned closure captures a shared state so that all callers (all
 * cable nodes, all 3 call sites) observing the same (or monotonically
 * increasing) time see a coherent turbulence trajectory. If t goes backwards
 * or stays the same, the last-computed value is returned without advancing
 * the state (idempotent during a single RK stage).
 */
inline std::function<Tddd(const Tddd&, double)>
makeAR1(const Tddd& U_mean, double TI, double T_L, unsigned seed) {
   struct State {
      Tddd U_mean;
      double sigma_u;   // = TI * |U_mean|
      double T_L;       // [s]
      std::mt19937 rng;
      std::normal_distribution<double> N01{0., 1.};
      Tddd x{0., 0., 0.};     // OU fluctuation state
      double last_t = -std::numeric_limits<double>::infinity();
      Tddd last_U;            // cached return value for repeated same-t queries
   };

   auto state = std::make_shared<State>();
   state->U_mean = U_mean;
   state->sigma_u = TI * Norm(U_mean);
   state->T_L = (T_L > 0.) ? T_L : 1.0;
   state->rng.seed(seed);
   state->last_U = U_mean;

   return [state](const Tddd& /*X*/, double t) -> Tddd {
      if (t <= state->last_t) return state->last_U;
      double dt = t - state->last_t;
      // First call (last_t = -inf): skip decay term, just evaluate U_mean.
      if (!std::isfinite(state->last_t)) {
         state->last_t = t;
         state->last_U = state->U_mean + state->x;
         return state->last_U;
      }
      double decay = std::exp(-dt / state->T_L);
      double diffusion = state->sigma_u * std::sqrt(std::max(0., 1. - decay * decay));
      Tddd xi{state->N01(state->rng), state->N01(state->rng), state->N01(state->rng)};
      state->x = state->x * decay + xi * diffusion;
      state->last_t = t;
      state->last_U = state->U_mean + state->x;
      return state->last_U;
   };
}

}  // namespace WindField
