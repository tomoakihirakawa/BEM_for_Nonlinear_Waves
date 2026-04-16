#pragma once

// -----------------------------------------------------------------------------
// Backward-compat shim.
//
// `MooringLine` was renamed to `LumpedCable` in Phase 0 of the cable-dynamics
// refactor (2026-04-12). The class logic now lives in LumpedCable.hpp. This
// header re-exposes `MooringLine` as a type alias so that existing BEM /
// cable / example code that still `#include "MooringLine.hpp"` and refers to
// `MooringLine` continues to compile unchanged during the staged migration.
//
// Removal plan: once all call sites have been migrated to `LumpedCable`
// (Phase 2 for BEM, Phase 3 for cable_solver), this shim will be deleted.
// -----------------------------------------------------------------------------

#include "LumpedCable.hpp"

using MooringLine = LumpedCable;
