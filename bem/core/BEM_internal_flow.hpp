#pragma once

#include "BEM_BoundaryValues.hpp"
#include "Network.hpp"

void prepareBEMVelocityAtCache(const std::vector<Network *> &fluidObjects);
Tddd getBEMVelocityAt_cached(const Tddd &a, const std::vector<Network *> &fluidObjects);
