#pragma once

#include <vector>

#include "BEM_time_domain_types.hpp"
#include "Network.hpp"

void applyInitialConditions(const std::vector<Network*>& FluidObject,
                            const std::vector<Network*>& RigidBodyObject,
                            const std::vector<Network*>& SoftBodyObject,
                            const std::vector<Network*>& AllObjects,
                            bool use_true_quadratic,
                            NodeRelocationSurface& node_relocation_surface_ref);
