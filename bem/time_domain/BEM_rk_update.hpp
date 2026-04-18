#pragma once

#include <unordered_set>
#include <vector>

#include "BEM_time_domain_types.hpp"
#include "Network.hpp"

void computeSignedDistances(const std::vector<BEM_DOF_Base*>& fluid_nodes);
double computeMeanPhi(const std::unordered_set<networkFace*>& fluid_faces);
void applyAbsorptionAndPush(const std::vector<BEM_DOF_Base*>& fluid_nodes,
                            double mean_phi,
                            bool use_ale = false);
void applyInterpolationRelocation(const std::vector<Network*>& FluidObject,
                                  bool use_true_quadratic,
                                  NodeRelocationSurface surface,
                                  InterpolationMidpointMode midpoint_mode);
void subtractMeanPhi(const std::unordered_set<networkFace*>& fluid_faces,
                     const std::vector<BEM_DOF_Base*>& fluid_nodes);
