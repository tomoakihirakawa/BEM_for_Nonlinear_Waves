#ifndef BEM_MESH_RELAXATION_HPP
#define BEM_MESH_RELAXATION_HPP

#include <functional>
#include <string>
#include <vector>

#include "BEM_time_domain_types.hpp"
#include "Network.hpp"

struct PreRelaxParams {
   int loop;
   double coef;
   int cycles = 1;
   std::string output_tag = "";
};

void preRelaxMesh(
    const std::vector<Network*>& FluidObject,
    const std::vector<Network*>& AllObjects,
    const std::vector<Network*>& contact_objects,
    bool use_true_quadratic,
    NodeRelocationSurface& node_relocation_surface_ref,
    double simulation_time,
    const PreRelaxParams& params,
    std::function<void(Network*, const std::string&)> write_mesh = nullptr);

#endif
