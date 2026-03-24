#ifndef BEM_MESH_RELAXATION_HPP
#define BEM_MESH_RELAXATION_HPP

#include "BEM_calculateVelocities.hpp"
#include "BEM_remesh_main.hpp"
#include "BEM_setBoundaryTypes.hpp"

/* ========================================================================== */
/*  Pre-relaxation: flip + smoothing to improve initial mesh quality           */
/*                                                                             */
/*  Called at t=0 to improve OBJ-loaded mesh before time stepping.             */
/*  When initial conditions are present (solitary wave etc.), called again      */
/*  after IC application to recover mesh quality from X_mid wave projection.   */
/* ========================================================================== */

struct PreRelaxParams {
   int loop;                     // smoothing iteration count
   double coef;                  // smoothing coefficient
   int cycles = 1;               // flip + smoothing cycle count
   std::string output_tag = "";  // output file tag (empty = no output)
};

inline void preRelaxMesh(
    const std::vector<Network*>& FluidObject,
    const std::vector<Network*>& AllObjects,
    const std::vector<Network*>& contact_objects,
    bool use_true_quadratic,
    NodeRelocationSurface& node_relocation_surface_ref,
    double simulation_time,
    const PreRelaxParams& params,
    std::function<void(Network*, const std::string&)> write_mesh = nullptr) {

   TimeWatch pre_watch;
   std::cout << Green << "[preRelaxMesh] begin (loop=" << params.loop
             << ", coef=" << params.coef
             << ", cycles=" << params.cycles
             << ")" << colorReset << std::endl;

   // Prepare spatial acceleration for contact queries.
   _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);

   for (auto& water : FluidObject) {
      constexpr double relax_dt = 1.0;
      const double rad = std::acos(-1.0) / 180.0;
      const Tdd flip_limitD = {20.0 * rad, 20.0 * rad};
      const Tdd flip_limitN = {5.0 * rad, 5.0 * rad};

      auto init_relax_rk = [&]() {
         for (const auto& p : water->getPoints()) {
            p->u_reloc.fill(0.0);
            p->RK_X.initialize(relax_dt, simulation_time, ToX(p), 1);
         }
         if (use_true_quadratic) {
            for (auto* l : water->getBoundaryLines()) {
               l->u_reloc.fill(0.0);
               l->RK_X.initialize(relax_dt, simulation_time, l->X_mid, 1);
            }
         }
      };

      // Always initialize X_mid to straight-edge midpoints.
      // On first call: OBJ leaves X_mid at {0,0,0}.
      // After IC: X_mid holds wave-projected values that may create tiny subfaces.
      // In both cases, resetting to linear midpoints is correct because:
      // - pre-relax moves vertices only, making old X_mid stale
      // - IC re-applies X_mid from theory after pre-relax completes
      if (use_true_quadratic) {
         for (auto* l : water->getBoundaryLines()) {
            auto [pA, pB] = l->getPoints();
            l->setXSingle(0.5 * (pA->X + pB->X));
         }
      }

      refreshBoundaryStatesAndTypes(water, contact_objects);
      if (use_true_quadratic)
         computeAllCornerMidpointOffsets(water);

      // Pre-relaxation uses pseudo_quadratic surface: X_mid may not have
      // curvature info (1st call) or may have wave-distorted values (2nd call).
      // pseudo_quadratic uses neighbor reconstruction for better approximation.
      auto saved_surface = node_relocation_surface_ref;
      if (node_relocation_surface_ref == NodeRelocationSurface::true_quadratic)
         node_relocation_surface_ref = NodeRelocationSurface::pseudo_quadratic;

      for (int c = 0; c < params.cycles; ++c) {
         init_relax_rk();
         water->setGeometricPropertiesForce();
         flipIfBatched(*water, flip_limitD, flip_limitN, "pre-relax");
         calculateVecToSurface(*water, params.loop, params.coef);
         for (const auto& p : water->getPoints())
            p->setXSingle(RK_with_Ubuff(p));
      }

      node_relocation_surface_ref = saved_surface;
      water->setGeometricPropertiesForce();

      if (write_mesh && !params.output_tag.empty())
         write_mesh(water, params.output_tag);
   }

   if (bem_verbose())
      PrintLap(pre_watch, "[preRelaxMesh] done");
}

#endif // BEM_MESH_RELAXATION_HPP
