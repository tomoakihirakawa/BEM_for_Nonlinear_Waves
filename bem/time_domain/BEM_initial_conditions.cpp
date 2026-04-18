#include "BEM_initial_conditions.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ranges>

#include "BEM_legacy_globals.hpp"
#include "BEM_node_face_state.hpp"
#include "BEM_setBoundaryTypes.hpp"
#include "BEM_calculateVelocities.hpp"

// Initial condition application at t=0 — extracted from main_time_domain.cpp.
// Applies analytical η(x) and φ(x) to the free surface.
// Called once at the start of simulation when ic_eta/ic_phi are provided.

void applyInitialConditions(const std::vector<Network*>& FluidObject,
                            const std::vector<Network*>& RigidBodyObject,
                            const std::vector<Network*>& SoftBodyObject,
                            const std::vector<Network*>& AllObjects,
                            bool use_true_quadratic,
                            NodeRelocationSurface& node_relocation_surface_ref) {
  bool has_any_ic = std::ranges::any_of(FluidObject, [](const Network* w) { return w->ic_eta && w->ic_phi; });
  if (!has_any_ic)
    return;

  TimeWatch ic_watch;

  // Ensure boundary types are set (needed for Dirichlet/Neumann classification).
  // Must always run here because remesh_for_main_loop may have invalidated
  // boundary types even if initial_mesh_pre_relax already called setBoundaryTypes.
  _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);
  for (auto& water : FluidObject) {
    refreshBoundaryStatesAndTypes(water, Join(RigidBodyObject, SoftBodyObject));
    if (use_true_quadratic)
      computeAllCornerMidpointOffsets(water);
  }

  for (auto& water : FluidObject) {
    if (!water->ic_eta || !water->ic_phi)
      continue;

    std::cout << Green << "[initial_condition] applying to " << water->getName() << colorReset << std::endl;

    // Step 1: Move Dirichlet (free surface) nodes to z = η(x, y, t=0)
    for (const auto& p : water->getPoints()) {
      if (p->Dirichlet || p->CORNER) {
        auto X = ToX(p);
        double eta = water->ic_eta(X, 0.0);
        p->setXSingle({X[0], X[1], eta});
      }
    }
    water->setGeometricPropertiesForce();

    // Step 2: Mesh smoothing to improve mesh quality after deformation
    {
      // Use pseudo_quadratic surface for IC smoothing (same reason as pre-relaxation:
      // X_mid has no curvature info at this point).
      auto saved_surface = node_relocation_surface_ref;
      if (node_relocation_surface_ref == NodeRelocationSurface::true_quadratic)
        node_relocation_surface_ref = NodeRelocationSurface::pseudo_quadratic;
      constexpr double relax_dt = 1.0;
      for (const auto& p : water->getPoints()) {
        p->u_reloc.fill(0.0);
        p->RK_X.initialize(relax_dt, 0.0, ToX(p), 1);
      }
      calculateVecToSurface(*water, 30, 0.05);
      for (const auto& p : water->getPoints()) {
        p->setXSingle(RK_with_Ubuff(p));
        p->vecToSurface.fill(0.);
      }
      water->setGeometricPropertiesForce();
      node_relocation_surface_ref = saved_surface;
    }

    // Step 3: Set φ from theory at the final (post-smoothing) node positions
    {
      int n_dirichlet = 0, n_corner = 0, n_neumann = 0, n_other = 0;
      double phi_min = 1e30, phi_max = -1e30;
      for (const auto& p : water->getPoints()) {
        const bool has_dirichlet = hasAnyDirichletBoundaryState(p);
        const bool has_neumann = hasAnyNeumannBoundaryState(p);
        if (has_dirichlet && !has_neumann)
          ++n_dirichlet;
        else if (has_dirichlet && has_neumann)
          ++n_corner;
        else if (has_neumann)
          ++n_neumann;
        else
          ++n_other;

        if (has_dirichlet) {
          auto X = ToX(p);
          double phi_val = water->ic_phi(X, 0.0);
          std::get<0>(p->phiphin) = phi_val;
          phi_min = std::min(phi_min, phi_val);
          phi_max = std::max(phi_max, phi_val);
        }
      }
      std::cout << Green << "[initial_condition] node classification: "
                << "Dirichlet=" << n_dirichlet << " CORNER=" << n_corner
                << " Neumann=" << n_neumann << " other=" << n_other
                << "\n  phi applied to " << (n_dirichlet + n_corner) << " / "
                << (n_dirichlet + n_corner + n_neumann + n_other) << " nodes"
                << "\n  phi range: [" << phi_min << ", " << phi_max << "]"
                << colorReset << std::endl;
    }

    // Step 4: Set phi_mid and X_mid from theory for Dirichlet-containing lines.
    // 要素タイプに依らず適用する — X_mid と phi_mid は可視化・接触判定でも
    // 参照されるため、linear/pseudo_quadratic でも初期化しておく方が整合する。
    // 前提: remesh_for_main_loop の safety により X_mid は既に 0.5*(pA+pB) で
    // 初期化済み。linear 要素で l->phiphin[0] は RK 中に端点平均で上書きされる。
    {
      int n_mid_set = 0;
      double phi_mid_min = 1e30, phi_mid_max = -1e30;
      double max_dev_from_avg = 0;
      for (auto* l : water->getBoundaryLines()) {
        auto [pA, pB] = l->getPoints();
        if (hasAnyDirichletBoundaryState(l)) {
          auto X_mid_projected = l->X_mid;
          X_mid_projected[2] = water->ic_eta(X_mid_projected, 0.0);
          l->setXSingle(X_mid_projected);

          double phi_val = water->ic_phi(l->X_mid, 0.0);
          double phi_avg = 0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
          l->phiphin[0] = phi_val;
          phi_mid_min = std::min(phi_mid_min, phi_val);
          phi_mid_max = std::max(phi_mid_max, phi_val);
          max_dev_from_avg = std::max(max_dev_from_avg, std::abs(phi_val - phi_avg));
          ++n_mid_set;
        }
      }
      water->setGeometricPropertiesForce();
      std::cout << Green << "[initial_condition] phi_mid set from theory: n=" << n_mid_set
                << " range=[" << phi_mid_min << ", " << phi_mid_max << "]"
                << " max|phi_mid - avg(phiA,phiB)|=" << max_dev_from_avg
                << colorReset << std::endl;
    }

    std::cout << Green << "[initial_condition] done for " << water->getName()
              << Blue << "\nElapsed time: " << Red << ic_watch() << colorReset << " s\n";
  }
}
