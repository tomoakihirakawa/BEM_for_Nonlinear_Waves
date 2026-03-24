#pragma once

// Absorption-zone damping, RK position/phi push, interpolation relocation,
// and mean-phi normalization — extracted from main_time_domain.cpp.
//
// All functions live in an anonymous namespace so they are file-local.

#include <cmath>
#include <iostream>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BEM_midpoint_debug.hpp"
#include "Network.hpp"
#include "BEM_node_face_state.hpp"
#include "BEM_calculateVelocities.hpp"

namespace {

// ============================================================================
//  NaN-guarded RK push helpers
// ============================================================================

Tddd pushPosition(RungeKutta<Tddd>& RK, const Tddd& u_reloc, const Tddd& fallback, std::size_t& nan_count) {
  RK.push(u_reloc);
  auto X = RK.getX();
  if (!isFinite(X)) {
    RK.repush(Tddd{0., 0., 0.});
    X = RK.getX();
    if (!isFinite(X))
      X = fallback;
    ++nan_count;
  }
  return X;
}

double pushPhi(RungeKutta<double>& RK, double dphi_dt, std::size_t& nan_count) {
  if (!std::isfinite(dphi_dt)) {
    dphi_dt = 0.;
    ++nan_count;
  }
  RK.push(dphi_dt);
  double phi = RK.getX();
  if (!std::isfinite(phi)) {
    RK.repush(0.);
    phi = RK.getX();
    if (!std::isfinite(phi))
      phi = 0.;
    ++nan_count;
  }
  return phi;
}

// ============================================================================
//  Absorption-zone gamma / eta / ref_phi computation
// ============================================================================

double computeAbsorption(const Network* absorber, double signed_distance,
                         RungeKutta<Tddd>& RK, Tddd& u_reloc, bool do_eta_phi,
                         double mean_phi,
                         std::size_t& nan_pos_count, std::size_t& nan_phi_count) {
  double gamma = absorber->absorb_gamma(signed_distance);
  if (!std::isfinite(gamma)) {
    gamma = 0.;
    ++nan_phi_count;
  }
  double ref_phi = 0.;
  if (do_eta_phi) {
    ref_phi = mean_phi;
    auto nextX = RK.getX(u_reloc); //次時刻の位置を計算
    const double nextT = RK.getTimeAtNextStep();
    if (isFinite(nextX)) {
      const double eta = absorber->absorb_eta(nextX, nextT);
      const double dt_rk = RK.getdt();
      const double to_eta_in_z = eta - nextX[2]; //次時刻の位置と吸収面までのz距離の差
      if (std::isfinite(dt_rk) && dt_rk > 0. && std::isfinite(to_eta_in_z)) {
        const double u2_new = u_reloc[2] + gamma * to_eta_in_z / dt_rk;
        if (std::isfinite(u2_new))
          u_reloc[2] = u2_new;
        else
          ++nan_pos_count;
      } else
        ++nan_pos_count;
      auto nextX_abs = RK.getX(u_reloc);
      const double phi_abs = isFinite(nextX_abs) ? absorber->absorb_phi(nextX_abs, nextT) : std::numeric_limits<double>::quiet_NaN();
      if (std::isfinite(phi_abs))
        ref_phi = phi_abs + mean_phi;
      else
        ++nan_phi_count;
    } else {
      ++nan_pos_count;
    }
  }
  return ref_phi;
}

// ============================================================================
//  Signed distance computation
// ============================================================================

void computeSignedDistances(const std::vector<BEM_DOF_Base*>& fluid_nodes) {
  _Pragma("omp parallel for") for (const auto& node : fluid_nodes) {
    if (node->absorbedBy != nullptr) {
      auto [f, X_nearest] = node->absorbedBy->Nearest(node->getPosition());
      node->signed_distance = Norm(node->getPosition() - X_nearest);
    } else
      node->signed_distance = 0;
  }
}

// ============================================================================
//  Mean phi over fluid faces
// ============================================================================

double computeMeanPhi(const std::unordered_set<networkFace*>& fluid_faces) {
  double mean_phi = 0., total_area = 0;
  for (const auto& f : fluid_faces) {
    auto [p0, p1, p2] = f->getPoints();
    mean_phi += (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3 * f->area;
    total_area += f->area;
  }
  if (!(total_area > 0.) || !std::isfinite(total_area))
    return 0.;
  return mean_phi / total_area;
}

// ============================================================================
//  Apply absorption + NaN-guarded RK push to all fluid points and midpoints
// ============================================================================

void applyAbsorptionAndPush(const std::vector<BEM_DOF_Base*>& fluid_nodes,
                            double mean_phi,
                            bool use_ale = false) {
  const bool debug_nan_guard = []() {
    if (const char* env = std::getenv("BEM_RELOCATION_NAN_DEBUG"))
      return std::string(env) != "0";
    return false;
  }();
  std::size_t nan_guard_point_count = 0;
  std::size_t nan_guard_phi_count = 0;

  for (auto* node : fluid_nodes) {
    // ALE: u_reloc (= u_total + smoothing) for both position and phi push.
    //      Position and phi evolve consistently in the relocated frame.
    // non-ALE (none/interpolation): u_total (pure Lagrangian) for push.
    //      For interpolation, position is moved to X_reloc and phi is
    //      re-interpolated separately after RK4 completes.
    Tddd u_push = use_ale ? node->u_reloc : node->u_total;

    node->u_absorbed.fill(0.);
    node->phi_absorbed = 0.;
    node->absorb_gamma = 0.;
    double ref_phi = 0;
    if (node->absorbedBy != nullptr) {
      bool has_dirichlet = hasAnyDirichletBoundaryState(node);
      auto uz_before = u_push[2];
      ref_phi = computeAbsorption(node->absorbedBy,
                                  node->signed_distance,
                                  node->RK_X,
                                  u_push,
                                  has_dirichlet,
                                  mean_phi,
                                  nan_guard_point_count, nan_guard_phi_count);
      node->u_absorbed[2] = u_push[2] - uz_before;
      node->absorb_gamma = node->absorbedBy->absorb_gamma(node->signed_distance);
      if (!hasAnyDirichletBoundaryState(node) && node->CORNER) {
        throw std::runtime_error("Error: Absorption for CORNER nodes without Dirichlet faces is not supported. Please ensure that all CORNER nodes have at least one Dirichlet face or remove the CORNER classification.");
      }
    }

    if (!isFinite(u_push)) {
      u_push = isFinite(node->u_total) ? node->u_total : Tddd{0., 0., 0.};
      ++nan_guard_point_count;
    }
    if (!isFinite(node->RK_X.getX(u_push))) {
      u_push = {0., 0., 0.};
      ++nan_guard_point_count;
    }

    if (hasAnyDirichletBoundaryState(node)) {
      double dphi_dt = node->DphiDt_damped({node->absorb_gamma, ref_phi}, u_push, 0.);
      double dphi_dt_undamped = node->DphiDt_damped({0., 0.}, u_push, 0.);
      node->phi_absorbed = (dphi_dt - dphi_dt_undamped) * node->RK_phi.getdt();
      std::get<0>(node->phiphin) = pushPhi(node->RK_phi, dphi_dt, nan_guard_phi_count);
    }
    node->setXSingle(pushPosition(node->RK_X, u_push, node->getPosition(), nan_guard_point_count));
    node->phi_tmp = 0;
  }
  if (debug_nan_guard && (nan_guard_point_count > 0 || nan_guard_phi_count > 0)) {
    std::cout << Magenta << "[relocation:nan-guard] corrected_u=" << nan_guard_point_count
              << " corrected_phi=" << nan_guard_phi_count << colorReset << std::endl;
  }
}

// ============================================================================
//  POST-RK interpolation relocation
// ============================================================================

// These enums are defined in main_time_domain.cpp before this header is included.
// Forward-declared here so the function signature is self-documenting.
// enum class NodeRelocationSurface { linear, pseudo_quadratic, true_quadratic };
// enum class InterpolationMidpointMode { nearest };

void applyInterpolationRelocation(const std::vector<Network*>& FluidObject,
                                  bool use_true_quadratic,
                                  NodeRelocationSurface surface,
                                  InterpolationMidpointMode midpoint_mode) {
  for (auto water : FluidObject) {

    // 1. Copy phiphin at the Lagrangian-push position (before relocation moves nodes).
    //    These values were computed with u_total (not u_reloc), so they are
    //    consistent with the current node positions.
    for (auto* p : water->getPoints())
      p->copy_phiphin();
    for (auto* l : water->getBoundaryLines())
      l->copy_phiphin();

    // 2. Set relocation_face / relocation_param for X_reloc on the current mesh.
    //    X_reloc is the smoothing-corrected target position. We find where it sits
    //    on the current (Lagrangian-pushed) Dirichlet surface, so that phi can be
    //    interpolated from phiphin_copy at that parametric location.
    auto setRelocationParam = [](auto* entity) {
      const bool has_dirichlet_state = std::ranges::any_of(entity->getBoundaryFaces(), [&](const auto* f) {
        return getNodeFaceBoundaryType(entity, f) == NodeFaceBoundaryType::Dirichlet;
      });
      if (!has_dirichlet_state) {
        entity->relocation_face = nullptr;
        entity->relocation_param = {0., 0.};
        return;
      }
      const Tddd target = entity->X_reloc;
      double best_dist = 1E+20;
      networkFace* best_face = nullptr;
      Tdd best_param = {0., 0.};
      for (auto* f : entity->getBoundaryFaces()) {
        if (f->penetratedBody || getNodeFaceBoundaryType(entity, f) != NodeFaceBoundaryType::Dirichlet)
          continue;
        Tdd param;
        // Use actual post-push coordinates (not RK_without_Ubuff) since RK4 is already complete.
        auto actualPos = [](const auto* q) -> Tddd { return getNodeX(q); };
        Tddd X_near = NearestOnDirichletFace(target, f, &param, actualPos);
        double dist = Norm(X_near - target);
        if (dist < best_dist) {
          best_dist = dist;
          best_face = f;
          best_param = param;
        }
      }
      if (best_face) {
        entity->relocation_face = best_face;
        entity->relocation_param = best_param;
      } else {
        entity->relocation_face = nullptr;
        entity->relocation_param = {0., 0.};
      }
    };
    for (const auto& p : water->getPoints())
      setRelocationParam(p);
    for (auto* l : water->getBoundaryLines())
      setRelocationParam(l);

    // 3. Interpolate phi from phiphin_copy using the face/param set above.
    auto interpolate_scalar = [&](networkFace* f, const double t0, const double t1, bool phi0_phin1 /*0 or 1*/) -> double {
      auto [p0, l0, p1, l1, p2, l2] = f->PLPLPL;
      if (surface == NodeRelocationSurface::true_quadratic) {
        const auto N = TriShape<6>(t0, t1);
        return N[0] * p0->phiphin_copy[phi0_phin1] + N[1] * p1->phiphin_copy[phi0_phin1] + N[2] * p2->phiphin_copy[phi0_phin1] + N[3] * l0->phiphin_copy[phi0_phin1] + N[4] * l1->phiphin_copy[phi0_phin1] + N[5] * l2->phiphin_copy[phi0_phin1];
      } else if (surface == NodeRelocationSurface::pseudo_quadratic && f->dodecaPoints[0]) {
        return f->dodecaPoints[0]->interpolate(t0, t1, [&](networkPoint* q) -> double { return q->phiphin_copy[phi0_phin1]; });
      } else
        return t0 * p0->phiphin_copy[phi0_phin1] + t1 * p1->phiphin_copy[phi0_phin1] + (1.0 - t0 - t1) * p2->phiphin_copy[phi0_phin1];
    };

    auto interpolateAndMove = [&](auto* entity) {
      // Interpolate phi at the parametric location of X_reloc on the current
      // (already RK-updated) mesh, then move the node to X_reloc.
      if (hasAnyDirichletBoundaryState(entity) && entity->relocation_face) {
        auto* f = entity->relocation_face;
        auto [t0, t1] = entity->relocation_param;
        const double phi_new = interpolate_scalar(f, t0, t1, 0);
        std::get<0>(entity->phiphin) = phi_new;
        for (auto& [face, d] : entity->dofs)
          if (isDirichletBieDofKey(entity, face))
            d.phi = phi_new;
      }
      // Move node to the relocation target position (computed by setRelocFromVecToSurface)
      entity->setXSingle(entity->X_reloc);
    };

    // 4. Apply to points
    for (auto* p : water->getPoints())
      interpolateAndMove(p);

    // 5. Apply to line midpoints (true_quadratic only)
    if (use_true_quadratic) {
      for (auto* l : water->getBoundaryLines())
        if (l->hasActiveBieDof())
          interpolateAndMove(l);
    }

    // 6. Recompute geometric properties after position change
    water->setGeometricPropertiesForce();

    dumpDebugMidpointLineState(water, "post-interpolation-relocation", -1, -1);
  }
}

// ============================================================================
//  Subtract mean phi (normalize velocity potential)
// ============================================================================

void subtractMeanPhi(const std::unordered_set<networkFace*>& fluid_faces,
                     const std::vector<BEM_DOF_Base*>& fluid_nodes) {
  double mean_phi = 0., total_area = 0;
  for (const auto& f : fluid_faces) {
    auto [p0, p1, p2] = f->getPoints();
    mean_phi += (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3 * f->area;
    total_area += f->area;
  }
  mean_phi /= total_area;
  for (auto* node : fluid_nodes) {
    std::get<0>(node->phiphin) -= mean_phi;
    // Sync dof.phi for Dirichlet DOFs to prevent stale values in next BVP
    for (auto& [face, d] : node->dofs)
      if (isDirichletBieDofKey(node, face))
        d.phi -= mean_phi;
  }
}

} // namespace
