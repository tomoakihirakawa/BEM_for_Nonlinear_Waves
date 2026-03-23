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

Tddd pushPosition(RungeKutta<Tddd>& RK, const Tddd& u_node, const Tddd& fallback, std::size_t& nan_count) {
  RK.push(u_node);
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
                         RungeKutta<Tddd>& RK, Tddd& u_node, bool do_eta_phi,
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
    auto nextX = RK.getX(u_node); //次時刻の位置を計算
    const double nextT = RK.getTimeAtNextStep();
    if (isFinite(nextX)) {
      const double eta = absorber->absorb_eta(nextX, nextT);
      const double dt_rk = RK.getdt();
      const double to_eta_in_z = eta - nextX[2]; //次時刻の位置と吸収面までのz距離の差
      if (std::isfinite(dt_rk) && dt_rk > 0. && std::isfinite(to_eta_in_z)) {
        const double u2_new = u_node[2] + gamma * to_eta_in_z / dt_rk;
        if (std::isfinite(u2_new))
          u_node[2] = u2_new;
        else
          ++nan_pos_count;
      } else
        ++nan_pos_count;
      auto nextX_abs = RK.getX(u_node);
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
                            double mean_phi) {
  const bool debug_nan_guard = []() {
    if (const char* env = std::getenv("BEM_RELOCATION_NAN_DEBUG"))
      return std::string(env) != "0";
    return false;
  }();
  std::size_t nan_guard_point_count = 0;
  std::size_t nan_guard_phi_count = 0;

  for (auto* node : fluid_nodes) {
    node->u_absorbed.fill(0.);
    node->phi_absorbed = 0.;
    node->absorb_gamma = 0.;
    double ref_phi = 0;
    if (node->absorbedBy != nullptr) {
      bool has_dirichlet = hasAnyDirichletBoundaryState(node);
      auto uz_before = node->u_node[2];
      ref_phi = computeAbsorption(node->absorbedBy,
                                  node->signed_distance,
                                  node->RK_X,
                                  node->u_node,
                                  has_dirichlet,
                                  mean_phi,
                                  nan_guard_point_count, nan_guard_phi_count);
      node->u_absorbed[2] = node->u_node[2] - uz_before;
      node->absorb_gamma = node->absorbedBy->absorb_gamma(node->signed_distance);
      if (!hasAnyDirichletBoundaryState(node) && node->CORNER) {
        throw std::runtime_error("Error: Absorption for CORNER nodes without Dirichlet faces is not supported. Please ensure that all CORNER nodes have at least one Dirichlet face or remove the CORNER classification.");
      }
    }

    if (!isFinite(node->u_node)) {
      node->u_node = isFinite(node->u_total) ? node->u_total : Tddd{0., 0., 0.};
      ++nan_guard_point_count;
    }
    if (!isFinite(node->RK_X.getX(node->u_node))) {
      node->u_node = {0., 0., 0.};
      ++nan_guard_point_count;
    }

    if (hasAnyDirichletBoundaryState(node)) {
      double dphi_dt = node->DphiDt_damped({node->absorb_gamma, ref_phi}, node->u_node, 0.);
      double dphi_dt_undamped = node->DphiDt_damped({0., 0.}, node->u_node, 0.);
      node->phi_absorbed = (dphi_dt - dphi_dt_undamped) * node->RK_phi.getdt();
      std::get<0>(node->phiphin) = pushPhi(node->RK_phi, dphi_dt, nan_guard_phi_count);
    }
    node->setXSingle(pushPosition(node->RK_X, node->u_node, node->getPosition(), nan_guard_point_count));
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

    for (auto* p : water->getPoints())
      p->copy_phiphin();
    for (auto* l : water->getBoundaryLines())
      l->copy_phiphin();

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

    auto func = [&](auto p) {
      if (hasAnyDirichletBoundaryState(p) && p->relocation_face) {
        auto* f = p->relocation_face;
        auto [t0, t1] = p->relocation_param;
        const double phi_new = interpolate_scalar(f, t0, t1, 0);
        // const double phi_new = std::get<0>(p->phiphin);
        std::get<0>(p->phiphin) = phi_new;
        for (auto& [face, d] : p->dofs)
          if (isDirichletBieDofKey(p, face))
            d.phi = phi_new;
      }
    };

    // 2. Points: interpolate phi/phin using cached (face, t0, t1)
    for (auto* p : water->getPoints())
      func(p);

    // 3. Lines: true_quadratic interpolation using cached (face, t0, t1)
    if (use_true_quadratic) {
      for (auto* l : water->getBoundaryLines())
        func(l);
    }

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
  }
}

} // namespace
