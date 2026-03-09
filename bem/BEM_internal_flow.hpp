#pragma once

#include "BEM_BoundaryValues.hpp"
#include "Network.hpp"

inline void prepareBEMVelocityAtCache(const std::vector<Network *> &fluidObjects) {
  for (const auto &water : fluidObjects) {
    water->setGeometricPropertiesForce();
    _Pragma("omp parallel for") for (const auto &integ_f : water->getBoundaryFaces()) integ_f->setIntegrationInfo();
  }
}

inline Tddd getBEMVelocityAt_cached(const Tddd &a, const std::vector<Network *> &fluidObjects) {
  constexpr double singular_eps = 1e-5;
  constexpr double near_ratio = 25.;

  struct NearFaceHit {
    bool hit = false;
    double dist2 = std::numeric_limits<double>::infinity();
    networkFace *face = nullptr;
    Tddd nearest = {0., 0., 0.};
    Tddd unitNormal = {0., 0., 0.};
    Tddd u_interp = {0., 0., 0.};
  };
  NearFaceHit best_hit;

  auto interpolateVelocityAtNearestPointOnFace = [&](networkFace *f) -> NearFaceHit {
    auto [p0, p1, p2] = f->getPoints();
    const auto [w0, w1, nearest, _n] = Nearest_(a, T3Tddd{p0->X, p1->X, p2->X}); // weights for p0, p1 (p2 is implicit)
    const double w2 = 1.0 - w0 - w1;
    Tddd unitNormal = _n;
    if (Norm(unitNormal) > 0.0 && Dot(unitNormal, f->normal) < 0.0)
      unitNormal = unitNormal * (-1.0);
    const auto diff = nearest - a;
    const double dist2 = Dot(diff, diff);
    return {
        true, dist2, f, nearest, unitNormal, p0->u_potential_BEM * w0 + p1->u_potential_BEM * w1 + p2->u_potential_BEM * w2,
    };
  };

  auto get_phi = [&](const networkPoint *p, networkFace *f) -> double {
    if (p->phiOnFace.count(f))
      return p->phiOnFace.at(f);
    if (p->phiOnFace.count(nullptr))
      return p->phiOnFace.at(nullptr);
    return std::get<0>(p->phiphin);
  };
  auto get_phin = [&](const networkPoint *p, networkFace *f) -> double {
    if (p->phinOnFace.count(f))
      return p->phinOnFace.at(f);
    if (p->phinOnFace.count(nullptr))
      return p->phinOnFace.at(nullptr);
    return std::get<1>(p->phiphin);
  };

  Tddd u = {0., 0., 0.};
  for (const auto &water : fluidObjects) {
    const double scale = water->getScale();
    for (const auto &integ_f : water->getBoundaryFaces()) {
      auto [q0, q1, q2] = integ_f->getPoints();

      networkPoint *closest_p = q0;
      {
        const double d0 = Norm(q0->X - a);
        const double d1 = Norm(q1->X - a);
        const double d2 = Norm(q2->X - a);
        if (d1 < d0 && d1 <= d2)
          closest_p = q1;
        else if (d2 < d0 && d2 < d1)
          closest_p = q2;
      }

      int how_far = 1;
      {
        const auto center = (q0->X + q1->X + q2->X) / 3.;
        if (Norm(center - a) < scale / near_ratio)
          how_far = 2; // 最も精度がいい__array_GW10xGW10__を近傍では使う
      }

      if (integ_f->isLinearElement) {
        if (integ_f->map_Point_LinearIntegrationInfo_vector.size() <= static_cast<std::size_t>(how_far) || !integ_f->map_Point_LinearIntegrationInfo_vector[how_far].count(closest_p)) {
          integ_f->setIntegrationInfo();
        }
        if (integ_f->map_Point_LinearIntegrationInfo_vector.size() <= static_cast<std::size_t>(how_far) || !integ_f->map_Point_LinearIntegrationInfo_vector[how_far].count(closest_p)) {
          continue;
        }

        const auto points3 = integ_f->getPoints(closest_p);
        const Tddd phi012 = {get_phi(points3[0], integ_f), get_phi(points3[1], integ_f), get_phi(points3[2], integ_f)};
        const Tddd phin012 = {get_phin(points3[0], integ_f), get_phin(points3[1], integ_f), get_phin(points3[2], integ_f)};

        for (const auto &[t0t1, ww, shape3, X, cross, J_det] : integ_f->map_Point_LinearIntegrationInfo_vector[how_far].at(closest_p)) {
          const auto r_vec = X - a;
          const double r = Norm(r_vec);
          if (r < scale * singular_eps) {
            // Very close to a quadrature point => treat `a` as (almost) on the boundary.
            // Do not early-return: keep integrating other faces, but remember the nearest face to fix the normal component later.
            const auto hit = interpolateVelocityAtNearestPointOnFace(integ_f);
            if (hit.hit && hit.dist2 < best_hit.dist2)
              best_hit = hit;
            continue;
          }

          const double r2 = r * r;
          const double r3 = r2 * r;

          const Tddd n = Normalize(cross);
          const double r_dot_n = Dot(r_vec, n);

          const Tddd gradG = r_vec;
          const Tddd grad_dG_dn = -n + 3. * r_dot_n * r_vec / r2;
          const double weight = J_det * ww;

          const double phi_x = Dot(shape3, phi012);
          const double phin_x = Dot(shape3, phin012);
          u += (phin_x * gradG + phi_x * grad_dG_dn) * weight / r3;
        }
      } else if (integ_f->isPseudoQuadraticElement) {
        if (integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector.size() <= static_cast<std::size_t>(how_far) || !integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector[how_far].count(closest_p) || !integ_f->map_Point_BEM_IGIGn_info_init.count(closest_p)) {
          integ_f->setIntegrationInfo();
        }
        if (integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector.size() <= static_cast<std::size_t>(how_far) || !integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector[how_far].count(closest_p) || !integ_f->map_Point_BEM_IGIGn_info_init.count(closest_p)) {
          continue;
        }

        const auto &key_list = integ_f->map_Point_BEM_IGIGn_info_init.at(closest_p); // 24 entries
        for (const auto &[t0t1, ww, N4, X, cross, J_det] : integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector[how_far].at(closest_p)) {
          const auto r_vec = X - a;
          const double r = Norm(r_vec);
          if (r < scale * singular_eps) {
            // Very close to a quadrature point => treat `a` as (almost) on the boundary.
            // Do not early-return: keep integrating other faces, but remember the nearest face to fix the normal component later.
            const auto hit = interpolateVelocityAtNearestPointOnFace(integ_f);
            if (hit.hit && hit.dist2 < best_hit.dist2)
              best_hit = hit;
            continue;
          }

          double phi_x = 0.;
          double phin_x = 0.;
          for (int group = 0; group < 4; ++group) {
            for (int i = 0; i < 6; ++i) {
              const auto &entry = key_list[static_cast<std::size_t>(group * 6 + i)];
              auto *p = std::get<0>(entry);
              auto *f_key = std::get<1>(entry);
              const double N = N4[static_cast<std::size_t>(group)][static_cast<std::size_t>(i)];
              phi_x += N * get_phi(p, f_key);
              phin_x += N * get_phin(p, f_key);
            }
          }

          const double r2 = r * r;
          const double r3 = r2 * r;

          const Tddd n = Normalize(cross);
          const double r_dot_n = Dot(r_vec, n);

          const Tddd gradG = r_vec;
          const Tddd grad_dG_dn = -n + 3. * r_dot_n * r_vec / r2;
          const double weight = J_det * ww;

          u += (phin_x * gradG + phi_x * grad_dG_dn) * weight / r3;
        }
      }
    }
  }

  // Base BEM integral result.
  Tddd u_int = u / (4. * M_PI);

  // If we detected a near-singular hit, enforce that the normal component matches the interpolated boundary velocity
  // on the nearest face, while preserving the tangential component from the integral (more stable near the boundary).
  if (best_hit.hit && Norm(best_hit.unitNormal) > 0.0) {
    const auto n = best_hit.unitNormal;
    const double un_target = Dot(best_hit.u_interp, n);
    const double un_now = Dot(u_int, n);
    u_int = u_int + n * (un_target - un_now);
  }

  return u_int;
}
