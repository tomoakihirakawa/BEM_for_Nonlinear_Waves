#pragma once

#include <algorithm>
#include <limits>
#include <vector>

#include "BEM_calculateVelocities.hpp"
#include "BEM_setBoundaryTypes.hpp"

namespace BEMMeshPipeline {

inline Tddd currentNormal(const networkPoint* p) {
  Tddd normal = {0., 0., 0.};
  for (const auto* f : p->getBoundaryFaces()) {
    if (!f)
      continue;
    normal += f->area * f->normal;
  }
  const double n = Norm(normal);
  return (n > 1e-20 && std::isfinite(n)) ? normal / n : Tddd{0., 0., 0.};
}

inline double currentShiftLimit(const networkPoint* p, double factor) {
  double length = std::numeric_limits<double>::infinity();
  for (const auto* t : p->getTetras()) {
    auto [p0, p1, p2, p3] = t->Points;
    const double avg = (Norm(p0->X - p1->X) + Norm(p1->X - p2->X) + Norm(p2->X - p0->X) +
                        Norm(p0->X - p3->X) + Norm(p1->X - p3->X) + Norm(p2->X - p3->X)) /
                       6.0;
    if (std::isfinite(avg) && avg > 0.)
      length = std::min(length, avg);
  }
  if (!std::isfinite(length)) {
    for (const auto* f : p->getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      const double avg = (Norm(p0->X - p1->X) + Norm(p1->X - p2->X) + Norm(p2->X - p0->X)) / 3.0;
      if (std::isfinite(avg) && avg > 0.)
        length = std::min(length, avg);
    }
  }
  return std::isfinite(length) ? factor * length : 0.;
}

inline double currentShiftLimit(const networkLine* l, double factor) {
  auto [pA, pB] = l->getPoints();
  const double length = Norm(pA->X - pB->X);
  return (std::isfinite(length) && length > 0.) ? factor * length : 0.;
}

template <class Entity>
inline Tddd vectorToCurrentSurface(Entity* entity, Tddd X_shifted) {
  const bool has_dirichlet_state = hasAnyDirichletBoundaryState(entity);
  const bool has_neumann_state = hasAnyNeumannBoundaryState(entity);
  const bool is_mixed_state = has_dirichlet_state && has_neumann_state;

  if (!has_dirichlet_state && !has_neumann_state)
    return {0., 0., 0.};

  const Tddd X_not_shifted = entity->getPosition();
  Tddd vecToDirichlet = {1E+20, 1E+20, 1E+20};
  networkFace* closest_face = nullptr;
  Tdd best_face_param = {0., 0.};

  auto faces = entity->getBoundaryFaces();
  if (faces.empty())
    return {0., 0., 0.};

  if (has_dirichlet_state) {
    for (const auto* f_const : faces) {
      auto* f = const_cast<networkFace*>(f_const);
      if (getNodeFaceBoundaryType(entity, f) != NodeFaceBoundaryType::Dirichlet)
        continue;
      Tdd face_param;
      const Tddd X = NearestOnDirichletFace(
          X_shifted, f, &face_param,
          [](const auto* q) -> Tddd { return q->getPosition(); });
      if (Norm(vecToDirichlet) >= Norm(X - X_shifted)) {
        vecToDirichlet = X - X_shifted;
        closest_face = f;
        best_face_param = face_param;
      }
    }
    if (closest_face != nullptr) {
      entity->relocation_face = closest_face;
      entity->relocation_param = best_face_param;
      if (!has_neumann_state)
        return vecToDirichlet;
      X_shifted += vecToDirichlet;
    } else {
      if (!has_neumann_state)
        return {0., 0., 0.};
      vecToDirichlet = {0., 0., 0.};
    }
  }

  struct DirectionInfo {
    double distance;
    Tddd direction;
  };
  std::vector<DirectionInfo> direction_infos;
  std::vector<T3Tddd> current_triangles;
  int isInContact_pass_count = 0;

  auto add_vector = [&](const Tddd& V, Tddd n) {
    if (!isFinite(V) || !isFinite(n) || !(Norm(n) > 0))
      return;
    n = Normalize(n);
    const double new_dist = Dot(V, n);
    std::erase_if(direction_infos, [&](const auto& info) {
      return Dot(n, info.direction) > std::cos(M_PI / 180. * 20.) &&
             std::abs(info.distance) > std::abs(new_dist);
    });
    if (std::ranges::none_of(direction_infos, [&](const auto& info) {
          return Dot(n, info.direction) > std::cos(M_PI / 180. * 20.);
        })) {
      direction_infos.push_back({new_dist, n});
    }
  };

  if (has_neumann_state) {
    const double short_range = 0.01;
    const bool has_contact = std::ranges::any_of(entity->dofs, [](const auto& kv) {
      return !kv.second.contact_opponent_faces.empty();
    });
    if (!has_contact && entity->penetratedBody != nullptr) {
      auto [f, X_nearest] = entity->penetratedBody->Nearest(X_shifted);
      if (f != nullptr)
        add_vector(X_nearest - X_shifted, Normalize(X_nearest - X_shifted));
    }

    auto body_faces = bfs(getEffectiveContactFaces(entity), 2);
    entity->debug_body_vertices_count = static_cast<int>(body_faces.size());
    for (const auto* fluid_face : faces) {
      if (getNodeFaceBoundaryType(entity, fluid_face) != NodeFaceBoundaryType::Neumann)
        continue;
      for (const auto* body_face : body_faces) {
        if (!body_face)
          continue;
        const T3Tddd tri = ToX(body_face);
        if (!isInContact(X_not_shifted, fluid_face->normal, tri, entity->contact_range))
          continue;
        ++isInContact_pass_count;
        current_triangles.emplace_back(tri);
        auto [t0, t1, X_nearest, n] = Nearest_(X_shifted, tri);
        (void)t0;
        (void)t1;
        const Tddd To_nearest = X_nearest - X_shifted;
        const double distance = Norm(To_nearest);
        if ((isFlat(n, To_nearest, contactAcceptanceAngle(distance, entity->contact_range)) ||
             isFlat(n, -To_nearest, contactAcceptanceAngle(distance, entity->contact_range))) &&
            (entity->contact_range >= distance))
          add_vector(To_nearest, n);
        else if (short_range * entity->contact_range >= distance)
          add_vector(To_nearest, n);
      }
    }
  }

  entity->debug_direction_info_count = static_cast<int>(direction_infos.size());
  entity->debug_contact_faces_count = static_cast<int>(std::ranges::count_if(entity->dofs, [](const auto& kv) {
    return !kv.second.contact_opponent_faces.empty();
  }));
  entity->debug_isInContact_pass_count = isInContact_pass_count;

  if (direction_infos.empty())
    return is_mixed_state ? vecToDirichlet : Tddd{0., 0., 0.};

  std::vector<double> distances;
  std::vector<Tddd> directions;
  distances.reserve(direction_infos.size());
  directions.reserve(direction_infos.size());
  for (const auto& info : direction_infos) {
    distances.push_back(info.distance);
    directions.push_back(info.direction);
  }

  Tddd V_opt = optimalVector(distances, directions, {0., 0., 0.});
  if (!current_triangles.empty()) {
    const Tddd X_target = X_shifted + V_opt;
    V_opt = Nearest(X_target, current_triangles) - X_shifted;
  }
  return V_opt + (is_mixed_state ? vecToDirichlet : Tddd{0., 0., 0.});
}

struct AdjustStats {
  int iterations_done = 0;
  double max_move_final = 0.0;
  std::size_t nodes_adjusted = 0;
  std::size_t lines_snapped = 0;
};

inline AdjustStats adjustMeshAtRemeshStage(const std::vector<Network*>& fluid_nets,
                                           int max_iter = 10,
                                           double tol = 1e-6,
                                           double move_limit_factor = 0.3,
                                           const std::filesystem::path* debug_dir = nullptr) {
  (void)debug_dir;
  AdjustStats stats;
  double global_mean = 0.;
  int global_count = 0;
  for (auto* net : fluid_nets) {
    for (auto* l : net->getBoundaryLines()) {
      auto [pA, pB] = l->getPoints();
      const double len = Norm(pA->X - pB->X);
      if (std::isfinite(len) && len > 0.) {
        global_mean += len;
        ++global_count;
      }
    }
  }
  global_mean = (global_count > 0) ? global_mean / static_cast<double>(global_count) : 1.0;

  for (int iter = 0; iter < max_iter; ++iter) {
    double max_move = 0.;
    for (auto* net : fluid_nets) {
      std::vector<std::pair<networkPoint*, Tddd>> moves;
      for (auto* p : net->getBoundaryPoints()) {
        if (!hasAnyNeumannBoundaryState(p))
          continue;
        Tddd vec = vectorToCurrentSurface(p, p->X);
        if (!isFinite(vec))
          continue;
        const double n = Norm(vec);
        const double limit = currentShiftLimit(p, move_limit_factor);
        if (!(n > 1e-12) || !(limit > 0.))
          continue;
        const Tddd delta = std::min(n, limit) * Normalize(vec);
        moves.emplace_back(p, delta);
      }
      for (auto& [p, delta] : moves) {
        p->setXSingle(p->X + delta);
        max_move = std::max(max_move, Norm(delta));
        ++stats.nodes_adjusted;
      }
      net->setGeometricPropertiesForce();
    }
    stats.iterations_done = iter + 1;
    stats.max_move_final = max_move;
    if (max_move < tol * global_mean) {
      std::cout << Green << "[mesh_pipeline] clung converged at iter=" << iter
                << " max_move=" << max_move
                << " nodes_adjusted=" << stats.nodes_adjusted << colorReset << std::endl;
      break;
    }
    if (iter + 1 == max_iter)
      std::cout << Yellow << "[mesh_pipeline] clung reached max_iter=" << max_iter
                << " max_move=" << max_move
                << " nodes_adjusted=" << stats.nodes_adjusted << colorReset << std::endl;
  }

  for (auto* net : fluid_nets) {
    int snapped = 0;
    for (auto* l : net->getBoundaryLines()) {
      if (!hasAnyDirichletBoundaryState(l) && !hasAnyNeumannBoundaryState(l))
        continue;
      auto [pA, pB] = l->getPoints();
      const Tddd linear_mid = 0.5 * (pA->X + pB->X);
      Tddd vec = vectorToCurrentSurface(l, linear_mid);
      if (!isFinite(vec))
        continue;
      const double n = Norm(vec);
      const double limit = currentShiftLimit(l, move_limit_factor);
      Tddd delta = {0., 0., 0.};
      if (n > 1e-12 && limit > 0.)
        delta = std::min(n, limit) * Normalize(vec);
      l->setXSingle(linear_mid + delta);
      ++snapped;
      ++stats.lines_snapped;
    }
    net->setGeometricPropertiesForce();
    std::cout << Green << "[mesh_pipeline] midpoint snap: " << net->getName()
              << " snapped=" << snapped << colorReset << std::endl;
  }
  return stats;
}

} // namespace BEMMeshPipeline
