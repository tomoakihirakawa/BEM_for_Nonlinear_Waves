#pragma once

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "BEM_geometry_projector.hpp"
#include "BEM_inputfile_reader.hpp"
#include "BEM_pre_bvp_consistency.hpp"
// Legacy Stage 2 v1 waterline path kept for default-compatible operation while
// geometry_projector_v2 remains opt-in.  New BCInterface projection logic lives in
// BEM_geometry_projector.hpp.
#include "BEM_waterline_geometry_provider.hpp"

namespace BEM_Penetration {
void throwIfStructurePenetrated(const std::vector<Network*>&, const std::vector<Network*>&,
                                int, const std::string&, bool);
} // namespace BEM_Penetration

namespace BEMMeshPipeline {
namespace GeometryRepair {

template <class Entity>
inline bool isRepairableGeometryEntity(const Entity* entity) {
  return entity && (hasAnyDirichletBoundaryState(entity) ||
                    hasAnyNeumannBoundaryState(entity) ||
                    GeometryProjectorDetail::isBCInterfaceEntity(entity));
}

struct GeometryRepairReport {
  AdjustStats stats;
  int repaired_points = 0;
  int repaired_lines = 0;
  int repaired_bcinterface_lines = 0;
  int quality_damaged_lines = 0;
  int quality_damaged_faces = 0;
  int penetration_rejected_before_repair = 0;
  std::vector<networkLine*> damaged_lines;
  std::vector<networkFace*> damaged_faces;

  bool hasDamagedTopology() const {
    return quality_damaged_lines > 0 || quality_damaged_faces > 0;
  }
};

template <class T>
inline bool addUniquePtr(std::vector<T*>& items, T* item) {
  if (!item)
    return false;
  if (std::find(items.begin(), items.end(), item) != items.end())
    return false;
  items.push_back(item);
  return true;
}

inline void recordDamagedLine(GeometryRepairReport& report, networkLine* l) {
  if (addUniquePtr(report.damaged_lines, l))
    report.quality_damaged_lines = static_cast<int>(report.damaged_lines.size());
  if (!l)
    return;
  for (auto* f : l->getBoundaryFaces())
    addUniquePtr(report.damaged_faces, f);
  report.quality_damaged_faces = static_cast<int>(report.damaged_faces.size());
}

inline bool addUniqueDofEntity(std::vector<BEM_DOF_Base*>& entities, BEM_DOF_Base* entity) {
  if (!entity)
    return false;
  if (std::find(entities.begin(), entities.end(), entity) != entities.end())
    return false;
  entities.push_back(entity);
  return true;
}

inline std::vector<BEM_DOF_Base*> collectRawContactCheckEntities(networkLine* l) {
  std::vector<BEM_DOF_Base*> entities;
  if (!l)
    return entities;
  addUniqueDofEntity(entities, l);
  auto [pA, pB] = l->getPoints();
  addUniqueDofEntity(entities, pA);
  addUniqueDofEntity(entities, pB);
  for (auto* f : l->getBoundaryFaces()) {
    if (!f)
      continue;
    auto [fp0, fp1, fp2] = f->getPoints();
    if (fp0 && fp0->BCInterface)
      addUniqueDofEntity(entities, fp0);
    if (fp1 && fp1->BCInterface)
      addUniqueDofEntity(entities, fp1);
    if (fp2 && fp2->BCInterface)
      addUniqueDofEntity(entities, fp2);
    auto [fl0, fl1, fl2] = f->getLines();
    if (fl0 && fl0->BCInterface)
      addUniqueDofEntity(entities, fl0);
    if (fl1 && fl1->BCInterface)
      addUniqueDofEntity(entities, fl1);
    if (fl2 && fl2->BCInterface)
      addUniqueDofEntity(entities, fl2);
  }
  return entities;
}

inline void refreshRawContactEntities(const std::vector<BEM_DOF_Base*>& entities,
                                      const std::vector<Network*>& contact_objects) {
  for (auto* entity : entities) {
    if (!entity)
      continue;
    entity->clearContactFaces();
    entity->addContactFaces(contact_objects, false);
  }
}

inline bool rawContactEntitiesValid(const std::vector<BEM_DOF_Base*>& entities) {
  for (const auto* entity : entities) {
    if (!entity || !entity->BCInterface)
      continue;
    if (entity->getContactFaces().empty())
      return false;
  }
  return true;
}

inline bool midpointQualityRegresses(const BEMPreBVP::WaterlineMidpointStats& pre,
                                     const BEMPreBVP::WaterlineMidpointStats& post) {
  constexpr double kMinAcceptedSubtriAngleDeg = 2.0;
  constexpr double kMaxAcceptedSubtriAspect = 60.0;
  if (pre.waterline_subtri_count <= 0)
    return false;
  if (post.waterline_subtri_count <= 0)
    return true;
  if (!std::isfinite(post.subtri_min_angle_deg) ||
      post.subtri_min_angle_deg < kMinAcceptedSubtriAngleDeg)
    return true;
  if (!std::isfinite(post.subtri_max_aspect) ||
      post.subtri_max_aspect > kMaxAcceptedSubtriAspect)
    return true;
  if (std::isfinite(pre.subtri_min_angle_deg) && pre.subtri_min_angle_deg > 0.0 &&
      std::isfinite(post.subtri_min_angle_deg) &&
      post.subtri_min_angle_deg < 0.95 * pre.subtri_min_angle_deg)
    return true;
  if (std::isfinite(pre.subtri_max_aspect) && pre.subtri_max_aspect > 0.0 &&
      std::isfinite(post.subtri_max_aspect) &&
      post.subtri_max_aspect > 1.05 * pre.subtri_max_aspect)
    return true;
  if (std::isfinite(pre.subtri_min_area) && pre.subtri_min_area > 0.0 &&
      std::isfinite(post.subtri_min_area) &&
      post.subtri_min_area < 0.5 * pre.subtri_min_area)
    return true;
  return false;
}

inline bool midpointQualityUnacceptable(const BEMPreBVP::WaterlineMidpointStats& post) {
  constexpr double kMinAcceptedSubtriAngleDeg = 2.0;
  constexpr double kMaxAcceptedSubtriAspect = 60.0;
  if (post.waterline_subtri_count <= 0)
    return true;
  if (!std::isfinite(post.subtri_min_angle_deg) ||
      post.subtri_min_angle_deg < kMinAcceptedSubtriAngleDeg)
    return true;
  if (!std::isfinite(post.subtri_max_aspect) ||
      post.subtri_max_aspect > kMaxAcceptedSubtriAspect)
    return true;
  if (!std::isfinite(post.subtri_min_area) ||
      post.subtri_min_area <= 0.0)
    return true;
  return false;
}

inline void countLegacyWaterlineFailure(BEMPreBVP::WaterlineRefreshSummary& stats,
                                        WaterlineStatus status) {
  switch (status) {
  case WaterlineStatus::no_dirichlet_face:
    ++stats.failed_no_dirichlet;
    break;
  case WaterlineStatus::no_body_candidate:
    ++stats.failed_no_body;
    break;
  case WaterlineStatus::body_gap_too_large:
    ++stats.failed_body_gap;
    break;
  case WaterlineStatus::no_move_possible:
    ++stats.failed_no_move;
    break;
  case WaterlineStatus::invalid_projection:
    ++stats.failed_invalid;
    break;
  case WaterlineStatus::skipped_not_corner:
    ++stats.skipped;
    break;
  default:
    break;
  }
}

inline void countProjectedWaterlineFailure(BEMPreBVP::WaterlineRefreshSummary& stats,
                                           GeometryTargetStatus status) {
  switch (status) {
  case GeometryTargetStatus::no_reference_surface:
    ++stats.failed_no_dirichlet;
    break;
  case GeometryTargetStatus::no_body_surface:
    ++stats.failed_no_body;
    break;
  case GeometryTargetStatus::body_gap_too_large:
    ++stats.failed_body_gap;
    break;
  case GeometryTargetStatus::no_move_possible:
    ++stats.failed_no_move;
    break;
  case GeometryTargetStatus::invalid_projection:
    ++stats.failed_invalid;
    break;
  case GeometryTargetStatus::not_applicable:
    ++stats.skipped;
    break;
  default:
    break;
  }
}

inline void printWaterlineSummary(const BEMPreBVP::WaterlineRefreshSummary& stats,
                                  const std::string& mode) {
  const int converged_total = stats.ok + stats.ok_via_fallback;
  const int usable_total = converged_total + stats.no_convergence;
  const double converged_ratio = (stats.lines_total > 0) ? static_cast<double>(converged_total) / stats.lines_total : 0.0;
  const double usable_ratio = (stats.lines_total > 0) ? static_cast<double>(usable_total) / stats.lines_total : 0.0;
  std::cout << Green << "[mesh:waterline]"
            << " mode=" << mode
            << " lines=" << stats.lines_total
            << " ok=" << stats.ok
            << " fallback=" << stats.ok_via_fallback
            << " no_convergence=" << stats.no_convergence
            << " converged_ratio=" << converged_ratio
            << " usable_ratio=" << usable_ratio
            << " quality_damaged=" << stats.quality_damaged_after_repair
            << " quality_rejected=" << stats.quality_rejected_before_repair
            << " raw_contact_rejected=" << stats.raw_contact_rejected_before_repair
            << " failed={no_dir:" << stats.failed_no_dirichlet
            << ",no_body:" << stats.failed_no_body
            << ",body_gap:" << stats.failed_body_gap
            << ",no_move:" << stats.failed_no_move
            << ",invalid:" << stats.failed_invalid
            << "} skipped=" << stats.skipped
            << " clamped=" << stats.clamped_count
            << " max_free_gap=" << stats.max_free_gap
            << " max_body_gap=" << stats.max_body_gap
            << " max_move_ratio=" << stats.max_move_ratio
            << colorReset << std::endl;
}

inline BEMPreBVP::WaterlineRefreshSummary repairWaterlineMidpointsLegacy(
    const std::vector<Network*>& fluid_nets,
    const SimulationSettings::TimeDomainSettings::MeshPreparationPipelineSettings& settings) {
  BEMPreBVP::WaterlineRefreshSummary stats;
  for (auto* water : fluid_nets) {
    if (!water)
      continue;
    for (auto* l : water->getBoundaryLines()) {
      if (!l)
        continue;
      if (!l->BCInterface) {
        ++stats.skipped;
        continue;
      }
      ++stats.lines_total;
      auto [pA, pB] = l->getPoints();
      if (!pA || !pB) {
        ++stats.failed_invalid;
        continue;
      }

      // Query X_mid from the current endpoint midpoint.  The free-surface side
      // of the projection is still the ReferenceState, not the already repaired
      // current surface; the body side is the current body geometry.
      const Tddd x_linear = 0.5 * (pA->X + pB->X);
      WaterlineQuery query;
      query.l = l;
      query.X_linear = x_linear;
      query.move_limit_factor = settings.clung_move_limit_factor;
      const auto result = queryWaterlineGeometry(query);

      if (!(result.status == WaterlineStatus::ok ||
            result.status == WaterlineStatus::ok_via_fallback ||
            result.status == WaterlineStatus::no_convergence)) {
        countLegacyWaterlineFailure(stats, result.status);
        continue;
      }
      if (!isFinite(result.target_X_clamped) || !isFinite(result.delta_clamped) ||
          !std::isfinite(result.body_gap) || !std::isfinite(result.body_gap_threshold) ||
          result.body_gap > result.body_gap_threshold) {
        ++stats.failed_invalid;
        continue;
      }

      const auto adjacent_faces = l->getBoundaryFaces();
      const auto pre_quality = BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces);
      const auto post_quality = BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces, l, result.target_X_clamped);
      if (midpointQualityRegresses(pre_quality, post_quality)) {
        ++stats.quality_rejected_before_repair;
        continue;
      }

      l->setXSingle(result.target_X_clamped);
      l->corner_midpoint_offset = l->X_mid - x_linear;

      if (result.status == WaterlineStatus::ok)
        ++stats.ok;
      else if (result.status == WaterlineStatus::ok_via_fallback)
        ++stats.ok_via_fallback;
      else
        ++stats.no_convergence;
      stats.max_free_gap = std::max(stats.max_free_gap, result.free_surface_gap);
      stats.max_body_gap = std::max(stats.max_body_gap, result.body_gap);
      stats.max_move_ratio = std::max(stats.max_move_ratio, result.move_ratio);
      if (result.move_ratio >= 1.0 - 1e-12)
        ++stats.clamped_count;
    }
  }
  printWaterlineSummary(stats, "legacy_waterline");
  BEMPreBVP::set_latest_waterline_refresh_summary(stats);
  return stats;
}

inline AdjustStats repairCurrentSurfaceGeometry(const std::vector<Network*>& fluid_nets,
                                                int max_iter = 10,
                                                double tol = 1e-6,
                                                double move_limit_factor = 0.3,
                                                const std::filesystem::path* debug_dir = nullptr) {
  return adjustMeshAtRemeshStage(fluid_nets, max_iter, tol, move_limit_factor, debug_dir);
}

inline bool isAcceptedRepairTarget(const GeometryTarget& target) {
  if (!(target.status == GeometryTargetStatus::ok ||
        target.status == GeometryTargetStatus::ok_via_fallback ||
        target.status == GeometryTargetStatus::no_convergence))
    return false;
  if (!isFinite(target.target_X_clamped) || !isFinite(target.delta_clamped) ||
      !std::isfinite(target.move_ratio))
    return false;
  if (std::isfinite(target.body_gap_threshold) && target.body_gap_threshold > 0.0 &&
      std::isfinite(target.body_gap) && target.body_gap > target.body_gap_threshold)
    return false;
  return true;
}

inline double meanBoundaryEdgeLength(const std::vector<Network*>& fluid_nets) {
  double sum = 0.0;
  int count = 0;
  for (auto* net : fluid_nets) {
    if (!net)
      continue;
    for (auto* l : net->getBoundaryLines()) {
      if (!l)
        continue;
      auto [pA, pB] = l->getPoints();
      if (!pA || !pB)
        continue;
      const double len = Norm(pA->X - pB->X);
      if (std::isfinite(len) && len > 0.0) {
        sum += len;
        ++count;
      }
    }
  }
  return count > 0 ? sum / static_cast<double>(count) : 1.0;
}

inline double repairPenetrationTolerance(const double h_local) {
  constexpr double c_rel = 0.25;
  constexpr double c_abs_eps = 1e-6;
  return (h_local > 0.0) ? c_rel * h_local + c_abs_eps : c_abs_eps;
}

inline std::vector<Network*> repairPenetrationCandidateSolids(const std::vector<Network*>& fluid_nets,
                                                              const std::vector<Network*>& body_objects) {
  std::vector<Network*> out;
  out.reserve(body_objects.size());
  for (auto* body : body_objects) {
    if (!body)
      continue;
    std::size_t inside_count = 0;
    std::size_t total_count = 0;
    for (const auto* water : fluid_nets) {
      if (!water)
        continue;
      for (const auto* p : water->getBoundaryPoints()) {
        ++total_count;
        if (body->InsideQ(p->X))
          ++inside_count;
      }
    }
    const double inside_fraction =
        total_count > 0 ? static_cast<double>(inside_count) / static_cast<double>(total_count) : 0.0;
    if (inside_fraction <= 0.2)
      out.push_back(body);
  }
  return out;
}

inline double pointPenetrationExcess(const Tddd& x,
                                     const std::vector<Network*>& active_bodies,
                                     const double h_local) {
  const double tolerance = repairPenetrationTolerance(h_local);
  double worst = 0.0;
  for (auto* body : active_bodies) {
    if (!body || !body->InsideQ(x))
      continue;
    auto [near_f, near_x] = body->Nearest(x);
    const double distance = near_f ? Norm(x - near_x) : 0.0;
    if (std::isfinite(distance))
      worst = std::max(worst, distance - tolerance);
  }
  return std::max(0.0, worst);
}

inline GeometryProjectorOptions makeProjectedGeometryRepairOptions(
    const SimulationSettings::TimeDomainSettings::MeshPreparationPipelineSettings& settings,
    GeometryProjectionMode projection_mode) {
  GeometryProjectorOptions options;
  options.mode = projection_mode;
  options.move_limit_factor = settings.clung_move_limit_factor;
  options.max_iter = 4;
  options.tol_relative = 1e-6;
  options.body_bfs_fallback_depth = 1;
  options.enable_global_body_fallback = true;
  options.global_body_fallback_range_factor = 1.0;
  return options;
}

struct ProjectedGeometryRepairContext {
  const std::vector<Network*>& fluid_nets;
  const ReferenceState& reference;
  const std::vector<Network*>& body_objects;
  const SimulationSettings::TimeDomainSettings::MeshPreparationPipelineSettings& settings;
  GeometryProjectionMode projection_mode;
  GeometryProjectorOptions options;
  BoundaryGeometryTargetProvider target_provider;
  std::vector<Network*> active_penetration_bodies;
  GeometryRepairReport report;
  BEMPreBVP::WaterlineRefreshSummary waterline_stats;
  int warn_committed_without_line_contact = 0;
  int warn_committed_global_body = 0;

  ProjectedGeometryRepairContext(
      const std::vector<Network*>& fluid_nets_,
      const ReferenceState& reference_,
      const std::vector<Network*>& body_objects_,
      const SimulationSettings::TimeDomainSettings::MeshPreparationPipelineSettings& settings_,
      GeometryProjectionMode projection_mode_)
      : fluid_nets(fluid_nets_),
        reference(reference_),
        body_objects(body_objects_),
        settings(settings_),
        projection_mode(projection_mode_),
        options(makeProjectedGeometryRepairOptions(settings_, projection_mode_)),
        target_provider(reference_, body_objects_, options),
        active_penetration_bodies(repairPenetrationCandidateSolids(fluid_nets_, body_objects_)) {
    waterline_stats.valid = true;
  }
};

inline void repairVertexNodePositions(ProjectedGeometryRepairContext& ctx) {
  auto& report = ctx.report;
  auto& stats = report.stats;
  const double global_mean = meanBoundaryEdgeLength(ctx.fluid_nets);
  const int max_iter = std::max(1, ctx.settings.clung_max_iter);
  for (int iter = 0; iter < max_iter; ++iter) {
    double max_move = 0.0;
    for (auto* water : ctx.fluid_nets) {
      if (!water)
        continue;
      std::vector<std::pair<networkPoint*, Tddd>> moves;
      for (auto* p : water->getBoundaryPoints()) {
        if (!isRepairableGeometryEntity(p))
          continue;
        const auto target = ctx.target_provider.queryEntityTarget(p, water, p->X);
        if (!isAcceptedRepairTarget(target))
          continue;
        const double target_move = Norm(target.delta_clamped);
        if (!(target_move > 1e-12))
          continue;
        const double h_local = localEdgeLength(p);
        const double pre_excess = pointPenetrationExcess(p->X, ctx.active_penetration_bodies, h_local);
        const auto adjacent_faces = p->getBoundaryFaces();
        const auto waterline_adjacent_faces = BEMPreBVP::affected_waterline_faces(p);
        const auto local_pre_quality =
            BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces);
        const bool waterline_sensitive = local_pre_quality.waterline_subtri_count > 0;
        const auto full_quality_baseline =
            waterline_sensitive ? BEMPreBVP::make_waterline_quality_baseline(water->getBoundaryFaces())
                                : BEMPreBVP::WaterlineQualityBaseline{};
        const auto full_pre_quality =
            waterline_sensitive ? full_quality_baseline.full_stats
                                : BEMPreBVP::WaterlineMidpointStats{};
        bool rejected_by_penetration = false;
        bool rejected_by_waterline_quality = false;
        bool accepted = false;
        Tddd accepted_delta = {0., 0., 0.};
        for (const double alpha : {1.0, 0.5, 0.25, 0.125, 0.0}) {
          const Tddd candidate_delta = alpha * target.delta_clamped;
          const Tddd candidate_x = p->X + candidate_delta;
          const double post_excess = pointPenetrationExcess(candidate_x, ctx.active_penetration_bodies, h_local);
          const bool penetration_ok = pre_excess > 0.0
                                          ? post_excess < pre_excess - 1e-12
                                          : post_excess <= pre_excess + 1e-12;
          if (!penetration_ok) {
            rejected_by_penetration = true;
            continue;
          }
          if (waterline_sensitive) {
            const auto local_post_quality =
                BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces, nullptr, std::nullopt, p, candidate_x);
            const auto full_post_quality =
                ctx.settings.waterline_fast_quality_enabled
                    ? BEMPreBVP::measure_waterline_midpoint_quality_delta(
                          full_quality_baseline, waterline_adjacent_faces, nullptr, std::nullopt, p, candidate_x)
                    : BEMPreBVP::measure_waterline_midpoint_quality(water->getBoundaryFaces(), nullptr, std::nullopt, p, candidate_x);
            const bool local_quality_ok =
                !midpointQualityUnacceptable(local_post_quality) &&
                !midpointQualityRegresses(local_pre_quality, local_post_quality);
            const bool full_quality_ok =
                !midpointQualityUnacceptable(full_post_quality) &&
                !midpointQualityRegresses(full_pre_quality, full_post_quality);
            if (!local_quality_ok || !full_quality_ok) {
              rejected_by_waterline_quality = true;
              continue;
            }
          }
          accepted = true;
          accepted_delta = candidate_delta;
          break;
        }
        if (rejected_by_penetration)
          ++report.penetration_rejected_before_repair;
        if (rejected_by_waterline_quality) {
          for (auto* f : adjacent_faces)
            addUniquePtr(report.damaged_faces, f);
          report.quality_damaged_faces = static_cast<int>(report.damaged_faces.size());
        }
        if (!accepted)
          continue;
        const double move = Norm(accepted_delta);
        if (!(move > 1e-12))
          continue;
        moves.emplace_back(p, accepted_delta);
      }
      for (auto& [p, delta] : moves) {
        p->setXSingle(p->X + delta);
        max_move = std::max(max_move, Norm(delta));
        ++stats.nodes_adjusted;
        ++report.repaired_points;
      }
      water->setGeometricPropertiesForce();
    }
    stats.iterations_done = iter + 1;
    stats.max_move_final = max_move;
    if (max_move < ctx.settings.clung_tol * global_mean)
      break;
  }
}

inline void rebuildContactAndBoundaryStateAfterVertexRepair(ProjectedGeometryRepairContext& ctx) {
  std::vector<Network*> all_objects;
  all_objects.reserve(ctx.fluid_nets.size() + ctx.body_objects.size());
  for (auto* water : ctx.fluid_nets) {
    if (water)
      all_objects.push_back(water);
  }
  for (auto* body : ctx.body_objects) {
    if (body && std::find(all_objects.begin(), all_objects.end(), body) == all_objects.end())
      all_objects.push_back(body);
  }

  for (auto* net : all_objects) {
    if (net)
      net->setGeometricPropertiesForce();
  }
  _Pragma("omp parallel for") for (const auto& net : all_objects) {
    if (net)
      net->makeBuckets(net->getScale() / 10.);
  }

  auto& mutable_fluid_nets = const_cast<std::vector<Network*>&>(ctx.fluid_nets);
  refreshBoundaryStatesAndTypes(mutable_fluid_nets, ctx.body_objects);
  for (auto* water : ctx.fluid_nets) {
    if (water)
      computeAllBCInterfaceMidpointOffsets(water);
  }
}

inline bool tryRepairOrdinaryEdgeNodePosition(ProjectedGeometryRepairContext& ctx,
                                              Network* water,
                                              networkLine* l,
                                              int& damaged_log_count) {
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB)
    return false;
  const Tddd x_linear = 0.5 * (pA->X + pB->X);
  const auto target = ctx.target_provider.queryEntityTarget(l, water, x_linear);
  if (!isAcceptedRepairTarget(target))
    return false;
  if (l->BCInterface) {
    const auto adjacent_faces = l->getBoundaryFaces();
    const auto pre_quality = BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces);
    const auto post_quality = BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces, l, target.target_X_clamped);
    if (midpointQualityRegresses(pre_quality, post_quality)) {
      recordDamagedLine(ctx.report, l);
      if (ctx.settings.verbose_debug && damaged_log_count < 16) {
        ++damaged_log_count;
        std::cout << "[mesh:geometry:damaged]"
                  << " line=" << l
                  << " status=" << toString(target.status)
                  << " fallback=" << toString(target.fallback)
                  << " body_gap=" << target.body_gap
                  << " body_gap_threshold=" << target.body_gap_threshold
                  << " dirichlet_gap=" << target.dirichlet_gap
                  << " move_ratio=" << target.move_ratio
                  << " pre={count:" << pre_quality.waterline_subtri_count
                  << ",min_angle:" << pre_quality.subtri_min_angle_deg
                  << ",max_aspect:" << pre_quality.subtri_max_aspect
                  << ",min_area:" << pre_quality.subtri_min_area
                  << "} post={count:" << post_quality.waterline_subtri_count
                  << ",min_angle:" << post_quality.subtri_min_angle_deg
                  << ",max_aspect:" << post_quality.subtri_max_aspect
                  << ",min_area:" << post_quality.subtri_min_area
                  << "} x_linear={" << x_linear[0] << "," << x_linear[1] << "," << x_linear[2] << "}"
                  << " target={" << target.target_X_clamped[0] << "," << target.target_X_clamped[1] << "," << target.target_X_clamped[2] << "}"
                  << std::endl;
      }
    }
  }
  Tddd accepted_midpoint = target.target_X_clamped;
  const auto waterline_adjacent_faces = l->getBoundaryFaces();
  const auto local_pre_quality =
      BEMPreBVP::measure_waterline_midpoint_quality(waterline_adjacent_faces);
  if (local_pre_quality.waterline_subtri_count > 0) {
    const auto affected_faces = BEMPreBVP::affected_waterline_faces(l);
    const auto full_quality_baseline =
        BEMPreBVP::make_waterline_quality_baseline(water->getBoundaryFaces());
    const auto full_pre_quality = full_quality_baseline.full_stats;
    bool accepted = false;
    for (const double alpha : {1.0, 0.5, 0.25, 0.125, 0.0}) {
      const Tddd candidate_midpoint = x_linear + alpha * (target.target_X_clamped - x_linear);
      const auto local_post_quality =
          BEMPreBVP::measure_waterline_midpoint_quality(waterline_adjacent_faces, l, candidate_midpoint);
      const auto full_post_quality =
          ctx.settings.waterline_fast_quality_enabled
              ? BEMPreBVP::measure_waterline_midpoint_quality_delta(
                    full_quality_baseline, affected_faces, l, candidate_midpoint)
              : BEMPreBVP::measure_waterline_midpoint_quality(water->getBoundaryFaces(), l, candidate_midpoint);
      const bool local_quality_ok =
          !midpointQualityUnacceptable(local_post_quality) &&
          !midpointQualityRegresses(local_pre_quality, local_post_quality);
      const bool full_quality_ok =
          !midpointQualityUnacceptable(full_post_quality) &&
          !midpointQualityRegresses(full_pre_quality, full_post_quality);
      if (!local_quality_ok || !full_quality_ok)
        continue;
      accepted_midpoint = candidate_midpoint;
      accepted = true;
      break;
    }
    if (!accepted) {
      recordDamagedLine(ctx.report, l);
      return false;
    }
  }
  l->setXSingle(accepted_midpoint);
  l->corner_midpoint_offset = l->X_mid - x_linear;
  ++ctx.report.stats.lines_snapped;
  ++ctx.report.repaired_lines;
  return true;
}

inline bool tryRepairBCInterfaceEdgeNodePosition(ProjectedGeometryRepairContext& ctx,
                                                 Network* water,
                                                 networkLine* l) {
  auto& stats = ctx.waterline_stats;
  ++stats.lines_total;
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB) {
    ++stats.failed_invalid;
    return false;
  }

  const Tddd x_linear = 0.5 * (pA->X + pB->X);
  const bool has_line_contact = !getEffectiveContactFaces(l).empty();
  if (!has_line_contact) {
    ++ctx.warn_committed_without_line_contact;
    if (ctx.settings.verbose_debug && ctx.warn_committed_without_line_contact <= 8) {
      std::cout << "[waterline_repair:no_line_contact_fallback] line=" << l
                << " x_linear={" << x_linear[0] << "," << x_linear[1] << "," << x_linear[2] << "}"
                << std::endl;
    }
  }
  const auto target = ctx.target_provider.queryEntityTarget(l, water, x_linear);

  if (!(target.status == GeometryTargetStatus::ok ||
        target.status == GeometryTargetStatus::ok_via_fallback ||
        target.status == GeometryTargetStatus::no_convergence)) {
    countProjectedWaterlineFailure(stats, target.status);
    return false;
  }
  if (!isFinite(target.target_X_clamped) || !isFinite(target.delta_clamped) ||
      !std::isfinite(target.body_gap) || !std::isfinite(target.body_gap_threshold) ||
      target.body_gap > target.body_gap_threshold) {
    ++stats.failed_invalid;
    return false;
  }

  const auto adjacent_faces = l->getBoundaryFaces();
  const auto waterline_adjacent_faces = BEMPreBVP::affected_waterline_faces(l);
  const auto pre_quality = BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces);
  const auto full_quality_baseline =
      BEMPreBVP::make_waterline_quality_baseline(water->getBoundaryFaces());
  const auto full_pre_quality = full_quality_baseline.full_stats;
  const auto original_midpoint = l->X_mid;
  const auto original_offset = l->corner_midpoint_offset;
  const auto contact_check_entities = collectRawContactCheckEntities(l);

  struct Candidate {
    double alpha = 0.0;
    Tddd midpoint = {0., 0., 0.};
    double move_ratio = 0.0;
    BEMPreBVP::WaterlineMidpointStats quality;
    BEMPreBVP::WaterlineMidpointStats full_quality;
    bool quality_ok = false;
    bool full_quality_ok = false;
    bool raw_contact_ok = false;
    bool accepted = false;
  };

  std::vector<Candidate> candidates;
  Candidate accepted_candidate;
  bool rejected_by_quality = false;
  bool rejected_by_full_quality = false;
  bool rejected_by_raw_contact = false;

  for (const double alpha : {1.0, 0.5, 0.25, 0.125, 0.0}) {
    Candidate candidate;
    candidate.alpha = alpha;
    candidate.midpoint = x_linear + alpha * (target.target_X_clamped - x_linear);
    candidate.move_ratio = alpha * target.move_ratio;
    candidate.quality =
        BEMPreBVP::measure_waterline_midpoint_quality(adjacent_faces, l, candidate.midpoint);
    candidate.full_quality =
        ctx.settings.waterline_fast_quality_enabled
            ? BEMPreBVP::measure_waterline_midpoint_quality_delta(
                  full_quality_baseline, waterline_adjacent_faces, l, candidate.midpoint)
            : BEMPreBVP::measure_waterline_midpoint_quality(water->getBoundaryFaces(), l, candidate.midpoint);
    candidate.quality_ok =
        !midpointQualityUnacceptable(candidate.quality) &&
        !midpointQualityRegresses(pre_quality, candidate.quality);
    candidate.full_quality_ok =
        !midpointQualityUnacceptable(candidate.full_quality) &&
        !midpointQualityRegresses(full_pre_quality, candidate.full_quality);
    if (!candidate.quality_ok || !candidate.full_quality_ok) {
      rejected_by_quality = true;
      if (!candidate.full_quality_ok)
        rejected_by_full_quality = true;
      candidates.push_back(candidate);
      continue;
    }

    l->setXSingle(candidate.midpoint);
    refreshRawContactEntities(contact_check_entities, ctx.body_objects);
    candidate.raw_contact_ok = rawContactEntitiesValid(contact_check_entities);
    if (!candidate.raw_contact_ok) {
      rejected_by_raw_contact = true;
      candidates.push_back(candidate);
      continue;
    }

    candidate.accepted = true;
    accepted_candidate = candidate;
    candidates.push_back(candidate);
    break;
  }

  if (!accepted_candidate.accepted) {
    l->setXSingle(original_midpoint);
    l->corner_midpoint_offset = original_offset;
    refreshRawContactEntities(contact_check_entities, ctx.body_objects);
    if (rejected_by_quality)
      ++stats.quality_rejected_before_repair;
    if (rejected_by_raw_contact)
      ++stats.raw_contact_rejected_before_repair;
    recordDamagedLine(ctx.report, l);
    const auto& last_quality = candidates.empty() ? pre_quality : candidates.back().quality;
    if (ctx.settings.verbose_debug) {
      std::cout << "[waterline_repair:rejected]"
                << " line=" << l
                << " status=" << toString(target.status)
                << " fallback=" << toString(target.fallback)
                << " has_line_contact=" << (has_line_contact ? 1 : 0)
                << " rejected_by_quality=" << (rejected_by_quality ? 1 : 0)
                << " rejected_by_full_quality=" << (rejected_by_full_quality ? 1 : 0)
                << " rejected_by_raw_contact=" << (rejected_by_raw_contact ? 1 : 0)
                << " body_gap=" << target.body_gap
                << " body_gap_threshold=" << target.body_gap_threshold
                << " dirichlet_gap=" << target.dirichlet_gap
                << " move_ratio=" << target.move_ratio
                << " pre={count:" << pre_quality.waterline_subtri_count
                << ",min_angle:" << pre_quality.subtri_min_angle_deg
                << ",max_aspect:" << pre_quality.subtri_max_aspect
                << ",min_area:" << pre_quality.subtri_min_area
                << "} last={count:" << last_quality.waterline_subtri_count
                << ",min_angle:" << last_quality.subtri_min_angle_deg
                << ",max_aspect:" << last_quality.subtri_max_aspect
                << ",min_area:" << last_quality.subtri_min_area
                << "} full_pre={count:" << full_pre_quality.waterline_subtri_count
                << ",min_angle:" << full_pre_quality.subtri_min_angle_deg
                << ",max_aspect:" << full_pre_quality.subtri_max_aspect
                << "} x_linear={" << x_linear[0] << "," << x_linear[1] << "," << x_linear[2] << "}"
                << " target={" << target.target_X_clamped[0] << "," << target.target_X_clamped[1] << "," << target.target_X_clamped[2] << "}"
                << std::endl;
    }
    return false;
  }

  auto accepted_midpoint = accepted_candidate.midpoint;
  double accepted_move_ratio = accepted_candidate.move_ratio;
  auto post_quality = accepted_candidate.quality;
  const auto full_post_quality = accepted_candidate.full_quality;
  if (accepted_candidate.alpha < 1.0) {
    if (rejected_by_quality)
      ++stats.quality_rejected_before_repair;
    if (rejected_by_raw_contact)
      ++stats.raw_contact_rejected_before_repair;
    if (ctx.settings.verbose_debug) {
      std::cout << "[waterline_repair:damped]"
                << " line=" << l
                << " alpha=" << accepted_candidate.alpha
                << " rejected_by_quality=" << (rejected_by_quality ? 1 : 0)
                << " rejected_by_full_quality=" << (rejected_by_full_quality ? 1 : 0)
                << " rejected_by_raw_contact=" << (rejected_by_raw_contact ? 1 : 0)
                << " pre={count:" << pre_quality.waterline_subtri_count
                << ",min_angle:" << pre_quality.subtri_min_angle_deg
                << ",max_aspect:" << pre_quality.subtri_max_aspect
                << ",min_area:" << pre_quality.subtri_min_area
                << "} full={count:" << full_post_quality.waterline_subtri_count
                << ",min_angle:" << full_post_quality.subtri_min_angle_deg
                << ",max_aspect:" << full_post_quality.subtri_max_aspect
                << ",min_area:" << full_post_quality.subtri_min_area
                << "}"
                << " move_ratio=" << target.move_ratio
                << std::endl;
    }
  }
  if (target.fallback == GeometryFallbackLevel::global_body) {
    ++ctx.warn_committed_global_body;
    if (ctx.settings.verbose_debug && ctx.warn_committed_global_body <= 8) {
      std::cout << "[waterline_repair:global_body_commit] line=" << l
                << " has_line_contact=" << (has_line_contact ? 1 : 0)
                << " body_gap=" << target.body_gap
                << " threshold=" << target.body_gap_threshold
                << std::endl;
    }
  }

  l->setXSingle(accepted_midpoint);
  l->corner_midpoint_offset = l->X_mid - x_linear;
  ++ctx.report.repaired_lines;
  ++ctx.report.repaired_bcinterface_lines;
  ++ctx.report.stats.lines_snapped;

  if (target.status == GeometryTargetStatus::ok)
    ++stats.ok;
  else if (target.status == GeometryTargetStatus::ok_via_fallback)
    ++stats.ok_via_fallback;
  else
    ++stats.no_convergence;
  stats.max_free_gap = std::max(stats.max_free_gap, target.dirichlet_gap);
  stats.max_body_gap = std::max(stats.max_body_gap, target.body_gap);
  stats.max_move_ratio = std::max(stats.max_move_ratio, accepted_move_ratio);
  if (accepted_move_ratio >= 1.0 - 1e-12) {
    ++stats.clamped_count;
    BEMPreBVP::WaterlineClampedEntity entity;
    entity.line_index = l->midpoint_index;
    entity.endpoint0_index = pA->index;
    entity.endpoint1_index = pB->index;
    entity.reference_face_index = -1;
    entity.body_face_index = target.body_face ? target.body_face->index : -1;
    entity.free_gap = target.dirichlet_gap;
    entity.body_gap = target.body_gap;
    entity.move_ratio = accepted_move_ratio;
    entity.local_min_angle_before = pre_quality.subtri_min_angle_deg;
    entity.local_min_angle_after = post_quality.subtri_min_angle_deg;
    entity.local_max_aspect_before = pre_quality.subtri_max_aspect;
    entity.local_max_aspect_after = post_quality.subtri_max_aspect;
    entity.midpoint = accepted_midpoint;
    stats.clamped_entities.emplace_back(entity);
  }
  return true;
}

inline void repairEdgeNodePositions(ProjectedGeometryRepairContext& ctx) {
  for (auto* water : ctx.fluid_nets) {
    if (!water)
      continue;
    int snapped = 0;
    int damaged_log_count = 0;
    for (auto* l : water->getBoundaryLines()) {
      if (!l)
        continue;
      if (ctx.settings.waterline_geometry && !l->BCInterface)
        ++ctx.waterline_stats.skipped;
      if (!isRepairableGeometryEntity(l))
        continue;

      const bool repaired = (ctx.settings.waterline_geometry && l->BCInterface)
                                ? tryRepairBCInterfaceEdgeNodePosition(ctx, water, l)
                                : tryRepairOrdinaryEdgeNodePosition(ctx, water, l, damaged_log_count);
      if (repaired)
        ++snapped;
    }
    water->setGeometricPropertiesForce();
    std::cout << Green << "[mesh:geometry] repair name=" << water->getName()
              << " snapped=" << snapped << colorReset << std::endl;
  }

  if (ctx.settings.waterline_geometry) {
    if (ctx.settings.verbose_debug &&
        (ctx.warn_committed_without_line_contact > 0 || ctx.warn_committed_global_body > 0)) {
      std::cout << "[waterline_repair:contact_summary]"
                << " committed_without_line_contact=" << ctx.warn_committed_without_line_contact
                << " committed_global_body=" << ctx.warn_committed_global_body
                << std::endl;
    }
    printWaterlineSummary(ctx.waterline_stats, "geometry_projector_v2");
    BEMPreBVP::set_latest_waterline_refresh_summary(ctx.waterline_stats);
  }
}

inline GeometryRepairReport finalizeProjectedGeometryRepair(
    std::vector<Network*>& fluid_nets,
    std::vector<Network*>& all_objects,
    const std::vector<Network*>& contact_objects,
    int step,
    const std::string& label,
    bool true_quadratic) {
  GeometryRepairReport report;
  BEM_Penetration::throwIfStructurePenetrated(fluid_nets, contact_objects, step, label, true_quadratic);
  for (auto* net : all_objects) {
    if (net)
      net->setGeometricPropertiesForce();
  }
  _Pragma("omp parallel for") for (const auto& net : all_objects) {
    if (net)
      net->makeBuckets(net->getScale() / 10.);
  }
  refreshBoundaryStatesAndTypes(fluid_nets, contact_objects);
  return report;
}

inline GeometryRepairReport repairCurrentSurfaceGeometryWithProjector(
    const std::vector<Network*>& fluid_nets,
    const ReferenceState& reference,
    const std::vector<Network*>& body_objects,
    const SimulationSettings::TimeDomainSettings::MeshPreparationPipelineSettings& settings,
    GeometryProjectionMode projection_mode) {
  ProjectedGeometryRepairContext ctx(fluid_nets, reference, body_objects, settings, projection_mode);
  repairVertexNodePositions(ctx);
  rebuildContactAndBoundaryStateAfterVertexRepair(ctx);
  repairEdgeNodePositions(ctx);
  return ctx.report;
}

inline GeometryRepairReport repairWaterlineMidpointsIfEnabled(
    std::vector<Network*>& fluid_nets,
    std::vector<Network*>& all_objects,
    const std::vector<Network*>& contact_objects,
    const ReferenceState& reference,
    const SimulationSettings::TimeDomainSettings::MeshPreparationPipelineSettings& settings,
    GeometryProjectionMode projection_mode,
    int step,
    const std::string& label,
    bool true_quadratic) {
  GeometryRepairReport report;
  if (!settings.waterline_geometry)
    return report;
  if (settings.geometry_projector_v2) {
    return finalizeProjectedGeometryRepair(
        fluid_nets, all_objects, contact_objects, step, label, true_quadratic);
  }
  for (auto* water : fluid_nets)
    computeAllBCInterfaceMidpointOffsets(water);
  (void)reference;
  (void)projection_mode;
  (void)repairWaterlineMidpointsLegacy(fluid_nets, settings);
  BEM_Penetration::throwIfStructurePenetrated(fluid_nets, contact_objects, step, label, true_quadratic);
  for (auto* net : all_objects)
    net->setGeometricPropertiesForce();
  _Pragma("omp parallel for") for (const auto& net : all_objects) net->makeBuckets(net->getScale() / 10.);
  refreshBoundaryStatesAndTypes(fluid_nets, contact_objects);
  return report;
}

} // namespace GeometryRepair
} // namespace BEMMeshPipeline
