// Separate translation unit for remesh_for_main_loop — the main remeshing function.
// Compiled independently to reduce main_time_domain.cpp build time.

#define BEM
#include "Network.hpp"
#include "BEM_remesh_main.hpp"

void remesh_for_main_loop(Network& water, int time_step, double min_edge_length, bool tetrahedralize, bool surface_flip,
                          const CollisionSettings& collision_settings,
                          bool surface_split, bool surface_collapse,
                          bool skip_post_remesh_quality_rejects) {
  const double rad = M_PI / 180.0;
  const double cos_3rad = std::cos(3.0 * rad);
  const double cos_rad = std::cos(rad);
  const double global_mean_len = Mean(extLength(water.getLines()));
  const double limit_len = (min_edge_length > 0.0) ? min_edge_length : global_mean_len * 0.1;
  constexpr double min_local_face_area_ratio = 0.05;

  // Safety: ensure X_mid is at least at the linear midpoint for all boundary lines.
  for (auto* l : water.getBoundaryLines()) {
    if (l->X_mid[0] == 0. && l->X_mid[1] == 0. && l->X_mid[2] == 0.) {
      auto [pA, pB] = l->getPoints();
      const Tddd mid = 0.5 * (pA->X + pB->X);
      if (mid[0] != 0. || mid[1] != 0. || mid[2] != 0.) {
        l->setXSingle(mid);
      }
    }
  }

  std::cout << "孤立した四面体の削除を開始" << std::endl;
  water.DeleteIsolatedTetras();
  std::cout << "内部要素の削除を開始" << std::endl;
  water.DeleteInteriorTetras();

  std::unordered_set<networkLine*> protected_lines;
  bool collision_ok = detectAndResolveCollisions(water, time_step, collision_settings, protected_lines);
  if (!collision_ok) {
    constexpr double max_fold_ratio = 0.15;
    constexpr double max_protected_ratio = 0.5;
    auto folded = detectFoldedFaces(water, collision_settings.normal_reversal_cos, global_mean_len);
    size_t n_boundary_faces = water.getBoundaryFaces().size();
    size_t n_boundary_lines = water.getBoundaryLines().size();
    double fold_ratio = (n_boundary_faces > 0) ? static_cast<double>(folded.size()) / n_boundary_faces : 0.0;
    double protected_ratio = (n_boundary_lines > 0) ? static_cast<double>(protected_lines.size()) / n_boundary_lines : 0.0;
    if (fold_ratio > max_fold_ratio)
      throw step_failure("collision unresolved + fold ratio " + std::to_string(fold_ratio) + " > " + std::to_string(max_fold_ratio) + " at time_step " + std::to_string(time_step));
    if (protected_ratio > max_protected_ratio)
      throw step_failure("collision unresolved + protected ratio " + std::to_string(protected_ratio) + " > " + std::to_string(max_protected_ratio) + " (" + std::to_string(protected_lines.size()) + " / " + std::to_string(n_boundary_lines) + " lines)" + " at time_step " + std::to_string(time_step));
    std::cout << Yellow << "[remesh] time_step " << time_step
              << ": collision unresolved (fold_ratio=" << fold_ratio
              << ", protected_ratio=" << protected_ratio << "), continuing"
              << colorReset << std::endl;
  }

  auto build_protection_halo = [&](const std::unordered_set<networkLine*>& seeds) {
    std::unordered_set<networkLine*> halo = seeds;
    std::vector<networkLine*> queue(seeds.begin(), seeds.end());
    for (auto* l : queue) {
      if (!l) continue;
      auto [p0, p1] = l->getPoints();
      auto add_from_point = [&](networkPoint* p) {
        if (!p) return;
        for (auto* nl : p->getBoundaryLines())
          if (nl) halo.insert(nl);
      };
      add_from_point(p0);
      add_from_point(p1);
      for (auto* f : l->getBoundaryFaces())
        if (f)
          for (auto* nl : f->getLines())
            if (nl) halo.insert(nl);
    }
    return halo;
  };
  const auto protected_halo_lines = build_protection_halo(protected_lines);

  water.setGeometricPropertiesForce();
  water.checkConnectivity();

  auto is_near_protected = [&](networkLine* l) {
    if (!l) return false;
    return protected_halo_lines.count(l) > 0;
  };

  const std::size_t boundary_line_count = water.getBoundaryLines().size();
  if (boundary_line_count > 0 && protected_lines.size() > 0) {
    std::cout << Yellow << "[remesh] time_step " << time_step
              << ": protected region is " << protected_lines.size() << " / " << boundary_line_count
              << " boundary lines (" << protected_halo_lines.size() << " incl. 1-ring halo)"
              << colorReset << std::endl;
  }
  const bool heavy_collision_protection = false;

  {
    int n_pre_topo_err = 0;
    for (auto* l : water.getBoundaryLines())
      if (!l->checkTopology()) n_pre_topo_err++;
    if (n_pre_topo_err > 0) {
      std::cerr << Red << "[remesh] time_step " << time_step
                << ": " << n_pre_topo_err << " pre-existing topology errors detected"
                << colorReset << std::endl;
      throw step_failure("pre-existing topology errors (" + std::to_string(n_pre_topo_err) + " lines) at time_step " + std::to_string(time_step));
    }
  }

  const int iter_divide_collapse = 3;
  for (auto i = 0; i < iter_divide_collapse; i++) {
    struct CollapseImpactMetrics {
      double length_ratio = 1E+100;
      double max_face_area_ratio = 1E+100;
      double volume_change_ratio = 1E+100;
      bool valid = false;
    };
    auto collapse_impact_metrics = [&](const networkLine* l) {
      CollapseImpactMetrics out;
      if (!l) return out;
      const auto faces = l->getBoundaryFaces();
      if (faces.size() != 2 || !faces[0] || !faces[1]) return out;
      const double local_mean_len = localEdgeLength(l);
      const double len = l->length();
      if (!(local_mean_len > 0.0) || !(len > 0.0) || !std::isfinite(local_mean_len) || !std::isfinite(len)) return out;
      const double area0 = boundaryFaceArea(faces[0]);
      const double area1 = boundaryFaceArea(faces[1]);
      const double mean_area0 = localMeanFaceArea(faces[0]);
      const double mean_area1 = localMeanFaceArea(faces[1]);
      if (!(area0 > 0.0) || !(area1 > 0.0) || !(mean_area0 > 0.0) || !(mean_area1 > 0.0)) return out;
      if (!std::isfinite(area0) || !std::isfinite(area1) || !std::isfinite(mean_area0) || !std::isfinite(mean_area1)) return out;
      out.length_ratio = len / local_mean_len;
      out.max_face_area_ratio = std::max(area0 / mean_area0, area1 / mean_area1);
      const double volume_before = signedPatchVolumeAroundLineAfterCollapse(l, false);
      const double volume_after = signedPatchVolumeAroundLineAfterCollapse(l, true);
      out.volume_change_ratio = std::abs(volume_after - volume_before) / std::pow(local_mean_len, 3);
      out.valid = std::isfinite(out.length_ratio) && std::isfinite(out.max_face_area_ratio) && std::isfinite(out.volume_change_ratio);
      return out;
    };
    constexpr double small_impact_length_ratio = 0.25;
    constexpr double small_impact_volume_change_ratio = 0.02;
    auto is_small_impact_collapse_candidate = [&](const networkLine* l) {
      const auto metrics = collapse_impact_metrics(l);
      if (!metrics.valid) return false;
      return metrics.length_ratio < small_impact_length_ratio && metrics.volume_change_ratio < small_impact_volume_change_ratio;
    };

    auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };
    auto collect_non_adjacent = [&](const std::vector<networkLine*>& candidates) {
      std::vector<networkLine*> batch;
      std::unordered_set<networkPoint*> used_points;
      batch.reserve(candidates.size());
      used_points.reserve(candidates.size() * 2);
      for (auto* l : candidates) {
        if (!l) continue;
        auto [p0, p1] = l->getPoints();
        if (!p0 || !p1) continue;
        if (used_points.contains(p0) || used_points.contains(p1)) continue;
        used_points.insert(p0);
        used_points.insert(p1);
        batch.emplace_back(l);
      }
      return batch;
    };
    auto gather_neighbor_lines = [&](const networkLine* l, std::unordered_set<networkLine*>& out) {
      if (!l) return;
      auto [p0, p1] = l->getPoints();
      auto add_from_point = [&](networkPoint* p) {
        if (!p) return;
        for (auto* nl : p->getBoundaryLines())
          if (nl) out.insert(nl);
      };
      add_from_point(p0);
      add_from_point(p1);
      for (auto* f : l->getBoundaryFaces())
        if (f)
          for (auto* nl : f->getLines())
            if (nl) out.insert(nl);
    };

    static const bool enable_corner_neumann_debug = (std::getenv("BEM_CORNER_DEBUG") && std::string(std::getenv("BEM_CORNER_DEBUG")) != "0");
    auto log_corner_connected_neumann_line = [&](const char* phase, const networkLine* l, const char* event = nullptr) {
      if (!enable_corner_neumann_debug) return;
      if (!isPostTypeCollapseTargetLine(l)) return;
      auto [p0, p1] = l->getPoints();
      const auto faces = l->getBoundaryFaces();
      const int n_faces = static_cast<int>(faces.size());
      const double len = l->length();
      const double local_mean_len_val = (n_faces == 2 && faces[0] && faces[1]) ? localEdgeLength(l) : -1.0;
      double alt0 = -1.0, alt1 = -1.0, alt_threshold = -1.0;
      double alt_ratio0 = -1.0, alt_ratio1 = -1.0;
      double aspect_ratio0 = -1.0, aspect_ratio1 = -1.0;
      double min_angle_deg0 = -1.0, min_angle_deg1 = -1.0;
      double area_ratio0 = -1.0, area_ratio1 = -1.0, normal_dot = -2.0;
      int common_points = -1;
      int p0_lines = p0 ? static_cast<int>(p0->getLines().size()) : -1;
      int p1_lines = p1 ? static_cast<int>(p1->getLines().size()) : -1;
      int opp0_lines = -1, opp1_lines = -1;
      bool face0_dirichlet = false, face0_neumann = false;
      bool face1_dirichlet = false, face1_neumann = false;
      if (n_faces == 2 && faces[0] && faces[1] && p0 && p1) {
        auto* f0 = faces[0]; auto* f1 = faces[1];
        face0_dirichlet = isDirichletBoundaryState(p0, f0) && isDirichletBoundaryState(p1, f0);
        face0_neumann = isNeumannBoundaryState(p0, f0) && isNeumannBoundaryState(p1, f0);
        face1_dirichlet = isDirichletBoundaryState(p0, f1) && isDirichletBoundaryState(p1, f1);
        face1_neumann = isNeumannBoundaryState(p0, f1) && isNeumannBoundaryState(p1, f1);
        auto [a, this0, b, l1_loc, p2, l2_loc] = f0->getPointsAndLines(const_cast<networkLine*>(l));
        auto [q0, this1, q1, e1, q2, e2] = f1->getPointsAndLines(const_cast<networkLine*>(l));
        if (this0 == l && this1 == l && a == q1 && b == q0) {
          const double area0_val = boundaryFaceArea(f0);
          const double area1_val = boundaryFaceArea(f1);
          const double max_edge0 = std::max({Norm(f0->getPoints()[0]->X - f0->getPoints()[1]->X), Norm(f0->getPoints()[1]->X - f0->getPoints()[2]->X), Norm(f0->getPoints()[2]->X - f0->getPoints()[0]->X)});
          const double max_edge1 = std::max({Norm(f1->getPoints()[0]->X - f1->getPoints()[1]->X), Norm(f1->getPoints()[1]->X - f1->getPoints()[2]->X), Norm(f1->getPoints()[2]->X - f1->getPoints()[0]->X)});
          alt0 = (max_edge0 > 1e-20) ? 2.0 * area0_val / max_edge0 : 0.0;
          alt1 = (max_edge1 > 1e-20) ? 2.0 * area1_val / max_edge1 : 0.0;
          alt_threshold = (local_mean_len_val > 0.0) ? 0.1 * local_mean_len_val : -1.0;
          alt_ratio0 = (local_mean_len_val > 0.0) ? alt0 / local_mean_len_val : -1.0;
          alt_ratio1 = (local_mean_len_val > 0.0) ? alt1 / local_mean_len_val : -1.0;
          aspect_ratio0 = faceAspectRatio(f0); aspect_ratio1 = faceAspectRatio(f1);
          min_angle_deg0 = minimumInteriorAngleDeg(f0); min_angle_deg1 = minimumInteriorAngleDeg(f1);
          const double mean_area0_val = localMeanFaceArea(f0);
          const double mean_area1_val = localMeanFaceArea(f1);
          area_ratio0 = (mean_area0_val > 0.0) ? area0_val / mean_area0_val : -1.0;
          area_ratio1 = (mean_area1_val > 0.0) ? area1_val / mean_area1_val : -1.0;
          normal_dot = Dot(f0->normal, f1->normal);
          common_points = static_cast<int>(Intersection(p0->getNeighborPointsOnSurfaces(), p1->getNeighborPointsOnSurfaces()).size());
          opp0_lines = p2 ? static_cast<int>(p2->getLines().size()) : -1;
          opp1_lines = q2 ? static_cast<int>(q2->getLines().size()) : -1;
        }
      }
      std::cout << Magenta << "[corner_neumann_debug] " << phase;
      if (event) std::cout << " event=" << event;
      std::cout << " p0=" << (p0 ? p0->index : -1) << " p1=" << (p1 ? p1->index : -1)
                << " len=" << len << " local_mean_len=" << local_mean_len_val << " faces=" << n_faces
                << " alt0=" << alt0 << " alt1=" << alt1 << " alt_threshold=" << alt_threshold
                << " alt_ratio0=" << alt_ratio0 << " alt_ratio1=" << alt_ratio1
                << " aspect_ratio0=" << aspect_ratio0 << " aspect_ratio1=" << aspect_ratio1
                << " min_angle_deg0=" << min_angle_deg0 << " min_angle_deg1=" << min_angle_deg1
                << " area_ratio0=" << area_ratio0 << " area_ratio1=" << area_ratio1
                << " normal_dot=" << normal_dot << " common_points=" << common_points
                << " line_flags={D:" << l->Dirichlet << ",N:" << l->Neumann << ",C:" << l->CORNER << "}"
                << " face0_flags={D:" << face0_dirichlet << ",N:" << face0_neumann << "}"
                << " face1_flags={D:" << face1_dirichlet << ",N:" << face1_neumann << "}"
                << " endpoint0_flags={D:" << (p0 ? p0->Dirichlet : false) << ",N:" << (p0 ? p0->Neumann : false) << ",C:" << (p0 ? p0->CORNER : false) << "}"
                << " endpoint1_flags={D:" << (p1 ? p1->Dirichlet : false) << ",N:" << (p1 ? p1->Neumann : false) << ",C:" << (p1 ? p1->CORNER : false) << "}"
                << " endpoint_lines={" << p0_lines << "," << p1_lines << "}"
                << " opposite_lines={" << opp0_lines << "," << opp1_lines << "}"
                << " x0=" << (p0 ? p0->X : Tddd{0., 0., 0.}) << " x1=" << (p1 ? p1->X : Tddd{0., 0., 0.})
                << colorReset << std::endl;
    };

    std::vector<networkLine*> corner_connected_neumann_shortlist;
    for (auto* l : water.getBoundaryLines())
      if (isPostTypeCollapseTargetLine(l)) corner_connected_neumann_shortlist.emplace_back(l);
    std::ranges::sort(corner_connected_neumann_shortlist, [](const auto* a, const auto* b) { return a->length() < b->length(); });
    if (corner_connected_neumann_shortlist.size() > 5) corner_connected_neumann_shortlist.resize(5);
    for (auto* l : corner_connected_neumann_shortlist) log_corner_connected_neumann_line("shortlist", l);

    auto should_split = [&](networkLine* l, double local_mean_len_val) {
      if (!l) return false;
      if (is_near_protected(l)) return false;
      if (local_mean_len_val <= 0.) return false;
      auto len = l->length();
      if (len <= 2.0 * limit_len || len <= 1.4 * local_mean_len_val) return false;
      auto surfaces = l->getBoundaryFaces();
      if (surfaces.size() != 2) return false;
      double min_neighbor_len = 1e20;
      for (const auto& f : surfaces)
        for (const auto& nl : f->getLines())
          if (nl != l && nl->length() > 0.) min_neighbor_len = std::min(min_neighbor_len, static_cast<double>(nl->length()));
      if (min_neighbor_len < 1e19 && len < 2.0 * min_neighbor_len) return false;
      return Dot(surfaces[0]->normal, surfaces[1]->normal) < cos_3rad;
    };

    if (surface_split) {
      int total_splits = 0;
      const int max_splits_per_step = std::clamp(static_cast<int>(water.getBoundaryLines().size()) / 50, 10, 100);
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines()) dirty.insert(l);
      for (auto iter = 0; iter < 20 && total_splits < max_splits_per_step; iter++) {
        std::vector<networkLine*> candidates;
        candidates.reserve(dirty.size());
        std::unordered_map<const networkLine*, double> local_mean_cache;
        local_mean_cache.reserve(dirty.size());
        for (auto* l : dirty) {
          if (!line_alive(l)) continue;
          auto it = local_mean_cache.find(l);
          auto lml = (it != local_mean_cache.end()) ? it->second : localEdgeLength(l);
          if (it == local_mean_cache.end()) local_mean_cache.emplace(l, lml);
          if (should_split(l, lml)) candidates.emplace_back(l);
        }
        if (candidates.empty()) break;
        auto batch = collect_non_adjacent(candidates);
        if (batch.empty()) break;
        bool divided_any = false;
        bool split_topo_error = false;
        std::unordered_set<networkLine*> touched;
        for (auto* l : batch) {
          if (!line_alive(l)) continue;
          auto len = l->length();
          auto lml = localEdgeLength(l);
          if (!should_split(l, lml)) continue;
          bool pre_ok = true;
          for (auto* f : l->getBoundaryFaces())
            for (auto* nl : f->getLines())
              if (nl && !nl->checkTopology()) { pre_ok = false; break; }
          if (!pre_ok) {
            std::cerr << Yellow << "[remesh] time_step " << time_step << ": skipping split (neighbor topology already bad)" << colorReset << std::endl;
            continue;
          }
          auto [sp0, sp1] = l->getPoints();
          bool near_collision = is_near_protected(l);
          std::cout << Red << "time_step " << time_step << ": splitting line"
                    << " len=" << len << " local_mean=" << lml
                    << " p0=" << sp0->X << " p1=" << sp1->X
                    << (near_collision ? " [NEAR_COLLISION_ZONE]" : "")
                    << colorReset << std::endl;
          l->Split();
          divided_any = true;
          ++total_splits;
          gather_neighbor_lines(l, touched);
          if (total_splits >= max_splits_per_step) break;
        }
        if (!divided_any) break;
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        for (const auto& l : water.getLines())
          if (!l->checkTopology()) { split_topo_error = true; break; }
        if (split_topo_error) {
          std::cerr << Red << "[remesh] time_step " << time_step << ": topology error after split batch, stopping further splits" << colorReset << std::endl;
          break;
        }
        dirty = std::move(touched);
        if (dirty.empty()) break;
      }
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    {
      bool post_split_ok = true;
      for (const auto& l : water.getLines())
        if (!l->checkTopology()) { post_split_ok = false; break; }
      if (!post_split_ok)
        throw step_failure("topology error after division at time_step " + std::to_string(time_step));
    }

    if (surface_flip && !heavy_collision_protection)
      flipIfBatched(water, {10 * rad, 10 * rad}, {3 * rad, 3 * rad}, "pre-collapse", &protected_halo_lines);

    if (surface_collapse)
      for (auto iter = 0; iter < 20; iter++) {
        bool found_small_line = false;
        if (surface_flip) {
          for (auto tiny_flip_iter = 0; tiny_flip_iter < 10; ++tiny_flip_iter) {
            std::vector<networkLine*> tiny_flip_candidates;
            std::unordered_set<networkLine*> seen;
            for (auto* face : water.getBoundaryFaces()) {
              if (!face) continue;
              const double area = boundaryFaceArea(face);
              if (!(area > 0.0) || !std::isfinite(area)) continue;
              const double mean_area = localMeanFaceArea(face);
              if (!(mean_area > 0.0) || !std::isfinite(mean_area)) continue;
              if (!(area / mean_area < min_local_face_area_ratio)) continue;
              for (auto* l : face->getLines()) {
                if (!l || !line_alive(l) || l->CORNER || is_near_protected(l)) continue;
                if (seen.insert(l).second) tiny_flip_candidates.emplace_back(l);
              }
            }
            if (tiny_flip_candidates.empty()) break;
            auto batch = collect_non_adjacent(RandomSample(tiny_flip_candidates));
            if (batch.empty()) break;
            bool changed_tiny_faces = false;
            std::unordered_set<networkLine*> touched;
            for (auto* l : batch) {
              if (!line_alive(l) || l->CORNER || is_near_protected(l)) continue;
              const double before_ratio = worstTinyFaceAreaRatioOnLine(l);
              const double predicted_ratio = predictedWorstTinyFaceAreaRatioAfterFlip(l);
              if (flipTinyFaceIfImproves(l, min_local_face_area_ratio)) {
                found_small_line = true;
                changed_tiny_faces = true;
                gather_neighbor_lines(l, touched);
                std::cout << "time_step " << time_step << ": line flipped due to tiny local-area face. ratio_before=" << before_ratio << " ratio_after_pred=" << predicted_ratio << std::endl;
              }
            }
            if (!changed_tiny_faces) break;
            water.setGeometricPropertiesForce();
            water.checkConnectivity();
          }
        }

        for (auto tiny_iter = 0; tiny_iter < 10; ++tiny_iter) {
          std::vector<networkLine*> tiny_candidates;
          std::unordered_set<networkLine*> seen;
          for (auto* l : water.getBoundaryLines()) {
            if (!l || !line_alive(l) || l->CORNER || is_near_protected(l)) continue;
            if (!is_small_impact_collapse_candidate(l)) continue;
            if (seen.insert(l).second) tiny_candidates.emplace_back(l);
          }
          for (auto* face : water.getBoundaryFaces()) {
            if (!face) continue;
            const double area = boundaryFaceArea(face);
            if (!(area > 0.0) || !std::isfinite(area)) continue;
            const double mean_area = localMeanFaceArea(face);
            if (!(mean_area > 0.0) || !std::isfinite(mean_area)) continue;
            if (!(area / mean_area < min_local_face_area_ratio)) continue;
            networkLine* shortest = nullptr;
            double shortest_len = 1E+100;
            for (auto* l : face->getLines()) {
              if (!l || !line_alive(l) || l->CORNER || is_near_protected(l)) continue;
              const double len = l->length();
              if (!(len > 0.0) || !std::isfinite(len)) continue;
              if (len < shortest_len) { shortest_len = len; shortest = l; }
            }
            if (shortest && seen.insert(shortest).second) tiny_candidates.emplace_back(shortest);
          }
          if (tiny_candidates.empty()) break;
          auto batch = collect_non_adjacent(tiny_candidates);
          if (batch.empty()) break;
          bool changed_tiny_faces = false;
          for (auto* l : batch) {
            if (!line_alive(l) || l->CORNER || is_near_protected(l)) continue;
            const auto attached_faces = l->getBoundaryFaces();
            double worst_ratio = 1E+100;
            for (auto* face : attached_faces) {
              const double area = boundaryFaceArea(face);
              const double mean_area = localMeanFaceArea(face);
              if (!(area > 0.0) || !(mean_area > 0.0)) continue;
              worst_ratio = std::min(worst_ratio, area / mean_area);
            }
            const auto metrics = collapse_impact_metrics(l);
            if (l->Collapse()) {
              found_small_line = true;
              changed_tiny_faces = true;
              std::cout << "time_step " << time_step << ": line merged due to tiny local-area face. worst_ratio=" << worst_ratio;
              if (metrics.valid)
                std::cout << " len_ratio=" << metrics.length_ratio << " max_face_ratio=" << metrics.max_face_area_ratio << " volume_change_ratio=" << metrics.volume_change_ratio;
              std::cout << std::endl;
            }
          }
          if (!changed_tiny_faces) break;
          water.setGeometricPropertiesForce();
          water.checkConnectivity();
        }

        std::unordered_set<networkLine*> dirty;
        for (auto* l : water.getBoundaryLines()) dirty.insert(l);
        while (!dirty.empty()) {
          std::vector<networkLine*> candidates;
          candidates.reserve(dirty.size());
          for (auto* l : dirty) { if (!line_alive(l)) continue; candidates.emplace_back(l); }
          if (candidates.empty()) break;
          auto batch = collect_non_adjacent(candidates);
          if (batch.empty()) break;
          bool changed_in_batch = false;
          std::unordered_set<networkLine*> touched;
          for (auto* l : batch) {
            if (!line_alive(l)) continue;
            if (l->CORNER) continue;
            if (is_near_protected(l)) continue;
            auto surfaces = l->getBoundaryFaces();
            if (surfaces.size() != 2) continue;
            auto f0 = surfaces[0]; auto f1 = surfaces[1];
            if (!f0 || !f1) continue;
            auto len = l->length();
            auto [a, b, c] = f0->getPoints(l);
            auto [_, __, d] = f1->getPoints(l);
            auto local_mean_len_val = localEdgeLength(l);
            {
              double area0 = TriangleArea(f0->getPoints()[0]->X, f0->getPoints()[1]->X, f0->getPoints()[2]->X);
              double area1 = TriangleArea(f1->getPoints()[0]->X, f1->getPoints()[1]->X, f1->getPoints()[2]->X);
              double max_edge0 = std::max({Norm(f0->getPoints()[0]->X - f0->getPoints()[1]->X), Norm(f0->getPoints()[1]->X - f0->getPoints()[2]->X), Norm(f0->getPoints()[2]->X - f0->getPoints()[0]->X)});
              double max_edge1 = std::max({Norm(f1->getPoints()[0]->X - f1->getPoints()[1]->X), Norm(f1->getPoints()[1]->X - f1->getPoints()[2]->X), Norm(f1->getPoints()[2]->X - f1->getPoints()[0]->X)});
              double alt0 = (max_edge0 > 1e-20) ? 2.0 * area0 / max_edge0 : 0.0;
              double alt1 = (max_edge1 > 1e-20) ? 2.0 * area1 / max_edge1 : 0.0;
              double alt_threshold = local_mean_len_val * 0.1;
              double mean_area0 = localMeanFaceArea(f0); double mean_area1 = localMeanFaceArea(f1);
              double area_ratio0 = (mean_area0 > 0.0) ? area0 / mean_area0 : 1E+100;
              double area_ratio1 = (mean_area1 > 0.0) ? area1 / mean_area1 : 1E+100;
              double alt_ratio0 = (local_mean_len_val > 0.0) ? alt0 / local_mean_len_val : 1E+100;
              double alt_ratio1 = (local_mean_len_val > 0.0) ? alt1 / local_mean_len_val : 1E+100;
              double min_angle_deg0 = minimumInteriorAngleDeg(f0);
              double min_angle_deg1 = minimumInteriorAngleDeg(f1);
              if (!l->CORNER && std::min(alt_ratio0, alt_ratio1) < 0.2 && std::min(area_ratio0, area_ratio1) < 0.2 && std::min(min_angle_deg0, min_angle_deg1) < 20.0) {
                gather_neighbor_lines(l, touched);
                if (std::ranges::find(corner_connected_neumann_shortlist, l) != corner_connected_neumann_shortlist.end())
                  log_corner_connected_neumann_line("collapse", l, "corner_neumann_area_alt_attempt");
                if (l->Collapse()) { changed_in_batch = true; continue; }
              }
              if (alt0 < alt_threshold || alt1 < alt_threshold) {
                gather_neighbor_lines(l, touched);
                if (std::ranges::find(corner_connected_neumann_shortlist, l) != corner_connected_neumann_shortlist.end())
                  log_corner_connected_neumann_line("collapse", l, "low_altitude_attempt");
                if (l->Collapse()) {
                  found_small_line = true; changed_in_batch = true;
                  std::cout << "time_step " << time_step << ": line merged due to low altitude. alt0=" << alt0 << " alt1=" << alt1 << " threshold=" << alt_threshold << std::endl;
                  continue;
                }
              }
            }
            if (l->Neumann && len > limit_len * 0.5) continue;
            if ((len < local_mean_len_val * 0.4 || len < limit_len)) {
              gather_neighbor_lines(l, touched);
              l->Collapse();
              found_small_line = true; changed_in_batch = true;
              std::cout << "time_step " << time_step << ": line merged due to small length. length = " << len << ", local mean length = " << local_mean_len_val << std::endl;
              continue;
            }
            double volume = TetrahedronVolume(a->X, b->X, c->X, d->X);
            double ref_volume = std::pow(local_mean_len_val, 3) * 0.5;
            (void)volume; (void)ref_volume;
            auto [p0, p1, p2] = f0->getPoints(l);
            auto ci0 = CircumradiusToInradius(p0->X, p1->X, p2->X);
            auto [q0, q1, q2] = f1->getPoints(l);
            auto ci1 = CircumradiusToInradius(q0->X, q1->X, q2->X);
            const double cos_5rad = std::cos(5.0 * rad);
            if (l->Neumann) { if (Dot(f0->normal, -f1->normal) < cos_5rad) continue; }
            if (l->length() < 0.1 * (Norm(p2->X - p1->X) + Norm(p2->X - p0->X) + Norm(q2->X - q1->X) + Norm(q2->X - q0->X)) / 4.)
              if (l->Collapse()) { found_small_line = true; changed_in_batch = true; std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl; continue; }
            if (ci0 > 50. || ci1 > 50.)
              if (l->Collapse()) { found_small_line = true; changed_in_batch = true; std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl; continue; }
            if (Norm(p2->X - q2->X) < std::sqrt(TriangleArea(p0->X, p1->X, p2->X) + TriangleArea(q0->X, q1->X, q2->X)) / 10.) {
              l->Flip(true); l->Collapse();
              found_small_line = true; changed_in_batch = true;
              std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl;
              continue;
            }
            if (Dot(f0->normal, -f1->normal) > cos_rad)
              if (l->Collapse()) { found_small_line = true; changed_in_batch = true; std::cout << "time_step " << time_step << ": line merged due to opposite normals. n0 = " << f0->normal << ", n1 = " << f1->normal << std::endl; continue; }
          }
          if (changed_in_batch) { water.setGeometricPropertiesForce(); water.checkConnectivity(); }
          if (!changed_in_batch) break;
          dirty = std::move(touched);
        }
        if (!found_small_line) break;
      }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    {
      bool post_merge_ok = true;
      for (const auto& l : water.getLines())
        if (!l->checkTopology()) { post_merge_ok = false; break; }
      if (!post_merge_ok)
        throw step_failure("topology error after merge at time_step " + std::to_string(time_step));
    }

    if (surface_flip && !heavy_collision_protection)
      flipIfBatched(water, {20 * rad, 20 * rad}, {5 * rad, 5 * rad}, "post-collapse", &protected_halo_lines);

    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    {
      bool post_flip_ok = true;
      for (const auto& l : water.getLines())
        if (!l->checkTopology()) { post_flip_ok = false; break; }
      if (!post_flip_ok)
        throw step_failure("topology error after flip at time_step " + std::to_string(time_step));
    }
  }

  {
    constexpr double max_fold_ratio = 0.15;
    auto folded = detectFoldedFaces(water, collision_settings.normal_reversal_cos, global_mean_len);
    size_t n_boundary_faces = water.getBoundaryFaces().size();
    double fold_ratio = (n_boundary_faces > 0) ? static_cast<double>(folded.size()) / n_boundary_faces : 0.0;
    if (!folded.empty())
      std::cout << Yellow << "[remesh] time_step " << time_step << ": post-remesh fold check: " << folded.size() << " / " << n_boundary_faces << " folded faces (ratio=" << fold_ratio << ")" << colorReset << std::endl;
    if (!skip_post_remesh_quality_rejects && fold_ratio > max_fold_ratio)
      throw step_failure("post-remesh fold ratio " + std::to_string(fold_ratio) + " > " + std::to_string(max_fold_ratio) + " at time_step " + std::to_string(time_step));
  }
  {
    const auto tiny = checkTinyFacesRelativeToLocalMean(water);
    if (tiny.worst_face && tiny.worst_area_ratio < min_local_face_area_ratio) {
      std::cout << Yellow << "[remesh] time_step " << time_step << ": tiny face check: face=" << tiny.worst_face->index << " area_ratio=" << tiny.worst_area_ratio << " area=" << tiny.worst_area << " local_mean_area=" << tiny.worst_local_mean_area << colorReset << std::endl;
      // throw する前に CORNER 辺を含めた collapse を試みる
      if (collapseFaceByIndexIfPossible(water, tiny.worst_face->index)) {
        std::cout << Green << "[remesh] time_step " << time_step << ": rescued tiny face " << tiny.worst_face->index << " by collapse (including CORNER edges)" << colorReset << std::endl;
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        // collapse 成功 → throw せずに続行（再チェックは次のループで行われる）
      } else if (!skip_post_remesh_quality_rejects) {
        throw step_failure("post-remesh tiny face area ratio " + std::to_string(tiny.worst_area_ratio) +
                           " < " + std::to_string(min_local_face_area_ratio) +
                           " on " + water.getName() +
                           " at time_step " + std::to_string(time_step) +
                           " (face_index=" + std::to_string(tiny.worst_face->index) +
                           ", subface_index=" + std::to_string(tiny.worst_subface_index) +
                           ", face_area=" + std::to_string(tiny.worst_area) +
                           ", local_mean_area=" + std::to_string(tiny.worst_local_mean_area) + ")");
      }
    }
  }

#ifdef USE_TETGEN
  if (tetrahedralize)
    water.tetrahedralize();
#endif
  water.setGeometricPropertiesForce();
  water.improveTetrahedraDelaunay();
}
