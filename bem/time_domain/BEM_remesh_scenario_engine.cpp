// BEM_remesh_scenario_engine.cpp
//
// Implementation of the per-edge multi-scenario trial engine.
// See BEM_remesh_scenario_engine.hpp for the public API and overview.
//
// Used by the Trial/local_patch remesh path for per-edge candidate patches.

#define BEM
#include "Network.hpp"
#include "basic_surface_geometry.hpp"
#include "BEM_remesh_scenario_engine.hpp"
#include "BEM_remesh_common.hpp"

#include <cmath>
#include <limits>
#include <unordered_set>
#include <utility>

namespace scenario_engine {

// ---------------------------------------------------------------------------
// Shared helpers (area_weighted_centroid lives in BEM_remesh_common.hpp)
// ---------------------------------------------------------------------------

namespace {

// Smallest dot product between adjacent face normal pairs in `net`.
// Returns 1.0 if no 2-face adjacency exists (flat / isolated edges).
inline double worst_normal_dot(const Network& net) {
  double min_dot = 1.0;
  for (auto* l : net.getLines()) {
    auto fs = l->getBoundaryFaces();
    if (fs.size() == 2 && fs[0] && fs[1]) {
      const double d = Dot(fs[0]->normal, fs[1]->normal);
      if (std::isfinite(d))
        min_dot = std::min(min_dot, d);
    }
  }
  return min_dot;
}

// Fold detection with source-aware relative threshold.
//
// Returns true if the post-op patch contains an adjacency worse than both:
//   - the absolute floor rs.quality_normal_flip_cos, AND
//   - src_min_dot - epsilon (pre-existing sharp features don't count as folds)
//
// Equivalent threshold = min(rs.quality_normal_flip_cos, src_min_dot - eps).
// Because `min` of cosines is the more permissive side, pre-existing sharp
// source patches get a looser threshold exactly matching their own sharpness.
inline bool patch_has_new_fold(const Network& post,
                               double src_min_dot,
                               double abs_floor_cos) {
  constexpr double eps = 0.05;
  const double thr = std::min(abs_floor_cos, src_min_dot - eps);
  for (auto* l : post.getLines()) {
    auto fs = l->getBoundaryFaces();
    if (fs.size() == 2 && fs[0] && fs[1] &&
        Dot(fs[0]->normal, fs[1]->normal) < thr)
      return true;
  }
  return false;
}

// patchQuality: higher is better. Returns kFlipSentinel on fold (relative).
// score = mean_ir + 1/(1+local_cv_mean)
//       - w * (fs_target_dist + 0.5 * theta_target_dist)
// rim (patch-border) points/edges/faces are excluded since they're pinned.
double patchQuality(const Network& net,
                    const Context& ctx,
                    double src_min_dot) {
  constexpr double kFlipSentinel = -1.0e18;
  if (patch_has_new_fold(net, src_min_dot, ctx.rs.quality_normal_flip_cos))
    return kFlipSentinel;

  const auto border_vec = net.getBorderPoints();
  std::unordered_set<const networkPoint*> rim(border_vec.begin(), border_vec.end());
  auto is_rim = [&](const networkPoint* p) { return p && rim.count(p); };

  double ir_sum = 0.0;
  double equilateral_dev_sum = 0.0;
  int ir_count = 0;
  int equilateral_dev_count = 0;
  for (auto* f : net.getBoundaryFaces()) {
    auto [p0, p1, p2] = f->getPoints();
    if (is_rim(p0) && is_rim(p1) && is_rim(p2))
      continue;
    ir_sum += InradiusToCircumradius(p0->X, p1->X, p2->X);
    const double equilateral_dev = equilateral_coordinate_deviation(p0->X, p1->X, p2->X);
    if (std::isfinite(equilateral_dev)) {
      equilateral_dev_sum += equilateral_dev;
      ++equilateral_dev_count;
    }
    ++ir_count;
  }
  const double mean_ir = (ir_count > 0) ? (ir_sum / static_cast<double>(ir_count)) : 1.0;
  const double equilateral_deviation = (equilateral_dev_count > 0)
                                           ? equilateral_dev_sum / static_cast<double>(equilateral_dev_count)
                                           : 0.0;

  double fs_dist_sum = 0.;
  int fs_count = 0;
  double th_dist_sum = 0.;
  int th_count = 0;
  for (auto* l : net.getBoundaryLines()) {
    auto [p0, p1] = l->getPoints();
    if (is_rim(p0) && is_rim(p1))
      continue;
    if ((l->Dirichlet || l->BCInterface) && ctx.free_surface_target_len > 1e-20) {
      const double Llen = l->length();
      if (Llen > 1e-20) {
        const double lr = std::log(Llen / ctx.free_surface_target_len);
        fs_dist_sum += lr * lr;
        ++fs_count;
      }
    }
    if (!l->geom_curvature.valid || !p0 || !p1)
      continue;
    const Tddd edge = p1->X - p0->X;
    const double theta = surface_geometry::edgeCurvatureAngle(
        edge, l->geom_curvature.k1, l->geom_curvature.k2,
        l->geom_curvature.PD1, l->geom_curvature.PD2);
    if (std::isfinite(theta) && ctx.theta_target > 1e-10 && theta > 1e-10) {
      const double lr = std::log(theta / ctx.theta_target);
      th_dist_sum += lr * lr;
      ++th_count;
    }
  }
  const double fs_target_dist = (fs_count > 0) ? fs_dist_sum / fs_count : 0.0;
  const double theta_target_dist = (th_count > 0) ? th_dist_sum / th_count : 0.0;

  double local_cv_mean = 0.;
  {
    double cv_sum = 0.;
    int cv_count = 0;
    for (auto* p : net.getBoundaryPoints()) {
      if (!p || is_rim(p))
        continue;
      cv_sum += p->localEdgeLengthCV();
      ++cv_count;
    }
    if (cv_count > 0)
      local_cv_mean = cv_sum / cv_count;
  }

  const double w = ctx.rs.quality_resolution_weight;
  return mean_ir + 1. / (1. + local_cv_mean) -
         w * (fs_target_dist + 0.5 * theta_target_dist + equilateral_deviation);
}

// Flip sub-op (valence 'v', Delaunay 'd', or 0 = AND both).
// Optional lower-case 's' suffix means "sharp-safe": skip sharp-edge flips.
void do_flip(Network& net, char variant, bool skip_sharp_edges, int valence_version, const Context& ctx) {
  const double feat_rad = ctx.feature_angle;
  const bool use_sharp_sector_valence = (variant == 'v' && valence_version == 2);
  const double cos_thr = ctx.rs.quality_normal_flip_cos;
  const double vertex_cos_thr =
      (cos_thr > std::cos(30.0 * M_PI / 180.0))
          ? cos_thr
          : std::cos(30.0 * M_PI / 180.0);
  const double min_ang = ctx.rs.quality_min_angle_deg * M_PI / 180.0;
  const auto border_pts_vec = net.getBorderPoints();
  std::unordered_set<const networkPoint*> border_set(
      border_pts_vec.begin(), border_pts_vec.end());

  auto is_outer_rim = [&](const networkLine* l) {
    auto [p0, p1] = l->getPoints();
    return p0 && p1 && border_set.count(p0) && border_set.count(p1);
  };
  auto eligible = [&](networkLine* pl, int& vg_out) -> bool {
    if (!pl || pl->getBoundaryFaces().size() != 2)
      return false;
    if (is_outer_rim(pl))
      return false;
    if (pl->BCInterface)
      return false;
    if (skip_sharp_edges && pl->SharpQ(feat_rad))
      return false;
    auto [p0, p1] = pl->getPoints();
    if (!p0 || !p1)
      return false;
    if ((int)p0->getLines().size() <= 3)
      return false;
    if ((int)p1->getLines().size() <= 3)
      return false;
    auto [vb, va] = use_sharp_sector_valence
                        ? sharpSectorValenceDeviationScore(pl, feat_rad)
                        : valenceDeviationScore(pl);
    if (use_sharp_sector_valence &&
        (vb == std::numeric_limits<int>::max() || va == std::numeric_limits<int>::max()))
      return false;
    vg_out = vb - va;
    const bool valence_gain = (vg_out > 0);
    const bool delaunay_violation = is_non_delaunay(pl);
    if (variant == 'v') {
      if (!valence_gain)
        return false;
    } else if (variant == 'd') {
      if (!delaunay_violation)
        return false;
    } else {
      if (!valence_gain)
        return false;
      if (!delaunay_violation)
        return false;
    }
    if (!flip_preserves_normals(pl, cos_thr, min_ang))
      return false;
    if (!flip_preserves_vertex_normals(pl, vertex_cos_thr))
      return false;
    if (!flip_dihedral_change_ok(pl))
      return false;
    if (!flip_midpoint_drift_ok(pl))
      return false;
    return true;
  };

  const int max_flip_iterations =
      std::max(8, 2 * static_cast<int>(net.Lines.size()));
  for (int flip_iter = 0; flip_iter < max_flip_iterations; ++flip_iter) {
    struct Cand {
      networkLine* l;
      int vg;
      double minang;
    };
    std::vector<Cand> cands;
    for (auto* pl : net.Lines) {
      int vg = 0;
      if (!eligible(pl, vg))
        continue;
      cands.push_back({pl, vg, flip_new_min_angle(pl)});
    }
    if (cands.empty())
      break;
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) {
                if (a.vg != b.vg)
                  return a.vg > b.vg;
                return a.minang > b.minang;
              });
    try {
      cands[0].l->Flip(true);
    } catch (...) {
    }
  }
}

// Smooth sub-op: selected kernel + feature-aware delta + step cap.
void do_smooth(Network& net,
               const Context& ctx,
               SimulationSettings::RemeshingSettings::SmoothMode smooth_mode,
               bool sharp_safe) {
  const double cos_thr = ctx.rs.quality_normal_flip_cos;
  const double min_ang = ctx.rs.quality_min_angle_deg * M_PI / 180.0;
  const SizingField sf = make_uniform_field(ctx.free_surface_target_len);
  auto touches_sharp_line = [&](const networkPoint* p) {
    if (!p)
      return false;
    for (auto* l : p->getBoundaryLines())
      if (l && l->SharpQ(ctx.feature_angle))
        return true;
    return false;
  };
  constexpr int sub_iterations = 3;
  for (int si = 0; si < sub_iterations; ++si) {
    for (auto* p : net.getPoints()) {
      if (p->BorderQ() || p->BCInterface)
        continue;
      if (sharp_safe && touches_sharp_line(p))
        continue;
      Tddd V;
      if (!remesh_smooth_delta(p, smooth_mode, sf, ctx.feature_angle, V))
        continue;
      const double disp = Norm(V);
      if (!(disp > 1e-14))
        continue;
      const double limit = remesh_smooth_step_limit_ratio(smooth_mode) * localEdgeLength(p);
      const Tddd V_capped = (disp <= limit) ? V : (limit / disp) * V;
      const double alpha = find_safe_smooth_step(p, V_capped, cos_thr, min_ang);
      if (alpha > 0.0)
        p->setXSingle(p->X + alpha * V_capped);
    }
    for (const auto& l : net.getLines())
      l->setBoundsSingle();
    for (const auto& p : net.getPoints())
      p->setFaces();
    for (const auto& f : net.getFaces())
      f->setGeometricProperties(ToX(f->getPoints()));
  }
}

// Apply one scenario string to a patch. Returns false on per-op failure.
bool run_patch_ops(const std::string& ops,
                   Network& net,
                   networkLine* cl,
                   const Context& ctx) {
  const double cos_thr = ctx.rs.quality_normal_flip_cos;
  const double min_ang = ctx.rs.quality_min_angle_deg * M_PI / 180.0;
  using RemeshSettings = SimulationSettings::RemeshingSettings;
  using ScenarioOp = RemeshSettings::ScenarioOp;
  using SmoothMode = RemeshSettings::SmoothMode;

  std::vector<RemeshSettings::ScenarioToken> tokens;
  std::string parse_error;
  if (!RemeshSettings::parse_scenario_tokens(
          ops, RemeshSettings::ScenarioScope::Patch, &tokens, &parse_error))
    return false;

  for (const auto& token : tokens) {
    if (token.op == ScenarioOp::Flip) {
      do_flip(net, token.variant, token.sharp_safe, token.valence_version, ctx);
    } else if (token.op == ScenarioOp::Collapse) {
      try {
        if (cl->getBoundaryFaces().size() != 2)
          return false;
        auto [cp0, cp1] = cl->getPoints();
        if (!cp0 || !cp1)
          return false;
        if (cp0->BorderQ() && cp1->BorderQ())
          return false;
        if (cp0->getBoundaryLines().size() <= 3)
          return false;
        if (cp1->getBoundaryLines().size() <= 3)
          return false;
        bool degen = false;
        for (auto* f : cl->getBoundaryFaces())
          if (f && (f->area < 1e-14 || !std::isfinite(f->area))) {
            degen = true;
            break;
          }
        if (degen)
          return false;
        const auto target_opt = feature_aware_collapse_target(cp0, cp1, ctx.feature_angle);
        if (!target_opt)
          return false;
        const Tddd target = *target_opt;
        if (!collapse_preserves_normals(cl, target, cos_thr, min_ang))
          return false;
        double kphi = 0.;
        double kphin = 0.;
        if (!ctx.rs.defer_scalar_interpolation) {
          kphi = 0.5 * (std::get<0>(cp0->phiphin) + std::get<0>(cp1->phiphin));
          kphin = 0.5 * (std::get<1>(cp0->phiphin) + std::get<1>(cp1->phiphin));
        }
        auto* kp = cl->Collapse(target);
        if (!kp)
          return false;
        if (!ctx.rs.defer_scalar_interpolation) {
          kp->phiphin = {kphi, kphin};
          kp->phiphin_t = {0., 0.};
        }
      } catch (...) {
        return false;
      }
    } else if (token.op == ScenarioOp::Smooth) {
      const SmoothMode smooth_mode =
          token.has_smooth_mode_suffix ? token.smooth_mode : ctx.rs.smooth_mode;
      do_smooth(net, ctx, smooth_mode, token.sharp_safe);
    } else if (token.op == ScenarioOp::Split) {
      auto [pA, pB] = cl->getPoints();
      if (!pA || !pB)
        return false;
      const Tddd mid = 0.5 * (pA->X + pB->X);
      auto* np = cl->Split(mid);
      if (!np)
        return false;
      if (!ctx.rs.defer_scalar_interpolation) {
        const double nphi = 0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
        const double nphin = 0.5 * (std::get<1>(pA->phiphin) + std::get<1>(pB->phiphin));
        np->phiphin = {nphi, nphin};
        np->phiphin_t = {0., 0.};
      }
    }
  }
  return true;
}

// Run one scenario on one patch, score it, check fold guard.
TrialResult run_trial(const std::string& ops,
                      std::shared_ptr<Network> patch,
                      networkLine* cl,
                      double src_min_dot,
                      const Context& ctx) {
  TrialResult r;
  r.ops = ops;
  if (!cl)
    return r;
  try {
    if (!run_patch_ops(ops, *patch, cl, ctx))
      return r;

    // rim preservation: all border points must still coincide with originals
    auto bset = patch->getBorderPoints();
    for (auto* bp : bset) {
      bool found = false;
      for (auto* orig_p : patch->copied_points) {
        if (Norm(bp->X - orig_p->X) < 1e-12) {
          found = true;
          break;
        }
      }
      if (!found)
        return r;
    }

    // 軽量幾何再計算 (setDodecaPoints は border-ring で OMP 不安全なため除外)
    for (const auto& l : patch->getLines())
      l->setBoundsSingle();
    for (const auto& p : patch->getPoints())
      p->setFaces();
    for (const auto& f : patch->getFaces())
      f->setGeometricProperties(ToX(f->getPoints()));
    patch->computePrincipalCurvatures();

    // Fold hard guard (relative to source).
    if (patch_has_new_fold(*patch, src_min_dot, ctx.rs.quality_normal_flip_cos))
      return r;

    r.score = patchQuality(*patch, ctx, src_min_dot);
    r.valid = true;
    r.patch = patch;
  } catch (...) {
  }
  return r;
}

} // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

TrialResult runTrials(Network& water,
                      networkLine* target_line,
                      const std::vector<std::string>& scenarios,
                      bool aggressive,
                      const Context& ctx) {
  TrialResult best;
  if (!target_line || scenarios.empty())
    return best;

  // Reference patch (for src_min_dot, score_before, debug source).
  auto patch0 = std::make_shared<Network>("file_name_is_not_given", "ref");
  const int ring = 2;

  water.copyLocalPatch(*patch0, target_line, ring);
  // IMPORTANT: copyLocalPatch does NOT always leave face normals in a state
  // usable for Dot comparisons. Recompute the same way run_trial does for
  // post-op patches — otherwise worst_normal_dot returns 1.0 on patches with
  // stale/zero normals and the relative fold threshold collapses back to the
  // absolute floor (rs.quality_normal_flip_cos). That's what causes split
  // trials on Neumann-sharp edges to be silently rejected by patch_has_new_fold.
  for (const auto& l : patch0->getLines())
    l->setBoundsSingle();
  for (const auto& p : patch0->getPoints())
    p->setFaces();
  for (const auto& f : patch0->getFaces())
    f->setGeometricProperties(ToX(f->getPoints()));
  patch0->computePrincipalCurvatures();
  const double src_min_dot = worst_normal_dot(*patch0);
  const double score_before = patchQuality(*patch0, ctx, src_min_dot);

  // copyLocalPatch は OMP-unsafe (race condition あり)。同期入力でも commit 数が
  // run ごとに変わる (例: 2929 → 2140 → 1024 splits at step=0) ことを実測で確認。
  // おそらく water 側の一部 mutable 状態 or new networkFace の link() が thread-safe でない。
  // → copyLocalPatch は serial、run_trial だけ OMP で回す。
  std::vector<std::shared_ptr<Network>> patches(scenarios.size());
  std::vector<networkLine*> copied_lines(scenarios.size(), nullptr);
  for (std::size_t i = 0; i < scenarios.size(); ++i) {
    patches[i] = std::make_shared<Network>("file_name_is_not_given", scenarios[i]);
    copied_lines[i] = water.copyLocalPatch(*patches[i], target_line, ring);
  }

  // run_trial 側は thread-local patch のみ触るので OMP 安全。
  // schedule(dynamic, 1): scenarios の重さが不均等 (PFvSFvSFvS は P より遥かに重い)。
  std::vector<TrialResult> results(scenarios.size());
  _Pragma("omp parallel for schedule(dynamic, 1)") for (std::size_t i = 0; i < scenarios.size(); ++i)
      results[i] = run_trial(scenarios[i], patches[i], copied_lines[i], src_min_dot, ctx);

  // Pick best valid trial.  Even when callers request an aggressive pass,
  // require improvement over the source patch; accepting a merely valid but
  // lower-scoring patch lets local remesh degrade nearby elements.
  (void)aggressive;
  double best_score = score_before + ctx.rs.quality_score_improve_margin;
  int n_valid = 0, n_improved = 0;
  for (auto& r : results) {
    if (!r.valid)
      continue;
    // Sentinel / NaN guard (defensive — run_trial's fold guard should have filtered).
    if (!std::isfinite(r.score) || r.score <= -1.0e17)
      continue;
    ++n_valid;
    if (r.score > best_score) {
      ++n_improved;
      best_score = r.score;
      best = std::move(r);
    }
  }
  if (n_valid == 0)
    best.reject_code = 1;
  else if (n_improved == 0)
    best.reject_code = 3;
  if (best.valid)
    best.reject_code = 0;
  best.source_patch = patch0;
  return best;
}

} // namespace scenario_engine
