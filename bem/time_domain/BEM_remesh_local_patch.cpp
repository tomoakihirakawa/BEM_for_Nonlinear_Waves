// BEM_remesh_local_patch.cpp
//
// Method: "local_patch_trial_multi_objective"
//
// User-built implementation of remeshing based on:
//   - Local 2-ring patch copy for each trigger edge
//   - Multi-scenario trial (P, PS, PF, PFS, PFSFS, PFSFSFS, ...) in parallel
//   - Selection by composite patchQuality score
//       (worst_ir - w·(area_deviation + equilateral_deviation))
//   - replacePatch for atomic application
//
// Geometry repair and scalar sampling are routed through
// BoundaryGeometryTargetProvider. Line nodes are remapped after topology edits
// are finished, instead of being treated as trial-time source of truth.
//
// Dispatch entry: BEM_remesh_main.cpp -> remesh_for_main_loop() routes here
// when rs.remesh_method == "local_patch_trial_multi_objective" (or "trial").

#define BEM
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <unordered_map>

#include "Network.hpp"
#include "basic_surface_geometry.hpp"
#include "BEM_remesh_main.hpp"
#include "BEM_remesh_common.hpp"
#include "BEM_remesh_debug_output.hpp"
#include "BEM_remesh_lightgbm.hpp"
#include "BEM_remesh_logging.hpp"
#include "BEM_remesh_trial_geometry_inputs.hpp"
#include "BEM_geometry_repair.hpp"
#include "BEM_pre_bvp_consistency.hpp"
#include "BEM_step_failure.hpp"
#include "BEM_legacy_globals.hpp"

// 前方宣言
bool collapseQualityGuardOK(const networkLine* l, double local_mean_len_val);
bool findBestCollapse(const networkLine* l, Tddd& bestX, double global_max_len);
bool shouldFlipByQuality(const networkLine* l);

void remesh_local_patch_trial_multi_objective(
    Network& water, int time_step,
    const SimulationSettings::RemeshingSettings& rs,
    bool skip_post_remesh_quality_rejects,
    const std::string& patch_output_directory,
    double simulation_time,
    PVDWriter* candidate_patches_pvd,
    PVDWriter* remeshed_patches_pvd,
    PVDWriter* trigger_edges_pvd,
    PVDWriter* edges_pvd,
    const BEMMeshPipeline::RemeshTrialGeometryInputs* trial_geometry_inputs,
    const std::string& phase_debug_tag,
    PVDWriter* remeshed_water_unrepaired_pvd,
    int step_retry,
    RemeshWaterSnapshotWriter remeshed_water_snapshot_writer) {

  // rs から個別フラグを展開
  const bool tetrahedralize = rs.tetrahedralize;
  const bool surface_flip = rs.surface_flip;
  const bool surface_split = rs.surface_split;
  const bool surface_collapse = rs.surface_collapse;
  const bool surface_smoothing = rs.surface_smoothing;
  const auto& collision_settings = rs.collision;
  const bool remesh_debug_output = rs.debug_output;
  const bool remesh_debug_verbose = rs.debug_verbose;
  const bool remesh_debug_edges = rs.debug_edges || rs.debug_output;
  const bool stratified_training_sweep =
      rs.training_data_enabled &&
      rs.training_data_mode == "stratified_op_scenario_sweep";
  const bool use_lightgbm_ranking =
      rs.lightgbm_ranking_enabled && !stratified_training_sweep;
  const bool use_lightgbm_recommender =
      use_lightgbm_ranking && rs.lightgbm_mode == "recommender";

  BEMMeshPipeline::RemeshAI::LightGBMTextModel lightgbm_scenario_model;
  if (use_lightgbm_ranking) {
    if (rs.lightgbm_mode != "ranker" && rs.lightgbm_mode != "recommender")
      throw step_failure("remesh_lightgbm mode must be ranker or recommender: " +
                         rs.lightgbm_mode);
    if (rs.lightgbm_safe_model_path.empty())
      throw step_failure("remesh_lightgbm enabled but remesh_lightgbm_safe_model_path is empty");
    if (!lightgbm_scenario_model.load(rs.lightgbm_safe_model_path))
      throw step_failure("failed to load remesh LightGBM model: " +
                         lightgbm_scenario_model.lastError());
    std::cout << "[remesh:lgbm] scenario_ranking=enabled"
              << " mode=" << rs.lightgbm_mode
              << " model=" << rs.lightgbm_safe_model_path
              << " top_k=" << (use_lightgbm_recommender ? rs.lightgbm_recommender_top_k : 0)
              << std::endl;
  }

  auto sanitize_phase_tag = [](std::string s) {
    for (char& c : s) {
      const bool ok = (c >= '0' && c <= '9') ||
                      (c >= 'A' && c <= 'Z') ||
                      (c >= 'a' && c <= 'z') ||
                      c == '_' || c == '-';
      if (!ok)
        c = '_';
    }
    return s.empty() ? std::string("phase2") : s;
  };
  auto dump_remesh_phase_mesh = [&](const std::string& stage) {
    if (!remesh_debug_output || patch_output_directory.empty())
      return;
    const std::filesystem::path out_dir(patch_output_directory);
    std::filesystem::create_directories(out_dir);
    const std::string tag = sanitize_phase_tag(phase_debug_tag);
    const std::string filename = water.getName() + "_" + stage + "_" + tag + "_" + std::to_string(time_step) + ".vtu";
    mk_vtu((out_dir / filename).string(), water.getBoundaryFaces());
    if (remesh_debug_verbose)
      std::cout << "[mesh:dump] wrote " << (out_dir / filename).string() << std::endl;
  };
  struct BCInterfaceContactStats {
    int point_total = 0;
    int point_no_contact = 0;
    int line_total = 0;
    int line_no_contact = 0;
    int total() const { return point_total + line_total; }
    int no_contact() const { return point_no_contact + line_no_contact; }
  };

  auto count_bcinterface_no_contact = [](const Network& net, bool ignore_patch_border_lines = false) {
    BCInterfaceContactStats stats;
    auto point_touches_patch_border = [](const networkPoint* p) {
      if (!p)
        return false;
      for (auto* l : p->getBoundaryLines())
        if (l && l->BorderQ())
          return true;
      return false;
    };
    for (auto* p : net.getBoundaryPoints()) {
      if (!p || !p->BCInterface)
        continue;
      if (ignore_patch_border_lines && point_touches_patch_border(p))
        continue;
      ++stats.point_total;
      if (p->getContactFaces().empty())
        ++stats.point_no_contact;
    }
    for (auto* l : net.getBoundaryLines()) {
      if (!l || !l->BCInterface)
        continue;
      if (ignore_patch_border_lines && l->BorderQ())
        continue;
      ++stats.line_total;
      if (l->getContactFaces().empty())
        ++stats.line_no_contact;
    }
    return stats;
  };

  auto refresh_contact_state_after_topology_change = [&](const char* reason,
                                                         bool hard_fail_on_no_contact = true) -> BCInterfaceContactStats {
    if (!trial_geometry_inputs || !trial_geometry_inputs->enabled() || !trial_geometry_inputs->body_objects)
      return {};
    water.setGeometricPropertiesForce();
    water.makeBuckets(water.getScale() / 10.);
    for (auto* body : *trial_geometry_inputs->body_objects) {
      if (!body)
        continue;
      body->setGeometricPropertiesForce();
      body->makeBuckets(body->getScale() / 10.);
    }
    refreshBoundaryStatesAndTypes(&water,
                                  *trial_geometry_inputs->body_objects,
                                  BoundaryStateRefreshOptions{.verbose = remesh_debug_verbose});
    computeAllBCInterfaceMidpointOffsets(&water);

    const auto stats = count_bcinterface_no_contact(water, /*ignore_patch_border_lines=*/false);
    if (remesh_debug_verbose) {
      std::cout << "[remesh:contact] refresh reason=" << reason
                << " bci_no_contact=" << stats.no_contact() << "/" << stats.total()
                << " point=" << stats.point_no_contact << "/" << stats.point_total
                << " line=" << stats.line_no_contact << "/" << stats.line_total
                << std::endl;
    }
    if (hard_fail_on_no_contact && stats.no_contact() > 0) {
      throw step_failure(std::string("remesh produced BCInterface point/line(s) without raw contact after ") +
                         reason + ": " + std::to_string(stats.no_contact()) + "/" +
                         std::to_string(stats.total()) +
                         " (points=" + std::to_string(stats.point_no_contact) + "/" +
                         std::to_string(stats.point_total) +
                         ", lines=" + std::to_string(stats.line_no_contact) + "/" +
                         std::to_string(stats.line_total) + ")");
    }
    return stats;
  };

  auto refresh_trial_patch_boundary_state = [&](Network& net) {
    if (!trial_geometry_inputs || !trial_geometry_inputs->enabled() || !trial_geometry_inputs->body_objects)
      return BCInterfaceContactStats{};

    refreshBoundaryStatesAndTypes(&net,
                                  *trial_geometry_inputs->body_objects,
                                  BoundaryStateRefreshOptions{.verbose = false});
    return count_bcinterface_no_contact(net, /*ignore_patch_border_lines=*/true);
  };

  // Remesh debug VTU output is intentionally kept behind remesh_debug_output.
  // The simulation path still records profiles and rejection counts, while
  // large candidate/patch/trigger VTU files are opt-in diagnostics.
  int attempt_counter = 0;

  enum RejectCode : int {
    kRejectSuccess = 0,
    kRejectNoValid = 1,
    kRejectHdExceeded = 2,
    kRejectScoreNoImprove = 3,
    kRejectReplaceFailed = 4,
    kRejectRepairFailed = 5,
    kRejectPatchOpsFailed = 6,
    kRejectBorderMoved = 7,
    kRejectSubsurfaceAltitude = 8,
    kRejectBCInterfaceNoContact = 9,
    kRejectTinyFace = 10,
    kRejectWaterlineQuality = 11,
    kRejectIrRegression = 12,
    kRejectException = 13,
    kRejectCodeCount = 14
  };

  auto rejectCodeName = [](int code) -> const char* {
    switch (code) {
    case kRejectSuccess:
      return "success";
    case kRejectNoValid:
      return "novalid";
    case kRejectHdExceeded:
      return "hd";
    case kRejectScoreNoImprove:
      return "score";
    case kRejectReplaceFailed:
      return "replace";
    case kRejectRepairFailed:
      return "repair";
    case kRejectPatchOpsFailed:
      return "ops";
    case kRejectBorderMoved:
      return "border";
    case kRejectSubsurfaceAltitude:
      return "subsurface";
    case kRejectBCInterfaceNoContact:
      return "bci_contact";
    case kRejectTinyFace:
      return "tiny_face";
    case kRejectWaterlineQuality:
      return "waterline_quality";
    case kRejectIrRegression:
      return "ir_regression";
    case kRejectException:
      return "exception";
    default:
      return "unknown";
    }
  };

  struct RemeshEdgeTrialStats {
    std::array<int, 4> attempts{};
    std::array<int, 4> successes{};
    std::array<std::array<int, kRejectCodeCount>, 4> rejects{};
    int since_report = 0;

    void reset() {
      attempts.fill(0);
      successes.fill(0);
      for (auto& row : rejects)
        row.fill(0);
      since_report = 0;
    }
  };

  auto edgeTrialCategory = [](const networkLine* l) -> int {
    if (l && l->BCInterface)
      return 2;
    if (l && l->Dirichlet)
      return 1;
    if (l && l->Neumann)
      return 0;
    return 3;
  };

  auto printEdgeTrialStats = [&](const char* op, RemeshEdgeTrialStats& stats,
                                 int scenario_count, bool partial = false,
                                 bool final = false) {
    if (stats.since_report <= 0)
      return;
    constexpr std::array<const char*, 4> labels = {"N", "D", "B", "O"};
    std::cout << "[remesh:edge_trials] step=" << time_step
              << " op=" << op
              << " scenarios=" << scenario_count
              << " window=" << stats.since_report
              << (partial ? " partial=1" : "")
              << (final ? " final=1" : "");
    for (int c = 0; c < 4; ++c) {
      if (c == 3 && stats.attempts[c] == 0)
        continue;
      std::cout << " " << labels[c] << ":" << stats.successes[c] << "/" << stats.attempts[c] << "(";
      bool first = true;
      for (int rc = 1; rc < kRejectCodeCount; ++rc) {
        if (stats.rejects[c][rc] == 0)
          continue;
        if (!first)
          std::cout << ",";
        std::cout << rejectCodeName(rc) << ":" << stats.rejects[c][rc];
        first = false;
      }
      if (stats.attempts[c] > 0 && stats.rejects[c][kRejectWaterlineQuality] == 0) {
        if (!first)
          std::cout << ",";
        std::cout << rejectCodeName(kRejectWaterlineQuality) << ":0";
      }
      std::cout << ")";
    }
    std::cout << std::endl;
    if (!final)
      stats.reset();
  };

  auto recordEdgeTrialStatsByCategory = [&](RemeshEdgeTrialStats& stats, int category,
                                            bool success, int reject_code) {
    const int c = (category >= 0 && category < 4) ? category : 3;
    ++stats.attempts[c];
    ++stats.since_report;
    if (success) {
      ++stats.successes[c];
      return;
    }
    if (reject_code > kRejectSuccess && reject_code < kRejectCodeCount)
      ++stats.rejects[c][reject_code];
    else
      ++stats.rejects[c][kRejectNoValid];
  };

  auto recordEdgeTrialStats = [&](RemeshEdgeTrialStats& stats, const networkLine* l,
                                  bool success, int reject_code) {
    recordEdgeTrialStatsByCategory(stats, edgeTrialCategory(l), success, reject_code);
  };

  auto recordBothEdgeTrialStats = [&](RemeshEdgeTrialStats& window_stats,
                                      RemeshEdgeTrialStats& total_stats,
                                      int category, bool success, int reject_code) {
    recordEdgeTrialStatsByCategory(window_stats, category, success, reject_code);
    recordEdgeTrialStatsByCategory(total_stats, category, success, reject_code);
  };

  BEMMeshPipeline::RemeshDebug::Recorder remesh_debug(
      remesh_debug_output,
      water,
      time_step,
      patch_output_directory,
      simulation_time,
      candidate_patches_pvd,
      remeshed_patches_pvd,
      trigger_edges_pvd);

  auto lineFlagIndex = [](const networkLine* l) -> int {
    return BEMMeshPipeline::RemeshDebug::lineFlagIndex(l);
  };
  auto edgeFaceSideIndex = [](const networkLine* l, const networkFace* f) -> int {
    const auto bc = getEdgeNodeFaceBoundaryType(l, f);
    if (bc == NodeFaceBoundaryType::Neumann)
      return 0;
    if (bc == NodeFaceBoundaryType::Dirichlet)
      return 1;
    return 2;
  };
  std::array<std::array<long long, 3>, 3> rt_subsurface_worst_side_by_trigger{};
  std::array<std::array<long long, 4>, 3> rt_subsurface_worst_line_flag_by_trigger{};
  std::array<long long, 3> rt_subsurface_worst_samples{};

  auto collectCandidatePatch = [&](const Network& patch, int op, int reason, int aid,
                                   int success, int reject_code) {
    remesh_debug.collectCandidatePatch(patch, op, reason, aid, success, reject_code);
  };
  auto collectRemeshedPatch = [&](const Network& patch, int op, int reason, int aid) {
    remesh_debug.collectRemeshedPatch(patch, op, reason, aid);
  };
  auto snapshotPatchTris = [&](const Network& patch) {
    return remesh_debug.snapshotPatchTris(patch);
  };
  auto commitRemeshedSnapshot = [&](const std::vector<T3Tddd>& tris, int op, int reason, int aid) {
    remesh_debug.commitRemeshedSnapshot(tris, op, reason, aid);
  };
  auto collectTriggerEdge = [&](const networkLine* l, int op, int reason, int aid,
                                int success, int reject_code) {
    remesh_debug.collectTriggerEdge(l, op, reason, aid, success, reject_code);
  };
  auto markAttemptSuccess = [&](int aid) {
    remesh_debug.markAttemptSuccess(aid, kRejectSuccess);
  };
  auto markAttemptReplaceFailed = [&](int aid) {
    remesh_debug.markAttemptReplaceFailed(aid, kRejectReplaceFailed);
  };
  auto snapshotFaceTri = [](const networkFace* f) {
    return BEMMeshPipeline::RemeshDebug::snapshotFaceTri(f);
  };

  // ========================================================================
  // [1] パラメータ設定
  // ========================================================================

  const double rad = M_PI / 180.0;
  const double global_mean_len = Mean(extLength(water.getLines()));
  // 自由表面の目標辺長は水平方向（x,y）のメッシュスケールから決める
  double free_surface_target_len;
  double horiz_diag_dbg = 0.0;
  {
    const auto bnds = water.getBounds();
    double dx = std::get<1>(std::get<0>(bnds)) - std::get<0>(std::get<0>(bnds));
    double dy = std::get<1>(std::get<1>(bnds)) - std::get<0>(std::get<1>(bnds));
    horiz_diag_dbg = std::sqrt(dx * dx + dy * dy);
    free_surface_target_len = horiz_diag_dbg / rs.len_target_divisor;
  }
  const double global_min_len = free_surface_target_len * rs.len_global_min_ratio;
  const double global_max_len = free_surface_target_len * rs.len_global_max_ratio;

  // Free surface / wall 辺長の内訳を診断ログに出す
  {
    std::size_t n_dirichlet = 0, n_corner = 0, n_neumann = 0;
    double fs_len_min = std::numeric_limits<double>::infinity();
    double fs_len_max = 0.0, fs_len_sum = 0.0;
    double wall_len_min = std::numeric_limits<double>::infinity();
    double wall_len_max = 0.0, wall_len_sum = 0.0;
    for (auto* l : water.getBoundaryLines()) {
      double L = l->length();
      if (l->BCInterface) {
        ++n_corner;
      } else if (l->Dirichlet) {
        ++n_dirichlet;
        fs_len_min = std::min(fs_len_min, L);
        fs_len_max = std::max(fs_len_max, L);
        fs_len_sum += L;
      } else {
        ++n_neumann;
        wall_len_min = std::min(wall_len_min, L);
        wall_len_max = std::max(wall_len_max, L);
        wall_len_sum += L;
      }
    }
    // 現在の needsSplit / needsCollapse は len > 4/3*target / len < 4/5*target
    // に統一。ログでは実際に使われている 4/3, 4/5 を表示する。
    const double fs_thr_split = (4.0 / 3.0) * free_surface_target_len;
    const double fs_thr_collapse = (4.0 / 5.0) * free_surface_target_len;
    if (remesh_debug_verbose) {
      std::cout << "[remesh:target] horiz_diag=" << horiz_diag_dbg
                << " divisor=" << rs.len_target_divisor
                << " fs_target=" << free_surface_target_len
                << " fs_split_thr=" << fs_thr_split << " (=4/3·fs_target)"
                << " fs_collapse_thr=" << fs_thr_collapse << " (=4/5·fs_target)"
                << " global_max=" << global_max_len << " (legacy, unused)"
                << " global_min=" << global_min_len << " (legacy, unused)"
                << std::endl;
      std::cout << "[remesh:edges] Dirichlet=" << n_dirichlet
                << " (len min=" << fs_len_min
                << " max=" << fs_len_max
                << " mean=" << (n_dirichlet ? fs_len_sum / n_dirichlet : 0.0) << ")"
                << "  Neumann=" << n_neumann
                << " (len min=" << wall_len_min
                << " max=" << wall_len_max
                << " mean=" << (n_neumann ? wall_len_sum / n_neumann : 0.0) << ")"
                << "  BCInterface=" << n_corner << std::endl;
    }
  }
  constexpr double min_local_face_area_ratio = 0.05;

  // 曲面忠実度（θ）: theta_target = 2π/N で円筒 N 分割が目標
  const bool curvature_remesh_enabled = rs.theta_enabled;
  const double theta_target = 2.0 * M_PI / rs.theta_target_N;
  const double theta_split = rs.theta_split_ratio * theta_target;
  const double theta_collapse = rs.theta_collapse_ratio * theta_target;

  // メッシュ均一性（grading）
  constexpr double a_grading = 0.5;
  constexpr double grading_split = 1.0 + a_grading; // g > 1.5 → split
  constexpr double grading_collapse = 1.0 - a_grading;

  // ========================================================================
  // [2] 前処理: 曲率計算、X_mid 初期化、四面体整理
  // ========================================================================

  if (curvature_remesh_enabled)
    water.computePrincipalCurvatures();

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

  // ========================================================================
  // [3] 衝突検出・解決 → 保護領域の設定
  // ========================================================================

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
    std::cout << Yellow << "[remesh:quality] step=" << time_step
              << " collision_unresolved=1 fold_ratio=" << fold_ratio
              << " protected_ratio=" << protected_ratio
              << colorReset << std::endl;
  }

  const auto protected_halo_lines = remesh_detail::build_protection_halo(protected_lines);

  water.setGeometricPropertiesForce();
  water.checkConnectivity();

  const std::size_t boundary_line_count = water.getBoundaryLines().size();
  if (boundary_line_count > 0 && protected_lines.size() > 0) {
    std::cout << Yellow << "[remesh:protect] step=" << time_step
              << " protected_lines=" << protected_lines.size() << "/" << boundary_line_count
              << " halo_lines=" << protected_halo_lines.size()
              << colorReset << std::endl;
  }
  const bool heavy_collision_protection = false;

  // ========================================================================
  // [4] topology 事前チェック
  // ========================================================================

  {
    int n_pre_topo_err = 0;
    for (auto* l : water.getBoundaryLines())
      if (!l->checkTopology())
        n_pre_topo_err++;
    if (n_pre_topo_err > 0) {
      std::cerr << Red << "[remesh:topology] step=" << time_step
                << " pre_existing_errors=" << n_pre_topo_err
                << colorReset << std::endl;
      throw step_failure("pre-existing topology errors (" + std::to_string(n_pre_topo_err) + " lines) at time_step " + std::to_string(time_step));
    }
  }

  // ========================================================================
  // [5] split / flip / collapse ループ（iter_split_collapse 回繰り返し）
  // ========================================================================

  const bool training_data_record_enabled =
      rs.training_data_enabled && rs.training_data_mode != "off";
  const int iter_split_collapse = training_data_record_enabled
                                      ? std::max(0, rs.training_data_max_iter)
                                      : std::max(0, rs.iter_split_collapse);

  auto joinScenarioList = [](const std::vector<std::string>& scenarios) {
    std::ostringstream os;
    os << "{";
    for (std::size_t i = 0; i < scenarios.size(); ++i) {
      if (i)
        os << ",";
      os << scenarios[i];
    }
    os << "}";
    return os.str();
  };
  auto scenarioListHas = [](const std::vector<std::string>& scenarios,
                            const std::string& needle) {
    return std::any_of(scenarios.begin(), scenarios.end(),
                       [&](const std::string& s) {
                         return s.find(needle) != std::string::npos;
                       });
  };
  const bool active_scenarios_have_fv2 =
      scenarioListHas(rs.split_scenarios, "Fv2") ||
      scenarioListHas(rs.collapse_scenarios, "Fv2") ||
      scenarioListHas(rs.flip_scenarios, "Fv2") ||
      scenarioListHas(rs.global_scenarios, "Fv2");
  std::cout << "[remesh:scenarios] step=" << time_step
            << " has_Fv2=" << (active_scenarios_have_fv2 ? 1 : 0)
            << " split=" << rs.split_scenarios.size() << joinScenarioList(rs.split_scenarios)
            << " collapse=" << rs.collapse_scenarios.size() << joinScenarioList(rs.collapse_scenarios)
            << " flip=" << rs.flip_scenarios.size() << joinScenarioList(rs.flip_scenarios)
            << " global=" << rs.global_scenarios.size() << joinScenarioList(rs.global_scenarios)
            << std::endl;

  const std::filesystem::path training_out_dir = patch_output_directory.empty()
                                                     ? std::filesystem::path(".")
                                                     : std::filesystem::path(patch_output_directory);
  const std::filesystem::path training_trial_csv =
      training_out_dir / (water.getName() + "_remesh_training_trials.csv");
  const std::filesystem::path training_iter_csv =
      training_out_dir / (water.getName() + "_remesh_training_iters.csv");
  const std::filesystem::path training_applied_csv =
      training_out_dir / (water.getName() + "_remesh_training_applied.csv");
  const std::filesystem::path training_metadata_json =
      training_out_dir / (water.getName() + "_remesh_training_metadata.json");

  auto csvString = [](const std::string& s) {
    std::string out = "\"";
    for (char c : s) {
      if (c == '"')
        out += "\"\"";
      else
        out += c;
    }
    out += "\"";
    return out;
  };
  auto appendWithHeader = [](const std::filesystem::path& path,
                             const std::string& header,
                             const std::string& row) {
    const bool need_header = !std::filesystem::exists(path) ||
                             std::filesystem::file_size(path) == 0;
    std::ofstream ofs(path, std::ios::app);
    if (need_header)
      ofs << header << '\n';
    ofs << row << '\n';
  };

  constexpr const char* training_trial_header =
      "time_step,iter,op,reason,edge_category,p0_index,p1_index,len,len_ratio,theta,theta_ratio,"
      "scenario,valid,reject_code,reject_name,score_before,score_after,score_gain,"
      "ir_mean_before,ir_mean_after,ir_mean_gain,ir_worst_before,ir_worst_after,"
      "area_deviation_before,area_deviation_after,area_deviation_gain,"
      "equilateral_deviation_before,equilateral_deviation_after,equilateral_deviation_gain,"
      "hd,hd_ref,hd_body,"
      "hd_limit,bcinterface_no_contact,subsurface_altitude_rel,saved_kind";
  constexpr const char* training_iter_header =
      "time_step,iter,candidate_edges,trial_count,safe_count,selected_count,applied_count,failed_apply_count";
  constexpr const char* training_applied_header =
      "time_step,iter,op,reason,edge_category,p0_index,p1_index,scenario,score_gain,weight,applied,reject_code,reject_name";

  if (training_data_record_enabled) {
    std::filesystem::create_directories(training_out_dir);
    auto writeScenarioArray = [&](std::ofstream& ofs,
                                  const char* name,
                                  const std::vector<std::string>& scenarios,
                                  bool trailing_comma) {
      ofs << "  \"" << name << "\": [";
      for (std::size_t i = 0; i < scenarios.size(); ++i) {
        if (i)
          ofs << ", ";
        ofs << csvString(scenarios[i]);
      }
      ofs << "]" << (trailing_comma ? "," : "") << "\n";
    };
    std::ofstream meta(training_metadata_json);
    meta << "{\n";
    meta << "  \"mode\": " << csvString(rs.training_data_mode) << ",\n";
    meta << "  \"training_data_mode\": " << csvString(rs.training_data_mode) << ",\n";
    meta << "  \"time_step\": " << time_step << ",\n";
    meta << "  \"max_iter\": " << iter_split_collapse << ",\n";
    meta << "  \"max_edges_per_iter\": " << rs.training_data_max_edges_per_iter << ",\n";
    meta << "  \"random_extra_edges_per_iter\": " << rs.training_data_random_extra_edges_per_iter << ",\n";
    meta << "  \"max_apply_per_iter\": " << rs.training_data_max_apply_per_iter << ",\n";
    meta << "  \"sweep_split_edges_per_iter\": " << rs.training_data_sweep_split_edges_per_iter << ",\n";
    meta << "  \"sweep_collapse_edges_per_iter\": " << rs.training_data_sweep_collapse_edges_per_iter << ",\n";
    meta << "  \"sweep_flip_edges_per_iter\": " << rs.training_data_sweep_flip_edges_per_iter << ",\n";
    meta << "  \"sweep_bcinterface_fraction\": " << rs.training_data_sweep_bcinterface_fraction << ",\n";
    meta << "  \"sweep_dirichlet_fraction\": " << rs.training_data_sweep_dirichlet_fraction << ",\n";
    meta << "  \"sweep_neumann_fraction\": " << rs.training_data_sweep_neumann_fraction << ",\n";
    meta << "  \"min_score_gain\": " << rs.training_data_min_score_gain << ",\n";
    meta << "  \"quality_hd_diag_ratio\": " << rs.quality_hd_diag_ratio << ",\n";
    meta << "  \"quality_hd_curv_ratio\": " << rs.quality_hd_curv_ratio << ",\n";
    meta << "  \"uses_target_aware_hd\": true,\n";
    meta << "  \"region_order\": [\"BCInterface\", \"Dirichlet\", \"Neumann\"],\n";
    writeScenarioArray(meta, "split_scenarios", rs.split_scenarios, true);
    writeScenarioArray(meta, "collapse_scenarios", rs.collapse_scenarios, true);
    writeScenarioArray(meta, "flip_scenarios", rs.flip_scenarios, false);
    meta << "}\n";
  }

  // [profile] 各セクションの累積時間
  double t_split = 0., t_collapse = 0., t_smooth_ir = 0., t_smooth_angle = 0., t_topo = 0., t_curv = 0.;
  // runTrials の内訳
  double t_rt_ref = 0., t_rt_copy = 0., t_rt_run = 0., t_rt_pick = 0.;
  int rt_calls = 0;
  long long rt_repair_aware_calls = 0;
  long long rt_repair_failed = 0;
  long long rt_repaired_points = 0;
  long long rt_repaired_lines = 0;
  long long rt_hd_checks = 0;
  long long rt_hd_early_exits = 0;
  long long rt_hd_samples = 0;
  long long rt_hd_nearest_calls = 0;
  long long rt_hd_target_missing = 0;
  double rt_hd_ref_max = 0.0;
  double rt_hd_body_max = 0.0;
  long long rt_altitude_rescue_attempts = 0;
  long long rt_altitude_rescue_successes = 0;
  long long rt_altitude_rescue_split_worst_line = 0;
  long long rt_altitude_rescue_split_longest_face_edge = 0;
  long long rt_subsurface_trigger_edges = 0;
  long long rt_best_rejected_altitude_cases = 0;
  double rt_best_rejected_altitude_before = 0.0;
  double rt_best_rejected_altitude_gain = -std::numeric_limits<double>::infinity();
  double rt_best_rejected_altitude_after = 0.0;
  double rt_max_free_gap = 0.0;
  double rt_max_body_gap = 0.0;
  double rt_max_move_ratio = 0.0;
  // collapse 内訳
  double t_c_needs = 0., t_c_runtrials = 0., t_c_replace = 0., t_c_postop = 0.;
  // split 内訳
  double t_s_needs = 0., t_s_runtrials = 0., t_s_replace = 0., t_s_postop = 0.;
  // flip scenario polish 内訳
  double t_f_runtrials = 0.;
  using clk = std::chrono::high_resolution_clock;
  auto sec_since = [](const clk::time_point& t0) {
    return std::chrono::duration<double>(clk::now() - t0).count();
  };

  // Operation limits: -1 = unlimited, 0 = disabled, >0 = cap.
  const int total_ops_limit = rs.total_ops_limit;
  const int split_ops_limit = rs.split_ops_limit;
  const int collapse_ops_limit = rs.collapse_ops_limit;
  const int flip_ops_limit = rs.flip_ops_limit;
  std::size_t split_ops_count = 0;
  std::size_t collapse_ops_count = 0;
  std::size_t flip_ops_count = 0;
  bool topology_changed = false;
  auto op_limit_reached = [](const std::size_t used, const int limit) {
    return limit >= 0 && used >= static_cast<std::size_t>(limit);
  };
  auto int_limit_reached = [](const int used, const int limit) {
    return limit >= 0 && used >= limit;
  };
  auto total_ops_count = [&]() {
    return split_ops_count + collapse_ops_count + flip_ops_count;
  };
  auto total_ops_reached = [&]() { return op_limit_reached(total_ops_count(), total_ops_limit); };
  auto split_ops_reached = [&]() { return op_limit_reached(split_ops_count, split_ops_limit); };
  auto collapse_ops_reached = [&]() { return op_limit_reached(collapse_ops_count, collapse_ops_limit); };
  auto flip_ops_reached = [&]() { return op_limit_reached(flip_ops_count, flip_ops_limit); };
  auto all_ops_reached = [&]() {
    return total_ops_reached() ||
           (split_ops_reached() && collapse_ops_reached() && flip_ops_reached());
  };
  struct FinalLineRemapStats {
    int total = 0;
    int linear = 0;
    int projected = 0;
    int failed = 0;
    int scalar_updated = 0;
    int scalar_failed = 0;
    double max_move = 0.0;
  };

  auto finalLineTargetAccepted = [](const BEMMeshPipeline::GeometryTarget& target) {
    using BEMMeshPipeline::GeometryTargetStatus;
    return target.status == GeometryTargetStatus::ok ||
           target.status == GeometryTargetStatus::ok_via_fallback ||
           target.status == GeometryTargetStatus::no_convergence;
  };
  auto finalLineConstructionTargetAccepted = [&](const BEMMeshPipeline::GeometryTarget& target,
                                                 const Tddd& x_seed) {
    if (!finalLineTargetAccepted(target) || !isFinite(target.target_X))
      return false;
    if (std::isfinite(target.shift_limit) && target.shift_limit > 0.0) {
      const double move = Norm(target.target_X - x_seed);
      if (!std::isfinite(move) || move > target.shift_limit)
        return false;
    }
    if (std::isfinite(target.body_gap_threshold) && target.body_gap_threshold > 0.0 &&
        std::isfinite(target.body_gap) && target.body_gap > target.body_gap_threshold)
      return false;
    return true;
  };

  auto finalRemapLineNodes = [&]() {
    FinalLineRemapStats stats;
    const bool true_quadratic = use_true_quadratic_element;
    auto commit_line = [&](networkLine* l,
                           const Tddd& x_linear,
                           const Tddd& x_final,
                           bool projected,
                           const BEMMeshPipeline::BoundaryGeometryTargetProvider* target_provider = nullptr,
                           const BEMMeshPipeline::GeometryTarget* target = nullptr) {
      if (!l || !isFinite(x_final))
        return;
      const double move = isFinite(l->X_mid) ? Norm(x_final - l->X_mid) : 0.0;
      l->setXSingle(x_final);
      l->corner_midpoint_offset = l->X_mid - x_linear;
      if (hasAnyDirichletBoundaryState(l)) {
        if (target_provider && target) {
          const auto phi =
              target_provider->sampleReferenceField(*target, BEMMeshPipeline::ProviderFieldKind::Phi);
          if (phi.ok) {
            std::get<0>(l->phiphin) = phi.value;
            const auto phi_t =
                target_provider->sampleReferenceField(*target, BEMMeshPipeline::ProviderFieldKind::PhiT);
            if (phi_t.ok)
              std::get<0>(l->phiphin_t) = phi_t.value;
            ++stats.scalar_updated;
          } else {
            ++stats.scalar_failed;
          }
        } else if (!projected) {
          auto [pA, pB] = l->getPoints();
          if (pA && pB) {
            std::get<0>(l->phiphin) =
                0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
            std::get<0>(l->phiphin_t) =
                0.5 * (std::get<0>(pA->phiphin_t) + std::get<0>(pB->phiphin_t));
            ++stats.scalar_updated;
          } else {
            ++stats.scalar_failed;
          }
        }
      }
      stats.max_move = std::max(stats.max_move, move);
      if (projected)
        ++stats.projected;
      else
        ++stats.linear;
    };
    auto is_repairable_line = [](const networkLine* l) {
      return l && (hasAnyDirichletBoundaryState(l) ||
                   hasAnyNeumannBoundaryState(l) ||
                   BEMMeshPipeline::GeometryProjectorDetail::isBCInterfaceEntity(l));
    };

    if (true_quadratic && trial_geometry_inputs && trial_geometry_inputs->enabled()) {
      BEMMeshPipeline::GeometryProjectorOptions options;
      options.mode = trial_geometry_inputs->projection_mode;
      options.move_limit_factor = trial_geometry_inputs->mesh_pipeline->clung_move_limit_factor;
      options.max_iter = 4;
      options.tol_relative = 1e-6;
      options.max_body_gap_factor = 5.0;
      options.contact_range_factor = 2.0;
      options.trust_radius_factor = 5.0;
      options.body_bfs_fallback_depth = 1;
      options.enable_global_body_fallback = true;
      options.global_body_fallback_range_factor = 1.0;
      options.feature_angle_rad = rs.feature_angle_deg * M_PI / 180.0;
      const BEMMeshPipeline::BoundaryGeometryTargetProvider target_provider{
          *trial_geometry_inputs->reference, *trial_geometry_inputs->body_objects, options,
          trial_geometry_inputs->provider_epoch};
      for (auto* l : water.getBoundaryLines()) {
        if (!l || l->getBoundaryFaces().size() < 2)
          continue;
        auto [pA, pB] = l->getPoints();
        if (!pA || !pB)
          continue;
        ++stats.total;
        const Tddd x_linear = 0.5 * (pA->X + pB->X);
        if (!is_repairable_line(l)) {
          commit_line(l, x_linear, x_linear, /*projected=*/false);
          continue;
        }
        const auto target = target_provider.queryLineSampleTarget(l, &water, 0.5, x_linear);
        if (!finalLineConstructionTargetAccepted(target, x_linear)) {
          ++stats.failed;
          continue;
        }
        commit_line(l, x_linear, target.target_X, /*projected=*/true, &target_provider, &target);
      }
    } else {
      for (auto* l : water.getBoundaryLines()) {
        if (!l || l->getBoundaryFaces().size() < 2)
          continue;
        auto [pA, pB] = l->getPoints();
        if (!pA || !pB)
          continue;
        ++stats.total;
        const Tddd x_linear = 0.5 * (pA->X + pB->X);
        commit_line(l, x_linear, x_linear, /*projected=*/false);
      }
    }
    water.setGeometricPropertiesForce();
    refresh_contact_state_after_topology_change("final-line-remap",
                                                /*hard_fail_on_no_contact=*/false);
    std::cout << "[remesh:final_line_remap] step=" << time_step
              << " true_quadratic=" << (true_quadratic ? 1 : 0)
              << " total=" << stats.total
              << " projected=" << stats.projected
              << " linear=" << stats.linear
              << " failed=" << stats.failed
              << " scalar_updated=" << stats.scalar_updated
              << " scalar_failed=" << stats.scalar_failed
              << " max_move=" << stats.max_move
              << std::endl;
    return stats;
  };

  for (auto i = 0; i < iter_split_collapse; i++) {
    if (all_ops_reached()) {
      std::cout << "[remesh:limit] per_op_limits"
                << " total=" << total_ops_count() << "/" << total_ops_limit
                << " split=" << split_ops_count << "/" << split_ops_limit
                << " collapse=" << collapse_ops_count << "/" << collapse_ops_limit
                << " flip=" << flip_ops_count << "/" << flip_ops_limit
                << " skip_remaining_from_iter=" << i
                << std::endl;
      break;
    }
    auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };
    struct TrainingIterationStats {
      int candidate_edges = 0;
      int trial_count = 0;
      int safe_count = 0;
      int selected_count = 0;
      int applied_count = 0;
      int failed_apply_count = 0;
      int saved_success = 0;
      int saved_failure = 0;
    };
    TrainingIterationStats training_iter_stats;
    auto trainingCandidateLimitReached = [&]() {
      return training_data_record_enabled &&
             int_limit_reached(training_iter_stats.candidate_edges,
                               rs.training_data_max_edges_per_iter);
    };
    auto trainingApplyLimitReached = [&]() {
      return training_data_record_enabled &&
             int_limit_reached(training_iter_stats.applied_count,
                               rs.training_data_max_apply_per_iter);
    };
    auto endpointIndex0 = [](const networkLine* l) {
      if (!l)
        return -1;
      auto [p0, p1] = l->getPoints();
      return p0 ? p0->index : -1;
    };
    auto endpointIndex1 = [](const networkLine* l) {
      if (!l)
        return -1;
      auto [p0, p1] = l->getPoints();
      return p1 ? p1->index : -1;
    };
    auto targetLengthFor = [&](const networkLine* l) {
      return (l && (l->Dirichlet || l->BCInterface)) ? free_surface_target_len : global_mean_len;
    };
    constexpr std::array<int, 3> remesh_region_order = {2, 1, 0}; // BCInterface, Dirichlet, Neumann
    auto remeshRegionName = [](int category) {
      switch (category) {
      case 2:
        return "BCInterface";
      case 1:
        return "Dirichlet";
      case 0:
        return "Neumann";
      default:
        return "Other";
      }
    };
    auto keepRegionCandidates = [&](std::vector<networkLine*>& candidates, int category) {
      std::erase_if(candidates, [&](const networkLine* l) {
        return edgeTrialCategory(l) != category;
      });
    };
    auto trainingShuffleTail = [&](std::vector<networkLine*>& candidates,
                                   int remesh_iter,
                                   int op_code,
                                   int category) {
      if (!training_data_record_enabled || candidates.size() <= 1)
        return;
      const int random_extra = std::max(0, rs.training_data_random_extra_edges_per_iter);
      if (random_extra == 0)
        return;
      std::size_t priority_keep = 0;
      if (rs.training_data_max_edges_per_iter > 0) {
        const int keep = std::max(0, rs.training_data_max_edges_per_iter - random_extra);
        priority_keep = std::min<std::size_t>(static_cast<std::size_t>(keep), candidates.size());
      }
      auto mix = [](std::uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
      };
      const std::uint64_t seed =
          mix(static_cast<std::uint64_t>(rs.training_data_random_seed)) ^
          mix(static_cast<std::uint64_t>(time_step)) ^
          mix(static_cast<std::uint64_t>(remesh_iter + 1) << 16) ^
          mix(static_cast<std::uint64_t>(op_code + 17) << 32) ^
          mix(static_cast<std::uint64_t>(category + 31) << 48);
      std::mt19937_64 rng(seed);
      std::shuffle(candidates.begin() + static_cast<std::ptrdiff_t>(priority_keep),
                   candidates.end(), rng);
    };
    auto is_moving_contact_body = [](const Network* net) {
      if (!net)
        return false;
      const auto mode = std::get<0>(net->velocity_name_start);
      if (!mode.empty() && mode != "fixed")
        return true;
      constexpr double eps = 1e-14;
      return Norm(net->velocityTranslational()) > eps ||
             Norm(net->velocityRotational()) > eps;
    };
    auto is_moving_body_neumann_line = [&](const networkLine* l) {
      if (!l || !l->Neumann)
        return false;
      for (auto* cf : getEffectiveContactFaces(l)) {
        if (cf && is_moving_contact_body(cf->getNetwork()))
          return true;
      }
      return false;
    };
    // feature 保護の閾値（これ以上の法線変化は feature として保護）
    const double feature_angle = rs.feature_angle_deg * M_PI / 180.0;

    auto makeRemeshTargetProvider =
        [&](const BEMMeshPipeline::RemeshTrialGeometryInputs& c) {
          BEMMeshPipeline::GeometryProjectorOptions options;
          options.mode = c.projection_mode;
          options.move_limit_factor = c.mesh_pipeline->clung_move_limit_factor;
          options.max_iter = 4;
          options.tol_relative = 1e-6;
          options.max_body_gap_factor = 5.0;
          options.contact_range_factor = 2.0;
          options.trust_radius_factor = 5.0;
          options.body_bfs_fallback_depth = 1;
          options.enable_global_body_fallback = true;
          options.global_body_fallback_range_factor = 1.0;
          options.feature_angle_rad = feature_angle;
          return BEMMeshPipeline::BoundaryGeometryTargetProvider{*c.reference, *c.body_objects, options, c.provider_epoch};
        };

    enum class ThetaSource {
      Invalid,
      Provider,
      CandidateFallback
    };
    struct EdgeThetaSample {
      bool valid = false;
      double theta = -1.0;
      ThetaSource source = ThetaSource::Invalid;
      bool used_reference = false;
      bool used_body = false;
    };
    struct ProviderThetaStats {
      long long all_provider = 0;
      long long all_invalid = 0;
      long long all_candidate_fallback = 0;
    } provider_theta_stats;
    auto candidateMeshTheta = [](const networkLine* l, const Tddd& edge) -> double {
      if (!l || !l->geom_curvature.valid)
        return -1.0;
      return surface_geometry::edgeCurvatureAngle(
          edge, l->geom_curvature.k1, l->geom_curvature.k2,
          l->geom_curvature.PD1, l->geom_curvature.PD2);
    };
    auto targetThetaForEdge = [&](Network& net,
                                  networkLine* l,
                                  bool allow_candidate_fallback,
                                  bool record_stats) {
      EdgeThetaSample out;
      if (!curvature_remesh_enabled || !l)
        return out;
      auto [p0, p1] = l->getPoints();
      if (!p0 || !p1)
        return out;
      const Tddd edge = p1->X - p0->X;
      if (!isFinite(edge) || Norm(edge) <= 1e-14)
        return out;

      if (trial_geometry_inputs && trial_geometry_inputs->enabled()) {
        const auto target_provider = makeRemeshTargetProvider(*trial_geometry_inputs);
        auto account = [&](const BEMMeshPipeline::TargetThetaResult& r) {
          if (r.valid && std::isfinite(r.theta) && r.theta > 0.0) {
            out.valid = true;
            out.theta = std::max(out.theta, r.theta);
            out.source = ThetaSource::Provider;
            out.used_reference = out.used_reference || r.used_reference;
            out.used_body = out.used_body || r.used_body;
          }
        };
        const Tddd x_mid = 0.5 * (p0->X + p1->X);
        account(target_provider.queryTargetTheta(l, &net, x_mid, edge));
        account(target_provider.queryTargetTheta(p0, &net, p0->X, edge));
        account(target_provider.queryTargetTheta(p1, &net, p1->X, edge));
      }
      if (!out.valid && allow_candidate_fallback) {
        const double theta = candidateMeshTheta(l, edge);
        if (std::isfinite(theta) && theta > 0.0) {
          out.valid = true;
          out.theta = theta;
          out.source = ThetaSource::CandidateFallback;
        }
      }
      if (record_stats) {
        if (out.source == ThetaSource::Provider) {
          ++provider_theta_stats.all_provider;
        } else if (out.source == ThetaSource::CandidateFallback) {
          ++provider_theta_stats.all_candidate_fallback;
        } else {
          ++provider_theta_stats.all_invalid;
        }
      }
      return out;
    };
    auto edgeThetaValue = [&](networkLine* l) -> double {
      if (!l)
        return -1.0;
      const auto theta = targetThetaForEdge(water, l,
                                            /*allow_candidate_fallback=*/true,
                                            /*record_stats=*/false);
      return theta.valid ? theta.theta : -1.0;
    };
	    auto selectLightGBMScenarios = [&](networkLine* target_line,
	                                       const std::vector<std::string>& scenarios,
	                                       int op_code,
	                                       int reason_code) {
      if (!use_lightgbm_ranking ||
          !lightgbm_scenario_model.loaded() ||
          scenarios.size() <= 1 ||
          (op_code != 0 && op_code != 1 && op_code != 2))
        return scenarios;

      const double len = target_line ? target_line->length() : -1.0;
      const double target_len = targetLengthFor(target_line);
      const double len_ratio = (target_len > 1e-20) ? len / target_len : -1.0;
      const double theta = edgeThetaValue(target_line);
      const double theta_ratio =
          (theta >= 0.0 && theta_split > 1e-20) ? theta / theta_split : -1.0;
      const int category = edgeTrialCategory(target_line);

      struct ScoredScenario {
        std::string scenario;
        double score;
        std::size_t original_index;
      };
      std::vector<ScoredScenario> scored;
      scored.reserve(scenarios.size());
      for (std::size_t i = 0; i < scenarios.size(); ++i) {
        BEMMeshPipeline::RemeshAI::RemeshScenarioFeatures f;
        f.op = op_code;
        f.reason = reason_code;
        f.edge_category = category;
        f.len = len;
        f.len_ratio = len_ratio;
        f.theta = theta;
        f.theta_ratio = theta_ratio;
        f.scenario = scenarios[i];
        scored.push_back({scenarios[i],
                          lightgbm_scenario_model.predictProbability(f),
                          i});
      }

      std::stable_sort(scored.begin(), scored.end(),
                       [](const ScoredScenario& a, const ScoredScenario& b) {
                         if (a.score != b.score)
                           return a.score > b.score;
                         return a.original_index < b.original_index;
                       });
      const std::size_t selected_count =
          use_lightgbm_recommender
              ? std::min<std::size_t>(static_cast<std::size_t>(rs.lightgbm_recommender_top_k),
                                      scored.size())
              : scored.size();
      std::vector<std::string> selected;
      selected.reserve(selected_count);
      for (std::size_t i = 0; i < selected_count; ++i)
        selected.push_back(scored[i].scenario);

      if (selected.empty() && !scored.empty())
        selected.push_back(scored.front().scenario);

      if (selected.empty())
        return scenarios;

      std::ostringstream oss;
      for (std::size_t i = 0; i < scored.size(); ++i) {
        if (i > 0)
          oss << ',';
        oss << scored[i].scenario << ':' << scored[i].score;
      }
      std::cout << "[remesh:lgbm_scenarios] step=" << time_step
                << " mode=" << rs.lightgbm_mode
                << " op=" << op_code
                << " reason=" << reason_code
                << " category=" << category
                << " p0=" << endpointIndex0(target_line)
                << " p1=" << endpointIndex1(target_line)
                << (use_lightgbm_recommender ? " recommended=" : " ranked=")
                << selected.size() << "/" << scenarios.size()
                << " scores={" << oss.str() << "}"
                << std::endl;

	      return selected;
	    };
	    auto scenarioCountForEdgeTrialLog = [&](const std::vector<std::string>& scenarios) {
	      if (!use_lightgbm_recommender)
	        return static_cast<int>(scenarios.size());
	      return static_cast<int>(std::min<std::size_t>(
	          static_cast<std::size_t>(rs.lightgbm_recommender_top_k),
	          scenarios.size()));
	    };
	    auto lightGBMEdgeScore = [&](networkLine* target_line,
	                                 const std::vector<std::string>& scenarios,
	                                 int op_code,
	                                 int reason_code) {
      if (!use_lightgbm_ranking ||
          !lightgbm_scenario_model.loaded() ||
          scenarios.empty() ||
          (op_code != 0 && op_code != 1 && op_code != 2))
        return -std::numeric_limits<double>::infinity();

      const double len = target_line ? target_line->length() : -1.0;
      const double target_len = targetLengthFor(target_line);
      const double len_ratio = (target_len > 1e-20) ? len / target_len : -1.0;
      const double theta = edgeThetaValue(target_line);
      const double theta_ratio =
          (theta >= 0.0 && theta_split > 1e-20) ? theta / theta_split : -1.0;
      const int category = edgeTrialCategory(target_line);

      double best = -std::numeric_limits<double>::infinity();
      for (const auto& scenario : scenarios) {
        BEMMeshPipeline::RemeshAI::RemeshScenarioFeatures f;
        f.op = op_code;
        f.reason = reason_code;
        f.edge_category = category;
        f.len = len;
        f.len_ratio = len_ratio;
        f.theta = theta;
        f.theta_ratio = theta_ratio;
        f.scenario = scenario;
        best = std::max(best, lightgbm_scenario_model.predictProbability(f));
      }
      return best;
    };
    auto appendProviderThetaStatsCsv = [&]() {
      if (patch_output_directory.empty())
        return;
      std::filesystem::create_directories(patch_output_directory);
      const std::filesystem::path path =
          std::filesystem::path(patch_output_directory) / "provider_theta_stats.csv";
      const bool write_header = !std::filesystem::exists(path) ||
                                std::filesystem::file_size(path) == 0;
      std::ofstream out(path, std::ios::app);
      if (!out)
        return;
      if (write_header) {
        out << "time_step,phase,all_provider,all_invalid,all_candidate_fallback\n";
      }
      out << time_step << ',' << sanitize_phase_tag(phase_debug_tag) << ','
          << provider_theta_stats.all_provider << ','
          << provider_theta_stats.all_invalid << ','
          << provider_theta_stats.all_candidate_fallback << '\n';
      std::cout << "[remesh:provider_theta] step=" << time_step
                << " phase=" << sanitize_phase_tag(phase_debug_tag)
                << " all_provider=" << provider_theta_stats.all_provider
                << " all_invalid=" << provider_theta_stats.all_invalid
                << " all_candidate_fallback=" << provider_theta_stats.all_candidate_fallback
                << std::endl;
    };

    // --- θ 判定（split/collapse 共通） ---
    // 辺中点の曲率を使い、端点曲率の max による過剰 split を回避
    auto edgeThetaVerdict = [&](const networkLine* l) -> EdgeThetaVerdict {
      if (!curvature_remesh_enabled)
        return EdgeThetaVerdict::Keep;
      auto* mutable_line = const_cast<networkLine*>(l);
      const auto theta_sample = targetThetaForEdge(water, mutable_line,
                                                   /*allow_candidate_fallback=*/true,
                                                   /*record_stats=*/true);
      if (!theta_sample.valid)
        return EdgeThetaVerdict::CurvatureInvalid;
      const double theta = theta_sample.theta;
      if (theta >= theta_split)
        return EdgeThetaVerdict::SplitCandidate;
      if (theta < theta_collapse)
        return EdgeThetaVerdict::CollapseCandidate;
      return EdgeThetaVerdict::Keep;
    };

    // --- needsSplit: 解像度的にこの辺は粗すぎるか？ ---
    // 優先順位:
    //   1. global 上限超え → 無条件 split
    //   2. θ（曲面忠実度）→ SplitCandidate なら split
    //   3. grading（メッシュ均一性）→ 周囲より突出して長い辺を split
    //   4. CurvatureInvalid → 二面角 fallback
    // 辺の隣接面に鈍角（閾値超え）の内角があるか
    auto hasObtuseAngle = [](const networkLine* l, double angle_threshold = M_PI / 2.) -> bool {
      for (auto* f : l->getBoundaryFaces()) {
        if (!f)
          continue;
        auto [p0, p1, p2] = f->getPoints();
        for (auto& [a, b, c] : std::vector<std::tuple<Tddd, Tddd, Tddd>>{{p0->X, p1->X, p2->X}, {p1->X, p2->X, p0->X}, {p2->X, p0->X, p1->X}}) {
          Tddd v1 = b - a, v2 = c - a;
          double len1 = Norm(v1), len2 = Norm(v2);
          if (len1 > 1e-20 && len2 > 1e-20) {
            double cos_angle = std::clamp(Dot(v1, v2) / (len1 * len2), -1.0, 1.0);
            if (std::acos(cos_angle) > angle_threshold)
              return true;
          }
        }
      }
      return false;
    };

    // --- needsSplit / needsCollapse: この辺を split / collapse すべきか ---
    // Trial では標準的な split / collapse 閾値を使う:
    //   split   : len > 4/3·target
    //   collapse: len < 4/5·target
    // target は自由表面 (Dirichlet/BCInterface) なら free_surface_target_len、
    // それ以外は global_mean_len。
    auto needsSplit = [&](networkLine* l) -> bool {
      const auto len = l->length();
      if (!std::isfinite(len) || len <= 0.)
        return false;
      const double target = (l->Dirichlet || l->BCInterface)
                                ? free_surface_target_len
                                : global_mean_len;
      if (!(target > 0.))
        return false;
      if (edgeThetaVerdict(l) == EdgeThetaVerdict::SplitCandidate)
        return true;
      return len > (4.0 / 3.0) * target;
    };

    auto needsCollapse = [&](networkLine* l) -> bool {
      if (!l || l->BCInterface)
        return false;
      const auto len = l->length();
      if (!std::isfinite(len) || len <= 0.)
        return false;
      const double target = (l->Dirichlet || l->BCInterface)
                                ? free_surface_target_len
                                : global_mean_len;
      if (!(target > 0.))
        return false;
      return len < (4.0 / 5.0) * target;
    };

    auto splitReasonCodeForLine = [&](networkLine* l,
                                      const std::unordered_set<networkLine*>& subsurface_trigger_set) {
      if (subsurface_trigger_set.contains(l))
        return 5;
      const double len = l ? l->length() : -1.0;
      if (len > global_max_len)
        return 0;
      const auto v = edgeThetaVerdict(l);
      if (v == EdgeThetaVerdict::SplitCandidate)
        return 1;
      if (v == EdgeThetaVerdict::CurvatureInvalid)
        return 2;
      if (hasObtuseAngle(l))
        return 3;
      return 4;
    };
    auto splitReasonName = [](int reason_code) {
      switch (reason_code) {
      case 0:
        return "global_max";
      case 1:
        return "theta";
      case 2:
        return "dihedral";
      case 3:
        return "obtuse";
      case 5:
        return "subsurface_altitude";
      default:
        return "length";
      }
    };
    auto collapseReasonCodeForLine = [&](networkLine* l) {
      const double len = l ? l->length() : -1.0;
      if (len < global_min_len)
        return 0;
      if (hasObtuseAngle(l))
        return 4;
      const auto v = edgeThetaVerdict(l);
      if (v == EdgeThetaVerdict::CollapseCandidate)
        return 1;
      if (v == EdgeThetaVerdict::CurvatureInvalid)
        return 2;
      return 3;
    };

    struct SubsurfaceTriggerSet {
      std::vector<networkLine*> lines;
      std::size_t bad_faces = 0;
      double min_altitude_rel = 1E+100;
    };

    auto collectSubsurfaceTriggerLines = [&]() {
      SubsurfaceTriggerSet out;
      if (!rs.subsurface_altitude_reject.enabled)
        return out;
      std::unordered_set<networkLine*> seen;
      auto add_line = [&](networkLine* l) {
        if (!line_alive(l) || seen.contains(l))
          return;
        seen.insert(l);
        out.lines.emplace_back(l);
      };
      auto add_face_trigger = [&](networkFace* f) {
        if (!f)
          return;
        ++out.bad_faces;
        bool added_interface = false;
        for (auto* fl : f->getLines()) {
          if (line_alive(fl) && fl->BCInterface) {
            add_line(fl);
            added_interface = true;
          }
        }
        if (added_interface)
          return;
        networkLine* longest = nullptr;
        double longest_len = -1.0;
        for (auto* fl : f->getLines()) {
          if (!line_alive(fl))
            continue;
          const double len = fl->length();
          if (std::isfinite(len) && len > longest_len) {
            longest = fl;
            longest_len = len;
          }
        }
        add_line(longest);
      };

      for (auto* line : water.getBoundaryLines()) {
        if (!line_alive(line))
          continue;
        auto faces = line->getBoundaryFaces();
        if (faces.size() != 2 || !faces[0] || !faces[1])
          continue;
        const double dot = std::clamp(Dot(faces[0]->normal, faces[1]->normal), -1., 1.);
        const double theta_deg = 180. - std::acos(dot) * 180. / M_PI;
        if (!(rs.subsurface_altitude_reject.min_edge_angle_deg < theta_deg &&
              theta_deg < rs.subsurface_altitude_reject.max_edge_angle_deg))
          continue;
        const auto d0 = faceAltitudeRelativeToSharedEdgeDetail(
            faces[0], line, /*use_quadratic_subfaces=*/false);
        const auto d1 = faceAltitudeRelativeToSharedEdgeDetail(
            faces[1], line, /*use_quadratic_subfaces=*/false);
        if (d0.altitude_rel < rs.subsurface_altitude_reject.min_face_altitude_rel) {
          out.min_altitude_rel = std::min(out.min_altitude_rel, d0.altitude_rel);
          add_face_trigger(faces[0]);
        }
        if (d1.altitude_rel < rs.subsurface_altitude_reject.min_face_altitude_rel) {
          out.min_altitude_rel = std::min(out.min_altitude_rel, d1.altitude_rel);
          add_face_trigger(faces[1]);
        }
      }
      return out;
    };

    // ------ 全辺の幾何/曲率情報を VTU 出力（i==0 のみ） ------
    auto writeEdgesVTU = [&]() {
      if (!edges_pvd && patch_output_directory.empty())
        return;
      const std::string vtu_name = water.getName() + "_edges_" + std::to_string(time_step) + ".vtu";
      const std::string filename = patch_output_directory + "/" + vtu_name;
      FILE* fp = fopen(filename.c_str(), "wb");
      if (!fp)
        return;
      const auto lines = water.getBoundaryLines();
      int n_cells = (int)lines.size();
      int n_points = n_cells * 2;
      fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
      fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
      fprintf(fp, "<UnstructuredGrid>\n");
      fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", n_cells, n_points);
      fprintf(fp, "<Points>\n");
      fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' format='ascii'>\n");
      for (auto* l : lines) {
        auto [p0, p1] = l->getPoints();
        fprintf(fp, "%f %f %f\n", (float)std::get<0>(p0->X), (float)std::get<1>(p0->X), (float)std::get<2>(p0->X));
        fprintf(fp, "%f %f %f\n", (float)std::get<0>(p1->X), (float)std::get<1>(p1->X), (float)std::get<2>(p1->X));
      }
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "</Points>\n");
      fprintf(fp, "<Cells>\n");
      fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
      for (int i = 0; i < n_cells; ++i)
        fprintf(fp, "%d %d\n", i * 2, i * 2 + 1);
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
      for (int i = 0; i < n_cells; ++i)
        fprintf(fp, "%d\n", (i + 1) * 2);
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
      for (int i = 0; i < n_cells; ++i)
        fprintf(fp, "3\n"); // VTK_LINE
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "</Cells>\n");
      fprintf(fp, "<CellData>\n");

      auto write_scalar = [&](const char* name, auto extractor) {
        fprintf(fp, "<DataArray type='Float32' Name='%s' NumberOfComponents='1' format='ascii'>\n", name);
        for (auto* l : lines)
          fprintf(fp, "%f\n", (float)extractor(l));
        fprintf(fp, "</DataArray>\n");
      };
      auto write_vec3 = [&](const char* name, auto extractor) {
        fprintf(fp, "<DataArray type='Float32' Name='%s' NumberOfComponents='3' format='ascii'>\n", name);
        for (auto* l : lines) {
          auto v = extractor(l);
          fprintf(fp, "%f %f %f\n", (float)std::get<0>(v), (float)std::get<1>(v), (float)std::get<2>(v));
        }
        fprintf(fp, "</DataArray>\n");
      };

      write_scalar("k1", [](networkLine* l) { return l->geom_curvature.valid ? l->geom_curvature.k1 : 0.0; });
      write_scalar("k2", [](networkLine* l) { return l->geom_curvature.valid ? l->geom_curvature.k2 : 0.0; });
      write_scalar("kmax", [](networkLine* l) { return l->geom_curvature.valid ? l->geom_curvature.kmax : 0.0; });
      write_scalar("len", [](networkLine* l) { return l->length(); });
      write_scalar("theta", [](networkLine* l) {
        if (!l->geom_curvature.valid)
          return 0.0;
        auto [p0, p1] = l->getPoints();
        Tddd edge = p1->X - p0->X;
        return surface_geometry::edgeCurvatureAngle(edge, l->geom_curvature.k1, l->geom_curvature.k2, l->geom_curvature.PD1, l->geom_curvature.PD2);
      });
      write_scalar("valid", [](networkLine* l) { return l->geom_curvature.valid ? 1.0 : 0.0; });
      write_scalar("is_sharp", [](networkLine* l) { return l->SharpQ() ? 1.0 : 0.0; });
      write_scalar("is_split_target", [&](networkLine* l) { return needsSplit(l) ? 1.0 : 0.0; });
      write_scalar("is_collapse_target", [&](networkLine* l) { return needsCollapse(l) ? 1.0 : 0.0; });
      write_vec3("PD1", [](networkLine* l) -> Tddd { return l->geom_curvature.valid ? l->geom_curvature.PD1 : Tddd{0., 0., 0.}; });
      write_vec3("PD2", [](networkLine* l) -> Tddd { return l->geom_curvature.valid ? l->geom_curvature.PD2 : Tddd{0., 0., 0.}; });
      write_vec3("n_avg", [](networkLine* l) -> Tddd {
        auto faces = l->getBoundaryFaces();
        Tddd n = {0., 0., 0.};
        for (auto* f : faces)
          if (f)
            n += f->area * f->normal;
        double nn = Norm(n);
        return (nn > 1e-20) ? (n / nn) : Tddd{0., 0., 0.};
      });

      fprintf(fp, "</CellData>\n");
      fprintf(fp, "</Piece>\n");
      fprintf(fp, "</UnstructuredGrid>\n");
      fprintf(fp, "</VTKFile>\n");
      fclose(fp);
      std::cout << "[edges debug] wrote " << filename << " (" << n_cells << " edges)" << std::endl;
      if (edges_pvd) {
        edges_pvd->push(vtu_name, simulation_time);
        edges_pvd->output();
      }
    };

    if (remesh_debug_edges && i == 0)
      writeEdgesVTU();

    // ------ 共通ヘルパー（split / collapse 共用） ------

    struct ProviderFieldAssignment {
      bool update_phi = false;
      double phi = 0.0;
      bool update_phi_t = false;
      double phi_t = 0.0;
    };

    auto sampleProviderReferenceFields =
        [&](const BEMMeshPipeline::BoundaryGeometryTargetProvider& target_provider,
            const BEMMeshPipeline::GeometryTarget& target,
            const networkPoint* p,
            ProviderFieldAssignment& fields) {
          fields = {};
          if (rs.defer_scalar_interpolation || !hasAnyDirichletBoundaryState(p))
            return true;
          const auto phi = target_provider.sampleReferenceField(
              target, BEMMeshPipeline::ProviderFieldKind::Phi);
          if (!phi.ok)
            return false;
          fields.update_phi = true;
          fields.phi = phi.value;
          const auto phi_t = target_provider.sampleReferenceField(
              target, BEMMeshPipeline::ProviderFieldKind::PhiT);
          if (phi_t.ok) {
            fields.update_phi_t = true;
            fields.phi_t = phi_t.value;
          }
          return true;
        };

    auto applyProviderFieldAssignment = [](networkPoint* p,
                                           const ProviderFieldAssignment& fields) {
      if (!p)
        return;
      if (fields.update_phi)
        std::get<0>(p->phiphin) = fields.phi;
      if (fields.update_phi_t)
        std::get<0>(p->phiphin_t) = fields.phi_t;
    };

    struct VirtualRepairStats {
      int repaired_points = 0;
      int repaired_lines = 0;
      int repair_failed = 0;
      double max_free_gap = 0.0;
      double max_body_gap = 0.0;
      double max_move_ratio = 0.0;
    };

    auto repairTargetAccepted = [](const BEMMeshPipeline::GeometryTarget& target) {
      using BEMMeshPipeline::GeometryTargetStatus;
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
    };

    auto repairStatusIsFailure = [](BEMMeshPipeline::GeometryTargetStatus status) {
      using BEMMeshPipeline::GeometryTargetStatus;
      return status == GeometryTargetStatus::no_reference_surface ||
             status == GeometryTargetStatus::body_gap_too_large ||
             status == GeometryTargetStatus::invalid_projection;
    };

    auto accumulateRepairStats = [](VirtualRepairStats& stats, const BEMMeshPipeline::GeometryTarget& target) {
      if (std::isfinite(target.dirichlet_gap))
        stats.max_free_gap = std::max(stats.max_free_gap, target.dirichlet_gap);
      if (std::isfinite(target.body_gap))
        stats.max_body_gap = std::max(stats.max_body_gap, target.body_gap);
      if (std::isfinite(target.move_ratio))
        stats.max_move_ratio = std::max(stats.max_move_ratio, target.move_ratio);
    };

    auto mergeVirtualRepairStats = [](VirtualRepairStats& dst, const VirtualRepairStats& src) {
      dst.repaired_points += src.repaired_points;
      dst.repaired_lines += src.repaired_lines;
      dst.repair_failed += src.repair_failed;
      dst.max_free_gap = std::max(dst.max_free_gap, src.max_free_gap);
      dst.max_body_gap = std::max(dst.max_body_gap, src.max_body_gap);
      dst.max_move_ratio = std::max(dst.max_move_ratio, src.max_move_ratio);
    };

    auto meanPatchBoundaryEdgeLength = [](const Network& net) {
      double sum = 0.0;
      int count = 0;
      for (auto* l : net.getBoundaryLines()) {
        if (!l)
          continue;
        const double len = l->length();
        if (std::isfinite(len) && len > 0.0) {
          sum += len;
          ++count;
        }
      }
      return count > 0 ? sum / static_cast<double>(count) : 1.0;
    };

    auto applyVirtualGeometryRepair = [&](Network& net) {
      VirtualRepairStats stats;
      if (!trial_geometry_inputs || !trial_geometry_inputs->enabled())
        return stats;

      const auto target_provider = makeRemeshTargetProvider(*trial_geometry_inputs);
      const auto border_vec = net.getBorderPoints();
      std::unordered_set<const networkPoint*> border_points(border_vec.begin(), border_vec.end());
      auto is_border_point = [&](const networkPoint* p) {
        return p && border_points.contains(p);
      };
      auto is_repairable_point = [](const networkPoint* p) {
        return p && (hasAnyDirichletBoundaryState(p) ||
                     hasAnyNeumannBoundaryState(p) ||
                     BEMMeshPipeline::GeometryProjectorDetail::isBCInterfaceEntity(p));
      };
      const double local_mean = meanPatchBoundaryEdgeLength(net);
      const double contact_zero_tol =
          std::max(10.0 * trial_geometry_inputs->mesh_pipeline->clung_tol * local_mean, 1e-10);
      auto all_contact_faces_near_zero = [&](const auto* entity, const Tddd& x) {
        const auto faces = getEffectiveContactFaces(entity);
        if (faces.empty())
          return false;
        for (auto* f : faces) {
          if (!f)
            return false;
          const double d = Norm(Nearest(x, f) - x);
          if (!std::isfinite(d) || d > contact_zero_tol)
            return false;
        }
        return true;
      };
      const int max_iter = std::max(1, trial_geometry_inputs->mesh_pipeline->clung_max_iter);
      struct PointRepairMove {
        networkPoint* p = nullptr;
        Tddd delta = {0., 0., 0.};
        ProviderFieldAssignment fields;
      };
      for (int iter = 0; iter < max_iter; ++iter) {
        double max_move = 0.0;
        std::vector<PointRepairMove> moves;
        for (auto* p : net.getBoundaryPoints()) {
          if (!p || is_border_point(p))
            continue;
          if (!is_repairable_point(p))
            continue;
          if (all_contact_faces_near_zero(p, p->X))
            continue;
          const auto target = target_provider.queryEntityTarget(p, &net, p->X);
          if (!repairTargetAccepted(target)) {
            if (target.status == BEMMeshPipeline::GeometryTargetStatus::no_body_surface)
              continue;
            if (repairStatusIsFailure(target.status))
              ++stats.repair_failed;
            continue;
          }
          accumulateRepairStats(stats, target);
          ProviderFieldAssignment fields;
          if (!sampleProviderReferenceFields(target_provider, target, p, fields)) {
            ++stats.repair_failed;
            continue;
          }
          const double move = Norm(target.delta_clamped);
          if (!(move > 1e-12) && !fields.update_phi && !fields.update_phi_t)
            continue;
          moves.push_back({p, target.delta_clamped, fields});
          max_move = std::max(max_move, move);
        }
        for (auto& move : moves) {
          auto* p = move.p;
          if (!p)
            continue;
          if (Norm(move.delta) > 1e-12)
            p->setXSingle(p->X + move.delta);
          applyProviderFieldAssignment(p, move.fields);
          ++stats.repaired_points;
        }
        if (!moves.empty()) {
          for (const auto& l : net.getLines())
            l->setBoundsSingle();
          for (const auto& p : net.getPoints())
            p->setFaces();
          for (const auto& f : net.getFaces())
            f->setGeometricProperties(ToX(f->getPoints()));
          // Point motion changes distances and normals used by the per
          // (node, face) contact states. Refresh before the next repair
          // iteration.
          refresh_trial_patch_boundary_state(net);
        }

        // Trial repair intentionally does not move line X_mid. Line nodes are
        // mapped once after vertex repair/topology changes are finished.
        if (max_move < trial_geometry_inputs->mesh_pipeline->clung_tol * local_mean)
          break;
      }

      for (const auto& l : net.getLines())
        l->setBoundsSingle();
      for (const auto& p : net.getPoints())
        p->setFaces();
      for (const auto& f : net.getFaces())
        f->setGeometricProperties(ToX(f->getPoints()));
      net.computePrincipalCurvatures(false);
      return stats;
    };

    // 三角形の内角の 60° からの最大偏差 [rad]（0 が最良、π/3 が最悪）
    auto maxAngleDeviation = [](const Tddd& a, const Tddd& b, const Tddd& c) -> double {
      constexpr double ideal = M_PI / 3.;
      double max_dev = 0.;
      for (auto& [p, q, r] : std::vector<std::tuple<Tddd, Tddd, Tddd>>{{a, b, c}, {b, c, a}, {c, a, b}}) {
        Tddd v1 = q - p, v2 = r - p;
        double l1 = Norm(v1), l2 = Norm(v2);
        if (l1 > 1e-20 && l2 > 1e-20) {
          double angle = std::acos(std::clamp(Dot(v1, v2) / (l1 * l2), -1.0, 1.0));
          max_dev = std::max(max_dev, std::abs(angle - ideal));
        }
      }
      return max_dev;
    };

    // ============================================================
    // patchQuality: パッチの品質を 1 つの数値で評価する複合スコア
    // ============================================================
    //
    // ★ 高いほど良い（HIGHER IS BETTER）
    //   - 最悪（メッシュ反転）= 非常に小さい値（-1e18）で即却下扱い
    //
    // スコア式:
    //   score = worst_ir - w * (area_deviation + equilateral_deviation)
    //
    //   fs_target_dist    = mean_{Dirichlet/BCInterface edges} ( log(len / fs_target) )^2
    //   theta_target_dist = mean_{valid curvature edges}  ( log(theta / theta_target) )^2
    //   area_deviation   = mean_{faces} ( log(area / local_mean_area) )^2
    //   equilateral_deviation = mean_{faces,vertices} |X_ideal - X|^2 / h_ideal^2
    //
    // 特徴:
    //   - 過大・過小の両方を penalize（log space で対称）
    //   - 大きな逸脱は 2 乗で強く penalize（旧 1/(1+x) の飽和を回避）
    //   - mean (not max) なので 1 本の悪い辺に支配されず、全体が target に寄る方向へ
    //   - w = rs.quality_resolution_weight（既定 1.0）で重みを調整可
    //
    // 注意:
    //   patch 外縁も score に含める。外縁を除外すると、境界近傍で
    //   細長い三角形を作る候補が core の改善だけで accept されるため。
    auto worstSubsurfaceAltitudeRel = [&](Network& net) {
      if (!rs.subsurface_altitude_reject.enabled)
        return 1E+100;
      const auto result = checkSubsurfaceFaceAltitude(
          net, rs.subsurface_altitude_reject, /*use_quadratic_subfaces=*/false);
      if (!result.worst_line || !std::isfinite(result.worst_altitude_rel))
        return 1E+100;
      return result.worst_altitude_rel;
    };

    struct PatchQualityResult {
      double score = -1.0e18;
      double mean_ir = 0.0;
      double worst_ir = 0.0;
      int face_count = 0;
      double local_cv_mean = 0.0;
      double fs_target_dist = 0.0;
      double theta_target_dist = 0.0;
      double area_deviation = 0.0;
      double equilateral_deviation = 0.0;
    };

    auto patchQuality = [&](Network& net) {
      constexpr double kFlipSentinel = -1.0e18;
      PatchQualityResult q;
      q.score = kFlipSentinel;

      // ----- 反転チェック: 隣接面の法線が逆を向いていれば即却下 -----
      for (auto* l : net.getLines()) {
        auto faces = l->getBoundaryFaces();
        if (faces.size() == 2 && faces[0] && faces[1]) {
          if (Dot(faces[0]->normal, faces[1]->normal) < rs.quality_normal_flip_cos)
            return q;
        }
      }

      // ----- 面ループ: 形状品質（patch 外縁も含める）-----
      double worst_ir = 1.0;
      double ir_sum = 0.0;
      double area_dev_sum = 0.0;
      int area_dev_count = 0;
      double equilateral_dev_sum = 0.0;
      int equilateral_dev_count = 0;
      int face_count = 0;
      for (auto* f : net.getBoundaryFaces()) {
        if (!f)
          continue;
        auto [p0, p1, p2] = f->getPoints();
        if (!p0 || !p1 || !p2)
          continue;
        const double ir = std::pow(InradiusToCircumradius(p0->X, p1->X, p2->X), 2.);
        worst_ir = std::min(worst_ir, ir);
        ir_sum += ir;
        const double area = TriangleArea(p0->X, p1->X, p2->X);
        const double local_mean_area = localMeanFaceArea(f);
        if (area > 0.0 && local_mean_area > 0.0 &&
            std::isfinite(area) && std::isfinite(local_mean_area)) {
          const double log_area_ratio = std::log(area / local_mean_area);
          area_dev_sum += log_area_ratio * log_area_ratio;
          ++area_dev_count;
        }
        const double equilateral_dev = equilateral_coordinate_deviation(p0->X, p1->X, p2->X);
        if (std::isfinite(equilateral_dev)) {
          equilateral_dev_sum += equilateral_dev;
          ++equilateral_dev_count;
        }
        ++face_count;
      }
      if (face_count == 0)
        return q;
      q.face_count = face_count;
      q.worst_ir = worst_ir;
      q.mean_ir = ir_sum / static_cast<double>(face_count);
      q.area_deviation = (area_dev_count > 0)
                             ? area_dev_sum / static_cast<double>(area_dev_count)
                             : 0.0;
      q.equilateral_deviation = (equilateral_dev_count > 0)
                                    ? equilateral_dev_sum / static_cast<double>(equilateral_dev_count)
                                    : std::numeric_limits<double>::infinity();

      auto targetSurfaceTheta = [&](networkLine* l, networkPoint* p0, networkPoint* p1, const Tddd& edge) -> double {
        if (!curvature_remesh_enabled)
          return -1.0;
        if (!l || !p0 || !p1)
          return -1.0;
        const bool allow_candidate_fallback =
            !(trial_geometry_inputs && trial_geometry_inputs->enabled());
        EdgeThetaSample theta;
        if (trial_geometry_inputs && trial_geometry_inputs->enabled()) {
          const auto target_provider = makeRemeshTargetProvider(*trial_geometry_inputs);
          auto account = [&](const BEMMeshPipeline::TargetThetaResult& r) {
            if (r.valid && std::isfinite(r.theta) && r.theta > 0.0) {
              theta.valid = true;
              theta.theta = std::max(theta.theta, r.theta);
              theta.source = ThetaSource::Provider;
            }
          };
          account(target_provider.queryTargetTheta(p0, &net, p0->X, edge));
          account(target_provider.queryTargetTheta(p1, &net, p1->X, edge));
        }
        if (!theta.valid && allow_candidate_fallback) {
          const double candidate_theta = candidateMeshTheta(l, edge);
          if (std::isfinite(candidate_theta) && candidate_theta > 0.0) {
            theta.valid = true;
            theta.theta = candidate_theta;
            theta.source = ThetaSource::CandidateFallback;
          }
        }
        return theta.valid ? theta.theta : -1.0;
      };

      // ----- 辺ループ: 目標距離 (log^2 ratio の平均、patch 外縁も含める)-----
      double fs_dist_sum = 0.;
      int fs_count = 0;
      double th_dist_sum = 0.;
      int th_count = 0;
      for (auto* l : net.getBoundaryLines()) {
        auto [p0, p1] = l->getPoints();

        // Dirichlet/BCInterface 辺 → free surface target との距離
        if ((l->Dirichlet || l->BCInterface) && free_surface_target_len > 1e-20) {
          double L = l->length();
          if (L > 1e-20) {
            double lr = std::log(L / free_surface_target_len); // 長いほど正、短いほど負
            fs_dist_sum += lr * lr;
            ++fs_count;
          }
        }
        if (!p0 || !p1)
          continue;
        Tddd edge = p1->X - p0->X;
        // 曲率 theta は candidate patch の edge 長さ・方向を使い、
        // 曲率だけを repair target surface から取得する。
        const double theta = targetSurfaceTheta(l, p0, p1, edge);
        if (curvature_remesh_enabled && std::isfinite(theta) && theta_target > 1e-10 && theta > 1e-10) {
          double lr = std::log(theta / theta_target);
          th_dist_sum += lr * lr;
          ++th_count;
        }
      }
      double fs_target_dist = (fs_count > 0) ? fs_dist_sum / fs_count : 0.0;
      double theta_target_dist = (th_count > 0) ? th_dist_sum / th_count : 0.0;
      q.fs_target_dist = fs_target_dist;
      q.theta_target_dist = theta_target_dist;

      // ----- 節点ループ: local CV (local grading、patch 外縁も含める)-----
      double local_cv_mean = 0.;
      {
        double cv_sum = 0.;
        int cv_count = 0;
        for (auto* p : net.getBoundaryPoints()) {
          if (!p)
            continue;
          cv_sum += p->localEdgeLengthCV();
          ++cv_count;
        }
        if (cv_count > 0)
          local_cv_mean = cv_sum / cv_count;
      }
      q.local_cv_mean = local_cv_mean;

      // ----- スコア合成 -----
      const double w = rs.quality_resolution_weight;
      // double score = q.mean_ir + 1. / (1. + local_cv_mean) - w * (fs_target_dist + 0.5 * theta_target_dist);
      // double score = q.mean_ir - w * (fs_target_dist + 0.5 * theta_target_dist);
      double score = q.worst_ir - w * (q.area_deviation + q.equilateral_deviation);
      q.score = score;
      return q;
    };

    // computeHdLimit はパッチ版を runTrials 直前に定義（重複を避けるためここからは削除）

    // --- TrialResult: 1シナリオの試行結果（best パッチを保持） ---
    struct TrialResult {
      std::string ops;
      std::string reject_detail;
      double score = -1.;
      double score_before = -1.;
      double ir_mean = 0.0;
      double ir_mean_before = 0.0;
      double ir_worst = 0.0;
      double ir_worst_before = 0.0;
      double area_deviation = 0.0;
      double area_deviation_before = 0.0;
      double equilateral_deviation = 0.0;
      double equilateral_deviation_before = 0.0;
      double hd = std::numeric_limits<double>::infinity();
      double hd_limit = std::numeric_limits<double>::infinity();
      double hd_ref_limit = std::numeric_limits<double>::infinity();
      double hd_body_limit = std::numeric_limits<double>::infinity();
      double hd_patch_diag = 0.0;
      double hd_curvature_radius = std::numeric_limits<double>::infinity();
      double hd_kmax = 0.0;
      double hd_ref_diag_limit = std::numeric_limits<double>::infinity();
      double hd_ref_curv_limit = std::numeric_limits<double>::infinity();
      double hd_body_diag_limit = std::numeric_limits<double>::infinity();
      double hd_body_curv_limit = std::numeric_limits<double>::infinity();
      bool valid = false;
      int reject_code = kRejectNoValid; // valid=false の時の reject reason
      VirtualRepairStats repair_stats;
      long long hd_checks = 0;
      long long hd_early_exits = 0;
      long long hd_samples = 0;
      long long hd_nearest_calls = 0;
      long long hd_target_missing_count = 0;
      double hd_ref = 0.0;
      double hd_body = 0.0;
      int bcinterface_no_contact = 0;
      int bcinterface_total = 0;
      bool altitude_rescue_attempted = false;
      bool altitude_rescue_succeeded = false;
      int altitude_rescue_split_worst_line = 0;
      int altitude_rescue_split_longest_face_edge = 0;
      double subsurface_altitude_rel = 1E+100;
      int subsurface_worst_side = 2;      // 0=Neumann side, 1=Dirichlet side, 2=invalid
      int subsurface_worst_line_flag = 3; // 0=Neumann line, 1=Dirichlet line, 2=BCInterface line, 3=invalid
      int subsurface_worst_face_index = -1;
      int subsurface_worst_line_index = -1;
      int subsurface_worst_subface_index = -1;
      double subsurface_worst_angle_deg = 0.0;
      bool has_subsurface_worst_tri = false;
      T3Tddd subsurface_worst_tri = T3Tddd{{Tddd{0., 0., 0.}, Tddd{0., 0., 0.}, Tddd{0., 0., 0.}}};
      std::shared_ptr<Network> patch;        // best AFTER パッチ
      std::shared_ptr<Network> source_patch; // BEFORE 参照パッチ (patch0)
    };

    auto refreshTrialPatchGeometry = [](Network& net, bool update_curvature = true) {
      for (const auto& l : net.getLines())
        l->setBoundsSingle();
      for (const auto& p : net.getPoints())
        p->setFaces();
      for (const auto& f : net.getFaces())
        f->setGeometricProperties(ToX(f->getPoints()));
      if (update_curvature)
        net.computePrincipalCurvatures(false);
    };

    enum class TrialSmoothKernel { Settings,
                                   CircumradiusToInradius,
                                   AreaLaplacian };

    struct TrialSmoothVariant {
      TrialSmoothKernel kernel = TrialSmoothKernel::Settings;
      bool sharp_safe = false;
    };

    auto doTrialSmooth = [&](Network& net,
                             TrialSmoothVariant variant = {}) {
      const double cos_thr = rs.quality_normal_flip_cos;
      const double min_ang_rad = rs.quality_min_angle_deg * M_PI / 180.0;
      const SizingField trial_sf = make_uniform_field(free_surface_target_len);
      auto trial_smooth_mode = [&]() {
        using SmoothMode = SimulationSettings::RemeshingSettings::SmoothMode;
        switch (variant.kernel) {
        case TrialSmoothKernel::CircumradiusToInradius:
          return SmoothMode::CircumradiusToInradius;
        case TrialSmoothKernel::AreaLaplacian:
          return SmoothMode::AreaLaplacian;
        case TrialSmoothKernel::Settings:
          return rs.smooth_mode;
        }
        return rs.smooth_mode;
      }();
      const bool use_projected_trial_smooth =
          trial_geometry_inputs && trial_geometry_inputs->enabled();
      std::optional<BEMMeshPipeline::BoundaryGeometryTargetProvider> target_provider;
      if (use_projected_trial_smooth)
        target_provider.emplace(makeRemeshTargetProvider(*trial_geometry_inputs));
      // Trial-only diagnostic: after virtual repair, try one smoothing block on
      // the candidate copy to see whether quality improvement must be interleaved.
      for (int si = 0; si < 20; ++si) {
        const auto border_vec = net.getBorderPoints();
        std::unordered_set<const networkPoint*> border_points(border_vec.begin(), border_vec.end());
        auto is_feature_line = [&](const networkLine* l) {
          return l && (l->BCInterface || l->SharpQ(feature_angle));
        };
        auto touches_feature_line = [&](const networkPoint* p) {
          if (!p)
            return false;
          if (p->BCInterface)
            return true;
          for (auto* l : p->getBoundaryLines())
            if (is_feature_line(l))
              return true;
          return false;
        };
        auto incident_feature_count = [&](const networkPoint* p) {
          int count = 0;
          if (!p)
            return count;
          std::unordered_set<const networkLine*> seen;
          for (auto* l : p->getBoundaryLines()) {
            if (!l || !seen.insert(l).second)
              continue;
            if (is_feature_line(l))
              ++count;
          }
          for (auto* nb : p->getNeighbors()) {
            if (!nb)
              continue;
            for (auto* l : nb->getBoundaryLines()) {
              if (!l || !seen.insert(l).second)
                continue;
              if (is_feature_line(l))
                ++count;
            }
          }
          return count;
        };
        auto feature_proximity_weight = [&](const networkPoint* p) {
          const int n = incident_feature_count(p);
          if (n <= 0)
            return 1.0;
          return std::min(2.5, 1.0 + 0.25 * static_cast<double>(n));
        };
        auto allows_projected_trial_smoothing = [&](const networkPoint* p) {
          if (!p || border_points.contains(p) || p->BorderQ())
            return false;
          return hasAnyDirichletBoundaryState(p) ||
                 hasAnyNeumannBoundaryState(p) ||
                 BEMMeshPipeline::GeometryProjectorDetail::isBCInterfaceEntity(p);
        };
        auto allows_unprojected_trial_smoothing = [&](const networkPoint* p) {
          if (!p || border_points.contains(p) || p->BorderQ())
            return false;
          if (BEMMeshPipeline::GeometryProjectorDetail::isBCInterfaceEntity(p))
            return true;
          return allows_free_surface_smoothing(p);
        };
        struct SmoothUpdate {
          networkPoint* p = nullptr;
          Tddd x = {0., 0., 0.};
          ProviderFieldAssignment fields;
        };
        std::vector<SmoothUpdate> updates;
        for (auto* p : net.getPoints()) {
          if (variant.sharp_safe && touches_feature_line(p))
            continue;
          const bool allowed = use_projected_trial_smooth
                                   ? allows_projected_trial_smoothing(p)
                                   : allows_unprojected_trial_smoothing(p);
          if (!allowed)
            continue;
          Tddd V;
          if (!remesh_smooth_delta(p, trial_smooth_mode, trial_sf, feature_angle, V))
            continue;
          if (use_projected_trial_smooth)
            V *= feature_proximity_weight(p);
          double vn = Norm(V);
          if (vn > 1e-12) {
            double limit = remesh_smooth_step_limit_ratio(trial_smooth_mode) * localEdgeLength(p);
            Tddd V_capped = std::min(vn, limit) * Normalize(V);
            const double alpha = find_safe_smooth_step(p, V_capped, cos_thr, min_ang_rad);
            if (alpha > 0.0) {
              Tddd next_x = p->X + alpha * V_capped;
              ProviderFieldAssignment fields;
              if (use_projected_trial_smooth) {
                const auto target = target_provider->queryEntityTarget(p, &net, next_x);
                if (!repairTargetAccepted(target) || !isFinite(target.target_X_clamped))
                  continue;
                if (!sampleProviderReferenceFields(*target_provider, target, p, fields))
                  continue;
                next_x = target.target_X_clamped;
              }
              updates.push_back({p, next_x, fields});
            }
          }
        }
        for (auto& update : updates) {
          if (!update.p)
            continue;
          update.p->setXSingle(update.x);
          applyProviderFieldAssignment(update.p, update.fields);
        }
        refreshTrialPatchGeometry(net, false);
        if (use_projected_trial_smooth && !updates.empty())
          refresh_trial_patch_boundary_state(net);
      }
    };

    auto borderPreserved = [](Network& net) {
      auto bset = net.getBorderPoints();
      for (auto* bp : bset) {
        bool found = false;
        for (auto* orig_p : net.copied_points) {
          if (orig_p && Norm(bp->X - orig_p->X) < 1e-12) {
            found = true;
            break;
          }
        }
        if (!found)
          return false;
      }
      return true;
    };

    auto isOuterRimLine = [](const std::unordered_set<const networkPoint*>& border_set,
                             const networkLine* l) {
      if (!l)
        return true;
      auto [p0, p1] = l->getPoints();
      return p0 && p1 && border_set.contains(p0) && border_set.contains(p1);
    };

    auto longestNonRimEdge = [&](const std::unordered_set<const networkPoint*>& border_set,
                                 networkFace* face) -> networkLine* {
      if (!face)
        return nullptr;
      networkLine* best = nullptr;
      double best_len = -1.0;
      for (auto* l : face->getLines()) {
        if (!l || isOuterRimLine(border_set, l))
          continue;
        const double len = l->length();
        if (std::isfinite(len) && len > best_len) {
          best = l;
          best_len = len;
        }
      }
      return best;
    };

    struct SplitPointCandidate {
      bool ok = false;
      Tddd x = {0., 0., 0.};
      bool projected = false;
      std::string reject_detail;
    };

    using CollapsePointCandidate = SplitPointCandidate;

    auto providerLineSamplePoint = [&](Network& net,
                                       networkLine* l,
                                       const Tddd& x_linear) -> SplitPointCandidate {
      SplitPointCandidate out;
      if (!l) {
        out.reject_detail = "line_null";
        return out;
      }
      if (!isFinite(x_linear)) {
        out.reject_detail = "midpoint_not_finite";
        return out;
      }
      out.ok = true;
      out.x = x_linear;

      if (!trial_geometry_inputs || !trial_geometry_inputs->enabled())
        return out;

      const bool repairable =
          hasAnyDirichletBoundaryState(l) ||
          hasAnyNeumannBoundaryState(l) ||
          BEMMeshPipeline::GeometryProjectorDetail::isBCInterfaceEntity(l);
      if (!repairable)
        return out;

      const auto target_provider = makeRemeshTargetProvider(*trial_geometry_inputs);
      const auto target = target_provider.queryLineSampleTarget(l, &net, 0.5, x_linear);
      if (!finalLineConstructionTargetAccepted(target, x_linear)) {
        out.ok = false;
        std::ostringstream oss;
        oss << "line_sample_target_rejected"
            << ":status=" << BEMMeshPipeline::toString(target.status)
            << ":fallback=" << BEMMeshPipeline::toString(target.fallback)
            << ":dir_gap=" << target.dirichlet_gap
            << ":body_gap=" << target.body_gap
            << ":body_gap_threshold=" << target.body_gap_threshold
            << ":move_ratio=" << target.move_ratio
            << ":shift_limit=" << target.shift_limit;
        out.reject_detail = oss.str();
        return out;
      }
      out.x = target.target_X;
      out.projected = true;
      return out;
    };

    auto splitPointForLine = [&](Network& net, networkLine* l) -> SplitPointCandidate {
      SplitPointCandidate out;
      if (!l)
        return out;
      auto [p0, p1] = l->getPoints();
      if (!p0 || !p1) {
        out.reject_detail = "line_endpoint_missing";
        return out;
      }
      const Tddd x_linear = 0.5 * (p0->X + p1->X);
      if (!isFinite(x_linear)) {
        out.reject_detail = "line_midpoint_not_finite";
        return out;
      }
      return providerLineSamplePoint(net, l, x_linear);
    };

    auto collapsePointForLine = [&](Network& net,
                                    networkLine* l,
                                    bool preserve_p0,
                                    bool preserve_p1) -> CollapsePointCandidate {
      CollapsePointCandidate out;
      if (!l)
        return out;
      auto [p0, p1] = l->getPoints();
      if (!p0 || !p1)
        return out;
      if (preserve_p0 && preserve_p1)
        return out;

      const bool preserve_bc0 = p0->BCInterface && !p1->BCInterface;
      const bool preserve_bc1 = p1->BCInterface && !p0->BCInterface;
      if ((preserve_p0 || preserve_bc0) && (preserve_p1 || preserve_bc1))
        return out;
      if (preserve_p0 || preserve_bc0) {
        if (!isFinite(p0->X))
          return out;
        out.ok = true;
        out.x = p0->X;
        return out;
      }
      if (preserve_p1 || preserve_bc1) {
        if (!isFinite(p1->X))
          return out;
        out.ok = true;
        out.x = p1->X;
        return out;
      }

      const Tddd x_linear = 0.5 * (p0->X + p1->X);
      if (!isFinite(x_linear))
        return out;
      return providerLineSamplePoint(net, l, x_linear);
    };

    enum class AltitudeRescueSplitKind { WorstLine,
                                         LongestFaceEdge };

    auto splitForAltitudeRescue = [&](Network& net,
                                      const SubsurfaceAltitudeCheckResult& issue,
                                      const AltitudeRescueSplitKind kind,
                                      int& split_worst_line,
                                      int& split_longest_face_edge) {
      const auto border_vec = net.getBorderPoints();
      const std::unordered_set<const networkPoint*> border_set(border_vec.begin(), border_vec.end());

      auto try_split = [&](networkLine* l) {
        if (!l || isOuterRimLine(border_set, l))
          return false;
        try {
          const auto split_point = splitPointForLine(net, l);
          return split_point.ok && l->Split(split_point.x, /*update_midpoints_after=*/false) != nullptr;
        } catch (...) {
          return false;
        }
      };

      if (kind == AltitudeRescueSplitKind::WorstLine && try_split(issue.worst_line)) {
        ++split_worst_line;
        return true;
      }

      if (kind == AltitudeRescueSplitKind::LongestFaceEdge &&
          try_split(longestNonRimEdge(border_set, issue.worst_face))) {
        ++split_longest_face_edge;
        return true;
      }

      return false;
    };

    // --- runPatchOps: パッチ上で ops を実行（trial 専用） ---
    auto runPatchOps = [&](const std::string& ops,
                           Network& net,
                           networkLine* cl,
                           std::string* fail_detail = nullptr) -> bool {
      auto failPatchOps = [&](std::string detail) {
        if (fail_detail)
          *fail_detail = std::move(detail);
        return false;
      };
      auto doFlip = [&](char variant, bool skip_sharp_edges, int valence_version) {
        // Patch outer rim = lines whose BOTH endpoints are patch-border
        // points. Flipping them alters valence of the un-updated mesh outside
        // the patch → skip to keep the patch-atomic contract.
        const double feat_rad = rs.feature_angle_deg * M_PI / 180.0;
        const bool use_sharp_sector_valence = (variant == 'v' && valence_version == 2);
        const double cos_thr = rs.quality_normal_flip_cos;
        const double vertex_cos_thr =
            (cos_thr > std::cos(30.0 * M_PI / 180.0))
                ? cos_thr
                : std::cos(30.0 * M_PI / 180.0);
        const double min_ang_rad = rs.quality_min_angle_deg * M_PI / 180.0;
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
            return false; // patch 外縁は除外
          if (pl->BCInterface)
            return false;
          if (skip_sharp_edges && pl->SharpQ(feat_rad))
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
            if (!valence_gain || !delaunay_violation)
              return false;
          }
          if (!flip_preserves_normals(pl, cos_thr, min_ang_rad))
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
          // Priority: valence gain DESC, tiebreak by post-flip min angle DESC.
          std::sort(cands.begin(), cands.end(),
                    [](const Cand& a, const Cand& b) {
                      if (a.vg != b.vg)
                        return a.vg > b.vg;
                      return a.minang > b.minang;
                    });
          // Apply the single best candidate per iteration (true best-first).
          // A flip changes which faces belong to each edge/node.  Refresh the
          // per-(node, face) contact states before the next flip candidate is
          // selected, otherwise later flips/smoothing can read stale BC flags.
          try {
            if (cands[0].l->Flip(true, /*update_midpoints_after=*/false))
              refresh_trial_patch_boundary_state(net);
          } catch (...) {
          }
        }
      };

      // ops 開始前の border を記憶（操作でトポロジが変わっても初期状態を参照）
      const auto initial_border = net.getBorderPoints();

      auto doSmooth = [&](TrialSmoothVariant variant) { doTrialSmooth(net, variant); };

      std::vector<SimulationSettings::RemeshingSettings::ScenarioToken> scenario_tokens;
      std::string scenario_parse_error;
      if (!SimulationSettings::RemeshingSettings::parse_scenario_tokens(
              ops, SimulationSettings::RemeshingSettings::ScenarioScope::Patch,
              &scenario_tokens, &scenario_parse_error)) {
        return failPatchOps("invalid_scenario:" + scenario_parse_error);
      }

      for (const auto& token : scenario_tokens) {
        using ScenarioOp = SimulationSettings::RemeshingSettings::ScenarioOp;
        if (token.op == ScenarioOp::Flip) {
          doFlip(token.variant, token.sharp_safe, token.valence_version);
        } else if (token.op == ScenarioOp::Collapse) {
          try {
            if (cl->getBoundaryFaces().size() != 2)
              return failPatchOps("collapse_trigger_not_two_face");
            if (cl->BCInterface)
              return failPatchOps("collapse_bcinterface");
            auto [cp0, cp1] = cl->getPoints();
            bool cp0_border = initial_border.count(cp0);
            bool cp1_border = initial_border.count(cp1);
            if (cp0_border && cp1_border)
              return failPatchOps("collapse_both_endpoints_border");
            const auto collapse_point = collapsePointForLine(net, cl, cp0_border, cp1_border);
            if (!collapse_point.ok)
              return failPatchOps(collapse_point.reject_detail.empty()
                                      ? "collapse_target_rejected"
                                      : ("collapse_target_rejected:" + collapse_point.reject_detail));
            if (!cl->Collapse(collapse_point.x, /*update_midpoints_after=*/false))
              return failPatchOps("collapse_topology_failed");
            refresh_trial_patch_boundary_state(net);
          } catch (...) {
            return failPatchOps("collapse_exception");
          }
        } else if (token.op == ScenarioOp::Smooth) {
          TrialSmoothVariant smooth_variant;
          if (token.has_smooth_mode_suffix) {
            using SmoothMode = SimulationSettings::RemeshingSettings::SmoothMode;
            if (token.smooth_mode == SmoothMode::CircumradiusToInradius)
              smooth_variant.kernel = TrialSmoothKernel::CircumradiusToInradius;
            else
              smooth_variant.kernel = TrialSmoothKernel::AreaLaplacian;
          }
          smooth_variant.sharp_safe = token.sharp_safe;
          // setDodecaPoints を除外した軽量版 setGeometricPropertiesForce
          // DodecaPoints はパッチ border 辺を通じてパッチ外面を参照するため OMP 並列で不安全
          refreshTrialPatchGeometry(net, true);
          doSmooth(smooth_variant);
          refresh_trial_patch_boundary_state(net);
          const auto smooth_repair = applyVirtualGeometryRepair(net);
          if (smooth_repair.repair_failed > 0)
            return false;
          refreshTrialPatchGeometry(net, true);
          refresh_trial_patch_boundary_state(net);
        } else if (token.op == ScenarioOp::Split) {
          const auto split_point = splitPointForLine(net, cl);
          if (!split_point.ok)
            return failPatchOps(split_point.reject_detail.empty()
                                    ? "split_target_rejected"
                                    : ("split_target_rejected:" + split_point.reject_detail));
          if (!cl->Split(split_point.x, /*update_midpoints_after=*/false))
            return failPatchOps("split_topology_failed");
          refresh_trial_patch_boundary_state(net);
        }
      }
      return true;
    };

    auto measureWaterlineFaceQuality = [](const V_netFp& faces) {
      BEMPreBVP::Stats stats;
      for (const auto* f : faces)
        BEMPreBVP::accumulate_face_quality(stats, f);
      return stats;
    };

    auto waterlineFaceQualityFailsLimit = [&](const BEMPreBVP::Stats& stats) {
      if (!rs.waterline_protection_enabled || stats.waterline_faces <= 0)
        return false;
      const bool bad_min_angle =
          !std::isfinite(stats.waterline_face_min_angle_deg) ||
          stats.waterline_face_min_angle_deg < rs.waterline_face_min_angle_deg;
      const bool bad_max_aspect =
          !std::isfinite(stats.waterline_face_aspect_max) ||
          stats.waterline_face_aspect_max > rs.waterline_face_max_aspect;
      return bad_min_angle || bad_max_aspect;
    };

    // --- runTrial: 事前コピー済みパッチ上で操作+評価（並列実行される） ---
    struct HdLimits {
      double ref = std::numeric_limits<double>::infinity();
      double body = std::numeric_limits<double>::infinity();
      double diag = 0.0;
      double curvature_radius = std::numeric_limits<double>::infinity();
      double kmax = 0.0;
      double ref_diag = std::numeric_limits<double>::infinity();
      double ref_curv = std::numeric_limits<double>::infinity();
      double body_diag = std::numeric_limits<double>::infinity();
      double body_curv = std::numeric_limits<double>::infinity();
    };

    auto runTrial = [&](const std::string& ops, std::shared_ptr<Network> patch, networkLine* cl,
                        const V_netFp& ref_faces,
                        const HdLimits hd_limits) -> TrialResult {
      TrialResult result;
      result.ops = ops;
      result.hd_ref_limit = hd_limits.ref;
      result.hd_body_limit = hd_limits.body;
      result.hd_limit = hd_limits.body; // legacy log/CSV column: previous single limit value
      result.hd_patch_diag = hd_limits.diag;
      result.hd_curvature_radius = hd_limits.curvature_radius;
      result.hd_kmax = hd_limits.kmax;
      result.hd_ref_diag_limit = hd_limits.ref_diag;
      result.hd_ref_curv_limit = hd_limits.ref_curv;
      result.hd_body_diag_limit = hd_limits.body_diag;
      result.hd_body_curv_limit = hd_limits.body_curv;
      auto rejectIfBCInterfaceNoContact = [&](const BCInterfaceContactStats& stats) {
        if (stats.no_contact() <= 0)
          return false;
        result.reject_code = kRejectBCInterfaceNoContact;
        result.bcinterface_no_contact = stats.no_contact();
        result.bcinterface_total = stats.total();
        return true;
      };
      if (!cl)
        return result;
      if (!trial_geometry_inputs || !trial_geometry_inputs->enabled()) {
        result.reject_code = kRejectRepairFailed;
        return result;
      }

      try {
        if (cl->BCInterface)
          refresh_trial_patch_boundary_state(*patch);
        std::string patch_ops_detail;
        if (!runPatchOps(ops, *patch, cl, &patch_ops_detail)) {
          result.reject_code = kRejectPatchOpsFailed;
          result.reject_detail = std::move(patch_ops_detail);
          return result;
        }
        if (rejectIfBCInterfaceNoContact(refresh_trial_patch_boundary_state(*patch)))
          return result;

        // 外周保存チェック
        if (!borderPreserved(*patch)) {
          result.reject_code = kRejectBorderMoved;
          return result;
        }

        if (rejectIfBCInterfaceNoContact(refresh_trial_patch_boundary_state(*patch)))
          return result;
        result.repair_stats = applyVirtualGeometryRepair(*patch);
        if (result.repair_stats.repair_failed > 0) {
          result.reject_code = kRejectRepairFailed;
          return result;
        }

        // setDodecaPoints を除外した軽量版（パッチ border 辺がパッチ外面を参照する問題を回避）
        refreshTrialPatchGeometry(*patch, true);
        if (rejectIfBCInterfaceNoContact(refresh_trial_patch_boundary_state(*patch)))
          return result;
        auto rejectIfWaterlineFaceQualityBad = [&](const V_netFp& faces) {
          if (!rs.waterline_protection_enabled)
            return false;
          const auto face_quality = measureWaterlineFaceQuality(faces);
          if (waterlineFaceQualityFailsLimit(face_quality)) {
            result.reject_code = kRejectWaterlineQuality;
            return true;
          }
          return false;
        };
        if (rejectIfWaterlineFaceQualityBad(patch->getBoundaryFaces()))
          return result;
        if (rs.waterline_protection_enabled) {
          V_netFp all_patch_faces;
          all_patch_faces.reserve(patch->getFaces().size());
          for (auto* f : patch->getFaces())
            all_patch_faces.emplace_back(f);
          if (rejectIfWaterlineFaceQualityBad(all_patch_faces))
            return result;
        }
        {
          const auto tiny = checkTinyFacesRelativeToLocalMean(
              *patch, /*use_quadratic_subfaces=*/false);
          if (tiny.worst_face && tiny.worst_area_ratio < min_local_face_area_ratio) {
            result.reject_code = kRejectTinyFace;
            return result;
          }
        }

        auto currentSubsurfaceIssue = [&]() {
          if (!rs.subsurface_altitude_reject.enabled)
            return SubsurfaceAltitudeCheckResult{};
          return checkSubsurfaceFaceAltitude(
              *patch, rs.subsurface_altitude_reject, /*use_quadratic_subfaces=*/false);
        };

        auto subsurfaceAltitudeRejected = [&](const SubsurfaceAltitudeCheckResult& issue) {
          if (!rs.subsurface_altitude_reject.enabled)
            return false;
          result.subsurface_altitude_rel = issue.worst_line ? issue.worst_altitude_rel : 1E+100;
          if (issue.worst_line && issue.worst_face) {
            result.subsurface_worst_side = edgeFaceSideIndex(issue.worst_line, issue.worst_face);
            result.subsurface_worst_line_flag = lineFlagIndex(issue.worst_line);
            result.subsurface_worst_face_index = issue.worst_face->index;
            result.subsurface_worst_subface_index = issue.worst_subface_index;
            result.subsurface_worst_angle_deg = issue.worst_angle_deg;
            result.subsurface_worst_tri = snapshotFaceTri(issue.worst_face);
            result.has_subsurface_worst_tri = true;
          }
          return std::isfinite(result.subsurface_altitude_rel) &&
                 result.subsurface_altitude_rel < rs.subsurface_altitude_reject.min_face_altitude_rel;
        };

        auto issue = currentSubsurfaceIssue();
        if (subsurfaceAltitudeRejected(issue)) {
          if (!trial_geometry_inputs || !trial_geometry_inputs->enabled()) {
            result.reject_code = kRejectSubsurfaceAltitude;
            return result;
          }

          result.altitude_rescue_attempted = true;
          bool rescued = false;
          auto applyRescuePostprocess = [&]() {
            refreshTrialPatchGeometry(*patch, true);
            doTrialSmooth(*patch);
            const auto rescue_repair = applyVirtualGeometryRepair(*patch);
            mergeVirtualRepairStats(result.repair_stats, rescue_repair);
            if (rescue_repair.repair_failed > 0) {
              result.reject_code = kRejectRepairFailed;
              return false;
            }
            refreshTrialPatchGeometry(*patch, true);
            if (rejectIfBCInterfaceNoContact(refresh_trial_patch_boundary_state(*patch)))
              return false;
            if (!borderPreserved(*patch)) {
              result.reject_code = kRejectBorderMoved;
              return false;
            }
            issue = currentSubsurfaceIssue();
            if (!subsurfaceAltitudeRejected(issue)) {
              rescued = true;
              return false;
            }
            return true;
          };

          for (int rescue_iter = 0; rescue_iter < 2; ++rescue_iter) {
            bool did_split = false;
            if (splitForAltitudeRescue(*patch, issue, AltitudeRescueSplitKind::WorstLine,
                                       result.altitude_rescue_split_worst_line,
                                       result.altitude_rescue_split_longest_face_edge)) {
              did_split = true;
              if (!applyRescuePostprocess())
                break;
            }
            if (rescued)
              break;
            if (splitForAltitudeRescue(*patch, issue, AltitudeRescueSplitKind::LongestFaceEdge,
                                       result.altitude_rescue_split_worst_line,
                                       result.altitude_rescue_split_longest_face_edge)) {
              did_split = true;
              if (!applyRescuePostprocess())
                break;
            }
            if (!did_split || rescued)
              break;
          }

          result.altitude_rescue_succeeded = rescued;
          if (result.reject_code == kRejectRepairFailed || result.reject_code == kRejectBorderMoved)
            return result;
          if (!rescued) {
            result.reject_code = kRejectSubsurfaceAltitude;
            return result;
          }
        }

        V_netFp modified_faces = patch->getBoundaryFaces();
        result.hd = 0.;
        auto rejectIfHdExceeded = [&]() {
          result.reject_code = kRejectHdExceeded;
          ++result.hd_early_exits;
          result.patch = patch; // debug output: keep failed AFTER patch
          return result;
        };
        auto accountTargetSurfaceHd = [&](double distance, bool body_target) {
          if (!std::isfinite(distance))
            return;
          result.hd = std::max(result.hd, distance);
          if (body_target)
            result.hd_body = std::max(result.hd_body, distance);
          else
            result.hd_ref = std::max(result.hd_ref, distance);
        };
        auto targetHdExceeded = [&]() {
          return result.hd_target_missing_count > 0 ||
                 (std::isfinite(hd_limits.ref) && result.hd_ref > hd_limits.ref) ||
                 (std::isfinite(hd_limits.body) && result.hd_body > hd_limits.body);
        };

        const bool target_aware_hd =
            trial_geometry_inputs && trial_geometry_inputs->enabled();
        if (!target_aware_hd) {
          result.reject_code = kRejectRepairFailed;
          return result;
        }
        {
          const auto hd_target_provider = makeRemeshTargetProvider(*trial_geometry_inputs);
          auto accountReferenceTarget = [&](const Tddd& x) {
            ++result.hd_samples;
            const auto projection = hd_target_provider.queryReferenceDistance(&water, x);
            if (!projection.ok || !std::isfinite(projection.distance)) {
              ++result.hd_target_missing_count;
              return false;
            }
            accountTargetSurfaceHd(projection.distance, /*body_target=*/false);
            return !targetHdExceeded();
          };
          auto accountFaceBodyTarget = [&](networkFace* f, const Tddd& x) {
            ++result.hd_samples;
            BEMMeshPipeline::GeometryProjectorDetail::BodyProjection best_body;
            bool found = false;
            auto visit_entity = [&](auto* entity) {
              if (!entity)
                return;
              if (getNodeFaceBoundaryType(entity, f) != NodeFaceBoundaryType::Neumann)
                return;
              const auto body = hd_target_provider.queryBodyDistanceForEntity(entity, x);
              ++result.hd_nearest_calls;
              if (!body.ok || !std::isfinite(body.distance))
                return;
              if (!found || body.distance < best_body.distance) {
                best_body = body;
                found = true;
              }
            };
            if (f) {
              auto [p0, p1, p2] = f->getPoints();
              auto [l0, l1, l2] = f->getLines();
              visit_entity(p0);
              visit_entity(p1);
              visit_entity(p2);
              visit_entity(l0);
              visit_entity(l1);
              visit_entity(l2);
            }
            if (!found) {
              ++result.hd_target_missing_count;
              return false;
            }
            accountTargetSurfaceHd(best_body.distance, /*body_target=*/true);
            return !targetHdExceeded();
          };
          auto accountEntityTarget = [&](auto* entity, const Tddd& x) {
            ++result.hd_checks;
            const auto target = hd_target_provider.queryEntityTarget(entity, &water, x);
            if (!repairTargetAccepted(target)) {
              ++result.hd_target_missing_count;
              return false;
            }
            ++result.hd_samples;
            const auto flags = hd_target_provider.targetFlags(entity);
            if (flags.need_reference)
              accountTargetSurfaceHd(target.dirichlet_gap, /*body_target=*/false);
            if (flags.need_body)
              accountTargetSurfaceHd(target.body_gap, /*body_target=*/true);
            return !targetHdExceeded();
          };
          auto evalBoundaryFaceQuad = [](networkFace* f, double t0, double t1) -> Tddd {
            auto [p0, p1, p2] = f->getPoints();
            auto N = TriShape<6>(t0, t1);
            const Tddd m0 = 0.5 * (p0->X + p1->X);
            const Tddd m1 = 0.5 * (p1->X + p2->X);
            const Tddd m2 = 0.5 * (p2->X + p0->X);
            return N[0] * p0->X + N[1] * p1->X + N[2] * p2->X +
                   N[3] * m0 + N[4] * m1 + N[5] * m2;
          };
          auto runTargetAwareFaceHd = [&]() {
            ++result.hd_checks;
            constexpr auto t0t1_samples = UniformPointsOnTriangle_00_10_01<3>();
            for (auto* f : modified_faces) {
              if (!f)
                continue;
              const auto flags = hd_target_provider.targetFlags(f);
              const bool need_ref = flags.need_reference;
              const bool need_body = flags.need_body;
              if (!need_ref && !need_body)
                continue;
              for (const auto& [t0, t1] : t0t1_samples) {
                const auto x = evalBoundaryFaceQuad(f, t0, t1);
                if (need_ref && !accountReferenceTarget(x))
                  return false;
                if (need_body && !accountFaceBodyTarget(f, x))
                  return false;
              }
            }
            return true;
          };
          auto runTargetAwareMidpointHd = [&]() {
            std::unordered_set<networkLine*> lines;
            for (auto* f : modified_faces)
              if (f)
                for (auto* l : f->getLines())
                  lines.insert(l);
            for (auto* l : lines) {
              if (!l || l->getBoundaryFaces().size() < 2)
                continue;
              auto [p0, p1] = l->getPoints();
              if (!p0 || !p1)
                continue;
              const Tddd x_linear = 0.5 * (p0->X + p1->X);
              if (!isFinite(x_linear))
                continue;
              if (!accountEntityTarget(l, x_linear))
                return false;
            }
            return true;
          };
          if (!runTargetAwareFaceHd() || !runTargetAwareMidpointHd())
            return rejectIfHdExceeded();
        }

        const auto quality_after = patchQuality(*patch);
        result.score = quality_after.score;
        result.ir_mean = quality_after.mean_ir;
        result.ir_worst = quality_after.worst_ir;
        result.area_deviation = quality_after.area_deviation;
        result.equilateral_deviation = quality_after.equilateral_deviation;
        result.valid = true;
        result.patch = patch;
      } catch (...) {
        result.reject_code = kRejectException;
        // OMP 並列内で例外を伝播させない
      }
      return result;
    };

    // --- computeHdLimit: パッチ上で許容 Hausdorff 距離を計算 ---
    // 形状忠実度の上限を決める。reference/body の意味別に別々の上限を持つ。
    // 各上限は 2 つの制約の min を取る:
    //   (a) bbox 対角 × rs.quality_hd_diag_ratio  — パッチ全体のスケール制限
    //   (b) 曲率半径 R × rs.quality_hd_curv_ratio — 角の近傍ではより厳しく
    // 呼び出し前提: patch の geometric properties と curvature が最新であること
    auto computeHdLimits = [&rs](const Network& patch) -> HdLimits {
      auto bnds = patch.getBounds();
      double dx = std::get<1>(std::get<0>(bnds)) - std::get<0>(std::get<0>(bnds));
      double dy = std::get<1>(std::get<1>(bnds)) - std::get<0>(std::get<1>(bnds));
      double dz = std::get<1>(std::get<2>(bnds)) - std::get<0>(std::get<2>(bnds));
      double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
      // 曲率半径制限: パッチ内最大 kmax から R = 1/kmax_max
      double kmax_max = 0.;
      for (auto* p : patch.getBoundaryPoints())
        if (p->geom_curvature.valid)
          kmax_max = std::max(kmax_max, p->geom_curvature.kmax);
      double R = (kmax_max > 1e-10) ? (1.0 / kmax_max) : std::numeric_limits<double>::infinity();
      HdLimits out;
      out.diag = diag;
      out.curvature_radius = R;
      out.kmax = kmax_max;
      out.ref_diag = 0.5 * rs.quality_hd_diag_ratio * diag;
      out.ref_curv = 0.5 * rs.quality_hd_curv_ratio * R;
      out.body_diag = rs.quality_hd_diag_ratio * diag;
      out.body_curv = rs.quality_hd_curv_ratio * R;
      out.ref = std::min(out.ref_diag, out.ref_curv);
      out.body = std::min(out.body_diag, out.body_curv);
      return out;
    };

    // --- runTrials: copy は逐次、操作+評価は並列 ---
    // 採用条件は常に「trial patch の score が before + margin を超えること」。
    // 以前は aggressive split で score 改善なしの valid trial も採用していたが、
    // 局所的に悪い修正が通る原因になるため、split/collapse/flip で統一する。
    auto runTrials = [&](networkLine* target_line,
                         const std::vector<std::string>& scenarios,
                         bool aggressive = false,
                         int remesh_iter = -1,
                         int op_code = -1,
                         int reason_code = -1) -> TrialResult {
      ++rt_calls;
      if (training_data_record_enabled)
        ++training_iter_stats.candidate_edges;
      auto accumulate_rt_repair_stats = [&](const VirtualRepairStats& s) {
        rt_repair_failed += s.repair_failed;
        rt_repaired_points += s.repaired_points;
        rt_repaired_lines += s.repaired_lines;
        rt_max_free_gap = std::max(rt_max_free_gap, s.max_free_gap);
        rt_max_body_gap = std::max(rt_max_body_gap, s.max_body_gap);
        rt_max_move_ratio = std::max(rt_max_move_ratio, s.max_move_ratio);
      };
      auto t_ref_begin = clk::now();
      auto patch0 = std::make_shared<Network>("file_name_is_not_given", "ref");
      water.copyLocalPatch(*patch0, target_line, 2);
      refresh_trial_patch_boundary_state(*patch0);
      const auto source_repair_stats = applyVirtualGeometryRepair(*patch0);
      refreshTrialPatchGeometry(*patch0, true);
      refresh_trial_patch_boundary_state(*patch0);
      if (trial_geometry_inputs && trial_geometry_inputs->enabled()) {
        ++rt_repair_aware_calls;
        accumulate_rt_repair_stats(source_repair_stats);
      }
      patch0->computePrincipalCurvatures();
      V_netFp ref_faces = patch0->getBoundaryFaces();
      const auto quality_before = patchQuality(*patch0);
      double score_before = quality_before.score;
      const double source_altitude_rel = worstSubsurfaceAltitudeRel(*patch0);

      // copyLocalPatch 内部で既に setGeometricPropertiesForce 済み (Network.cpp:2889)、
      // computePrincipalCurvatures は geom_curvature のみ書くので face props は最新。
      const HdLimits hd_limits = computeHdLimits(*patch0);
      t_rt_ref += sec_since(t_ref_begin);

      // 各シナリオ用のパッチをコピー
      // copyLocalPatch は OMP-unsafe (race condition で malloc double-free クラッシュ)。
      // scenario_engine.cpp と同じ方針で **serial にする**。run_trial のみ OMP で回す。
      // 参考: bem/time_domain/BEM_remesh_scenario_engine.cpp の runTrials 内コメント。
      auto t_copy_begin = clk::now();
      std::vector<std::shared_ptr<Network>> patches(scenarios.size());
      std::vector<networkLine*> copied_lines(scenarios.size(), nullptr);
      for (std::size_t i = 0; i < scenarios.size(); ++i) {
        patches[i] = std::make_shared<Network>("file_name_is_not_given", scenarios[i]);
        copied_lines[i] = water.copyLocalPatch(*patches[i], target_line, 2);
      }
      t_rt_copy += sec_since(t_copy_begin);

      // 並列: 操作 + 評価
      auto t_run_begin = clk::now();
      std::vector<TrialResult> results(scenarios.size());
      _Pragma("omp parallel for") for (std::size_t i = 0; i < scenarios.size(); ++i)
          results[i] = runTrial(scenarios[i], patches[i], copied_lines[i], ref_faces, hd_limits);
      t_rt_run += sec_since(t_run_begin);
      for (auto& r : results) {
        r.score_before = score_before;
        r.ir_mean_before = quality_before.mean_ir;
        r.ir_worst_before = quality_before.worst_ir;
        r.area_deviation_before = quality_before.area_deviation;
        r.equilateral_deviation_before = quality_before.equilateral_deviation;
        r.hd_ref_limit = hd_limits.ref;
        r.hd_body_limit = hd_limits.body;
        r.hd_limit = hd_limits.body;
        r.hd_patch_diag = hd_limits.diag;
        r.hd_curvature_radius = hd_limits.curvature_radius;
        r.hd_kmax = hd_limits.kmax;
        r.hd_ref_diag_limit = hd_limits.ref_diag;
        r.hd_ref_curv_limit = hd_limits.ref_curv;
        r.hd_body_diag_limit = hd_limits.body_diag;
        r.hd_body_curv_limit = hd_limits.body_curv;
      }
      auto hasUnacceptableIrRegression = [](const TrialResult& r) {
        constexpr double kWorstIrRegressionRatio = 0.95;
        constexpr double kWorstIrRegressionAbsTol = 1e-4;
        if (!r.valid)
          return false;
        if (!std::isfinite(r.ir_worst) || !std::isfinite(r.ir_worst_before))
          return true;
        if (r.ir_worst_before <= 0.0)
          return false;
        return r.ir_worst + kWorstIrRegressionAbsTol <
               kWorstIrRegressionRatio * r.ir_worst_before;
      };
      for (auto& r : results) {
        if (!hasUnacceptableIrRegression(r))
          continue;
        r.valid = false;
        r.reject_code = kRejectIrRegression;
        r.reject_detail = "ir_worst_regression";
      }

      if (target_line && target_line->BCInterface) {
        std::unordered_map<std::string, int> patch_ops_reject_details;
        int patch_ops_rejects = 0;
        for (const auto& r : results) {
          if (r.reject_code != kRejectPatchOpsFailed)
            continue;
          ++patch_ops_rejects;
          const std::string key = r.reject_detail.empty() ? "unspecified" : r.reject_detail;
          ++patch_ops_reject_details[key];
        }
        if (patch_ops_rejects > 0) {
          std::ostringstream oss;
          bool first = true;
          for (const auto& [detail, count] : patch_ops_reject_details) {
            if (!first)
              oss << ',';
            first = false;
            oss << detail << ':' << count;
          }
          std::cout << Yellow
                    << "[remesh:bci_ops_reject] step=" << time_step
                    << " phase=" << sanitize_phase_tag(phase_debug_tag)
                    << " op_code=" << op_code
                    << " reason_code=" << reason_code
                    << " p0=" << endpointIndex0(target_line)
                    << " p1=" << endpointIndex1(target_line)
                    << " len=" << target_line->length()
                    << " target_len=" << targetLengthFor(target_line)
                    << " scenarios=" << scenarios.size()
                    << " patch_ops_rejects=" << patch_ops_rejects
                    << " details={" << oss.str() << "}"
                    << colorReset << std::endl;
        }
      }

      if (training_data_record_enabled) {
        training_iter_stats.trial_count += static_cast<int>(results.size());
        const int category = edgeTrialCategory(target_line);
        const int p0_index = endpointIndex0(target_line);
        const int p1_index = endpointIndex1(target_line);
        const double len = target_line ? target_line->length() : -1.0;
        const double target_len = targetLengthFor(target_line);
        const double len_ratio = (target_len > 1e-20) ? len / target_len : -1.0;
        const double theta = edgeThetaValue(target_line);
        const double theta_ratio = (theta >= 0.0 && theta_split > 1e-20) ? theta / theta_split : -1.0;
        for (std::size_t si = 0; si < results.size(); ++si) {
          auto& r = results[si];
          const double gain = r.score - r.score_before;
          const bool safe = r.valid && gain > rs.training_data_min_score_gain;
          if (safe)
            ++training_iter_stats.safe_count;
          const bool save_row =
              safe ? !int_limit_reached(training_iter_stats.saved_success,
                                        rs.training_data_max_saved_success_per_iter)
                   : !int_limit_reached(training_iter_stats.saved_failure,
                                        rs.training_data_max_saved_failure_per_iter);
          if (!save_row)
            continue;
          if (safe)
            ++training_iter_stats.saved_success;
          else
            ++training_iter_stats.saved_failure;
          std::ostringstream row;
          row << time_step << ',' << remesh_iter << ','
              << op_code << ',' << reason_code << ',' << category << ','
              << p0_index << ',' << p1_index << ','
              << len << ',' << len_ratio << ',' << theta << ',' << theta_ratio << ','
              << csvString(si < scenarios.size() ? scenarios[si] : r.ops) << ','
              << (r.valid ? 1 : 0) << ',' << r.reject_code << ','
              << csvString(rejectCodeName(r.reject_code)) << ','
              << r.score_before << ',' << r.score << ',' << gain << ','
              << r.ir_mean_before << ',' << r.ir_mean << ','
              << (r.ir_mean - r.ir_mean_before) << ','
              << r.ir_worst_before << ',' << r.ir_worst << ','
              << r.area_deviation_before << ',' << r.area_deviation << ','
              << (r.area_deviation - r.area_deviation_before) << ','
              << r.equilateral_deviation_before << ',' << r.equilateral_deviation << ','
              << (r.equilateral_deviation - r.equilateral_deviation_before) << ','
              << r.hd << ',' << r.hd_ref << ',' << r.hd_body << ','
              << r.hd_limit << ',' << r.bcinterface_no_contact << ','
              << r.subsurface_altitude_rel << ','
              << (safe ? "success" : "failure");
          appendWithHeader(training_trial_csv, training_trial_header, row.str());
        }
      }

      for (const auto& r : results) {
        rt_hd_checks += r.hd_checks;
        rt_hd_early_exits += r.hd_early_exits;
        rt_hd_samples += r.hd_samples;
        rt_hd_nearest_calls += r.hd_nearest_calls;
        rt_hd_target_missing += r.hd_target_missing_count;
        rt_hd_ref_max = std::max(rt_hd_ref_max, r.hd_ref);
        rt_hd_body_max = std::max(rt_hd_body_max, r.hd_body);
      }

      if (trial_geometry_inputs && trial_geometry_inputs->enabled())
        for (const auto& r : results) {
          accumulate_rt_repair_stats(r.repair_stats);
          if (r.altitude_rescue_attempted)
            ++rt_altitude_rescue_attempts;
          if (r.altitude_rescue_succeeded)
            ++rt_altitude_rescue_successes;
          rt_altitude_rescue_split_worst_line += r.altitude_rescue_split_worst_line;
          rt_altitude_rescue_split_longest_face_edge += r.altitude_rescue_split_longest_face_edge;
          const int trigger_bc = lineFlagIndex(target_line);
          if (trigger_bc >= 0 && trigger_bc < 3 &&
              r.reject_code == kRejectSubsurfaceAltitude) {
            if (r.subsurface_worst_side >= 0 && r.subsurface_worst_side < 3)
              ++rt_subsurface_worst_side_by_trigger[trigger_bc][r.subsurface_worst_side];
            if (r.subsurface_worst_line_flag >= 0 && r.subsurface_worst_line_flag < 4)
              ++rt_subsurface_worst_line_flag_by_trigger[trigger_bc][r.subsurface_worst_line_flag];
            ++rt_subsurface_worst_samples[trigger_bc];
            if (r.has_subsurface_worst_tri)
              remesh_debug.recordSubsurfaceRejectFace(
                  r.subsurface_worst_tri,
                  trigger_bc,
                  r.subsurface_worst_side,
                  r.subsurface_worst_line_flag,
                  r.subsurface_altitude_rel,
                  r.subsurface_worst_angle_deg);
          }
        }

      // スコア表出力（デバッグ時のみ有効化）
      // {
      //   std::cout << "  [trials] before=" << std::setprecision(3) << score_before
      //             << " hd_lim=" << hd_limit << std::endl;
      //   for (std::size_t i = 0; i < results.size(); ++i) {
      //     auto& r = results[i];
      //     std::cout << "    " << std::setw(8) << std::left << scenarios[i]
      //               << " valid=" << r.valid
      //               << " score=" << std::setprecision(3) << std::setw(7) << r.score
      //               << " hd=" << std::setprecision(3) << std::setw(9) << r.hd
      //               << (r.valid && r.hd <= hd_limit && r.score > score_before && r.score >= 0.05 ? " *" : "")
      //               << std::endl;
      //   }
      // }

      auto t_pick_begin = clk::now();
      TrialResult best;
      (void)aggressive;
      double score_margin = rs.quality_score_improve_margin;
      const double best_score_threshold = score_before + score_margin;
      double best_score = best_score_threshold;
      int n_valid = 0, n_within_hd = 0, n_improved = 0;
      double best_subsurface_rejected_altitude = -std::numeric_limits<double>::infinity();
      std::array<int, kRejectCodeCount> invalid_counts{};
      TrialResult best_hd_reject;
      bool has_best_hd_reject = false;
      TrialResult best_score_reject;
      bool has_best_score_reject = false;
      for (auto& r : results) {
        if (!r.valid) {
          if (r.reject_code >= 0 && r.reject_code < kRejectCodeCount)
            ++invalid_counts[r.reject_code];
          if (r.reject_code == kRejectHdExceeded && r.patch && !has_best_hd_reject) {
            best_hd_reject = r;
            has_best_hd_reject = true;
          }
          if (r.reject_code == kRejectSubsurfaceAltitude &&
              std::isfinite(r.subsurface_altitude_rel))
            best_subsurface_rejected_altitude =
                std::max(best_subsurface_rejected_altitude, r.subsurface_altitude_rel);
          continue;
        }
        ++n_valid;
        if ((std::isfinite(hd_limits.ref) && r.hd_ref > hd_limits.ref) ||
            (std::isfinite(hd_limits.body) && r.hd_body > hd_limits.body))
          continue;
        ++n_within_hd;
        if (!has_best_score_reject || r.score > best_score_reject.score) {
          best_score_reject = r;
          has_best_score_reject = true;
        }
        if (r.score > best_score) {
          ++n_improved;
          best_score = r.score;
          best = std::move(r);
        }
      }
      // 失敗原因の分類 (valid=true になった場合は 0 に上書き)
      if (n_valid == 0) {
        if (invalid_counts[kRejectHdExceeded] > 0) {
          if (has_best_hd_reject)
            best = best_hd_reject;
          best.reject_code = kRejectHdExceeded;
        } else {
          int dominant_code = kRejectNoValid;
          int dominant_count = 0;
          for (int rc = 0; rc < kRejectCodeCount; ++rc) {
            if (invalid_counts[rc] > dominant_count) {
              dominant_count = invalid_counts[rc];
              dominant_code = rc;
            }
          }
          best.reject_code = dominant_count > 0 ? dominant_code : kRejectNoValid;
        }
      } else if (n_within_hd == 0)
        best.reject_code = kRejectHdExceeded;
      else if (n_improved == 0) {
        if (has_best_score_reject)
          best = best_score_reject;
        best.valid = false;
        best.reject_code = kRejectScoreNoImprove;
      }
      if (best.valid)
        best.reject_code = kRejectSuccess; // 後で replacePatch 失敗なら 4 に上書きする
      if (std::isfinite(best_subsurface_rejected_altitude)) {
        ++rt_best_rejected_altitude_cases;
        if (std::isfinite(source_altitude_rel))
          rt_best_rejected_altitude_before =
              std::max(rt_best_rejected_altitude_before, source_altitude_rel);
        rt_best_rejected_altitude_after =
            std::max(rt_best_rejected_altitude_after, best_subsurface_rejected_altitude);
        if (std::isfinite(source_altitude_rel))
          rt_best_rejected_altitude_gain =
              std::max(rt_best_rejected_altitude_gain,
                       best_subsurface_rejected_altitude - source_altitude_rel);
      }
      if (target_line && target_line->BCInterface) {
        std::array<int, kRejectCodeCount> scenario_reject_counts{};
        for (const auto& r : results) {
          const int code = r.valid ? kRejectSuccess : r.reject_code;
          if (code >= 0 && code < kRejectCodeCount)
            ++scenario_reject_counts[code];
        }
        std::ostringstream oss;
        bool first = true;
        for (int rc = 0; rc < kRejectCodeCount; ++rc) {
          if (scenario_reject_counts[rc] <= 0)
            continue;
          if (!first)
            oss << ',';
          first = false;
          oss << rejectCodeName(rc) << ':' << scenario_reject_counts[rc];
        }
        std::cout << Yellow
                  << "[remesh:bci_trial] step=" << time_step
                  << " phase=" << sanitize_phase_tag(phase_debug_tag)
                  << " op_code=" << op_code
                  << " reason_code=" << reason_code
                  << " p0=" << endpointIndex0(target_line)
                  << " p1=" << endpointIndex1(target_line)
                  << " len=" << target_line->length()
                  << " target_len=" << targetLengthFor(target_line)
                  << " scenarios=" << scenarios.size()
                  << " valid=" << n_valid
                  << " within_hd=" << n_within_hd
                  << " improved=" << n_improved
                  << " score_before=" << score_before
                  << " score_threshold=" << best_score_threshold
                  << " best_score=" << best_score
                  << " hd_ref_limit=" << hd_limits.ref
                  << " hd_body_limit=" << hd_limits.body
                  << " hd_limit_source_ref=" << (hd_limits.ref_diag <= hd_limits.ref_curv ? "diag" : "curv")
                  << " hd_limit_source_body=" << (hd_limits.body_diag <= hd_limits.body_curv ? "diag" : "curv")
                  << " hd_patch_diag=" << hd_limits.diag
                  << " hd_R=" << hd_limits.curvature_radius
                  << " hd_kmax=" << hd_limits.kmax
                  << " hd_ref_diag_limit=" << hd_limits.ref_diag
                  << " hd_ref_curv_limit=" << hd_limits.ref_curv
                  << " hd_body_diag_limit=" << hd_limits.body_diag
                  << " hd_body_curv_limit=" << hd_limits.body_curv
                  << " hd_reject_ref=" << (has_best_hd_reject ? best_hd_reject.hd_ref : -1.0)
                  << " hd_reject_body=" << (has_best_hd_reject ? best_hd_reject.hd_body : -1.0)
                  << " hd_reject_side="
                  << (has_best_hd_reject
                          ? ((best_hd_reject.hd_target_missing_count > 0)
                                 ? "missing"
                                 : ((std::isfinite(hd_limits.ref) &&
                                     best_hd_reject.hd_ref > hd_limits.ref)
                                        ? "ref"
                                        : ((std::isfinite(hd_limits.body) &&
                                            best_hd_reject.hd_body > hd_limits.body)
                                               ? "body"
                                               : "none")))
                          : "none")
                  << " final=" << rejectCodeName(best.reject_code)
                  << " scenario_results={" << oss.str() << "}"
                  << colorReset << std::endl;
      }
      best.source_patch = patch0;
      t_rt_pick += sec_since(t_pick_begin);
      return best;
    };

    auto recordTrainingApplied = [&](int remesh_iter,
                                     int op_code,
                                     int reason_code,
                                     int category,
                                     int p0_index,
                                     int p1_index,
                                     const TrialResult& trial,
                                     bool applied,
                                     int reject_code) {
      if (!training_data_record_enabled)
        return;
      const double gain = trial.score - trial.score_before;
      const double base_weight = std::max(gain - rs.training_data_min_score_gain, 1e-12);
      const double weight = rs.training_data_weight_alpha > 0.0
                                ? std::pow(base_weight, rs.training_data_weight_alpha)
                                : 1.0;
      std::ostringstream row;
      row << time_step << ',' << remesh_iter << ','
          << op_code << ',' << reason_code << ',' << category << ','
          << p0_index << ',' << p1_index << ','
          << csvString(trial.ops) << ','
          << gain << ',' << weight << ','
          << (applied ? 1 : 0) << ',' << reject_code << ','
          << csvString(rejectCodeName(reject_code));
      appendWithHeader(training_applied_csv, training_applied_header, row.str());
      if (applied)
        ++training_iter_stats.applied_count;
      else
        ++training_iter_stats.failed_apply_count;
    };

    auto emitTrainingIterationSummary = [&]() {
      if (!training_data_record_enabled)
        return;
      std::ostringstream row;
      row << time_step << ',' << i << ','
          << training_iter_stats.candidate_edges << ','
          << training_iter_stats.trial_count << ','
          << training_iter_stats.safe_count << ','
          << training_iter_stats.selected_count << ','
          << training_iter_stats.applied_count << ','
          << training_iter_stats.failed_apply_count;
      appendWithHeader(training_iter_csv, training_iter_header, row.str());
      std::cout << "[remesh:training] step=" << time_step
                << " iter=" << i
                << " mode=" << rs.training_data_mode
                << " candidate_edges=" << training_iter_stats.candidate_edges
                << " trials=" << training_iter_stats.trial_count
                << " safe=" << training_iter_stats.safe_count
                << " selected=" << training_iter_stats.selected_count
                << " applied=" << training_iter_stats.applied_count
                << " failed_apply=" << training_iter_stats.failed_apply_count
                << " out=" << training_trial_csv.string()
                << std::endl;
    };

    auto splitCandidatePriority = [&](const networkLine* l,
                                      const std::unordered_set<networkLine*>& subsurface_trigger_set) {
      if (subsurface_trigger_set.contains(const_cast<networkLine*>(l)))
        return 0;
      if (l->BCInterface)
        return 1;
      if (is_moving_body_neumann_line(l))
        return 2;
      if (l->Dirichlet)
        return 3;
      return 4;
    };

    auto collectSplitCandidates = [&](const std::unordered_set<networkLine*>& dirty,
                                      const std::unordered_set<networkLine*>& subsurface_trigger_set) {
      std::vector<networkLine*> candidates;
      candidates.reserve(dirty.size());
      for (auto* l : dirty) {
        if (line_alive(l) && (needsSplit(l) || subsurface_trigger_set.contains(l)))
          candidates.emplace_back(l);
      }
      // Moving-body Neumann edges are processed before Dirichlet edges.  With a
      // small total_ops_limit, Dirichlet-only priority left the body-side mesh
      // too coarse near moving structures.
      std::ranges::sort(candidates, [&](const networkLine* a, const networkLine* b) {
        const int ka = splitCandidatePriority(a, subsurface_trigger_set);
        const int kb = splitCandidatePriority(b, subsurface_trigger_set);
        if (ka != kb)
          return ka < kb;
        return a->length() > b->length();
      });
      return candidates;
    };

    auto collectCollapseCandidates = [&](const std::unordered_set<networkLine*>& dirty,
                                         const std::unordered_set<networkLine*>& rejected_collapse) {
      std::vector<networkLine*> candidates;
      candidates.reserve(dirty.size());
      for (auto* l : dirty) {
        if (!line_alive(l) || rejected_collapse.count(l))
          continue;
        if (needsCollapse(l))
          candidates.emplace_back(l);
      }
      std::ranges::sort(candidates, [](const auto* a, const auto* b) {
        return a->length() < b->length();
      });
      return candidates;
    };

    struct AcceptedPatchUndoElements {
      std::unordered_set<networkPoint*> points;
      std::unordered_set<networkLine*> lines;
      std::unordered_set<networkFace*> faces;
    };

    auto collectAcceptedPatchUndoElements = [](const Network& patch) {
      AcceptedPatchUndoElements undo;
      std::unordered_set<networkPoint*> border_points;
      std::unordered_set<networkLine*> border_lines;
      for (auto* l : patch.getLines()) {
        if (!l)
          continue;
        if (l->getBoundaryFaces().size() < 2) {
          border_lines.insert(l);
          auto [p0, p1] = l->getPoints();
          if (p0)
            border_points.insert(p0);
          if (p1)
            border_points.insert(p1);
        }
      }
      for (auto* p : patch.getPoints())
        if (p && !border_points.count(p))
          undo.points.insert(p);
      for (auto* l : patch.getLines())
        if (l && !border_lines.count(l))
          undo.lines.insert(l);
      for (auto* f : patch.getFaces())
        if (f)
          undo.faces.insert(f);
      return undo;
    };

    auto isWaterlineGuardRejectCode = [](int code) {
      return code == kRejectWaterlineQuality ||
             code == kRejectBCInterfaceNoContact;
    };

    int phase2_waterline_guard_rejects_this_iter = 0;
    int phase2_waterline_trial_rejects_this_iter = 0;
    auto waterlineGuardRejectsThisIter = [&]() {
      return phase2_waterline_guard_rejects_this_iter + phase2_waterline_trial_rejects_this_iter;
    };

    struct WaterlineGuardSnapshot {
      BEMPreBVP::WaterlineMidpointStats midpoint;
      BEMPreBVP::Stats face;
      BCInterfaceContactStats raw_contact;
      BEMPreBVP::WaterlineRefreshSummary refresh;
      bool face_bad = false;
      bool raw_contact_bad = false;
      bool clamped_bad = false;
      bool hard_ok = false;
    };

    auto currentWaterlineGuardSnapshot = [&]() {
      WaterlineGuardSnapshot s;
      s.midpoint = BEMPreBVP::measure_waterline_midpoint_quality(water.getBoundaryFaces());
      s.face = measureWaterlineFaceQuality(water.getBoundaryFaces());
      s.face_bad = waterlineFaceQualityFailsLimit(s.face);
      s.raw_contact = count_bcinterface_no_contact(water, /*ignore_patch_border_lines=*/false);
      s.raw_contact_bad = s.raw_contact.no_contact() > 0;
      s.refresh = BEMPreBVP::latest_waterline_refresh_summary();
      s.clamped_bad = s.refresh.valid && s.refresh.clamped_count > 0 &&
                      s.refresh.max_move_ratio > 0.8;
      s.hard_ok = !s.face_bad && !s.raw_contact_bad &&
                  !s.clamped_bad;
      return s;
    };

    auto describeWaterlineGuardSnapshot = [](const WaterlineGuardSnapshot& s) {
      std::ostringstream oss;
      oss << "midpoint_diag_subtri=" << s.midpoint.waterline_subtri_count
          << " midpoint_diag_min_angle_deg="
          << (std::isfinite(s.midpoint.subtri_min_angle_deg) ? s.midpoint.subtri_min_angle_deg : 0.0)
          << " midpoint_diag_max_aspect=" << s.midpoint.subtri_max_aspect
          << " face_count=" << s.face.waterline_faces
          << " face_min_angle_deg="
          << (std::isfinite(s.face.waterline_face_min_angle_deg) ? s.face.waterline_face_min_angle_deg : 0.0)
          << " face_max_aspect=" << s.face.waterline_face_aspect_max
          << " raw_bci_no_contact=" << s.raw_contact.no_contact() << "/" << s.raw_contact.total()
          << " raw_point=" << s.raw_contact.point_no_contact << "/" << s.raw_contact.point_total
          << " raw_line=" << s.raw_contact.line_no_contact << "/" << s.raw_contact.line_total
          << " clamped=" << (s.refresh.valid ? s.refresh.clamped_count : 0)
          << " max_move_ratio=" << (s.refresh.valid ? s.refresh.max_move_ratio : 0.0);
      return oss.str();
    };

    auto requireCurrentWaterlineGuard = [&](const char* context) {
      const auto snapshot = currentWaterlineGuardSnapshot();
      if (snapshot.hard_ok)
        return snapshot;
      std::ostringstream oss;
      oss << context << " waterline guard failed on " << water.getName()
          << " at time_step " << time_step
          << " (" << describeWaterlineGuardSnapshot(snapshot) << ")";
      throw step_failure(oss.str());
      return snapshot;
    };

    auto rollbackAcceptedPatch = [&](const TrialResult& trial,
                                     const AcceptedPatchUndoElements& undo,
                                     const char* op_name) {
      if (!trial.source_patch)
        return false;
      trial.source_patch->copied_points = undo.points;
      trial.source_patch->copied_lines = undo.lines;
      trial.source_patch->copied_faces = undo.faces;
      const bool restored = water.replacePatch(*trial.source_patch);
      if (restored) {
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        if (curvature_remesh_enabled)
          water.computePrincipalCurvatures(false);
        const std::string rollback_reason = std::string(op_name) + "-rollback";
        refresh_contact_state_after_topology_change(rollback_reason.c_str());
      }
      return restored;
    };

    auto waterlineFaceQualityRegressed = [&](const BEMPreBVP::Stats& before,
                                             const BEMPreBVP::Stats& after) {
      if (!waterlineFaceQualityFailsLimit(after))
        return false;
      if (before.waterline_faces <= 0)
        return true;
      constexpr double eps = 1e-12;
      const bool min_angle_regressed =
          !std::isfinite(after.waterline_face_min_angle_deg) ||
          (std::isfinite(before.waterline_face_min_angle_deg) &&
           after.waterline_face_min_angle_deg < before.waterline_face_min_angle_deg - eps);
      const bool aspect_regressed =
          !std::isfinite(after.waterline_face_aspect_max) ||
          (std::isfinite(before.waterline_face_aspect_max) &&
           after.waterline_face_aspect_max > before.waterline_face_aspect_max + eps);
      return min_angle_regressed || aspect_regressed;
    };

    auto rejectAcceptedPatchIfWaterlineBad = [&](const TrialResult& trial,
                                                 const AcceptedPatchUndoElements& undo,
                                                 const char* op_name,
                                                 int aid,
                                                 int& accepted_reject_code,
                                                 const BEMPreBVP::Stats& face_quality_before) {
      accepted_reject_code = kRejectSuccess;
      if (!rs.waterline_protection_enabled)
        return false;
      const auto q = BEMPreBVP::measure_waterline_midpoint_quality(water.getBoundaryFaces());
      const auto face_quality = measureWaterlineFaceQuality(water.getBoundaryFaces());
      const bool waterline_face_quality_regressed =
          waterlineFaceQualityRegressed(face_quality_before, face_quality);
      const auto raw_contact_stats = count_bcinterface_no_contact(water, /*ignore_patch_border_lines=*/false);
      const bool raw_contact_bad = raw_contact_stats.no_contact() > 0;
      if (!waterline_face_quality_regressed && !raw_contact_bad)
        return false;
      accepted_reject_code = waterline_face_quality_regressed
                                 ? kRejectWaterlineQuality
                                 : kRejectBCInterfaceNoContact;
      std::cout << Yellow << "[remesh:accepted_guard] step=" << time_step
                << " op=" << op_name
                << " reject=" << (waterline_face_quality_regressed && raw_contact_bad ? "waterline_face_quality_regressed+raw_bci_contact" : (waterline_face_quality_regressed ? "waterline_face_quality_regressed" : "raw_bci_contact"))
                << " midpoint_diag_subtri_count=" << q.waterline_subtri_count
                << " midpoint_diag_min_angle_deg=" << (std::isfinite(q.subtri_min_angle_deg) ? q.subtri_min_angle_deg : 0.0)
                << " midpoint_diag_max_aspect=" << q.subtri_max_aspect
                << " face_count=" << face_quality.waterline_faces
                << " face_min_angle_deg=" << (std::isfinite(face_quality_before.waterline_face_min_angle_deg) ? face_quality_before.waterline_face_min_angle_deg : 0.0)
                << "->" << (std::isfinite(face_quality.waterline_face_min_angle_deg) ? face_quality.waterline_face_min_angle_deg : 0.0)
                << " face_max_aspect=" << face_quality_before.waterline_face_aspect_max
                << "->" << face_quality.waterline_face_aspect_max
                << " face_min_angle_limit=" << rs.waterline_face_min_angle_deg
                << " face_max_aspect_limit=" << rs.waterline_face_max_aspect
                << " raw_bci_no_contact=" << raw_contact_stats.no_contact() << "/" << raw_contact_stats.total()
                << " raw_point=" << raw_contact_stats.point_no_contact << "/" << raw_contact_stats.point_total
                << " raw_line=" << raw_contact_stats.line_no_contact << "/" << raw_contact_stats.line_total
                << colorReset << std::endl;
      const bool restored = rollbackAcceptedPatch(trial, undo, op_name);
      std::cout << (restored ? Green : Red)
                << "[remesh:accepted_guard] step=" << time_step
                << " op=" << op_name
                << " rollback=" << (restored ? "success" : "failed")
                << colorReset << std::endl;
      if (!restored)
        throw step_failure(std::string("accepted ") + op_name +
                           " waterline/contact guard rollback failed at time_step " +
                           std::to_string(time_step));
      ++phase2_waterline_guard_rejects_this_iter;
      requireCurrentWaterlineGuard((std::string("accepted ") + op_name +
                                    " rollback")
                                       .c_str());
      markAttemptReplaceFailed(aid);
      return true;
    };

    const auto phase2_guard_baseline = requireCurrentWaterlineGuard("phase2 baseline");
    (void)phase2_guard_baseline;

    if (stratified_training_sweep) {
      auto makeSweepRng = [&](int remesh_iter, int op_code, int category) {
        auto mix = [](std::uint64_t x) {
          x += 0x9e3779b97f4a7c15ULL;
          x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
          x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
          return x ^ (x >> 31);
        };
        const std::uint64_t seed =
            mix(static_cast<std::uint64_t>(rs.training_data_random_seed)) ^
            mix(static_cast<std::uint64_t>(time_step + 101)) ^
            mix(static_cast<std::uint64_t>(remesh_iter + 1) << 16) ^
            mix(static_cast<std::uint64_t>(op_code + 17) << 32) ^
            mix(static_cast<std::uint64_t>(category + 31) << 48);
        return std::mt19937_64(seed);
      };

      auto categoryQuotas = [&](int total_quota) {
        std::array<int, 3> q{0, 0, 0}; // Neumann, Dirichlet, BCInterface
        if (total_quota <= 0)
          return q;
        const double f_b = std::max(0.0, rs.training_data_sweep_bcinterface_fraction);
        const double f_d = std::max(0.0, rs.training_data_sweep_dirichlet_fraction);
        const double f_n = std::max(0.0, rs.training_data_sweep_neumann_fraction);
        const double sum = f_b + f_d + f_n;
        const double nb = sum > 0.0 ? f_b / sum : 0.4;
        const double nd = sum > 0.0 ? f_d / sum : 0.3;
        int q_b = static_cast<int>(std::round(total_quota * nb));
        int q_d = static_cast<int>(std::round(total_quota * nd));
        q_b = std::clamp(q_b, 0, total_quota);
        q_d = std::clamp(q_d, 0, total_quota - q_b);
        q[2] = q_b;
        q[1] = q_d;
        q[0] = std::max(0, total_quota - q_b - q_d);
        return q;
      };

      auto selectStratifiedEdges = [&](std::vector<networkLine*> candidates,
                                       int total_quota,
                                       int op_code) {
        std::vector<networkLine*> selected;
        if (total_quota <= 0 || candidates.empty())
          return selected;
        selected.reserve(std::min<int>(total_quota, static_cast<int>(candidates.size())));
        std::array<std::vector<networkLine*>, 4> by_category;
        for (auto* l : candidates) {
          if (!line_alive(l))
            continue;
          const int c = edgeTrialCategory(l);
          by_category[(0 <= c && c <= 2) ? c : 3].push_back(l);
        }
        const auto quotas = categoryQuotas(total_quota);
        std::unordered_set<networkLine*> used;
        for (int category : remesh_region_order) {
          auto& bucket = by_category[category];
          auto rng = makeSweepRng(i, op_code, category);
          std::shuffle(bucket.begin(), bucket.end(), rng);
          const int take = std::min<int>(quotas[category], static_cast<int>(bucket.size()));
          for (int k = 0; k < take; ++k) {
            selected.push_back(bucket[k]);
            used.insert(bucket[k]);
          }
        }
        std::vector<networkLine*> leftovers;
        leftovers.reserve(candidates.size());
        for (auto* l : candidates)
          if (line_alive(l) && !used.contains(l))
            leftovers.push_back(l);
        auto rng = makeSweepRng(i, op_code, 99);
        std::shuffle(leftovers.begin(), leftovers.end(), rng);
        for (auto* l : leftovers) {
          if (static_cast<int>(selected.size()) >= total_quota)
            break;
          selected.push_back(l);
        }
        return selected;
      };

      auto faceMinAngle = [](const networkFace* f) {
        if (!f)
          return std::numeric_limits<double>::infinity();
        auto [p0, p1, p2] = f->getPoints();
        if (!p0 || !p1 || !p2)
          return 0.0;
        return triangle_min_angle(p0->X, p1->X, p2->X);
      };

      auto collectFlipSweepCandidates = [&]() {
        struct LineAngle {
          networkLine* l = nullptr;
          double min_angle = std::numeric_limits<double>::infinity();
        };
        std::unordered_map<networkLine*, double> all_line_min_angle;
        std::unordered_map<networkLine*, double> poor_line_min_angle;
        const double poor_face_angle = std::max(2.0 * rs.quality_min_angle_deg * M_PI / 180.0,
                                                15.0 * M_PI / 180.0);
        auto account_line = [](std::unordered_map<networkLine*, double>& map,
                               networkLine* l,
                               double min_angle) {
          auto it = map.find(l);
          if (it == map.end())
            map.emplace(l, min_angle);
          else
            it->second = std::min(it->second, min_angle);
        };
        for (auto* f : water.getBoundaryFaces()) {
          const double min_angle = faceMinAngle(f);
          if (!std::isfinite(min_angle))
            continue;
          const bool poor = min_angle < poor_face_angle;
          for (auto* l : f->getLines()) {
            if (!line_alive(l))
              continue;
            account_line(all_line_min_angle, l, min_angle);
            if (poor)
              account_line(poor_line_min_angle, l, min_angle);
          }
        }
        const auto& source = poor_line_min_angle.empty() ? all_line_min_angle : poor_line_min_angle;
        std::unordered_map<networkLine*, double> fallback_line_min_angle;
        if (source.empty()) {
          for (auto* l : water.getBoundaryLines())
            if (line_alive(l))
              fallback_line_min_angle.emplace(l, std::numeric_limits<double>::infinity());
        }
        const auto& final_source = source.empty() ? fallback_line_min_angle : source;
        std::vector<LineAngle> tmp;
        tmp.reserve(final_source.size());
        for (const auto& [l, a] : final_source)
          tmp.push_back({l, a});
        std::ranges::sort(tmp, [](const LineAngle& a, const LineAngle& b) {
          if (a.min_angle != b.min_angle)
            return a.min_angle < b.min_angle;
          return a.l->length() > b.l->length();
        });
        std::vector<networkLine*> out;
        out.reserve(tmp.size());
        for (const auto& x : tmp)
          out.push_back(x.l);
        return out;
      };

      struct SweepBestTrial {
        TrialResult trial;
        networkLine* line = nullptr;
        int op_code = -1;
        int reason_code = -1;
        int category = 3;
        int aid = -1;
      };

      auto opName = [](int op_code) {
        switch (op_code) {
        case 0:
          return "split";
        case 1:
          return "collapse";
        case 2:
          return "flip";
        default:
          return "unknown";
        }
      };

      auto applySweepBest = [&](const SweepBestTrial& best) {
        if (!best.trial.valid || !best.trial.patch || !line_alive(best.line) ||
            trainingApplyLimitReached() || total_ops_reached())
          return false;
        if ((best.op_code == 0 && split_ops_reached()) ||
            (best.op_code == 1 && collapse_ops_reached()) ||
            (best.op_code == 2 && flip_ops_reached()))
          return false;
        auto t_rep_begin = clk::now();
        auto trial_tris_snapshot = snapshotPatchTris(*best.trial.patch);
        const auto accepted_undo = collectAcceptedPatchUndoElements(*best.trial.patch);
        const auto face_quality_before = measureWaterlineFaceQuality(water.getBoundaryFaces());
        const int p0_index = endpointIndex0(best.line);
        const int p1_index = endpointIndex1(best.line);
        ++training_iter_stats.selected_count;
        const bool replaced = water.replacePatch(*best.trial.patch);
        if (best.op_code == 0)
          t_s_replace += sec_since(t_rep_begin);
        else if (best.op_code == 1)
          t_c_replace += sec_since(t_rep_begin);
        if (!replaced) {
          recordTrainingApplied(i, best.op_code, best.reason_code, best.category,
                                p0_index, p1_index, best.trial,
                                /*applied=*/false, kRejectReplaceFailed);
          markAttemptReplaceFailed(best.aid);
          return false;
        }
        refresh_contact_state_after_topology_change(
            (std::string(opName(best.op_code)) + "-training-sweep-replacePatch").c_str(),
            /*hard_fail_on_no_contact=*/false);
        int accepted_reject_code = kRejectSuccess;
        if (rejectAcceptedPatchIfWaterlineBad(best.trial, accepted_undo,
                                              opName(best.op_code), best.aid,
                                              accepted_reject_code,
                                              face_quality_before)) {
          recordTrainingApplied(i, best.op_code, best.reason_code, best.category,
                                p0_index, p1_index, best.trial,
                                /*applied=*/false, accepted_reject_code);
          return false;
        }
        recordTrainingApplied(i, best.op_code, best.reason_code, best.category,
                              p0_index, p1_index, best.trial,
                              /*applied=*/true, kRejectSuccess);
        commitRemeshedSnapshot(trial_tris_snapshot, best.op_code, best.reason_code, best.aid);
        markAttemptSuccess(best.aid);
        topology_changed = true;
        if (best.op_code == 0)
          ++split_ops_count;
        else if (best.op_code == 1)
          ++collapse_ops_count;
        else if (best.op_code == 2)
          ++flip_ops_count;
        std::cout << "[remesh:training_sweep_apply] step=" << time_step
                  << " iter=" << i
                  << " op=" << opName(best.op_code)
                  << " category=" << remeshRegionName(best.category)
                  << " reason=" << best.reason_code
                  << " scenario=" << best.trial.ops
                  << " score_gain=" << (best.trial.score - best.trial.score_before)
                  << " split=" << split_ops_count << "/" << split_ops_limit
                  << " collapse=" << collapse_ops_count << "/" << collapse_ops_limit
                  << " flip=" << flip_ops_count << "/" << flip_ops_limit
                  << std::endl;
        return true;
      };

      auto runSweepForOp = [&](int op_code,
                               std::vector<networkLine*> candidates,
                               const std::vector<std::string>& scenarios,
                               int quota) {
        std::optional<SweepBestTrial> best;
        if (quota <= 0 || scenarios.empty() || candidates.empty())
          return best;
        candidates = selectStratifiedEdges(std::move(candidates), quota, op_code);
        std::array<int, 4> selected_by_category{};
        for (auto* l : candidates) {
          const int c = edgeTrialCategory(l);
          ++selected_by_category[(0 <= c && c <= 2) ? c : 3];
        }
        std::cout << "[remesh:training_sweep] step=" << time_step
                  << " iter=" << i
                  << " op=" << opName(op_code)
                  << " selected=" << candidates.size()
                  << " scenarios=" << scenarios.size()
                  << " N=" << selected_by_category[0]
                  << " D=" << selected_by_category[1]
                  << " B=" << selected_by_category[2]
                  << " O=" << selected_by_category[3]
                  << std::endl;
        for (auto* l : candidates) {
          if (trainingCandidateLimitReached() || !line_alive(l))
            break;
          int reason_code = 6;
          if (op_code == 0) {
            std::unordered_set<networkLine*> empty_subsurface;
            reason_code = splitReasonCodeForLine(l, empty_subsurface);
          } else if (op_code == 1) {
            reason_code = collapseReasonCodeForLine(l);
          }
          const int aid = attempt_counter++;
          const int category = edgeTrialCategory(l);
          auto t_rt_begin = clk::now();
          auto trial = runTrials(l, scenarios, false, i, op_code, reason_code);
          const double elapsed = sec_since(t_rt_begin);
          if (op_code == 0)
            t_s_runtrials += elapsed;
          else if (op_code == 1)
            t_c_runtrials += elapsed;
          else
            t_f_runtrials += elapsed;
          const int rej = trial.valid ? 0 : trial.reject_code;
          collectTriggerEdge(l, op_code, reason_code, aid, /*success=*/0, rej);
          if (!trial.valid && trial.patch) {
            if (trial.source_patch)
              collectCandidatePatch(*trial.source_patch, op_code, reason_code, aid, /*success=*/-1, rej);
            collectCandidatePatch(*trial.patch, op_code, reason_code, aid, /*success=*/0, rej);
          } else if (trial.source_patch) {
            collectCandidatePatch(*trial.source_patch, op_code, reason_code, aid, /*success=*/0, rej);
          }
          if (!trial.valid)
            continue;
          const double gain = trial.score - trial.score_before;
          if (!best || gain > best->trial.score - best->trial.score_before)
            best = SweepBestTrial{std::move(trial), l, op_code, reason_code, category, aid};
        }
        return best;
      };

      const auto t_sweep_begin = clk::now();
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);

      if (surface_split && !split_ops_reached() && !total_ops_reached()) {
        const auto subsurface_triggers = collectSubsurfaceTriggerLines();
        std::unordered_set<networkLine*> subsurface_trigger_set(
            subsurface_triggers.lines.begin(), subsurface_triggers.lines.end());
        rt_subsurface_trigger_edges += static_cast<long long>(subsurface_triggers.lines.size());
        auto best = runSweepForOp(0,
                                  collectSplitCandidates(dirty, subsurface_trigger_set),
                                  rs.split_scenarios,
                                  rs.training_data_sweep_split_edges_per_iter);
        if (best)
          applySweepBest(*best);
      }

      dirty.clear();
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);
      if (surface_collapse && !collapse_ops_reached() && !total_ops_reached()) {
        std::unordered_set<networkLine*> rejected_collapse;
        auto best = runSweepForOp(1,
                                  collectCollapseCandidates(dirty, rejected_collapse),
                                  rs.collapse_scenarios,
                                  rs.training_data_sweep_collapse_edges_per_iter);
        if (best)
          applySweepBest(*best);
      }

      if (!rs.flip_scenarios.empty() && !flip_ops_reached() && !total_ops_reached()) {
        auto best = runSweepForOp(2,
                                  collectFlipSweepCandidates(),
                                  rs.flip_scenarios,
                                  rs.training_data_sweep_flip_edges_per_iter);
        if (best)
          applySweepBest(*best);
      }
      t_split += sec_since(t_sweep_begin);

      appendProviderThetaStatsCsv();
      auto t_topo_begin = clk::now();
      water.setGeometricPropertiesForce();
      water.checkConnectivity();
      for (const auto& l : water.getLines())
        if (!l->checkTopology())
          throw step_failure("topology error after stratified training sweep at time_step " +
                             std::to_string(time_step));
      if (curvature_remesh_enabled)
        water.computePrincipalCurvatures(false);
      t_topo += sec_since(t_topo_begin);
      emitTrainingIterationSummary();
      continue;
    }

    // ------ [5a] split パス ------
    // ループ構造: 1 回の iteration で最大 1 split（batch 内最初の成功で break）。
    // ops_limit まで、閾値違反辺を非隣接 batch で処理する。break 条件は
    // "候補なし" / "batch 空" / "成功 0" / "topology エラー" のみ。
    auto t_split_begin = clk::now();
    if (surface_split && !split_ops_reached() && !total_ops_reached()) {
      int total_splits = 0;
      bool stop_splits_due_to_tiny_face = false;
      RemeshEdgeTrialStats split_stats;
      RemeshEdgeTrialStats split_total_stats;
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);

      for (const int active_category : remesh_region_order) {
        if (stop_splits_due_to_tiny_face || split_ops_reached() || total_ops_reached())
          break;
        while (!stop_splits_due_to_tiny_face &&
               !split_ops_reached() && !total_ops_reached() &&
               !trainingApplyLimitReached() && !trainingCandidateLimitReached()) {
          auto t_needs_begin = clk::now();
          const auto subsurface_triggers = collectSubsurfaceTriggerLines();
          std::unordered_set<networkLine*> subsurface_trigger_set(
              subsurface_triggers.lines.begin(), subsurface_triggers.lines.end());
          rt_subsurface_trigger_edges += static_cast<long long>(subsurface_triggers.lines.size());
          if (remesh_debug_verbose && !subsurface_triggers.lines.empty()) {
            std::cout << "[subsurface trigger] faces=" << subsurface_triggers.bad_faces
                      << " edges=" << subsurface_triggers.lines.size()
                      << " min_altitude_rel=" << subsurface_triggers.min_altitude_rel
                      << std::endl;
          }
          auto candidates = collectSplitCandidates(dirty, subsurface_trigger_set);
          keepRegionCandidates(candidates, active_category);
          trainingShuffleTail(candidates, i, 0, active_category);
          if (use_lightgbm_ranking) {
            std::stable_sort(candidates.begin(), candidates.end(),
                             [&](networkLine* a, networkLine* b) {
                               const double sa = lightGBMEdgeScore(
                                   a, rs.split_scenarios, 0,
                                   splitReasonCodeForLine(a, subsurface_trigger_set));
                               const double sb = lightGBMEdgeScore(
                                   b, rs.split_scenarios, 0,
                                   splitReasonCodeForLine(b, subsurface_trigger_set));
                               if (sa != sb)
                                 return sa > sb;
                               return a->length() > b->length();
                             });
          }
          if (candidates.empty()) {
            t_s_needs += sec_since(t_needs_begin);
            break;
          }
          auto batch = remesh_detail::collect_non_adjacent(candidates);
          t_s_needs += sec_since(t_needs_begin);
          if (batch.empty())
            break;

          // batch を走査し、最初に replacePatch が成功した辺で抜ける
          bool divided_any = false;
          for (auto* l : batch) {
            if (trainingCandidateLimitReached())
              break;
            if (!line_alive(l) || edgeTrialCategory(l) != active_category ||
                !(needsSplit(l) || subsurface_trigger_set.contains(l)))
              continue;
            auto len = l->length();

            // split の reason コード分類: 0=global_max, 1=theta, 2=dihedral, 3=obtuse, 4=length
            int reason_code = 4;
            const char* split_reason = "length";
            if (subsurface_trigger_set.contains(l)) {
              reason_code = 5;
              split_reason = "subsurface_altitude";
            } else if (len > global_max_len) {
              reason_code = 0;
              split_reason = "global_max";
            } else {
              auto v = edgeThetaVerdict(l);
              if (v == EdgeThetaVerdict::SplitCandidate) {
                reason_code = 1;
                split_reason = "theta";
              } else if (v == EdgeThetaVerdict::CurvatureInvalid) {
                reason_code = 2;
                split_reason = "dihedral";
              } else if (hasObtuseAngle(l)) {
                reason_code = 3;
                split_reason = "obtuse";
              }
            }

            int aid = attempt_counter++;
            const int trial_category = edgeTrialCategory(l);
            const auto split_scenarios =
                selectLightGBMScenarios(l, rs.split_scenarios, 0, reason_code);
            auto t_rt_begin = clk::now();
            auto trial = runTrials(l, split_scenarios, rs.aggressive_split, i, 0, reason_code);
            t_s_runtrials += sec_since(t_rt_begin);

            int rej = trial.valid ? 0 : trial.reject_code;
            collectTriggerEdge(l, /*op=*/0, reason_code, aid, /*success=*/0, rej);
            if (!trial.valid && trial.patch) {
              if (trial.source_patch)
                collectCandidatePatch(*trial.source_patch, /*op=*/0, reason_code, aid, /*success=*/-1, rej);
              collectCandidatePatch(*trial.patch, /*op=*/0, reason_code, aid, /*success=*/0, rej);
            } else if (trial.source_patch)
              collectCandidatePatch(*trial.source_patch, /*op=*/0, reason_code, aid, /*success=*/0, rej);

            if (!trial.valid) {
              if (isWaterlineGuardRejectCode(rej))
                ++phase2_waterline_trial_rejects_this_iter;
              recordBothEdgeTrialStats(split_stats, split_total_stats, trial_category, /*success=*/false, rej);
              if (split_stats.since_report >= 10)
                printEdgeTrialStats("split", split_stats, static_cast<int>(split_scenarios.size()));
              continue;
            }

            std::cout << Red << "time_step " << time_step << ": split reason=" << split_reason
                      << " len=" << len
                      << " best_ops=" << trial.ops
                      << " score_before=" << trial.score_before
                      << " score_after=" << trial.score
                      << " score_gain=" << (trial.score - trial.score_before)
                      << " ir_mean_before=" << trial.ir_mean_before
                      << " ir_mean_after=" << trial.ir_mean
                      << " ir_mean_gain=" << (trial.ir_mean - trial.ir_mean_before)
                      << " ir_worst_before=" << trial.ir_worst_before
                      << " ir_worst_after=" << trial.ir_worst
                      << " area_deviation_before=" << trial.area_deviation_before
                      << " area_deviation_after=" << trial.area_deviation
                      << " area_deviation_gain=" << (trial.area_deviation - trial.area_deviation_before)
                      << " equilateral_deviation_before=" << trial.equilateral_deviation_before
                      << " equilateral_deviation_after=" << trial.equilateral_deviation
                      << " equilateral_deviation_gain=" << (trial.equilateral_deviation - trial.equilateral_deviation_before)
                      << " hd=" << trial.hd
                      << " hd_ref=" << trial.hd_ref
                      << " hd_ref_limit=" << trial.hd_ref_limit
                      << " hd_body=" << trial.hd_body
                      << " hd_body_limit=" << trial.hd_body_limit
                      << " hd_limit_source_ref=" << (trial.hd_ref_diag_limit <= trial.hd_ref_curv_limit ? "diag" : "curv")
                      << " hd_limit_source_body=" << (trial.hd_body_diag_limit <= trial.hd_body_curv_limit ? "diag" : "curv")
                      << " hd_patch_diag=" << trial.hd_patch_diag
                      << " hd_R=" << trial.hd_curvature_radius
                      << " hd_kmax=" << trial.hd_kmax
                      << " hd_ref_diag_limit=" << trial.hd_ref_diag_limit
                      << " hd_ref_curv_limit=" << trial.hd_ref_curv_limit
                      << " hd_body_diag_limit=" << trial.hd_body_diag_limit
                      << " hd_body_curv_limit=" << trial.hd_body_curv_limit
                      << " hd_limit=" << trial.hd_limit
                      << colorReset << std::endl;
            auto t_rep_begin = clk::now();
            auto trial_tris_snapshot = snapshotPatchTris(*trial.patch);
            const auto accepted_undo = collectAcceptedPatchUndoElements(*trial.patch);
            const auto face_quality_before = measureWaterlineFaceQuality(water.getBoundaryFaces());
            const int p0_index = endpointIndex0(l);
            const int p1_index = endpointIndex1(l);
            if (training_data_record_enabled)
              ++training_iter_stats.selected_count;
            bool replaced = water.replacePatch(*trial.patch);
            t_s_replace += sec_since(t_rep_begin);
            if (replaced) {
              refresh_contact_state_after_topology_change("split-replacePatch",
                                                          /*hard_fail_on_no_contact=*/false);
              int accepted_reject_code = kRejectSuccess;
              if (rejectAcceptedPatchIfWaterlineBad(trial, accepted_undo, "split", aid,
                                                    accepted_reject_code,
                                                    face_quality_before)) {
                recordTrainingApplied(i, 0, reason_code, trial_category, p0_index, p1_index,
                                      trial, /*applied=*/false, accepted_reject_code);
                recordBothEdgeTrialStats(split_stats, split_total_stats, trial_category, /*success=*/false, accepted_reject_code);
                printEdgeTrialStats("split", split_stats, static_cast<int>(split_scenarios.size()), /*partial=*/split_stats.since_report < 10);
                dirty.erase(l);
                break;
              }
              recordTrainingApplied(i, 0, reason_code, trial_category, p0_index, p1_index,
                                    trial, /*applied=*/true, kRejectSuccess);
              recordBothEdgeTrialStats(split_stats, split_total_stats, trial_category, /*success=*/true, kRejectSuccess);
              printEdgeTrialStats("split", split_stats, static_cast<int>(split_scenarios.size()), /*partial=*/split_stats.since_report < 10);
              commitRemeshedSnapshot(trial_tris_snapshot, /*op=*/0, reason_code, aid);
              markAttemptSuccess(aid);
              divided_any = true;
              topology_changed = true;
              ++total_splits;
              ++split_ops_count;
              std::cout << "[remesh:ops_limit] step=" << time_step
                        << " split=" << split_ops_count << "/" << split_ops_limit
                        << " collapse=" << collapse_ops_count << "/" << collapse_ops_limit
                        << " flip=" << flip_ops_count << "/" << flip_ops_limit
                        << " last=split" << std::endl;
              const auto tiny_after_split = checkTinyFacesRelativeToLocalMean(
                  water, /*use_quadratic_subfaces=*/false);
              if (tiny_after_split.worst_face &&
                  tiny_after_split.worst_area_ratio < min_local_face_area_ratio) {
                const int tiny_face_index = tiny_after_split.worst_face->index;
                std::cout << Yellow << "[remesh:quality] step=" << time_step
                          << " tiny_face_after_split=" << tiny_face_index
                          << " area_ratio=" << tiny_after_split.worst_area_ratio
                          << " action=try_immediate_collapse_or_stop_splits"
                          << colorReset << std::endl;
                  if (!collapse_ops_reached() && !total_ops_reached() &&
                      collapseFaceByIndexIfPossible(water, tiny_face_index)) {
                  water.setGeometricPropertiesForce();
                  water.checkConnectivity();
                  if (curvature_remesh_enabled)
                    water.computePrincipalCurvatures(false);
                  ++collapse_ops_count;
                  std::cout << "[remesh:ops_limit] step=" << time_step
                            << " split=" << split_ops_count << "/" << split_ops_limit
                            << " collapse=" << collapse_ops_count << "/" << collapse_ops_limit
                            << " flip=" << flip_ops_count << "/" << flip_ops_limit
                            << " last=tiny_face_collapse" << std::endl;
                  std::cout << Green << "[remesh:quality] step=" << time_step
                            << " rescued_tiny_face=" << tiny_face_index
                            << " action=immediate_collapse_after_split"
                            << colorReset << std::endl;
                } else {
                  stop_splits_due_to_tiny_face = true;
		                  std::cout << Yellow << "[remesh:quality] step=" << time_step
		                            << " tiny_face=" << tiny_face_index
		                            << " action=stop_remaining_splits"
		                            << colorReset << std::endl;
                }
              }
              break;
            } else {
              recordTrainingApplied(i, 0, reason_code, trial_category, p0_index, p1_index,
                                    trial, /*applied=*/false, kRejectReplaceFailed);
              recordBothEdgeTrialStats(split_stats, split_total_stats, trial_category, /*success=*/false, kRejectReplaceFailed);
              if (split_stats.since_report >= 10)
                printEdgeTrialStats("split", split_stats, static_cast<int>(split_scenarios.size()));
              markAttemptReplaceFailed(aid);
            }
          }
          if (!divided_any) {
            break;
          }

          // 後処理: 幾何再計算、topology 検査、dirty set リフレッシュ
          auto t_post_begin = clk::now();
          water.setGeometricPropertiesForce();
          water.checkConnectivity();
          if (curvature_remesh_enabled)
            water.computePrincipalCurvatures(false);
          bool split_topo_error = false;
          for (const auto& l : water.getLines())
            if (!l->checkTopology()) {
              split_topo_error = true;
              break;
            }
          t_s_postop += sec_since(t_post_begin);
          if (split_topo_error) {
            std::cerr << Red << "[remesh:topology] step=" << time_step
                      << " error_after_split_batch=1 action=stop_splits"
                      << colorReset << std::endl;
            break;
          }

          dirty.clear();
          for (auto* ll : water.getBoundaryLines())
            dirty.insert(ll);
        }
        if (remesh_debug_verbose)
          std::cout << "[remesh:region] step=" << time_step
                    << " op=split region=" << remeshRegionName(active_category)
                    << " total_splits=" << total_splits << std::endl;
      }
      printEdgeTrialStats("split", split_stats, scenarioCountForEdgeTrialLog(rs.split_scenarios), /*partial=*/true);
      printEdgeTrialStats("split", split_total_stats, scenarioCountForEdgeTrialLog(rs.split_scenarios), /*partial=*/false, /*final=*/true);
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    {
      bool post_split_ok = true;
      for (const auto& l : water.getLines())
        if (!l->checkTopology()) {
          post_split_ok = false;
          break;
        }
      if (!post_split_ok)
        throw step_failure("topology error after division at time_step " + std::to_string(time_step));
    }

    t_split += sec_since(t_split_begin);

    // [5b] 削除: flipIfBatched は line-node remap 前の X_mid と競合するため廃止。
    // flip はパッチ内の ops（PF, PFS, CF, CFS 等）で実行され、Provider repair で戻す。
    auto t_curv_begin = clk::now();
    if (curvature_remesh_enabled)
      water.computePrincipalCurvatures(false);
    t_curv += sec_since(t_curv_begin);

    // ------ [5c] collapse パス ------
    auto t_collapse_begin = clk::now();
    // split と同じ構造: 1 iter = 最大 1 collapse。ops_limit まで候補を処理する。
    // rejected_collapse: 試行失敗 or 真の collapse ではなかった (smooth のみ) 辺を記憶して再試行しない。
    std::unordered_set<networkLine*> rejected_collapse;
    if (surface_collapse && !collapse_ops_reached() && !total_ops_reached()) {
      int total_collapses = 0;
      RemeshEdgeTrialStats collapse_stats;
      RemeshEdgeTrialStats collapse_total_stats;
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);

      for (const int active_category : remesh_region_order) {
        if (collapse_ops_reached() || total_ops_reached())
          break;
        while (!collapse_ops_reached() && !total_ops_reached() &&
               !trainingApplyLimitReached() && !trainingCandidateLimitReached()) {
          auto t_needs_begin = clk::now();
          auto candidates = collectCollapseCandidates(dirty, rejected_collapse);
          keepRegionCandidates(candidates, active_category);
          trainingShuffleTail(candidates, i, 1, active_category);
          if (use_lightgbm_ranking) {
            std::stable_sort(candidates.begin(), candidates.end(),
                             [&](networkLine* a, networkLine* b) {
                               const double sa = lightGBMEdgeScore(
                                   a, rs.collapse_scenarios, 1,
                                   collapseReasonCodeForLine(a));
                               const double sb = lightGBMEdgeScore(
                                   b, rs.collapse_scenarios, 1,
                                   collapseReasonCodeForLine(b));
                               if (sa != sb)
                                 return sa > sb;
                               return a->length() < b->length();
                             });
          }
          if (candidates.empty()) {
            t_c_needs += sec_since(t_needs_begin);
            break;
          }
          auto batch = remesh_detail::collect_non_adjacent(candidates);
          t_c_needs += sec_since(t_needs_begin);
          if (batch.empty())
            break;

          bool changed_in_batch = false;
          for (auto* l : batch) {
            if (trainingCandidateLimitReached())
              break;
            if (!line_alive(l) || edgeTrialCategory(l) != active_category || !needsCollapse(l))
              continue;
            auto [p0_orig, p1_orig] = l->getPoints();
            if (!p0_orig || !p1_orig)
              continue;

            double len = l->length();
            double theta = edgeThetaValue(l);

            // collapse の reason コード分類:
            // 0=global_min, 1=theta+grading, 2=dihedral+grading, 3=degenerate, 4=obtuse
            int reason_code = 1;
            if (len < global_min_len)
              reason_code = 0;
            else if (hasObtuseAngle(l))
              reason_code = 4;
            else {
              auto v = edgeThetaVerdict(l);
              if (v == EdgeThetaVerdict::CollapseCandidate)
                reason_code = 1;
              else if (v == EdgeThetaVerdict::CurvatureInvalid)
                reason_code = 2;
              else
                reason_code = 3;
            }

            int aid = attempt_counter++;
            const int trial_category = edgeTrialCategory(l);
            const auto collapse_scenarios =
                selectLightGBMScenarios(l, rs.collapse_scenarios, 1, reason_code);
            auto t_rt_begin = clk::now();
            auto trial = runTrials(l, collapse_scenarios, false, i, 1, reason_code);
            t_c_runtrials += sec_since(t_rt_begin);

            int rej = trial.valid ? 0 : trial.reject_code;
            collectTriggerEdge(l, /*op=*/1, reason_code, aid, /*success=*/0, rej);
            if (!trial.valid && trial.patch) {
              if (trial.source_patch)
                collectCandidatePatch(*trial.source_patch, /*op=*/1, reason_code, aid, /*success=*/-1, rej);
              collectCandidatePatch(*trial.patch, /*op=*/1, reason_code, aid, /*success=*/0, rej);
            } else if (trial.source_patch)
              collectCandidatePatch(*trial.source_patch, /*op=*/1, reason_code, aid, /*success=*/0, rej);

            if (!trial.valid) {
              if (isWaterlineGuardRejectCode(rej))
                ++phase2_waterline_trial_rejects_this_iter;
              recordBothEdgeTrialStats(collapse_stats, collapse_total_stats, trial_category, /*success=*/false, rej);
              if (collapse_stats.since_report >= 10)
                printEdgeTrialStats("collapse", collapse_stats, static_cast<int>(collapse_scenarios.size()));
              rejected_collapse.insert(l);
              continue;
            }

            std::cout << "time_step " << time_step << ": collapse len=" << len
                      << " theta=" << theta
                      << " best_ops=" << trial.ops
                      << " score_before=" << trial.score_before
                      << " score_after=" << trial.score
                      << " score_gain=" << (trial.score - trial.score_before)
                      << " ir_mean_before=" << trial.ir_mean_before
                      << " ir_mean_after=" << trial.ir_mean
                      << " ir_mean_gain=" << (trial.ir_mean - trial.ir_mean_before)
                      << " ir_worst_before=" << trial.ir_worst_before
                      << " ir_worst_after=" << trial.ir_worst
                      << " area_deviation_before=" << trial.area_deviation_before
                      << " area_deviation_after=" << trial.area_deviation
                      << " area_deviation_gain=" << (trial.area_deviation - trial.area_deviation_before)
                      << " equilateral_deviation_before=" << trial.equilateral_deviation_before
                      << " equilateral_deviation_after=" << trial.equilateral_deviation
                      << " equilateral_deviation_gain=" << (trial.equilateral_deviation - trial.equilateral_deviation_before)
                      << " hd=" << trial.hd
                      << " hd_ref=" << trial.hd_ref
                      << " hd_ref_limit=" << trial.hd_ref_limit
                      << " hd_body=" << trial.hd_body
                      << " hd_body_limit=" << trial.hd_body_limit
                      << " hd_limit_source_ref=" << (trial.hd_ref_diag_limit <= trial.hd_ref_curv_limit ? "diag" : "curv")
                      << " hd_limit_source_body=" << (trial.hd_body_diag_limit <= trial.hd_body_curv_limit ? "diag" : "curv")
                      << " hd_patch_diag=" << trial.hd_patch_diag
                      << " hd_R=" << trial.hd_curvature_radius
                      << " hd_kmax=" << trial.hd_kmax
                      << " hd_ref_diag_limit=" << trial.hd_ref_diag_limit
                      << " hd_ref_curv_limit=" << trial.hd_ref_curv_limit
                      << " hd_body_diag_limit=" << trial.hd_body_diag_limit
                      << " hd_body_curv_limit=" << trial.hd_body_curv_limit
                      << " hd_limit=" << trial.hd_limit
                      << std::endl;

            auto t_rep_begin = clk::now();
            auto trial_tris_snapshot = snapshotPatchTris(*trial.patch);
            const auto accepted_undo = collectAcceptedPatchUndoElements(*trial.patch);
            const auto face_quality_before = measureWaterlineFaceQuality(water.getBoundaryFaces());
            const int p0_index = endpointIndex0(l);
            const int p1_index = endpointIndex1(l);
            if (training_data_record_enabled)
              ++training_iter_stats.selected_count;
            bool replaced = water.replacePatch(*trial.patch);
            t_c_replace += sec_since(t_rep_begin);
            if (replaced) {
              refresh_contact_state_after_topology_change("collapse-replacePatch",
                                                          /*hard_fail_on_no_contact=*/false);
              int accepted_reject_code = kRejectSuccess;
              if (rejectAcceptedPatchIfWaterlineBad(trial, accepted_undo, "collapse", aid,
                                                    accepted_reject_code,
                                                    face_quality_before)) {
                recordTrainingApplied(i, 1, reason_code, trial_category, p0_index, p1_index,
                                      trial, /*applied=*/false, accepted_reject_code);
                recordBothEdgeTrialStats(collapse_stats, collapse_total_stats, trial_category, /*success=*/false, accepted_reject_code);
                printEdgeTrialStats("collapse", collapse_stats, static_cast<int>(collapse_scenarios.size()), /*partial=*/collapse_stats.since_report < 10);
                rejected_collapse.insert(l);
                break;
              }
              recordTrainingApplied(i, 1, reason_code, trial_category, p0_index, p1_index,
                                    trial, /*applied=*/true, kRejectSuccess);
              recordBothEdgeTrialStats(collapse_stats, collapse_total_stats, trial_category, /*success=*/true, kRejectSuccess);
              printEdgeTrialStats("collapse", collapse_stats, static_cast<int>(collapse_scenarios.size()), /*partial=*/collapse_stats.since_report < 10);
              commitRemeshedSnapshot(trial_tris_snapshot, /*op=*/1, reason_code, aid);
              markAttemptSuccess(aid);
              changed_in_batch = true;
              topology_changed = true;
              ++total_collapses;
              ++collapse_ops_count;
              std::cout << "[remesh:ops_limit] step=" << time_step
                        << " split=" << split_ops_count << "/" << split_ops_limit
                        << " collapse=" << collapse_ops_count << "/" << collapse_ops_limit
                        << " flip=" << flip_ops_count << "/" << flip_ops_limit
                        << " last=collapse" << std::endl;
              // trial.ops に 'C' が含まれない = smooth のみで collapse されなかった辺は再試行対象外
              if (trial.ops.find('C') == std::string::npos)
                rejected_collapse.insert(l);
              break;
            } else {
              recordTrainingApplied(i, 1, reason_code, trial_category, p0_index, p1_index,
                                    trial, /*applied=*/false, kRejectReplaceFailed);
              recordBothEdgeTrialStats(collapse_stats, collapse_total_stats, trial_category, /*success=*/false, kRejectReplaceFailed);
              if (collapse_stats.since_report >= 10)
                printEdgeTrialStats("collapse", collapse_stats, static_cast<int>(collapse_scenarios.size()));
              markAttemptReplaceFailed(aid);
              rejected_collapse.insert(l);
            }
          }
          if (!changed_in_batch) {
            break;
          }

          auto t_post_begin = clk::now();
          water.setGeometricPropertiesForce();
          water.checkConnectivity();
          if (curvature_remesh_enabled)
            water.computePrincipalCurvatures(false);
          t_c_postop += sec_since(t_post_begin);

          dirty.clear();
          for (auto* ll : water.getBoundaryLines())
            dirty.insert(ll);
        }
        if (remesh_debug_verbose)
          std::cout << "[remesh:region] step=" << time_step
                    << " op=collapse region=" << remeshRegionName(active_category)
                    << " total_collapses=" << total_collapses << std::endl;
      }
      printEdgeTrialStats("collapse", collapse_stats, scenarioCountForEdgeTrialLog(rs.collapse_scenarios), /*partial=*/true);
      printEdgeTrialStats("collapse", collapse_total_stats, scenarioCountForEdgeTrialLog(rs.collapse_scenarios), /*partial=*/false, /*final=*/true);
    }

    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    {
      bool post_merge_ok = true;
      for (const auto& l : water.getLines())
        if (!l->checkTopology()) {
          post_merge_ok = false;
          break;
        }
      if (!post_merge_ok)
        throw step_failure("topology error after merge at time_step " + std::to_string(time_step));
    }
    t_collapse += sec_since(t_collapse_begin);

    // ------ [5c.5] flip scenario polish ------
    // Split/collapse scenarios stay anchored on P/C.  Pure flip scenarios are run
    // here once per outer iteration, before the global flip/smooth pass.
    if (!rs.flip_scenarios.empty() && !flip_ops_reached() && !total_ops_reached()) {
      RemeshEdgeTrialStats flip_stats;
      RemeshEdgeTrialStats flip_total_stats;
      auto faceMinAngle = [](const networkFace* f) {
        if (!f)
          return std::numeric_limits<double>::infinity();
        auto [p0, p1, p2] = f->getPoints();
        if (!p0 || !p1 || !p2)
          return 0.0;
        return triangle_min_angle(p0->X, p1->X, p2->X);
      };
      const double poor_face_angle = std::max(2.0 * rs.quality_min_angle_deg * M_PI / 180.0,
                                              15.0 * M_PI / 180.0);
      struct PoorFaceCandidate {
        networkFace* f = nullptr;
        double min_angle = std::numeric_limits<double>::infinity();
      };
      std::vector<PoorFaceCandidate> poor_faces;
      poor_faces.reserve(water.getBoundaryFaces().size());
      for (auto* f : water.getBoundaryFaces()) {
        const double min_angle = faceMinAngle(f);
        if (min_angle >= poor_face_angle)
          continue;
        poor_faces.push_back({f, min_angle});
      }
      std::ranges::sort(poor_faces, [](const PoorFaceCandidate& a, const PoorFaceCandidate& b) {
        return a.min_angle < b.min_angle;
      });
      for (const int active_category : remesh_region_order) {
        if (flip_ops_reached() || total_ops_reached() ||
            trainingApplyLimitReached() || trainingCandidateLimitReached())
          break;
        for (const auto& face_entry : poor_faces) {
          if (flip_ops_reached() || total_ops_reached() ||
              trainingApplyLimitReached() || trainingCandidateLimitReached())
            break;
          auto* f = face_entry.f;
          if (!f)
            continue;
          struct FaceTrial {
            TrialResult trial;
            networkLine* l = nullptr;
            int aid = -1;
            int category = 3;
          };
          std::optional<FaceTrial> best_face_trial;
          std::unordered_set<networkLine*> seen_lines;
          int tried_edges = 0;
          int valid_edges = 0;
          for (auto* l : f->getLines()) {
            if (trainingCandidateLimitReached())
              break;
            if (!line_alive(l) ||
                edgeTrialCategory(l) != active_category || !seen_lines.insert(l).second)
              continue;
            ++tried_edges;
            const int aid = attempt_counter++;
            const int trial_category = edgeTrialCategory(l);
            const auto flip_scenarios =
                selectLightGBMScenarios(l, rs.flip_scenarios, 2, 6);
            auto t_rt_begin = clk::now();
            auto trial = runTrials(l, flip_scenarios, false, i, 2, 6);
            t_f_runtrials += sec_since(t_rt_begin);

            const int rej = trial.valid ? 0 : trial.reject_code;
            collectTriggerEdge(l, /*op=*/2, /*reason=*/6, aid, /*success=*/0, rej);
            if (!trial.valid && trial.patch) {
              if (trial.source_patch)
                collectCandidatePatch(*trial.source_patch, /*op=*/2, /*reason=*/6, aid, /*success=*/-1, rej);
              collectCandidatePatch(*trial.patch, /*op=*/2, /*reason=*/6, aid, /*success=*/0, rej);
            } else if (trial.source_patch)
              collectCandidatePatch(*trial.source_patch, /*op=*/2, /*reason=*/6, aid, /*success=*/0, rej);

            if (!trial.valid) {
              if (isWaterlineGuardRejectCode(rej))
                ++phase2_waterline_trial_rejects_this_iter;
              recordBothEdgeTrialStats(flip_stats, flip_total_stats, trial_category, /*success=*/false, rej);
              if (flip_stats.since_report >= 10)
                printEdgeTrialStats("flip", flip_stats, static_cast<int>(flip_scenarios.size()));
              continue;
            }
            ++valid_edges;
            if (!best_face_trial || trial.score > best_face_trial->trial.score)
              best_face_trial = FaceTrial{std::move(trial), l, aid, trial_category};
          }

          if (!best_face_trial)
            continue;

          std::cout << "[remesh:flip_face] step=" << time_step
                    << " poor_faces=" << poor_faces.size()
                    << " face_min_angle_deg=" << face_entry.min_angle * 180.0 / M_PI
                    << " tried_edges=" << tried_edges
                    << " valid_edges=" << valid_edges
                    << " best_ops=" << best_face_trial->trial.ops
                    << " best_score=" << best_face_trial->trial.score
                    << " best_score_gain=" << best_face_trial->trial.score - best_face_trial->trial.score_before
                    << std::endl;
          std::cout << "time_step " << time_step << ": flip face_min_angle_deg="
                    << face_entry.min_angle * 180.0 / M_PI
                    << " best_ops=" << best_face_trial->trial.ops
                    << " len=" << best_face_trial->l->length()
                    << " score_before=" << best_face_trial->trial.score_before
                    << " score_after=" << best_face_trial->trial.score
                    << " score_gain=" << (best_face_trial->trial.score - best_face_trial->trial.score_before)
                    << " ir_mean_before=" << best_face_trial->trial.ir_mean_before
                    << " ir_mean_after=" << best_face_trial->trial.ir_mean
                    << " ir_mean_gain=" << (best_face_trial->trial.ir_mean - best_face_trial->trial.ir_mean_before)
                    << " ir_worst_before=" << best_face_trial->trial.ir_worst_before
                    << " ir_worst_after=" << best_face_trial->trial.ir_worst
                    << " area_deviation_before=" << best_face_trial->trial.area_deviation_before
                    << " area_deviation_after=" << best_face_trial->trial.area_deviation
                    << " area_deviation_gain=" << (best_face_trial->trial.area_deviation - best_face_trial->trial.area_deviation_before)
                    << " equilateral_deviation_before=" << best_face_trial->trial.equilateral_deviation_before
                    << " equilateral_deviation_after=" << best_face_trial->trial.equilateral_deviation
                    << " equilateral_deviation_gain=" << (best_face_trial->trial.equilateral_deviation - best_face_trial->trial.equilateral_deviation_before)
                    << " hd=" << best_face_trial->trial.hd
                    << " hd_ref=" << best_face_trial->trial.hd_ref
                    << " hd_ref_limit=" << best_face_trial->trial.hd_ref_limit
                    << " hd_body=" << best_face_trial->trial.hd_body
                    << " hd_body_limit=" << best_face_trial->trial.hd_body_limit
                    << " hd_limit_source_ref=" << (best_face_trial->trial.hd_ref_diag_limit <= best_face_trial->trial.hd_ref_curv_limit ? "diag" : "curv")
                    << " hd_limit_source_body=" << (best_face_trial->trial.hd_body_diag_limit <= best_face_trial->trial.hd_body_curv_limit ? "diag" : "curv")
                    << " hd_patch_diag=" << best_face_trial->trial.hd_patch_diag
                    << " hd_R=" << best_face_trial->trial.hd_curvature_radius
                    << " hd_kmax=" << best_face_trial->trial.hd_kmax
                    << " hd_ref_diag_limit=" << best_face_trial->trial.hd_ref_diag_limit
                    << " hd_ref_curv_limit=" << best_face_trial->trial.hd_ref_curv_limit
                    << " hd_body_diag_limit=" << best_face_trial->trial.hd_body_diag_limit
                    << " hd_body_curv_limit=" << best_face_trial->trial.hd_body_curv_limit
                    << " hd_limit=" << best_face_trial->trial.hd_limit
                    << std::endl;
          auto trial_tris_snapshot = snapshotPatchTris(*best_face_trial->trial.patch);
          const auto accepted_undo = collectAcceptedPatchUndoElements(*best_face_trial->trial.patch);
          const auto face_quality_before = measureWaterlineFaceQuality(water.getBoundaryFaces());
          const int p0_index = endpointIndex0(best_face_trial->l);
          const int p1_index = endpointIndex1(best_face_trial->l);
          if (training_data_record_enabled)
            ++training_iter_stats.selected_count;
          if (water.replacePatch(*best_face_trial->trial.patch)) {
            refresh_contact_state_after_topology_change("flip-scenario-replacePatch",
                                                        /*hard_fail_on_no_contact=*/false);
            int accepted_reject_code = kRejectSuccess;
            if (rejectAcceptedPatchIfWaterlineBad(best_face_trial->trial, accepted_undo,
                                                  "flip", best_face_trial->aid,
                                                  accepted_reject_code,
                                                  face_quality_before)) {
              recordTrainingApplied(i, 2, 6, best_face_trial->category, p0_index, p1_index,
                                    best_face_trial->trial, /*applied=*/false, accepted_reject_code);
              recordBothEdgeTrialStats(flip_stats, flip_total_stats, best_face_trial->category, /*success=*/false, accepted_reject_code);
              printEdgeTrialStats("flip", flip_stats, scenarioCountForEdgeTrialLog(rs.flip_scenarios), /*partial=*/flip_stats.since_report < 10);
              break;
            }
            recordTrainingApplied(i, 2, 6, best_face_trial->category, p0_index, p1_index,
                                  best_face_trial->trial, /*applied=*/true, kRejectSuccess);
            recordBothEdgeTrialStats(flip_stats, flip_total_stats, best_face_trial->category, /*success=*/true, kRejectSuccess);
            printEdgeTrialStats("flip", flip_stats, scenarioCountForEdgeTrialLog(rs.flip_scenarios), /*partial=*/flip_stats.since_report < 10);
            commitRemeshedSnapshot(trial_tris_snapshot, /*op=*/2, /*reason=*/6, best_face_trial->aid);
            markAttemptSuccess(best_face_trial->aid);
            topology_changed = true;
            ++flip_ops_count;
            std::cout << "[remesh:ops_limit] step=" << time_step
                      << " split=" << split_ops_count << "/" << split_ops_limit
                      << " collapse=" << collapse_ops_count << "/" << collapse_ops_limit
                      << " flip=" << flip_ops_count << "/" << flip_ops_limit
                      << " last=flip" << std::endl;
            break;
          }

          recordTrainingApplied(i, 2, 6, best_face_trial->category, p0_index, p1_index,
                                best_face_trial->trial, /*applied=*/false, kRejectReplaceFailed);
          recordBothEdgeTrialStats(flip_stats, flip_total_stats, best_face_trial->category, /*success=*/false, kRejectReplaceFailed);
          markAttemptReplaceFailed(best_face_trial->aid);
          if (flip_stats.since_report >= 10)
            printEdgeTrialStats("flip", flip_stats, scenarioCountForEdgeTrialLog(rs.flip_scenarios));
        }
        if (remesh_debug_verbose)
          std::cout << "[remesh:region] step=" << time_step
                    << " op=flip region=" << remeshRegionName(active_category)
                    << std::endl;
      }
      printEdgeTrialStats("flip", flip_stats, scenarioCountForEdgeTrialLog(rs.flip_scenarios), /*partial=*/true);
      printEdgeTrialStats("flip", flip_total_stats, scenarioCountForEdgeTrialLog(rs.flip_scenarios), /*partial=*/false, /*final=*/true);
      water.setGeometricPropertiesForce();
      water.checkConnectivity();
      if (curvature_remesh_enabled)
        water.computePrincipalCurvatures(false);
    }

    // smoothing 共通: 面の最長辺を target line として返す
    auto longestEdgeOf = [](networkFace* f) -> networkLine* {
      networkLine* best = nullptr;
      double best_len = 0.;
      for (auto* l : f->getLines()) {
        if (!l)
          continue;
        double len = l->length();
        if (len > best_len) {
          best_len = len;
          best = l;
        }
      }
      return best;
    };

    // Historical scenario-based smoothing passes are disabled. Global
    // pass_smooth() handles smoothing uniformly at the end of the iter loop.
    // Keep this no-op lambda only so older local references remain harmless.
    [[maybe_unused]] auto runSmoothingPass = [&](const std::string& tag, std::function<bool(networkFace*)> face_predicate) {
      std::unordered_set<networkLine*> seen;
      std::vector<networkLine*> candidates;
      for (auto* f : water.getBoundaryFaces()) {
        if (!f || !face_predicate(f))
          continue;
        auto* l = longestEdgeOf(f);
        if (l && seen.insert(l).second)
          candidates.push_back(l);
      }
      if (candidates.empty())
        return;
      // Scenario-based smoothing is replaced by global pass_smooth().
      // This lambda is retained for reference but intentionally no-ops.
      (void)tag;
      (void)face_predicate;
    };

    // Profile compatibility: old smoothing-stage timers remain zero because
    // smoothing now runs in the global polish block.
    auto t_smooth_ir_begin = clk::now();
    t_smooth_ir += sec_since(t_smooth_ir_begin);
    auto t_smooth_angle_begin = clk::now();
    t_smooth_angle += sec_since(t_smooth_angle_begin);

    // [5d] 削除: flipIfBatched は廃止（5b と同じ理由）。
    auto t_curv2_begin = clk::now();
    if (curvature_remesh_enabled)
      water.computePrincipalCurvatures(false);
    t_curv += sec_since(t_curv2_begin);

    // ------ [5d] Global scenarios + repair ------
    // Patch scenarios run on copied 2-ring Networks.  Global scenarios run the
    // same F/S vocabulary on the main water Network after patch trials settle.
    // P/C are deliberately patch-only because global split/collapse needs a
    // separate candidate-edge policy.
    const bool polish_after_guard_no_progress =
        !topology_changed && waterlineGuardRejectsThisIter() > 0;
    const bool run_global_scenarios =
        !rs.global_scenarios.empty() && (topology_changed || polish_after_guard_no_progress);
    if (rs.global_scenarios.empty()) {
      if (remesh_debug_verbose) {
        std::cout << "[remesh:global] skipped reason=no_global_scenarios" << std::endl;
      }
    } else if (!run_global_scenarios) {
      if (remesh_debug_verbose) {
        std::cout << "[remesh:global] skipped reason=no_topology_change" << std::endl;
      }
    } else {
      if (polish_after_guard_no_progress) {
        std::cout << Yellow << "[remesh:phase2_no_progress] step=" << time_step
                  << " action=run_global_polish_after_waterline_guard"
                  << " waterline_guard_rejects=" << phase2_waterline_guard_rejects_this_iter
                  << " trial_guard_rejects=" << phase2_waterline_trial_rejects_this_iter
                  << colorReset << std::endl;
      }
      dump_remesh_phase_mesh("phase2_after_topology_before_global_smoothing");
      auto t_global_polish_begin = clk::now();
      const double feat_rad = rs.feature_angle_deg * M_PI / 180.0;
      const double cos_thr = rs.quality_normal_flip_cos;
      const double min_ang = rs.quality_min_angle_deg * M_PI / 180.0;
      const SizingField sf = make_uniform_field(free_surface_target_len);
      using RemeshSettings = SimulationSettings::RemeshingSettings;
      using ScenarioOp = RemeshSettings::ScenarioOp;
      using ScenarioToken = RemeshSettings::ScenarioToken;
      using FlipMode = RemeshSettings::FlipMode;
      using SmoothMode = RemeshSettings::SmoothMode;
      auto flip_mode_label = [](FlipMode m) {
        switch (m) {
        case FlipMode::ValenceStrict:
          return "valence_strict";
        case FlipMode::ValenceRelaxed:
          return "valence_relaxed";
        case FlipMode::ValenceOnly:
          return "valence_only";
        case FlipMode::DelaunayOnly:
          return "delaunay_only";
        }
        return "?";
      };
      auto smooth_mode_label = [](SmoothMode m) {
        switch (m) {
        case SmoothMode::AreaLaplacian:
          return "area_laplacian";
        case SmoothMode::OdtLaplacian:
          return "odt_laplacian";
        case SmoothMode::OdtCircumcenter:
          return "odt_circumcenter";
        case SmoothMode::CircumradiusToInradius:
          return "circumradius_to_inradius";
        }
        return "?";
      };
      auto flip_mode_for_token = [&](const ScenarioToken& token) {
        if (token.variant == 'v')
          return FlipMode::ValenceOnly;
        if (token.variant == 'd')
          return FlipMode::DelaunayOnly;
        return rs.flip_mode;
      };
      auto smooth_mode_for_token = [&](const ScenarioToken& token) {
        return token.has_smooth_mode_suffix ? token.smooth_mode : rs.smooth_mode;
      };
      auto smooth_iter_for_mode = [&](SmoothMode mode) {
        return mode == SmoothMode::OdtCircumcenter
                   ? std::max(0, rs.odt_smooth_iterations)
                   : std::max(0, rs.global_smooth_iterations);
      };
      struct GlobalSmoothStats {
        std::size_t moved = 0;
        std::size_t projected = 0;
        std::size_t dirichlet_projected = 0;
        std::size_t neumann_projected = 0;
        std::size_t projection_failed = 0;
        std::size_t boundary_points = 0;
        std::size_t sharp_points = 0;
        std::size_t bcinterface_points = 0;
        std::size_t candidate_points = 0;
        std::size_t candidate_sharp_points = 0;
        std::size_t sharp_safe_skipped = 0;
        std::size_t weighted_moved = 0;
        std::size_t ring0_moved = 0;
        std::size_t ring1_moved = 0;
        std::size_t ring2_moved = 0;
        std::size_t far_moved = 0;
        double max_weight = 1.0;
      };
      auto projected_boundary_smooth = [&](SmoothMode smooth_mode,
                                           int smooth_iter,
                                           bool sharp_safe) -> GlobalSmoothStats {
        GlobalSmoothStats stats;
        if (!trial_geometry_inputs || !trial_geometry_inputs->enabled())
          return stats;

        const auto target_provider = makeRemeshTargetProvider(*trial_geometry_inputs);
        auto touches_sharp_line = [&](const networkPoint* p) {
          if (!p)
            return false;
          for (auto* l : p->getBoundaryLines())
            if (l && l->SharpQ(feat_rad))
              return true;
          return false;
        };
        auto allows_projected_smoothing = [&](const networkPoint* p) {
          if (!p || p->BorderQ() || p->BCInterface)
            return false;
          if (!hasAnyDirichletBoundaryState(p) && !hasAnyNeumannBoundaryState(p))
            return false;
          // Keep the actual waterline endpoints pinned.  They are handled by
          // the BCInterface midpoint/geometry repair passes.
          for (auto* l : p->getBoundaryLines())
            if (l && l->BCInterface)
              return false;
          return true;
        };
        for (auto* p : water.getBoundaryPoints()) {
          if (!p)
            continue;
          ++stats.boundary_points;
          if (p->SharpQ(feat_rad))
            ++stats.sharp_points;
          if (p->BCInterface)
            ++stats.bcinterface_points;
          if (allows_projected_smoothing(p)) {
            ++stats.candidate_points;
            if (p->SharpQ(feat_rad))
              ++stats.candidate_sharp_points;
            if (sharp_safe && touches_sharp_line(p))
              ++stats.sharp_safe_skipped;
          }
        }
        auto is_feature_line = [&](const networkLine* l) {
          return l && (l->BCInterface || l->SharpQ(feat_rad));
        };
        auto incident_feature_count = [&](const networkPoint* p) {
          int count = 0;
          if (!p)
            return count;
          std::unordered_set<const networkLine*> seen;
          for (auto* l : p->getBoundaryLines()) {
            if (!l || !seen.insert(l).second)
              continue;
            if (is_feature_line(l))
              ++count;
          }
          for (auto* nb : p->getNeighbors()) {
            if (!nb)
              continue;
            for (auto* l : nb->getBoundaryLines()) {
              if (!l || !seen.insert(l).second)
                continue;
              if (is_feature_line(l))
                ++count;
            }
          }
          return count;
        };
        auto feature_proximity_weight = [&](const networkPoint* p) {
          const int n = incident_feature_count(p);
          if (n <= 0)
            return 1.0;
          return std::min(2.5, 1.0 + 0.25 * static_cast<double>(n));
        };
        auto count_moved_weight_bucket = [&](double weight) {
          if (weight >= 2.0)
            ++stats.ring0_moved;
          else if (weight >= 1.5)
            ++stats.ring1_moved;
          else if (weight > 1.0)
            ++stats.ring2_moved;
          else
            ++stats.far_moved;
        };
        auto smooth_vector = [&](networkPoint* p, Tddd& V_out) {
          return remesh_smooth_delta(p, smooth_mode, sf, feat_rad, V_out);
        };

        for (int si = 0; si < smooth_iter; ++si) {
          struct PendingSmoothUpdate {
            networkPoint* p = nullptr;
            Tddd x = {0., 0., 0.};
            double weight = 1.0;
            bool neumann = false;
          };
          std::vector<PendingSmoothUpdate> updates;
          for (auto* p : water.getBoundaryPoints()) {
            if (!allows_projected_smoothing(p))
              continue;
            if (sharp_safe && touches_sharp_line(p))
              continue;
            Tddd V;
            if (!smooth_vector(p, V) || !isFinite(V))
              continue;
            const double weight = feature_proximity_weight(p);
            V *= weight;
            stats.max_weight = std::max(stats.max_weight, weight);
            const double disp = Norm(V);
            const double limit = remesh_smooth_step_limit_ratio(smooth_mode) * localEdgeLength(p);
            if (!(disp > 1e-14) || !(limit > 0.0))
              continue;
            const Tddd V_capped = (disp <= limit) ? V : ((limit / disp) * V);
            const double alpha = find_safe_smooth_step(p, V_capped, cos_thr, min_ang);
            if (alpha <= 0.0)
              continue;

            const Tddd proposed = p->X + alpha * V_capped;
            const auto target = target_provider.queryEntityTarget(p, &water, proposed);
            if (!repairTargetAccepted(target) || !isFinite(target.target_X_clamped)) {
              ++stats.projection_failed;
              continue;
            }
            const Tddd final_x = target.target_X_clamped;
            if (!smooth_preserves_normals(p, final_x, cos_thr, min_ang)) {
              ++stats.projection_failed;
              continue;
            }
            updates.push_back({p, final_x, weight, hasAnyNeumannBoundaryState(p)});
          }
          for (const auto& update : updates) {
            update.p->setXSingle(update.x);
            ++stats.moved;
            ++stats.projected;
            if (update.weight > 1.0)
              ++stats.weighted_moved;
            count_moved_weight_bucket(update.weight);
            if (update.neumann)
              ++stats.neumann_projected;
            else
              ++stats.dirichlet_projected;
          }
          refresh_contact_state_after_topology_change("global-smoothing-iteration");
        }
        return stats;
      };
      auto run_global_repair_after_polish = [&]() {
        if (!trial_geometry_inputs || !trial_geometry_inputs->enabled()) {
          if (remesh_debug_verbose)
            std::cout << "[remesh:global_repair] skipped reason=provider_disabled" << std::endl;
          return;
        }
        std::vector<Network*> fluid_nets{&water};
        std::vector<Network*> all_objects{&water};
        for (auto* body : *trial_geometry_inputs->body_objects) {
          if (!body)
            continue;
          if (std::find(all_objects.begin(), all_objects.end(), body) == all_objects.end())
            all_objects.push_back(body);
        }

        const auto repair_report =
            BEMMeshPipeline::GeometryRepair::repairCurrentSurfaceGeometryWithProjector(
                fluid_nets,
                *trial_geometry_inputs->reference,
                *trial_geometry_inputs->body_objects,
                *trial_geometry_inputs->mesh_pipeline,
                trial_geometry_inputs->projection_mode);
        BEMMeshPipeline::GeometryRepair::finalizeProjectedGeometryRepair(
            fluid_nets,
            all_objects,
            *trial_geometry_inputs->body_objects,
            time_step,
            "post-global-scenario-repair-phase2",
            use_true_quadratic_element);
        refresh_contact_state_after_topology_change("global-repair");
        std::cout << Green << "[remesh:global_repair] step=" << time_step
                  << " repaired_points=" << repair_report.repaired_points
                  << " repaired_lines=" << repair_report.repaired_lines
                  << " damaged_faces=" << repair_report.quality_damaged_faces
                  << " waterline_repaired_points=0"
                  << " waterline_repaired_lines=" << repair_report.repaired_bcinterface_lines
                  << colorReset << std::endl;
      };
      bool global_polish_started = false;
      try {
        for (const auto& scenario : rs.global_scenarios) {
          std::vector<ScenarioToken> tokens;
          std::string parse_error;
          if (!RemeshSettings::parse_scenario_tokens(
                  scenario, RemeshSettings::ScenarioScope::Global, &tokens, &parse_error)) {
            throw step_failure("invalid normalized remesh_global_scenarios entry \"" +
                               scenario + "\": " + parse_error);
          }
          std::cout << "[remesh:global_scenario] step=" << time_step
                    << " scenario=" << scenario
                    << " tokens=" << tokens.size() << std::endl;
          for (const auto& token : tokens) {
            global_polish_started = true;
            if (token.op == ScenarioOp::Flip) {
              const auto flip_mode = flip_mode_for_token(token);
              const bool use_sharp_sector_valence =
                  token.variant == 'v' && token.valence_version == 2;
              const auto flipped = pass_flip(
                  water, flip_mode, feat_rad, cos_thr, min_ang,
                  use_sharp_sector_valence);
              std::cout << "[remesh:global_scenario] step=" << time_step
                        << " scenario=" << scenario
                        << " token=" << token.spelling
                        << " op=flip"
                        << " mode=" << flip_mode_label(flip_mode)
                        << " valence_version=" << token.valence_version
                        << " sharp_safe=1"
                        << " flipped=" << flipped << std::endl;
              dump_remesh_phase_mesh("phase2_after_global_flip");
            } else if (token.op == ScenarioOp::Smooth) {
              const auto smooth_mode = smooth_mode_for_token(token);
              const int smooth_iter = smooth_iter_for_mode(smooth_mode);
              const bool use_projected_smooth =
                  smooth_iter > 0 && trial_geometry_inputs && trial_geometry_inputs->enabled();
              const auto projected_stats =
                  use_projected_smooth
                      ? projected_boundary_smooth(smooth_mode, smooth_iter, token.sharp_safe)
                      : GlobalSmoothStats{};
              const std::size_t smooth_moved =
                  smooth_iter <= 0
                      ? 0
                      : (use_projected_smooth
                             ? projected_stats.moved
                             : pass_smooth(water, smooth_mode, sf, smooth_iter, feat_rad,
                                           cos_thr, min_ang, token.sharp_safe));
              std::cout << "[remesh:global_scenario] step=" << time_step
                        << " scenario=" << scenario
                        << " token=" << token.spelling
                        << " op=smooth"
                        << " method=" << smooth_mode_label(smooth_mode)
                        << " sharp_safe=" << (token.sharp_safe ? 1 : 0)
                        << " iterations=" << smooth_iter
                        << (use_projected_smooth ? " projected_boundary=1" : " pure_dirichlet_only=1")
                        << " moved_points=" << smooth_moved;
              if (use_projected_smooth) {
                std::cout << " projected=" << projected_stats.projected
                          << " dirichlet_projected=" << projected_stats.dirichlet_projected
                          << " neumann_projected=" << projected_stats.neumann_projected
                          << " projection_failed=" << projected_stats.projection_failed
                          << " sharp_points=" << projected_stats.sharp_points
                          << "/" << projected_stats.boundary_points
                          << " bcinterface_points=" << projected_stats.bcinterface_points
                          << " candidate_sharp_points=" << projected_stats.candidate_sharp_points
                          << "/" << projected_stats.candidate_points
                          << " sharp_safe_skipped=" << projected_stats.sharp_safe_skipped
                          << " feature_weighted=" << projected_stats.weighted_moved
                          << " feature_weight_bucket={ge2:" << projected_stats.ring0_moved
                          << ",ge1p5:" << projected_stats.ring1_moved
                          << ",gt1:" << projected_stats.ring2_moved
                          << ",far:" << projected_stats.far_moved
                          << "} max_weight=" << projected_stats.max_weight;
              }
              std::cout << std::endl;
              dump_remesh_phase_mesh("phase2_after_global_smoothing");
            } else {
              throw step_failure("global scenario unexpectedly contained patch-only token \"" +
                                 token.spelling + "\"");
            }
          }
        }
      } catch (const step_failure&) {
        throw;
      } catch (const std::exception& e) {
        std::cerr << "[remesh local_patch] global scenario polish aborted: "
                  << e.what() << std::endl;
      }
      if (global_polish_started) {
        run_global_repair_after_polish();
        dump_remesh_phase_mesh("phase2_after_global_repair");
      }
      t_smooth_angle += sec_since(t_global_polish_begin);
    }
    appendProviderThetaStatsCsv();

    // ------ [5e] topology チェック ------
    auto t_topo_begin = clk::now();
    water.setGeometricPropertiesForce();
    water.checkConnectivity();
    {
      bool post_flip_ok = true;
      for (const auto& l : water.getLines())
        if (!l->checkTopology()) {
          post_flip_ok = false;
          break;
        }
      if (!post_flip_ok)
        throw step_failure("topology error after flip at time_step " + std::to_string(time_step));
    }
    t_topo += sec_since(t_topo_begin);

    emitTrainingIterationSummary();
  }

  // Topology edits and vertex repair/smoothing are complete. Now map line
  // nodes once using the final endpoint positions.
  const auto final_line_remap_stats = finalRemapLineNodes();
  (void)final_line_remap_stats;

  // [profile] セクション別累積時間サマリ
  {
    BEMMeshPipeline::RemeshLog::StageTiming stage{
        t_split, t_collapse, t_smooth_ir, t_smooth_angle, t_curv, t_topo};
    BEMMeshPipeline::RemeshLog::RunTrialTiming trial{
        rt_calls, t_rt_ref, t_rt_copy, t_rt_run, t_rt_pick};
    BEMMeshPipeline::RemeshLog::HdLimitStats hd{
        rt_hd_checks, rt_hd_early_exits, rt_hd_samples, rt_hd_nearest_calls,
        rt_hd_target_missing, rt_hd_ref_max, rt_hd_body_max};
    BEMMeshPipeline::RemeshLog::RepairAwareTrialStats repair;
    repair.enabled = trial_geometry_inputs && trial_geometry_inputs->enabled();
    repair.source_calls = rt_repair_aware_calls;
    repair.repaired_points = rt_repaired_points;
    repair.repaired_lines = rt_repaired_lines;
    repair.repair_failed = rt_repair_failed;
    repair.altitude_rescue_attempts = rt_altitude_rescue_attempts;
    repair.altitude_rescue_successes = rt_altitude_rescue_successes;
    repair.altitude_rescue_split_worst_line = rt_altitude_rescue_split_worst_line;
    repair.altitude_rescue_split_longest_face_edge = rt_altitude_rescue_split_longest_face_edge;
    repair.subsurface_trigger_edges = rt_subsurface_trigger_edges;
    repair.best_rejected_altitude_cases = rt_best_rejected_altitude_cases;
    repair.best_rejected_altitude_before = rt_best_rejected_altitude_before;
    repair.best_rejected_altitude_after = rt_best_rejected_altitude_after;
    repair.best_rejected_altitude_gain =
        std::isfinite(rt_best_rejected_altitude_gain) ? rt_best_rejected_altitude_gain : 0.0;
    repair.max_free_gap = rt_max_free_gap;
    repair.max_body_gap = rt_max_body_gap;
    repair.max_move_ratio = rt_max_move_ratio;
    BEMMeshPipeline::RemeshLog::PhaseParts split_parts{
        t_s_needs, t_s_runtrials, t_s_replace, t_s_postop};
    BEMMeshPipeline::RemeshLog::PhaseParts collapse_parts{
        t_c_needs, t_c_runtrials, t_c_replace, t_c_postop};
    BEMMeshPipeline::RemeshLog::printLocalPatchProfile(
        time_step, stage, trial, hd, repair, split_parts, collapse_parts);
    if (t_f_runtrials > 0.0)
      std::cout << "[remesh:flip_profile] step=" << time_step
                << " runtrials=" << t_f_runtrials << "s" << std::endl;

    if (repair.enabled && remesh_debug_verbose) {
      const char* trigger_names[3] = {"Neumann", "Dirichlet", "BCInterface"};
      const char* side_names[3] = {"NeumannSide", "DirichletSide", "InvalidSide"};
      const char* line_names[4] = {"NeumannLine", "DirichletLine", "BCInterfaceLine", "InvalidLine"};
      for (int bc = 0; bc < 3; ++bc) {
        if (rt_subsurface_worst_samples[bc] == 0)
          continue;
        std::cout << "[subsurface reject detail] trigger=" << trigger_names[bc]
                  << " samples=" << rt_subsurface_worst_samples[bc]
                  << " worst_face_side={";
        for (int i = 0; i < 3; ++i) {
          if (i)
            std::cout << ",";
          std::cout << side_names[i] << ":" << rt_subsurface_worst_side_by_trigger[bc][i];
        }
        std::cout << "} worst_line_flag={";
        for (int i = 0; i < 4; ++i) {
          if (i)
            std::cout << ",";
          std::cout << line_names[i] << ":" << rt_subsurface_worst_line_flag_by_trigger[bc][i];
        }
        std::cout << "}" << std::endl;
      }
    }
  }

  if (remeshed_water_unrepaired_pvd && remeshed_water_snapshot_writer)
    remeshed_water_snapshot_writer(&water, step_retry, "remeshed", *remeshed_water_unrepaired_pvd);

  // ========================================================================
  // [6] 事後チェック: folded face、tiny face
  // ========================================================================

  {
    constexpr double max_fold_ratio = 0.15;
    auto folded = detectFoldedFaces(water, collision_settings.normal_reversal_cos, global_mean_len);
    size_t n_boundary_faces = water.getBoundaryFaces().size();
    double fold_ratio = (n_boundary_faces > 0) ? static_cast<double>(folded.size()) / n_boundary_faces : 0.0;
    if (!folded.empty())
      std::cout << Yellow << "[remesh:quality] step=" << time_step
                << " folded_faces=" << folded.size() << "/" << n_boundary_faces
                << " ratio=" << fold_ratio << colorReset << std::endl;
    if (!skip_post_remesh_quality_rejects && fold_ratio > max_fold_ratio)
      throw step_failure("post-remesh fold ratio " + std::to_string(fold_ratio) + " > " + std::to_string(max_fold_ratio) + " at time_step " + std::to_string(time_step));
  }
  {
    const auto tiny = checkTinyFacesRelativeToLocalMean(water);
    if (tiny.worst_face && tiny.worst_area_ratio < min_local_face_area_ratio) {
      std::cout << Yellow << "[remesh:quality] step=" << time_step
                << " tiny_face=" << tiny.worst_face->index
                << " area_ratio=" << tiny.worst_area_ratio
                << " area=" << tiny.worst_area
                << " local_mean_area=" << tiny.worst_local_mean_area
                << colorReset << std::endl;
      if (collapseFaceByIndexIfPossible(water, tiny.worst_face->index)) {
        std::cout << Green << "[remesh:quality] step=" << time_step
                  << " rescued_tiny_face=" << tiny.worst_face->index
                  << " action=collapse_bcinterface_allowed"
                  << colorReset << std::endl;
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
  if (rs.waterline_protection_enabled && !skip_post_remesh_quality_rejects) {
    const auto waterline_quality =
        BEMPreBVP::measure_waterline_midpoint_quality(water.getBoundaryFaces());
    if (waterline_quality.waterline_subtri_count > 0) {
      const bool min_angle_bad =
          !std::isfinite(waterline_quality.subtri_min_angle_deg) ||
          waterline_quality.subtri_min_angle_deg < rs.waterline_midpoint_min_angle_deg;
      const bool aspect_bad =
          !std::isfinite(waterline_quality.subtri_max_aspect) ||
          waterline_quality.subtri_max_aspect > rs.waterline_midpoint_max_aspect;
      if (min_angle_bad || aspect_bad) {
        std::ostringstream oss;
        oss << "post-remesh waterline midpoint quality failed on " << water.getName()
            << " at time_step " << time_step
            << " (subtri_count=" << waterline_quality.waterline_subtri_count
            << ", min_angle_deg=" << waterline_quality.subtri_min_angle_deg
            << ", max_aspect=" << waterline_quality.subtri_max_aspect
            << ", min_angle_limit=" << rs.waterline_midpoint_min_angle_deg
            << ", max_aspect_limit=" << rs.waterline_midpoint_max_aspect << ")";
        throw step_failure(oss.str());
      }
    }
    BEMPreBVP::Stats waterline_face_quality;
    for (const auto* f : water.getBoundaryFaces())
      BEMPreBVP::accumulate_face_quality(waterline_face_quality, f);
    const bool waterline_face_min_angle_bad =
        waterline_face_quality.waterline_faces > 0 &&
        (!std::isfinite(waterline_face_quality.waterline_face_min_angle_deg) ||
         waterline_face_quality.waterline_face_min_angle_deg < rs.waterline_face_min_angle_deg);
    const bool waterline_face_aspect_bad =
        waterline_face_quality.waterline_faces > 0 &&
        (!std::isfinite(waterline_face_quality.waterline_face_aspect_max) ||
         waterline_face_quality.waterline_face_aspect_max > rs.waterline_face_max_aspect);
    if (waterline_face_min_angle_bad || waterline_face_aspect_bad) {
      std::ostringstream oss;
      oss << "post-remesh waterline face quality failed on " << water.getName()
          << " at time_step " << time_step
          << " (faces=" << waterline_face_quality.waterline_faces
          << ", min_angle_deg=" << waterline_face_quality.waterline_face_min_angle_deg
          << ", max_aspect=" << waterline_face_quality.waterline_face_aspect_max
          << ", min_angle_limit=" << rs.waterline_face_min_angle_deg
          << ", max_aspect_limit=" << rs.waterline_face_max_aspect << ")";
      throw step_failure(oss.str());
    }
    const auto raw_contact_stats = count_bcinterface_no_contact(water, /*ignore_patch_border_lines=*/false);
    if (raw_contact_stats.no_contact() > 0) {
      std::ostringstream oss;
      oss << "post-remesh raw BCInterface contact failed on " << water.getName()
          << " at time_step " << time_step
          << " (no_contact=" << raw_contact_stats.no_contact() << "/" << raw_contact_stats.total()
          << ", points=" << raw_contact_stats.point_no_contact << "/" << raw_contact_stats.point_total
          << ", lines=" << raw_contact_stats.line_no_contact << "/" << raw_contact_stats.line_total
          << ")";
      throw step_failure(oss.str());
    }
    const auto& refresh = BEMPreBVP::latest_waterline_refresh_summary();
    if (refresh.valid && refresh.clamped_count > 0 && refresh.max_move_ratio > 0.8) {
      std::cout << Yellow
                << "[mesh:waterline] post-remesh repair clamped but contact/quality remained acceptable"
                << " name=" << water.getName()
                << " time_step=" << time_step
                << " clamped=" << refresh.clamped_count
                << " max_move_ratio=" << refresh.max_move_ratio
                << " midpoint_min_angle=" << waterline_quality.subtri_min_angle_deg
                << " midpoint_max_aspect=" << waterline_quality.subtri_max_aspect
                << colorReset << std::endl;
    }
  }

  // ========================================================================
  // [7] 四面体化 + 幾何プロパティ更新
  // ========================================================================

#ifdef USE_TETGEN
  if (tetrahedralize)
    water.tetrahedralize();
#endif
  water.setGeometricPropertiesForce();
  water.improveTetrahedraDelaunay();

  // パッチ VTU は remesh_debug_output=true の時だけ書く。
  remesh_debug.flushAll();

  // reject_code / bc_type 内訳サマリ (success=0 のみ集計)
  if (remesh_debug.enabled()) {
    BEMMeshPipeline::RemeshLog::printRejectSummary(
        remesh_debug.triggerEdgeRecords(), kRejectCodeCount);
  }
}

// computePrincipalCurvatures is now Network::computePrincipalCurvatures() in Network.cpp

bool collapseQualityGuardOK(const networkLine* l, const double local_mean_len_val) {
  auto surfaces = l->getBoundaryFaces();
  if (surfaces.size() != 2)
    return false;
  auto* f0 = surfaces[0];
  auto* f1 = surfaces[1];
  if (!f0 || !f1)
    return false;

  double area0 = TriangleArea(f0->getPoints()[0]->X, f0->getPoints()[1]->X, f0->getPoints()[2]->X);
  double area1 = TriangleArea(f1->getPoints()[0]->X, f1->getPoints()[1]->X, f1->getPoints()[2]->X);
  double max_edge0 = std::max({Norm(f0->getPoints()[0]->X - f0->getPoints()[1]->X),
                               Norm(f0->getPoints()[1]->X - f0->getPoints()[2]->X),
                               Norm(f0->getPoints()[2]->X - f0->getPoints()[0]->X)});
  double max_edge1 = std::max({Norm(f1->getPoints()[0]->X - f1->getPoints()[1]->X),
                               Norm(f1->getPoints()[1]->X - f1->getPoints()[2]->X),
                               Norm(f1->getPoints()[2]->X - f1->getPoints()[0]->X)});
  double alt0 = (max_edge0 > 1e-20) ? 2.0 * area0 / max_edge0 : 0.0;
  double alt1 = (max_edge1 > 1e-20) ? 2.0 * area1 / max_edge1 : 0.0;

  // Minimum altitude guard
  double alt_threshold = local_mean_len_val * 0.05;
  if (alt0 < alt_threshold || alt1 < alt_threshold)
    return false;

  // Minimum angle guard
  double min_angle0 = minimumInteriorAngleDeg(f0);
  double min_angle1 = minimumInteriorAngleDeg(f1);
  if (std::min(min_angle0, min_angle1) < 5.0)
    return false;

  // Normal direction guard (same as existing Neumann check at line 583)
  if (l->Neumann) {
    const double cos_5deg = std::cos(5.0 * M_PI / 180.0);
    if (Dot(f0->normal, -f1->normal) < cos_5deg)
      return false;
  }

  // Predict collapse outcome: check that moving endpoints to midpoint
  // does not flip normals of surrounding faces
  {
    auto [p0, p1] = l->getPoints();
    Tddd targetX = 0.5 * (p0->X + p1->X);
    if (p0->BCInterface && !p1->BCInterface)
      targetX = p0->X;
    else if (p1->BCInterface && !p0->BCInterface)
      targetX = p1->X;

    // Check all boundary faces connected to p0 or p1 (except f0, f1 which will be deleted)
    constexpr double cos_max_normal_change = std::cos(25.0 * M_PI / 180.0);
    for (auto* p : {p0, p1}) {
      for (auto* f : p->getBoundaryFaces()) {
        if (f == f0 || f == f1)
          continue;
        auto [fp0, fp1, fp2] = f->getPoints();
        int moving_idx = (fp0 == p0 || fp0 == p1) ? 0 : (fp1 == p0 || fp1 == p1) ? 1
                                                                                 : 2;
        double cos_change = surface_geometry::predictNormalChangeAfterCollapse(
            fp0->X, fp1->X, fp2->X, f->normal, targetX, moving_idx);
        if (cos_change < cos_max_normal_change)
          return false;
      }
    }
  }

  return true;
}

// ---------------------------------------------------------------------------
// findBestCollapse: collapse 前後のスコア比較で最適位置を決定
// ---------------------------------------------------------------------------
//
// surfaces_before: collapse 前の全面（位置探索の参照曲面）
// surfaces_after: collapse 後の残る面（固定2頂点のみ、candidateX で三角形を構成）
//
// 各 surfaces_before の面のパラメトリック空間上で candidateX を探索し、
// surfaces_after のスコア min(InradiusToCircumradius) が最大になる位置を返す。

bool findBestCollapse(const networkLine* l, Tddd& bestX, double global_max_len) {
  auto [p0, p1] = l->getPoints();
  if (!p0 || !p1)
    return false;

  // hard stop
  if (l->BCInterface)
    return false;
  auto adjacent = l->getBoundaryFaces();
  if (adjacent.size() != 2 || !adjacent[0] || !adjacent[1])
    return false;
  auto* f0 = adjacent[0];
  auto* f1 = adjacent[1];

  // --- collapse 前のスコア: 残る面それぞれの実際の品質の最悪値 ---
  double score_before = 1.0;
  for (auto* p : {p0, p1})
    for (auto* f : p->getBoundaryFaces()) {
      auto pts = f->getPoints(p); // pts[0]=p, pts[1],pts[2]=固定
      score_before = std::min(score_before, InradiusToCircumradius(pts[0]->X, pts[1]->X, pts[2]->X));
    }

  // --- collapse 後に残る面の情報を収集 ---
  // getPoints(p) で p を先頭にした頂点配列を返す → [0]=可動, [1],[2]=固定
  struct RemainingFace {
    Tddd current_vertex; // collapse 前の可動頂点（p0 or p1）の位置
    Tddd fixed0, fixed1; // 固定2頂点
  };
  std::vector<RemainingFace> remaining;
  for (auto* p : {p0, p1}) {
    for (auto* f : p->getBoundaryFaces()) {
      if (f == f0 || f == f1)
        continue;
      auto pts = f->getPoints(p); // pts[0]=p, pts[1],pts[2]=固定
      remaining.push_back({pts[0]->X, pts[1]->X, pts[2]->X});
    }
  }

  if (remaining.empty())
    return false;

  // --- surfaces_before: collapse 前の全面（位置探索の参照曲面）---
  struct BeforeFace {
    Tddd a, b, c;
  };
  std::vector<BeforeFace> surfaces_before;
  {
    std::unordered_set<networkFace*> face_set;
    for (auto* f : p0->getBoundaryFaces())
      face_set.insert(f);
    for (auto* f : p1->getBoundaryFaces())
      face_set.insert(f);
    for (auto* f : face_set) {
      auto [fa, fb, fc] = f->getPoints();
      surfaces_before.push_back({fa->X, fb->X, fc->X});
    }
  }

  // --- スコア評価関数（candidateX での collapse 後の品質）---
  auto evaluate_score = [&](const Tddd& candidateX) -> double {
    double min_q = 1.0;
    for (auto& rf : remaining) {
      // 辺長チェック
      if (Norm(candidateX - rf.fixed0) > global_max_len ||
          Norm(candidateX - rf.fixed1) > global_max_len)
        return 0.0;
      double q = InradiusToCircumradius(candidateX, rf.fixed0, rf.fixed1);
      min_q = std::min(min_q, q);
    }
    return min_q;
  };

  // --- 各 surfaces_before 上で格子探索 ---
  constexpr double min_acceptable_score = 0.05; // 絶対下限
  double best_score = 0.0;
  bestX = 0.5 * (p0->X + p1->X);

  for (auto& bf : surfaces_before) {
    constexpr int N = 5;
    for (int i = 0; i <= N; ++i) {
      for (int j = 0; j <= N - i; ++j) {
        double t0 = static_cast<double>(i) / N;
        double t1 = static_cast<double>(j) / N;
        Tddd X = (1.0 - t0 - t1) * bf.a + t0 * bf.b + t1 * bf.c;
        double s = evaluate_score(X);
        if (s > best_score) {
          best_score = s;
          bestX = X;
        }
      }
    }
  }

  // 改善判定:
  // 1. collapse 後のスコアが collapse 前より良い
  // 2. collapse 後のスコアが絶対下限以上
  return best_score > score_before && best_score >= min_acceptable_score;
}

// ---------------------------------------------------------------------------
// shouldFlip: flip で隣接2面の品質スコアが改善するか判定
// ---------------------------------------------------------------------------
//
// flip 前: f0=(p0,p1,p2), f1=(p0,p1,q2)
// flip 後: fA=(p0,p2,q2), fB=(p1,q2,p2)
//
// min(q) が改善するなら true。頂点は動かないので位置探索は不要。

bool shouldFlipByQuality(const networkLine* l) {
  // if (!l || l->BCInterface)
  //   return false; 2026/05/07外しした．BCInterface でも削除して問題ない．
  auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1])
    return false;

  auto [p0, p1] = l->getPoints();
  if (!p0 || !p1)
    return false;

  // f0, f1 の p0,p1 以外の頂点を取得
  auto pts0 = faces[0]->getPoints(l); // [p0_side, p1_side, p2]
  auto pts1 = faces[1]->getPoints(l); // [p1_side, p0_side, q2]
  // getPoints(line) は [a, _, b, _, c, _] 形式なので getPoints(point) を使う
  networkPoint* p2 = nullptr;
  networkPoint* q2 = nullptr;
  {
    auto [a0, a1, a2] = faces[0]->getPoints();
    p2 = (a0 != p0 && a0 != p1) ? a0 : (a1 != p0 && a1 != p1) ? a1
                                                              : a2;
    auto [b0, b1, b2] = faces[1]->getPoints();
    q2 = (b0 != p0 && b0 != p1) ? b0 : (b1 != p0 && b1 != p1) ? b1
                                                              : b2;
  }
  if (!p2 || !q2 || p2 == q2)
    return false;

  // flip 前のスコア
  double q_before_0 = InradiusToCircumradius(p0->X, p1->X, p2->X);
  double q_before_1 = InradiusToCircumradius(p0->X, p1->X, q2->X);
  double score_before = std::min(q_before_0, q_before_1);

  // flip 後のスコア
  double q_after_0 = InradiusToCircumradius(p0->X, p2->X, q2->X);
  double q_after_1 = InradiusToCircumradius(p1->X, q2->X, p2->X);
  double score_after = std::min(q_after_0, q_after_1);

  return score_after > score_before;
}
