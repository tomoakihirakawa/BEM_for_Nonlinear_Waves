// Separate translation unit for remesh_for_main_loop — the main remeshing function.
// Compiled independently to reduce main_time_domain.cpp build time.

#define BEM
#include "Network.hpp"
#include "basic_surface_geometry.hpp"
#include "BEM_remesh_main.hpp"

// 前方宣言
bool collapseQualityGuardOK(const networkLine* l, double local_mean_len_val);
bool findBestCollapse(const networkLine* l, Tddd& bestX, double global_max_len);
bool shouldFlipByQuality(const networkLine* l);
std::pair<int, int> valenceDeviationScore(const networkLine* l, int s_mean = 6);

void remesh_for_main_loop(Network& water, int time_step,
                          const SimulationSettings::RemeshingSettings& rs,
                          bool skip_post_remesh_quality_rejects,
                          const std::string& patch_output_directory,
                          double simulation_time,
                          PVDWriter* candidate_patches_pvd,
                          PVDWriter* remeshed_patches_pvd,
                          PVDWriter* trigger_edges_pvd,
                          PVDWriter* edges_pvd) {

  // rs から個別フラグを展開
  const bool tetrahedralize = rs.tetrahedralize;
  const bool surface_flip = rs.surface_flip;
  const bool surface_split = rs.surface_split;
  const bool surface_collapse = rs.surface_collapse;
  const bool surface_smoothing = rs.surface_smoothing;
  const auto& collision_settings = rs.collision;

  // ========================================================================
  // Remesh debug VTU 出力: 3 レイヤー (candidate_patches / remeshed_patches / trigger_edges) を
  // 共通 attempt_id で紐付けて蓄積、time_step 終了時にまとめて出力。
  //
  // op: 0=split, 1=collapse, 2=smooth
  // reason (op=split):    0=global_max, 1=theta, 2=dihedral, 3=obtuse
  // reason (op=collapse): 0=global_min, 1=theta+grading, 2=dihedral+grading, 3=degenerate, 4=obtuse
  // reason (op=smooth):   0 (分類なし)
  // success: 0=失敗 / 1=replacePatch 成功
  //
  // スレッド安全性: runTrials 内部の omp parallel for はレコードに触れない。
  // 外側の for (batch) ループは逐次なので record vectors / attempt_counter に mutex は不要。
  // ========================================================================
  int attempt_counter = 0;

  struct CandidatePatchRecord {
    T3Tddd tri;
    int attempt_id;
    int op;
    int reason;
    int success;
    int reject_code;  // 0=success, 1=no_valid, 2=hd_exceeded, 3=no_improve, 4=replace_failed
  };
  struct RemeshedPatchRecord {
    T3Tddd tri;
    int attempt_id;
    int op;
    int reason;
  };
  struct TriggerEdgeRecord {
    Tddd p0, p1;
    int attempt_id;
    int op;
    int reason;
    int success;
    int reject_code;
    int bc_type;  // 0=Neumann, 1=Dirichlet, 2=CORNER
  };
  std::vector<CandidatePatchRecord> candidate_patch_records;
  std::vector<RemeshedPatchRecord>  remeshed_patch_records;
  std::vector<TriggerEdgeRecord>    trigger_edge_records;

  auto collectCandidatePatch = [&](const Network& patch, int op, int reason, int aid,
                                    int success, int reject_code) {
    for (auto* f : patch.getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      candidate_patch_records.push_back({{p0->X, p1->X, p2->X}, aid, op, reason, success, reject_code});
    }
  };
  auto collectRemeshedPatch = [&](const Network& patch, int op, int reason, int aid) {
    for (auto* f : patch.getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      remeshed_patch_records.push_back({{p0->X, p1->X, p2->X}, aid, op, reason});
    }
  };
  // trial.patch は replacePatch 成功時に water へ face を移管するため、
  // replacePatch 呼び出し後は getBoundaryFaces() が空になる。
  // そのため座標を事前にスナップショットし、成功時のみコミットする。
  auto snapshotPatchTris = [](const Network& patch) -> std::vector<T3Tddd> {
    std::vector<T3Tddd> tris;
    for (auto* f : patch.getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      tris.push_back({p0->X, p1->X, p2->X});
    }
    return tris;
  };
  auto commitRemeshedSnapshot = [&](const std::vector<T3Tddd>& tris, int op, int reason, int aid) {
    for (const auto& t : tris)
      remeshed_patch_records.push_back({t, aid, op, reason});
  };
  auto collectTriggerEdge = [&](const networkLine* l, int op, int reason, int aid,
                                 int success, int reject_code) {
    auto [p0, p1] = l->getPoints();
    if (p0 && p1) {
      int bc = l->CORNER ? 2 : (l->Dirichlet ? 1 : 0);
      trigger_edge_records.push_back({p0->X, p1->X, aid, op, reason, success, reject_code, bc});
    }
  };
  // 試行成功後: success=1 / reject_code=0 に更新
  auto markAttemptSuccess = [&](int aid) {
    for (auto& r : trigger_edge_records)
      if (r.attempt_id == aid) { r.success = 1; r.reject_code = 0; }
    for (auto& r : candidate_patch_records)
      if (r.attempt_id == aid) { r.success = 1; r.reject_code = 0; }
  };
  // replacePatch 失敗時: reject_code = 4 に更新
  auto markAttemptReplaceFailed = [&](int aid) {
    for (auto& r : trigger_edge_records)
      if (r.attempt_id == aid) r.reject_code = 4;
    for (auto& r : candidate_patch_records)
      if (r.attempt_id == aid) r.reject_code = 4;
  };

  // 三角 VTU 共通書き出し: CellData は Int32 (attempt_id, op, reason, [success])
  auto writeTrianglePatchVTU = [&](const std::string& tag,
                                   const auto& records,
                                   bool has_success,
                                   PVDWriter* pvd) {
    if (!pvd && patch_output_directory.empty())
      return;
    const std::string vtu_name = water.getName() + "_" + tag + "_" + std::to_string(time_step) + ".vtu";
    const std::string filename = patch_output_directory.empty()
                                     ? "/tmp/" + vtu_name
                                     : patch_output_directory + "/" + vtu_name;
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp)
      return;
    int n_cells = records.size();
    int n_points = n_cells * 3;
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", n_cells, n_points);
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' format='ascii'>\n");
    for (auto& r : records) {
      auto& [a, b, c] = r.tri;
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(a), (float)std::get<1>(a), (float)std::get<2>(a));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(b), (float)std::get<1>(b), (float)std::get<2>(b));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(c), (float)std::get<1>(c), (float)std::get<2>(c));
    }
    fprintf(fp, "</DataArray>\n</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d %d %d\n", i*3, i*3+1, i*3+2);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d\n", (i+1)*3);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "5\n"); // VTK_TRIANGLE
    fprintf(fp, "</DataArray>\n</Cells>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type='Int32' Name='attempt_id' format='ascii'>\n");
    for (auto& r : records) fprintf(fp, "%d\n", r.attempt_id);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='op' format='ascii'>\n");
    for (auto& r : records) fprintf(fp, "%d\n", r.op);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='reason' format='ascii'>\n");
    for (auto& r : records) fprintf(fp, "%d\n", r.reason);
    fprintf(fp, "</DataArray>\n");
    if (has_success) {
      fprintf(fp, "<DataArray type='Int32' Name='success' format='ascii'>\n");
      for (auto& r : records) {
        if constexpr (requires { r.success; })
          fprintf(fp, "%d\n", r.success);
      }
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "<DataArray type='Int32' Name='reject_code' format='ascii'>\n");
      for (auto& r : records) {
        if constexpr (requires { r.reject_code; })
          fprintf(fp, "%d\n", r.reject_code);
      }
      fprintf(fp, "</DataArray>\n");
    }
    fprintf(fp, "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(fp);
    std::cout << "[patch debug] wrote " << filename << " (" << n_cells << " faces)" << std::endl;
    if (pvd) {
      pvd->push(vtu_name, simulation_time);
      pvd->output();
    }
  };

  auto flushCandidatePatchVTU = [&]() {
    writeTrianglePatchVTU("candidate_patches", candidate_patch_records, /*has_success=*/true, candidate_patches_pvd);
  };
  auto flushRemeshedPatchVTU = [&]() {
    writeTrianglePatchVTU("remeshed_patches", remeshed_patch_records, /*has_success=*/false, remeshed_patches_pvd);
  };

  auto flushTriggerEdgesVTU = [&]() {
    if (!trigger_edges_pvd && patch_output_directory.empty())
      return;
    const std::string vtu_name = water.getName() + "_trigger_edges_" + std::to_string(time_step) + ".vtu";
    const std::string filename = patch_output_directory.empty()
                                     ? "/tmp/" + vtu_name
                                     : patch_output_directory + "/" + vtu_name;
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp)
      return;
    int n_cells = trigger_edge_records.size();
    int n_points = n_cells * 2;
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", n_cells, n_points);
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' format='ascii'>\n");
    for (auto& r : trigger_edge_records) {
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(r.p0), (float)std::get<1>(r.p0), (float)std::get<2>(r.p0));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(r.p1), (float)std::get<1>(r.p1), (float)std::get<2>(r.p1));
    }
    fprintf(fp, "</DataArray>\n</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d %d\n", i*2, i*2+1);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d\n", (i+1)*2);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "3\n"); // VTK_LINE
    fprintf(fp, "</DataArray>\n</Cells>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type='Int32' Name='attempt_id' format='ascii'>\n");
    for (auto& r : trigger_edge_records) fprintf(fp, "%d\n", r.attempt_id);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='op' format='ascii'>\n");
    for (auto& r : trigger_edge_records) fprintf(fp, "%d\n", r.op);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='reason' format='ascii'>\n");
    for (auto& r : trigger_edge_records) fprintf(fp, "%d\n", r.reason);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='success' format='ascii'>\n");
    for (auto& r : trigger_edge_records) fprintf(fp, "%d\n", r.success);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='reject_code' format='ascii'>\n");
    for (auto& r : trigger_edge_records) fprintf(fp, "%d\n", r.reject_code);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='bc_type' format='ascii'>\n");
    for (auto& r : trigger_edge_records) fprintf(fp, "%d\n", r.bc_type);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(fp);
    std::cout << "[patch debug] wrote " << filename << " (" << n_cells << " edges)" << std::endl;
    if (trigger_edges_pvd) {
      trigger_edges_pvd->push(vtu_name, simulation_time);
      trigger_edges_pvd->output();
    }
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
      if (l->CORNER) { ++n_corner; }
      else if (l->Dirichlet) { ++n_dirichlet; fs_len_min = std::min(fs_len_min, L); fs_len_max = std::max(fs_len_max, L); fs_len_sum += L; }
      else { ++n_neumann; wall_len_min = std::min(wall_len_min, L); wall_len_max = std::max(wall_len_max, L); wall_len_sum += L; }
    }
    const double fs_thr_split    = rs.len_fs_split_ratio    * free_surface_target_len;
    const double fs_thr_collapse = rs.len_fs_collapse_ratio * free_surface_target_len;
    std::cout << "[remesh target] horiz_diag=" << horiz_diag_dbg
              << " divisor=" << rs.len_target_divisor
              << " fs_target=" << free_surface_target_len
              << " fs_split_thr=" << fs_thr_split
              << " fs_collapse_thr=" << fs_thr_collapse
              << " global_max=" << global_max_len
              << " global_min=" << global_min_len << std::endl;
    std::cout << "[remesh edges] Dirichlet=" << n_dirichlet
              << " (len min=" << fs_len_min
              << " max=" << fs_len_max
              << " mean=" << (n_dirichlet ? fs_len_sum / n_dirichlet : 0.0) << ")"
              << "  Neumann=" << n_neumann
              << " (len min=" << wall_len_min
              << " max=" << wall_len_max
              << " mean=" << (n_neumann ? wall_len_sum / n_neumann : 0.0) << ")"
              << "  CORNER=" << n_corner << std::endl;
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
    std::cout << Yellow << "[remesh] time_step " << time_step
              << ": collision unresolved (fold_ratio=" << fold_ratio
              << ", protected_ratio=" << protected_ratio << "), continuing"
              << colorReset << std::endl;
  }

  const auto protected_halo_lines = remesh_detail::build_protection_halo(protected_lines);

  water.setGeometricPropertiesForce();
  water.checkConnectivity();

  const std::size_t boundary_line_count = water.getBoundaryLines().size();
  if (boundary_line_count > 0 && protected_lines.size() > 0) {
    std::cout << Yellow << "[remesh] time_step " << time_step
              << ": protected region is " << protected_lines.size() << " / " << boundary_line_count
              << " boundary lines (" << protected_halo_lines.size() << " incl. 1-ring halo)"
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
      std::cerr << Red << "[remesh] time_step " << time_step
                << ": " << n_pre_topo_err << " pre-existing topology errors detected"
                << colorReset << std::endl;
      throw step_failure("pre-existing topology errors (" + std::to_string(n_pre_topo_err) + " lines) at time_step " + std::to_string(time_step));
    }
  }

  // ========================================================================
  // [5] split / flip / collapse ループ（iter_split_collapse 回繰り返し）
  // ========================================================================

  const int iter_split_collapse = std::max(1, rs.iter_split_collapse);

  // [profile] 各セクションの累積時間
  double t_split = 0., t_collapse = 0., t_smooth_ir = 0., t_smooth_angle = 0., t_topo = 0., t_curv = 0.;
  // runTrials の内訳
  double t_rt_ref = 0., t_rt_copy = 0., t_rt_run = 0., t_rt_pick = 0.;
  int rt_calls = 0;
  // collapse 内訳
  double t_c_needs = 0., t_c_runtrials = 0., t_c_replace = 0., t_c_postop = 0.;
  // split 内訳
  double t_s_needs = 0., t_s_runtrials = 0., t_s_replace = 0., t_s_postop = 0.;
  using clk = std::chrono::high_resolution_clock;
  auto sec_since = [](const clk::time_point& t0) {
    return std::chrono::duration<double>(clk::now() - t0).count();
  };

  for (auto i = 0; i < iter_split_collapse; i++) {
    auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };

    // feature 保護の閾値（これ以上の法線変化は feature として保護）
    const double feature_angle = rs.feature_angle_deg * M_PI / 180.0;

    // --- θ 判定（split/collapse 共通） ---
    // 辺中点の曲率を使い、端点曲率の max による過剰 split を回避
    auto edgeThetaVerdict = [&](const networkLine* l) -> EdgeThetaVerdict {
      if (!l->geom_curvature.valid)
        return EdgeThetaVerdict::CurvatureInvalid;
      auto [p0, p1] = l->getPoints();
      Tddd edge = p1->X - p0->X;
      double theta = surface_geometry::edgeCurvatureAngle(edge, l->geom_curvature.k1, l->geom_curvature.k2, l->geom_curvature.PD1, l->geom_curvature.PD2);
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
        if (!f) continue;
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

    // --- needsSplit: この辺を split すべきか ---
    // 自由表面 (Dirichlet/CORNER) は free_surface_target_len のヒステリシス上限を満たさない限り却下
    // それ以外の辺は曲率忠実度 / 絶対長さで判定
    auto needsSplit = [&](networkLine* l) -> bool {
      auto len = l->length();

      // 自由表面辺: len > len_fs_split_ratio * fs_target を満たさないと split しない
      if (l->Dirichlet || l->CORNER) {
        return len > rs.len_fs_split_ratio * free_surface_target_len;
      }

      // 一般辺: 絶対長さ上限超過
      if (len > global_max_len)
        return true;

      // 一般辺: 曲率忠実度
      auto verdict = edgeThetaVerdict(l);
      if (verdict == EdgeThetaVerdict::SplitCandidate)
        return true;

      return false;
    };

    // --- needsCollapse: この辺を collapse すべきか ---
    // 自由表面 (Dirichlet/CORNER) は free_surface_target_len のヒステリシス下限を下回らない限り却下
    // それ以外の辺は曲率忠実度 / 絶対長さで判定
    auto needsCollapse = [&](networkLine* l) -> bool {
      auto len = l->length();

      // 自由表面辺: len < len_fs_collapse_ratio * fs_target を満たさないと collapse しない
      if (l->Dirichlet || l->CORNER) {
        return len < rs.len_fs_collapse_ratio * free_surface_target_len;
      }

      // 一般辺: 絶対長さ下限未満
      if (len < global_min_len)
        return true;

      // 一般辺: 曲率忠実度
      auto verdict = edgeThetaVerdict(l);
      if (verdict == EdgeThetaVerdict::CollapseCandidate)
        return true;

      return false;
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
      for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d %d\n", i * 2, i * 2 + 1);
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
      for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d\n", (i + 1) * 2);
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
      for (int i = 0; i < n_cells; ++i) fprintf(fp, "3\n"); // VTK_LINE
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "</Cells>\n");
      fprintf(fp, "<CellData>\n");

      auto write_scalar = [&](const char* name, auto extractor) {
        fprintf(fp, "<DataArray type='Float32' Name='%s' NumberOfComponents='1' format='ascii'>\n", name);
        for (auto* l : lines) fprintf(fp, "%f\n", (float)extractor(l));
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
        if (!l->geom_curvature.valid) return 0.0;
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
        for (auto* f : faces) if (f) n += f->area * f->normal;
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

    if (i == 0) writeEdgesVTU();

    // ------ 共通ヘルパー（split / collapse 共用） ------

    // --- OldFaceData: コピー元の面の6節点データ ---
    struct OldFaceData {
      T3Tddd tri;                   // 3頂点座標 {p0, p1, p2}
      std::array<Tddd, 3> X_mids;   // 3辺の X_mid {mid01, mid12, mid20}
      std::array<double, 6> phi6;   // 6節点 phi
      std::array<double, 6> phi_t6; // 6節点 phi_t
    };

    // ops 実行前の patch 面データを保存し、投影先候補として使う
    // linear 要素では辺節点 DOF を持たない前提なので、辺節点の幾何/phi/phi_t は端点平均で与える
    // quadratic 要素では辺節点の値をそのまま使う
    auto collectOldFaceData = [](const Network& net) -> std::vector<OldFaceData> {
      std::vector<OldFaceData> result;
      for (auto* f : net.getBoundaryFaces()) {
        auto [fp0, fp1, fp2] = f->getPoints();
        auto [fl0, fl1, fl2] = f->getLines();
        // 辺節点の値: midpoint_index >= 0（quadratic）ならそのまま使う
        //             midpoint_index < 0（linear）なら端点平均
        auto edgeMidX = [](const networkLine* l) -> Tddd {
          if (l->midpoint_index >= 0)
            return l->X_mid;
          auto [pA, pB] = l->getPoints();
          return 0.5 * (pA->X + pB->X);
        };
        auto edgePhi = [](const networkLine* l) -> double {
          if (l->midpoint_index >= 0)
            return std::get<0>(l->phiphin);
          auto [pA, pB] = l->getPoints();
          return 0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
        };
        auto edgePhiT = [](const networkLine* l) -> double {
          if (l->midpoint_index >= 0)
            return std::get<0>(l->phiphin_t);
          auto [pA, pB] = l->getPoints();
          return 0.5 * (std::get<0>(pA->phiphin_t) + std::get<0>(pB->phiphin_t));
        };
        result.push_back({{fp0->X, fp1->X, fp2->X},
                          {edgeMidX(fl0), edgeMidX(fl1), edgeMidX(fl2)},
                          {std::get<0>(fp0->phiphin), std::get<0>(fp1->phiphin), std::get<0>(fp2->phiphin),
                           edgePhi(fl0), edgePhi(fl1), edgePhi(fl2)},
                          {std::get<0>(fp0->phiphin_t), std::get<0>(fp1->phiphin_t), std::get<0>(fp2->phiphin_t),
                           edgePhiT(fl0), edgePhiT(fl1), edgePhiT(fl2)}});
      }
      return result;
    };

    // old_faces 全体から2次曲面上の最近点パラメータ (t0,t1) を求める
    // 1. linear Nearest_ で t0t1 の初期値を取得
    // 2. Newton 反復で2次曲面上の正確な t0t1 を取得（node relocation と同じ手法）
    auto projectOntoQuadSurface = [](const Tddd& X, const std::vector<OldFaceData>& old_faces)
        -> std::tuple<Tddd, double, double, int> {
      // 6節点座標から TriShape<6> で曲面上の点を計算
      auto evalQuad = [](double t0, double t1, const OldFaceData& face) -> Tddd {
        auto N = TriShape<6>(t0, t1);
        auto& [v0, v1, v2] = face.tri;
        auto& [m01, m12, m20] = face.X_mids;
        return {
            N[0] * std::get<0>(v0) + N[1] * std::get<0>(v1) + N[2] * std::get<0>(v2) + N[3] * std::get<0>(m01) + N[4] * std::get<0>(m12) + N[5] * std::get<0>(m20),
            N[0] * std::get<1>(v0) + N[1] * std::get<1>(v1) + N[2] * std::get<1>(v2) + N[3] * std::get<1>(m01) + N[4] * std::get<1>(m12) + N[5] * std::get<1>(m20),
            N[0] * std::get<2>(v0) + N[1] * std::get<2>(v1) + N[2] * std::get<2>(v2) + N[3] * std::get<2>(m01) + N[4] * std::get<2>(m12) + N[5] * std::get<2>(m20)};
      };
      // D_TriShape<6> 微分で曲面の接線ベクトルを計算
      auto evalQuadDeriv = [](double t0, double t1, int di, int dj, const OldFaceData& face) -> Tddd {
        auto& [v0, v1, v2] = face.tri;
        auto& [m01, m12, m20] = face.X_mids;
        std::array<double, 6> DN;
        if (di == 1 && dj == 0)
          DN = D_TriShape<6, 1, 0>(t0, t1);
        else if (di == 0 && dj == 1)
          DN = D_TriShape<6, 0, 1>(t0, t1);
        else if (di == 2 && dj == 0)
          DN = D_TriShape<6, 2, 0>(t0, t1);
        else if (di == 0 && dj == 2)
          DN = D_TriShape<6, 0, 2>(t0, t1);
        else
          DN = D_TriShape<6, 1, 1>(t0, t1);
        return {
            DN[0] * std::get<0>(v0) + DN[1] * std::get<0>(v1) + DN[2] * std::get<0>(v2) + DN[3] * std::get<0>(m01) + DN[4] * std::get<0>(m12) + DN[5] * std::get<0>(m20),
            DN[0] * std::get<1>(v0) + DN[1] * std::get<1>(v1) + DN[2] * std::get<1>(v2) + DN[3] * std::get<1>(m01) + DN[4] * std::get<1>(m12) + DN[5] * std::get<1>(m20),
            DN[0] * std::get<2>(v0) + DN[1] * std::get<2>(v1) + DN[2] * std::get<2>(v2) + DN[3] * std::get<2>(m01) + DN[4] * std::get<2>(m12) + DN[5] * std::get<2>(m20)};
      };

      double min_dist = 1e+20;
      Tddd best_near = X;
      double best_t0 = 0, best_t1 = 0;
      int best_idx = -1;
      for (int i = 0; i < (int)old_faces.size(); ++i) {
        // 1. linear Nearest_ で初期 t0t1
        auto [t0_init, t1_init, near_linear, normal] = Nearest_(X, old_faces[i].tri);
        // 2. refineNearestParam で曲面上の正確な t0t1（共通 Newton ヘルパ）
        auto [t0, t1] = refineNearestParam(
            X, Tdd{t0_init, t1_init},
            [&](double a, double b) { return evalQuad(a, b, old_faces[i]); },
            [&](double a, double b, int di, int dj) { return evalQuadDeriv(a, b, di, dj, old_faces[i]); });
        Tddd near_quad = evalQuad(t0, t1, old_faces[i]);
        double dist = Norm(near_quad - X);
        if (dist < min_dist) {
          min_dist = dist;
          best_near = near_quad;
          best_t0 = t0;
          best_t1 = t1;
          best_idx = i;
        }
      }
      return {best_near, best_t0, best_t1, best_idx};
    };

    // 内部頂点を先に投影し、その更新後の頂点位置から内部辺節点候補を作って投影する
    // 投影先は orig_faces 全体から距離最小の候補面を選ぶ
    auto projectOntoOldSurface = [&](Network& net, const std::vector<OldFaceData>& orig_faces) {
      if (orig_faces.empty())
        return;
      // 外周点を特定（面が < 2 の辺の端点）
      std::unordered_set<networkPoint*> bset;
      for (auto* ll : net.getLines())
        if (ll->getBoundaryFaces().size() < 2) {
          auto [p0, p1] = ll->getPoints();
          if (p0)
            bset.insert(p0);
          if (p1)
            bset.insert(p1);
        }

      // 内部頂点を投影 + phi/phi_t 補間（外周点のみ固定）
      for (auto* p : net.getPoints()) {
        if (bset.count(p))
          continue;
        auto [near, t0, t1, idx] = projectOntoQuadSurface(p->X, orig_faces);
        if (idx >= 0) {
          p->setXSingle(near);
          auto N = TriShape<6>(t0, t1);
          double phi = 0., phi_t = 0.;
          for (int k = 0; k < 6; ++k) {
            phi += N[k] * orig_faces[idx].phi6[k];
            phi_t += N[k] * orig_faces[idx].phi_t6[k];
          }
          std::get<0>(p->phiphin) = phi;
          std::get<0>(p->phiphin_t) = phi_t;
        }
      }

      // 内部辺節点（X_mid）を、移動後の端点平均から投影
      for (auto* l : net.getLines()) {
        if (l->getBoundaryFaces().size() < 2)
          continue;
        auto [pA, pB] = l->getPoints();
        Tddd X_mid_candidate = 0.5 * (pA->X + pB->X);
        auto [near, t0, t1, idx] = projectOntoQuadSurface(X_mid_candidate, orig_faces);
        if (idx >= 0) {
          l->setXSingle(near);
          auto N = TriShape<6>(t0, t1);
          double phi = 0., phi_t = 0.;
          for (int k = 0; k < 6; ++k) {
            phi += N[k] * orig_faces[idx].phi6[k];
            phi_t += N[k] * orig_faces[idx].phi_t6[k];
          }
          std::get<0>(l->phiphin) = phi;
          std::get<0>(l->phiphin_t) = phi_t;
        }
      }
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
    // 新スコア式（target-distance ベース）:
    //   score = worst_ir
    //         + 1 / (1 + local_cv_mean)
    //         - w * (fs_target_dist + 0.5 * theta_target_dist)
    //
    //   fs_target_dist    = mean_{Dirichlet/CORNER edges} ( log(len / fs_target) )^2
    //   theta_target_dist = mean_{valid curvature edges}  ( log(theta / theta_target) )^2
    //
    // 特徴:
    //   - 過大・過小の両方を penalize（log space で対称）
    //   - 大きな逸脱は 2 乗で強く penalize（旧 1/(1+x) の飽和を回避）
    //   - mean (not max) なので 1 本の悪い辺に支配されず、全体が target に寄る方向へ
    //   - w = rs.quality_resolution_weight（既定 1.0）で重みを調整可
    auto patchQuality = [&](const Network& net) {
      constexpr double kFlipSentinel = -1.0e18;

      // ----- 反転チェック: 隣接面の法線が逆を向いていれば即却下 -----
      for (auto* l : net.getLines()) {
        auto faces = l->getBoundaryFaces();
        if (faces.size() == 2 && faces[0] && faces[1]) {
          if (Dot(faces[0]->normal, faces[1]->normal) < rs.quality_normal_flip_cos)
            return kFlipSentinel;
        }
      }

      // ----- 面ループ: 形状品質 -----
      double worst_ir = 1.0;
      for (auto* f : net.getBoundaryFaces()) {
        auto [p0, p1, p2] = f->getPoints();
        worst_ir = std::min(worst_ir, InradiusToCircumradius(p0->X, p1->X, p2->X));
      }

      // ----- 辺ループ: 目標距離 (log^2 ratio の平均) -----
      double fs_dist_sum = 0.;   int fs_count = 0;
      double th_dist_sum = 0.;   int th_count = 0;
      for (auto* l : net.getBoundaryLines()) {
        // Dirichlet/CORNER 辺 → free surface target との距離
        if ((l->Dirichlet || l->CORNER) && free_surface_target_len > 1e-20) {
          double L = l->length();
          if (L > 1e-20) {
            double lr = std::log(L / free_surface_target_len);  // 長いほど正、短いほど負
            fs_dist_sum += lr * lr;
            ++fs_count;
          }
        }
        // 曲率 valid 辺 → theta_target との距離
        if (!l->geom_curvature.valid)
          continue;
        auto [p0, p1] = l->getPoints();
        if (!p0 || !p1)
          continue;
        Tddd edge = p1->X - p0->X;
        double theta = surface_geometry::edgeCurvatureAngle(
            edge, l->geom_curvature.k1, l->geom_curvature.k2,
            l->geom_curvature.PD1, l->geom_curvature.PD2);
        if (std::isfinite(theta) && theta_target > 1e-10 && theta > 1e-10) {
          double lr = std::log(theta / theta_target);
          th_dist_sum += lr * lr;
          ++th_count;
        }
      }
      double fs_target_dist    = (fs_count > 0) ? fs_dist_sum / fs_count : 0.0;
      double theta_target_dist = (th_count > 0) ? th_dist_sum / th_count : 0.0;

      // ----- 節点ループ: local CV (local grading) -----
      double local_cv_mean = 0.;
      {
        double cv_sum = 0.;
        int cv_count = 0;
        for (auto* p : net.getBoundaryPoints()) {
          if (!p) continue;
          cv_sum += p->localEdgeLengthCV();
          ++cv_count;
        }
        if (cv_count > 0)
          local_cv_mean = cv_sum / cv_count;
      }

      // ----- スコア合成 -----
      const double w = rs.quality_resolution_weight;
      double score = worst_ir
                   + 1. / (1. + local_cv_mean)
                   - w * (fs_target_dist + 0.5 * theta_target_dist);
      return score;
    };



    // computeHdLimit はパッチ版を runTrials 直前に定義（重複を避けるためここからは削除）

    // --- TrialResult: 1シナリオの試行結果（best パッチを保持） ---
    // reject_code: 0=success, 1=no valid trial, 2=hd exceeded, 3=score not improved,
    //              4=replacePatch failed (外側で設定)
    struct TrialResult {
      std::string ops;
      double score = -1.;
      double hd = std::numeric_limits<double>::infinity();
      bool valid = false;
      int reject_code = 1;                    // valid=false の時の reject reason
      std::shared_ptr<Network> patch;         // best AFTER パッチ
      std::shared_ptr<Network> source_patch;  // BEFORE 参照パッチ (patch0)
    };

    // --- runPatchOps: パッチ上で ops を実行（trial 専用） ---
    auto runPatchOps = [&](const std::string& ops, Network& net, networkLine* cl) -> bool {
      auto doFlip = [&]() {
        while (true) {
          networkLine* best_line = nullptr;
          int best_improvement = 0;
          double best_len = -1.;
          for (auto* pl : net.Lines) {
            if (pl->getBoundaryFaces().size() != 2)
              continue;
            auto [vb, va] = valenceDeviationScore(pl);
            int improvement = vb - va;
            if (improvement <= 0)
              continue;
            double len = pl->length();
            if (improvement > best_improvement || (improvement == best_improvement && len > best_len)) {
              best_line = pl;
              best_improvement = improvement;
              best_len = len;
            }
          }
          if (!best_line)
            break;
          best_line->Flip(true);
        }
      };

      // ops 開始前の border を記憶（操作でトポロジが変わっても初期状態を参照）
      const auto initial_border = net.getBorderPoints();

      auto doSmooth = [&]() {
        auto pos_of = [](const networkPoint* q) -> Tddd { return q->X; };
        // 反復回数を 3 → 6 に増やし、各点の局所平衡により近づける
        for (int si = 0; si < 20; ++si) {
          for (auto* p : net.getPoints()) {
            if (p->BorderQ())
              continue;
            Tddd V = qualitySmoothingVector(p, p->X, pos_of);
            if (!isFlat(p, feature_angle)) {
              // feature edge 上の点: crease 方向（PD2）にのみスムージングを許可
              if (p->geom_curvature.valid) {
                double pd2_norm = Norm(p->geom_curvature.PD2);
                if (pd2_norm < 1e-12)
                  continue;
                Tddd dir = p->geom_curvature.PD2 / pd2_norm;
                V = Dot(V, dir) * dir;
              } else {
                // curvature 無効時: 隣接面法線の Cross でフォールバック
                auto faces = p->getBoundaryFaces();
                std::vector<Tddd> crease_dirs;
                for (std::size_t i = 0; i < faces.size(); ++i)
                  for (std::size_t j = i + 1; j < faces.size(); ++j) {
                    double cos_angle = Dot(faces[i]->normal, faces[j]->normal);
                    if (cos_angle < std::cos(feature_angle)) {
                      Tddd c = Cross(faces[i]->normal, faces[j]->normal);
                      double cn = Norm(c);
                      if (cn > 1e-12)
                        crease_dirs.push_back(c / cn);
                    }
                  }
                if (crease_dirs.empty())
                  continue;
                if (crease_dirs.size() == 1) {
                  V = Dot(V, crease_dirs[0]) * crease_dirs[0];
                } else {
                  continue; // 角点 → 固定
                }
              }
            } else if (p->Neumann) {
              Tddd n = p->getNormalAreaAveraged();
              double nn = Norm(n);
              if (nn > 1e-12) {
                n = n / nn;
                V = V - Dot(V, n) * n;
              }
            }
            double vn = Norm(V);
            if (vn > 1e-12) {
              // 変位クランプを 0.3→0.5 倍 localEdgeLength に緩和
              double limit = 0.05 * localEdgeLength(p);
              p->setXSingle(p->X + std::min(vn, limit) * Normalize(V));
            }
          }
          // setDodecaPoints を除外した軽量版
          for (const auto& l : net.getLines()) l->setBoundsSingle();
          for (const auto& p : net.getPoints()) p->setFaces();
          for (const auto& f : net.getFaces()) f->setGeometricProperties(ToX(f->getPoints()));
        }
      };

      for (char op : ops) {
        if (op == 'F') {
          doFlip();
        } else if (op == 'C') {
          try {
            if (cl->getBoundaryFaces().size() != 2)
              return false;
            auto [cp0, cp1] = cl->getPoints();
            bool cp0_border = initial_border.count(cp0);
            bool cp1_border = initial_border.count(cp1);
            if (cp0_border && cp1_border)
              return false; // 両端が border → collapse 不可
            Tddd targetX = 0.5 * (cp0->X + cp1->X);
            if (cp0_border)
              targetX = cp0->X; // border 点は移動させない
            else if (cp1_border)
              targetX = cp1->X;
            if (!cl->Collapse(targetX))
              return false;
          } catch (...) {
            return false;
          }
        } else if (op == 'S') {
          // setDodecaPoints を除外した軽量版 setGeometricPropertiesForce
          // DodecaPoints はパッチ border 辺を通じてパッチ外面を参照するため OMP 並列で不安全
          for (const auto& l : net.getLines()) l->setBoundsSingle();
          for (const auto& p : net.getPoints()) p->setFaces();
          for (const auto& f : net.getFaces()) f->setGeometricProperties(ToX(f->getPoints()));
          net.computePrincipalCurvatures(false);
          doSmooth();
        } else if (op == 'P') {
          if (!cl->Split(cl->X_mid))
            return false;
        }
      }
      return true;
    };

    // --- runTrial: 事前コピー済みパッチ上で操作+評価（並列実行される） ---
    auto runTrial = [&](const std::string& ops, std::shared_ptr<Network> patch, networkLine* cl,
                        const std::vector<OldFaceData>& orig_faces,
                        const V_netFp& ref_faces) -> TrialResult {
      TrialResult result;
      result.ops = ops;
      if (!cl)
        return result;

      try {
      if (!runPatchOps(ops, *patch, cl))
        return result;

      // 外周保存チェック
      {
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
            return result;
        }
      }

      projectOntoOldSurface(*patch, orig_faces);

      // setDodecaPoints を除外した軽量版（パッチ border 辺がパッチ外面を参照する問題を回避）
      for (const auto& l : patch->getLines()) l->setBoundsSingle();
      for (const auto& p : patch->getPoints()) p->setFaces();
      for (const auto& f : patch->getFaces()) f->setGeometricProperties(ToX(f->getPoints()));
      patch->computePrincipalCurvatures();

      V_netFp modified_faces = patch->getBoundaryFaces();
      result.hd = std::max(DirectedHausdorffDistance(ref_faces, modified_faces, 3),
                           DirectedHausdorffDistance(modified_faces, ref_faces, 3));

      {
        auto midpointMaxDist = [](const V_netFp& node_faces, const V_netFp& target_faces) -> double {
          double max_d = 0.;
          std::unordered_set<networkLine*> lines;
          for (auto* f : node_faces)
            if (f)
              for (auto* l : f->getLines())
                lines.insert(l);
          for (auto* l : lines) {
            if (l->getBoundaryFaces().size() < 2)
              continue;
            double min_d = 1e+20;
            for (auto* f : target_faces) {
              if (!f)
                continue;
              auto [w0, w1, nearest, normal] = Nearest_(l->X_mid, ToX(f));
              min_d = std::min(min_d, Norm(l->X_mid - nearest));
            }
            if (min_d < 1e+20)
              max_d = std::max(max_d, min_d);
          }
          return max_d;
        };
        result.hd = std::max(result.hd, midpointMaxDist(ref_faces, modified_faces));
        result.hd = std::max(result.hd, midpointMaxDist(modified_faces, ref_faces));
      }

      result.score = patchQuality(*patch);
      result.valid = true;
      result.patch = patch;
      } catch (...) {
        // OMP 並列内で例外を伝播させない
      }
      return result;
    };



    // --- computeHdLimit: パッチ上で許容 Hausdorff 距離を計算 ---
    // 形状忠実度の上限を決める。2 つの制約の min を取る:
    //   (a) bbox 対角 × rs.quality_hd_diag_ratio  — パッチ全体のスケール制限
    //   (b) 曲率半径 R × rs.quality_hd_curv_ratio — 角の近傍ではより厳しく
    // 呼び出し前提: patch の geometric properties と curvature が最新であること
    auto computeHdLimit = [&rs](const Network& patch) -> double {
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
      return std::min(rs.quality_hd_diag_ratio * diag, rs.quality_hd_curv_ratio * R);
    };

    // --- runTrials: copy は逐次、操作+評価は並列 ---
    auto runTrials = [&](networkLine* target_line, const std::vector<std::string>& scenarios) -> TrialResult {
      ++rt_calls;
      auto t_ref_begin = clk::now();
      auto patch0 = std::make_shared<Network>("file_name_is_not_given", "ref");
      water.copyLocalPatch(*patch0, target_line, 2);
      patch0->computePrincipalCurvatures();
      V_netFp ref_faces = patch0->getBoundaryFaces();
      double score_before = patchQuality(*patch0);

      // copyLocalPatch 内部で既に setGeometricPropertiesForce 済み (Network.cpp:2889)、
      // computePrincipalCurvatures は geom_curvature のみ書くので face props は最新。
      double hd_limit = computeHdLimit(*patch0);
      // ops 適用前のパッチ状態は全シナリオで同一なので orig_faces は 1 回だけ計算して共有
      auto orig_faces_shared = collectOldFaceData(*patch0);
      t_rt_ref += sec_since(t_ref_begin);

      // 各シナリオ用のパッチをコピー
      auto t_copy_begin = clk::now();
      std::vector<std::shared_ptr<Network>> patches(scenarios.size());
      std::vector<networkLine*> copied_lines(scenarios.size(), nullptr);
      _Pragma("omp parallel for") for (std::size_t i = 0; i < scenarios.size(); ++i) {
        patches[i] = std::make_shared<Network>("file_name_is_not_given", scenarios[i]);
        copied_lines[i] = water.copyLocalPatch(*patches[i], target_line, 2);
      }
      t_rt_copy += sec_since(t_copy_begin);

      // 並列: 操作 + 評価
      auto t_run_begin = clk::now();
      std::vector<TrialResult> results(scenarios.size());
      _Pragma("omp parallel for") for (std::size_t i = 0; i < scenarios.size(); ++i)
          results[i] = runTrial(scenarios[i], patches[i], copied_lines[i], orig_faces_shared, ref_faces);
      t_rt_run += sec_since(t_run_begin);

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
      double best_score = score_before + rs.quality_score_improve_margin;
      int n_valid = 0, n_within_hd = 0, n_improved = 0;
      for (auto& r : results) {
        if (!r.valid)
          continue;
        ++n_valid;
        if (r.hd > hd_limit)
          continue;
        ++n_within_hd;
        if (r.score > best_score) {
          ++n_improved;
          best_score = r.score;
          best = std::move(r);
        }
      }
      // 失敗原因の分類 (valid=true になった場合は 0 に上書き)
      if (n_valid == 0)
        best.reject_code = 1;              // NO_VALID_TRIAL
      else if (n_within_hd == 0)
        best.reject_code = 2;              // HD_EXCEEDED
      else if (n_improved == 0)
        best.reject_code = 3;              // SCORE_NOT_IMPROVED
      if (best.valid)
        best.reject_code = 0;              // 後で replacePatch 失敗なら 4 に上書きする
      best.source_patch = patch0;
      t_rt_pick += sec_since(t_pick_begin);
      return best;
    };

    // ------ [5a] split パス ------
    // ループ構造: 1 回の iteration で最大 1 split（batch 内最初の成功で break）。
    // 外側の終了条件は単一: total_splits >= max_splits_per_step。
    // 他の break 条件は "候補なし" / "batch 空" / "成功 0" / "topology エラー"。
    auto t_split_begin = clk::now();
    if (surface_split) {
      int total_splits = 0;
      const int max_splits_per_step = rs.max_splits_per_step;
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);

      while (total_splits < max_splits_per_step) {
        auto t_needs_begin = clk::now();
        std::vector<networkLine*> candidates;
        candidates.reserve(dirty.size());
        for (auto* l : dirty) {
          if (line_alive(l) && needsSplit(l))
            candidates.emplace_back(l);
        }
        if (candidates.empty()) { t_s_needs += sec_since(t_needs_begin); break; }
        // CORNER 辺を優先 → 次に Dirichlet → 最後に Neumann。同種内は辺長の長い順
        std::ranges::sort(candidates, [](const networkLine* a, const networkLine* b) {
          int ka = a->CORNER ? 0 : (a->Dirichlet ? 1 : 2);
          int kb = b->CORNER ? 0 : (b->Dirichlet ? 1 : 2);
          if (ka != kb) return ka < kb;
          return a->length() > b->length();
        });
        auto batch = remesh_detail::collect_non_adjacent(candidates);
        t_s_needs += sec_since(t_needs_begin);
        if (batch.empty())
          break;

        // batch を走査し、最初に replacePatch が成功した辺で抜ける
        bool divided_any = false;
        for (auto* l : batch) {
          if (!line_alive(l) || !needsSplit(l))
            continue;
          auto len = l->length();

          // split の reason コード分類: 0=global_max, 1=theta, 2=dihedral, 3=obtuse
          int reason_code = 1;
          const char* split_reason = "theta";
          if (len > global_max_len) { reason_code = 0; split_reason = "global_max"; }
          else {
            auto v = edgeThetaVerdict(l);
            if (v == EdgeThetaVerdict::SplitCandidate)      { reason_code = 1; split_reason = "theta"; }
            else if (v == EdgeThetaVerdict::CurvatureInvalid){ reason_code = 2; split_reason = "dihedral"; }
            else if (hasObtuseAngle(l))                      { reason_code = 3; split_reason = "obtuse"; }
          }

          int aid = attempt_counter++;
          auto t_rt_begin = clk::now();
          auto trial = runTrials(l, rs.split_scenarios);
          t_s_runtrials += sec_since(t_rt_begin);

          int rej = trial.valid ? 0 : trial.reject_code;
          collectTriggerEdge(l, /*op=*/0, reason_code, aid, /*success=*/0, rej);
          if (trial.source_patch)
            collectCandidatePatch(*trial.source_patch, /*op=*/0, reason_code, aid, /*success=*/0, rej);

          if (!trial.valid)
            continue;

          std::cout << Red << "time_step " << time_step << ": split reason=" << split_reason
                    << " len=" << len << " best_ops=" << trial.ops << colorReset << std::endl;
          auto t_rep_begin = clk::now();
          auto trial_tris_snapshot = snapshotPatchTris(*trial.patch);
          bool replaced = water.replacePatch(*trial.patch);
          t_s_replace += sec_since(t_rep_begin);
          if (replaced) {
            commitRemeshedSnapshot(trial_tris_snapshot, /*op=*/0, reason_code, aid);
            markAttemptSuccess(aid);
            divided_any = true;
            ++total_splits;
            break;
          } else {
            markAttemptReplaceFailed(aid);
          }
        }
        if (!divided_any)
          break;

        // 後処理: 幾何再計算、topology 検査、dirty set リフレッシュ
        auto t_post_begin = clk::now();
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        if (curvature_remesh_enabled)
          water.computePrincipalCurvatures(false);
        bool split_topo_error = false;
        for (const auto& l : water.getLines())
          if (!l->checkTopology()) { split_topo_error = true; break; }
        t_s_postop += sec_since(t_post_begin);
        if (split_topo_error) {
          std::cerr << Red << "[remesh] time_step " << time_step
                    << ": topology error after split batch, stopping further splits"
                    << colorReset << std::endl;
          break;
        }

        dirty.clear();
        for (auto* ll : water.getBoundaryLines())
          dirty.insert(ll);
      }
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

    // [5b] 削除: flipIfBatched は2次曲面投影なしで X_mid がずれるため廃止。
    // flip はパッチ内の ops（PF, PFS, CF, CFS 等）で実行され、projectOntoOldSurface で投影される。
    auto t_curv_begin = clk::now();
    if (curvature_remesh_enabled)
      water.computePrincipalCurvatures(false);
    t_curv += sec_since(t_curv_begin);

    // ------ [5c] collapse パス ------
    auto t_collapse_begin = clk::now();
    // ------ [5c] collapse パス ------
    // split と同じ構造: 1 iter = 最大 1 collapse。外側終了条件は total_collapses >= max。
    // rejected_collapse: 試行失敗 or 真の collapse ではなかった (smooth のみ) 辺を記憶して再試行しない。
    std::unordered_set<networkLine*> rejected_collapse;
    if (surface_collapse) {
      int total_collapses = 0;
      const int max_collapses_per_step = rs.max_collapses_per_step;
      // ログ用 theta 値取得
      auto getEdgeThetaForLog = [](const networkLine* l) -> double {
        if (!l->geom_curvature.valid)
          return -1.0;
        auto [p0, p1] = l->getPoints();
        Tddd edge = p1->X - p0->X;
        return surface_geometry::edgeCurvatureAngle(edge, l->geom_curvature.k1, l->geom_curvature.k2, l->geom_curvature.PD1, l->geom_curvature.PD2);
      };

      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);

      while (total_collapses < max_collapses_per_step) {
        auto t_needs_begin = clk::now();
        std::vector<networkLine*> candidates;
        candidates.reserve(dirty.size());
        for (auto* l : dirty) {
          if (!line_alive(l) || rejected_collapse.count(l))
            continue;
          if (needsCollapse(l))
            candidates.emplace_back(l);
        }
        if (candidates.empty()) { t_c_needs += sec_since(t_needs_begin); break; }
        // Sort by length ascending — collapse shortest first
        std::ranges::sort(candidates, [](const auto* a, const auto* b) { return a->length() < b->length(); });
        auto batch = remesh_detail::collect_non_adjacent(candidates);
        t_c_needs += sec_since(t_needs_begin);
        if (batch.empty())
          break;

        bool changed_in_batch = false;
        for (auto* l : batch) {
          if (!line_alive(l) || !needsCollapse(l))
            continue;
          auto [p0_orig, p1_orig] = l->getPoints();
          if (!p0_orig || !p1_orig)
            continue;

          double len = l->length();
          double theta = getEdgeThetaForLog(l);

          // collapse の reason コード分類:
          // 0=global_min, 1=theta+grading, 2=dihedral+grading, 3=degenerate, 4=obtuse
          int reason_code = 1;
          if (len < global_min_len) reason_code = 0;
          else if (hasObtuseAngle(l)) reason_code = 4;
          else {
            auto v = edgeThetaVerdict(l);
            if (v == EdgeThetaVerdict::CollapseCandidate) reason_code = 1;
            else if (v == EdgeThetaVerdict::CurvatureInvalid) reason_code = 2;
            else reason_code = 3;
          }

          int aid = attempt_counter++;
          auto t_rt_begin = clk::now();
          auto trial = runTrials(l, rs.collapse_scenarios);
          t_c_runtrials += sec_since(t_rt_begin);

          int rej = trial.valid ? 0 : trial.reject_code;
          collectTriggerEdge(l, /*op=*/1, reason_code, aid, /*success=*/0, rej);
          if (trial.source_patch)
            collectCandidatePatch(*trial.source_patch, /*op=*/1, reason_code, aid, /*success=*/0, rej);

          if (!trial.valid) {
            rejected_collapse.insert(l);
            continue;
          }

          std::cout << "time_step " << time_step << ": collapse len=" << len
                    << " theta=" << theta << " best_ops=" << trial.ops << std::endl;

          auto t_rep_begin = clk::now();
          auto trial_tris_snapshot = snapshotPatchTris(*trial.patch);
          bool replaced = water.replacePatch(*trial.patch);
          t_c_replace += sec_since(t_rep_begin);
          if (replaced) {
            commitRemeshedSnapshot(trial_tris_snapshot, /*op=*/1, reason_code, aid);
            markAttemptSuccess(aid);
            changed_in_batch = true;
            ++total_collapses;
            // trial.ops に 'C' が含まれない = smooth のみで collapse されなかった辺は再試行対象外
            if (trial.ops.find('C') == std::string::npos)
              rejected_collapse.insert(l);
            break;
          } else {
            markAttemptReplaceFailed(aid);
            rejected_collapse.insert(l);
          }
        }
        if (!changed_in_batch)
          break;

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


    // smoothing 共通: 面の最長辺を target line として返す
    auto longestEdgeOf = [](networkFace* f) -> networkLine* {
      networkLine* best = nullptr;
      double best_len = 0.;
      for (auto* l : f->getLines()) {
        if (!l) continue;
        double len = l->length();
        if (len > best_len) { best_len = len; best = l; }
      }
      return best;
    };

    auto runSmoothingPass = [&](const std::string& tag, std::function<bool(networkFace*)> face_predicate) {
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
      auto batch = remesh_detail::collect_non_adjacent(candidates);
      int total = 0;
      const int max_per_step = rs.max_smoothing_per_step;
      for (auto* l : batch) {
        if (!line_alive(l))
          continue;
        int aid = attempt_counter++;
        if (rs.smoothing_scenarios.empty())
          continue; // smoothing_scenarios が空なら smoothing pass スキップ
        auto trial = runTrials(l, rs.smoothing_scenarios);
        int rej = trial.valid ? 0 : trial.reject_code;
        collectTriggerEdge(l, /*op=*/2, /*reason=*/0, aid, /*success=*/0, rej);
        if (trial.source_patch)
          collectCandidatePatch(*trial.source_patch, /*op=*/2, /*reason=*/0, aid, /*success=*/0, rej);
        if (!trial.valid)
          continue;
        std::cout << Magenta << "time_step " << time_step << ": smooth(" << tag << ") best_ops=" << trial.ops << colorReset << std::endl;
        auto trial_tris_snapshot = snapshotPatchTris(*trial.patch);
        if (water.replacePatch(*trial.patch)) {
          commitRemeshedSnapshot(trial_tris_snapshot, /*op=*/2, /*reason=*/0, aid);
          markAttemptSuccess(aid);
          ++total;
          if (total >= max_per_step) break;
        } else {
          markAttemptReplaceFailed(aid);
        }
      }
      if (total > 0) {
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        if (curvature_remesh_enabled)
          water.computePrincipalCurvatures(false);
      }
    };

    // InradiusToCircumradius が低い面を対象: FS / FSFS / FSFSFS のみ試行
    auto t_smooth_ir_begin = clk::now();
    if (surface_smoothing) {
      constexpr double ir_threshold = 0.5;
      runSmoothingPass("ir", [](networkFace* f) {
        auto [p0, p1, p2] = f->getPoints();
        return InradiusToCircumradius(p0->X, p1->X, p2->X) < ir_threshold;
      });
    }
    t_smooth_ir += sec_since(t_smooth_ir_begin);

    // 60° からのズレが大きい面を対象: FS / FSFS / FSFSFS のみ試行
    auto t_smooth_angle_begin = clk::now();
    if (surface_smoothing) {
      constexpr double angle_dev_threshold = M_PI / 6.; // 30°
      runSmoothingPass("angle", [&](networkFace* f) {
        auto [p0, p1, p2] = f->getPoints();
        return maxAngleDeviation(p0->X, p1->X, p2->X) > angle_dev_threshold;
      });
    }
    t_smooth_angle += sec_since(t_smooth_angle_begin);

    // CORNER 近傍の面を強制対象: 水線(CORNER edge)に接する面、およびその 1-ring 隣接面。
    // CORNER trigger の patch は中心点の制約が厳しく F/S で quality score が改善しにくいため、
    // 独立 smoothing pass で周辺 edge を広めに拾って flip+smooth を適用する。
    if (surface_smoothing) {
      std::unordered_set<networkFace*> corner_adjacent;
      // 1-hop: CORNER edge を持つ face
      for (auto* l : water.getBoundaryLines()) {
        if (!l || !l->CORNER) continue;
        for (auto* f : l->getBoundaryFaces())
          if (f) corner_adjacent.insert(f);
      }
      // 2-hop: さらに隣接する face も対象に
      std::unordered_set<networkFace*> corner_adjacent_2ring = corner_adjacent;
      for (auto* f : corner_adjacent) {
        for (auto* l : f->getLines())
          if (l)
            for (auto* f2 : l->getBoundaryFaces())
              if (f2) corner_adjacent_2ring.insert(f2);
      }
      runSmoothingPass("corner_halo", [&](networkFace* f) {
        return corner_adjacent_2ring.count(f) > 0;
      });
    }


    // [5d] 削除: flipIfBatched は廃止（5b と同じ理由）。
    auto t_curv2_begin = clk::now();
    if (curvature_remesh_enabled)
      water.computePrincipalCurvatures(false);
    t_curv += sec_since(t_curv2_begin);

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
  }

  // [profile] セクション別累積時間サマリ
  {
    double t_total = t_split + t_collapse + t_smooth_ir + t_smooth_angle + t_topo + t_curv;
    std::cout << "[remesh profile] step=" << time_step
              << " split=" << t_split << "s"
              << " collapse=" << t_collapse << "s"
              << " smooth_ir=" << t_smooth_ir << "s"
              << " smooth_angle=" << t_smooth_angle << "s"
              << " curv=" << t_curv << "s"
              << " topo=" << t_topo << "s"
              << " total=" << t_total << "s" << std::endl;
    double t_rt_total = t_rt_ref + t_rt_copy + t_rt_run + t_rt_pick;
    std::cout << "[remesh profile] runTrials calls=" << rt_calls
              << " ref=" << t_rt_ref << "s"
              << " copy=" << t_rt_copy << "s"
              << " run=" << t_rt_run << "s"
              << " pick=" << t_rt_pick << "s"
              << " total=" << t_rt_total << "s" << std::endl;
    std::cout << "[remesh profile] split parts:"
              << " needs=" << t_s_needs << "s"
              << " runTrials=" << t_s_runtrials << "s"
              << " replace=" << t_s_replace << "s"
              << " postop=" << t_s_postop << "s" << std::endl;
    std::cout << "[remesh profile] collapse parts:"
              << " needs=" << t_c_needs << "s"
              << " runTrials=" << t_c_runtrials << "s"
              << " replace=" << t_c_replace << "s"
              << " postop=" << t_c_postop << "s" << std::endl;
  }

  // ========================================================================
  // [6] 事後チェック: folded face、tiny face
  // ========================================================================

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

  // ========================================================================
  // [7] 四面体化 + 幾何プロパティ更新
  // ========================================================================

#ifdef USE_TETGEN
  if (tetrahedralize)
    water.tetrahedralize();
#endif
  water.setGeometricPropertiesForce();
  water.improveTetrahedraDelaunay();

  // パッチ VTU をまとめて出力 (3 レイヤー)
  flushCandidatePatchVTU();
  flushRemeshedPatchVTU();
  flushTriggerEdgesVTU();

  // reject_code / bc_type 内訳サマリ (success=0 のみ集計)
  {
    std::array<int, 3> bc_count{0, 0, 0};                 // 0=Neumann, 1=Dirichlet, 2=CORNER
    std::array<std::array<int, 5>, 3> reject_by_bc{};     // reject_by_bc[bc][reject_code]
    std::array<int, 3> succ_by_bc{0, 0, 0};
    for (const auto& r : trigger_edge_records) {
      if (r.bc_type < 0 || r.bc_type > 2)
        continue;
      ++bc_count[r.bc_type];
      if (r.success == 1)
        ++succ_by_bc[r.bc_type];
      else if (r.reject_code >= 0 && r.reject_code <= 4)
        ++reject_by_bc[r.bc_type][r.reject_code];
    }
    const char* bc_names[3] = {"Neumann", "Dirichlet", "CORNER"};
    const char* rc_names[5] = {"success", "no_valid", "hd_exceeded", "score_no_improve", "replace_failed"};
    std::cout << "[remesh reject-summary] (trigger_edges, per bc_type)" << std::endl;
    for (int bc = 0; bc < 3; ++bc) {
      if (bc_count[bc] == 0) continue;
      std::cout << "  " << bc_names[bc] << " total=" << bc_count[bc]
                << " success=" << succ_by_bc[bc];
      for (int rc = 1; rc <= 4; ++rc)
        std::cout << " " << rc_names[rc] << "=" << reject_by_bc[bc][rc];
      std::cout << std::endl;
    }
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
    if (p0->CORNER && !p1->CORNER)
      targetX = p0->X;
    else if (p1->CORNER && !p0->CORNER)
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
  if (l->CORNER)
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

// valenceDeviationScore: flip 前後の頂点次数バランスを評価
// 返り値: {v_before, v_after} — 理想次数 s_mean からの偏差の二乗和（整数）
// v_after < v_before ならトポロジ的に改善
std::pair<int, int> valenceDeviationScore(const networkLine* l, int s_mean) {
  constexpr int INF = std::numeric_limits<int>::max();
  if (!l)
    return {INF, INF};
  auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1])
    return {INF, INF};

  auto [p0, p1] = l->getPoints();
  if (!p0 || !p1)
    return {INF, INF};

  auto* p2 = faces[0]->getPointOpposite(l);
  auto* p3 = faces[1]->getPointOpposite(l);
  if (!p2 || !p3 || p2 == p3)
    return {INF, INF};

  int s0 = p0->getLines().size();
  int s1 = p1->getLines().size();
  int s2 = p2->getLines().size();
  int s3 = p3->getLines().size();

  // flip 後 p0, p1 の次数が -1 されるので、3以下なら flip 不可
  if (s0 <= 3 || s1 <= 3)
    return {INF, INF};

  auto sq = [](int x) { return x * x; };
  int v_before = sq(s0 - s_mean) + sq(s1 - s_mean) + sq(s2 - s_mean) + sq(s3 - s_mean);
  int v_after = sq(s0 - 1 - s_mean) + sq(s1 - 1 - s_mean) + sq(s2 + 1 - s_mean) + sq(s3 + 1 - s_mean);
  return {v_before, v_after};
}

bool shouldFlipByQuality(const networkLine* l) {
  if (!l || l->CORNER)
    return false;
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
