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

void remesh_for_main_loop(Network& water, int time_step, double /*min_edge_length*/, bool tetrahedralize, bool surface_flip,
                          const CollisionSettings& collision_settings,
                          bool surface_split, bool surface_collapse,
                          bool surface_smoothing,
                          bool skip_post_remesh_quality_rejects,
                          const std::string& patch_output_directory,
                          double simulation_time,
                          PVDWriter* patch_pvd) {

  // パッチ VTU デバッグ出力: 全パッチを蓄積し、time_step 終了時にまとめて出力
  struct PatchFaceRecord {
    T3Tddd tri;   // 3頂点の座標
    double stage; // 0.0 = split, 1.0 = collapse
  };
  std::vector<PatchFaceRecord> patch_records;

  auto collectPatchFaces = [&](const Network& patch, double stage_val) {
    for (auto* f : patch.getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      patch_records.push_back({{p0->X, p1->X, p2->X}, stage_val});
    }
  };

  auto flushPatchVTU = [&]() {
    if (!patch_pvd && patch_output_directory.empty())
      return;
    const std::string vtu_name = water.getName() + "_patch_" + std::to_string(time_step) + ".vtu";
    const std::string filename = patch_output_directory.empty()
                                     ? "/tmp/" + vtu_name
                                     : patch_output_directory + "/" + vtu_name;
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp)
      return;
    int n_cells = patch_records.size();
    int n_points = n_cells * 3;
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", n_cells, n_points);
    // Points
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' format='ascii'>\n");
    for (auto& r : patch_records) {
      auto& [a, b, c] = r.tri;
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(a), (float)std::get<1>(a), (float)std::get<2>(a));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(b), (float)std::get<1>(b), (float)std::get<2>(b));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(c), (float)std::get<1>(c), (float)std::get<2>(c));
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    // Cells
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i)
      fprintf(fp, "%d %d %d\n", i * 3, i * 3 + 1, i * 3 + 2);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i)
      fprintf(fp, "%d\n", (i + 1) * 3);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i)
      fprintf(fp, "5\n"); // VTK_TRIANGLE
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    // CellData: stage (0=split, 1=collapse)
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type='Float32' Name='stage' NumberOfComponents='1' format='ascii'>\n");
    for (auto& r : patch_records)
      fprintf(fp, "%f\n", (float)r.stage);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    std::cout << "[patch debug] wrote " << filename << " (" << n_cells << " faces)" << std::endl;
    if (patch_pvd) {
      patch_pvd->push(vtu_name, simulation_time);
      patch_pvd->output();
    }
  };

  // ========================================================================
  // [1] パラメータ設定
  // ========================================================================

  const double rad = M_PI / 180.0;
  const double global_mean_len = Mean(extLength(water.getLines()));
  const double global_min_len = global_mean_len * 0.3; // これ以下の辺は存在させない
  const double global_max_len = global_mean_len * 3.0; // これ以上の辺は作らない
  constexpr double min_local_face_area_ratio = 0.05;

  // 曲面忠実度（θ）: theta_target = 2π/N で円筒 N 分割が目標
  constexpr bool curvature_remesh_enabled = true;                   // TODO: move to settings
  constexpr double theta_target = 2.0 * M_PI / 40.0;                // 円筒40分割
  constexpr double a_theta = 0.25;                                  // ヒステリシス帯（a < 1/3 で振動回避）
  constexpr double theta_split = (1.0 + a_theta) * theta_target;    // 32分割以下で split
  constexpr double theta_collapse = (1.0 - a_theta) * theta_target; // 53分割以上で collapse

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
  // [TEST] copyLocalPatch のテスト（初回のみ実行）
  // ========================================================================
  {
    auto* test_line = *water.getBoundaryLines().begin();
    if (test_line) {
      Network patch("file_name_is_not_given", "test_patch");
      auto* copied_line = water.copyLocalPatch(patch, test_line, 2);
      std::cout << "[copyLocalPatch test] original: Points=" << water.getPoints().size()
                << " Lines=" << water.getLines().size()
                << " Faces=" << water.getFaces().size() << std::endl;
      std::cout << "[copyLocalPatch test] patch:    Points=" << patch.getPoints().size()
                << " Lines=" << patch.getLines().size()
                << " Faces=" << patch.getFaces().size() << std::endl;
      if (copied_line) {
        std::cout << "[copyLocalPatch test] copied_line length=" << copied_line->length() << std::endl;
        // パッチ外周点（面が1つしかない辺の端点）
        std::unordered_set<networkPoint*> patch_boundary;
        for (auto* ll : patch.getLines()) {
          if (ll->getBoundaryFaces().size() < 2) {
            auto [pp0, pp1] = ll->getPoints();
            if (pp0)
              patch_boundary.insert(pp0);
            if (pp1)
              patch_boundary.insert(pp1);
          }
        }
        std::cout << "[copyLocalPatch test] patch boundary=" << patch_boundary.size()
                  << " internal=" << (patch.getPoints().size() - patch_boundary.size()) << std::endl;
        // 品質スコアのテスト
        double worst = 1.0;
        for (auto* f : patch.getFaces()) {
          if (!f->BoundaryQ())
            continue;
          auto [p0, p1, p2] = f->getPoints();
          double q = InradiusToCircumradius(p0->X, p1->X, p2->X);
          worst = std::min(worst, q);
        }
        std::cout << "[copyLocalPatch test] worst InradiusToCircumradius=" << worst << std::endl;
        // Collapse テスト（commonPoints チェック）
        auto [tp0, tp1] = copied_line->getPoints();
        auto commonPts = Intersection(tp0->getNeighborPointsOnSurfaces(), tp1->getNeighborPointsOnSurfaces());
        std::cout << "[copyLocalPatch test] commonPoints=" << commonPts.size() << " (should be 2)" << std::endl;
      } else {
        std::cout << "[copyLocalPatch test] ERROR: copied_line is null" << std::endl;
      }
    }
  }

  // ========================================================================
  // [5] split / flip / collapse ループ（iter_divide_collapse 回繰り返し）
  // ========================================================================

  const int iter_divide_collapse = 10;
  for (auto i = 0; i < iter_divide_collapse; i++) {
    auto line_alive = [&](const networkLine* l) { return l && (water.Lines.find(const_cast<networkLine*>(l)) != water.Lines.end()); };

    // feature 保護の閾値（30° — これ以上の法線変化は feature として保護）
    constexpr double feature_angle = 30.0 * M_PI / 180.0;

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
    auto needsSplit = [&](networkLine* l) -> bool {
      auto len = l->length();

      // global 上限を超える長辺は無条件で split 候補
      if (len > global_max_len)
        return true;

      auto verdict = edgeThetaVerdict(l);

      // 曲面忠実度
      if (verdict == EdgeThetaVerdict::SplitCandidate)
        return true;
      if (verdict == EdgeThetaVerdict::CollapseCandidate)
        return false;

      // CurvatureInvalid → 二面角 fallback（非 CORNER のみ）
      if (verdict == EdgeThetaVerdict::CurvatureInvalid && !l->CORNER) {
        auto surfaces = l->getBoundaryFaces();
        if (surfaces.size() == 2 && surfaces[0] && surfaces[1]) {
          double alpha = std::acos(std::clamp(
              Dot(surfaces[0]->normal, surfaces[1]->normal), -1.0, 1.0));
          if (alpha >= theta_split)
            return true;
        }
      }
      return false;
    };

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
        auto [t0, t1, near_linear, normal] = Nearest_(X, old_faces[i].tri);
        // 2. Newton 反復で2次曲面上の正確な t0t1
        for (int iter = 0; iter < 5; ++iter) {
          Tddd Xq = evalQuad(t0, t1, old_faces[i]);
          Tddd diff = Xq - X;
          Tddd dXdt0 = evalQuadDeriv(t0, t1, 1, 0, old_faces[i]);
          Tddd dXdt1 = evalQuadDeriv(t0, t1, 0, 1, old_faces[i]);
          double g0 = Dot(diff, dXdt0), g1 = Dot(diff, dXdt1);
          if (std::abs(g0) < 1e-10 && std::abs(g1) < 1e-10)
            break;
          Tddd d2Xdt00 = evalQuadDeriv(t0, t1, 2, 0, old_faces[i]);
          Tddd d2Xdt11 = evalQuadDeriv(t0, t1, 0, 2, old_faces[i]);
          Tddd d2Xdt01 = evalQuadDeriv(t0, t1, 1, 1, old_faces[i]);
          double H00 = Dot(dXdt0, dXdt0) + Dot(diff, d2Xdt00);
          double H11 = Dot(dXdt1, dXdt1) + Dot(diff, d2Xdt11);
          double H01 = Dot(dXdt0, dXdt1) + Dot(diff, d2Xdt01);
          double det = H00 * H11 - H01 * H01;
          if (std::abs(det) < 1e-20)
            break;
          t0 += -(H11 * g0 - H01 * g1) / det;
          t1 += -(H00 * g1 - H01 * g0) / det;
          t0 = std::max(0.0, t0);
          t1 = std::max(0.0, t1);
          if (t0 + t1 > 1.0) {
            double s = t0 + t1;
            t0 /= s;
            t1 /= s;
          }
        }
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

    auto worstScore = [](const Network& net) {
      double worst = 1.0;
      for (auto* f : net.getBoundaryFaces()) {
        auto [p0, p1, p2] = f->getPoints();
        double q = InradiusToCircumradius(p0->X, p1->X, p2->X);
        worst = std::min(worst, q);
      }
      // 反転チェック: 隣接面の法線が反転していれば -1 を返す
      for (auto* l : net.getLines()) {
        auto faces = l->getBoundaryFaces();
        if (faces.size() == 2 && faces[0] && faces[1]) {
          if (Dot(faces[0]->normal, faces[1]->normal) < -0.3)
            return -1.0;
        }
      }
      return worst;
    };

    auto patchBoundaryPoints = [](const Network& net) {
      std::unordered_set<networkPoint*> bpts;
      for (auto* l : net.getLines()) {
        if (l->getBoundaryFaces().size() < 2) {
          auto [p0, p1] = l->getPoints();
          if (p0)
            bpts.insert(p0);
          if (p1)
            bpts.insert(p1);
        }
      }
      return bpts;
    };

    auto computeHdLimit = [](networkLine* l) -> double {
      double kmax_sum = 0.;
      int kmax_count = 0;
      std::unordered_set<networkPoint*> face_verts;
      for (auto* f : l->getFaces()) {
        auto [fp0, fp1, fp2] = f->getPoints();
        face_verts.insert(fp0);
        face_verts.insert(fp1);
        face_verts.insert(fp2);
      }
      for (auto* p : face_verts) {
        if (p->geom_curvature.valid) {
          kmax_sum += p->geom_curvature.kmax;
          ++kmax_count;
        }
      }
      double kmax_avg = (kmax_count > 0) ? (kmax_sum / kmax_count) : 0.;
      double curvature_radius_R = (kmax_avg > 1e-10) ? (1.0 / kmax_avg) : std::numeric_limits<double>::infinity();
      return std::min(0.1 * curvature_radius_R, 0.1 * l->length());
    };

    // --- TrialResult: 1シナリオの試行結果（best パッチを保持） ---
    struct TrialResult {
      std::string ops;
      double score = -1.;
      double hd = std::numeric_limits<double>::infinity();
      bool valid = false;
      std::shared_ptr<Network> patch; // best パッチを保持
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

      auto doSmooth = [&]() {
        auto bset = patchBoundaryPoints(net);
        auto pos_of = [](const networkPoint* q) -> Tddd { return q->X; };
        for (int si = 0; si < 3; ++si) {
          for (auto* p : net.getPoints()) {
            if (bset.count(p))
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
              double limit = 0.3 * localEdgeLength(p);
              p->setXSingle(p->X + std::min(vn, limit) * Normalize(V));
            }
          }
          net.setGeometricPropertiesForce();
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
            if (!cl->Collapse(0.5 * (cp0->X + cp1->X)))
              return false;
          } catch (...) {
            return false;
          }
        } else if (op == 'S') {
          net.setGeometricPropertiesForce();
          net.computePrincipalCurvatures(false);
          doSmooth();
        } else if (op == 'P') {
          if (!cl->Split())
            return false;
        }
      }
      return true;
    };

    // --- runTrial: パッチ上で1シナリオを試行、パッチを保持して返す ---
    auto runTrial = [&](const std::string& ops, networkLine* target_line,
                        const V_netFp& ref_faces) -> TrialResult {
      TrialResult result;
      result.ops = ops;
      auto patch = std::make_shared<Network>("file_name_is_not_given", ops);
      auto* cl = water.copyLocalPatch(*patch, target_line, 1);
      if (!cl)
        return result;
      patch->computePrincipalCurvatures(true);

      // ops 実行前の patch 面データを保存
      auto orig_faces = collectOldFaceData(*patch);

      if (!runPatchOps(ops, *patch, cl))
        return result;

      // 頂点→辺節点の順に投影 + phi/phi_t 補間
      projectOntoOldSurface(*patch, orig_faces);

      patch->setGeometricPropertiesForce();
      V_netFp modified_faces = patch->getBoundaryFaces();
      result.hd = std::max(DirectedHausdorffDistance(ref_faces, modified_faces, 3),
                           DirectedHausdorffDistance(modified_faces, ref_faces, 3));

      // 辺中点 X_mid を面上サンプルとは独立に HD 評価に加える。
      //
      // 面上の均等サンプル点による HD では、flip で頂点位置が変わらない場合に
      // HD ≈ 0 となり検出できない。しかし flip で面の接続が変わると、
      // 面の頂点は元の曲面に接しているが辺の中点が曲面から浮くケースがある。
      // 特に CORNER（構造物の角）を跨ぐ flip では、新辺の X_mid が角を跨いだ
      // 位置になり、元の曲面のどの面からも離れる。
      //
      // 各 X_mid を独立なサンプル点として target 面集合への最小距離を求め、
      // その最大を HD に加算することで、このような異常パッチを却下できる。
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
      result.score = worstScore(*patch);
      result.valid = true;
      result.patch = patch; // パッチを保持
      return result;
    };

    // --- runTrials: 複数シナリオを並列実行し、best の TrialResult を返す ---
    auto runTrials = [&](networkLine* target_line,
                         const std::vector<std::string>& scenarios,
                         double hd_limit) -> TrialResult {
      Network patch0("file_name_is_not_given", "ref");
      water.copyLocalPatch(patch0, target_line, 1);
      V_netFp ref_faces = patch0.getBoundaryFaces();
      double score_before = worstScore(patch0);

      std::vector<TrialResult> results(scenarios.size());
#pragma omp parallel for schedule(dynamic)
      for (std::size_t i = 0; i < scenarios.size(); ++i)
        results[i] = runTrial(scenarios[i], target_line, ref_faces);

      TrialResult best;
      double best_score = score_before;
      for (auto& r : results) {
        if (!r.valid || r.hd > hd_limit)
          continue;
        if (r.score > best_score && r.score >= 0.05) {
          best_score = r.score;
          best = std::move(r);
        }
      }
      return best; // .valid == false なら no-op
    };

    // ------ [5a] split パス ------
    if (surface_split) {
      int total_splits = 0;
      const int max_splits_per_step = std::clamp(static_cast<int>(water.getBoundaryLines().size()) / 50, 10, 10);
      std::unordered_set<networkLine*> dirty;
      for (auto* l : water.getBoundaryLines())
        dirty.insert(l);
      for (auto iter = 0; iter < 20 && total_splits < max_splits_per_step; iter++) {
        std::vector<networkLine*> candidates;
        candidates.reserve(dirty.size());
        for (auto* l : dirty) {
          if (!line_alive(l))
            continue;
          if (needsSplit(l))
            candidates.emplace_back(l);
        }
        if (candidates.empty())
          break;
        auto batch = remesh_detail::collect_non_adjacent(candidates);
        if (batch.empty())
          break;
        bool divided_any = false;
        bool split_topo_error = false;
        for (auto* l : batch) {
          if (!line_alive(l) || !needsSplit(l))
            continue;
          auto len = l->length();

          // split の理由をログ
          const char* split_reason = "unknown";
          if (len > global_max_len)
            split_reason = "global_max";
          else {
            auto v = edgeThetaVerdict(l);
            if (v == EdgeThetaVerdict::SplitCandidate)
              split_reason = "theta";
            else if (v == EdgeThetaVerdict::CurvatureInvalid)
              split_reason = "dihedral";
          }

          // trial 評価（並列）→ best パッチを replacePatch で適用
          const double hd_limit = computeHdLimit(l);
          auto trial = runTrials(l, {"P", "PF", "PS", "PFS", "PSF", "PSFS", "PSFSFS"}, hd_limit);

          if (!trial.valid)
            continue;

          std::cout << Red << "time_step " << time_step << ": split" << " reason=" << split_reason << " len=" << len << " best_ops=" << trial.ops << colorReset << std::endl;
          collectPatchFaces(*trial.patch, 0.0);
          if (water.replacePatch(*trial.patch)) {
            divided_any = true;
            ++total_splits;
            break; // メッシュが変わったのでバッチを捨てて再収集
          }
          if (total_splits >= max_splits_per_step)
            break;
        }
        if (!divided_any)
          break;
        water.setGeometricPropertiesForce();
        water.checkConnectivity();
        if (curvature_remesh_enabled)
          water.computePrincipalCurvatures(false);
        for (const auto& l : water.getLines())
          if (!l->checkTopology()) {
            split_topo_error = true;
            break;
          }
        if (split_topo_error) {
          std::cerr << Red << "[remesh] time_step " << time_step << ": topology error after split batch, stopping further splits" << colorReset << std::endl;
          break;
        }
        // replacePatch 後は dirty を全辺に再設定
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

    // [5b] 削除: flipIfBatched は2次曲面投影なしで X_mid がずれるため廃止。
    // flip はパッチ内の ops（PF, PFS, CF, CFS 等）で実行され、projectOntoOldSurface で投影される。
    if (curvature_remesh_enabled)
      water.computePrincipalCurvatures(false);

    // ------ [5c] collapse パス ------
    std::unordered_set<networkLine*> rejected_collapse;
    if (surface_collapse)
      for (auto iter = 0; iter < 20; iter++) {
        bool found_small_line = false;
        // tiny face flip and tiny face collapse are disabled
        // theta based split/collapse control should handle these cases
        // if (surface_flip) { ... tiny_flip ... }
        // for (auto tiny_iter ...) { ... tiny_collapse ... }

        // --- needsCollapse: この辺を collapse する理由があるか？ ---
        // 以下のいずれかに該当すれば collapse 候補（OR 条件）:
        //   A. 極短辺（dt 保護）
        //   B. 曲面忠実度 AND grading の両方が要求
        //   C. 品質トリガー: 隣接面が既に潰れている（既存の悪い面を修復）
        auto needsCollapse = [&](networkLine* l) -> bool {
          auto len = l->length();

          // A. 極端に短い辺は無条件で collapse（dt 保護）
          if (len < global_min_len)
            return true;

          // B. 曲面忠実度 AND grading
          {
            auto verdict = edgeThetaVerdict(l);

            // SplitCandidate → 曲率的に粗い → collapse しない
            if (verdict != EdgeThetaVerdict::SplitCandidate) {
              bool curvature_wants_collapse = false;
              if (verdict == EdgeThetaVerdict::CollapseCandidate) {
                curvature_wants_collapse = true;
              } else if (verdict == EdgeThetaVerdict::CurvatureInvalid) {
                auto surfaces = l->getBoundaryFaces();
                if (surfaces.size() == 2 && surfaces[0] && surfaces[1]) {
                  double alpha = std::acos(std::clamp(
                      Dot(surfaces[0]->normal, surfaces[1]->normal), -1.0, 1.0));
                  curvature_wants_collapse = (alpha < theta_collapse);
                }
              }

              bool grading_wants_collapse = false;
              auto local_mean = localEdgeLength(l);
              if (local_mean > 0. && len < grading_collapse * local_mean)
                grading_wants_collapse = true;

              if (curvature_wants_collapse && grading_wants_collapse)
                return true;
            }
          }

          // C. 品質トリガー: 隣接面が既に潰れている → 修復のために collapse
          {
            auto surfaces = l->getBoundaryFaces();
            if (surfaces.size() == 2 && surfaces[0] && surfaces[1]) {
              auto local_mean_len_val = localEdgeLength(l);
              if (local_mean_len_val > 0.) {
                for (auto* f : surfaces) {
                  auto pts = f->getPoints();
                  double area = TriangleArea(pts[0]->X, pts[1]->X, pts[2]->X);
                  double max_edge = std::max({Norm(pts[0]->X - pts[1]->X),
                                              Norm(pts[1]->X - pts[2]->X),
                                              Norm(pts[2]->X - pts[0]->X)});
                  double alt = (max_edge > 1e-20) ? 2.0 * area / max_edge : 0.0;
                  // 低高度: 隣接面の高度が local_mean の 3% 未満（本当に潰れた面のみ）
                  if (alt < 0.03 * local_mean_len_val)
                    return true;
                }
              }
            }
          }

          return false;
        };

        // ログ用 theta 値取得
        auto getEdgeThetaForLog = [](const networkLine* l) -> double {
          if (!l->geom_curvature.valid)
            return -1.0;
          auto [p0, p1] = l->getPoints();
          Tddd edge = p1->X - p0->X;
          return surface_geometry::edgeCurvatureAngle(
              edge, l->geom_curvature.k1, l->geom_curvature.k2,
              l->geom_curvature.PD1, l->geom_curvature.PD2);
        };

        std::unordered_set<networkLine*> dirty;
        for (auto* l : water.getBoundaryLines())
          dirty.insert(l);
        while (!dirty.empty()) {
          std::vector<networkLine*> candidates;
          candidates.reserve(dirty.size());
          for (auto* l : dirty) {
            if (!line_alive(l))
              continue;
            if (rejected_collapse.count(l))
              continue;
            if (needsCollapse(l))
              candidates.emplace_back(l);
          }
          if (candidates.empty())
            break;
          // Sort by length ascending — collapse shortest first
          std::ranges::sort(candidates, [](const auto* a, const auto* b) { return a->length() < b->length(); });
          auto batch = remesh_detail::collect_non_adjacent(candidates);
          if (batch.empty())
            break;
          bool changed_in_batch = false;
          for (auto* l : batch) {
            if (!line_alive(l) || !needsCollapse(l))
              continue;
            auto [p0_orig, p1_orig] = l->getPoints();
            if (!p0_orig || !p1_orig)
              continue;

            const double hd_limit = computeHdLimit(l);
            double len = l->length();
            double theta = getEdgeThetaForLog(l);

            auto trial = runTrials(l, {"C", "CS", "CF", "CFS", "FC", "FCS", "F", "FS", "CSFS", "CSFSFS"}, hd_limit);

            if (!trial.valid) {
              std::cout << Yellow << "time_step " << time_step
                        << ": collapse REJECTED len=" << len
                        << " hausdorff_distance_limit=" << hd_limit
                        << colorReset << std::endl;
              rejected_collapse.insert(l);
              continue;
            }

            std::cout << "time_step " << time_step
                      << ": collapse len=" << len
                      << " theta=" << theta
                      << " best_ops=" << trial.ops << std::endl;

            collectPatchFaces(*trial.patch, 1.0);
            if (water.replacePatch(*trial.patch)) {
              changed_in_batch = true;
              if (trial.ops.find('C') != std::string::npos)
                found_small_line = true;
              else
                rejected_collapse.insert(l);
              break; // メッシュが変わったのでバッチを捨てて再収集
            } else {
              rejected_collapse.insert(l);
            }
          }
          if (changed_in_batch) {
            water.setGeometricPropertiesForce();
            water.checkConnectivity();
            if (curvature_remesh_enabled)
              water.computePrincipalCurvatures(false);
          }
          if (!changed_in_batch)
            break;
          dirty.clear();
          for (auto* ll : water.getBoundaryLines())
            dirty.insert(ll);
        }
        if (!found_small_line)
          break;
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

    // [5d] 削除: flipIfBatched は廃止（5b と同じ理由）。
    if (curvature_remesh_enabled)
      water.computePrincipalCurvatures(false);

    // ------ [5e] topology チェック ------
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

  // パッチ VTU をまとめて出力
  flushPatchVTU();
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
