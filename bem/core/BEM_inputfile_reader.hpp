#pragma once

#include <filesystem>
#include <optional>
#include <regex>

#include "BEM_collision_settings.hpp"
#include "BEM_output_info.hpp"
#include "BEM_time_domain_types.hpp"
#include "Network.hpp"
#include "basic.hpp"
#include "TanakaSolitaryWave.hpp"

std::string toLowerCase(const std::string& str);

struct SimulationSettings {
  enum class DomainMode { Time,
                          Frequency };
  // Canonical settings file name (legacy: setting.json).
  const std::string settings_filename = "settings.json";
  const std::string legacy_setting_filename = "setting.json";

  DomainMode domain_mode = DomainMode::Time;

  // Raw JSON (kept for debugging/experimentation; should not be used for simulation logic in main.cpp).
  std::filesystem::path settings_file_path;
  JSON settingJSON;

  /* ----------------------------- common settings ---------------------------- */
  struct CommonSettings {
    std::filesystem::path input_directory;
    std::filesystem::path output_directory;
    // Also store physical constants (they are also applied to global variables for legacy code paths).
    double water_density = _WATER_DENSITY_;
    double gravity = _GRAVITY_;
  } common;

  /* --------------------------- time-domain settings ------------------------- */
  struct TimeDomainSettings {
    int end_time_step = 0;
    double end_time = 0.0;
    double max_dt = 0.0;
    int rk_order = 4; // Runge-Kutta order (1=Euler, 2=Heun, 3=RK3, 4=RK4)
    struct NodeRelocationSettings {
      using Method = ::NodeRelocationMethod;
      using Surface = ::NodeRelocationSurface;
      using MidpointMode = ::InterpolationMidpointMode;
      Method method = Method::none;
      int period = 0;
      Surface surface = Surface::pseudo_quadratic;
      MidpointMode midpoint_mode = MidpointMode::nearest;
      bool surface_explicitly_set = false;
      bool debug_output = false;
    } node_relocation;
  } time;

  /* ------------------------- frequency-domain settings ---------------------- */
  struct FrequencyDomainSettings {
    std::vector<double> omegas;
    std::vector<int> dofs;
  } frequency;

  /* ------------------------------ BEM settings ------------------------------ */
  struct BEMSettings {
    struct ElementSettings {
      bool linear = true;
      bool pseudo_quadratic = false;
      bool true_quadratic = false;
      bool quadratic_linear_hybrid = false;
    } element;

    struct SolverSettings {
      std::string solver_type = "LU";
      std::string preconditioner_type = "NONE";
      // ILU/MILU/ILUT sparsity / neighborhood settings
      std::string ilu_neighborhood_type = "BUCKETS"; // "BUCKETS" (default) or "K-RING"
      int ilu_kring_num = 1;                         // 0..20 (only used when type == "K-RING")
      // MILU params (used when preconditioner_type == "MILU")
      double milu_omega = 1.0; // diagonal compensation weight (MILU(0) -> omega=1)
      // ILUT params (used when preconditioner_type == "ILUT")
      double ilut_drop_tol = 1e-3;
      int ilut_max_entries_per_row = 50;
      double ilut_pivot_min = 1e-12;
      // Schwarz / block-Jacobi params (used when preconditioner_type == "SCHWARZ")
      int schwarz_core_k = 1;           // graph BFS depth for non-overlapping "cores"
      int schwarz_overlap_k = 1;        // additional overlap BFS depth
      int schwarz_max_core_size = 64;   // cap for core block size
      int schwarz_max_block_size = 128; // cap for full block size (core+overlap)
      double schwarz_pivot_min = 1e-12; // fallback diagonal clamp for singular/near-singular blocks
      double schwarz_diag_shift = 0.0;  // optional diagonal regularization added to each dense block
      std::string coupling_type = "NONE";
      double coupling_tol = 1e-10;
      std::vector<double> coupling_params;
      double solver_tol = 1e-9;
      int solver_max_iter = 500;
      int solver_restart = 100;
      // Metal M2L GPU acceleration settings (requires GMRES + FMM)
      bool use_metal_m2l = false;         // Enable Metal GPU acceleration for M2L
      bool metal_m2l_threadgroup = false; // true: use threadgroup parallelization
      bool metal_m2l_sort_terms = false;  // Sort terms for improved memory locality
      // Nearfield integration mode: "scalar" (default), "simd" (NEON float 4-target), "simd_double" (NEON double 2-target), "metal" (GPU)
      std::string nearfield_mode = "cell_scalar";
      // P2M quadrature: number of Dunavant points (1, 3, 6, 7). Default=1 for backward compatibility.
      int p2m_quadrature_points = 6;
      // MAC criterion parameter (smaller = stricter separation = faster but less accurate)
      double mac_theta = 0.25;
      // FMM tree structure settings (requires GMRES + FMM)
      int fmm_max_level = 7;          // Maximum tree depth
      int fmm_bucket_max_points = 50; // Bucket split threshold (number of points)
      // Pressure-based detachment
      bool enable_pressure_detachment = false;
      double detachment_pressure_threshold = 0.0;
      int detachment_consecutive_steps = 3;
    } solver;
  } bem;

  /* ---------------------------- remeshing settings -------------------------- */
  struct RemeshingSettings {
    // min_edge_length removed — limit_len is now always global_mean_len * 0.1
    double min_edge_length = 0.0;  // kept for backward compatibility but unused

    // Meshing options
    bool tetrahedralize = false;
    bool surface_flip = false;
    bool surface_split = false;
    bool surface_collapse = false;
    bool surface_smoothing = false;
    bool improve_tetrahedra = false;

    // Optional simulation-time windows for remeshing logic (currently consumed in main.cpp).
    double stop_remesh_time = 1E+10;
    double force_remesh_time = 0.0;
    int grid_refinement = 0;

    // Initial mesh pre-relaxation (no physical time advance)
    struct InitialMeshPreRelax {
      // Optional; controlled by settings.json "initial_mesh_pre_relax".
      bool enabled = true;
      int loop = 20;
      double coef = 0.005;
    } initial_mesh_pre_relax;

    // Surface collision detection and resolution
    // CollisionSettings is defined in BEM_collision.hpp (separate compilation unit)
    CollisionSettings collision;

    // SubsurfaceAltitudeRejectSettings is defined in BEM_collision_settings.hpp
    ::SubsurfaceAltitudeRejectSettings subsurface_altitude_reject;

    bool shell_visualization = false;
    bool front_advancing_debug = false;

    // Curvature-based mesh density control
    double max_edge_length = 0.0;  // 0 なら無効。曲率ベース h_target の上限に使う。

    // ---- Remesh numerical parameters (default = current hardcoded values) ----

    // A. 辺長制御 (len_*)
    double len_target_divisor = 40.0;       // free_surface_target_len = sqrt(dx^2+dy^2) / この値
    double len_fs_split_ratio = 1.1;        // 自由表面辺: len > fs_target * この値 で split
    double len_fs_collapse_ratio = 0.55;    // 自由表面辺: len < fs_target * この値 で collapse
    double len_global_max_ratio = 3.0;      // 全辺共通: len > fs_target * この値 で無条件 split
    double len_global_min_ratio = 0.2;      // 全辺共通: len < fs_target * この値 で無条件 collapse

    // B. 曲率忠実度 (theta_*)
    bool theta_enabled = true;              // 曲率ベース split/collapse の有効化
    double theta_target_N = 60.0;           // 円筒 N 分割目標 → theta_target = 2π/N
    double theta_split_ratio = 1.2;         // θ > theta_target * この値 で split
    double theta_collapse_ratio = 0.8;      // θ < theta_target * この値 で collapse

    // C. スムージング
    double feature_angle_deg = 30.0;        // feature edge 保護角度 [deg]

    // D. ループ制御
    //   各パスで試行する ops シナリオを直接指定できる。runTrials は各シナリオを
    //   並列に試行し、patchQuality が最大のものを採択する。
    //   ops の文字:  P=Split (trigger 辺), C=Collapse (trigger 辺), F=Flip (patch 内全 line valence 最適化),
    //              S=Smooth (patch 内全 point の品質向上変位)
    //   例:  split_scenarios = {"P","PS","PF","PFS","PFSFS","PFSFSFS"}
    //        smoothing_scenarios = {"FS","FSFS","FSFSFS"}
    std::vector<std::string> split_scenarios     = {"P","PS","PF","PFS","PFSFS","PFSFSFS"};
    std::vector<std::string> collapse_scenarios  = {"C","CS","CF","CFS","CFSFS","CFSFSFS"};
    std::vector<std::string> smoothing_scenarios = {"FS","FSFS","FSFSFS"};
    // 各パスの 1 ステップあたり成功上限
    int max_splits_per_step    = 10;
    int max_collapses_per_step = 10;
    int max_smoothing_per_step = 30;  // smoothing 独立パスの per-BC_predicate 上限（ir/angle/corner_halo 各々）
    // [5] split / smoothing / collapse ループを何回繰り返すか（外側ループ回数）
    int iter_split_collapse = 3;

    // リメッシュ方式の選択
    //   "trial"    : 既存。patchQuality ベースの local trial-and-error 方式（重いが物理保存強い）
    //   "explicit" : MeshLab 風 Isotropic Explicit。split→collapse→flip→smooth を直接適用（速い）
    std::string remesh_method = "trial";

    // E. Trial acceptance / quality gate (緩めれば採用されやすくなる)
    double quality_hd_diag_ratio = 0.1;        // hd_limit = min(diag * this, R * quality_hd_curv_ratio)
    double quality_hd_curv_ratio = 0.2;
    double quality_normal_flip_cos = -0.3;     // 隣接面法線内積 < この値 で反転判定（致命的不適合）
    double quality_score_improve_margin = 0.0; // trial.score > score_before + この値 で採用
    // F. patchQuality スコア設計
    //   score = worst_ir + 1/(1+local_cv_mean)
    //         - resolution_weight * (fs_target_dist + 0.5 * theta_target_dist)
    //   fs_target_dist = mean_{Dirichlet辺} (log(len/fs_target))^2   ← 長短両方向に対称
    //   theta_target_dist = mean_{曲率valid辺} (log(theta/theta_target))^2
    // resolution_weight を大きくすると target 達成が優先される (local quality より)
    double quality_resolution_weight = 1.0;
  } remeshing;

  /* ------------------------------- VPM settings ----------------------------- */
  struct VPMSettings {
    bool enabled = false;
    std::size_t wall_min_absorb_receivers = 1;
    double wall_min_absorb_total_weight = 1e-5;
    double sigma_factor = 1.0;
    std::string stretching_scheme = "transpose";
    std::string PSE_correction = "curvature";
  } vpm;

  /* ----------------------------- checkpoint settings ----------------------- */
  struct CheckpointSettings {
    int interval = 0;         // 0=disabled, N=save every N steps
    std::string restart_from; // checkpoint file path (empty=fresh start)
    int max_keep = 3;         // number of checkpoint files to keep (0=keep all)
  } checkpoint;

  /* ------------------------------- objects/data ----------------------------- */
  std::map<std::string, outputInfo> NetOutputInfo;
  std::vector<Network*> FluidObject, RigidBodyObject, SoftBodyObject, AbsorberObject;
  std::vector<JSON> MeasurementJSONs;

  SimulationSettings(std::filesystem::path input_directory, DomainMode mode_in = DomainMode::Time);
  /* -------------------------------------------------------------------------- */

  void setOutputInfo(auto output_name, const std::filesystem::path& output_directory) {
    std::cout << "setOutputInfo" << std::endl;
    this->NetOutputInfo[output_name].pvd_file_name = output_name;
    this->NetOutputInfo[output_name].vtu_file_name = output_name + "_";
    this->NetOutputInfo[output_name].PVD = new PVDWriter(output_directory / (output_name + ".pvd"));
  };

  /* -------------------------------------------------------------------------- */

  void setTypes(Network* net) {
    std::cout << "setTypes" << std::endl;
    //$ set type
    auto type = net->inputJSON.at("type")[0];
    if (type.contains("RigidBody")) {
      std::cout << "RigidBody" << std::endl;
      RigidBodyObject.emplace_back(net);
      net->isRigidBody = true;
      net->isSoftBody = net->isFluid = false;
    } else if (type.contains("SoftBody") || type.contains("FixedBody")) {
      std::cout << "SoftBody" << std::endl;
      SoftBodyObject.emplace_back(net);
      net->isSoftBody = true;
      net->isRigidBody = net->isFluid = false;
    } else if (type.contains("Fluid")) {
      std::cout << "Fluid" << std::endl;
      FluidObject.emplace_back(net);
      net->isFluid = true;
      net->isRigidBody = net->isSoftBody = false;

      // Parse initial condition for Fluid objects
      net->inputJSON.find("initial_condition", [&](auto STR_VEC) {
        if (STR_VEC.empty())
          return;
        std::string ic_type = toLowerCase(STR_VEC[0]);
        if (ic_type == "solitary_wave_theory") {
          auto vec = stod(std::vector<std::string>(STR_VEC.begin() + 1, STR_VEC.end()));
          if (vec.size() < 2)
            throw std::runtime_error("initial_condition solitary_wave_theory requires at least [H/h, depth]");
          double Hh = vec[0], h = vec[1];
          double bz = (vec.size() > 2) ? vec[2] : 0.0;
          double x0 = (vec.size() > 3) ? vec[3] : 0.0;
          auto sw = std::make_shared<TanakaSolitaryWave>();
          sw->solve(Hh, h, bz, x0);
          net->ic_phi = [sw](const Tddd& X, double t) { return sw->phi(X, t); };
          net->ic_eta = [sw](const Tddd& X, double t) { return sw->eta(X, t); };
          std::cout << "  [initial_condition] solitary_wave_theory: H/h=" << Hh << ", h=" << h << ", bz=" << bz << ", x0=" << x0 << std::endl;
        } else if (ic_type == "wave_theory") {
          auto vec = stod(std::vector<std::string>(STR_VEC.begin() + 1, STR_VEC.end()));
          if (vec.size() < 4)
            throw std::runtime_error("initial_condition wave_theory requires [A, T, h, bottom_z]");
          auto wt = std::make_shared<WaterWaveTheory>();
          wt->A = vec[0];
          wt->set_T_h(vec[1], vec[2]);
          wt->bottom_z = vec[3];
          if (vec.size() > 4)
            wt->theta = vec[4] / 180. * M_PI;
          net->ic_phi = [wt](const Tddd& X, double t) { return wt->phi(X, t); };
          net->ic_eta = [wt](const Tddd& X, double t) { return wt->eta(X, t); };
          std::cout << "  [initial_condition] wave_theory: A=" << wt->A << ", T=" << vec[1] << ", h=" << vec[2] << std::endl;
        } else if (ic_type == "wave_theory_l") {
          auto vec = stod(std::vector<std::string>(STR_VEC.begin() + 1, STR_VEC.end()));
          if (vec.size() < 4)
            throw std::runtime_error("initial_condition wave_theory_L requires [A, L, h, bottom_z]");
          auto wt = std::make_shared<WaterWaveTheory>();
          wt->A = vec[0];
          wt->set_L_h(vec[1], vec[2]);
          wt->bottom_z = vec[3];
          if (vec.size() > 4)
            wt->theta = vec[4] / 180. * M_PI;
          net->ic_phi = [wt](const Tddd& X, double t) { return wt->phi(X, t); };
          net->ic_eta = [wt](const Tddd& X, double t) { return wt->eta(X, t); };
          std::cout << "  [initial_condition] wave_theory_L: A=" << wt->A << ", L=" << vec[1] << ", h=" << vec[2] << std::endl;
        } else {
          throw std::runtime_error("Unknown initial_condition type: " + ic_type);
        }
      });

    } else if (type.contains("Absorber") || type.contains("absorb") || type.contains("damping")) {
      std::cout << "Absorber" << std::endl;
      AbsorberObject.emplace_back(net);
      net->isAbsorber = true;
      net->isRigidBody = net->isSoftBody = net->isFluid = false;
      net->inputJSON.find("wave_theory", [&](auto STR_VEC) {
        net->water_wave_theory = WaterWaveTheory();
        auto vec = stod(STR_VEC);
        net->water_wave_theory.A = vec[0];
        net->water_wave_theory.set_T_h(vec[1], vec[2]);
        net->water_wave_theory.bottom_z = vec[3];
        if (vec.size() > 4)
          net->water_wave_theory.theta = vec[4] / 180. * M_PI;
        net->absorb_velocity = [net](const Tddd& X, double t) { return net->water_wave_theory.gradPhi(X, t); };
        net->absorb_gradPhi_t = [net](const Tddd& X, double t) { return net->water_wave_theory.gradPhi_t(X, t); };
        net->absorb_phi = [net](const Tddd& X, const double t) { return net->water_wave_theory.phi(X, t); };
        net->absorb_eta = [net](const Tddd& X, const double t) { return net->water_wave_theory.eta(X, t); };
        net->absorb_gamma = [net](double sd) { return std::clamp(sd / (3. * net->water_wave_theory.L), 0., 1.); };
      });
      net->inputJSON.find("random_wave_theory", [&](auto STR_VEC) {
        double H13, T13, h, bottom_z;
        std::string tmp = STR_VEC[0];
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](char c) { return std::tolower(c); });
        if (tmp == "jonswap") {
          double height = std::stod(STR_VEC[1]);
          double period = std::stod(STR_VEC[2]);
          h = std::stod(STR_VEC[3]);
          bottom_z = std::stod(STR_VEC[4]);
          double gamma = (STR_VEC.size() >= 6) ? std::stod(STR_VEC[5]) : 3.3;
          // Check for parameter mode (7th element from GUI): default is Hm0_Tp for JONSWAP
          std::string param_mode_str = (STR_VEC.size() >= 7) ? STR_VEC[6] : "Hm0_Tp";
          WaveParamMode param_mode = WaveParamMode::HM0_TP;
          if (param_mode_str == "H13_T13")
            param_mode = WaveParamMode::H13_T13;
          std::cout << "Using JONSWAP (" << param_mode << ")" << std::endl;
          net->random_water_wave_theory =
              RandomWaterWaveTheory::create(SpectrumType::JONSWAP, param_mode,
                                            height, period, gamma, h, bottom_z);
        } else {

          if (STR_VEC.size() < 4)
            throw std::runtime_error("Default random_wave_theory requires 4 elements: [H13, T13, h, bottom_z]");

          auto vec = stod(STR_VEC);
          H13 = vec[0];
          T13 = vec[1];
          h = vec[2];
          bottom_z = vec[3];
          net->random_water_wave_theory =
              RandomWaterWaveTheory::Bretschneider(H13, T13, h, bottom_z);
        }

        std::cout << net->random_water_wave_theory;

        net->absorb_velocity = [net](const Tddd& X, double t) { return net->random_water_wave_theory.gradPhi(X, t); };
        net->absorb_gradPhi_t = [net](const Tddd& X, double t) { return net->random_water_wave_theory.gradPhi_t(X, t); };
        net->absorb_phi = [net](const Tddd& X, const double t) { return net->random_water_wave_theory.phi(X, t); };
        net->absorb_eta = [net](const Tddd& X, const double t) { return net->random_water_wave_theory.eta(X, t); };
        net->absorb_gamma = [net](double sd) { return std::clamp(sd / (3. * net->random_water_wave_theory.reference_wavelength), 0., 1.); };
      });
      net->inputJSON.find("solitary_wave_theory", [&](auto STR_VEC) {
        // Format: [H/h, depth, bottom_z, x_crest]
        // Example: [0.3, 1.0, 0.0, 0.0]
        auto vec = stod(STR_VEC);
        if (vec.size() < 2)
          throw std::runtime_error("solitary_wave_theory requires at least [H/h, depth]");
        double Hh = vec[0], h = vec[1];
        double bz = (vec.size() > 2) ? vec[2] : 0.0;
        double x0 = (vec.size() > 3) ? vec[3] : 0.0;
        auto sw = std::make_shared<TanakaSolitaryWave>();
        sw->solve(Hh, h, bz, x0);
        net->absorb_velocity = [sw](const Tddd& X, double t) { return sw->gradPhi(X, t); };
        net->absorb_gradPhi_t = [sw](const Tddd& X, double t) { return sw->gradPhi_t(X, t); };
        net->absorb_phi = [sw](const Tddd& X, const double t) { return sw->phi(X, t); };
        net->absorb_eta = [sw](const Tddd& X, const double t) { return sw->eta(X, t); };
        net->absorb_gamma = [sw](double sd) { return std::clamp(sd / (3. * sw->L), 0., 1.); };
      });
      net->inputJSON.find("wave_theory_L", [&](auto STR_VEC) {
        net->water_wave_theory = WaterWaveTheory();
        auto vec = stod(STR_VEC);
        net->water_wave_theory.A = vec[0];
        net->water_wave_theory.set_L_h(vec[1], vec[2]);
        net->water_wave_theory.bottom_z = vec[3];
        if (vec.size() > 4)
          net->water_wave_theory.theta = vec[4] / 180. * M_PI;
        net->absorb_velocity = [net](const Tddd& X, double t) { return net->water_wave_theory.gradPhi(X, t); };
        net->absorb_gradPhi_t = [net](const Tddd& X, double t) { return net->water_wave_theory.gradPhi_t(X, t); };
        net->absorb_phi = [net](const Tddd& X, const double t) { return net->water_wave_theory.phi(X, t); };
        net->absorb_eta = [net](const Tddd& X, const double t) { return net->water_wave_theory.eta(X, t); };
        net->absorb_gamma = [net](double sd) { return std::clamp(sd / (3. * net->water_wave_theory.L), 0., 1.); };
      });
    }

    //! isFixedはdefaultでfalse．指定された分だけ順に置き換わる
    //! ただし，指定が１つだけなら，それを全てのに適用する．
    //! 後方互換: "fixed" キーも受け入れる（"isFixed" を優先）
    const char* fixed_key = net->inputJSON.find("isFixed") ? "isFixed" : (net->inputJSON.find("fixed") ? "fixed" : nullptr);
    if (fixed_key) {
      if (net->inputJSON.find("fixed") && net->inputJSON.find("isFixed"))
        std::cerr << "Warning: both \"isFixed\" and \"fixed\" found in " << net->getName() << ".json. Using \"isFixed\"." << std::endl;
      auto isFixed = stob(net->inputJSON.at(fixed_key));
      if (isFixed.size() == 1)
        net->isFixed.fill(isFixed[0]);
      else
        for (auto i = 0; i < isFixed.size(); ++i)
          net->isFixed[i] = isFixed[i];
    }
    // もし，velocityがfixedなら，isFixedをtrueにする．
    if (net->inputJSON.find("velocity")) {
      auto velocity = net->inputJSON.at("velocity");
      if (velocity.size() == 1 && velocity[0] == "fixed") {
        net->isFixed.fill(true);
      } else if (velocity.size() > 1 && velocity[0] == "fixed") {
        for (auto i = 0; i < velocity.size() - 1 && i < net->isFixed.size(); ++i)
          net->isFixed[i] = true;
      }
    }
    // for (auto i = 0; i < 10; ++i)
    //    AreaWeightedSmoothingPreserveShape(net->getPoints(), 0.1);
    //$ set velocity
    std::cout << "set velocity" << std::endl;
    net->isFloatingBody = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating");
    // velocityにfileが指定されている場合は，そのファイルを読み込み，
    // interpolationBsplineであるnet->intpMotionRigidBodyをsetする．
    /*

    ##　物体の動きを値として与える場合

    以下は魚の体の動き{t,x,y,z,q0,q1,q2}をファイルに保存し，それを読み込む例

    ```python
    bodyA = {"name": "bodyA",
       "type": "RigidBody",
       "COM": [0., 0, 0.25],
       "mass": 10**10,
       "MOI": [10**10, 10**10, 10**10],
       "output": "json",
       #  "velocity": ["sin", 0, 0.1, 5, 0, 0, 0, 0, 0, 1],
       "velocity": ["file", "./study_fish/bodyA.dat"],
       "objfile": objfolder + "/bodyA20.obj"}
    ```

    この方法で，関数で与えられない複雑な動きも与えることができる．
    また，境界条件は，速度として与える必要があるが，ここでは位置を与えている．

    */

    net->inputJSON.find("velocity", [&](auto STR_VEC) {
      if (STR_VEC[0].contains("file")) {
        if (STR_VEC.size() == 1)
          throw std::runtime_error("Failed to open the input file.");
        std::ifstream file(STR_VEC[1]);
        std::cout << "displacement file : " << STR_VEC[1] << std::endl;
        if (!file.is_open())
          throw std::runtime_error("Failed to open the input file.");
        std::vector<double> T;
        std::vector<std::array<double, 6>> XYZ_Angles;
        std::string line;
        while (std::getline(file, line)) {
          if (line[0] == '#')
            continue;
          std::replace(line.begin(), line.end(), ',', ' ');
          line = std::regex_replace(line, std::regex("\\s+"), " ");
          std::istringstream iss(line);
          double t = 0.0, x = 0.0, y = 0.0, z = 0.0, q0 = 0.0, q1 = 0.0, q2 = 0.0;
          iss >> t >> x >> y >> z >> q0 >> q1 >> q2;
          T.push_back(t);
          XYZ_Angles.push_back({x, y, z, q0, q1, q2});
        }
        file.close();
        net->intpMotionRigidBody.set(3, T, XYZ_Angles);
      }
    });
    net->isFloatingBody = net->isFloatingBody || (net->inputJSON.find("acceleration") && net->inputJSON.at("acceleration")[0] == "floating");

    net->resetInitialX();
    net->setGeometricPropertiesForce();
  };
};
