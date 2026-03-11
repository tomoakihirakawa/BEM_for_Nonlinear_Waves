#pragma once

#include <filesystem>
#include <optional>

#include "BEM.hpp"
#include "basic.hpp"
#include "TanakaSolitaryWave.hpp"

std::string toLowerCase(const std::string &str) {
  std::string strLower = str;
  std::transform(strLower.begin(), strLower.end(), strLower.begin(), [](unsigned char c) { return std::tolower(c); });
  return strLower;
}

struct SimulationSettings {
  enum class DomainMode { Time, Frequency };
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
    struct NodeRelocationSettings {
      enum class Method { none, ALE, interpolation };
      enum class Surface { linear, pseudo_quadratic, true_quadratic };
      enum class MidpointMode { nearest, cubic_hermite };
      Method method = Method::none;
      int period = 0;
      Surface surface = Surface::pseudo_quadratic;
      MidpointMode midpoint_mode = MidpointMode::nearest;
      bool surface_explicitly_set = false;
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
      bool use_metal_m2l = false;           // Enable Metal GPU acceleration for M2L
      bool metal_m2l_threadgroup = false;   // true: use threadgroup parallelization
      bool metal_m2l_sort_terms = false;    // Sort terms for improved memory locality
      // Nearfield integration mode: "scalar" (default), "simd" (NEON float 4-target), "simd_double" (NEON double 2-target), "metal" (GPU)
      std::string nearfield_mode = "cell_scalar";
      // P2M quadrature: number of Dunavant points (1, 3, 6, 7). Default=1 for backward compatibility.
      int p2m_quadrature_points = 6;
      // MAC criterion parameter (smaller = stricter separation = faster but less accurate)
      double mac_theta = 0.25;
      // FMM tree structure settings (requires GMRES + FMM)
      int fmm_max_level = 7;               // Maximum tree depth
      int fmm_bucket_max_points = 50;      // Bucket split threshold (number of points)
    } solver;
  } bem;

  /* ---------------------------- remeshing settings -------------------------- */
  struct RemeshingSettings {
    // Used by remesh_for_main_loop
    double min_edge_length = 0.0;

    // Meshing options
    bool tetrahedralize = false;
    bool surface_flip = false;
    bool surface_split = false;
    bool surface_collapse = false;
    bool improve_tetrahedra = false;

    // Optional simulation-time windows for remeshing logic (currently consumed in main.cpp).
    double stop_remesh_time = 1E+10;
    double force_remesh_time = 0.0;
    int grid_refinement = 0;

    // Initial ALE-style mesh pre-relaxation (no physical time advance)
    struct InitialMeshPreRelax {
      // Optional; controlled by settings.json "initial_mesh_pre_relax".
      bool enabled = true;
      int loop = 50;
      double coef = 0.01;
    } initial_mesh_pre_relax;

    // Surface collision detection and resolution
    struct CollisionSettings {
      bool enabled = false;                  // Master enable/disable
      double proximity_factor = 1.5;         // mean_edge_length × this factor = detection threshold
      double normal_reversal_cos = -0.3;     // Adjacent face normal dot < this → folding detected
      int min_zone_faces = 3;                // Minimum faces in zone to attempt resolution
      bool detect_folding = true;            // Problem B: adjacent face folding detection
      bool detect_non_adjacent = true;       // Problem A: non-adjacent collision detection
      bool resolve_with_tetgen = true;       // Use TetGen local repair for non-adjacent collision zones
    } collision;

    struct SubsurfaceAltitudeRejectSettings {
      bool enabled = true;
      double min_face_altitude_rel = 0.05;
      double min_edge_angle_deg = 1.0;
      double max_edge_angle_deg = 120.0;
    } subsurface_altitude_reject;

    bool shell_visualization = false;
    bool front_advancing_debug = false;
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
    int interval = 0;           // 0=disabled, N=save every N steps
    std::string restart_from;   // checkpoint file path (empty=fresh start)
    int max_keep = 3;           // number of checkpoint files to keep (0=keep all)
  } checkpoint;

  /* ------------------------------- objects/data ----------------------------- */
  std::map<std::string, outputInfo> NetOutputInfo;
  std::vector<Network *> FluidObject, RigidBodyObject, SoftBodyObject, AbsorberObject;
  std::vector<JSON> MeasurementJSONs;

  SimulationSettings(std::filesystem::path input_directory, DomainMode mode_in = DomainMode::Time)
      : settings_file_path([&] {
          const auto canonical = input_directory / settings_filename;
          if (std::filesystem::exists(canonical))
            return canonical;
          const auto legacy = input_directory / legacy_setting_filename;
          if (std::filesystem::exists(legacy)) {
            std::cout << Yellow << "Warning: using legacy settings file name \"" << legacy_setting_filename
                      << "\". Please rename it to \"" << settings_filename << "\"." << colorReset << std::endl;
            return legacy;
          }
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                              "Missing settings file: expected \"" + settings_filename + "\" (or legacy \"" + legacy_setting_filename + "\") in: " + input_directory.string());
        }()),
        settingJSON(settings_file_path) {
    domain_mode = mode_in;
    common.input_directory = std::move(input_directory);
    const auto mode_label = (domain_mode == DomainMode::Time) ? std::string("time_domain") : std::string("frequency_domain");

    auto pick_key = [&](std::initializer_list<const char *> keys) -> std::optional<std::string> {
      for (const auto *k : keys) {
        if (settingJSON.find(k))
          return std::string(k);
      }
      return std::nullopt;
    };
    auto require_key = [&](std::initializer_list<const char *> keys, const std::string &label) -> std::string {
      auto hit = pick_key(keys);
      if (!hit) {
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                            "Missing key in settings.json (" + mode_label + "): " + label);
      }
      return *hit;
    };

    const auto output_key = require_key({"output_directory"}, "output_directory");
    const auto input_files_key = require_key({"input_files"}, "input_files");
    common.output_directory = settingJSON.at(output_key)[0];
    if (common.output_directory.is_relative())
      common.output_directory = common.input_directory / common.output_directory;
    // Allow overriding output directory without editing the input case.
    // Useful for sandboxed runs / short benchmarks.
    if (const char *out_env = std::getenv("BEM_OUTPUT_DIR")) {
      if (out_env[0] != '\0') {
        common.output_directory = std::filesystem::path(out_env);
        std::cout << Yellow << "Override output_directory via BEM_OUTPUT_DIR: " << common.output_directory << colorReset << std::endl;
        std::filesystem::create_directories(common.output_directory);
      }
    }
    (void)input_files_key;

    // Helper: parse node_relocation from a string array
    // ["method", period, "surface", "midpoint_mode"]
    auto parse_node_relocation = [&](const std::vector<std::string> &arr) {
      using M = TimeDomainSettings::NodeRelocationSettings::Method;
      using S = TimeDomainSettings::NodeRelocationSettings::Surface;
      using MM = TimeDomainSettings::NodeRelocationSettings::MidpointMode;
      auto &nr = time.node_relocation;
      auto method_str = toLowerCase(arr[0]);
      if (method_str == "none")
        nr.method = M::none;
      else if (method_str == "ale")
        nr.method = M::ALE;
      else if (method_str == "interpolation")
        nr.method = M::interpolation;
      else
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                            "node_relocation method must be none, ALE, or interpolation (got: " + arr[0] + ")");
      if (arr.size() >= 2)
        nr.period = std::stoi(arr[1]);
      if (arr.size() >= 3) {
        auto surf_str = toLowerCase(arr[2]);
        nr.surface_explicitly_set = true;
        if (surf_str == "linear")
          nr.surface = S::linear;
        else if (surf_str.contains("true"))
          nr.surface = S::true_quadratic;
        else if (surf_str.contains("pseudo") || surf_str.contains("quad"))
          nr.surface = S::pseudo_quadratic;
        else
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                              "node_relocation surface must be linear, pseudo_quadratic, or true_quadratic (got: " + arr[2] + ")");
      }
      if (arr.size() >= 4) {
        auto midpoint_mode_str = toLowerCase(arr[3]);
        if (midpoint_mode_str == "nearest")
          nr.midpoint_mode = MM::nearest;
        else if (midpoint_mode_str == "cubic_hermite" || midpoint_mode_str == "hermite")
          nr.midpoint_mode = MM::cubic_hermite;
        else
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                              "node_relocation midpoint_mode must be nearest or cubic_hermite (got: " + arr[3] + ")");
      }
    };

    // Helper: parse legacy ALE/ALEPERIOD keys into node_relocation
    auto parse_legacy_ale = [&](const std::string &ale_value, int period) {
      using M = TimeDomainSettings::NodeRelocationSettings::Method;
      using S = TimeDomainSettings::NodeRelocationSettings::Surface;
      auto &nr = time.node_relocation;
      nr.period = period;
      auto val = toLowerCase(ale_value);
      if (val == "none") {
        nr.method = M::none;
      } else if (val == "remap") {
        nr.method = M::interpolation;
      } else {
        nr.method = M::ALE;
        if (val == "linear") {
          nr.surface = S::linear;
          nr.surface_explicitly_set = true;
        }
        // "pseudo_quad", "true_quadratic" etc. → surface auto-resolved from element type
      }
      std::cout << Yellow << "Note: 'ALE'/'ALEPERIOD' keys are deprecated. Use 'node_relocation' instead. "
                << "(e.g., \"node_relocation\": [\"ALE\", " << period << "])" << colorReset << std::endl;
    };

    if (domain_mode == DomainMode::Time) {
      const auto end_step_key = require_key({"time_end_time_step", "end_time_step"}, "end_time_step");
      const auto end_time_key = require_key({"time_end_time", "end_time"}, "end_time");
      const auto max_dt_key = require_key({"time_max_dt", "max_dt"}, "max_dt");

      time.end_time_step = stoi(settingJSON.at(end_step_key))[0];
      time.end_time = stod(settingJSON.at(end_time_key))[0];
      time.max_dt = stod(settingJSON.at(max_dt_key))[0];

      if (const auto nr_key = pick_key({"node_relocation"})) {
        parse_node_relocation(settingJSON.at(*nr_key));
      } else {
        // Legacy fallback
        const auto ale_key = require_key({"time_ALE", "ALE"}, "ALE");
        const auto ale_period_key = require_key({"time_ALEPERIOD", "ALEPERIOD"}, "ALEPERIOD");
        parse_legacy_ale(settingJSON.at(ale_key)[0], stoi(settingJSON.at(ale_period_key))[0]);
      }
    } else {
      if (const auto end_step_key = pick_key({"time_end_time_step", "end_time_step"}))
        time.end_time_step = stoi(settingJSON.at(*end_step_key))[0];
      if (const auto end_time_key = pick_key({"time_end_time", "end_time"}))
        time.end_time = stod(settingJSON.at(*end_time_key))[0];
      if (const auto max_dt_key = pick_key({"time_max_dt", "max_dt"}))
        time.max_dt = stod(settingJSON.at(*max_dt_key))[0];

      if (const auto nr_key = pick_key({"node_relocation"})) {
        parse_node_relocation(settingJSON.at(*nr_key));
      } else {
        int period = 0;
        if (const auto ale_period_key = pick_key({"time_ALEPERIOD", "ALEPERIOD"}))
          period = stoi(settingJSON.at(*ale_period_key))[0];
        if (const auto ale_key = pick_key({"time_ALE", "ALE"}))
          parse_legacy_ale(settingJSON.at(*ale_key)[0], period);
      }
    }

    bool coupling_explicitly_set = false;
    {

      std::cout << "input_directory : " << common.input_directory << std::endl;

      /* -------------------------------------------------------------------------- */
      /*                           read settings.json                               */
      /* -------------------------------------------------------------------------- */

      for (auto &[key, value] : settingJSON())
        std::cout << key << ": " << value << std::endl;

      // Required key validation already handled above (mode-aware).

      /* ----------------------------- meshing options ---------------------------- */

      settingJSON.find("meshing_options", [&](const auto &STR_VEC) {
        for (const auto &STR_VEC : STR_VEC) {
          if (STR_VEC == "tetrahedralize")
            remeshing.tetrahedralize = true;
          if (STR_VEC == "surface_flip")
            remeshing.surface_flip = true;
          if (STR_VEC == "surface_split")
            remeshing.surface_split = true;
          if (STR_VEC == "surface_collapse")
            remeshing.surface_collapse = true;
          if (STR_VEC == "improve_tetrahedra")
            remeshing.improve_tetrahedra = true;
        }
      });

      /* -------------------------------------------------------------------------- */

      if (settingJSON.find("element")) {
        auto elem_str = toLowerCase(settingJSON.at("element")[0]);
        if (elem_str == "linear") {
          bem.element.linear = true;
        } else if (elem_str == "true_quadratic" || elem_str == "true_quad" || elem_str == "quadratic") {
          bem.element.linear = false;
          bem.element.true_quadratic = true;
        } else if (elem_str.contains("quad") && elem_str.contains("pseudo")) {
          bem.element.linear = false;
          bem.element.pseudo_quadratic = true;
        } else {
          bem.element.linear = true;
        }
      } else {
        bem.element.linear = true;
      }

      if (bem.element.linear)
        std::cout << "LINEAR_ELEMENT" << std::endl;
      if (bem.element.pseudo_quadratic)
        std::cout << "PSEUDO_QUADRATIC_ELEMENT" << std::endl;
      if (bem.element.true_quadratic)
        std::cout << "TRUE_QUADRATIC_ELEMENT" << std::endl;

      {
        using M = TimeDomainSettings::NodeRelocationSettings::Method;
        using S = TimeDomainSettings::NodeRelocationSettings::Surface;
        using MM = TimeDomainSettings::NodeRelocationSettings::MidpointMode;
        auto &nr = time.node_relocation;
        if (nr.method == M::none)
          std::cout << "NODE_RELOCATION: none (pure Lagrangian)" << std::endl;
        else {
          std::cout << "NODE_RELOCATION: " << (nr.method == M::ALE ? "ALE" : "interpolation")
                    << ", period=" << nr.period;
          if (nr.method == M::interpolation)
            std::cout << ", midpoint_mode=" << (nr.midpoint_mode == MM::cubic_hermite ? "cubic_hermite" : "nearest");
          std::cout << std::endl;
        }
      }

      std::cout << "" << std::endl;

      if (settingJSON.find("WATER_DENSITY", [&](auto STR_VEC) {
            common.water_density = stod(STR_VEC[0]);
            _WATER_DENSITY_ = common.water_density;
          }))
        std::cout << "WATER_DENSITY: " << common.water_density << std::endl;

      if (settingJSON.find("GRAVITY", [&](auto STR_VEC) {
            common.gravity = stod(STR_VEC[0]);
            _GRAVITY_ = common.gravity;
          }))
        std::cout << "GRAVITY: " << common.gravity << std::endl;

      auto set_preconditioner = [&](std::string v) {
        v = toLowerCase(StringTrim(v, {" "}));
        if (v.empty() || v == "none" || v == "off" || v == "no" || v == "false" || v == "0") {
          bem.solver.preconditioner_type = "NONE";
          return;
        }
        if (v == "ilu") {
          bem.solver.preconditioner_type = "ILU";
          return;
        }
        if (v == "milu" || v == "milu0" || v == "milu(0)") {
          bem.solver.preconditioner_type = "MILU";
          return;
        }
        if (v == "ilut") {
          bem.solver.preconditioner_type = "ILUT";
          return;
        }
        if (v == "schwarz" || v == "additive_schwarz" || v == "as" || v == "block-jacobi" || v == "block_jacobi" || v == "bj") {
          bem.solver.preconditioner_type = "SCHWARZ";
          return;
        }
        if (v == "diagonal" || v == "jacobi") {
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "preconditioner \"" + v + "\" is not supported (removed). Use \"ILU\", \"MILU\", \"ILUT\", \"SCHWARZ\", or \"NONE\".");
        }
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Unknown preconditioner \"" + v + "\". Supported: \"ILU\", \"MILU\", \"ILUT\", \"SCHWARZ\", \"NONE\".");
      };

      auto set_ilu_neighborhood = [&](std::string type, const std::optional<std::string> &num_str = std::nullopt) {
        type = toLowerCase(StringTrim(type, {" "}));
        if (type == "buckets" || type == "bucket") {
          bem.solver.ilu_neighborhood_type = "BUCKETS";
          return;
        }
        if (type == "k-ring" || type == "kring" || type == "k_ring" || type == "k-ring ") {
          if (!num_str.has_value())
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ILU neighborhood \"k-ring\" requires an integer num (0..20).");
          const int num = std::stoi(StringTrim(*num_str, {" "}));
          if (num < 0 || num > 20)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ILU k-ring num must be in [0,20].");
          bem.solver.ilu_neighborhood_type = "K-RING";
          bem.solver.ilu_kring_num = num;
          return;
        }
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Unknown ILU neighborhood \"" + type + "\". Supported: \"k-ring\" (with num 0..20), \"buckets\".");
      };

      auto set_coupling = [&](std::string v) {
        v = toLowerCase(StringTrim(v, {" "}));
        if (v.empty() || v == "none" || v == "off") {
          bem.solver.coupling_type = "NONE";
          return;
        }
        if (v == "anderson" || v == "aa") {
          bem.solver.coupling_type = "ANDERSON";
          return;
        }
        if (v == "aitken") {
          bem.solver.coupling_type = "AITKEN";
          return;
        }
        if (v == "broyden") {
          bem.solver.coupling_type = "BROYDEN";
          return;
        }
        if (v == "newton") {
          bem.solver.coupling_type = "NEWTON";
          return;
        }
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Unknown coupling \"" + v + "\". Supported: \"NONE\", \"ANDERSON\", \"AITKEN\", \"BROYDEN\".");
      };

      // Solver settings
      settingJSON.find("solver", [&](const auto &v) {
        if (v.size() >= 1)
          bem.solver.solver_type = StringTrim(v[0], {" "});
        if (v.size() >= 2)
          bem.solver.solver_tol = std::stod(v[1]);
        if (v.size() >= 3)
          bem.solver.solver_max_iter = std::stoi(v[2]);
        if (v.size() >= 4)
          bem.solver.solver_restart = std::stoi(v[3]);
        if (v.size() >= 5)
          set_preconditioner(v[4]);
        if (v.size() >= 6) {
          if (bem.solver.preconditioner_type == "ILU" || bem.solver.preconditioner_type == "MILU" || bem.solver.preconditioner_type == "ILUT" || bem.solver.preconditioner_type == "SCHWARZ") {
            set_ilu_neighborhood(v[5], (v.size() >= 7) ? std::optional<std::string>(v[6]) : std::nullopt);
            if (bem.solver.preconditioner_type == "ILUT") {
              // Optional ILUT params:
              // solver (examples):
              // - [GMRES, tol, max_iter, restart, ILUT, buckets, <drop_tol>, <max_entries_per_row>, <pivot_min>]
              // - [GMRES, tol, max_iter, restart, ILUT, k-ring, <num>, <drop_tol>, <max_entries_per_row>, <pivot_min>]
              std::size_t idx = 6;
              if (bem.solver.ilu_neighborhood_type == "K-RING")
                idx = 7;
              if (v.size() > idx)
                bem.solver.ilut_drop_tol = std::stod(v[idx]);
              if (v.size() > idx + 1)
                bem.solver.ilut_max_entries_per_row = std::stoi(v[idx + 1]);
              if (v.size() > idx + 2)
                bem.solver.ilut_pivot_min = std::stod(v[idx + 2]);
            } else if (bem.solver.preconditioner_type == "MILU") {
              // Optional MILU params:
              // - [GMRES, tol, max_iter, restart, MILU, buckets, <omega>]
              // - [GMRES, tol, max_iter, restart, MILU, k-ring, <num>, <omega>]
              std::size_t idx = 6;
              if (bem.solver.ilu_neighborhood_type == "K-RING")
                idx = 7;
              if (v.size() > idx)
                bem.solver.milu_omega = std::stod(v[idx]);
            } else if (bem.solver.preconditioner_type == "SCHWARZ") {
              if (bem.solver.ilu_neighborhood_type != "K-RING") {
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "SCHWARZ preconditioner currently requires ILU neighborhood type \"k-ring\".");
              }
              // Optional SCHWARZ params:
              // - [GMRES, tol, max_iter, restart, SCHWARZ, k-ring, <num>, <core_k>, <overlap_k>, <max_core_size>, <max_block_size>, <pivot_min>, <diag_shift>]
              std::size_t idx = 7; // after k-ring num
              if (v.size() > idx)
                bem.solver.schwarz_core_k = std::stoi(v[idx]);
              if (v.size() > idx + 1)
                bem.solver.schwarz_overlap_k = std::stoi(v[idx + 1]);
              if (v.size() > idx + 2)
                bem.solver.schwarz_max_core_size = std::stoi(v[idx + 2]);
              if (v.size() > idx + 3)
                bem.solver.schwarz_max_block_size = std::stoi(v[idx + 3]);
              if (v.size() > idx + 4)
                bem.solver.schwarz_pivot_min = std::stod(v[idx + 4]);
              if (v.size() > idx + 5)
                bem.solver.schwarz_diag_shift = std::stod(v[idx + 5]);
            }
          } else {
            std::cout << Yellow << "Warning: extra solver parameters after preconditioner are ignored because preconditioner=\"" << bem.solver.preconditioner_type << "\"." << colorReset << std::endl;
          }
        }
      });
      settingJSON.find("solver_tol", [&](const auto &v) { bem.solver.solver_tol = std::stod(v[0]); });
      settingJSON.find("solver_max_iter", [&](const auto &v) { bem.solver.solver_max_iter = std::stoi(v[0]); });
      settingJSON.find("solver_restart", [&](const auto &v) { bem.solver.solver_restart = std::stoi(v[0]); });
      settingJSON.find("preconditioner", [&](const auto &v) { set_preconditioner(v[0]); });
      settingJSON.find("milu_omega", [&](const auto &v) { bem.solver.milu_omega = std::stod(v[0]); });
      settingJSON.find("ilut_drop_tol", [&](const auto &v) { bem.solver.ilut_drop_tol = std::stod(v[0]); });
      settingJSON.find("ilut_max_entries_per_row", [&](const auto &v) { bem.solver.ilut_max_entries_per_row = std::stoi(v[0]); });
      settingJSON.find("ilut_pivot_min", [&](const auto &v) { bem.solver.ilut_pivot_min = std::stod(v[0]); });
      settingJSON.find("schwarz_core_k", [&](const auto &v) { bem.solver.schwarz_core_k = std::stoi(v[0]); });
      settingJSON.find("schwarz_overlap_k", [&](const auto &v) { bem.solver.schwarz_overlap_k = std::stoi(v[0]); });
      settingJSON.find("schwarz_max_core_size", [&](const auto &v) { bem.solver.schwarz_max_core_size = std::stoi(v[0]); });
      settingJSON.find("schwarz_max_block_size", [&](const auto &v) { bem.solver.schwarz_max_block_size = std::stoi(v[0]); });
      settingJSON.find("schwarz_pivot_min", [&](const auto &v) { bem.solver.schwarz_pivot_min = std::stod(v[0]); });
      settingJSON.find("schwarz_diag_shift", [&](const auto &v) { bem.solver.schwarz_diag_shift = std::stod(v[0]); });
      settingJSON.find("coupling", [&](const auto &v) {
        coupling_explicitly_set = true;
        if (v.size() >= 1)
          set_coupling(v[0]);
        if (v.size() >= 2)
          bem.solver.coupling_tol = std::stod(v[1]);
        if (v.size() > 2) {
          bem.solver.coupling_params.clear();
          for (size_t i = 2; i < v.size(); ++i)
            bem.solver.coupling_params.push_back(std::stod(v[i]));
        }
      });

      // Metal M2L GPU acceleration settings (requires GMRES + FMM)
      settingJSON.find("use_metal_m2l", [&](const auto &v) { bem.solver.use_metal_m2l = stob(v)[0]; });
      settingJSON.find("metal_m2l_threadgroup", [&](const auto &v) { bem.solver.metal_m2l_threadgroup = stob(v)[0]; });
      settingJSON.find("metal_m2l_sort_terms", [&](const auto &v) { bem.solver.metal_m2l_sort_terms = stob(v)[0]; });

      // P2M quadrature points (Dunavant rule)
      settingJSON.find("p2m_quadrature_points", [&](const auto &v) { bem.solver.p2m_quadrature_points = std::stoi(v[0]); });

      // MAC criterion parameter
      settingJSON.find("mac_theta", [&](const auto &v) { bem.solver.mac_theta = std::stod(v[0]); });

      // Nearfield integration mode
      settingJSON.find("nearfield_mode", [&](const auto &v) { bem.solver.nearfield_mode = v[0]; });

      // FMM tree structure settings (requires GMRES + FMM)
      settingJSON.find("fmm_max_level", [&](const auto &v) { bem.solver.fmm_max_level = stoi(v)[0]; });
      settingJSON.find("fmm_bucket_max_points", [&](const auto &v) { bem.solver.fmm_bucket_max_points = stoi(v)[0]; });

      // Remeshing schedule / controls
      settingJSON.find("stop_remesh_time", [&](auto STR_VEC) { remeshing.stop_remesh_time = stod(STR_VEC[0]); });
      settingJSON.find("force_remesh_time", [&](auto STR_VEC) { remeshing.force_remesh_time = stod(STR_VEC[0]); });
      settingJSON.find("grid_refinement", [&](auto STR_VEC) { remeshing.grid_refinement = stoi(STR_VEC[0]); });
      settingJSON.find("min_edge_length", [&](const auto &v) { remeshing.min_edge_length = std::stod(v[0]); });

      // Initial mesh pre-relaxation (ALE-style), before the first remesh/collapse.
      // settings.json spec: "initial_mesh_pre_relax": ["true", "<loop:int>", "<coef:double>"]
      settingJSON.find("initial_mesh_pre_relax", [&](const auto &v) {
        if (v.empty())
          return;
        remeshing.initial_mesh_pre_relax.enabled = stob(v)[0];
        if (!remeshing.initial_mesh_pre_relax.enabled)
          return;
        if (v.size() < 3) {
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "initial_mesh_pre_relax must be: [true, loop(int), coef(double)]");
        }
        remeshing.initial_mesh_pre_relax.loop = std::stoi(v[1]);
        remeshing.initial_mesh_pre_relax.coef = std::stod(v[2]);
      });

      // Surface collision detection and resolution
      settingJSON.find("collision_enabled", [&](const auto &v) { remeshing.collision.enabled = stob(v)[0]; });
      settingJSON.find("collision_proximity_factor", [&](const auto &v) { remeshing.collision.proximity_factor = std::stod(v[0]); });
      settingJSON.find("collision_normal_reversal_cos", [&](const auto &v) { remeshing.collision.normal_reversal_cos = std::stod(v[0]); });
      settingJSON.find("collision_min_zone_faces", [&](const auto &v) { remeshing.collision.min_zone_faces = std::stoi(v[0]); });
      settingJSON.find("collision_detect_folding", [&](const auto &v) { remeshing.collision.detect_folding = stob(v)[0]; });
      settingJSON.find("collision_detect_non_adjacent", [&](const auto &v) { remeshing.collision.detect_non_adjacent = stob(v)[0]; });
      settingJSON.find("collision_resolve_with_tetgen", [&](const auto &v) { remeshing.collision.resolve_with_tetgen = stob(v)[0]; });
      settingJSON.find("subsurface_altitude_reject", [&](const auto &v) { remeshing.subsurface_altitude_reject.enabled = stob(v)[0]; });
      settingJSON.find("subsurface_min_face_altitude_rel", [&](const auto &v) { remeshing.subsurface_altitude_reject.min_face_altitude_rel = std::stod(v[0]); });
      settingJSON.find("subsurface_min_edge_angle_deg", [&](const auto &v) { remeshing.subsurface_altitude_reject.min_edge_angle_deg = std::stod(v[0]); });
      settingJSON.find("subsurface_max_edge_angle_deg", [&](const auto &v) { remeshing.subsurface_altitude_reject.max_edge_angle_deg = std::stod(v[0]); });
      settingJSON.find("shell_visualization", [&](const auto &v) { remeshing.shell_visualization = stob(v)[0]; });
      settingJSON.find("front_advancing_debug", [&](const auto &v) { remeshing.front_advancing_debug = stob(v)[0]; });

      // VPM (optional)
      settingJSON.find("use_VPM", [&](const auto &v) { vpm.enabled = stob(v)[0]; });
      settingJSON.find("VPM_wall_min_absorb_receivers", [&](const auto &v) { vpm.wall_min_absorb_receivers = static_cast<std::size_t>(std::stoll(v[0])); });
      settingJSON.find("VPM_wall_min_absorb_total_weight", [&](const auto &v) { vpm.wall_min_absorb_total_weight = std::stod(v[0]); });
      settingJSON.find("VPM_sigma_factor", [&](const auto &v) { vpm.sigma_factor = std::stod(v[0]); });
      settingJSON.find("VPM_stretching_scheme", [&](const auto &v) {
        if (!v.empty())
          vpm.stretching_scheme = v[0];
      });
      settingJSON.find("VPM_PSE_correction", [&](const auto &v) {
        if (!v.empty())
          vpm.PSE_correction = v[0];
      });

      // Checkpoint settings
      settingJSON.find("checkpoint_interval", [&](const auto &v) { checkpoint.interval = std::stoi(v[0]); });
      settingJSON.find("restart_from_checkpoint", [&](const auto &v) { if (!v.empty()) checkpoint.restart_from = v[0]; });
      settingJSON.find("checkpoint_max_keep", [&](const auto &v) { checkpoint.max_keep = std::stoi(v[0]); });

      // Frequency-domain optional settings (parsed even in time mode; used by main_freq_domain).
      auto parse_double_list = [&](const std::vector<std::string> &v, std::vector<double> &out) {
        out.clear();
        out.reserve(v.size());
        for (const auto &s : v)
          out.push_back(std::stod(s));
      };
      auto parse_int_list = [&](const std::vector<std::string> &v, std::vector<int> &out) {
        out.clear();
        out.reserve(v.size());
        for (const auto &s : v)
          out.push_back(std::stoi(s));
      };
      auto set_omegas = [&](const std::vector<std::string> &v) { parse_double_list(v, frequency.omegas); };
      auto set_dofs = [&](const std::vector<std::string> &v) { parse_int_list(v, frequency.dofs); };
      if (!settingJSON.find("freq_omegas", set_omegas)) {
        if (!settingJSON.find("frequency_omegas", set_omegas)) {
          settingJSON.find("omegas", set_omegas);
        }
      }
      if (!settingJSON.find("freq_dofs", set_dofs)) {
        if (!settingJSON.find("frequency_dofs", set_dofs)) {
          settingJSON.find("dofs", set_dofs);
        }
      }

      /* -------------------------------------------------------------------------- */
      /*                read JSON files, input_files in settings.json               */
      /* -------------------------------------------------------------------------- */

      for (auto input_file_name : settingJSON["input_files"]) {
        std::cout << Green << common.input_directory / input_file_name << colorReset << std::endl;
        JSON injson(common.input_directory / input_file_name);

        /* -------------------- display contents of the JSON file ------------------- */
        for (auto &[key, value] : injson())
          std::cout << Green << key << colorReset << ": " << value << std::endl;
        auto object_name = injson.at("name")[0];

        /* --------------------------- check required keys -------------------------- */
        if (!injson.find("name"))
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"name\" key");
        if (!injson.find("type"))
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"type\" key");
        auto type = injson.at("type")[0];
        if ((type.contains("Fluid") || type.contains("Body")) && !injson.find("objfile"))
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"objfile\" key");
        else if (type.contains("Measurement") || type.contains("gauge")) {
          MeasurementJSONs.emplace_back(injson);
          std::cout << "type = " << type << std::endl;
          std::cout << "skipped" << std::endl;
          continue;
        }

        /* ------------------------ create Network object --------------------------- */

        if (!injson.find("ignore") || !stob(injson["ignore"])[0]) {

          std::string objfile_str = injson.at("objfile")[0];
          std::filesystem::path objpath(objfile_str);
          if (objpath.is_relative())
            objpath = common.input_directory / objpath;
          auto net = new Network(objpath.string(), object_name);
          net->inputJSON = injson;
          net->applyTransformations(injson);
          setOutputInfo(injson.at("name")[0], common.output_directory);            //* for boundary surface */
          setOutputInfo(injson.at("name")[0] + "_tetra", common.output_directory); //* for volume mesh */
          if (type.contains("Fluid")) {
            setOutputInfo(injson.at("name")[0] + "_shell_faces", common.output_directory);
            setOutputInfo(injson.at("name")[0] + "_shell_tets", common.output_directory);
            setOutputInfo(injson.at("name")[0] + "_shell_inner_tets", common.output_directory);
            setOutputInfo(injson.at("name")[0] + "_shell_outer_tets", common.output_directory);
            setOutputInfo(injson.at("name")[0] + "_aft_candidates", common.output_directory);
          }

          /* ----------------------------- 境界面データの回転と平行移動 ----------------------------- */
          if (injson.find("rotation")) {
            const auto angle_theta = stod(injson.at("rotation"));
            std::array<double, 3> axis = {angle_theta[0], angle_theta[1], angle_theta[2]};
            const auto angle = angle_theta[3];
            net->rotate(angle, axis);
          }
          if (injson.find("translation")) {
            const auto translation = stod(injson.at("translation"));
            net->translate(std::array<double, 3>{translation[0], translation[1], translation[2]});
          }

          setTypes(net);
          std::filesystem::copy_file(common.input_directory / input_file_name, common.output_directory / input_file_name, std::filesystem::copy_options::overwrite_existing);
          mk_vtu(common.output_directory / (object_name + "_init.vtu"), {net->getFaces()});
          /* -------------------------------------------------------------------------- */
          for (auto &[key, value] : injson()) {
            if (key.contains("mooring")) {
              //$ ------------------------ MOORING ---------------------- */
              //$ extract key the contains "mooring" as key name. extract them as vector
              std::cout << "/* ------------------------ MOORING ---------------------- */" << std::endl;
              std::cout << "initialize mooring" << std::endl;

              //! check if the value contains 12 elements
              if (value.size() != 13) //! show contents of the value nicely
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mooring line must have 12 elements");

              const auto name = value[0];
              const int n_points = std::stoi(value[8]); //$ number of points
              const std::array<double, 3> X_begin = {std::stod(value[1]), std::stod(value[2]), std::stod(value[3])};
              const std::array<double, 3> X_end = {std::stod(value[4]), std::stod(value[5]), std::stod(value[6])};
              const double total_length = std::stod(value[7]); //! [m]
              const double w = std::stod(value[9]);            //! [kg/m]
              const double stiffness = std::stod(value[10]);   //! [N/m]
              const double damp = std::stod(value[11]);        //! [N/(m/s^2)]
              const double diam = std::stod(value[12]);        //! [m]

              std::cout << std::right << std::setw(15) << "name : " << name << std::endl;
              std::cout << std::right << std::setw(15) << "X_begin : " << X_begin << std::endl;
              std::cout << std::right << std::setw(15) << "X_end : " << X_end << std::endl;
              std::cout << std::right << std::setw(15) << "total_length : " << total_length << std::endl;
              std::cout << std::right << std::setw(15) << "n_points : " << n_points << std::endl;
              std::cout << std::right << std::setw(15) << "mass per unit length : " << w << std::endl;
              std::cout << std::right << std::setw(15) << "stiffness : " << stiffness << std::endl;
              std::cout << std::right << std::setw(15) << "damp : " << damp << std::endl;
              std::cout << std::right << std::setw(15) << "diam : " << diam << std::endl;
              std::cout << std::right << std::setw(15) << "total mass" << w * total_length << std::endl;

              std::cout << std::right << std::setw(15) << "initialize MooringLine" << std::endl;
              auto mooring_net = new MooringLine(X_begin, X_end, total_length, n_points);
              mooring_net->setName(name);

              std::cout << std::right << std::setw(15) << "MooringLine->getPoints().size() = " << mooring_net->getPoints().size() << std::endl;
              std::cout << std::right << std::setw(15) << "MooringLine->getName() = " << mooring_net->getName() << std::endl;

              std::cout << std::right << std::setw(15) << "MooringLine->setDensityStiffnessDampingDiameter" << std::endl;
              mooring_net->setDensityStiffnessDampingDiameter(w, stiffness, damp, diam);

              net->mooringLines.emplace_back(mooring_net);
              setOutputInfo(name, common.output_directory);

              std::cout << std::right << std::setw(15) << "MooringLine->setEquilibriumState" << std::endl;
              auto boundary_condition = [&](networkPoint *p) {
                if (p == mooring_net->lastPoint || p == mooring_net->firstPoint) {
                  p->acceleration.fill(0.);
                  p->velocity.fill(0.);
                }
              };

              mooring_net->setEquilibriumState(boundary_condition);

              std::cout << "/* ------------------------------------------------------- */" << std::endl;
            }
          }
          //$ ------------------------------------------------------- */
        } else
          Print("skipped");
      }
    }

    // If coupling is not explicitly specified and there are no movable bodies, disable coupling.
    if (!coupling_explicitly_set) {
      bool has_movable_body = false;
      for (auto *net : Join(RigidBodyObject, SoftBodyObject)) {
        for (const auto &fixed : net->isFixed) {
          if (!fixed) {
            has_movable_body = true;
            break;
          }
        }
        if (has_movable_body)
          break;
      }
      if (!has_movable_body) {
        bem.solver.coupling_type = "NONE";
        bem.solver.coupling_params.clear();
      }
    }
  }

  /* -------------------------------------------------------------------------- */

  void setOutputInfo(auto output_name, const std::filesystem::path &output_directory) {
    std::cout << "setOutputInfo" << std::endl;
    this->NetOutputInfo[output_name].pvd_file_name = output_name;
    this->NetOutputInfo[output_name].vtu_file_name = output_name + "_";
    this->NetOutputInfo[output_name].PVD = new PVDWriter(output_directory / (output_name + ".pvd"));
  };

  /* -------------------------------------------------------------------------- */

  void setTypes(Network *net) {
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
          net->ic_phi = [sw](const Tddd &X, double t) { return sw->phi(X, t); };
          net->ic_eta = [sw](const Tddd &X, double t) { return sw->eta(X, t); };
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
          net->ic_phi = [wt](const Tddd &X, double t) { return wt->phi(X, t); };
          net->ic_eta = [wt](const Tddd &X, double t) { return wt->eta(X, t); };
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
          net->ic_phi = [wt](const Tddd &X, double t) { return wt->phi(X, t); };
          net->ic_eta = [wt](const Tddd &X, double t) { return wt->eta(X, t); };
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
        net->absorb_velocity = [net](const Tddd &X, double t) { return net->water_wave_theory.gradPhi(X, t); };
        net->absorb_gradPhi_t = [net](const Tddd &X, double t) { return net->water_wave_theory.gradPhi_t(X, t); };
        net->absorb_phi = [net](const Tddd &X, const double t) { return net->water_wave_theory.phi(X, t); };
        net->absorb_eta = [net](const Tddd &X, const double t) { return net->water_wave_theory.eta(X, t); };
        net->absorb_gamma = [net](double sd) { return std::clamp(sd / (3. * net->water_wave_theory.L), 0., 1.); };
      });
      net->inputJSON.find("random_wave_theory", [&](auto STR_VEC) {
        double H13, T13, h, bottom_z;
        std::string tmp = STR_VEC[0];
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](char c) { return std::tolower(c); });
        if (tmp == "jonswap") {
          std::cout << "Using JONSWAP" << std::endl;
          H13 = std::stod(STR_VEC[1]);
          T13 = std::stod(STR_VEC[2]);
          h = std::stod(STR_VEC[3]);
          bottom_z = std::stod(STR_VEC[4]);
          net->random_water_wave_theory = RandomWaterWaveTheory(H13, T13, h, bottom_z);
          if (STR_VEC.size() == 6)
            net->random_water_wave_theory.gamma = std::stod(STR_VEC[5]);
          else
            net->random_water_wave_theory.gamma = 3.3;

          net->random_water_wave_theory.setSpectrumType(SpectrumType::JONSWAP);
        } else {

          if (STR_VEC.size() < 4)
            throw std::runtime_error("Default random_wave_theory requires 4 elements: [H13, T13, h, bottom_z]");

          auto vec = stod(STR_VEC);
          H13 = vec[0];
          T13 = vec[1];
          h = vec[2];
          bottom_z = vec[3];
          net->random_water_wave_theory = RandomWaterWaveTheory(H13, T13, h, bottom_z);
        }

        std::cout << net->random_water_wave_theory;

        net->absorb_velocity = [net](const Tddd &X, double t) { return net->random_water_wave_theory.gradPhi(X, t); };
        net->absorb_gradPhi_t = [net](const Tddd &X, double t) { return net->random_water_wave_theory.gradPhi_t(X, t); };
        net->absorb_phi = [net](const Tddd &X, const double t) { return net->random_water_wave_theory.phi(X, t); };
        net->absorb_eta = [net](const Tddd &X, const double t) { return net->random_water_wave_theory.eta(X, t); };
        net->absorb_gamma = [net](double sd) { return std::clamp(sd / (3. * net->random_water_wave_theory.L13), 0., 1.); };
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
        net->absorb_velocity = [sw](const Tddd &X, double t) { return sw->gradPhi(X, t); };
        net->absorb_gradPhi_t = [sw](const Tddd &X, double t) { return sw->gradPhi_t(X, t); };
        net->absorb_phi = [sw](const Tddd &X, const double t) { return sw->phi(X, t); };
        net->absorb_eta = [sw](const Tddd &X, const double t) { return sw->eta(X, t); };
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
        net->absorb_velocity = [net](const Tddd &X, double t) { return net->water_wave_theory.gradPhi(X, t); };
        net->absorb_gradPhi_t = [net](const Tddd &X, double t) { return net->water_wave_theory.gradPhi_t(X, t); };
        net->absorb_phi = [net](const Tddd &X, const double t) { return net->water_wave_theory.phi(X, t); };
        net->absorb_eta = [net](const Tddd &X, const double t) { return net->water_wave_theory.eta(X, t); };
        net->absorb_gamma = [net](double sd) { return std::clamp(sd / (3. * net->water_wave_theory.L), 0., 1.); };
      });
    }

    //! isFixedはdefaultでfalse．指定された分だけ順に置き換わる
    //! ただし，指定が１つだけなら，それを全てのに適用する．
    if (net->inputJSON.find("isFixed")) {
      auto isFixed = stob(net->inputJSON.at("isFixed"));
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
          // std::cout << t << " " << x << " " << y << " " << z << " " << q0 << " " << q1 << " " << q2 << std::endl;
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
