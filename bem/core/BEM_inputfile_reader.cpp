#include "BEM_inputfile_reader.hpp"

#include <algorithm>
#include <cctype>

std::string toLowerCase(const std::string& str) {
  std::string strLower = str;
  std::transform(strLower.begin(), strLower.end(), strLower.begin(), [](unsigned char c) { return std::tolower(c); });
  return strLower;
}

SimulationSettings::SimulationSettings(std::filesystem::path input_directory, SimulationSettings::DomainMode mode_in)
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

    auto pick_key = [&](std::initializer_list<const char*> keys) -> std::optional<std::string> {
      for (const auto* k : keys) {
        if (settingJSON.find(k))
          return std::string(k);
      }
      return std::nullopt;
    };
    auto require_key = [&](std::initializer_list<const char*> keys, const std::string& label) -> std::string {
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
    if (const char* out_env = std::getenv("BEM_OUTPUT_DIR")) {
      if (out_env[0] != '\0') {
        common.output_directory = std::filesystem::path(out_env);
        std::cout << Yellow << "Override output_directory via BEM_OUTPUT_DIR: " << common.output_directory << colorReset << std::endl;
        std::filesystem::create_directories(common.output_directory);
      }
    }
    (void)input_files_key;

    // Helper: parse node_relocation from a string array
    // ["method", period, "surface", "midpoint_mode"]
    auto parse_node_relocation = [&](const std::vector<std::string>& arr) {
      using M = TimeDomainSettings::NodeRelocationSettings::Method;
      using S = TimeDomainSettings::NodeRelocationSettings::Surface;
      using MM = TimeDomainSettings::NodeRelocationSettings::MidpointMode;
      auto& nr = time.node_relocation;
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
        else
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                              "node_relocation midpoint_mode must be nearest (got: " + arr[3] + ")");
      }
    };

    // Helper: parse legacy ALE/ALEPERIOD keys into node_relocation
    auto parse_legacy_ale = [&](const std::string& ale_value, int period) {
      using M = TimeDomainSettings::NodeRelocationSettings::Method;
      using S = TimeDomainSettings::NodeRelocationSettings::Surface;
      auto& nr = time.node_relocation;
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

      for (auto& [key, value] : settingJSON())
        std::cout << key << ": " << value << std::endl;

      // Required key validation already handled above (mode-aware).

      /* ----------------------------- meshing options ---------------------------- */

      settingJSON.find("meshing_options", [&](const auto& STR_VEC) {
        for (const auto& STR_VEC : STR_VEC) {
          if (STR_VEC == "tetrahedralize")
            remeshing.tetrahedralize = true;
          if (STR_VEC == "surface_flip")
            remeshing.surface_flip = true;
          if (STR_VEC == "surface_split")
            remeshing.surface_split = true;
          if (STR_VEC == "surface_collapse")
            remeshing.surface_collapse = true;
          if (STR_VEC == "surface_smoothing")
            remeshing.surface_smoothing = true;
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
        } else if (elem_str.contains("hybrid")) {
          bem.element.linear = false;
          bem.element.quadratic_linear_hybrid = true;
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
      if (bem.element.quadratic_linear_hybrid)
        std::cout << "QUADRATIC_LINEAR_HYBRID_ELEMENT" << std::endl;

      {
        using M = TimeDomainSettings::NodeRelocationSettings::Method;
        using S = TimeDomainSettings::NodeRelocationSettings::Surface;
        using MM = TimeDomainSettings::NodeRelocationSettings::MidpointMode;
        auto& nr = time.node_relocation;
        if (nr.method == M::none)
          std::cout << "NODE_RELOCATION: none (pure Lagrangian)" << std::endl;
        else {
          std::cout << "NODE_RELOCATION: " << (nr.method == M::ALE ? "ALE" : "interpolation")
                    << ", period=" << nr.period;
          if (nr.method == M::interpolation)
            std::cout << ", midpoint_mode=nearest";
          std::cout << std::endl;
        }
      }

      settingJSON.find("node_relocation_debug", [&](const auto& v) {
        time.node_relocation.debug_output = stob(v)[0];
      });
      std::cout << "NODE_RELOCATION_DEBUG: " << (time.node_relocation.debug_output ? "true" : "false") << std::endl;

      // Runge-Kutta order (1=Euler, 2=Heun, 3=RK3, 4=RK4)
      settingJSON.find("rk_order", [&](const auto& v) {
        int order = std::stoi(v[0]);
        if (order < 1 || order > 4)
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "rk_order must be 1, 2, 3, or 4 (got " + std::to_string(order) + ")");
        time.rk_order = order;
      });
      std::cout << "RK_ORDER: " << time.rk_order << std::endl;

      std::cout << "" << std::endl;

      // Accept both uppercase (canonical) and lowercase (legacy) keys
      {
        const char* wd_key = settingJSON.find("WATER_DENSITY") ? "WATER_DENSITY" : (settingJSON.find("water_density") ? "water_density" : nullptr);
        if (wd_key && settingJSON.find(wd_key, [&](auto STR_VEC) {
              common.water_density = stod(STR_VEC[0]);
              _WATER_DENSITY_ = common.water_density;
            }))
          std::cout << "WATER_DENSITY: " << common.water_density << std::endl;
      }
      {
        const char* g_key = settingJSON.find("GRAVITY") ? "GRAVITY" : (settingJSON.find("gravity") ? "gravity" : nullptr);
        if (g_key && settingJSON.find(g_key, [&](auto STR_VEC) {
              common.gravity = stod(STR_VEC[0]);
              _GRAVITY_ = common.gravity;
            }))
          std::cout << "GRAVITY: " << common.gravity << std::endl;
      }

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

      auto set_ilu_neighborhood = [&](std::string type, const std::optional<std::string>& num_str = std::nullopt) {
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
      settingJSON.find("solver", [&](const auto& v) {
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
      settingJSON.find("solver_tol", [&](const auto& v) { bem.solver.solver_tol = std::stod(v[0]); });
      settingJSON.find("solver_max_iter", [&](const auto& v) { bem.solver.solver_max_iter = std::stoi(v[0]); });
      settingJSON.find("solver_restart", [&](const auto& v) { bem.solver.solver_restart = std::stoi(v[0]); });
      settingJSON.find("preconditioner", [&](const auto& v) { set_preconditioner(v[0]); });
      settingJSON.find("milu_omega", [&](const auto& v) { bem.solver.milu_omega = std::stod(v[0]); });
      settingJSON.find("ilut_drop_tol", [&](const auto& v) { bem.solver.ilut_drop_tol = std::stod(v[0]); });
      settingJSON.find("ilut_max_entries_per_row", [&](const auto& v) { bem.solver.ilut_max_entries_per_row = std::stoi(v[0]); });
      settingJSON.find("ilut_pivot_min", [&](const auto& v) { bem.solver.ilut_pivot_min = std::stod(v[0]); });
      settingJSON.find("schwarz_core_k", [&](const auto& v) { bem.solver.schwarz_core_k = std::stoi(v[0]); });
      settingJSON.find("schwarz_overlap_k", [&](const auto& v) { bem.solver.schwarz_overlap_k = std::stoi(v[0]); });
      settingJSON.find("schwarz_max_core_size", [&](const auto& v) { bem.solver.schwarz_max_core_size = std::stoi(v[0]); });
      settingJSON.find("schwarz_max_block_size", [&](const auto& v) { bem.solver.schwarz_max_block_size = std::stoi(v[0]); });
      settingJSON.find("schwarz_pivot_min", [&](const auto& v) { bem.solver.schwarz_pivot_min = std::stod(v[0]); });
      settingJSON.find("schwarz_diag_shift", [&](const auto& v) { bem.solver.schwarz_diag_shift = std::stod(v[0]); });
      settingJSON.find("coupling", [&](const auto& v) {
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
      settingJSON.find("use_metal_m2l", [&](const auto& v) { bem.solver.use_metal_m2l = stob(v)[0]; });
      settingJSON.find("metal_m2l_threadgroup", [&](const auto& v) { bem.solver.metal_m2l_threadgroup = stob(v)[0]; });
      settingJSON.find("metal_m2l_sort_terms", [&](const auto& v) { bem.solver.metal_m2l_sort_terms = stob(v)[0]; });

      // P2M quadrature points (Dunavant rule)
      settingJSON.find("p2m_quadrature_points", [&](const auto& v) { bem.solver.p2m_quadrature_points = std::stoi(v[0]); });

      // MAC criterion parameter
      settingJSON.find("mac_theta", [&](const auto& v) { bem.solver.mac_theta = std::stod(v[0]); });

      // Nearfield integration mode
      settingJSON.find("nearfield_mode", [&](const auto& v) { bem.solver.nearfield_mode = v[0]; });

      // FMM tree structure settings (requires GMRES + FMM)
      settingJSON.find("fmm_max_level", [&](const auto& v) { bem.solver.fmm_max_level = stoi(v)[0]; });
      settingJSON.find("fmm_bucket_max_points", [&](const auto& v) { bem.solver.fmm_bucket_max_points = stoi(v)[0]; });

      // Pressure-based detachment
      settingJSON.find("pressure_detachment_enabled", [&](const auto& v) { bem.solver.enable_pressure_detachment = stob(v)[0]; });
      settingJSON.find("pressure_detachment_threshold", [&](const auto& v) { bem.solver.detachment_pressure_threshold = std::stod(v[0]); });
      settingJSON.find("pressure_detachment_consecutive_steps", [&](const auto& v) { bem.solver.detachment_consecutive_steps = std::stoi(v[0]); });

      // Remeshing schedule / controls
      settingJSON.find("stop_remesh_time", [&](auto STR_VEC) { remeshing.stop_remesh_time = stod(STR_VEC[0]); });
      settingJSON.find("force_remesh_time", [&](auto STR_VEC) { remeshing.force_remesh_time = stod(STR_VEC[0]); });
      settingJSON.find("grid_refinement", [&](auto STR_VEC) { remeshing.grid_refinement = stoi(STR_VEC[0]); });
      settingJSON.find("min_edge_length", [&](const auto& v) { remeshing.min_edge_length = std::stod(v[0]); });
      settingJSON.find("max_edge_length", [&](const auto& v) { remeshing.max_edge_length = std::stod(v[0]); });

      // Initial mesh pre-relaxation (ALE-style), before the first remesh/collapse.
      // settings.json spec: "initial_mesh_pre_relax": ["true", "<loop:int>", "<coef:double>"]
      settingJSON.find("initial_mesh_pre_relax", [&](const auto& v) {
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
      settingJSON.find("collision_enabled", [&](const auto& v) { remeshing.collision.enabled = stob(v)[0]; });
      settingJSON.find("collision_proximity_factor", [&](const auto& v) { remeshing.collision.proximity_factor = std::stod(v[0]); });
      settingJSON.find("collision_normal_reversal_cos", [&](const auto& v) { remeshing.collision.normal_reversal_cos = std::stod(v[0]); });
      settingJSON.find("collision_min_zone_faces", [&](const auto& v) { remeshing.collision.min_zone_faces = std::stoi(v[0]); });
      settingJSON.find("collision_detect_folding", [&](const auto& v) { remeshing.collision.detect_folding = stob(v)[0]; });
      settingJSON.find("collision_detect_non_adjacent", [&](const auto& v) { remeshing.collision.detect_non_adjacent = stob(v)[0]; });
      settingJSON.find("collision_resolve_with_tetgen", [&](const auto& v) { remeshing.collision.resolve_with_tetgen = stob(v)[0]; });
      settingJSON.find("subsurface_altitude_reject", [&](const auto& v) { remeshing.subsurface_altitude_reject.enabled = stob(v)[0]; });
      settingJSON.find("subsurface_min_face_altitude_rel", [&](const auto& v) { remeshing.subsurface_altitude_reject.min_face_altitude_rel = std::stod(v[0]); });
      settingJSON.find("subsurface_min_edge_angle_deg", [&](const auto& v) { remeshing.subsurface_altitude_reject.min_edge_angle_deg = std::stod(v[0]); });
      settingJSON.find("subsurface_max_edge_angle_deg", [&](const auto& v) { remeshing.subsurface_altitude_reject.max_edge_angle_deg = std::stod(v[0]); });
      settingJSON.find("shell_visualization", [&](const auto& v) { remeshing.shell_visualization = stob(v)[0]; });
      settingJSON.find("front_advancing_debug", [&](const auto& v) { remeshing.front_advancing_debug = stob(v)[0]; });

      // Remesh numerical parameters
      // A. 辺長制御
      settingJSON.find("remesh_len_target_divisor",     [&](const auto& v) { remeshing.len_target_divisor = std::stod(v[0]); });
      settingJSON.find("remesh_len_fs_split_ratio",     [&](const auto& v) { remeshing.len_fs_split_ratio = std::stod(v[0]); });
      settingJSON.find("remesh_len_fs_collapse_ratio",  [&](const auto& v) { remeshing.len_fs_collapse_ratio = std::stod(v[0]); });
      settingJSON.find("remesh_len_global_max_ratio",   [&](const auto& v) { remeshing.len_global_max_ratio = std::stod(v[0]); });
      settingJSON.find("remesh_len_global_min_ratio",   [&](const auto& v) { remeshing.len_global_min_ratio = std::stod(v[0]); });
      // B. 曲率忠実度
      settingJSON.find("remesh_theta_enabled",          [&](const auto& v) { remeshing.theta_enabled = stob(v)[0]; });
      settingJSON.find("remesh_theta_target_N",         [&](const auto& v) { remeshing.theta_target_N = std::stod(v[0]); });
      settingJSON.find("remesh_theta_split_ratio",      [&](const auto& v) { remeshing.theta_split_ratio = std::stod(v[0]); });
      settingJSON.find("remesh_theta_collapse_ratio",   [&](const auto& v) { remeshing.theta_collapse_ratio = std::stod(v[0]); });
      // C. スムージング
      settingJSON.find("remesh_feature_angle_deg",      [&](const auto& v) { remeshing.feature_angle_deg = std::stod(v[0]); });
      // D. ループ制御
      settingJSON.find("remesh_max_splits",             [&](const auto& v) { remeshing.max_splits_per_step = std::stoi(v[0]); });
      settingJSON.find("remesh_max_collapses",          [&](const auto& v) { remeshing.max_collapses_per_step = std::stoi(v[0]); });
      settingJSON.find("remesh_max_smoothing",          [&](const auto& v) { remeshing.max_smoothing_per_step = std::stoi(v[0]); });
      settingJSON.find("remesh_iter_split_collapse",    [&](const auto& v) { remeshing.iter_split_collapse = std::stoi(v[0]); });
      // シナリオリスト: JSON 配列（["P","PS"]）またはカンマ区切り文字列（"P,PS"）のどちらでも受理
      auto parse_scenarios = [](const auto& v, std::vector<std::string>& out) {
        out.clear();
        for (const auto& s : v) {
          std::string buf;
          for (char c : s) {
            if (c == ',') {
              auto t = StringTrim(buf, {" "});
              if (!t.empty()) out.push_back(t);
              buf.clear();
            } else {
              buf += c;
            }
          }
          auto t = StringTrim(buf, {" "});
          if (!t.empty()) out.push_back(t);
        }
      };
      settingJSON.find("remesh_split_scenarios",        [&](const auto& v) { parse_scenarios(v, remeshing.split_scenarios); });
      settingJSON.find("remesh_collapse_scenarios",     [&](const auto& v) { parse_scenarios(v, remeshing.collapse_scenarios); });
      settingJSON.find("remesh_smoothing_scenarios",    [&](const auto& v) { parse_scenarios(v, remeshing.smoothing_scenarios); });
      // Diagnostics
      {
        auto print_list = [](const char* label, const std::vector<std::string>& v) {
          std::cout << "[remesh] " << label << " =";
          for (const auto& s : v) std::cout << " " << s;
          std::cout << std::endl;
        };
        print_list("split_scenarios",     remeshing.split_scenarios);
        print_list("collapse_scenarios",  remeshing.collapse_scenarios);
        print_list("smoothing_scenarios", remeshing.smoothing_scenarios);
      }
      // E. 採用判定 (quality gate)
      settingJSON.find("remesh_quality_hd_diag_ratio",         [&](const auto& v) { remeshing.quality_hd_diag_ratio = std::stod(v[0]); });
      settingJSON.find("remesh_quality_hd_curv_ratio",         [&](const auto& v) { remeshing.quality_hd_curv_ratio = std::stod(v[0]); });
      settingJSON.find("remesh_quality_normal_flip_cos",       [&](const auto& v) { remeshing.quality_normal_flip_cos = std::stod(v[0]); });
      settingJSON.find("remesh_quality_score_improve_margin",  [&](const auto& v) { remeshing.quality_score_improve_margin = std::stod(v[0]); });
      settingJSON.find("remesh_quality_resolution_weight",     [&](const auto& v) { remeshing.quality_resolution_weight = std::stod(v[0]); });

      // VPM (optional)
      settingJSON.find("use_VPM", [&](const auto& v) { vpm.enabled = stob(v)[0]; });
      settingJSON.find("VPM_wall_min_absorb_receivers", [&](const auto& v) { vpm.wall_min_absorb_receivers = static_cast<std::size_t>(std::stoll(v[0])); });
      settingJSON.find("VPM_wall_min_absorb_total_weight", [&](const auto& v) { vpm.wall_min_absorb_total_weight = std::stod(v[0]); });
      settingJSON.find("VPM_sigma_factor", [&](const auto& v) { vpm.sigma_factor = std::stod(v[0]); });
      settingJSON.find("VPM_stretching_scheme", [&](const auto& v) {
        if (!v.empty())
          vpm.stretching_scheme = v[0];
      });
      settingJSON.find("VPM_PSE_correction", [&](const auto& v) {
        if (!v.empty())
          vpm.PSE_correction = v[0];
      });

      // Checkpoint settings
      settingJSON.find("checkpoint_interval", [&](const auto& v) { checkpoint.interval = std::stoi(v[0]); });
      settingJSON.find("restart_from_checkpoint", [&](const auto& v) { if (!v.empty()) checkpoint.restart_from = v[0]; });
      settingJSON.find("checkpoint_max_keep", [&](const auto& v) { checkpoint.max_keep = std::stoi(v[0]); });

      // Frequency-domain optional settings (parsed even in time mode; used by main_freq_domain).
      auto parse_double_list = [&](const std::vector<std::string>& v, std::vector<double>& out) {
        out.clear();
        out.reserve(v.size());
        for (const auto& s : v)
          out.push_back(std::stod(s));
      };
      auto parse_int_list = [&](const std::vector<std::string>& v, std::vector<int>& out) {
        out.clear();
        out.reserve(v.size());
        for (const auto& s : v)
          out.push_back(std::stoi(s));
      };
      auto set_omegas = [&](const std::vector<std::string>& v) { parse_double_list(v, frequency.omegas); };
      auto set_dofs = [&](const std::vector<std::string>& v) { parse_int_list(v, frequency.dofs); };
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
        for (auto& [key, value] : injson())
          std::cout << Green << key << colorReset << ": " << value << std::endl;
        auto object_name = injson.at("name")[0];

        /* --------------------------- check required keys -------------------------- */
        if (!injson.find("name"))
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"name\" key");
        if (!injson.find("type"))
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"type\" key");
        auto type = injson.at("type")[0];
        {
          bool has_objfile = injson.find("objfile");
          bool has_objfiles = injson.find("objfiles");
          if (has_objfile && has_objfiles)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + ": cannot specify both 'objfile' and 'objfiles'");
          if (has_objfiles && !type.contains("RigidBody"))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + ": 'objfiles' is only supported for RigidBody type");
          if ((type.contains("Fluid") || type.contains("Body")) && !has_objfile && !has_objfiles)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, (input_directory / input_file_name).string() + " does not have \"objfile\" or \"objfiles\" key");
        }
        if (type.contains("Measurement") || type.contains("gauge")) {
          MeasurementJSONs.emplace_back(injson);
          std::cout << "type = " << type << std::endl;
          std::cout << "skipped" << std::endl;
          continue;
        }

        /* ------------------------ create Network object --------------------------- */

        if (!injson.find("ignore") || !stob(injson["ignore"])[0]) {

          Network* net;
          if (injson.find("objfiles")) {
            // Multi-part rigid body: merge multiple obj files into one Network
            const auto& obj_paths = injson.at("objfiles");
            if (obj_paths.empty())
              throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'objfiles' is empty");

            // Component names: use component_names if given, otherwise derive from filename stems
            std::vector<std::string> comp_names;
            if (injson.find("component_names")) {
              comp_names = injson.at("component_names");
              if (comp_names.size() != obj_paths.size())
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                                    "'component_names' size (" + std::to_string(comp_names.size()) +
                                        ") must match 'objfiles' size (" + std::to_string(obj_paths.size()) + ")");
            } else {
              for (const auto& p : obj_paths)
                comp_names.push_back(std::filesystem::path(p).stem().string());
            }

            // Validate: no duplicate names
            std::unordered_set<std::string> seen_names;
            for (const auto& n : comp_names) {
              if (seen_names.contains(n))
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "duplicate component name: " + n);
              seen_names.insert(n);
            }

            net = new Network("file_name_is_not_given", object_name);
            net->component_names = comp_names;
            for (std::size_t i = 0; i < obj_paths.size(); i++) {
              std::filesystem::path objpath(obj_paths[i]);
              if (objpath.is_relative())
                objpath = common.input_directory / objpath;
              if (!std::filesystem::exists(objpath))
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "obj file not found: " + objpath.string());
              net->mergeObjFile(objpath.string(), static_cast<int>(i));
            }
          } else {
            // Single obj file (existing behavior)
            std::string objfile_str = injson.at("objfile")[0];
            std::filesystem::path objpath(objfile_str);
            if (objpath.is_relative())
              objpath = common.input_directory / objpath;
            net = new Network(objpath.string(), object_name);
          }
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

          // rotation/translation は applyTransformations (Network.cpp) で処理済み

          setTypes(net);
          std::filesystem::copy_file(common.input_directory / input_file_name, common.output_directory / input_file_name, std::filesystem::copy_options::overwrite_existing);
          mk_vtu(common.output_directory / (object_name + "_init.vtu"), {net->getFaces()});
          /* -------------------------------------------------------------------------- */
          // Phase 2 (2026-04-12): mooring input migrated from per-line raw
          // `new MooringLine` + `setEquilibriumState` to `LumpedCableSystem`.
          // The system is created lazily on the first mooring_* key, then
          // every subsequent mooring_* key adds another cable. After the loop,
          // a single `solveEquilibrium()` call relaxes all cables together
          // with the fast cable_solver-style RK4 driver (no dt warmup ramp).
          for (auto& [key, value] : injson()) {
            if (key.contains("mooring")) {
              //$ ------------------------ MOORING ---------------------- */
              std::cout << "/* ------------------------ MOORING ---------------------- */" << std::endl;
              std::cout << "initialize mooring" << std::endl;

              //! check if the value contains 13 elements
              if (value.size() != 13)
                throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "mooring line must have 13 elements");

              const auto name = value[0];
              const int n_points = std::stoi(value[8]);
              const std::array<double, 3> X_begin = {std::stod(value[1]), std::stod(value[2]), std::stod(value[3])};
              const std::array<double, 3> X_end = {std::stod(value[4]), std::stod(value[5]), std::stod(value[6])};
              const double total_length = std::stod(value[7]);
              const double w = std::stod(value[9]);
              const double stiffness = std::stod(value[10]);
              const double damp = std::stod(value[11]);
              const double diam = std::stod(value[12]);

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

              if (!net->cable_system)
                net->cable_system = std::make_unique<LumpedCableSystem>();

              CableProperties props{w, stiffness, damp, diam};

              // X_begin = world-fixed anchor (seabed/tower).
              // X_end   = fairlead, attached to this floating body. The
              //          initial world position seeds the straight-line
              //          discretisation; once attached the node's initialX
              //          will be used by CableAttachment::currentWorldPosition()
              //          for the rigid-body transform.
              CableAttachment end_a = CableAttachment::worldFixed(X_begin);
              CableAttachment end_b = CableAttachment::onBody(net, X_end);

              auto* cable = net->cable_system->addCable(name, end_a, end_b,
                                                         total_length, n_points,
                                                         props);
              setOutputInfo(name, common.output_directory);

              std::cout << std::right << std::setw(15)
                        << "LumpedCable->getPoints().size() = " << cable->getPoints().size() << std::endl;
              std::cout << std::right << std::setw(15)
                        << "LumpedCable->getName() = " << cable->getName() << std::endl;
              std::cout << "/* ------------------------------------------------------- */" << std::endl;
            }
          }
          //$ -- After all mooring_* keys parsed, run a single batch equilibrium --
          if (net->cable_system && net->cable_system->size() > 0) {
            std::cout << "[mooring] solving initial equilibrium for "
                      << net->cable_system->size() << " cable(s) on body "
                      << net->getName() << std::endl;
            net->cable_system->solveEquilibrium(/*tol=*/0.01,
                                                 /*max_steps=*/500000,
                                                 /*snapshot_interval=*/10000);
          }
          //$ ------------------------------------------------------- */
        } else
          Print("skipped");
      }
    }

    // If coupling is not explicitly specified and there are no movable bodies, disable coupling.
    if (!coupling_explicitly_set) {
      bool has_movable_body = false;
      for (auto* net : Join(RigidBodyObject, SoftBodyObject)) {
        for (const auto& fixed : net->isFixed) {
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
