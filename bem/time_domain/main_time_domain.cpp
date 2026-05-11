/*

cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=main_time_domain.cpp

*/
#include "pch.hpp"
#include <array>
#include <cctype>
#include <cstdint>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <execinfo.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include "BEM_time_domain_types.hpp"

/*DOC_EXTRACT 0_0_BEM

# BEM-MEL

<img src="./sample_Goring1979.gif">

*/

// #define _debugging_

//# ============================================================================
//# SECTION 1: Global Variables
//#   Element type flags, solver/FMM parameters, simulation state
//#   Referenced by BEM_*.hpp via extern declarations
//# ============================================================================

bool use_linear_element = false;
bool use_pseudo_quadratic_element = false;
bool use_true_quadratic_element = false;
bool use_quadratic_linear_hybrid = false;
NodeRelocationMethod node_relocation_method = NodeRelocationMethod::none;
NodeRelocationSurface node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
InterpolationMidpointMode interpolation_midpoint_mode = InterpolationMidpointMode::nearest;
std::string solver_type = "LU";
std::string coupling_type = "NONE";
double coupling_tol = 1e-10;
std::vector<double> coupling_params;
std::string preconditioner_type = "NONE";
std::string ilu_neighborhood_type = "BUCKETS";
int ilu_kring_num = 1;
double milu_omega = 1.0;
double ilut_drop_tol = 1e-3;
int ilut_max_entries_per_row = 50;
double ilut_pivot_min = 1e-12;
int schwarz_core_k = 1;
int schwarz_overlap_k = 1;
int schwarz_max_core_size = 64;
int schwarz_max_block_size = 128;
double schwarz_pivot_min = 1e-12;
double schwarz_diag_shift = 0.0;
double solver_tol = 1e-9;
int solver_max_iter = 500;
int solver_restart = 100;
int fmm_max_level = 7;
int fmm_bucket_max_points = 50;
#if defined(USE_METAL_M2L)
bool use_metal_m2l = false;
bool metal_m2l_threadgroup = false;
bool metal_m2l_sort_terms = false;
#endif
std::string nearfield_mode = "scalar";
int g_p2m_quadrature_points = 6;
int g_lu_far_dunavant_points = 0;
int g_lu_near_dunavant_points = 0;
double g_mac_theta = 0.25;

int time_step;
double simulation_time = 0;
bool enable_pressure_detachment = false;
double detachment_pressure_threshold = 0.0;
int detachment_consecutive_steps = 3;

//# ============================================================================
//# SECTION 2: Includes
//# ============================================================================

#define BEM
#define simulation
#include "Network.hpp"

JSONoutput jsonout;

#include "BEM.hpp"                  // Core BEM: BoundaryValues, calculateVelocities, setBoundaryTypes, solveBVP
#include "BEM_inputfile_reader.hpp" // JSON settings parser
#include "OutputCommon.hpp"         // Output context
#include "OutputJSON.hpp"           // JSON output writer
#include "OutputParaview.hpp"       // ParaView VTU/PVD output

#include "BEM_internal_flow.hpp" // Internal flow handling
#include "VPM.hpp"               // Vortex Particle Method

// --- Modules compiled as separate translation units (see CMakeLists.txt) ---
// BEM_checkpoint.cpp: binary checkpoint I/O
namespace BEM_Checkpoint {
void writeCheckpointToPath(const std::filesystem::path&, int, double, double,
                           const std::vector<Network*>&, const std::vector<Network*>&,
                           const std::vector<Network*>&, const VortexMethod&, bool);
void writeCheckpoint(const std::filesystem::path&, int, double, double,
                     const std::vector<Network*>&, const std::vector<Network*>&,
                     const std::vector<Network*>&, const VortexMethod&, bool);
std::tuple<int, double, double> readCheckpoint(const std::filesystem::path&,
                                               std::vector<Network*>&, const std::vector<Network*>&,
                                               const std::vector<Network*>&, VortexMethod&, bool);
void pruneCheckpoints(const std::filesystem::path&, int);
} // namespace BEM_Checkpoint

// BEM_penetration.cpp: structure penetration detection
namespace BEM_Penetration {
void throwIfStructurePenetrated(const std::vector<Network*>&, const std::vector<Network*>&,
                                int, const std::string&, bool);
} // namespace BEM_Penetration

// --- End separate TU forward declarations ---

#include "BEM_remesh_main.hpp"        // Remeshing (split/collapse/flip)
#include "BEM_remesh_logging.hpp"     // Shared remesh / mesh-preparation logs
#include "BEM_mesh_relaxation.hpp"    // Pre-relaxation (flip + smoothing at t=0)
#include "BEM_debug_helpers.hpp"      // Debug corner-point inspection, crash backtrace
#include "BEM_rk_update.hpp"          // Absorption, RK push, interpolation relocation, mean-phi
#include "BEM_initial_conditions.hpp" // Initial condition (IC) application at t=0
#include "BEM_time_step_control.hpp"  // Step retry, dt control, body state sync
#include "BEM_snapshot.hpp"           // Mesh preparation pipeline reference state
#include "BEM_field_mapper.hpp"       // Mesh preparation pipeline field mapping
#include "BEM_geometry_repair.hpp"    // Mesh preparation pipeline geometry repair
#include "BEM_topology_modifier.hpp"  // Mesh preparation pipeline topology repair

//# ============================================================================
//# SECTION 3: main()
//#
//#   Flow overview:
//#     3a. Settings & initialization
//#     3b. Initial mesh pre-relaxation (t=0 only)
//#     3c. Time loop
//#         3c-1. Output lambdas (write_step_outputs, write_failure_snapshot)
//#         3c-2. Step checkpoint & retry loop
//#             3c-2a. Remesh + quality checks
//#             3c-2b. Initial conditions (t=0 only)
//#             3c-2c. dt determination
//#             3c-2d. RK initialization
//#             3c-2e. RK loop (do-while)
//#                    - setBoundaryTypes (RK_step==1)
//#                    - Absorber assignment
//#                    - BVP: setPhiPhinOnFace -> BVP.solve -> velocity
//#                    - VPM advection (if enabled)
//#                    - Body update (pre-solve read + post-solve push)
//#                    - Absorption + RK push (-> BEM_rk_update.hpp)
//#                    - Quality checks
//#             3c-2f. Post-RK: interpolation relocation, VPM step, mean-phi
//#             3c-2g. Fixed-body pressure approximation
//#         3c-3. Step output
//# ============================================================================

int main(int argc, char** argv) {
  installCrashBacktraceIfRequested();
  enableRealFieldFmmDefaultForTimeDomain();
  std::clock_t cpu_clock_start = std::clock();
  auto wall_clock_start = std::chrono::high_resolution_clock::now();

  /*DOC_EXTRACT 0_1_BEM

  ## 入力ファイルの読み込み

  1. 境界条件の設定
  2. 境界値問題（BIE）を解き，$\phi$と$\phi_n$を求める
  3. 三角形の線形補間を使って節点の流速を計算する

  */

  //@ --- 3a. Settings & initialization ---

  if (!initializeLogFile("log.txt", argc, argv))
    return 1;
  if (argc <= 1)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "argv <= 1. write input json file directory!\nex.\n$ ./main ./input");
  SimulationSettings setting(argv[1], SimulationSettings::DomainMode::Time);
  std::filesystem::path input_directory = setting.common.input_directory;
  const double max_dt = setting.time.max_dt;
  const int node_relocation_period = setting.time.node_relocation.period;
  int end_time_step = setting.time.end_time_step;
  // BEM_BENCH_STEPS: override end_time_step for short profiling runs.
  if (const char* bench_steps_env = std::getenv("BEM_BENCH_STEPS")) {
    const int bench_steps = std::atoi(bench_steps_env);
    if (bench_steps > 0 && bench_steps < end_time_step)
      end_time_step = bench_steps;
  }
  const double end_time = setting.time.end_time;
  const double stop_remesh_time = setting.remeshing.stop_remesh_time;
  const double force_remesh_time = setting.remeshing.force_remesh_time;
  const bool tetrahedralize = setting.remeshing.tetrahedralize;
  const bool surface_flip = setting.remeshing.surface_flip;
  const int grid_refinement = setting.remeshing.grid_refinement;
  const auto mesh_pipeline = setting.time.mesh_preparation_pipeline;
  const bool use_mesh_preparation_pipeline = mesh_pipeline.enabled;
  const std::filesystem::path output_directory = setting.common.output_directory;
  use_linear_element = setting.bem.element.linear;
  use_pseudo_quadratic_element = setting.bem.element.pseudo_quadratic;
  use_quadratic_linear_hybrid = setting.bem.element.quadratic_linear_hybrid;
  // hybrid時もtrue_quadraticをtrueに — ソルバーの中点処理ゲートを有効化
  use_true_quadratic_element = setting.bem.element.true_quadratic || use_quadratic_linear_hybrid;
  // Node relocation: method
  {
    using M = SimulationSettings::TimeDomainSettings::NodeRelocationSettings::Method;
    switch (setting.time.node_relocation.method) {
    case M::ALE:
      node_relocation_method = NodeRelocationMethod::ALE;
      break;
    case M::interpolation:
      node_relocation_method = NodeRelocationMethod::interpolation;
      break;
    default:
      node_relocation_method = NodeRelocationMethod::none;
      break;
    }
  }
  if (use_mesh_preparation_pipeline && node_relocation_method != NodeRelocationMethod::none) {
    std::cout << Yellow
              << "[mesh:config] mesh preparation pipeline is enabled; legacy node relocation is disabled"
              << colorReset << std::endl;
    node_relocation_method = NodeRelocationMethod::none;
  }
  {
    interpolation_midpoint_mode = InterpolationMidpointMode::nearest;
  }
  // Node relocation: surface precision (auto-resolve from element type unless explicitly set)
  if (setting.time.node_relocation.surface_explicitly_set) {
    using S = SimulationSettings::TimeDomainSettings::NodeRelocationSettings::Surface;
    switch (setting.time.node_relocation.surface) {
    case S::linear:
      node_relocation_surface = NodeRelocationSurface::linear;
      break;
    case S::true_quadratic:
      node_relocation_surface = NodeRelocationSurface::true_quadratic;
      break;
    default:
      node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
      break;
    }
  } else {
    if (use_true_quadratic_element)
      node_relocation_surface = NodeRelocationSurface::true_quadratic;
    else if (use_mesh_preparation_pipeline)
      node_relocation_surface = NodeRelocationSurface::linear;
    else
      node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
  }
  if (use_mesh_preparation_pipeline && !setting.time.node_relocation.surface_explicitly_set)
    std::cout << Green << "[mesh:config] node_relocation_surface=auto phase3_interpolation_surface="
              << (node_relocation_surface == NodeRelocationSurface::true_quadratic ? "true_quadratic" : "linear")
              << colorReset << std::endl;
  solver_type = setting.bem.solver.solver_type;
  coupling_type = setting.bem.solver.coupling_type;
  coupling_tol = setting.bem.solver.coupling_tol;
  coupling_params = setting.bem.solver.coupling_params;
  preconditioner_type = setting.bem.solver.preconditioner_type;
  ilu_neighborhood_type = setting.bem.solver.ilu_neighborhood_type;
  ilu_kring_num = setting.bem.solver.ilu_kring_num;
  milu_omega = setting.bem.solver.milu_omega;
  ilut_drop_tol = setting.bem.solver.ilut_drop_tol;
  ilut_max_entries_per_row = setting.bem.solver.ilut_max_entries_per_row;
  ilut_pivot_min = setting.bem.solver.ilut_pivot_min;
  schwarz_core_k = setting.bem.solver.schwarz_core_k;
  schwarz_overlap_k = setting.bem.solver.schwarz_overlap_k;
  schwarz_max_core_size = setting.bem.solver.schwarz_max_core_size;
  schwarz_max_block_size = setting.bem.solver.schwarz_max_block_size;
  schwarz_pivot_min = setting.bem.solver.schwarz_pivot_min;
  schwarz_diag_shift = setting.bem.solver.schwarz_diag_shift;
  solver_tol = setting.bem.solver.solver_tol;
  solver_max_iter = setting.bem.solver.solver_max_iter;
  solver_restart = setting.bem.solver.solver_restart;
  fmm_max_level = setting.bem.solver.fmm_max_level;
  fmm_bucket_max_points = setting.bem.solver.fmm_bucket_max_points;
#if defined(USE_METAL_M2L)
  use_metal_m2l = setting.bem.solver.use_metal_m2l;
  metal_m2l_threadgroup = setting.bem.solver.metal_m2l_threadgroup;
  metal_m2l_sort_terms = setting.bem.solver.metal_m2l_sort_terms;
#endif
  nearfield_mode = setting.bem.solver.nearfield_mode;
  g_p2m_quadrature_points = setting.bem.solver.p2m_quadrature_points;
  g_mac_theta = setting.bem.solver.mac_theta;
  enable_pressure_detachment = setting.bem.solver.enable_pressure_detachment;
  detachment_pressure_threshold = setting.bem.solver.detachment_pressure_threshold;
  detachment_consecutive_steps = setting.bem.solver.detachment_consecutive_steps;

  std::map<std::string, outputInfo> NetOutputInfo = setting.NetOutputInfo;
  const bool shell_visualization = setting.remeshing.shell_visualization;
  const bool front_advancing_debug = setting.remeshing.front_advancing_debug;
  const bool node_relocation_debug_output = setting.time.node_relocation.debug_output;
  const auto subsurface_altitude_reject = setting.remeshing.subsurface_altitude_reject;
  const bool initial_mesh_pre_relax = setting.remeshing.initial_mesh_pre_relax.enabled;
  const int initial_mesh_pre_relax_loop = setting.remeshing.initial_mesh_pre_relax.loop;
  const double initial_mesh_pre_relax_coef = setting.remeshing.initial_mesh_pre_relax.coef;
  const double min_edge_length = setting.remeshing.min_edge_length;
  const int checkpoint_interval = setting.checkpoint.interval;
  const std::string restart_from_checkpoint = setting.checkpoint.restart_from;
  const int checkpoint_max_keep = setting.checkpoint.max_keep;
  std::vector<Network*> FluidObject = setting.FluidObject;
  std::vector<Network*> RigidBodyObject = setting.RigidBodyObject;
  std::vector<Network*> SoftBodyObject = setting.SoftBodyObject;
  std::vector<Network*> AbsorberObject = setting.AbsorberObject;
  std::vector<JSON> MeasurementJSONs = setting.MeasurementJSONs;
#ifdef USE_TETGEN
  if (tetrahedralize)
    for (auto& network : FluidObject)
      network->tetrahedralize();
#endif

  //$ --- Create output directory and copy source files ---
  std::filesystem::create_directories(output_directory);
  std::filesystem::copy_file(setting.settings_file_path, output_directory / "settings.json", std::filesystem::copy_options::overwrite_existing);

  // Copy source files to output directory for reproducibility.
  // Resolve source directory from __FILE__ (works regardless of working directory).
  const std::filesystem::path source_dir = std::filesystem::path(__FILE__).parent_path();
  auto safe_copy = [&](const std::filesystem::path& src) {
    try {
      if (std::filesystem::exists(src))
        std::filesystem::copy_file(src, output_directory / src.filename(), std::filesystem::copy_options::overwrite_existing);
    } catch (const std::filesystem::filesystem_error& e) {
      std::cerr << "[source-copy] skip " << src << ": " << e.what() << std::endl;
    }
  };
  safe_copy(source_dir / "main_time_domain.cpp");
  safe_copy(source_dir / "main.cpp");
  std::regex pattern("^BEM.*\\.hpp$");
  try {
    for (auto& entry : std::filesystem::directory_iterator(source_dir))
      if (std::regex_match(entry.path().filename().string(), pattern))
        safe_copy(entry.path());
  } catch (const std::filesystem::filesystem_error& e) {
    std::cerr << "[source-copy] skip directory " << source_dir << ": " << e.what() << std::endl;
  }

  PVDWriter cornerPointsPVD(output_directory / "cornerPointsPVD.pvd");
  PVDWriter DirichletSurfacePVD(output_directory / "DirichletSurface.pvd");
  PVDWriter vpm_pvd(output_directory / "vpm.pvd");
  const std::string water_name = FluidObject.empty() ? "" : FluidObject[0]->getName() + "_";
  PVDWriter candidate_patches_pvd((output_directory / (water_name + "candidate_patches.pvd")).string());
  PVDWriter remeshed_patches_pvd((output_directory / (water_name + "remeshed_patches.pvd")).string());
  PVDWriter trigger_edges_pvd((output_directory / (water_name + "trigger_edges.pvd")).string());
  PVDWriter edges_pvd((output_directory / (water_name + "edges.pvd")).string());
  PVDWriter provider_source_of_truth_pvd(output_directory / "provider_source_of_truth_quadratic.pvd");
  const std::string phase_water_name = FluidObject.empty() ? "water" : FluidObject[0]->getName();
  std::filesystem::remove(output_directory / (phase_water_name + "_remeshed.pvd"));
  std::filesystem::remove(output_directory / (phase_water_name + "_remeshed_repaired.pvd"));
  PVDWriter water_after_phase0_pvd(output_directory / (phase_water_name + "_after_phase0.pvd"));
  PVDWriter water_after_phase1_pvd(output_directory / (phase_water_name + "_after_phase1.pvd"));
  PVDWriter water_after_phase2_pvd(output_directory / (phase_water_name + "_after_phase2.pvd"));
  Print("setting done");

  /*DOC_EXTRACT 0_1_BEM

  ## 計算プログラムの概要

  | 項目 | 詳細|
  |---:|:---|
  | 要素 | 線形三角要素 |
  | 時間発展方法 | 4次のルンゲクッタ |
  | 解析領域 | 時間領域 |
  | 境界条件 | 水面の境界条件は非線形であるが，非線形のまま解く |

  ### 計算の流れ

  1. 境界条件の設定
  2. 境界値問題（BIE）を解き，$\phi$と$\phi_n$を求める
  3. 三角形の線形補間を使って節点の流速を計算する
  4. 次時刻の$\Omega(t+\Delta t)$がわかるので，修正流速を計算する
  5.
  浮体の加速度を計算する．境界値問題（BIE）を解き，$\phi_t$と$\phi_{nt}$を求め，浮体面上の圧力$p$を計算する必要がある
  6. 全境界面の節点の位置を更新．ディリクレ境界では$\phi$を次時刻の値へ更新

  */

  try {
    //# --- 3c. Time loop (main loop) ---
    Buckets<networkFace*> Buckets_Fluid_Faces;
    Buckets<networkPoint*> Buckets_Fluid_Points;
    std::vector<BEM_DOF_Base*> fluid_nodes;

    VortexMethod vpm;
    vpm.configure(setting.vpm.enabled, setting.vpm.stretching_scheme, setting.vpm.PSE_correction,
                  setting.vpm.wall_min_absorb_receivers, setting.vpm.wall_min_absorb_total_weight,
                  setting.vpm.sigma_factor);
    double time_setup = 0.0, time_solve = 0.0, unknownsizse = 0.0;

    //$ --- Checkpoint restart ---
    int start_time_step = 0;
    double last_successful_dt = max_dt;
    if (!restart_from_checkpoint.empty()) {
      auto [rstep, rtime, rdt] =
          BEM_Checkpoint::readCheckpoint(restart_from_checkpoint, FluidObject, RigidBodyObject, SoftBodyObject, vpm, use_true_quadratic_element);
      start_time_step = rstep + 1; // resume from the NEXT step
      simulation_time = rtime;
      last_successful_dt = rdt;
      std::cout << Green << "[restart] Loaded checkpoint: step=" << rstep
                << ", time=" << rtime << ", dt=" << rdt << colorReset << std::endl;

      // Output initial state from checkpoint so it appears in ParaView/JSON output
      {
        OutputContext ctx{.dt = rdt, .time_step = rstep, .simulation_time = rtime, .output_directory = output_directory, .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};
        OutputParaView::write_step(ctx, NetOutputInfo, FluidObject, RigidBodyObject, SoftBodyObject, {});
        if (shell_visualization)
          OutputParaView::write_shell_step(ctx, NetOutputInfo, FluidObject);
        if (front_advancing_debug)
          OutputParaView::write_aft_candidates_step(ctx, NetOutputInfo, FluidObject);
        std::cout << Green << "[restart] Wrote initial state (step=" << rstep << ") to output" << colorReset << std::endl;
      }
    }

    BEM_BVP BVP(FluidObject); // もともとあってバグった場所
    BVP.output_directory = output_directory;

    for (time_step = start_time_step; time_step < end_time_step; time_step++) {
      if (end_time < simulation_time)
        break;

      double dt = 1E+20;
      int RK_order = setting.time.rk_order;

      double spacing = 0.;
      for (auto& water : FluidObject)
        spacing += Mean(extLength(water->getLines())) * 10;
      spacing /= FluidObject.size();
      std::vector<Network*> AllObjects = Join(FluidObject, RigidBodyObject, SoftBodyObject);

      //@ --- 3b. Initial mesh pre-relaxation (t=0 only) ---
      if (time_step == 0 && start_time_step == 0 && initial_mesh_pre_relax && initial_mesh_pre_relax_loop > 0 && initial_mesh_pre_relax_coef > 0.0) {
        const auto contact_objects = Join(RigidBodyObject, SoftBodyObject);
        auto write_pre_relax_mesh = [&](Network* net, const std::string& tag) {
          const std::string suffix = tag.empty() ? "" : ("_" + tag);
          if (auto it = NetOutputInfo.find(net->getName()); it != NetOutputInfo.end()) {
            const auto& info = it->second;
            std::filesystem::path filename = info.vtu_file_name + std::string("pre_relax") + suffix + ".vtu";
            OutputParaView::mk_vtu_quadratic((output_directory / filename).string(), net->getBoundaryFaces(), dataForOutput(net, /*dt=*/0.0));
            if (info.PVD) {
              info.PVD->push(filename, simulation_time);
              info.PVD->output();
            }
          } else {
            std::filesystem::path filename = net->getName() + std::string("_pre_relax") + suffix + ".vtu";
            OutputParaView::mk_vtu_quadratic((output_directory / filename).string(), net->getBoundaryFaces(), dataForOutput(net, /*dt=*/0.0));
          }
          {
            std::filesystem::path filename = net->getName() + std::string("_pre_relax") + suffix + ".obj";
            std::ofstream ofs(output_directory / filename);
            if (ofs)
              createOBJ(ofs, *net);
          }
        };
        preRelaxMesh(FluidObject, AllObjects, contact_objects,
                     use_true_quadratic_element, node_relocation_surface, simulation_time,
                     {.loop = initial_mesh_pre_relax_loop, .coef = initial_mesh_pre_relax_coef, .output_tag = "after_relax"},
                     write_pre_relax_mesh);
      }

      if (time_step == 0 && start_time_step == 0 && shell_visualization) {
        OutputContext ctx{.dt = 0.0, .time_step = 0, .simulation_time = simulation_time, .output_directory = output_directory, .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};
        OutputParaView::write_shell_step(ctx, NetOutputInfo, FluidObject);
      }
      if (time_step == 0 && start_time_step == 0 && front_advancing_debug) {
        OutputContext ctx{.dt = 0.0, .time_step = 0, .simulation_time = simulation_time, .output_directory = output_directory, .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};
        OutputParaView::write_aft_candidates_step(ctx, NetOutputInfo, FluidObject);
      }

      //$ --- 3c-1. Output lambdas & per-step state ---
      std::vector<double> convergence;
      std::vector<double> ElapsedTimeSetup, ElapsedTimeSolve;
      double ElapsedTimeNodeRelocation = 0, ElapsedTimeTotal = 0;
      StepRetryState retry_state;

      auto sync_bodies = [&]() {
        sync_body_states_from_rk(RigidBodyObject, SoftBodyObject);
      };

      auto rebuild_output_buckets = [&]() -> double {
        if (FluidObject.empty())
          return simulation_time;

        CoordinateBounds bounds_org(FluidObject[0]->bounds);
        for (auto& water : FluidObject)
          bounds_org += water->bounds;
        CoordinateBounds bounds = bounds_org.scaledBounds(1.5);
        Buckets_Fluid_Faces.initialize(bounds, bounds.getScale() / 10.);
        Buckets_Fluid_Points.initialize(bounds, bounds.getScale() / 10.);
        for (auto& water : FluidObject) {
          for (const auto& f : water->getBoundaryFaces())
            Buckets_Fluid_Faces.add(f->getXtuple(), f);
          for (const auto& p : water->getPoints())
            Buckets_Fluid_Points.add(ToX(p), p);
        }

        double snapshot_time = simulation_time;
        if (!Buckets_Fluid_Points.data1D.empty())
          snapshot_time = (*(Buckets_Fluid_Points.data1D).begin())->RK_X.gett();
        return snapshot_time;
      };

      auto write_step_outputs = [&](const bool force_checkpoint) {
        if (Network::scalingDebugEnabled()) {
          std::cout << Magenta << "[scale] " << Cyan << "before write_step_outputs" << colorReset << std::endl;
          for (auto* net : FluidObject)
            net->debugPrintScalingState("time_loop write_step_outputs");
        }
        sync_bodies();
        const double snapshot_time = rebuild_output_buckets();

        for (auto water : FluidObject) {
          auto name = water->getName();
          std::ofstream ofs(output_directory / (name + ".obj"));
          createOBJ(ofs, *water);
          ofs.close();
        }

        OutputContext ctx{.dt = dt, .time_step = time_step, .simulation_time = snapshot_time, .output_directory = output_directory, .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};

        OutputJSON::write_step(ctx, jsonout, FluidObject, RigidBodyObject, SoftBodyObject, MeasurementJSONs, Buckets_Fluid_Faces.data1D, convergence, unknownsizse, time_setup, time_solve, BVP.last_ilu_build_time, BVP.last_ilu_apply_time_sum, BVP.last_gmres_iter_time_sum, BVP.last_A_sparse_nnz, BVP.last_A_sparse_avg_nnz);
        OutputParaView::write_step(ctx, NetOutputInfo, FluidObject, RigidBodyObject, SoftBodyObject, Buckets_Fluid_Faces.data1D);
        if (shell_visualization)
          OutputParaView::write_shell_step(ctx, NetOutputInfo, FluidObject);
        if (front_advancing_debug)
          OutputParaView::write_aft_candidates_step(ctx, NetOutputInfo, FluidObject);
        if (vpm.isEnabled())
          OutputParaView::write_vpm(ctx, vpm, vpm_pvd);

        if (force_checkpoint || (checkpoint_interval > 0 && time_step % checkpoint_interval == 0)) {
          BEM_Checkpoint::writeCheckpoint(output_directory, time_step, snapshot_time, dt,
                                          FluidObject, RigidBodyObject, SoftBodyObject, vpm,
                                          use_true_quadratic_element);
          if (!force_checkpoint && checkpoint_max_keep > 0)
            BEM_Checkpoint::pruneCheckpoints(output_directory, checkpoint_max_keep);
        }
      };

      auto write_failure_snapshot = [&](const int rk_step, const int step_retry) {
        if (Network::scalingDebugEnabled()) {
          std::cout << Magenta << "[scale] " << Cyan << "before write_failure_snapshot" << colorReset << std::endl;
          for (auto* net : FluidObject)
            net->debugPrintScalingState("time_loop write_failure_snapshot");
        }
        sync_bodies();
        const double snapshot_time = rebuild_output_buckets();
        const std::string tag = "failure_t" + std::to_string(time_step) +
                                "_rk" + std::to_string(rk_step) +
                                "_retry" + std::to_string(step_retry);

        for (auto* net : FluidObject) {
          auto it = NetOutputInfo.find(net->getName());
          const std::string prefix = (it != NetOutputInfo.end()) ? it->second.vtu_file_name : (net->getName() + "_");
          const auto filename = output_directory / (prefix + tag + ".vtu");
          OutputParaView::mk_vtu_quadratic(filename.string(), net->getBoundaryFaces(), dataForOutput(net, dt));

          auto tet_it = NetOutputInfo.find(net->getName() + "_tetra");
          if (tet_it != NetOutputInfo.end()) {
            const auto tet_filename = output_directory / (tet_it->second.vtu_file_name + tag + ".vtu");
            std::ofstream ofs(tet_filename);
            vtkUnstructuredGridWrite(ofs, net->getTetras());
          }
        }

        for (auto* net : RigidBodyObject) {
          auto it = NetOutputInfo.find(net->getName());
          const std::string prefix = (it != NetOutputInfo.end()) ? it->second.vtu_file_name : (net->getName() + "_");
          const auto filename = output_directory / (prefix + tag + ".vtu");
          OutputParaView::mk_vtu_quadratic(filename.string(), net->getBoundaryFaces(), dataForOutput(net, dt));
        }

        for (auto* net : SoftBodyObject) {
          auto it = NetOutputInfo.find(net->getName());
          const std::string prefix = (it != NetOutputInfo.end()) ? it->second.vtu_file_name : (net->getName() + "_");
          const auto filename = output_directory / (prefix + tag + ".vtu");
          OutputParaView::mk_vtu_quadratic(filename.string(), net->getBoundaryFaces(), dataForOutput(net, dt));
        }

        if (vpm.isEnabled())
          OutputParaView::write_vpm_vtp(output_directory / ("vpm_" + tag + ".vtp"), vpm);

        const auto checkpoint_path = output_directory / ("checkpoint_" + tag + ".bin");
        BEM_Checkpoint::writeCheckpointToPath(checkpoint_path, time_step, snapshot_time, dt,
                                              FluidObject, RigidBodyObject, SoftBodyObject, vpm,
                                              use_true_quadratic_element);
      };

      //@ --- 3c-2. Step checkpoint & retry loop ---
      auto step_checkpoint_path = output_directory / "_step_retry.bin";

      auto save_step_state = [&]() {
        sync_bodies();
        BEM_Checkpoint::writeCheckpointToPath(step_checkpoint_path, time_step, simulation_time, last_successful_dt,
                                              FluidObject, RigidBodyObject, SoftBodyObject, vpm,
                                              use_true_quadratic_element);
      };

      auto restore_step_state = [&](int step_retry) {
        std::cout << Red << "[step_reject] retry " << step_retry << "/" << StepRetryState::max_step_retries
                  << " with dt_override=" << retry_state.dt_override << colorReset << std::endl;
        auto [rstep, rtime, rdt] =
            BEM_Checkpoint::readCheckpoint(step_checkpoint_path, FluidObject, RigidBodyObject, SoftBodyObject,
                                           vpm, use_true_quadratic_element);
        simulation_time = rtime;
        BVP.WATERS = FluidObject;
        BVP.ilu_preconditioner.reset();
        BVP.ilut_preconditioner.reset();
        BVP.schwarz_preconditioner.reset();
        BVP.ilu_topology_hash_cached = 0;
        BVP.output_directory = output_directory;
        AllObjects = Join(FluidObject, RigidBodyObject, SoftBodyObject);
        if (Network::scalingDebugEnabled()) {
          std::cout << Magenta << "[scale] " << Cyan << "after restore_step_state" << colorReset << std::endl;
          for (auto* net : FluidObject)
            net->debugPrintScalingState("time_loop restore_step_state");
        }
      };

      // Overwrite step checkpoint with post-remesh state.
      // On BVP/RK failure, retry will restore this state and skip remesh,
      // avoiding topology changes and interpolation noise from re-remeshing.
      // Note: do NOT call sync_bodies() here — RK is not yet
      // initialized for this step, so reading RK state would corrupt body positions.
      auto save_post_remesh_state = [&]() {
        BEM_Checkpoint::writeCheckpointToPath(step_checkpoint_path, time_step, simulation_time, last_successful_dt,
                                              FluidObject, RigidBodyObject, SoftBodyObject, vpm,
                                              use_true_quadratic_element);
      };

      auto dump_pipeline_mesh = [&](const std::string& tag, bool enabled) {
        if (!enabled)
          return;
        for (auto* water : FluidObject) {
          if (!water)
            continue;
          const auto filename = output_directory / (water->getName() + "_" + tag + "_" + std::to_string(time_step) + ".vtu");
          OutputParaView::mk_vtu_quadratic(filename.string(), water->getBoundaryFaces(), dataForOutput(water, /*dt=*/0.0));
        }
      };

      auto write_water_phase_snapshot = [&](Network* water, int step_retry,
                                            const std::string& phase_tag,
                                            PVDWriter& pvd) {
        if (!water)
          return;
        const std::string filename =
            water->getName() + "_" + phase_tag + "_t" + std::to_string(time_step) +
            "_r" + std::to_string(step_retry) + ".vtu";
        OutputParaView::mk_vtu_quadratic(
            (output_directory / filename).string(),
            water->getBoundaryFaces(),
            dataForOutput(water, /*dt=*/0.0));
        pvd.push(filename, simulation_time);
        pvd.output();
      };
      auto write_all_water_phase_snapshots = [&](const std::string& phase_tag,
                                                 int step_retry,
                                                 PVDWriter& pvd) {
        for (auto* water : FluidObject)
          write_water_phase_snapshot(water, step_retry, phase_tag, pvd);
      };

      struct ProviderSourceOfTruthQuadraticCell {
        std::array<Tddd, 6> nodes;
        int reference_face_index = -1;
        int body_object_index = -1;
        int body_face_index = -1;
        int source_type = 0; // 1: snapshot reference, 2: initial-condition reference, 3: body target
        int bc_type = static_cast<int>(NodeFaceBoundaryType::Dirichlet);
      };

      auto provider_cell_midpoint_deviation = [](const std::array<Tddd, 6>& nodes) {
        return std::max({Norm(nodes[3] - 0.5 * (nodes[0] + nodes[1])),
                         Norm(nodes[4] - 0.5 * (nodes[1] + nodes[2])),
                         Norm(nodes[5] - 0.5 * (nodes[2] + nodes[0]))});
      };

      auto append_provider_source_cell =
          [](std::vector<ProviderSourceOfTruthQuadraticCell>& cells,
             const std::array<Tddd, 6>& nodes,
             int reference_face_index,
             int body_object_index,
             int body_face_index,
             int source_type,
             int bc_type) {
            for (const auto& X : nodes)
              if (!isFinite(X))
                return;
            cells.push_back({nodes, reference_face_index, body_object_index, body_face_index,
                             source_type, bc_type});
          };

      auto write_provider_source_of_truth_quadratic_vtu =
          [&](const std::filesystem::path& path,
              const std::vector<ProviderSourceOfTruthQuadraticCell>& cells) {
            std::ofstream ofs(path);
            if (!ofs)
              throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, path.string() + " can not be opened");

            constexpr std::array<std::array<int, 3>, 4> subfaces = {{{0, 3, 5},
                                                                      {3, 1, 4},
                                                                      {5, 4, 2},
                                                                      {3, 4, 5}}};
            ofs << std::setprecision(17);
            ofs << "<?xml version=\"1.0\"?>\n";
            ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            ofs << "  <UnstructuredGrid>\n";
            ofs << "    <Piece NumberOfPoints=\"" << cells.size() * 6
                << "\" NumberOfCells=\"" << cells.size() * subfaces.size() << "\">\n";
            ofs << "      <Points>\n";
            ofs << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (const auto& cell : cells) {
              for (const auto& X : cell.nodes) {
                for (int i = 0; i < 3; ++i) {
                  const double value = std::abs(X[i]) < 1E-14 ? 0.0 : X[i];
                  ofs << value << ' ';
                }
              }
            }
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "      </Points>\n";
            ofs << "      <CellData Scalars=\"cell_scalars\">\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"reference_face_index\" format=\"ascii\">\n";
            for (const auto& cell : cells)
              for (std::size_t i = 0; i < subfaces.size(); ++i)
                ofs << cell.reference_face_index << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"body_object_index\" format=\"ascii\">\n";
            for (const auto& cell : cells)
              for (std::size_t i = 0; i < subfaces.size(); ++i)
                ofs << cell.body_object_index << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"body_face_index\" format=\"ascii\">\n";
            for (const auto& cell : cells)
              for (std::size_t i = 0; i < subfaces.size(); ++i)
                ofs << cell.body_face_index << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"source_type\" format=\"ascii\">\n";
            for (const auto& cell : cells)
              for (std::size_t i = 0; i < subfaces.size(); ++i)
                ofs << cell.source_type << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"bc_type\" format=\"ascii\">\n";
            for (const auto& cell : cells)
              for (std::size_t i = 0; i < subfaces.size(); ++i)
                ofs << cell.bc_type << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"subface_index\" format=\"ascii\">\n";
            for (std::size_t c = 0; c < cells.size(); ++c)
              for (std::size_t i = 0; i < subfaces.size(); ++i)
                ofs << i << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Float64\" Name=\"midpoint_deviation\" format=\"ascii\">\n";
            for (const auto& cell : cells)
              for (std::size_t i = 0; i < subfaces.size(); ++i)
                ofs << provider_cell_midpoint_deviation(cell.nodes) << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "      </CellData>\n";
            ofs << "      <Cells>\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
            for (std::size_t i = 0; i < cells.size(); ++i) {
              const std::size_t base = i * 6;
              for (const auto& sf : subfaces)
                ofs << base + sf[0] << ' ' << base + sf[1] << ' ' << base + sf[2] << ' ';
            }
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
            int offset = 0;
            for (std::size_t i = 0; i < cells.size() * subfaces.size(); ++i) {
              offset += 3;
              ofs << offset << ' ';
            }
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
            for (std::size_t i = 0; i < cells.size() * subfaces.size(); ++i)
              ofs << 5 << ' ';
            ofs << "\n";
            ofs << "        </DataArray>\n";
            ofs << "      </Cells>\n";
            ofs << "    </Piece>\n";
            ofs << "  </UnstructuredGrid>\n";
            ofs << "</VTKFile>\n";
          };

      auto dump_provider_source_of_truth_quadratic =
          [&](const BEMMeshPipeline::ReferenceState& reference, int step_retry) {
            std::vector<ProviderSourceOfTruthQuadraticCell> cells;
            std::string source_tag;

            if (const auto* snapshot = std::get_if<BEMMeshPipeline::SnapshotReferenceState>(&reference)) {
              source_tag = "snapshot";
              cells.reserve(snapshot->snap.faces.size());
              for (std::size_t i = 0; i < snapshot->snap.faces.size(); ++i) {
                const auto& face = snapshot->snap.faces[i];
                append_provider_source_cell(cells, BEMMeshPipeline::snapshotFaceNodes(face),
                                            static_cast<int>(i), -1, -1, 1,
                                            static_cast<int>(face.bc_type));
              }
            } else if (const auto* ic = std::get_if<BEMMeshPipeline::InitialConditionReferenceState>(&reference)) {
              source_tag = "ic";
              std::size_t reserve_count = 0;
              for (const auto* water : FluidObject)
                if (water)
                  reserve_count += water->getBoundaryFaces().size();
              cells.reserve(reserve_count);

              for (auto* water : FluidObject) {
                if (!water || !water->ic_eta)
                  continue;
                auto project_ic = [&](Tddd X) {
                  X[2] = water->ic_eta(X, ic->simulation_time);
                  return X;
                };
                for (auto* f : water->getBoundaryFaces()) {
                  if (!BEMMeshPipeline::faceHasDirichletSource(f))
                    continue;
                  auto [p0, p1, p2] = f->getPoints();
                  auto [l0, l1, l2] = f->getLines();
                  const std::array<Tddd, 6> nodes = {
                      project_ic(p0->X), project_ic(p1->X), project_ic(p2->X),
                      project_ic(BEMMeshPipeline::snapshotMidPosition(l0)),
                      project_ic(BEMMeshPipeline::snapshotMidPosition(l1)),
                      project_ic(BEMMeshPipeline::snapshotMidPosition(l2))};
                  append_provider_source_cell(cells, nodes, f ? f->index : -1, -1, -1, 2,
                                              static_cast<int>(NodeFaceBoundaryType::Dirichlet));
                }
              }
            }

            int body_object_index = 0;
            auto append_body_target_surface = [&](Network* body) {
              if (!body)
                return;
              for (auto* f : body->getBoundaryFaces()) {
                if (!f)
                  continue;
                auto [p0, p1, p2] = f->getPoints();
                auto [l0, l1, l2] = f->getLines();
                if (!p0 || !p1 || !p2 || !l0 || !l1 || !l2)
                  continue;
                const std::array<Tddd, 6> nodes = {
                    p0->X, p1->X, p2->X,
                    BEMMeshPipeline::snapshotMidPosition(l0),
                    BEMMeshPipeline::snapshotMidPosition(l1),
                    BEMMeshPipeline::snapshotMidPosition(l2)};
                append_provider_source_cell(cells, nodes, -1, body_object_index,
                                            f->index, 3,
                                            static_cast<int>(NodeFaceBoundaryType::Neumann));
              }
              ++body_object_index;
            };
            for (auto* body : RigidBodyObject)
              append_body_target_surface(body);
            for (auto* body : SoftBodyObject)
              append_body_target_surface(body);

            if (cells.empty()) {
              std::cout << Yellow << "[mesh:provider] skip source-of-truth quadratic dump: no cells"
                        << colorReset << std::endl;
              return;
            }

            const std::filesystem::path filename =
                "provider_source_of_truth_quadratic_" + source_tag +
                "_t" + std::to_string(time_step) +
                "_r" + std::to_string(step_retry) + ".vtu";
            write_provider_source_of_truth_quadratic_vtu(output_directory / filename, cells);
            provider_source_of_truth_pvd.push(filename.string(), simulation_time);
            provider_source_of_truth_pvd.output();
          };

      auto refresh_pipeline_boundary_state = [&](const char* reason) {
        for (auto* net : AllObjects)
          net->setGeometricPropertiesForce();
        _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);
        refreshBoundaryStatesAndTypes(FluidObject, Join(RigidBodyObject, SoftBodyObject));
        for (auto* water : FluidObject) {
          computeAllBCInterfaceMidpointOffsets(water);
        }
        if (mesh_pipeline.verbose_debug) {
          int total_interface = 0;
          int no_contact = 0;
          for (auto* water : FluidObject) {
            if (!water)
              continue;
            for (auto* l : water->getBoundaryLines()) {
              if (!l || !l->BCInterface)
                continue;
              ++total_interface;
              if (getEffectiveContactFaces(l).empty())
                ++no_contact;
            }
          }
          std::cout << Green << "[mesh:boundary] refresh reason=" << reason
                    << " bci_no_contact=" << no_contact << "/" << total_interface
                    << colorReset << std::endl;
        }
      };

      struct Phase1RawContactStats {
        int point_total = 0;
        int point_no_contact = 0;
        int line_total = 0;
        int line_no_contact = 0;

        bool has_missing() const {
          return point_no_contact > 0 || line_no_contact > 0;
        }
      };

	      auto inspect_phase1_raw_contact = [&](int retry_index, bool write_csv) {
	        Phase1RawContactStats stats;
        std::ostringstream rows;

        auto adjacent_faces_text = [](const std::vector<networkFace*>& faces) {
          std::ostringstream oss;
          bool first = true;
          for (const auto* f : faces) {
            if (!f)
              continue;
            if (!first)
              oss << ';';
            first = false;
            oss << f->index;
          }
          return oss.str();
        };

        struct NearestBodyHit {
          std::string body_name;
          int face_index = -1;
          double distance = std::numeric_limits<double>::infinity();
        };

        auto nearest_body_hit = [&](const Tddd& X) {
          NearestBodyHit best;
          for (auto* body : Join(RigidBodyObject, SoftBodyObject)) {
            if (!body)
              continue;
            auto [f, near_x] = body->Nearest(X);
            if (!f)
              continue;
            const double distance = Norm(X - near_x);
            if (std::isfinite(distance) && distance < best.distance) {
              best.body_name = body->getName();
              best.face_index = f->index;
              best.distance = distance;
            }
          }
          return best;
        };

        auto append_common = [&](const char* kind, int index, const Tddd& X, double contact_range,
                                 int raw_contact_count, int effective_contact_count,
                                 Network* penetrated_body, const std::string& adjacent_faces,
                                 const std::string& endpoints, bool dirichlet,
                                 bool neumann, bool bcinterface, bool multiple) {
          const auto nearest = nearest_body_hit(X);
          rows << time_step << ',' << retry_index << ',' << kind << ',' << index << ','
               << X[0] << ',' << X[1] << ',' << X[2] << ','
               << contact_range << ',' << raw_contact_count << ','
               << effective_contact_count << ','
               << (penetrated_body ? penetrated_body->getName() : "") << ','
               << nearest.body_name << ',' << nearest.face_index << ','
               << (std::isfinite(nearest.distance) ? nearest.distance : -1.0) << ','
               << adjacent_faces << ',' << endpoints << ','
               << (dirichlet ? 1 : 0) << ',' << (neumann ? 1 : 0) << ','
               << (bcinterface ? 1 : 0) << ',' << (multiple ? 1 : 0) << '\n';
        };

        for (auto* water : FluidObject) {
          if (!water)
            continue;
          for (auto* p : water->getBoundaryPoints()) {
            if (!p || !p->BCInterface)
              continue;
            ++stats.point_total;
            const int raw_count = static_cast<int>(p->getContactFaces().size());
            if (raw_count > 0)
              continue;
            ++stats.point_no_contact;
            append_common("point", p->index, p->getPosition(), p->contact_range, raw_count,
                          static_cast<int>(getEffectiveContactFaces(p).size()), p->penetratedBody,
                          adjacent_faces_text(p->getBoundaryFaces()), "",
                          p->Dirichlet, p->Neumann, p->BCInterface, p->isMultipleNode);
          }

          for (auto* l : water->getBoundaryLines()) {
            if (!l || !l->BCInterface)
              continue;
            ++stats.line_total;
            const int raw_count = static_cast<int>(l->getContactFaces().size());
            if (raw_count > 0)
              continue;
            ++stats.line_no_contact;
            auto [pA, pB] = l->getPoints();
            std::ostringstream endpoints;
            endpoints << (pA ? pA->index : -1) << ';' << (pB ? pB->index : -1);
            append_common("line", l->midpoint_index, l->getPosition(), l->contact_range, raw_count,
                          static_cast<int>(getEffectiveContactFaces(l).size()), l->penetratedBody,
                          adjacent_faces_text(l->getBoundaryFaces()), endpoints.str(),
                          l->Dirichlet, l->Neumann, l->BCInterface, l->isMultipleNode);
          }
        }

        if (write_csv && stats.has_missing()) {
          const auto path = output_directory /
                            ("phase1_raw_contact_missing_t" + std::to_string(time_step) +
                             "_retry" + std::to_string(retry_index) + ".csv");
          std::ofstream ofs(path);
          ofs << "time_step,retry,kind,index,x,y,z,contact_range,raw_contact_count,"
                 "effective_contact_count,penetrated_body,nearest_body,nearest_face,"
                 "nearest_distance,adjacent_faces,endpoints,dirichlet,neumann,bcinterface,multiple\n";
          ofs << rows.str();
          std::cout << Yellow << "[mesh:waterline] wrote raw contact diagnostic " << path
                    << colorReset << std::endl;
        }

	        return stats;
	      };

	      auto write_phase1_waterline_quality_witness =
	          [&](int retry_index, const BEMPreBVP::WaterlineMidpointStats& q) {
	            const auto path = output_directory /
	                              ("phase1_waterline_quality_witness_t" + std::to_string(time_step) +
	                               "_retry" + std::to_string(retry_index) + ".csv");
	            std::ofstream ofs(path);
	            ofs << "time_step,retry,kind,face_index,subtri_index,value,limit,subtri_count,"
	                   "min_angle_deg,max_aspect\n";
	            ofs << time_step << ',' << retry_index << ",min_angle,"
	                << q.min_angle_witness.face_index << ','
	                << q.min_angle_witness.subtri_index << ','
	                << q.min_angle_witness.value << ','
	                << setting.remeshing.waterline_midpoint_min_angle_deg << ','
	                << q.waterline_subtri_count << ','
	                << q.subtri_min_angle_deg << ','
	                << q.subtri_max_aspect << '\n';
	            ofs << time_step << ',' << retry_index << ",max_aspect,"
	                << q.max_aspect_witness.face_index << ','
	                << q.max_aspect_witness.subtri_index << ','
	                << q.max_aspect_witness.value << ','
	                << setting.remeshing.waterline_midpoint_max_aspect << ','
	                << q.waterline_subtri_count << ','
	                << q.subtri_min_angle_deg << ','
	                << q.subtri_max_aspect << '\n';
	            std::cout << Yellow << "[mesh:waterline] wrote quality witness " << path
	                      << colorReset << std::endl;
	          };

	      auto write_waterline_clamped_entities =
	          [&](const std::string& tag, int retry_index,
	              const BEMPreBVP::WaterlineRefreshSummary& refresh) {
	            if (refresh.clamped_entities.empty())
	              return;
	            const auto path = output_directory /
	                              ("clamped_entities_" + tag + "_t" +
	                               std::to_string(time_step) + "_retry" +
	                               std::to_string(retry_index) + ".csv");
	            std::ofstream ofs(path);
	            ofs << "time_step,retry,tag,line_index,endpoint0_index,endpoint1_index,"
	                   "reference_face_index,body_face_index,free_gap,body_gap,move_ratio,"
	                   "local_min_angle_before,local_min_angle_after,"
	                   "local_max_aspect_before,local_max_aspect_after,mid_x,mid_y,mid_z\n";
	            for (const auto& e : refresh.clamped_entities) {
	              ofs << time_step << ',' << retry_index << ',' << tag << ','
	                  << e.line_index << ','
	                  << e.endpoint0_index << ','
	                  << e.endpoint1_index << ','
	                  << e.reference_face_index << ','
	                  << e.body_face_index << ','
	                  << e.free_gap << ','
	                  << e.body_gap << ','
	                  << e.move_ratio << ','
	                  << e.local_min_angle_before << ','
	                  << e.local_min_angle_after << ','
	                  << e.local_max_aspect_before << ','
	                  << e.local_max_aspect_after << ','
	                  << e.midpoint[0] << ','
	                  << e.midpoint[1] << ','
	                  << e.midpoint[2] << '\n';
	            }
	            std::cout << Yellow << "[mesh:waterline] wrote clamped entities " << path
	                      << colorReset << std::endl;
	          };

	      auto enforce_pre_bvp_guard = [&](const BEMPreBVP::Stats& guard_stats, int rk_step) {
        if (guard_stats.has_numeric_failure()) {
          throw step_failure("pre-BVP consistency guard found nonfinite/overflow active DOF values at RK_step " +
                             std::to_string(rk_step));
        }
        if (guard_stats.bcinterface_points_no_contact > 0 ||
            guard_stats.bcinterface_lines_no_contact > 0) {
          std::ostringstream oss;
          oss << "pre-BVP consistency guard found BCInterface entity without contact"
              << " at RK_step " << rk_step
              << " (points=" << guard_stats.bcinterface_points_no_contact
              << ", lines=" << guard_stats.bcinterface_lines_no_contact << ")";
          throw step_failure(oss.str());
        }
        if (setting.remeshing.waterline_protection_enabled && mesh_pipeline.waterline_geometry) {
          const auto& q = guard_stats.waterline_midpoint_quality;
          if (q.waterline_subtri_count > 0) {
            const bool min_angle_bad =
                !std::isfinite(q.subtri_min_angle_deg) ||
                q.subtri_min_angle_deg < setting.remeshing.waterline_midpoint_min_angle_deg;
            const bool aspect_bad =
                !std::isfinite(q.subtri_max_aspect) ||
                q.subtri_max_aspect > setting.remeshing.waterline_midpoint_max_aspect;
	            if (min_angle_bad || aspect_bad) {
	              std::ostringstream oss;
	              oss << "pre-BVP waterline midpoint quality failed at RK_step " << rk_step
	                  << " (subtri_count=" << q.waterline_subtri_count
	                  << ", min_angle_deg=" << q.subtri_min_angle_deg
	                  << ", max_aspect=" << q.subtri_max_aspect
	                  << ", min_angle_limit=" << setting.remeshing.waterline_midpoint_min_angle_deg
	                  << ", max_aspect_limit=" << setting.remeshing.waterline_midpoint_max_aspect
	                  << "," << BEMPreBVP::waterline_midpoint_witness_summary(q) << ")";
	              throw step_failure(oss.str());
	            }
          }
          if (guard_stats.waterline_faces > 0) {
            const bool face_min_angle_bad =
                !std::isfinite(guard_stats.waterline_face_min_angle_deg) ||
                guard_stats.waterline_face_min_angle_deg < setting.remeshing.waterline_face_min_angle_deg;
            const bool face_aspect_bad =
                !std::isfinite(guard_stats.waterline_face_aspect_max) ||
                guard_stats.waterline_face_aspect_max > setting.remeshing.waterline_face_max_aspect;
            if (face_min_angle_bad || face_aspect_bad) {
              std::ostringstream oss;
              oss << "pre-BVP waterline face quality failed at RK_step " << rk_step
                  << " (faces=" << guard_stats.waterline_faces
                  << ", min_angle_deg=" << guard_stats.waterline_face_min_angle_deg
                  << ", max_aspect=" << guard_stats.waterline_face_aspect_max
                  << ", min_angle_limit=" << setting.remeshing.waterline_face_min_angle_deg
                  << ", max_aspect_limit=" << setting.remeshing.waterline_face_max_aspect
                  << ")";
              throw step_failure(oss.str());
            }
          }
	          const auto& refresh = guard_stats.waterline_refresh;
	          if (refresh.valid && refresh.clamped_count > 0 && refresh.max_move_ratio > 0.8) {
	            write_waterline_clamped_entities("pre_bvp_rk" + std::to_string(rk_step), -1, refresh);
	            std::cout << Yellow
	                      << "[mesh:waterline] pre-BVP repair clamped but contact/quality remained acceptable"
	                      << " rk_step=" << rk_step
	                      << " clamped=" << refresh.clamped_count
	                      << " max_move_ratio=" << refresh.max_move_ratio
	                      << " midpoint_min_angle=" << q.subtri_min_angle_deg
	                      << " midpoint_max_aspect=" << q.subtri_max_aspect
	                      << colorReset << std::endl;
	          }
	        }
	      };

      bool mesh_prepared_for_rk = false;
      save_step_state();
      for (int step_retry = 0; step_retry <= StepRetryState::max_step_retries + (retry_state.degraded_mode ? 1 : 0); ++step_retry) {
        if (step_retry > 0)
          restore_step_state(step_retry);

        try {
          //@ --- 3c-2a. Remesh + quality checks ---
          std::optional<BEMMeshPipeline::ReferenceState> mesh_reference;
          bool initial_conditions_applied_by_pipeline = false;
          const bool run_mesh_preparation_this_retry = !use_mesh_preparation_pipeline || !mesh_prepared_for_rk;

          if (run_mesh_preparation_this_retry) {
            const bool is_initial_step = (time_step == 0 && start_time_step == 0);
            const auto contact_objects = Join(RigidBodyObject, SoftBodyObject);
            const auto geometry_projection_mode = use_true_quadratic_element
                                                      ? BEMMeshPipeline::GeometryProjectionMode::TrueQuadraticFace
                                                      : BEMMeshPipeline::GeometryProjectionMode::ReferenceQuadraticFace;

            if (use_mesh_preparation_pipeline) {
              _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);
              refreshBoundaryStatesAndTypes(FluidObject, contact_objects);

              if (is_initial_step) {
                applyInitialConditions(FluidObject, RigidBodyObject, SoftBodyObject, AllObjects, use_true_quadratic_element, node_relocation_surface);
                initial_conditions_applied_by_pipeline = true;
                if (initial_mesh_pre_relax && initial_mesh_pre_relax_loop > 0 && initial_mesh_pre_relax_coef > 0.0) {
                  const bool has_ic = std::ranges::any_of(FluidObject, [](const Network* w) { return w->ic_eta && w->ic_phi; });
                  if (has_ic) {
                    preRelaxMesh(FluidObject, AllObjects, contact_objects,
                                 use_true_quadratic_element, node_relocation_surface, simulation_time,
                                 {.loop = initial_mesh_pre_relax_loop, .coef = initial_mesh_pre_relax_coef, .output_tag = ""});
                    applyInitialConditions(FluidObject, RigidBodyObject, SoftBodyObject, AllObjects, use_true_quadratic_element, node_relocation_surface);
                  }
                }
              }

              refresh_pipeline_boundary_state("phase0-reference");
              mesh_reference.emplace(BEMMeshPipeline::makeReferenceState(FluidObject, time_step, start_time_step, simulation_time));
              dump_provider_source_of_truth_quadratic(*mesh_reference, step_retry);
              write_all_water_phase_snapshots("after_phase0", step_retry, water_after_phase0_pvd);
              dump_pipeline_mesh("phase1_before_geometry_repair", mesh_pipeline.debug_dump_all_phases);
              if (mesh_pipeline.dump_after_phase0 || mesh_pipeline.debug_dump_all_phases)
                dump_pipeline_mesh(BEMMeshPipeline::isInitialConditionReference(*mesh_reference) ? "phase0_ic_source" : "phase0_snapshot_source", true);

              const bool skip_initial_geometry_repair = is_initial_step;
              if (mesh_pipeline.phase1_clung && skip_initial_geometry_repair) {
                std::cout << Green << "[mesh:geometry] skip phase=phase1 reason=t0_initial_reference" << colorReset << std::endl;
              } else if (mesh_pipeline.phase1_clung) {
                BEMMeshPipeline::AdjustStats adjust_stats;
                BEMMeshPipeline::RemeshLog::GeometryRepairSummary phase1_summary;
                phase1_summary.phase = "phase1";
                if (mesh_pipeline.geometry_projector_v2) {
                  const auto repair_report = BEMMeshPipeline::GeometryRepair::repairCurrentSurfaceGeometryWithProjector(
                      FluidObject, *mesh_reference, contact_objects, mesh_pipeline, geometry_projection_mode);
                  adjust_stats = repair_report.stats;
                  phase1_summary.repaired_points = repair_report.repaired_points;
	                  phase1_summary.repaired_lines = repair_report.repaired_lines;
	                  phase1_summary.quality_damaged = repair_report.quality_damaged_lines;
	                  phase1_summary.damaged_faces = repair_report.quality_damaged_faces;
	                  phase1_summary.penetration_rejected = repair_report.penetration_rejected_before_repair;
                  BEMMeshPipeline::GeometryRepair::finalizeProjectedGeometryRepair(
                      FluidObject, AllObjects, contact_objects, time_step,
                      "post-mesh-pipeline-waterline-phase1", use_true_quadratic_element);
                } else {
                  adjust_stats = BEMMeshPipeline::GeometryRepair::repairCurrentSurfaceGeometry(
                      FluidObject,
                      mesh_pipeline.clung_max_iter,
                      mesh_pipeline.clung_tol,
                      mesh_pipeline.clung_move_limit_factor,
                      node_relocation_debug_output ? &output_directory : nullptr);
                }
                phase1_summary.iterations = adjust_stats.iterations_done;
                phase1_summary.max_move = adjust_stats.max_move_final;
                phase1_summary.nodes_adjusted = adjust_stats.nodes_adjusted;
                phase1_summary.lines_snapped = adjust_stats.lines_snapped;
                BEMMeshPipeline::RemeshLog::printPipelineGeometryRepair(phase1_summary);
                refresh_pipeline_boundary_state("after-phase1-geometry-repair");
                dump_pipeline_mesh("phase1_after_geometry_repair", mesh_pipeline.dump_after_phase1 || mesh_pipeline.debug_dump_all_phases);
              }

              if (!skip_initial_geometry_repair) {
                if (!mesh_pipeline.geometry_projector_v2) {
                  BEMMeshPipeline::GeometryRepair::repairWaterlineMidpointsIfEnabled(
                      FluidObject, AllObjects, contact_objects, *mesh_reference, mesh_pipeline,
                      geometry_projection_mode, time_step, "post-mesh-pipeline-waterline-phase1",
                      use_true_quadratic_element);
                  refresh_pipeline_boundary_state("after-phase1-waterline-repair");
                }
	                if (setting.remeshing.waterline_protection_enabled) {
	                  const auto& refresh = BEMPreBVP::latest_waterline_refresh_summary();
                  const auto raw_contact_stats = inspect_phase1_raw_contact(step_retry, true);
                  BEMPreBVP::WaterlineMidpointStats phase1_waterline_midpoint_quality;
                  for (auto* water : FluidObject) {
                    if (!water)
                      continue;
		                    BEMPreBVP::merge_waterline_midpoint_quality(
		                        phase1_waterline_midpoint_quality,
		                        BEMPreBVP::measure_waterline_midpoint_quality(water->getBoundaryFaces()));
		                  }
	                  if (raw_contact_stats.has_missing()) {
	                    std::ostringstream oss;
	                    oss << "phase1 waterline repair left BCInterface entity(s) without raw contact before remesh at time_step "
		                        << time_step
		                        << " (points=" << raw_contact_stats.point_no_contact << "/" << raw_contact_stats.point_total
		                        << ", lines=" << raw_contact_stats.line_no_contact << "/" << raw_contact_stats.line_total << ")";
		                    throw step_failure(oss.str());
		                  }
	                  const bool phase1_midpoint_quality_bad =
	                      phase1_waterline_midpoint_quality.waterline_subtri_count > 0 &&
	                      ((std::isfinite(phase1_waterline_midpoint_quality.subtri_min_angle_deg) &&
	                        phase1_waterline_midpoint_quality.subtri_min_angle_deg <
	                            setting.remeshing.waterline_midpoint_min_angle_deg) ||
	                       (std::isfinite(phase1_waterline_midpoint_quality.subtri_max_aspect) &&
	                        phase1_waterline_midpoint_quality.subtri_max_aspect >
	                            setting.remeshing.waterline_midpoint_max_aspect));
		                  if (phase1_midpoint_quality_bad) {
		                    write_phase1_waterline_quality_witness(step_retry, phase1_waterline_midpoint_quality);
		                    std::ostringstream oss;
		                    oss << "phase1 waterline midpoint quality failed before remesh at time_step "
		                        << time_step
	                        << " (subtri_count=" << phase1_waterline_midpoint_quality.waterline_subtri_count
	                        << ", min_angle_deg="
	                        << phase1_waterline_midpoint_quality.subtri_min_angle_deg
	                        << ", max_aspect="
	                        << phase1_waterline_midpoint_quality.subtri_max_aspect
	                        << ", min_angle_limit="
	                        << setting.remeshing.waterline_midpoint_min_angle_deg
	                        << ", max_aspect_limit="
	                        << setting.remeshing.waterline_midpoint_max_aspect
	                        << ", clamped=" << (refresh.valid ? refresh.clamped_count : 0)
		                        << ", max_move_ratio=" << (refresh.valid ? refresh.max_move_ratio : 0.0)
		                        << ", quality_damaged="
		                        << (refresh.valid ? refresh.quality_damaged_after_repair : 0)
		                        << "," << BEMPreBVP::waterline_midpoint_witness_summary(phase1_waterline_midpoint_quality)
		                        << ")";
		                    if (!mesh_pipeline.phase2_remesh) {
		                      throw step_failure(oss.str());
		                    }
		                    std::cout << Yellow << "[mesh:waterline] " << oss.str()
		                              << " action=defer_to_phase2_remesh"
		                              << colorReset << std::endl;
		                  }
                  if (refresh.valid && refresh.clamped_count > 0 && refresh.max_move_ratio > 0.8) {
                    write_waterline_clamped_entities("phase1", step_retry, refresh);
                    if (refresh.quality_damaged_after_repair > 0) {
                      std::ostringstream oss;
                      oss << "phase1 waterline repair clamped and degraded quality before remesh at time_step "
                          << time_step
                          << " (clamped=" << refresh.clamped_count
                          << ", max_move_ratio=" << refresh.max_move_ratio
                          << ", quality_damaged=" << refresh.quality_damaged_after_repair
                          << ", midpoint_min_angle="
                          << phase1_waterline_midpoint_quality.subtri_min_angle_deg
                          << ", midpoint_max_aspect="
                          << phase1_waterline_midpoint_quality.subtri_max_aspect << ")";
                      throw step_failure(oss.str());
                    }
	                    std::cout << Yellow
	                              << "[mesh:waterline] phase1 repair clamped but contact/quality remained acceptable"
	                              << " time_step=" << time_step
	                              << " clamped=" << refresh.clamped_count
	                              << " max_move_ratio=" << refresh.max_move_ratio
	                              << " midpoint_min_angle="
	                              << phase1_waterline_midpoint_quality.subtri_min_angle_deg
	                              << " midpoint_max_aspect="
	                              << phase1_waterline_midpoint_quality.subtri_max_aspect
	                              << colorReset << std::endl;
	                  }
                }
                dump_pipeline_mesh("phase1_after_waterline_repair", mesh_pipeline.dump_after_phase1 || mesh_pipeline.debug_dump_all_phases);
              }
              write_all_water_phase_snapshots("after_phase1", step_retry, water_after_phase1_pvd);
            }

            if (use_mesh_preparation_pipeline && mesh_pipeline.phase2_remesh && is_initial_step) {
              std::cout << Green << "[mesh:topology] skip reason=t0_initial_reference" << colorReset << std::endl;
              write_all_water_phase_snapshots("after_phase2", step_retry, water_after_phase2_pvd);
            } else if (!use_mesh_preparation_pipeline || mesh_pipeline.phase2_remesh) {
              _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);
              refreshBoundaryStatesAndTypes(FluidObject, contact_objects);
              auto remesh_settings = setting.remeshing;
              remesh_settings.defer_scalar_interpolation = use_mesh_preparation_pipeline && mesh_pipeline.phase3_interpolate;
              // Trial remesh が BoundaryGeometryTargetProvider を作るための材料。
              // reference は Dirichlet/free-surface 側の参照形状、
              // topology_repair_body_objects は Neumann/body 側の現在形状。
              const auto topology_repair_body_objects = contact_objects;
              const BEMMeshPipeline::RemeshTrialGeometryInputs trial_geometry_inputs{
                  mesh_reference ? &*mesh_reference : nullptr,
                  &topology_repair_body_objects,
                  &mesh_pipeline,
                  geometry_projection_mode,
                  static_cast<std::uint64_t>(std::max(time_step, 0) + 1)};
              const auto* trial_geometry_inputs_ptr =
                  (use_mesh_preparation_pipeline && mesh_reference && mesh_pipeline.geometry_projector_v2)
                      ? &trial_geometry_inputs
                      : nullptr;
              if (use_mesh_preparation_pipeline && mesh_pipeline.phase2_remesh)
                BEMMeshPipeline::RemeshLog::printPipelineTopologyRun("phase2");
              for (auto& water : FluidObject) {
                BEMMeshPipeline::TopologyModifier::applyLocalRemesh(
                    *water, time_step, remesh_settings,
                    retry_state.degraded_mode,
                    output_directory.string(), simulation_time,
                    &candidate_patches_pvd, &remeshed_patches_pvd,
                    &trigger_edges_pvd, &edges_pvd,
                    trial_geometry_inputs_ptr,
                    use_mesh_preparation_pipeline ? "phase2" : "");
                retry_state.collapse_repeatedly_rejected_faces(*water, step_retry);
                write_water_phase_snapshot(water, step_retry, "after_phase2", water_after_phase2_pvd);
                if (!retry_state.degraded_mode) {
                  refreshFaceBadQualityHistory(*water, time_step, std::nullopt, 0.1, subsurface_altitude_reject);
                  monitorTinyFaceRelativeToLocalMean(*water, time_step, std::nullopt);
                  throwIfSubsurfaceFaceAltitudeTooSmall(*water, time_step, std::nullopt, subsurface_altitude_reject);
                  throwIfTinyFaceRelativeToLocalMean(*water, time_step, std::nullopt, 0.01);
                }
              }
            }

            refreshBoundaryStatesAndTypes(FluidObject, contact_objects);
            BEM_Penetration::throwIfStructurePenetrated(FluidObject, contact_objects, time_step, use_mesh_preparation_pipeline ? "post-mesh-preparation" : "post-remesh", use_true_quadratic_element);

            if (use_mesh_preparation_pipeline && mesh_reference && mesh_pipeline.phase3_interpolate) {
              BEMMeshPipeline::FieldMapper::mapFieldsFromReferenceState(
                  FluidObject, *mesh_reference, node_relocation_surface,
                  mesh_pipeline.pseudo_quadratic_policy);
              if (mesh_pipeline.dump_after_phase3)
                dump_pipeline_mesh("phase3_interpolate", true);
            }

            if (use_mesh_preparation_pipeline) {
              BEMPreBVP::Options guard_options;
              guard_options.print_ok_summary = true;
              guard_options.print_waterline_midpoint_quality = mesh_pipeline.waterline_geometry;
              const auto guard_stats = BEMPreBVP::inspect(FluidObject, time_step, 0, guard_options);
              enforce_pre_bvp_guard(guard_stats, 0);
            }

            save_post_remesh_state();
            mesh_prepared_for_rk = use_mesh_preparation_pipeline;
          }

          //@ --- 3c-2b. Initial conditions (t=0 only) ---
          if (!use_mesh_preparation_pipeline && time_step == 0 && start_time_step == 0 && !initial_conditions_applied_by_pipeline) {
            applyInitialConditions(FluidObject, RigidBodyObject, SoftBodyObject, AllObjects, use_true_quadratic_element, node_relocation_surface);

            // Post-IC mesh quality recovery: IC projects X_mid to wave surface,
            // which can create tiny/collapsed subfaces with true_quadratic elements.
            // Re-run pre-relax (with flip + smoothing) to recover mesh quality,
            // then re-apply IC to correct any eta/phi drift from smoothing.
            if (initial_mesh_pre_relax && initial_mesh_pre_relax_loop > 0 && initial_mesh_pre_relax_coef > 0.0) {
              bool has_ic = std::ranges::any_of(FluidObject, [](const Network* w) { return w->ic_eta && w->ic_phi; });
              if (has_ic) {
                const auto contact_objects = Join(RigidBodyObject, SoftBodyObject);
                preRelaxMesh(FluidObject, AllObjects, contact_objects,
                             use_true_quadratic_element, node_relocation_surface, simulation_time,
                             {.loop = initial_mesh_pre_relax_loop, .coef = initial_mesh_pre_relax_coef, .output_tag = ""});
                applyInitialConditions(FluidObject, RigidBodyObject, SoftBodyObject, AllObjects, use_true_quadratic_element, node_relocation_surface);
              }
            }
          }

          CoordinateBounds bounds_org(FluidObject[0]->bounds);
          for (auto& water : FluidObject)
            bounds_org += water->bounds;

          CoordinateBounds bounds = bounds_org.scaledBounds(1.5);

          Buckets_Fluid_Faces.initialize(bounds, bounds.getScale() / 10.);
          Buckets_Fluid_Points.initialize(bounds, bounds.getScale() / 10.);

          //& --- 3c-2c. dt determination ---
          double dt_motion_min = max_dt, dt_geom_min = max_dt;
          for (auto& water : FluidObject) {
            show_info(*water);

            const double dt_motion = dt_mesh_motion(*water, max_dt, 0.5);
            const double dt_geom = dt_geometry(*water, max_dt, 0.2);
            dt_motion_min = std::min(dt_motion_min, dt_motion);
            dt_geom_min = std::min(dt_geom_min, dt_geom);
            const double dt_constrained = std::min(dt_motion, dt_geom);
            if (dt > dt_constrained)
              dt = dt_constrained;
            if (time_step <= 2)
              dt = dt / 10.;
          }
          std::cout << "[dt] dt_motion=" << dt_motion_min << " dt_geom=" << dt_geom_min;
          if (dt_geom_min < dt_motion_min)
            std::cout << Yellow << " (limited by geometry)" << colorReset;
          std::cout << std::endl;

          vpm.refreshKinematicsForDt([&](const std::array<double, 3>& x) { return getBEMVelocityAt_cached(x, FluidObject); });
          vpm.constrainDt(dt);
          if (dt < 1E-13)
            dt = 1E-13;
          const double dt_growth_limit = 1.25;
          if (last_successful_dt > 0.0) {
            const double dt_capped = last_successful_dt * dt_growth_limit;
            if (dt > dt_capped) {
              std::cout << Yellow << "[dt] growth limited: " << dt << " -> " << dt_capped << " (last_successful_dt=" << last_successful_dt << ")" << colorReset << std::endl;
              dt = dt_capped;
            }
          }
          // Apply dt override from step rejection (smaller dt for retry)
          if (retry_state.dt_override > 0 && retry_state.dt_override < dt) {
            std::cout << Yellow << "[step_reject] dt reduced: " << dt << " -> " << retry_state.dt_override << colorReset << std::endl;
            dt = retry_state.dt_override;
          }
          Print("===========================================================================");
          Print("       dt :", Red, std::setprecision(10), dt, colorReset);
          Print("time_step :", Red, time_step, colorReset);
          Print("real time :", Red, simulation_time, colorReset);
          Print("---------------------------------------------------------------------------");

          //^ --- VPM Strang splitting: diffusion half-step (beginning) ---
          vpm.applyStrangDiffusionHalf(dt, RigidBodyObject);

          // b@ -------------------------------------------------------------- */

          auto rebuild_fluid_global_bounds = [&]() {
            CoordinateBounds bounds_org(FluidObject[0]->bounds);
            for (auto& water : FluidObject)
              bounds_org += water->bounds;
            return CoordinateBounds(bounds_org.scaledBounds(1.5));
          };

          auto rebuild_Buckets_Fluid_Points = [&]() {
            const auto bounds = rebuild_fluid_global_bounds();
            Buckets_Fluid_Points.initialize(bounds, bounds.getScale() / 10.);
            for (auto& water : FluidObject)
              for (const auto& p : water->getPoints())
                Buckets_Fluid_Points.add(ToX(p), p);
          };

          auto rebuild_Buckets_Fluid_Faces = [&]() {
            const auto bounds = rebuild_fluid_global_bounds();
            Buckets_Fluid_Faces.initialize(bounds, bounds.getScale() / 10.);
            for (auto& water : FluidObject)
              for (const auto& f : water->getBoundaryFaces())
                if (f)
                  Buckets_Fluid_Faces.add(f->getXtuple(), f);
          };

          auto rebuild_fluid_nodes = [&]() {
            fluid_nodes.clear();
            for (auto* water : FluidObject) {
              for (auto* p : water->getPoints())
                fluid_nodes.push_back(p);
              for (auto* l : water->getBoundaryLines())
                if (l->hasActiveBieDof())
                  fluid_nodes.push_back(l);
            }
          };

          // b@ -------------------------------------------------------------- */

          auto RK_step = 0;
          TimeWatch watch;

          //& --- 3c-2e. RK loop (do-while) ---
          do {
            retry_state.current_rk_step = ++RK_step;
            //
            if (RK_step == 1) {
              //% --- setBoundaryTypes (RK_step==1 only) ---
              _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);
              if (bem_verbose())
                PrintLap(watch, "makeBuckets");
              refreshBoundaryStatesAndTypes(FluidObject, Join(RigidBodyObject, SoftBodyObject));
              for (auto& water : FluidObject) {
                water->setMinDepthFromBCInterface();
                computeAllBCInterfaceMidpointOffsets(water);
                const auto collapse_result = collapseCornerConnectedNeumannLinesAfterBoundaryTypes(*water, time_step);
                if (collapse_result.changed) {
                  refreshBoundaryStatesAndTypes(water, Join(RigidBodyObject, SoftBodyObject));
                  water->setMinDepthFromBCInterface();
                  computeAllBCInterfaceMidpointOffsets(water);
                  std::cout << Magenta << "[corner_neumann_post_type] recategorized after " << collapse_result.collapsed << " collapse(s)" << colorReset << std::endl;
                } else
                  logCornerConnectedNeumannLinesAfterBoundaryTypes(water, "post-setBoundaryTypes");
                dumpDebugCornerPointState(water, "post-setBoundaryTypes", time_step, RK_step);
              }
              BVP.matrix_size = setNodeFaceIndices(FluidObject);
              rebuild_Buckets_Fluid_Points();
              rebuild_Buckets_Fluid_Faces();
              rebuild_fluid_nodes();
              if (bem_verbose())
                PrintLap(watch, "setBoundaryTypes");

              /* -------------------------------------------------------------------------- */

              //& --- 3c-2d. RK initialization ---
              auto RK_init = [&](auto entity) {
                if (entity->hasActiveBieDof()) {
                  entity->RK_phi.initialize(dt, simulation_time, std::get<0>(entity->phiphin), RK_order);
                  entity->RK_X.initialize(dt, simulation_time, entity->getPosition(), RK_order);
                }
              };
              for (auto& water : FluidObject) {
                for (auto* p : water->getPoints())
                  RK_init(p);
                for (auto* l : water->getBoundaryLines())
                  RK_init(l);
              }

              for (const auto& net : RigidBodyObject) {
                net->RK_COM.initialize(dt, simulation_time, net->COM, RK_order);
                net->RK_Q.initialize(dt, simulation_time, net->Q(), RK_order);
                net->RK_Velocity.initialize(dt, simulation_time, net->velocity, RK_order);
                //
                if (net->interp_accel.size() > 10)
                  net->interp_accel.pop();
                net->interp_accel.push(simulation_time, net->acceleration);
              }
              for (const auto& net : SoftBodyObject) {
                // SoftBody でも（重心・姿勢・速度）を RK で進めるコードパスがあるため初期化が必要
                net->RK_COM.initialize(dt, simulation_time, net->COM, RK_order);
                net->RK_Q.initialize(dt, simulation_time, net->Q(), RK_order);
                net->RK_Velocity.initialize(dt, simulation_time, net->velocity, RK_order);
                for (const auto& p : net->getPoints())
                  p->RK_X.initialize(dt, simulation_time, ToX(p), RK_order);
              }

              // VPM particle update:
              // - diffusion is split (dt/2 at beginning/end of time-step)
              // - advection + stretching are advanced in the RK loop (no diffusion term)
              vpm.initializeRKIfEnabled(dt, simulation_time, RK_order);
            }

            auto RK_time = (*FluidObject[0]->getPoints().begin())->RK_X.gett(); //%各ルンゲクッタの時刻を使う
            std::cout << "RK_step = " << RK_step << "/" << RK_order << ", RK_time = " << RK_time << ", simulation_time = " << simulation_time << std::endl;

            //% --- Absorber assignment (before BVP) ---
            for (auto* node : fluid_nodes) {
              node->absorbedBy = nullptr;
              for (const auto& net : AbsorberObject)
                if (net->InsideQ(node->getPosition()))
                  node->absorbedBy = net;
            }

            //% --- BVP: setPhiPhinOnFace -> BVP.solve -> velocity ---
            setBodyVelocity(Join(RigidBodyObject, SoftBodyObject));
            if (bem_verbose())
              PrintLap(watch, "setBodyVelocity");

            setPhiPhinOnFace(FluidObject);
            for (auto water : FluidObject)
              dumpDebugCornerPointState(water, "post-setPhiPhinOnFace", time_step, RK_step);
            if (bem_verbose())
              PrintLap(watch, "setPhiPhinOnFace");

            //^ --- VPM Coupling: subtract vortex-induced velocity from Neumann BC ---
            vpm.subtractFromNeumannBC(FluidObject);

            {
              BEMPreBVP::Options guard_options;
              guard_options.print_ok_summary = use_mesh_preparation_pipeline && RK_step == 1;
              guard_options.print_waterline_midpoint_quality = mesh_pipeline.waterline_geometry;
              const auto guard_stats = BEMPreBVP::inspect(FluidObject, time_step, RK_step, guard_options);
              enforce_pre_bvp_guard(guard_stats, RK_step);
            }

            for (auto water : FluidObject)
              dumpDebugCornerPointState(water, "pre-BVP.solve", time_step, RK_step);
            auto [time_setup_, time_solve_, unknownsize_] = BVP.solve();
            if (solver_type == "GMRES" && !BVP.last_gmres_converged) {
              const double gmres_retry_threshold = std::max(1e-3, solver_tol * 1e6);
              std::cerr << Yellow << "[GMRES] time_step " << time_step << " RK_step " << RK_step << ": residual=" << BVP.last_gmres_residual_norm << " iter=" << BVP.last_gmres_total_iter << colorReset << std::endl;
              BVP.logLastGmresHotspots(std::cerr, 8);
              for (auto water : FluidObject)
                dumpDebugCornerPointState(water, "post-BVP.solve-gmres-failed", time_step, RK_step);
              if (BVP.last_gmres_residual_norm > gmres_retry_threshold) {
                throw step_failure("GMRES not converged at RK_step " + std::to_string(RK_step) + ", residual=" + std::to_string(BVP.last_gmres_residual_norm));
              }
              std::cerr << Yellow << "[GMRES] RK_step " << RK_step << ": residual " << BVP.last_gmres_residual_norm << " above solver_tol=" << solver_tol << " but below retry threshold=" << gmres_retry_threshold << ", continuing" << colorReset << std::endl;
            }
            time_setup = time_setup_;
            time_solve = time_solve_;
            ElapsedTimeSetup.push_back(time_setup);
            ElapsedTimeSolve.push_back(time_solve);
            unknownsizse = unknownsize_;
            if (bem_verbose())
              PrintLap(watch, "BVP.solve -> {Φ,Φn}が決まる");

            set_u_potential_BEM(FluidObject);
            if (bem_verbose())
              PrintLap(watch, "set_u_potential_BEM");

            set_u_total(FluidObject);
            if (bem_verbose())
              PrintLap(watch, "set_u_total");

            const bool do_node_relocation = (node_relocation_method != NodeRelocationMethod::none) && (RK_step == RK_order && time_step % node_relocation_period == 0);
            const auto* relocation_debug_dir = node_relocation_debug_output ? &output_directory : nullptr;
            setNodeVelocity(FluidObject, do_node_relocation ? 30 : 0, 0.05, relocation_debug_dir); //!
            if (bem_verbose())
              PrintLap(watch, "setNodeVelocity");

            //^ --- VPM Coupling: advection + stretching (no diffusion) ---
            vpm.advectAndStretchRK(RigidBodyObject, [&](const std::array<double, 3>& x) { return getBEMVelocityAt_cached(x, FluidObject); });

            //% --- Body update: pre-solve read (COM/Q sync, mooring) ---
            _Pragma("omp parallel for") for (const auto& net : Join(RigidBodyObject, SoftBodyObject)) {
              if (net->inputJSON.find("velocity")) {
                //! 時間まで静止させる
                if ((net->inputJSON.at("velocity")[0].contains("floating") && net->inputJSON.at("velocity").size() >= 2 && std::stod(net->inputJSON.at("velocity")[1]) > simulation_time)) {
                  net->velocity.fill(0);
                  net->acceleration.fill(0);
                }

                //! 静止していても，出力や力計算のためにCOM/Qは常に最新(RK)で同期しておく
                if (!net->inputJSON.at("velocity")[0].contains("fixed"))
                  std::cout << "updating " << net->getName() << "'s (Rigid/SoftBodyObject) velocity" << std::endl;
                net->COM = net->RK_COM.getX();
                net->Q = Normalize(net->RK_Q.getX());

                std::cout << "name = " << net->getName() << std::endl;
                std::cout << "net->velocityTranslational() = " << net->velocityTranslational() << std::endl;

                for (auto& mooring : net->mooringLines) {
                  //! mooring->lastPoint は，浮体ともに動く．
                  auto Xcurrent = mooring->lastPoint->X;
                  mooring->lastPoint->X_last = Xcurrent;
                  Tddd V = (nextPositionOnBody(net, mooring->lastPoint) - Xcurrent) / (net->RK_COM.getTimeAtNextStep() - simulation_time);
                  mooring->simulate(simulation_time, net->RK_COM.getTimeAtNextStep() - simulation_time, [&](networkPoint* p) {
                    if (p == mooring->firstPoint) {
                      p->acceleration.fill(0);
                      p->velocity.fill(0);
                    } else if (p == mooring->lastPoint) {
                      p->acceleration.fill(0);
                      p->velocity[0] = V[0];
                      p->velocity[1] = V[1];
                      p->velocity[2] = V[2];
                    }
                  });

                  if (net->RK_Q.finished)
                    mooring->applyMooringSimulationResult();

                  for (const auto& p : mooring->getPoints())
                    mooring->lastPoint->setX(nextPositionOnBody(net, mooring->lastPoint));
                }
              }

              // 人工的な粘性で結果が一致するようになるかどうかチェックする．

              bool use_given_velocity = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] != "update" && net->inputJSON.at("velocity")[0] != "floating");
              bool update_velocity_using_predetermined_accel = (net->inputJSON.find("velocity") && net->inputJSON.find("acceleration") && net->inputJSON.at("velocity")[0] == "update");
              bool update_velocity_using_solved_accel = (net->inputJSON.find("acceleration") && net->inputJSON.at("acceleration")[0] == "floating");
              bool update_velocity_using_solved_accel2 = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating");
              if (use_given_velocity) {
                std::cout << "use " << net->getName() << "'s (Rigid/SoftBodyObject) predetermiend velocity" << std::endl;
              } else if (update_velocity_using_solved_accel || update_velocity_using_predetermined_accel || update_velocity_using_solved_accel2) {
                std::cout << "updating " << net->getName() << "'s (Rigid/SoftBodyObject) velocity from acceleration" << std::endl;
                net->velocity = net->RK_Velocity.getX();
                std::cout << "acceleration = " << net->acceleration << std::endl;
                std::cout << "velocity = " << net->velocity << std::endl;
              } else {
                std::cout << net->getName() << "'s (Rigid/SoftBodyObject) velocity is not updated" << std::endl;
              }
            }

            //% --- phi_t BVP (floating body coupling) ---
            std::vector<Network*> movableObjects;
            for (auto& net : Join(RigidBodyObject, SoftBodyObject))
              if (std::ranges::any_of(net->isFixed, [](const auto& v) { return v == false; })) {
                movableObjects.push_back(net);
              }

            if (coupling_type == "NONE" || movableObjects.empty()) {
              convergence.clear();
              if (bem_verbose())
                PrintLap(watch, "Skipping BVP.solveForPhiPhin_t (no movable bodies or coupling=NONE).");
            } else {
              convergence = BVP.solveForPhiPhin_t(movableObjects);
              if (bem_verbose())
                PrintLap(watch, "BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる");
            }

            /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

            ### 浮体の重心位置・姿勢・速度の更新

            浮体の重心位置は，重心に関する運動方程式を解くことで求める．
            姿勢は，角運動量に関する運動方程式などを使って，各加速度を求める．姿勢はクオータニオンを使って表現する．

            */

            //% --- Body update: post-solve push (RK advance, geometry update) ---
            _Pragma("omp parallel for") for (const auto& net : Join(RigidBodyObject, SoftBodyObject)) {
              if (net->inputJSON.find("velocity")) {
                //! 時間まで静止させる
                if ((net->inputJSON.at("velocity")[0].contains("floating") && net->inputJSON.at("velocity").size() >= 2 && std::stod(net->inputJSON.at("velocity")[1]) > simulation_time)) {
                  net->velocity.fill(0);
                  net->acceleration.fill(0);
                }

                //! 静止した物体は速度ゼロが与えられているので，下を実行しても動かない
                if (!net->inputJSON.at("velocity")[0].contains("fixed")) {
                  std::cout << "updating " << net->getName() << "'s (RigidBodyObject) velocity" << std::endl;
                  auto V = net->velocityTranslational();
                  auto W = net->Q.AngularVelocityTodQdt(net->velocityRotational());
                  for (auto i = 0; i < 3; i++) {
                    if (net->isFixed[i])
                      V[i] = 0;
                    if (net->isFixed[i + 3])
                      W[i] = 0;
                  }
                  net->RK_COM.push(V);
                  net->RK_Q.push(W);
                }
                std::cout << "name = " << net->getName() << std::endl;
                std::cout << "net->velocityTranslational() = " << net->velocityTranslational() << std::endl;

                for (auto& mooring : net->mooringLines) {
                  //! mooring->lastPoint は，浮体ともに動く．
                  auto Xcurrent = mooring->lastPoint->X;
                  mooring->lastPoint->X_last = Xcurrent;
                  Tddd V = (nextPositionOnBody(net, mooring->lastPoint) - Xcurrent) / (net->RK_COM.getTimeAtNextStep() - simulation_time);
                  mooring->simulate(simulation_time, net->RK_COM.getTimeAtNextStep() - simulation_time, [&](networkPoint* p) {
                    if (p == mooring->firstPoint) {
                      p->acceleration.fill(0);
                      p->velocity.fill(0);
                    } else if (p == mooring->lastPoint) {
                      p->acceleration.fill(0);
                      p->velocity[0] = V[0];
                      p->velocity[1] = V[1];
                      p->velocity[2] = V[2];
                    }
                  });

                  if (net->RK_Q.finished)
                    mooring->applyMooringSimulationResult();

                  for (const auto& p : mooring->getPoints())
                    mooring->lastPoint->setX(nextPositionOnBody(net, mooring->lastPoint));
                }
              }
              // 人工的な粘性で結果が一致するようになるかどうかチェックする．
              bool use_given_velocity = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] != "update" && net->inputJSON.at("velocity")[0] != "floating");
              bool update_velocity_using_predetermined_accel = (net->inputJSON.find("velocity") && net->inputJSON.find("acceleration") && net->inputJSON.at("velocity")[0] == "update");
              bool update_velocity_using_solved_accel = (net->inputJSON.find("acceleration") && net->inputJSON.at("acceleration")[0] == "floating");
              bool update_velocity_using_solved_accel2 = (net->inputJSON.find("velocity") && net->inputJSON.at("velocity")[0] == "floating");
              if (use_given_velocity) {
                std::cout << "use " << net->getName() << "'s (RigidBodyObject) predetermiend velocity" << std::endl;
              } else if (update_velocity_using_solved_accel || update_velocity_using_predetermined_accel || update_velocity_using_solved_accel2) {
                std::cout << "updating " << net->getName() << "'s (RigidBodyObject) velocity from acceleration" << std::endl;
                net->RK_Velocity.push(net->acceleration);
                std::cout << "acceleration = " << net->acceleration << std::endl;
                std::cout << "velocity = " << net->velocity << std::endl;
              } else {
                std::cout << net->getName() << "'s (RigidBodyObject) velocity is not updated" << std::endl;
              }

              // Commit the newly advanced rigid-body RK pose before updating boundary points.
              if (net->RK_COM.steps > 0)
                net->COM = net->RK_COM.getX();
              if (net->RK_Q.steps > 0)
                net->Q = Normalize(net->RK_Q.getX());
              if (net->RK_Velocity.steps > 0)
                net->velocity = net->RK_Velocity.getX();

              if (net->isRigidBody) {
                for (const auto& p : net->getPoints())
                  p->setXSingle(net->rigidTransformation(p->initialX));
              } else {
                for (const auto& p : net->getPoints()) {
                  p->RK_X.push(p->velocityTranslational()); // setBodyVelocityで設定された速度を使う
                  p->setXSingle(p->RK_X.getX());
                }
              }

              net->setGeometricPropertiesForce();
            }

            //& --- Absorption + RK push ---
            computeSignedDistances(fluid_nodes);
            const double mean_phi = computeMeanPhi(Buckets_Fluid_Faces.data1D);
            applyAbsorptionAndPush(fluid_nodes, mean_phi,
                                   node_relocation_method == NodeRelocationMethod::ALE);

            //& --- Inactive line: 頂点の線形補間で値を維持 ---
            for (auto* water : FluidObject)
              for (auto* l : water->getBoundaryLines())
                if (!l->hasActiveBieDof()) {
                  auto [pA, pB] = l->getPoints();
                  l->setXSingle(0.5 * (pA->X + pB->X));
                  l->phiphin[0] = 0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
                }

            //! --- Quality checks (face altitude, tiny faces) ---
            for (auto net : AllObjects)
              net->setGeometricPropertiesForce();
            if (!retry_state.degraded_mode) {
              for (auto water : FluidObject)
                refreshFaceBadQualityHistory(*water, time_step, RK_step, 0.1, subsurface_altitude_reject);
              for (auto water : FluidObject)
                monitorTinyFaceRelativeToLocalMean(*water, time_step, RK_step);
              for (auto water : FluidObject)
                throwIfSubsurfaceFaceAltitudeTooSmall(*water, time_step, RK_step, subsurface_altitude_reject);
              for (auto water : FluidObject)
                throwIfTinyFaceRelativeToLocalMean(*water, time_step, RK_step, 0.01);
            }

            for (auto water : FluidObject) {
              auto name = water->getName();
              std::ofstream ofs(output_directory / (name + std::to_string(RK_step) + ".obj"));
              createOBJ(ofs, *water);
              ofs.close();
            }

            if (bem_verbose())
              std::cout << Blue << "Elapsed time: " << Red << watch() << colorReset << " s\n";
          } while (!((*FluidObject[0]->getPoints().begin())->RK_X.finished));

          //@ --- 3c-2f. Post-RK: interpolation relocation, VPM step, mean-phi ---
          sync_bodies();
          for (auto net : AllObjects)
            net->setGeometricPropertiesForce();
          BEM_Penetration::throwIfStructurePenetrated(FluidObject, Join(RigidBodyObject, SoftBodyObject), time_step, "post-RK", use_true_quadratic_element);

          // POST-RK4 interpolation relocation:
          // use X_reloc computed by setNodeVelocity(), project it onto the
          // current RK-updated mesh, interpolate phi there, then move nodes.
          if (!use_mesh_preparation_pipeline &&
              node_relocation_method == NodeRelocationMethod::interpolation &&
              time_step % node_relocation_period == 0) {
            applyInterpolationRelocation(FluidObject, use_true_quadratic_element, node_relocation_surface, interpolation_midpoint_mode);
            if (bem_verbose())
              PrintLap(watch, "Interpolation relocation completed");
          }

          std::cout << Red << "Total elapsed time: " << Red << watch()[1] << colorReset << " s\n";
          ElapsedTimeTotal = watch()[1];

          simulation_time = (*FluidObject[0]->getPoints().begin())->RK_X.gett();

          //^ --- VPM post-RK step (shedding, diffusion half-step, particle cleanup) ---
          // Contract: setBodyVelocity() must be called before postRKStep().
          setBodyVelocity(Join(RigidBodyObject, SoftBodyObject));
          vpm.postRKStep(dt, FluidObject, RigidBodyObject);
          vpm.refreshKinematicsForDt([&](const std::array<double, 3>& x) { return getBEMVelocityAt_cached(x, FluidObject); });

          subtractMeanPhi(Buckets_Fluid_Faces.data1D, fluid_nodes);

          for (const auto& net : RigidBodyObject)
            net->Q.normalize();

          //@ --- 3c-2g. Fixed-body pressure approximation (finite difference phi_t) ---
          {
            bool has_movable_body = false;
            for (const auto* net : Join(RigidBodyObject, SoftBodyObject))
              if (std::ranges::any_of(net->isFixed, [](const auto& v) { return v == false; })) {
                has_movable_body = true;
                break;
              }

            if (!has_movable_body) {
              static std::unordered_map<networkPoint*, double> prev_phi;
              for (auto* water : FluidObject) {
                for (auto* p : water->getBoundaryPoints()) {
                  const double phi = std::get<0>(p->phiphin);
                  double phi_t = 0.0;
                  if (dt > 0.0) {
                    if (auto it = prev_phi.find(p); it != prev_phi.end())
                      phi_t = (phi - it->second) / dt;
                  }
                  prev_phi[p] = phi;

                  std::get<0>(p->phiphin_t) = phi_t;
                  std::get<1>(p->phiphin_t) = 0.0;

                  if (p->Dirichlet)
                    p->pressure = p->pressure_BEM = 0.0;
                  else
                    p->pressure = p->pressure_BEM = -_WATER_DENSITY_ * (phi_t + 0.5 * Dot(p->u_total, p->u_total) + _GRAVITY_ * p->height());
                }
              }
            }
          }

          //* Step succeeded
          last_successful_dt = dt;
          break;

        } catch (step_failure& e) {
          if (Network::scalingDebugEnabled()) {
            std::cout << Magenta << "[scale] " << Cyan << "catch step_failure" << colorReset << std::endl;
            for (auto* net : FluidObject)
              net->debugPrintScalingState("time_loop catch step_failure");
          }
          if (retry_state.handle_step_failure(e, step_retry, dt, last_successful_dt, time_step, write_failure_snapshot) == RetryAction::BreakStep)
            break;
          continue;
        } catch (const std::exception& e) {
          if (Network::scalingDebugEnabled()) {
            std::cout << Magenta << "[scale] " << Cyan << "catch std::exception" << colorReset << std::endl;
            for (auto* net : FluidObject)
              net->debugPrintScalingState("time_loop catch std::exception");
          }
          step_failure wrapped(std::string("[firewall] unexpected exception in time-step body: ") + e.what());
          if (retry_state.handle_step_failure(wrapped, step_retry, dt, last_successful_dt, time_step, write_failure_snapshot) == RetryAction::BreakStep)
            break;
          continue;
        } catch (...) {
          if (Network::scalingDebugEnabled()) {
            std::cout << Magenta << "[scale] " << Cyan << "catch unknown exception" << colorReset << std::endl;
            for (auto* net : FluidObject)
              net->debugPrintScalingState("time_loop catch unknown exception");
          }
          step_failure wrapped("[firewall] unknown exception in time-step body");
          if (retry_state.handle_step_failure(wrapped, step_retry, dt, last_successful_dt, time_step, write_failure_snapshot) == RetryAction::BreakStep)
            break;
          continue;
        }
      } //! end retry loop
      if (retry_state.degraded_mode) {
        std::cerr << Yellow << "[step_reject] Exiting DEGRADED MODE for next time step." << colorReset << std::endl;
        retry_state.degraded_mode = false;
      }

      //$ --- 3c-3. Step output ---
      std::filesystem::remove(step_checkpoint_path);
      write_step_outputs(false);
    } //# end time loop
  } catch (const error_message& e) {
    e.print(); // 色付き表示 to cerr
    return 1;
  } catch (const std::exception& e) {
    std::cerr << Red << e.what() << colorReset << std::endl;
    return 1;
  }
  return 0;
};

/*DOC_EXTRACT 2_0_0_HOW_TO_RUN

# 実行方法

## ファイルのダウンロード

上書きされるので注意．ダウンロードしたら，`build_bem`ディレクトリに移動．

```sh
git clone https://github.com/tomoakihirakawa/cpp.git
cd ./cpp/builds/build_bem
```

## 入力ファイルの生成．

```sh
python3 input_generator.py
```

例えば，`./input_files/Hadzic2005`が生成される．

## プログラムのコンパイルと実行

`clean`でCMake関連のファイルを削除して（ゴミがあるかもしれないので），
`cmake`で`Makefile`を生成して，`make`でコンパイルする．

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

実行

```sh
./main ./input_files/Hadzic2005
```

*/

/*DOC_EXTRACT 3_0_EXAMPLES

# Examples

**[See the Examples here!](EXAMPLES.md)**

*/
