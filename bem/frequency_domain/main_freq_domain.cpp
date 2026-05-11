#include "pch.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BEM_time_domain_types.hpp"

// Global variables expected by legacy headers.
bool use_linear_element = false;
bool use_pseudo_quadratic_element = false;
bool use_true_quadratic_element = false;
bool use_quadratic_linear_hybrid = false;
NodeRelocationMethod node_relocation_method = NodeRelocationMethod::none;
NodeRelocationSurface node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
std::string solver_type = "GMRES";
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
bool enable_pressure_detachment = false;
int detachment_consecutive_steps = 3;
std::string nearfield_mode = "scalar";
int g_p2m_quadrature_points = 6;
int g_lu_far_dunavant_points = 0;
int g_lu_near_dunavant_points = 0;
double g_mac_theta = 0.25;

int fmm_max_level = 7;
int fmm_bucket_max_points = 50;

int time_step = 0;
double simulation_time = 0.0;

#define BEM
#include "Network.hpp"

#include "BEM_freqency_domain.hpp"
#include "BEM_inputfile_reader.hpp"
#include "BEM_pre_bvp_consistency.hpp"
#include "BEM_qtf.hpp"
#include "BEM_setBoundaryTypes.hpp"
#include "BEM_solveBVP.hpp"

namespace {

using Complex = std::complex<double>;
using Matrix6 = std::array<std::array<double, 6>, 6>;
using Matrix6Mask = std::array<std::array<bool, 6>, 6>;
using CVector6 = std::array<Complex, 6>;

struct BBox {
  Tddd min{{1e100, 1e100, 1e100}};
  Tddd max{{-1e100, -1e100, -1e100}};
};

BBox compute_bbox(const Network& net) {
  BBox b;
  for (auto* p : net.getPoints()) {
    b.min[0] = std::min(b.min[0], p->X[0]);
    b.min[1] = std::min(b.min[1], p->X[1]);
    b.min[2] = std::min(b.min[2], p->X[2]);
    b.max[0] = std::max(b.max[0], p->X[0]);
    b.max[1] = std::max(b.max[1], p->X[1]);
    b.max[2] = std::max(b.max[2], p->X[2]);
  }
  return b;
}

bool is_close(double a, double b, double eps) { return std::abs(a - b) <= eps; }

bool has_flag(int argc, char** argv, const std::string& flag) {
  for (int i = 2; i < argc; ++i) {
    if (argv[i] == flag)
      return true;
  }
  return false;
}

struct FrequencySpongeRuntime {
  enum class Mode {
    Radial,
    WallDistance,
    AbsorberWavelength,
    AbsorberWallMax
  };
  bool enabled = false;
  Mode mode = Mode::Radial;
  double r_start = 0.0;
  double length = 0.0;
  double mu_max = 0.0;
  int order = 2;
  double n_wavelengths = 3.0;
  std::string profile = "power";
  double floor_fraction = 0.0;
  double wall_mu_max = 0.0;
  double wall_width = 0.0;
  int wall_order = 2;
};

std::string frequency_sponge_mode_name(FrequencySpongeRuntime::Mode mode) {
  switch (mode) {
  case FrequencySpongeRuntime::Mode::Radial:
    return "radial";
  case FrequencySpongeRuntime::Mode::WallDistance:
    return "wall_distance";
  case FrequencySpongeRuntime::Mode::AbsorberWavelength:
    return "absorber_wavelength";
  case FrequencySpongeRuntime::Mode::AbsorberWallMax:
    return "absorber_wall_max";
  default:
    return "unknown";
  }
}

double finite_depth_wave_number(double omega, double depth, double gravity) {
  if (!(omega > 0.0) || !(gravity > 0.0))
    throw std::runtime_error("finite_depth_wave_number: omega and gravity must be positive");
  if (!(depth > 0.0) || !std::isfinite(depth))
    throw std::runtime_error("finite_depth_wave_number: depth must be positive and finite");
  auto f = [&](double k) {
    return gravity * k * std::tanh(k * depth) - omega * omega;
  };
  double lo = 0.0;
  double hi = std::max(omega * omega / gravity, omega / std::sqrt(gravity * depth));
  hi = std::max(hi, 1e-9);
  while (f(hi) < 0.0)
    hi *= 2.0;
  for (int iter = 0; iter < 100; ++iter) {
    const double mid = 0.5 * (lo + hi);
    if (f(mid) < 0.0)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

double smoothstep01(double s) {
  s = std::clamp(s, 0.0, 1.0);
  return s * s * (3.0 - 2.0 * s);
}

double sponge_profile_value_with_order(double s, const FrequencySpongeRuntime& sponge, int order_in) {
  s = std::clamp(s, 0.0, 1.0);
  const double order = static_cast<double>(order_in);
  if (sponge.profile == "power")
    return std::pow(s, order);
  if (sponge.profile == "smoothstep")
    return smoothstep01(s);
  if (sponge.profile == "smoothstep_power")
    return std::pow(smoothstep01(s), order);
  if (sponge.profile == "floor_power") {
    const double floor = std::clamp(sponge.floor_fraction, 0.0, 1.0);
    return floor + (1.0 - floor) * std::pow(s, order);
  }
  if (sponge.profile == "floor_smoothstep_power") {
    const double floor = std::clamp(sponge.floor_fraction, 0.0, 1.0);
    return floor + (1.0 - floor) * std::pow(smoothstep01(s), order);
  }
  throw std::runtime_error("unknown frequency sponge profile: " + sponge.profile);
}

double sponge_profile_value(double s, const FrequencySpongeRuntime& sponge) {
  return sponge_profile_value_with_order(s, sponge, sponge.order);
}

struct FaceSets {
  std::unordered_set<networkFace*> free_surface;  // Robin
  std::unordered_set<networkFace*> float_surface; // body boundary (radiation forcing)
  std::unordered_set<networkFace*> outer_wall;    // finite tank side wall
  std::unordered_set<networkFace*> bottom;        // tank bottom
};

std::string node_face_bc_name(NodeFaceBoundaryType t) {
  switch (t) {
  case NodeFaceBoundaryType::Dirichlet:
    return "dirichlet";
  case NodeFaceBoundaryType::Neumann:
    return "neumann";
  case NodeFaceBoundaryType::Undefined:
  default:
    return "undefined";
  }
}

std::string safe_number_label(double v) {
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(6) << v;
  std::string s = ss.str();
  for (char& c : s) {
    if (c == '-')
      c = 'm';
    else if (c == '.')
      c = 'p';
  }
  return s;
}

std::string summarized_bc_name(SummarizedNodeBoundaryType t) {
  switch (t) {
  case SummarizedNodeBoundaryType::Dirichlet:
    return "dirichlet";
  case SummarizedNodeBoundaryType::Neumann:
    return "neumann";
  case SummarizedNodeBoundaryType::Interface:
    return "interface";
  case SummarizedNodeBoundaryType::Invalid:
  default:
    return "invalid";
  }
}

std::string frequency_expected_bc_name(const networkFace* f, const FaceSets& face_sets) {
  if (!f)
    return "undefined";
  return face_sets.free_surface.contains(const_cast<networkFace*>(f)) ? "dirichlet" : "neumann";
}

std::string face_kind_name(const networkFace* f, const FaceSets& face_sets) {
  if (!f)
    return "none";
  if (face_sets.float_surface.contains(const_cast<networkFace*>(f)))
    return "float";
  if (face_sets.free_surface.contains(const_cast<networkFace*>(f)))
    return "free";
  if (face_sets.outer_wall.contains(const_cast<networkFace*>(f)))
    return "outer_wall";
  if (face_sets.bottom.contains(const_cast<networkFace*>(f)))
    return "bottom";
  return "wall";
}

FaceSets classify_faces_deepcwind(Network& water, const Network& float_body) {
  FaceSets sets;

  water.setGeometricPropertiesForce();

  const auto water_bounds = water.bounds;
  const double x_min = std::get<0>(water_bounds[0]);
  const double x_max = std::get<1>(water_bounds[0]);
  const double y_min = std::get<0>(water_bounds[1]);
  const double y_max = std::get<1>(water_bounds[1]);
  const double z_min = std::get<0>(water_bounds[2]);
  const double z_max = std::get<1>(water_bounds[2]);

  for (auto* f : water.getBoundaryFaces()) {
    const auto& c = f->centroid;
    const auto& n = f->normal;

    const bool on_free_surface = is_close(c[2], z_max, 1e-6) && (n[2] > 0.9);
    if (on_free_surface) {
      sets.free_surface.emplace(f);
      continue;
    }

    const bool on_outer_x = is_close(std::abs(c[0]), std::max(std::abs(x_min), std::abs(x_max)), 1e-6);
    const bool on_outer_y = is_close(std::abs(c[1]), std::max(std::abs(y_min), std::abs(y_max)), 1e-6);
    const bool on_bottom = is_close(c[2], z_min, 1e-5) && (n[2] < -0.9);
    if (on_outer_x || on_outer_y) {
      sets.outer_wall.emplace(f);
      continue;
    }
    if (on_bottom) {
      sets.bottom.emplace(f);
      continue; // tank walls/bottom
    }

    // Any remaining face (not free surface, not tank wall/bottom) is float body surface.
    sets.float_surface.emplace(f);
  }

  if (sets.free_surface.empty())
    throw std::runtime_error("classify_faces_deepcwind: free_surface faces not found");
  if (sets.float_surface.empty())
    throw std::runtime_error("classify_faces_deepcwind: float_surface faces not found");

  return sets;
}

template <class Entity>
std::string contact_body_names_for_node_face(const Entity* entity, const networkFace* f) {
  if (!entity || !f)
    return "";
  const auto* d = entity->findContactState(f);
  if (!d)
    return "";
  std::vector<std::string> names;
  for (auto* cf : d->contact_opponent_faces) {
    if (!cf || !cf->getNetwork())
      continue;
    const std::string name = cf->getNetwork()->getName();
    if (std::ranges::find(names, name) == names.end())
      names.push_back(name);
  }
  std::sort(names.begin(), names.end());
  std::ostringstream oss;
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (i)
      oss << "|";
    oss << names[i];
  }
  return oss.str();
}

struct FrequencyBcCheckSummary {
  std::size_t time_active_ids = 0;
  std::size_t freq_active_ids = 0;
  std::size_t point_face_pairs = 0;
  std::size_t line_face_pairs = 0;
  std::size_t node_face_mismatches = 0;
  std::size_t face_count = 0;
  std::size_t face_mismatches = 0;
  std::size_t float_contact_faces = 0;
  std::size_t free_faces = 0;
  std::size_t wall_faces = 0;
  BEMPreBVP::Stats pre_bvp_stats;
};

FrequencyBcCheckSummary write_frequency_contact_bc_check(const std::filesystem::path& outdir,
                                                         SimulationSettings& setting,
                                                         Network* water,
                                                         Network* float_body,
                                                         const FaceSets& face_sets,
                                                         const bem_frequency_domain::RobinUnknownPolicy radiation_robin_unknown_policy,
                                                         const bool radiation_face_neumann) {
  if (!water || !float_body)
    throw std::runtime_error("write_frequency_contact_bc_check: water/float body is null");

  std::filesystem::create_directories(outdir);
  FrequencyBcCheckSummary summary;
  summary.free_faces = face_sets.free_surface.size();
  summary.float_contact_faces = face_sets.float_surface.size();
  for (auto* f : water->getBoundaryFaces()) {
    if (!face_sets.free_surface.contains(f) && !face_sets.float_surface.contains(f))
      ++summary.wall_faces;
  }

  const std::vector<Network*> contact_objects = Join(setting.RigidBodyObject, setting.SoftBodyObject);
  const std::vector<Network*> bucket_objects = Join(setting.FluidObject, setting.RigidBodyObject, setting.SoftBodyObject, setting.AbsorberObject);

  auto prepare_buckets_like_time_domain = [&]() {
    for (auto* net : bucket_objects) {
      if (!net)
        continue;
      net->setGeometricPropertiesForce();
      net->makeBuckets(net->getScale() / 10.);
    }
  };

  auto expected_freq_bc = [&](const networkFace* f) {
    if (face_sets.free_surface.contains(const_cast<networkFace*>(f))) {
      return radiation_robin_unknown_policy == bem_frequency_domain::RobinUnknownPolicy::Phi
                 ? NodeFaceBoundaryType::Neumann
                 : NodeFaceBoundaryType::Dirichlet;
    }
    return NodeFaceBoundaryType::Neumann;
  };
  auto count_active_line_midpoint_dofs = [&]() {
    std::size_t count = 0;
    for (auto* l : water->getBoundaryLines()) {
      if (!l)
        continue;
      if (l->midpoint_index >= 0)
        ++count;
    }
    return count;
  };

  {
    bem_frequency_domain::BoundaryStateGuard guard(setting.FluidObject);
    prepare_buckets_like_time_domain();
    refreshBoundaryStatesAndTypes(setting.FluidObject, contact_objects);
    summary.time_active_ids = setNodeFaceIndices(setting.FluidObject);
    const auto time_active_line_midpoint_ids = count_active_line_midpoint_dofs();

    BEMPreBVP::Options options;
    options.print_ok_summary = true;
    options.print_waterline_midpoint_quality = use_true_quadratic_element;
    std::ostringstream pre_bvp_log;
    summary.pre_bvp_stats = BEMPreBVP::inspect(setting.FluidObject, 0, 0, options, pre_bvp_log);

    {
      std::ofstream fs(outdir / "frequency_bc_contact_summary.csv");
      fs << "key,value\n";
      fs << "contact_objects";
      fs << ",";
      for (std::size_t i = 0; i < contact_objects.size(); ++i) {
        if (i)
          fs << "|";
        fs << (contact_objects[i] ? contact_objects[i]->getName() : "null");
      }
      fs << "\n";
      fs << "bucket_objects";
      fs << ",";
      for (std::size_t i = 0; i < bucket_objects.size(); ++i) {
        if (i)
          fs << "|";
        fs << (bucket_objects[i] ? bucket_objects[i]->getName() : "null");
      }
      fs << "\n";
      fs << "time_active_ids," << summary.time_active_ids << "\n";
      fs << "time_active_line_midpoint_ids," << time_active_line_midpoint_ids << "\n";
      fs << "pre_bvp_unclassified_points," << summary.pre_bvp_stats.unclassified_points << "\n";
      fs << "pre_bvp_unclassified_lines," << summary.pre_bvp_stats.unclassified_lines << "\n";
      fs << "pre_bvp_bcinterface_points_no_contact," << summary.pre_bvp_stats.bcinterface_points_no_contact << "\n";
      fs << "pre_bvp_bcinterface_points," << summary.pre_bvp_stats.bcinterface_points << "\n";
      fs << "pre_bvp_bcinterface_lines_no_contact," << summary.pre_bvp_stats.bcinterface_lines_no_contact << "\n";
      fs << "pre_bvp_bcinterface_lines," << summary.pre_bvp_stats.bcinterface_lines << "\n";
      fs << "pre_bvp_nonfinite_dofs," << summary.pre_bvp_stats.nonfinite_dofs << "\n";
      fs << "pre_bvp_overflow_dofs," << summary.pre_bvp_stats.overflow_dofs << "\n";
      fs << "pre_bvp_has_warning," << (summary.pre_bvp_stats.has_warning() ? 1 : 0) << "\n";
      fs << "free_faces," << summary.free_faces << "\n";
      fs << "float_faces," << summary.float_contact_faces << "\n";
      fs << "wall_faces," << summary.wall_faces << "\n";
    }

    {
      std::ofstream fs(outdir / "frequency_bc_pre_bvp_guard.log");
      fs << pre_bvp_log.str();
      for (const auto& ex : summary.pre_bvp_stats.examples)
        fs << ex << "\n";
    }

    {
      std::ofstream fs(outdir / "frequency_bc_node_face_comparison.csv");
      fs << "entity_type,entity_id,face_index,x,y,z,face_kind,time_node_face_bc,freq_expected_bc,"
            "mismatch,contact_bodies,is_bc_interface,is_multiple_node,active_index\n";
      auto write_entity = [&](const char* entity_type, const auto* entity, int entity_id) {
        if (!entity)
          return;
        for (auto* f : entity->getBoundaryFaces()) {
          if (!f)
            continue;
          const auto time_bc = getNodeFaceBoundaryType(entity, f);
          const auto freq_bc = expected_freq_bc(f);
          const bool mismatch = time_bc != freq_bc;
          if (mismatch)
            ++summary.node_face_mismatches;
          const auto* d = entity->findActiveBieDofOrDefault(f);
          fs << entity_type << "," << entity_id << "," << f->index << ","
             << std::scientific << std::setprecision(12)
             << entity->getPosition()[0] << "," << entity->getPosition()[1] << "," << entity->getPosition()[2] << ","
             << face_kind_name(f, face_sets) << ","
             << node_face_bc_name(time_bc) << "," << node_face_bc_name(freq_bc) << ","
             << (mismatch ? 1 : 0) << ","
             << contact_body_names_for_node_face(entity, f) << ","
             << (entity->BCInterface ? 1 : 0) << ","
             << (entity->isMultipleNode ? 1 : 0) << ","
             << (d ? d->index : -1) << "\n";
        }
      };
      for (auto* p : water->getBoundaryPoints()) {
        summary.point_face_pairs += p ? p->getBoundaryFaces().size() : 0;
        write_entity("point", p, p ? p->index : -1);
      }
      for (auto* l : water->getBoundaryLines()) {
        summary.line_face_pairs += l ? l->getBoundaryFaces().size() : 0;
        auto [p0, p1] = l->getPoints();
        const int id = (p0 ? p0->index : -1);
        write_entity("line", l, id);
        (void)p1;
      }
    }

    {
      std::ofstream fs(outdir / "frequency_bc_face_comparison.csv");
      fs << "face_index,face_kind,cx,cy,cz,nx,ny,nz,area,time_point_neumann,time_point_dirichlet,"
            "time_line_neumann,time_line_dirichlet,freq_expected_bc,mismatch\n";
      for (auto* f : water->getBoundaryFaces()) {
        if (!f)
          continue;
        ++summary.face_count;
        int p_neu = 0, p_dir = 0, l_neu = 0, l_dir = 0;
        for (auto* p : f->getPoints()) {
          const auto bc = getNodeFaceBoundaryType(p, f);
          p_neu += (bc == NodeFaceBoundaryType::Neumann) ? 1 : 0;
          p_dir += (bc == NodeFaceBoundaryType::Dirichlet) ? 1 : 0;
        }
        for (auto* l : f->Lines) {
          const auto bc = getNodeFaceBoundaryType(l, f);
          l_neu += (bc == NodeFaceBoundaryType::Neumann) ? 1 : 0;
          l_dir += (bc == NodeFaceBoundaryType::Dirichlet) ? 1 : 0;
        }
        const auto freq_bc = expected_freq_bc(f);
        const bool all_nodes_match =
            (freq_bc == NodeFaceBoundaryType::Neumann) ? (p_dir == 0 && l_dir == 0) : (p_neu == 0 && l_neu == 0);
        if (!all_nodes_match)
          ++summary.face_mismatches;
        fs << f->index << "," << face_kind_name(f, face_sets) << ","
           << std::scientific << std::setprecision(12)
           << f->centroid[0] << "," << f->centroid[1] << "," << f->centroid[2] << ","
           << f->normal[0] << "," << f->normal[1] << "," << f->normal[2] << ","
           << f->area << ","
           << p_neu << "," << p_dir << "," << l_neu << "," << l_dir << ","
           << node_face_bc_name(freq_bc) << "," << (all_nodes_match ? 0 : 1) << "\n";
      }
    }

    guard.restore();
  }

  {
    bem_frequency_domain::BoundaryStateGuard guard(setting.FluidObject);
    bem_frequency_domain::BoundaryData freq_bc_data;
    freq_bc_data.face_bc = [&](const networkFace& f) -> bem_frequency_domain::FaceBC {
      return face_sets.free_surface.contains(const_cast<networkFace*>(&f)) ? bem_frequency_domain::FaceBC::Robin
                                                                           : bem_frequency_domain::FaceBC::Neumann;
    };
    freq_bc_data.robin_unknown_policy = radiation_robin_unknown_policy;
    if (radiation_face_neumann) {
      freq_bc_data.force_neumann_multiple = [&](const BEM_DOF_Base& node) -> bool {
        return std::ranges::any_of(node.getBoundaryFaces(), [&](const networkFace* f) {
          return face_sets.float_surface.contains(const_cast<networkFace*>(f));
        });
      };
    }
    std::unordered_set<networkFace*> robin_faces;
    bem_frequency_domain::apply_face_bc_to_mesh(setting.FluidObject, freq_bc_data, robin_faces);
    summary.freq_active_ids = setNodeFaceIndices(setting.FluidObject, [&](const BEM_DOF_Base* node) {
      return node && freq_bc_data.force_neumann_multiple && freq_bc_data.force_neumann_multiple(*node);
    });
    const auto freq_active_line_midpoint_ids = count_active_line_midpoint_dofs();
    std::ofstream fs(outdir / "frequency_bc_contact_summary.csv", std::ios::app);
    fs << "freq_active_line_midpoint_ids," << freq_active_line_midpoint_ids << "\n";
    guard.restore();
  }

  {
    std::ofstream fs(outdir / "frequency_bc_contact_summary.csv", std::ios::app);
    fs << "freq_active_ids," << summary.freq_active_ids << "\n";
    fs << "node_face_point_entities," << summary.point_face_pairs << "\n";
    fs << "node_face_line_entities," << summary.line_face_pairs << "\n";
    fs << "node_face_mismatches," << summary.node_face_mismatches << "\n";
    fs << "face_count," << summary.face_count << "\n";
    fs << "face_mismatches," << summary.face_mismatches << "\n";
    fs << "active_id_match," << (summary.time_active_ids == summary.freq_active_ids ? 1 : 0) << "\n";
    fs << "bc_match," << (summary.node_face_mismatches == 0 && summary.face_mismatches == 0 ? 1 : 0) << "\n";
  }

  return summary;
}

struct WettedBodyMeshExportSummary {
  std::size_t vertices = 0;
  std::size_t faces = 0;
  std::size_t faces_with_float_contact = 0;
  std::size_t faces_without_float_contact = 0;
  double area = 0.0;
  BBox bbox;
};

template <class Entity>
bool node_face_touches_network(const Entity* entity, const networkFace* f, const Network* target) {
  if (!entity || !f || !target)
    return false;
  const auto* state = entity->findContactState(f);
  if (!state)
    return false;
  for (auto* cf : state->contact_opponent_faces) {
    if (cf && cf->getNetwork() == target)
      return true;
  }
  return false;
}

bool face_touches_network_by_boundary_nodes(const networkFace* f, const Network* target) {
  if (!f || !target)
    return false;
  for (auto* p : f->getPoints()) {
    if (node_face_touches_network(p, f, target))
      return true;
  }
  for (auto* l : f->Lines) {
    if (node_face_touches_network(l, f, target))
      return true;
  }
  return false;
}

WettedBodyMeshExportSummary export_wetted_body_mesh_from_frequency_bc(const std::filesystem::path& outdir,
                                                                       SimulationSettings& setting,
                                                                       Network* water,
                                                                       Network* float_body,
                                                                       const FaceSets& face_sets) {
  if (!water || !float_body)
    throw std::runtime_error("export_wetted_body_mesh_from_frequency_bc: water/float body is null");

  std::filesystem::create_directories(outdir);

  const std::vector<Network*> contact_objects = Join(setting.RigidBodyObject, setting.SoftBodyObject);
  const std::vector<Network*> bucket_objects = Join(setting.FluidObject, setting.RigidBodyObject, setting.SoftBodyObject, setting.AbsorberObject);

  auto prepare_buckets_like_time_domain = [&]() {
    for (auto* net : bucket_objects) {
      if (!net)
        continue;
      net->setGeometricPropertiesForce();
      net->makeBuckets(net->getScale() / 10.);
    }
  };

  std::vector<networkFace*> wetted_faces;
  WettedBodyMeshExportSummary summary;

  bem_frequency_domain::BoundaryStateGuard guard(setting.FluidObject);
  prepare_buckets_like_time_domain();
  refreshBoundaryStatesAndTypes(setting.FluidObject, contact_objects);

  for (auto* f : water->getBoundaryFaces()) {
    if (!f || !face_sets.float_surface.contains(f))
      continue;
    wetted_faces.push_back(f);
    summary.area += f->area;
    if (face_touches_network_by_boundary_nodes(f, float_body))
      ++summary.faces_with_float_contact;
    else
      ++summary.faces_without_float_contact;
    for (auto* p : f->getPoints()) {
      if (!p)
        continue;
      summary.bbox.min[0] = std::min(summary.bbox.min[0], p->X[0]);
      summary.bbox.min[1] = std::min(summary.bbox.min[1], p->X[1]);
      summary.bbox.min[2] = std::min(summary.bbox.min[2], p->X[2]);
      summary.bbox.max[0] = std::max(summary.bbox.max[0], p->X[0]);
      summary.bbox.max[1] = std::max(summary.bbox.max[1], p->X[1]);
      summary.bbox.max[2] = std::max(summary.bbox.max[2], p->X[2]);
    }
  }

  std::unordered_map<const networkPoint*, int> vertex_id;
  std::vector<const networkPoint*> vertices;
  auto get_vertex_id = [&](const networkPoint* p) {
    auto [it, inserted] = vertex_id.emplace(p, static_cast<int>(vertices.size()) + 1);
    if (inserted)
      vertices.push_back(p);
    return it->second;
  };
  for (auto* f : wetted_faces) {
    for (auto* p : f->getPoints())
      get_vertex_id(p);
  }
  summary.vertices = vertices.size();
  summary.faces = wetted_faces.size();

  auto write_obj = [&](const std::filesystem::path& path, const bool flip_order) {
    std::ofstream fs(path);
    fs << "# Wetted body mesh exported from frequency-domain BC/contact classification.\n";
    fs << "# Source: water-side float contact faces. flip_order=" << (flip_order ? 1 : 0) << "\n";
    fs << "g float_wetted_from_bc\n";
    fs << std::scientific << std::setprecision(12);
    for (const auto* p : vertices)
      fs << "v " << p->X[0] << " " << p->X[1] << " " << p->X[2] << "\n";
    for (auto* f : wetted_faces) {
      std::vector<int> ids;
      for (auto* p : f->getPoints())
        ids.push_back(vertex_id.at(p));
      if (flip_order)
        std::reverse(ids.begin(), ids.end());
      fs << "f";
      for (int id : ids)
        fs << " " << id;
      fs << "\n";
    }
  };

  write_obj(outdir / "float_wetted_from_bc_body_normals.obj", true);
  write_obj(outdir / "float_wetted_from_bc_water_normals.obj", false);

  {
    std::ofstream fs(outdir / "float_wetted_from_bc_nemoh.dat");
    fs << "\t2\t0\n";
    fs << std::scientific << std::setprecision(12);
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      const auto* p = vertices[i];
      fs << (i + 1) << "\t" << p->X[0] << "\t" << p->X[1] << "\t" << p->X[2] << "\n";
    }
    fs << "0\t0\t0\t0\n";
    for (auto* f : wetted_faces) {
      std::vector<int> ids;
      for (auto* p : f->getPoints())
        ids.push_back(vertex_id.at(p));
      std::reverse(ids.begin(), ids.end()); // body outward normal for external solvers.
      if (ids.size() == 3)
        fs << ids[0] << "\t" << ids[1] << "\t" << ids[2] << "\t" << ids[2] << "\n";
      else if (ids.size() >= 4)
        fs << ids[0] << "\t" << ids[1] << "\t" << ids[2] << "\t" << ids[3] << "\n";
    }
    fs << "0\t0\t0\t0\n";
  }

  {
    std::ofstream fs(outdir / "wetted_body_mesh_summary.csv");
    fs << "key,value\n";
    fs << "source,frequency_domain_bc_contact\n";
    fs << "source_network,water\n";
    fs << "body_network," << float_body->getName() << "\n";
    fs << "body_normal_obj,float_wetted_from_bc_body_normals.obj\n";
    fs << "water_normal_obj,float_wetted_from_bc_water_normals.obj\n";
    fs << "nemoh_dat,float_wetted_from_bc_nemoh.dat\n";
    fs << "vertices," << summary.vertices << "\n";
    fs << "faces," << summary.faces << "\n";
    fs << "area," << std::setprecision(17) << summary.area << "\n";
    fs << "faces_with_float_contact," << summary.faces_with_float_contact << "\n";
    fs << "faces_without_float_contact," << summary.faces_without_float_contact << "\n";
    fs << "bbox_min_x," << summary.bbox.min[0] << "\n";
    fs << "bbox_max_x," << summary.bbox.max[0] << "\n";
    fs << "bbox_min_y," << summary.bbox.min[1] << "\n";
    fs << "bbox_max_y," << summary.bbox.max[1] << "\n";
    fs << "bbox_min_z," << summary.bbox.min[2] << "\n";
    fs << "bbox_max_z," << summary.bbox.max[2] << "\n";
    fs << "body_normal_face_order,flipped_from_water_surface\n";
    fs << "external_solver_note,NEMOH_Capytaine_use_body_wetted_surface_without_tank_water_absorber_sponge\n";
  }

  guard.restore();
  return summary;
}

std::vector<double> parse_omegas(int argc, char** argv, const std::vector<double>& settings_omegas) {
  // Default: a small sweep for a quick run.
  // Override: `--omega w1 w2 ...`
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == "--omega") {
      std::vector<double> w;
      for (int j = i + 1; j < argc; ++j) {
        if (std::string(argv[j]).starts_with("--"))
          break;
        w.push_back(std::stod(argv[j]));
      }
      if (w.empty())
        throw std::runtime_error("--omega provided but no values");
      return w;
    }
  }
  if (!settings_omegas.empty())
    return settings_omegas;
  return {0.5, 0.8, 1.1}; // rad/s
}

std::vector<int> parse_dofs(int argc, char** argv, const std::vector<int>& settings_dofs) {
  // Default: all 6.
  // Override: `--dofs 2 4` etc.
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == "--dofs") {
      std::vector<int> dofs;
      for (int j = i + 1; j < argc; ++j) {
        if (std::string(argv[j]).starts_with("--"))
          break;
        dofs.push_back(std::stoi(argv[j]));
      }
      if (dofs.empty())
        throw std::runtime_error("--dofs provided but no values");
      for (int& d : dofs) {
        if (d < 0 || d > 5)
          throw std::runtime_error("--dofs values must be in [0,5]");
      }
      std::sort(dofs.begin(), dofs.end());
      dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
      return dofs;
    }
  }
  if (!settings_dofs.empty()) {
    std::vector<int> dofs = settings_dofs;
    for (int d : dofs) {
      if (d < 0 || d > 5)
        throw std::runtime_error("settings dofs values must be in [0,5]");
    }
    std::sort(dofs.begin(), dofs.end());
    dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
    return dofs;
  }
  return {0, 1, 2, 3, 4, 5};
}

double parse_option_double(int argc, char** argv, const std::string& flag, double default_value) {
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == flag) {
      if (i + 1 >= argc || std::string(argv[i + 1]).starts_with("--"))
        throw std::runtime_error(flag + " requires one numeric value");
      return std::stod(argv[i + 1]);
    }
  }
  return default_value;
}

std::optional<double> parse_optional_double(int argc, char** argv, const std::string& flag) {
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == flag) {
      if (i + 1 >= argc || std::string(argv[i + 1]).starts_with("--"))
        throw std::runtime_error(flag + " requires one numeric value");
      return std::stod(argv[i + 1]);
    }
  }
  return std::nullopt;
}

std::string parse_option_string(int argc, char** argv, const std::string& flag, const std::string& default_value) {
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == flag) {
      if (i + 1 >= argc || std::string(argv[i + 1]).starts_with("--"))
        throw std::runtime_error(flag + " requires one value");
      return argv[i + 1];
    }
  }
  return default_value;
}

int parse_option_int(int argc, char** argv, const std::string& flag, int default_value) {
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == flag) {
      if (i + 1 >= argc || std::string(argv[i + 1]).starts_with("--"))
        throw std::runtime_error(flag + " requires one integer value");
      return std::stoi(argv[i + 1]);
    }
  }
  return default_value;
}

Tddd parse_option_tddd(int argc, char** argv, const std::string& flag, const Tddd& default_value) {
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == flag) {
      if (i + 3 >= argc)
        throw std::runtime_error(flag + " requires x y z");
      for (int k = 1; k <= 3; ++k) {
        if (std::string(argv[i + k]).starts_with("--"))
          throw std::runtime_error(flag + " requires x y z");
      }
      return {std::stod(argv[i + 1]), std::stod(argv[i + 2]), std::stod(argv[i + 3])};
    }
  }
  return default_value;
}

void require_sign_option(double value, const std::string& flag) {
  if (!(value == 1.0 || value == -1.0))
    throw std::runtime_error(flag + " must be +1 or -1");
}

Matrix6 make_zero_matrix6() {
  Matrix6 m{};
  for (auto& row : m)
    row.fill(0.0);
  return m;
}

Matrix6Mask make_false_matrix6_mask() {
  Matrix6Mask m{};
  for (auto& row : m)
    row.fill(false);
  return m;
}

double relative_pair_error(double a, double b) {
  const double scale = std::max({std::abs(a), std::abs(b), 1e-300});
  return std::abs(a - b) / scale;
}

Matrix6 symmetrize_radiation_matrix6(const Matrix6& raw, const Matrix6Mask& has_value) {
  Matrix6 sym = make_zero_matrix6();
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      const bool has_ij = has_value[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      const bool has_ji = has_value[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
      if (has_ij && has_ji)
        sym[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            0.5 * (raw[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                   raw[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)]);
      else if (has_ij)
        sym[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            raw[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      else if (has_ji)
        sym[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            raw[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
    }
  }
  return sym;
}

Matrix6 parse_row_major_matrix6(const std::vector<double>& values, const std::string& name) {
  if (values.size() != 36)
    throw std::runtime_error(name + " must contain 36 row-major values");
  Matrix6 m = make_zero_matrix6();
  for (int r = 0; r < 6; ++r)
    for (int c = 0; c < 6; ++c)
      m[r][c] = values[static_cast<std::size_t>(6 * r + c)];
  return m;
}

std::vector<int> nonfixed_dofs(const Network& body) {
  std::vector<int> out;
  for (int i = 0; i < 6; ++i) {
    if (!body.isFixed[static_cast<std::size_t>(i)])
      out.push_back(i);
  }
  return out;
}

std::array<double, 6> rigid_body_mass_diagonal(const Network& body) {
  std::array<double, 6> m{};
  m[0] = body.mass;
  m[1] = body.mass;
  m[2] = body.mass;
  for (int i = 3; i < 6; ++i)
    m[static_cast<std::size_t>(i)] = body.inertia[static_cast<std::size_t>(i)];
  return m;
}

Tddd velocity_unit_dof(int dof, const Tddd& x, const Tddd& ref, double rotation_sign = 1.0) {
  if (dof < 0 || dof > 5)
    return {0.0, 0.0, 0.0};
  if (dof <= 2) {
    Tddd v{0.0, 0.0, 0.0};
    v[dof] = 1.0;
    return v;
  }
  Tddd omega{0.0, 0.0, 0.0};
  omega[dof - 3] = 1.0;
  return rotation_sign * Cross(omega, x - ref);
}

double get_phi_on_face(const networkPoint* p, const networkFace* f) {
  if (!p)
    return 0.0;
  if (const auto* d = p->findActiveBieDofOrDefault(f))
    return d->phi;
  return std::get<0>(p->phiphin);
}

std::array<Complex, 6> integrate_complex_pressure_force(const Network& water, const std::unordered_set<networkFace*>& float_faces, const Tddd& ref_point, double rho, double omega, double moment_sign = 1.0) {
  const Complex I(0.0, 1.0);
  const Complex coef = I * omega * rho;

  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});

  // Degree-2 exact quadrature for products of linear fields on triangle.
  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};

  for (auto* f : float_faces) {
    auto [p0, p1, p2] = f->getPoints();
    const double phi0 = get_phi_on_face(p0, f);
    const double phi1 = get_phi_on_face(p1, f);
    const double phi2 = get_phi_on_face(p2, f);

    const Complex p_hat0 = coef * phi0;
    const Complex p_hat1 = coef * phi1;
    const Complex p_hat2 = coef * phi2;

    const auto& x0 = p0->X;
    const auto& x1 = p1->X;
    const auto& x2 = p2->X;
    const auto& n = f->normal;

    for (const auto& l : bary) {
      const double l0 = l[0], l1 = l[1], l2 = l[2];
      const Complex p_hat = l0 * p_hat0 + l1 * p_hat1 + l2 * p_hat2;
      const Tddd xq = l0 * x0 + l1 * x1 + l2 * x2;
      const std::array<Complex, 3> fq = {p_hat * n[0], p_hat * n[1], p_hat * n[2]};

      out[0] += w * fq[0] * f->area;
      out[1] += w * fq[1] * f->area;
      out[2] += w * fq[2] * f->area;

      const Tddd rq = xq - ref_point;
      const std::array<Complex, 3> tq = {
          rq[1] * fq[2] - rq[2] * fq[1],
          rq[2] * fq[0] - rq[0] * fq[2],
          rq[0] * fq[1] - rq[1] * fq[0],
      };
      out[3] += moment_sign * w * tq[0] * f->area;
      out[4] += moment_sign * w * tq[1] * f->area;
      out[5] += moment_sign * w * tq[2] * f->area;
    }
  }

  (void)water;
  return out;
}

enum class RadiationPostprocessVariant {
  LinearGeometryVertexField,
  LinearGeometryQuadraticField,
  StraightQuadraticGeometry6Node,
  TrueQuadraticGeometry6Node,
};

const char* radiation_postprocess_variant_name(RadiationPostprocessVariant variant) {
  switch (variant) {
  case RadiationPostprocessVariant::LinearGeometryVertexField:
    return "linear_geometry_vertex_field";
  case RadiationPostprocessVariant::LinearGeometryQuadraticField:
    return "linear_geometry_quadratic_field";
  case RadiationPostprocessVariant::StraightQuadraticGeometry6Node:
    return "straight_quadratic_geometry_6node";
  case RadiationPostprocessVariant::TrueQuadraticGeometry6Node:
    return "true_quadratic_geometry_6node";
  }
  return "unknown";
}

T6Tddd face_x6_for_postprocess(const networkFace* f, bool straight_midpoints) {
  auto [p0, p1, p2] = f->getPoints();
  auto [l0, l1, l2] = f->getLines();
  auto midpoint = [&](const networkLine* l) -> Tddd {
    if (straight_midpoints) {
      auto [pa, pb] = l->getPoints();
      if (pa && pb)
        return 0.5 * (pa->X + pb->X);
    }
    return l->X_mid;
  };
  return {p0->X, p1->X, p2->X, midpoint(l0), midpoint(l1), midpoint(l2)};
}

std::array<Complex, 6> integrate_pressure_force_postprocess_variant(const bem_frequency_domain::LinearSolution& sol,
                                                                    const std::unordered_set<networkFace*>& faces,
                                                                    const Tddd& ref_point,
                                                                    double rho,
                                                                    double moment_sign,
                                                                    RadiationPostprocessVariant variant) {
  const Complex I(0.0, 1.0);
  const Complex coef = I * sol.omega * rho;
  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});

  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary3 = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};
  constexpr std::array<bool, 3> all_true{true, true, true};

  auto add_force = [&](const Tddd& xq, const std::array<Complex, 3>& df) {
    out[0] += df[0];
    out[1] += df[1];
    out[2] += df[2];
    const Tddd r = xq - ref_point;
    out[3] += moment_sign * (r[1] * df[2] - r[2] * df[1]);
    out[4] += moment_sign * (r[2] * df[0] - r[0] * df[2]);
    out[5] += moment_sign * (r[0] * df[1] - r[1] * df[0]);
  };

  for (auto* f : faces) {
    const auto it = sol.face_field.find(f);
    if (it == sol.face_field.end())
      continue;
    const auto& field = it->second;
    auto [p0, p1, p2] = f->getPoints();
    const Tddd x0 = p0->X;
    const Tddd x1 = p1->X;
    const Tddd x2 = p2->X;

    const bool use_quadratic_field = bem_frequency_domain::face_field_is_quadratic(f, field);
    if ((variant == RadiationPostprocessVariant::StraightQuadraticGeometry6Node ||
         variant == RadiationPostprocessVariant::TrueQuadraticGeometry6Node) &&
        use_quadratic_field) {
      const auto p6 = bem_frequency_domain::phi6(field);
      const auto X6 = face_x6_for_postprocess(f, variant == RadiationPostprocessVariant::StraightQuadraticGeometry6Node);
      for (const auto& [x0q, x1q, w0w1] : __GWGW10__Tuple) {
        const auto bary = ModTriShape<3>(x0q, x1q);
        const double b0 = bary[0], b1 = bary[1];
        const auto N6 = f->trueQuadN6(b0, b1);
        const auto N6_geo = TriShape<6>(b0, b1, all_true);
        const auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
        const auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
        const Complex p_hat = coef * bem_frequency_domain::dot_shape6(N6, p6);
        const Tddd xq = Dot(N6_geo, X6);
        const Tddd area_vec = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
        const double weight = w0w1 * (1.0 - x0q);
        add_force(xq, {p_hat * area_vec[0] * weight,
                       p_hat * area_vec[1] * weight,
                       p_hat * area_vec[2] * weight});
      }
      continue;
    }

    const auto& n = f->normal;
    for (const auto& l : bary3) {
      const double l0 = l[0], l1 = l[1], l2 = l[2];
      Complex phi_q = l0 * field.phi[0] + l1 * field.phi[1] + l2 * field.phi[2];
      if (variant == RadiationPostprocessVariant::LinearGeometryQuadraticField && use_quadratic_field) {
        const auto N6 = f->trueQuadN6(l0, l1);
        phi_q = bem_frequency_domain::dot_shape6(N6, bem_frequency_domain::phi6(field));
      }
      const Complex p_hat = coef * phi_q;
      const Tddd xq = l0 * x0 + l1 * x1 + l2 * x2;
      add_force(xq, {w * p_hat * n[0] * f->area,
                     w * p_hat * n[1] * f->area,
                     w * p_hat * n[2] * f->area});
    }
  }

  return out;
}

void append_radiation_postprocess_probe(const std::filesystem::path& outdir,
                                        double omega,
                                        int dof_col,
                                        const char* solver,
                                        const bem_frequency_domain::LinearSolution& sol,
                                        const std::unordered_set<networkFace*>& float_faces,
                                        const Tddd& com,
                                        double rho,
                                        double moment_sign) {
  const std::vector<std::pair<std::string, Tddd>> refs = {{"COM", com}, {"origin", Tddd{0.0, 0.0, 0.0}}};
  const std::array<RadiationPostprocessVariant, 4> variants = {
      RadiationPostprocessVariant::LinearGeometryVertexField,
      RadiationPostprocessVariant::LinearGeometryQuadraticField,
      RadiationPostprocessVariant::StraightQuadraticGeometry6Node,
      RadiationPostprocessVariant::TrueQuadraticGeometry6Node,
  };

  struct ProbeValue {
    std::string variant;
    std::string ref_name;
    Tddd ref{};
    std::array<Complex, 6> z{};
  };
  std::vector<ProbeValue> values;
  values.reserve(refs.size() * variants.size());

  const auto probe_path = outdir / "radiation_postprocess_probe.csv";
  const bool probe_exists = std::filesystem::exists(probe_path) && std::filesystem::file_size(probe_path) > 0;
  std::ofstream pfs(probe_path, std::ios::app);
  if (!probe_exists)
    pfs << "omega,dof_col,dof_row,solver,variant,ref_name,ref_x,ref_y,ref_z,Z_re,Z_im,A,B\n";

  for (const auto& [ref_name, ref] : refs) {
    for (const auto variant : variants) {
      ProbeValue value;
      value.variant = radiation_postprocess_variant_name(variant);
      value.ref_name = ref_name;
      value.ref = ref;
      value.z = integrate_pressure_force_postprocess_variant(sol, float_faces, ref, rho, moment_sign, variant);
      values.push_back(value);
      for (int dof_row = 0; dof_row < 6; ++dof_row) {
        const Complex z = value.z[static_cast<std::size_t>(dof_row)];
        pfs << std::scientific << std::setprecision(12)
            << omega << "," << dof_col << "," << dof_row << "," << solver << ","
            << value.variant << "," << value.ref_name << ","
            << ref[0] << "," << ref[1] << "," << ref[2] << ","
            << z.real() << "," << z.imag() << ","
            << z.imag() / omega << "," << -z.real() << "\n";
      }
    }
  }

  const auto moment_path = outdir / "radiation_moment_reference_probe.csv";
  const bool moment_exists = std::filesystem::exists(moment_path) && std::filesystem::file_size(moment_path) > 0;
  std::ofstream mfs(moment_path, std::ios::app);
  if (!moment_exists)
    mfs << "omega,dof_col,dof_row,solver,variant,ref_a,ref_b,A_a,A_b,A_diff,B_a,B_b,B_diff\n";
  for (const auto& va : values) {
    if (va.ref_name != "COM" || va.variant != "true_quadratic_geometry_6node")
      continue;
    for (const auto& vb : values) {
      if (vb.ref_name != "origin" || vb.variant != va.variant)
        continue;
      for (int dof_row = 0; dof_row < 6; ++dof_row) {
        const Complex za = va.z[static_cast<std::size_t>(dof_row)];
        const Complex zb = vb.z[static_cast<std::size_t>(dof_row)];
        const double Aa = za.imag() / omega, Ab = zb.imag() / omega;
        const double Ba = -za.real(), Bb = -zb.real();
        mfs << std::scientific << std::setprecision(12)
            << omega << "," << dof_col << "," << dof_row << "," << solver << ","
            << va.variant << ",COM,origin,"
            << Aa << "," << Ab << "," << (Ab - Aa) << ","
            << Ba << "," << Bb << "," << (Bb - Ba) << "\n";
      }
    }
  }

  const auto diff_path = outdir / "radiation_quadratic_vs_linear_integral.csv";
  const bool diff_exists = std::filesystem::exists(diff_path) && std::filesystem::file_size(diff_path) > 0;
  std::ofstream dfs(diff_path, std::ios::app);
  if (!diff_exists)
    dfs << "omega,dof_col,dof_row,solver,ref_name,variant,baseline,A_variant,A_baseline,A_diff,A_rel_diff,B_variant,B_baseline,B_diff,B_rel_diff\n";
  for (const auto& base : values) {
    if (base.variant != "linear_geometry_vertex_field")
      continue;
    for (const auto& v : values) {
      if (v.ref_name != base.ref_name || v.variant == base.variant)
        continue;
      for (int dof_row = 0; dof_row < 6; ++dof_row) {
        const Complex zv = v.z[static_cast<std::size_t>(dof_row)];
        const Complex zb = base.z[static_cast<std::size_t>(dof_row)];
        const double Av = zv.imag() / omega, Ab = zb.imag() / omega;
        const double Bv = -zv.real(), Bb = -zb.real();
        dfs << std::scientific << std::setprecision(12)
            << omega << "," << dof_col << "," << dof_row << "," << solver << ","
            << v.ref_name << "," << v.variant << "," << base.variant << ","
            << Av << "," << Ab << "," << (Av - Ab) << "," << relative_pair_error(Av, Ab) << ","
            << Bv << "," << Bb << "," << (Bv - Bb) << "," << relative_pair_error(Bv, Bb) << "\n";
      }
    }
  }
}

double rigid_mode_normal_component(int dof, const Tddd& x, const Tddd& normal, const Tddd& ref_point, double moment_sign);

std::string pitch_surface_group(const networkFace* f) {
  if (!f)
    return "unknown";
  const double z = f->centroid[2];
  const double nz = std::abs(f->normal[2]);
  if (z > -1.5)
    return "waterline_near";
  if (z < -17.0 && nz > 0.45)
    return "heave_plate";
  if (z < -12.0 && nz > 0.45)
    return "bottom_plate";
  if (std::abs(f->normal[2]) < 0.35)
    return "column";
  return "other_float";
}

std::string wave_surface_group(const networkFace* f) {
  if (!f)
    return "unknown";
  const double z = f->centroid[2];
  const double nz = std::abs(f->normal[2]);
  if (z > -1.5)
    return "waterline_near";
  if (z < -17.0 && nz > 0.45)
    return "heave_plate";
  if (z < -8.0 && nz > 0.45)
    return "lower_horizontal";
  if (z < -8.0)
    return "pontoon_brace";
  if (nz < 0.35)
    return "column_side";
  return "other_float";
}

void append_radiation_rotation_mode_probe(const std::filesystem::path& path,
                                          double omega,
                                          int dof_col,
                                          const char* solver,
                                          const std::unordered_set<networkFace*>& float_faces,
                                          const Tddd& com,
                                          double rotation_sign,
                                          double moment_sign) {
  const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
  std::ofstream fs(path, std::ios::app);
  if (!exists)
    fs << "omega,dof_col,solver,ref_name,ref_x,ref_y,ref_z,face_index,group,area,"
          "q_l2,m_l2,q_minus_m_l2,q_plus_m_l2,max_abs_q_minus_m,max_abs_q_plus_m,q_area,m_area\n";

  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};
  const std::vector<std::pair<std::string, Tddd>> refs = {{"COM", com}, {"origin", Tddd{0.0, 0.0, 0.0}}};
  for (const auto& [ref_name, ref] : refs) {
    for (auto* f : float_faces) {
      auto [p0, p1, p2] = f->getPoints();
      long double q2 = 0.0L, m2 = 0.0L, diff2 = 0.0L, sum2 = 0.0L;
      long double q_area = 0.0L, m_area = 0.0L;
      double max_diff = 0.0, max_sum = 0.0;
      for (const auto& l : bary) {
        const Tddd xq = l[0] * p0->X + l[1] * p1->X + l[2] * p2->X;
        const double q = Dot(velocity_unit_dof(dof_col, xq, ref, rotation_sign), f->normal);
        const double m = rigid_mode_normal_component(dof_col, xq, f->normal, ref, moment_sign);
        const double wt = w * f->area;
        q2 += static_cast<long double>(wt) * q * q;
        m2 += static_cast<long double>(wt) * m * m;
        diff2 += static_cast<long double>(wt) * (q - m) * (q - m);
        sum2 += static_cast<long double>(wt) * (q + m) * (q + m);
        q_area += static_cast<long double>(wt) * q;
        m_area += static_cast<long double>(wt) * m;
        max_diff = std::max(max_diff, std::abs(q - m));
        max_sum = std::max(max_sum, std::abs(q + m));
      }
      fs << std::scientific << std::setprecision(12)
         << omega << "," << dof_col << "," << solver << ","
         << ref_name << "," << ref[0] << "," << ref[1] << "," << ref[2] << ","
         << f->index << "," << pitch_surface_group(f) << "," << f->area << ","
         << std::sqrt(static_cast<double>(q2)) << ","
         << std::sqrt(static_cast<double>(m2)) << ","
         << std::sqrt(static_cast<double>(diff2)) << ","
         << std::sqrt(static_cast<double>(sum2)) << ","
         << max_diff << "," << max_sum << ","
         << static_cast<double>(q_area) << "," << static_cast<double>(m_area) << "\n";
    }
  }
}

void append_radiation_pitch_bc_residual(const std::filesystem::path& path,
                                        double omega,
                                        int dof_col,
                                        const char* solver,
                                        const bem_frequency_domain::LinearSolution& sol,
                                        const std::unordered_set<networkFace*>& float_faces,
                                        const Tddd& ref_point,
                                        double rotation_sign) {
  struct Accum {
    std::size_t samples = 0;
    long double res2 = 0.0L;
    long double expected2 = 0.0L;
    long double phin2 = 0.0L;
    double max_abs_res = 0.0;
  };
  std::map<std::pair<std::string, std::string>, Accum> acc;
  auto add = [&](const std::string& group, const std::string& node_kind, const Complex& phin, double expected) {
    auto& a = acc[{group, node_kind}];
    const Complex expc{expected, 0.0};
    const Complex res = phin - expc;
    ++a.samples;
    a.res2 += std::norm(res);
    a.expected2 += expected * expected;
    a.phin2 += std::norm(phin);
    a.max_abs_res = std::max(a.max_abs_res, std::abs(res));
  };

  for (auto* f : float_faces) {
    const auto it = sol.face_field.find(f);
    if (it == sol.face_field.end())
      continue;
    const auto& field = it->second;
    const std::string group = pitch_surface_group(f);
    auto [p0, p1, p2] = f->getPoints();
    const std::array<networkPoint*, 3> pts{p0, p1, p2};
    for (std::size_t i = 0; i < pts.size(); ++i) {
      const double expected = Dot(velocity_unit_dof(dof_col, pts[i]->X, ref_point, rotation_sign), f->normal);
      add(group, "vertex", field.phin[i], expected);
    }
    if (bem_frequency_domain::face_field_is_quadratic(f, field)) {
      auto [l0, l1, l2] = f->getLines();
      const std::array<networkLine*, 3> lines{l0, l1, l2};
      for (std::size_t i = 0; i < lines.size(); ++i) {
        const double expected = Dot(velocity_unit_dof(dof_col, lines[i]->X_mid, ref_point, rotation_sign), f->normal);
        add(group, "midpoint", field.phin_mid[i], expected);
      }
    }
  }

  const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
  std::ofstream fs(path, std::ios::app);
  if (!exists)
    fs << "omega,dof_col,solver,group,node_kind,samples,res_l2,max_abs_res,expected_l2,phin_l2,ref_x,ref_y,ref_z,rotation_sign\n";
  for (const auto& [key, a] : acc) {
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver << ","
       << key.first << "," << key.second << "," << a.samples << ","
       << std::sqrt(static_cast<double>(a.res2)) << ","
       << a.max_abs_res << ","
       << std::sqrt(static_cast<double>(a.expected2)) << ","
       << std::sqrt(static_cast<double>(a.phin2)) << ","
       << ref_point[0] << "," << ref_point[1] << "," << ref_point[2] << ","
       << rotation_sign << "\n";
  }
}

std::string bem_entity_type_name(const BEM_DOF_Base* node) {
  if (!node)
    return "null";
  if (dynamic_cast<const networkPoint*>(node))
    return "point";
  if (dynamic_cast<const networkLine*>(node))
    return "line";
  return "unknown";
}

int bem_entity_id(const BEM_DOF_Base* node) {
  if (const auto* p = dynamic_cast<const networkPoint*>(node))
    return p->index;
  if (const auto* l = dynamic_cast<const networkLine*>(node))
    return l->midpoint_index;
  return -1;
}

std::string bie_key_kind_name(const BEM_DOF_Base* node, const networkFace* f) {
  if (!node)
    return "null";
  if (isDirichletBieDofKey(node, const_cast<networkFace*>(f)))
    return "dirichlet_type";
  if (isNeumannBieDofKey(node, const_cast<networkFace*>(f)))
    return "neumann_type";
  return "inactive_or_unknown";
}

std::string radiation_id_surface_name(const BEM_DOF_Base* node,
                                      const networkFace* f,
                                      const FaceSets& face_sets) {
  if (f)
    return face_kind_name(f, face_sets);
  if (!node)
    return "none";
  std::array<std::size_t, 3> counts{0, 0, 0}; // float, free, wall
  for (auto* adj : node->getBoundaryFaces()) {
    if (face_sets.float_surface.contains(adj))
      ++counts[0];
    else if (face_sets.free_surface.contains(adj))
      ++counts[1];
    else
      ++counts[2];
  }
  if (counts[0] >= counts[1] && counts[0] >= counts[2] && counts[0] > 0)
    return "float_default";
  if (counts[1] >= counts[0] && counts[1] >= counts[2] && counts[1] > 0)
    return "free_default";
  if (counts[2] > 0)
    return "wall_default";
  return "none";
}

std::string radiation_unknown_kind_name(const BEM_DOF_Base* node,
                                        const networkFace* f,
                                        const std::unordered_set<networkFace*>& robin_faces,
                                        const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin) {
  if (!node)
    return "null";
  bool robin = false;
  if (f) {
    robin = robin_faces.contains(const_cast<networkFace*>(f));
  } else {
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    robin = it != node_is_robin.end() && it->second;
  }
  if (robin)
    return isDirichletBieDofKey(node, const_cast<networkFace*>(f)) ? "robin_phin" : "robin_phi";
  if (isDirichletBieDofKey(node, const_cast<networkFace*>(f)))
    return "dirichlet_phin";
  if (isNeumannBieDofKey(node, const_cast<networkFace*>(f)))
    return "neumann_phi";
  return "inactive_or_unknown";
}

void append_radiation_fmm_matvec_compare(
    const std::filesystem::path& outdir,
    double omega,
    int dof_col,
    const char* solver_name,
    const std::vector<Complex>& dense_A,
    const std::vector<Complex>& dense_rhs,
    const std::function<std::vector<Complex>(const std::vector<Complex>&)>& fmm_matvec,
    const std::vector<Complex>& fmm_rhs,
    const std::vector<Complex>& u,
    const std::vector<bem_frequency_domain::Id>& id_by_index,
    const std::unordered_set<networkFace*>& robin_faces,
    const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
    const FaceSets& face_sets) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || dense_A.size() != n * n || dense_rhs.size() != n || fmm_rhs.size() != n || u.size() != n)
    return;
  auto aidx = [&](std::size_t row, std::size_t col) { return row + col * n; };
  auto dense_matvec = [&](const std::vector<Complex>& x) {
    std::vector<Complex> y(n, Complex{0.0, 0.0});
    for (std::size_t col = 0; col < n; ++col) {
      const Complex xc = x[col];
      if (xc == Complex{0.0, 0.0})
        continue;
      for (std::size_t row = 0; row < n; ++row)
        y[row] += dense_A[aidx(row, col)] * xc;
    }
    return y;
  };

  struct RowDiff {
    std::size_t row = 0;
    Complex dense{0.0, 0.0};
    Complex fmm{0.0, 0.0};
    double diff_abs = 0.0;
    double rel = 0.0;
  };
  auto write_probe = [&](const std::string& probe_name,
                         const std::vector<Complex>& dense_y,
                         const std::vector<Complex>& fmm_y) {
    if (dense_y.size() != n || fmm_y.size() != n)
      return;
    std::vector<RowDiff> rows;
    rows.reserve(n);
    long double dense_norm2 = 0.0L;
    long double fmm_norm2 = 0.0L;
    long double diff_norm2 = 0.0L;
    std::size_t worst_row = 0;
    double worst_abs = 0.0;
    double worst_rel = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const Complex d = dense_y[i];
      const Complex f = fmm_y[i];
      const double diff_abs = std::abs(d - f);
      const double denom = std::max({std::abs(d), std::abs(f), 1e-300});
      const double rel = diff_abs / denom;
      dense_norm2 += static_cast<long double>(std::norm(d));
      fmm_norm2 += static_cast<long double>(std::norm(f));
      diff_norm2 += static_cast<long double>(diff_abs) * static_cast<long double>(diff_abs);
      if (diff_abs > worst_abs) {
        worst_abs = diff_abs;
        worst_rel = rel;
        worst_row = i;
      }
      rows.push_back(RowDiff{i, d, f, diff_abs, rel});
    }

    const double dense_norm = std::sqrt(static_cast<double>(dense_norm2));
    const double fmm_norm = std::sqrt(static_cast<double>(fmm_norm2));
    const double diff_norm = std::sqrt(static_cast<double>(diff_norm2));
    const double rel_l2 = diff_norm / std::max({dense_norm, fmm_norm, 1e-300});
    const auto* worst_node = std::get<0>(id_by_index[worst_row]);
    const auto* worst_face = std::get<1>(id_by_index[worst_row]);
    const Tddd worst_x = worst_node ? worst_node->getPosition() : Tddd{0.0, 0.0, 0.0};

    const auto summary_path = outdir / "radiation_fmm_matvec_compare_summary.csv";
    const bool summary_exists = std::filesystem::exists(summary_path) && std::filesystem::file_size(summary_path) > 0;
    std::ofstream summary(summary_path, std::ios::app);
    if (!summary_exists) {
      summary << "omega,dof_col,solver,probe,n,dense_norm,fmm_norm,diff_l2,rel_l2,"
                 "max_abs_diff,max_rel_diff,worst_row,worst_entity,worst_entity_id,worst_face_index,"
                 "worst_surface,worst_key_kind,worst_unknown,is_bc_interface,worst_x,worst_y,worst_z\n";
    }
    summary << std::scientific << std::setprecision(12)
            << omega << "," << dof_col << "," << solver_name << "," << probe_name << ","
            << n << "," << dense_norm << "," << fmm_norm << "," << diff_norm << "," << rel_l2 << ","
            << worst_abs << "," << worst_rel << "," << worst_row << ","
            << bem_entity_type_name(worst_node) << "," << bem_entity_id(worst_node) << ","
            << (worst_face ? worst_face->index : -1) << ","
            << radiation_id_surface_name(worst_node, worst_face, face_sets) << ","
            << bie_key_kind_name(worst_node, worst_face) << ","
            << radiation_unknown_kind_name(worst_node, worst_face, robin_faces, node_is_robin) << ","
            << (worst_node && worst_node->BCInterface ? 1 : 0) << ","
            << worst_x[0] << "," << worst_x[1] << "," << worst_x[2] << "\n";

    std::ranges::sort(rows, [](const RowDiff& a, const RowDiff& b) {
      return a.diff_abs > b.diff_abs;
    });
    const auto rows_path = outdir / "radiation_fmm_matvec_compare_rows.csv";
    const bool rows_exists = std::filesystem::exists(rows_path) && std::filesystem::file_size(rows_path) > 0;
    std::ofstream out(rows_path, std::ios::app);
    if (!rows_exists) {
      out << "omega,dof_col,solver,probe,rank,row,entity_type,entity_id,entity_ptr,face_index,"
             "surface,key_kind,unknown,is_bc_interface,is_robin_node,x,y,z,"
             "dense_re,dense_im,fmm_re,fmm_im,diff_abs,rel_diff\n";
    }
    const std::size_t top_k = std::min<std::size_t>(200, rows.size());
    for (std::size_t rank = 0; rank < top_k; ++rank) {
      const auto& r = rows[rank];
      const auto* node = std::get<0>(id_by_index[r.row]);
      const auto* face = std::get<1>(id_by_index[r.row]);
      const Tddd x = node ? node->getPosition() : Tddd{0.0, 0.0, 0.0};
      const auto it_robin = node ? node_is_robin.find(const_cast<BEM_DOF_Base*>(node)) : node_is_robin.end();
      const bool node_robin = it_robin != node_is_robin.end() && it_robin->second;
      out << std::scientific << std::setprecision(12)
          << omega << "," << dof_col << "," << solver_name << "," << probe_name << ","
          << rank << "," << r.row << ","
          << bem_entity_type_name(node) << "," << bem_entity_id(node) << ","
          << reinterpret_cast<std::uintptr_t>(node) << ","
          << (face ? face->index : -1) << ","
          << radiation_id_surface_name(node, face, face_sets) << ","
          << bie_key_kind_name(node, face) << ","
          << radiation_unknown_kind_name(node, face, robin_faces, node_is_robin) << ","
          << (node && node->BCInterface ? 1 : 0) << ","
          << (node_robin ? 1 : 0) << ","
          << x[0] << "," << x[1] << "," << x[2] << ","
          << r.dense.real() << "," << r.dense.imag() << ","
          << r.fmm.real() << "," << r.fmm.imag() << ","
          << r.diff_abs << "," << r.rel << "\n";
    }
  };

  const auto dense_solution = dense_matvec(u);
  const auto fmm_solution = fmm_matvec(u);
  write_probe("solution_matvec", dense_solution, fmm_solution);
  write_probe("rhs", dense_rhs, fmm_rhs);
}

double expected_radiation_phin_on_key(const BEM_DOF_Base& node,
                                      const networkFace* f,
                                      const std::unordered_set<networkFace*>& float_faces,
                                      int dof_col,
                                      const Tddd& ref_point,
                                      double rotation_sign) {
  bool on_float = false;
  Tddd normal = {0.0, 0.0, 0.0};
  if (f) {
    on_float = float_faces.contains(const_cast<networkFace*>(f));
    normal = f->normal;
  } else {
    for (auto* adj : node.getBoundaryFaces()) {
      if (float_faces.contains(adj)) {
        on_float = true;
        normal += adj->normal;
      }
    }
    const double nm = Norm(normal);
    if (nm > 1e-15)
      normal /= nm;
  }
  if (!on_float)
    return 0.0;
  return Dot(velocity_unit_dof(dof_col, node.getPosition(), ref_point, rotation_sign), normal);
}

struct SolutionLookupProbe {
  std::string status = "miss";
  int index = -1;
  Complex value{0.0, 0.0};
  bool key_is_dirichlet_type = false;
  bool key_is_neumann_type = false;
};

std::unordered_map<bem_frequency_domain::Id, std::size_t, bem_frequency_domain::IdHash, bem_frequency_domain::IdEq>
build_solution_id_index(const bem_frequency_domain::Solution& sol) {
  std::unordered_map<bem_frequency_domain::Id, std::size_t, bem_frequency_domain::IdHash, bem_frequency_domain::IdEq> index;
  index.reserve(sol.n);
  for (std::size_t i = 0; i < sol.n; ++i)
    index.emplace(sol.id_by_index[i], i);
  return index;
}

SolutionLookupProbe lookup_solution_value_probe(
    const std::unordered_map<bem_frequency_domain::Id, std::size_t, bem_frequency_domain::IdHash, bem_frequency_domain::IdEq>& index,
    const bem_frequency_domain::Solution& raw,
    const BEM_DOF_Base* node,
    const networkFace* f,
    const std::vector<Complex>& values) {
  SolutionLookupProbe out;
  if (!node)
    return out;
  bem_frequency_domain::Id key{const_cast<BEM_DOF_Base*>(node), const_cast<networkFace*>(f)};
  auto it = index.find(key);
  if (it != index.end()) {
    out.status = "exact_face_key";
    out.index = static_cast<int>(it->second);
    out.value = values[it->second];
    out.key_is_dirichlet_type = isDirichletBieDofKey(node, const_cast<networkFace*>(f));
    out.key_is_neumann_type = isNeumannBieDofKey(node, const_cast<networkFace*>(f));
    return out;
  }
  key = {const_cast<BEM_DOF_Base*>(node), nullptr};
  it = index.find(key);
  if (it != index.end()) {
    out.status = "default_null_face";
    out.index = static_cast<int>(it->second);
    out.value = values[it->second];
    out.key_is_dirichlet_type = isDirichletBieDofKey(node, nullptr);
    out.key_is_neumann_type = isNeumannBieDofKey(node, nullptr);
    return out;
  }
  return out;
}

void append_radiation_solution_id_probe(const std::filesystem::path& path,
                                        double omega,
                                        int dof_col,
                                        const char* solver,
                                        const bem_frequency_domain::Solution& raw,
                                        const std::unordered_set<networkFace*>& float_faces,
                                        const Tddd& ref_point,
                                        double rotation_sign) {
  const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
  std::ofstream fs(path, std::ios::app);
  if (!exists)
    fs << "omega,dof_col,solver,solution_index,entity_type,entity_id,entity_ptr,face_index,group,"
          "key_kind,is_bc_interface,position_x,position_y,position_z,expected_phin,"
          "solution_phin_re,solution_phin_im,res_abs,solution_phi_re,solution_phi_im,u_re,u_im\n";
  for (std::size_t i = 0; i < raw.n; ++i) {
    auto [node, f] = raw.id_by_index[i];
    if (!node)
      continue;
    const double expected = expected_radiation_phin_on_key(*node, f, float_faces, dof_col, ref_point, rotation_sign);
    const Complex res = raw.phin[i] - Complex{expected, 0.0};
    const Tddd& x = node->getPosition();
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver << "," << i << ","
       << bem_entity_type_name(node) << "," << bem_entity_id(node) << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << (f ? f->index : -1) << "," << (f ? pitch_surface_group(f) : "default_id") << ","
       << bie_key_kind_name(node, f) << "," << (node->BCInterface ? 1 : 0) << ","
       << x[0] << "," << x[1] << "," << x[2] << ","
       << expected << ","
       << raw.phin[i].real() << "," << raw.phin[i].imag() << "," << std::abs(res) << ","
       << raw.phi[i].real() << "," << raw.phi[i].imag() << ","
       << raw.u[i].real() << "," << raw.u[i].imag() << "\n";
  }
}

void append_radiation_facefield_capture_probe(const std::filesystem::path& path,
                                              double omega,
                                              int dof_col,
                                              const char* solver,
                                              const bem_frequency_domain::Solution& raw,
                                              const bem_frequency_domain::LinearSolution& captured,
                                              const std::unordered_set<networkFace*>& float_faces,
                                              const Tddd& ref_point,
                                              double rotation_sign) {
  const auto index = build_solution_id_index(raw);
  const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
  std::ofstream fs(path, std::ios::app);
  if (!exists)
    fs << "omega,dof_col,solver,face_index,group,node_kind,local_node,entity_type,entity_id,entity_ptr,"
          "lookup_status,solution_index,key_is_dirichlet_type,key_is_neumann_type,expected_phin,"
          "solution_phin_re,solution_phin_im,facefield_phin_re,facefield_phin_im,"
          "solution_minus_facefield_abs,solution_res_abs,facefield_res_abs,"
          "solution_phi_re,solution_phi_im,facefield_phi_re,facefield_phi_im\n";

  auto write_sample = [&](networkFace* f,
                          const std::string& group,
                          const std::string& node_kind,
                          int local_node,
                          const BEM_DOF_Base* node,
                          const Complex& field_phi,
                          const Complex& field_phin) {
    const auto phin_lookup = lookup_solution_value_probe(index, raw, node, f, raw.phin);
    const auto phi_lookup = lookup_solution_value_probe(index, raw, node, f, raw.phi);
    const double expected = expected_radiation_phin_on_key(*node, f, float_faces, dof_col, ref_point, rotation_sign);
    const Complex expected_c{expected, 0.0};
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver << ","
       << (f ? f->index : -1) << "," << group << ","
       << node_kind << "," << local_node << ","
       << bem_entity_type_name(node) << "," << bem_entity_id(node) << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << phin_lookup.status << "," << phin_lookup.index << ","
       << (phin_lookup.key_is_dirichlet_type ? 1 : 0) << ","
       << (phin_lookup.key_is_neumann_type ? 1 : 0) << ","
       << expected << ","
       << phin_lookup.value.real() << "," << phin_lookup.value.imag() << ","
       << field_phin.real() << "," << field_phin.imag() << ","
       << std::abs(phin_lookup.value - field_phin) << ","
       << std::abs(phin_lookup.value - expected_c) << ","
       << std::abs(field_phin - expected_c) << ","
       << phi_lookup.value.real() << "," << phi_lookup.value.imag() << ","
       << field_phi.real() << "," << field_phi.imag() << "\n";
  };

  for (auto* f : float_faces) {
    const auto it = captured.face_field.find(f);
    if (it == captured.face_field.end())
      continue;
    const std::string group = pitch_surface_group(f);
    const auto& field = it->second;
    auto [p0, p1, p2] = f->getPoints();
    const std::array<networkPoint*, 3> pts{p0, p1, p2};
    for (std::size_t i = 0; i < pts.size(); ++i)
      write_sample(f, group, "vertex", static_cast<int>(i), pts[i], field.phi[i], field.phin[i]);
    if (bem_frequency_domain::face_field_is_quadratic(f, field)) {
      auto [l0, l1, l2] = f->getLines();
      const std::array<networkLine*, 3> lines{l0, l1, l2};
      for (std::size_t i = 0; i < lines.size(); ++i)
        write_sample(f, group, "midpoint", static_cast<int>(i), lines[i], field.phi_mid[i], field.phin_mid[i]);
    }
  }
}

Complex integrate_pitch_face_impedance_contribution(const bem_frequency_domain::LinearSolution& sol,
                                                    networkFace* f,
                                                    const Tddd& ref_point,
                                                    double rho,
                                                    double moment_sign) {
  const auto it = sol.face_field.find(f);
  if (it == sol.face_field.end())
    return Complex{0.0, 0.0};
  const Complex I{0.0, 1.0};
  const Complex coef = I * sol.omega * rho;
  const auto& field = it->second;
  Complex out{0.0, 0.0};
  auto add_force = [&](const Tddd& xq, const std::array<Complex, 3>& df) {
    const Tddd r = xq - ref_point;
    out += moment_sign * (r[2] * df[0] - r[0] * df[2]);
  };
  if (bem_frequency_domain::face_field_is_quadratic(f, field)) {
    constexpr std::array<bool, 3> all_true{true, true, true};
    const auto p6 = bem_frequency_domain::phi6(field);
    const auto X6 = bem_frequency_domain::face_x6(f);
    for (const auto& [x0q, x1q, w0w1] : __GWGW10__Tuple) {
      const auto bary = ModTriShape<3>(x0q, x1q);
      const double b0 = bary[0], b1 = bary[1];
      const auto N6 = f->trueQuadN6(b0, b1);
      const auto N6_geo = TriShape<6>(b0, b1, all_true);
      const auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
      const auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
      const Complex p_hat = coef * bem_frequency_domain::dot_shape6(N6, p6);
      const Tddd xq = Dot(N6_geo, X6);
      const Tddd area_vec = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
      const double weight = w0w1 * (1.0 - x0q);
      add_force(xq, {p_hat * area_vec[0] * weight, p_hat * area_vec[1] * weight, p_hat * area_vec[2] * weight});
    }
    return out;
  }

  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary3 = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};
  auto [p0, p1, p2] = f->getPoints();
  for (const auto& l : bary3) {
    const Complex phi_q = l[0] * field.phi[0] + l[1] * field.phi[1] + l[2] * field.phi[2];
    const Complex p_hat = coef * phi_q;
    const Tddd xq = l[0] * p0->X + l[1] * p1->X + l[2] * p2->X;
    add_force(xq, {w * p_hat * f->normal[0] * f->area,
                   w * p_hat * f->normal[1] * f->area,
                   w * p_hat * f->normal[2] * f->area});
  }
  return out;
}

void append_radiation_pitch_surface_contribution(const std::filesystem::path& path,
                                                 double omega,
                                                 int dof_col,
                                                 const char* solver,
                                                 const bem_frequency_domain::LinearSolution& sol,
                                                 const std::unordered_set<networkFace*>& float_faces,
                                                 const Tddd& ref_point,
                                                 double rho,
                                                 double moment_sign) {
  struct Sum {
    std::size_t faces = 0;
    double area = 0.0;
    Complex z{0.0, 0.0};
  };
  std::map<std::string, Sum> sums;
  const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
  std::ofstream fs(path, std::ios::app);
  if (!exists)
    fs << "omega,dof_col,solver,row_kind,group,face_index,area,cx,cy,cz,nx,ny,nz,Z_re,Z_im,A44,B44,missing_geometry_note\n";
  for (auto* f : float_faces) {
    const std::string group = pitch_surface_group(f);
    const Complex z = integrate_pitch_face_impedance_contribution(sol, f, ref_point, rho, moment_sign);
    auto& s = sums[group];
    ++s.faces;
    s.area += f->area;
    s.z += z;
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver << ",face,"
       << group << "," << f->index << "," << f->area << ","
       << f->centroid[0] << "," << f->centroid[1] << "," << f->centroid[2] << ","
       << f->normal[0] << "," << f->normal[1] << "," << f->normal[2] << ","
       << z.real() << "," << z.imag() << "," << z.imag() / omega << "," << -z.real()
       << ",black_braces_not_in_current_float_obj\n";
  }
  for (const auto& [group, s] : sums) {
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver << ",summary,"
       << group << ",-1," << s.area << ",0,0,0,0,0,0,"
       << s.z.real() << "," << s.z.imag() << "," << s.z.imag() / omega << "," << -s.z.real()
       << ",black_braces_not_in_current_float_obj\n";
	  }
	}

	bool radiation_coupling_probe_row_selected(int dof_row, int dof_col) {
	  const bool row_major = (dof_row == 0 || dof_row == 2 || dof_row == 4);
	  const bool col_major = (dof_col == 0 || dof_col == 2 || dof_col == 4);
	  return row_major && col_major;
	}

		const char* radiation_coupling_term_kind(int dof_row, int dof_col) {
		  if (dof_row == 4 && dof_col == 4)
		    return "pitch_diag";
		  if (dof_row == dof_col && (dof_row == 0 || dof_row == 2 || dof_row == 4))
		    return "low_frequency_damping_diag";
		  if ((dof_row == 4 && (dof_col == 0 || dof_col == 2)) ||
		      (dof_col == 4 && (dof_row == 0 || dof_row == 2)))
		    return "pitch_added_mass_coupling";
		  return "other_major_024";
		}

		std::array<Complex, 6> integrate_face_impedance_contribution(const bem_frequency_domain::LinearSolution& sol,
		                                                             networkFace* f,
		                                                             const Tddd& ref_point,
		                                                             double rho,
		                                                             double moment_sign) {
	  std::array<Complex, 6> out{};
	  out.fill(Complex{0.0, 0.0});
	  const auto it = sol.face_field.find(f);
	  if (it == sol.face_field.end())
	    return out;
	  const Complex I{0.0, 1.0};
	  const Complex coef = I * sol.omega * rho;
	  const auto& field = it->second;
	  auto add_force = [&](const Tddd& xq, const std::array<Complex, 3>& df) {
	    out[0] += df[0];
	    out[1] += df[1];
	    out[2] += df[2];
	    const Tddd r = xq - ref_point;
	    out[3] += moment_sign * (r[1] * df[2] - r[2] * df[1]);
	    out[4] += moment_sign * (r[2] * df[0] - r[0] * df[2]);
	    out[5] += moment_sign * (r[0] * df[1] - r[1] * df[0]);
	  };
	  if (bem_frequency_domain::face_field_is_quadratic(f, field)) {
	    constexpr std::array<bool, 3> all_true{true, true, true};
	    const auto p6 = bem_frequency_domain::phi6(field);
	    const auto X6 = bem_frequency_domain::face_x6(f);
	    for (const auto& [x0q, x1q, w0w1] : __GWGW10__Tuple) {
	      const auto bary = ModTriShape<3>(x0q, x1q);
	      const double b0 = bary[0], b1 = bary[1];
	      const auto N6 = f->trueQuadN6(b0, b1);
	      const auto N6_geo = TriShape<6>(b0, b1, all_true);
	      const auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
	      const auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
	      const Complex p_hat = coef * bem_frequency_domain::dot_shape6(N6, p6);
	      const Tddd xq = Dot(N6_geo, X6);
	      const Tddd area_vec = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
	      const double weight = w0w1 * (1.0 - x0q);
	      add_force(xq, {p_hat * area_vec[0] * weight,
	                     p_hat * area_vec[1] * weight,
	                     p_hat * area_vec[2] * weight});
	    }
	    return out;
	  }

	  constexpr double w = 1.0 / 3.0;
	  constexpr std::array<std::array<double, 3>, 3> bary3 = {{
	      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
	      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
	      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
	  }};
	  auto [p0, p1, p2] = f->getPoints();
	  for (const auto& l : bary3) {
	    const double l0 = l[0], l1 = l[1], l2 = l[2];
	    const Complex phi_q = l0 * field.phi[0] + l1 * field.phi[1] + l2 * field.phi[2];
	    const Complex p_hat = coef * phi_q;
	    const Tddd xq = l0 * p0->X + l1 * p1->X + l2 * p2->X;
	    add_force(xq, {w * p_hat * f->normal[0] * f->area,
	                   w * p_hat * f->normal[1] * f->area,
	                   w * p_hat * f->normal[2] * f->area});
	  }
		  return out;
		}

void append_wave_excitation_group_probe(const std::filesystem::path& outdir,
                                        double omega,
                                        const char* solver,
                                        const bem_frequency_domain::LinearSolution& incident,
                                        const bem_frequency_domain::LinearSolution& diffraction,
                                        const bem_frequency_domain::LinearSolution& total,
                                        const std::unordered_set<networkFace*>& float_faces,
                                        const Tddd& ref_point,
                                        double rho) {
  struct GroupSum {
    std::size_t faces = 0;
    double area = 0.0;
    std::array<Complex, 6> force{};
  };
  struct Component {
    const char* name = "";
    const bem_frequency_domain::LinearSolution* sol = nullptr;
    std::array<Complex, 6> all{};
  };
  std::array<Component, 3> components{{
      {"incident", &incident, bem_frequency_domain::integrate_linear_pressure_force(incident, float_faces, ref_point, rho)},
      {"diffraction", &diffraction, bem_frequency_domain::integrate_linear_pressure_force(diffraction, float_faces, ref_point, rho)},
      {"total", &total, bem_frequency_domain::integrate_linear_pressure_force(total, float_faces, ref_point, rho)},
  }};

  const auto group_path = outdir / "wave_excitation_face_group_contrib.csv";
  const bool group_exists = std::filesystem::exists(group_path) && std::filesystem::file_size(group_path) > 0;
  std::ofstream gfs(group_path, std::ios::app);
  if (!group_exists)
    gfs << "omega,component,solver,group,dof,faces,area,value_re,value_im,amp,phase,total_re,total_im,total_amp,amp_fraction,complex_fraction_re,complex_fraction_im\n";

  const auto face_path = outdir / "wave_excitation_face_contrib.csv";
  const bool face_exists = std::filesystem::exists(face_path) && std::filesystem::file_size(face_path) > 0;
  std::ofstream ffs(face_path, std::ios::app);
  if (!face_exists)
    ffs << "omega,component,solver,group,face_index,area,cx,cy,cz,nx,ny,nz,dof,value_re,value_im,amp,phase\n";

  for (const auto& comp : components) {
    std::map<std::string, GroupSum> sums;
    for (auto* f : float_faces) {
      if (!f || !comp.sol)
        continue;
      const std::string group = wave_surface_group(f);
      const auto z_face = integrate_face_impedance_contribution(*comp.sol, f, ref_point, rho, 1.0);
      auto& s = sums[group];
      ++s.faces;
      s.area += f->area;
      for (int dof = 0; dof < 6; ++dof)
        s.force[static_cast<std::size_t>(dof)] += z_face[static_cast<std::size_t>(dof)];

      for (int dof : {0, 2, 4}) {
        const Complex z = z_face[static_cast<std::size_t>(dof)];
        ffs << std::scientific << std::setprecision(12)
            << omega << "," << comp.name << "," << solver << ","
            << group << "," << f->index << "," << f->area << ","
            << f->centroid[0] << "," << f->centroid[1] << "," << f->centroid[2] << ","
            << f->normal[0] << "," << f->normal[1] << "," << f->normal[2] << ","
            << dof << "," << z.real() << "," << z.imag() << ","
            << std::abs(z) << "," << std::atan2(z.imag(), z.real()) << "\n";
      }
    }

    for (const auto& [group, s] : sums) {
      for (int dof = 0; dof < 6; ++dof) {
        const Complex z = s.force[static_cast<std::size_t>(dof)];
        const Complex z_total = comp.all[static_cast<std::size_t>(dof)];
        const double total_amp = std::abs(z_total);
        const double scale = std::max(total_amp, 1e-300);
        const Complex complex_fraction = total_amp > 1e-300 ? z / z_total : Complex{0.0, 0.0};
        gfs << std::scientific << std::setprecision(12)
            << omega << "," << comp.name << "," << solver << ","
            << group << "," << dof << "," << s.faces << "," << s.area << ","
            << z.real() << "," << z.imag() << ","
            << std::abs(z) << "," << std::atan2(z.imag(), z.real()) << ","
            << z_total.real() << "," << z_total.imag() << "," << total_amp << ","
            << std::abs(z) / scale << ","
            << complex_fraction.real() << "," << complex_fraction.imag() << "\n";
      }
    }
  }
}

		void append_radiation_coupling_probe(const std::filesystem::path& outdir,
		                                     double omega,
	                                     int dof_col,
	                                     const char* solver,
	                                     const bem_frequency_domain::LinearSolution& sol,
	                                     const std::unordered_set<networkFace*>& float_faces,
	                                     const Tddd& radiation_ref_point,
	                                     const Tddd& com,
	                                     double rho,
	                                     double moment_sign) {
	  const std::vector<std::pair<std::string, Tddd>> refs = {
	      {"solver_ref", radiation_ref_point},
	      {"COM", com},
	      {"origin", Tddd{0.0, 0.0, 0.0}},
	  };

	  const auto term_path = outdir / "radiation_coupling_terms.csv";
	  const bool term_exists = std::filesystem::exists(term_path) && std::filesystem::file_size(term_path) > 0;
	  std::ofstream tfs(term_path, std::ios::app);
	  if (!term_exists)
	    tfs << "omega,dof_row,dof_col,term_kind,solver,ref_name,ref_x,ref_y,ref_z,Z_re,Z_im,A,B\n";

	  const auto shift_path = outdir / "radiation_reference_shift_probe.csv";
	  const bool shift_exists = std::filesystem::exists(shift_path) && std::filesystem::file_size(shift_path) > 0;
	  std::ofstream rfs(shift_path, std::ios::app);
	  if (!shift_exists)
	    rfs << "omega,dof_row,dof_col,term_kind,solver,from_ref,to_ref,from_A,from_B,to_A,to_B,A_diff,B_diff,note\n";

	  std::map<std::string, std::array<Complex, 6>> z_by_ref;
	  for (const auto& [ref_name, ref] : refs) {
	    const auto z = bem_frequency_domain::integrate_linear_pressure_force(sol, float_faces, ref, rho, moment_sign);
	    z_by_ref[ref_name] = z;
	    for (int dof_row : {0, 2, 4}) {
	      if (!radiation_coupling_probe_row_selected(dof_row, dof_col))
	        continue;
	      const Complex zij = z[static_cast<std::size_t>(dof_row)];
	      tfs << std::scientific << std::setprecision(12)
	          << omega << "," << dof_row << "," << dof_col << ","
	          << radiation_coupling_term_kind(dof_row, dof_col) << ","
	          << solver << "," << ref_name << ","
	          << ref[0] << "," << ref[1] << "," << ref[2] << ","
	          << zij.real() << "," << zij.imag() << ","
	          << zij.imag() / omega << "," << -zij.real() << "\n";
	    }
	  }

	  const auto& z_from = z_by_ref["solver_ref"];
	  for (const std::string to_ref : {"COM", "origin"}) {
	    const auto it = z_by_ref.find(to_ref);
	    if (it == z_by_ref.end())
	      continue;
	    for (int dof_row : {0, 2, 4}) {
	      if (!radiation_coupling_probe_row_selected(dof_row, dof_col))
	        continue;
	      const Complex za = z_from[static_cast<std::size_t>(dof_row)];
	      const Complex zb = it->second[static_cast<std::size_t>(dof_row)];
	      const double Aa = za.imag() / omega, Ab = zb.imag() / omega;
	      const double Ba = -za.real(), Bb = -zb.real();
	      rfs << std::scientific << std::setprecision(12)
	          << omega << "," << dof_row << "," << dof_col << ","
	          << radiation_coupling_term_kind(dof_row, dof_col) << ","
	          << solver << ",solver_ref," << to_ref << ","
	          << Aa << "," << Ba << "," << Ab << "," << Bb << ","
	          << (Ab - Aa) << "," << (Bb - Ba)
	          << ",direct_reintegration_same_solution\n";
	    }
	  }

	  const auto contrib_path = outdir / "radiation_surface_contribution_by_term.csv";
	  const bool contrib_exists = std::filesystem::exists(contrib_path) && std::filesystem::file_size(contrib_path) > 0;
	  std::ofstream cfs(contrib_path, std::ios::app);
	  if (!contrib_exists)
	    cfs << "omega,dof_row,dof_col,term_kind,solver,ref_name,row_kind,group,face_index,area,cx,cy,cz,nx,ny,nz,Z_re,Z_im,A,B,missing_geometry_note\n";
	  for (const auto& [ref_name, ref] : refs) {
	    struct Sum {
	      std::size_t faces = 0;
	      double area = 0.0;
	      Complex z{0.0, 0.0};
	    };
	    std::map<std::pair<int, std::string>, Sum> sums;
	    for (auto* f : float_faces) {
	      const auto z_face = integrate_face_impedance_contribution(sol, f, ref, rho, moment_sign);
	      const std::string group = pitch_surface_group(f);
	      for (int dof_row : {0, 2, 4}) {
	        if (!radiation_coupling_probe_row_selected(dof_row, dof_col))
	          continue;
	        const Complex z = z_face[static_cast<std::size_t>(dof_row)];
	        auto& s = sums[{dof_row, group}];
	        ++s.faces;
	        s.area += f->area;
	        s.z += z;
	      }
	    }
	    for (const auto& [key, s] : sums) {
	      const int dof_row = key.first;
	      const std::string& group = key.second;
	      cfs << std::scientific << std::setprecision(12)
	          << omega << "," << dof_row << "," << dof_col << ","
	          << radiation_coupling_term_kind(dof_row, dof_col) << ","
	          << solver << "," << ref_name << ",summary,"
	          << group << ",-1," << s.area << ",0,0,0,0,0,0,"
	          << s.z.real() << "," << s.z.imag() << "," << s.z.imag() / omega << "," << -s.z.real()
	          << ",black_braces_not_in_current_float_obj\n";
	    }
	  }
	}

	void append_radiation_pitch_probe(const std::filesystem::path& outdir,
	                                  double omega,
	                                  int dof_col,
	                                  const char* solver,
                                  const bem_frequency_domain::Solution& raw,
                                  const bem_frequency_domain::LinearSolution& sol,
                                  const std::unordered_set<networkFace*>& float_faces,
                                  const Tddd& com,
                                  const Tddd& radiation_ref_point,
                                  double rho,
                                  double rotation_sign,
                                  double moment_sign) {
  append_radiation_solution_id_probe(outdir / "radiation_solution_id_probe.csv",
                                     omega,
                                     dof_col,
                                     solver,
                                     raw,
                                     float_faces,
                                     radiation_ref_point,
                                     rotation_sign);
  append_radiation_facefield_capture_probe(outdir / "radiation_facefield_capture_probe.csv",
                                           omega,
                                           dof_col,
                                           solver,
                                           raw,
                                           sol,
                                           float_faces,
                                           radiation_ref_point,
                                           rotation_sign);
  append_radiation_rotation_mode_probe(outdir / "radiation_rotation_mode_probe.csv",
                                       omega,
                                       dof_col,
                                       solver,
                                       float_faces,
                                       com,
                                       rotation_sign,
                                       moment_sign);
  append_radiation_pitch_bc_residual(outdir / "radiation_pitch_bc_residual.csv",
                                     omega,
                                     dof_col,
                                     solver,
                                     sol,
                                     float_faces,
                                     radiation_ref_point,
                                     rotation_sign);
  if (dof_col != 4)
    return;
  append_radiation_pitch_surface_contribution(outdir / "radiation_pitch_surface_contribution.csv",
                                              omega,
                                              dof_col,
                                              solver,
                                              sol,
                                              float_faces,
                                              radiation_ref_point,
                                              rho,
                                              moment_sign);
}

double rigid_mode_normal_component(int dof, const Tddd& x, const Tddd& normal, const Tddd& ref_point, double moment_sign = 1.0) {
  if (dof < 0 || dof > 5)
    return 0.0;
  if (dof <= 2)
    return normal[dof];
  const Tddd r = x - ref_point;
  const Tddd moment_mode = Cross(r, normal);
  return moment_sign * moment_mode[dof - 3];
}

std::array<Complex, 6> integrate_reciprocity_from_face_phi(const std::unordered_set<networkFace*>& float_faces,
                                                           const Tddd& ref_point,
                                                           double moment_sign,
                                                           const std::function<std::array<Complex, 3>(const networkFace*)>& phi_at_vertices) {
  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});

  // Bilinear diagnostic before pressure scaling:
  // S_ij = int_S phi_j m_i dS, where force modes m_i are n_i and (r x n)_(i-3).
  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};

  for (auto* f : float_faces) {
    const auto phi = phi_at_vertices(f);
    auto [p0, p1, p2] = f->getPoints();
    const auto& x0 = p0->X;
    const auto& x1 = p1->X;
    const auto& x2 = p2->X;
    const auto& n = f->normal;

    for (const auto& l : bary) {
      const double l0 = l[0], l1 = l[1], l2 = l[2];
      const Complex phi_q = l0 * phi[0] + l1 * phi[1] + l2 * phi[2];
      const Tddd xq = l0 * x0 + l1 * x1 + l2 * x2;
      for (int dof_row = 0; dof_row < 6; ++dof_row)
        out[static_cast<std::size_t>(dof_row)] += w * phi_q * rigid_mode_normal_component(dof_row, xq, n, ref_point, moment_sign) * f->area;
    }
  }

  return out;
}

std::array<Complex, 6> integrate_reciprocity_from_mesh(const std::unordered_set<networkFace*>& float_faces, const Tddd& ref_point, double moment_sign = 1.0) {
  return integrate_reciprocity_from_face_phi(float_faces, ref_point, moment_sign, [](const networkFace* f) {
    auto [p0, p1, p2] = f->getPoints();
    return std::array<Complex, 3>{
        Complex{get_phi_on_face(p0, f), 0.0},
        Complex{get_phi_on_face(p1, f), 0.0},
        Complex{get_phi_on_face(p2, f), 0.0},
    };
  });
}

std::array<Complex, 6> integrate_reciprocity_from_solution(const bem_frequency_domain::LinearSolution& sol,
                                                           const std::unordered_set<networkFace*>& float_faces,
                                                           const Tddd& ref_point,
                                                           double moment_sign = 1.0) {
  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});
  constexpr std::array<bool, 3> all_true{true, true, true};
  for (auto* f : float_faces) {
    const auto it = sol.face_field.find(f);
    if (it == sol.face_field.end())
      continue;
    const auto& field = it->second;
    if (!bem_frequency_domain::face_field_is_quadratic(f, field))
      continue;
    const auto p6 = bem_frequency_domain::phi6(field);
    const auto X6 = bem_frequency_domain::face_x6(f);
    for (const auto& [x0q, x1q, w0w1] : __GWGW10__Tuple) {
      const auto bary = ModTriShape<3>(x0q, x1q);
      const double b0 = bary[0], b1 = bary[1];
      const auto N6 = f->trueQuadN6(b0, b1);
      const auto N6_geo = TriShape<6>(b0, b1, all_true);
      const auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
      const auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
      const Complex phi_q = bem_frequency_domain::dot_shape6(N6, p6);
      const Tddd xq = Dot(N6_geo, X6);
      const Tddd area_vec = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
      const double weight = w0w1 * (1.0 - x0q);
      out[0] += phi_q * area_vec[0] * weight;
      out[1] += phi_q * area_vec[1] * weight;
      out[2] += phi_q * area_vec[2] * weight;
      const Tddd r = xq - ref_point;
      out[3] += moment_sign * phi_q * (r[1] * area_vec[2] - r[2] * area_vec[1]) * weight;
      out[4] += moment_sign * phi_q * (r[2] * area_vec[0] - r[0] * area_vec[2]) * weight;
      out[5] += moment_sign * phi_q * (r[0] * area_vec[1] - r[1] * area_vec[0]) * weight;
    }
  }
  bool used_quadratic = false;
  for (auto* f : float_faces) {
    const auto it = sol.face_field.find(f);
    if (it != sol.face_field.end() && bem_frequency_domain::face_field_is_quadratic(f, it->second)) {
      used_quadratic = true;
      break;
    }
  }
  if (used_quadratic) {
    const auto linear_part = integrate_reciprocity_from_face_phi(float_faces, ref_point, moment_sign, [&](const networkFace* f) {
      const auto it = sol.face_field.find(f);
      if (it == sol.face_field.end() || bem_frequency_domain::face_field_is_quadratic(f, it->second))
        return std::array<Complex, 3>{Complex{0.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}};
      return it->second.phi;
    });
    for (std::size_t i = 0; i < out.size(); ++i)
      out[i] += linear_part[i];
    return out;
  }

  return integrate_reciprocity_from_face_phi(float_faces, ref_point, moment_sign, [&](const networkFace* f) {
    const auto it = sol.face_field.find(f);
    if (it == sol.face_field.end())
      return std::array<Complex, 3>{Complex{0.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}};
    return it->second.phi;
  });
}

std::array<Complex, 6> integrate_reciprocity_from_face_phi_lumped(const std::unordered_set<networkFace*>& float_faces,
                                                                  const Tddd& ref_point,
                                                                  double moment_sign,
                                                                  const std::function<std::array<Complex, 3>(const networkFace*)>& phi_at_vertices) {
  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});

  // Mass-lumped diagnostic matching a nodal/collocation view:
  // S_i = sum_faces area/3 * sum_vertices phi(vertex) * m_i(vertex).
  for (auto* f : float_faces) {
    const auto phi = phi_at_vertices(f);
    auto [p0, p1, p2] = f->getPoints();
    const std::array<networkPoint*, 3> pts{p0, p1, p2};
    const double w = f->area / 3.0;
    for (std::size_t a = 0; a < pts.size(); ++a) {
      for (int dof_row = 0; dof_row < 6; ++dof_row) {
        out[static_cast<std::size_t>(dof_row)] +=
            w * phi[a] * rigid_mode_normal_component(dof_row, pts[a]->X, f->normal, ref_point, moment_sign);
      }
    }
  }

  return out;
}

std::array<Complex, 6> integrate_reciprocity_from_solution_lumped(const bem_frequency_domain::LinearSolution& sol,
                                                                  const std::unordered_set<networkFace*>& float_faces,
                                                                  const Tddd& ref_point,
                                                                  double moment_sign = 1.0) {
  bool used_quadratic = false;
  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});
  for (auto* f : float_faces) {
    const auto it = sol.face_field.find(f);
    if (it == sol.face_field.end() || !bem_frequency_domain::face_field_is_quadratic(f, it->second))
      continue;
    used_quadratic = true;
    const auto& field = it->second;
    auto [p0, p1, p2] = f->getPoints();
    auto [l0, l1, l2] = f->getLines();
    const std::array<Tddd, 6> xs{p0->X, p1->X, p2->X, l0->X_mid, l1->X_mid, l2->X_mid};
    const auto phi = bem_frequency_domain::phi6(field);
    // Area-lumped diagnostic for six true-quadratic nodes.  This is not used for coefficients.
    const double w = f->area / 6.0;
    for (std::size_t a = 0; a < xs.size(); ++a) {
      for (int dof_row = 0; dof_row < 6; ++dof_row) {
        out[static_cast<std::size_t>(dof_row)] +=
            w * phi[a] * rigid_mode_normal_component(dof_row, xs[a], f->normal, ref_point, moment_sign);
      }
    }
  }
  if (used_quadratic) {
    const auto linear_part = integrate_reciprocity_from_face_phi_lumped(float_faces, ref_point, moment_sign, [&](const networkFace* f) {
      const auto it = sol.face_field.find(f);
      if (it == sol.face_field.end() || bem_frequency_domain::face_field_is_quadratic(f, it->second))
        return std::array<Complex, 3>{Complex{0.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}};
      return it->second.phi;
    });
    for (std::size_t i = 0; i < out.size(); ++i)
      out[i] += linear_part[i];
    return out;
  }

  return integrate_reciprocity_from_face_phi_lumped(float_faces, ref_point, moment_sign, [&](const networkFace* f) {
    const auto it = sol.face_field.find(f);
    if (it == sol.face_field.end())
      return std::array<Complex, 3>{Complex{0.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}};
    return it->second.phi;
  });
}

Complex integrate_green_identity_pair(const bem_frequency_domain::LinearSolution& sol_i,
                                      const bem_frequency_domain::LinearSolution& sol_j,
                                      const std::unordered_set<networkFace*>& faces) {
  Complex out{0.0, 0.0};

  // Green's second identity diagnostic:
  // int_Gamma (phi_j * phin_i - phi_i * phin_j) dS.
  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};

  for (auto* f : faces) {
    const auto iti = sol_i.face_field.find(f);
    const auto itj = sol_j.face_field.find(f);
    if (iti == sol_i.face_field.end() || itj == sol_j.face_field.end())
      continue;
    const auto& fi = iti->second;
    const auto& fj = itj->second;
    if (bem_frequency_domain::face_field_is_quadratic(f, fi) && bem_frequency_domain::face_field_is_quadratic(f, fj)) {
      constexpr std::array<bool, 3> all_true{true, true, true};
      const auto X6 = bem_frequency_domain::face_x6(f);
      const auto phi_i6 = bem_frequency_domain::phi6(fi);
      const auto phi_j6 = bem_frequency_domain::phi6(fj);
      const auto phin_i6 = bem_frequency_domain::phin6(fi);
      const auto phin_j6 = bem_frequency_domain::phin6(fj);
      for (const auto& [x0q, x1q, w0w1] : __GWGW10__Tuple) {
        const auto bary = ModTriShape<3>(x0q, x1q);
        const double b0 = bary[0], b1 = bary[1];
        const auto N6 = f->trueQuadN6(b0, b1);
        const auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
        const auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
        const double dS = Norm(Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6))) * w0w1 * (1.0 - x0q);
        const Complex phi_i = bem_frequency_domain::dot_shape6(N6, phi_i6);
        const Complex phi_j = bem_frequency_domain::dot_shape6(N6, phi_j6);
        const Complex phin_i = bem_frequency_domain::dot_shape6(N6, phin_i6);
        const Complex phin_j = bem_frequency_domain::dot_shape6(N6, phin_j6);
        out += (phi_j * phin_i - phi_i * phin_j) * dS;
      }
      continue;
    }
    for (const auto& l : bary) {
      const double l0 = l[0], l1 = l[1], l2 = l[2];
      const Complex phi_i = l0 * fi.phi[0] + l1 * fi.phi[1] + l2 * fi.phi[2];
      const Complex phi_j = l0 * fj.phi[0] + l1 * fj.phi[1] + l2 * fj.phi[2];
      const Complex phin_i = l0 * fi.phin[0] + l1 * fi.phin[1] + l2 * fi.phin[2];
      const Complex phin_j = l0 * fj.phin[0] + l1 * fj.phin[1] + l2 * fj.phin[2];
      out += w * (phi_j * phin_i - phi_i * phin_j) * f->area;
    }
  }

  return out;
}

void append_radiation_green_identity(const std::filesystem::path& path,
                                     double omega,
                                     const char* solver_name,
                                     const std::unordered_map<int, bem_frequency_domain::LinearSolution>& solutions,
                                     const std::unordered_set<networkFace*>& all_faces,
                                     const FaceSets& face_sets) {
  if (solutions.size() < 2)
    return;

  std::unordered_set<networkFace*> wall_faces;
  wall_faces.reserve(all_faces.size());
  for (auto* f : all_faces) {
    if (!face_sets.float_surface.contains(f) && !face_sets.free_surface.contains(f))
      wall_faces.emplace(f);
  }

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header)
    fs << "omega,solver,dof_i,dof_j,surface,value_re,value_im,abs\n";

  auto write_surface = [&](int i, int j, const std::string& surface, const std::unordered_set<networkFace*>& faces) {
    const auto value = integrate_green_identity_pair(solutions.at(i), solutions.at(j), faces);
    fs << std::scientific << std::setprecision(12)
       << omega << "," << solver_name << "," << i << "," << j << "," << surface << ","
       << value.real() << "," << value.imag() << "," << std::abs(value) << "\n";
  };

  std::vector<int> keys;
  keys.reserve(solutions.size());
  for (const auto& [dof, _] : solutions)
    keys.push_back(dof);
  std::sort(keys.begin(), keys.end());

  for (std::size_t a = 0; a < keys.size(); ++a) {
    for (std::size_t b = a + 1; b < keys.size(); ++b) {
      const int i = keys[a];
      const int j = keys[b];
      write_surface(i, j, "float", face_sets.float_surface);
      write_surface(i, j, "free", face_sets.free_surface);
      write_surface(i, j, "wall", wall_faces);
      write_surface(i, j, "all", all_faces);
    }
  }
}

void append_radiation_effective_ntd_map(const std::filesystem::path& path,
                                        double omega,
                                        const char* solver_name,
                                        const std::unordered_map<int, bem_frequency_domain::LinearSolution>& solutions,
                                        const std::unordered_set<networkFace*>& all_faces,
                                        const FaceSets& face_sets,
                                        const Tddd& ref_point,
                                        double moment_sign) {
  if (solutions.size() < 2)
    return;

  std::unordered_set<networkFace*> wall_faces;
  wall_faces.reserve(all_faces.size());
  for (auto* f : all_faces) {
    if (!face_sets.float_surface.contains(f) && !face_sets.free_surface.contains(f))
      wall_faces.emplace(f);
  }

  std::vector<int> keys;
  keys.reserve(solutions.size());
  for (const auto& [dof, _] : solutions)
    keys.push_back(dof);
  std::sort(keys.begin(), keys.end());

  std::unordered_map<int, std::array<Complex, 6>> s_cols;
  std::unordered_map<int, std::array<Complex, 6>> s_lumped_cols;
  for (int dof : keys) {
    const auto& sol = solutions.at(dof);
    s_cols.emplace(dof, integrate_reciprocity_from_solution(sol, face_sets.float_surface, ref_point, moment_sign));
    s_lumped_cols.emplace(dof, integrate_reciprocity_from_solution_lumped(sol, face_sets.float_surface, ref_point, moment_sign));
  }

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header)
    fs << "omega,solver,dof_i,dof_j,"
          "Sij_re,Sij_im,Sji_re,Sji_im,asym_re,asym_im,asym_abs,asym_rel,"
          "asym_lumped_re,asym_lumped_im,asym_lumped_abs,asym_lumped_rel,"
          "green_float_re,green_float_im,green_float_abs,"
          "green_free_re,green_free_im,green_free_abs,"
          "green_wall_re,green_wall_im,green_wall_abs,"
          "green_all_re,green_all_im,green_all_abs,"
          "closure_float_abs,closure_all_sum_abs\n";

  auto rel = [](Complex a, Complex b, Complex diff) {
    return std::abs(diff) / std::max({std::abs(a), std::abs(b), 1e-300});
  };

  for (std::size_t a = 0; a < keys.size(); ++a) {
    for (std::size_t b = a + 1; b < keys.size(); ++b) {
      const int i = keys[a];
      const int j = keys[b];
      const Complex Sij = s_cols.at(j)[static_cast<std::size_t>(i)];
      const Complex Sji = s_cols.at(i)[static_cast<std::size_t>(j)];
      const Complex asym = Sij - Sji;
      const Complex Sij_lumped = s_lumped_cols.at(j)[static_cast<std::size_t>(i)];
      const Complex Sji_lumped = s_lumped_cols.at(i)[static_cast<std::size_t>(j)];
      const Complex asym_lumped = Sij_lumped - Sji_lumped;

      const auto green_float = integrate_green_identity_pair(solutions.at(i), solutions.at(j), face_sets.float_surface);
      const auto green_free = integrate_green_identity_pair(solutions.at(i), solutions.at(j), face_sets.free_surface);
      const auto green_wall = integrate_green_identity_pair(solutions.at(i), solutions.at(j), wall_faces);
      const auto green_all = integrate_green_identity_pair(solutions.at(i), solutions.at(j), all_faces);
      const Complex green_sum = green_float + green_free + green_wall;

      fs << std::scientific << std::setprecision(12)
         << omega << "," << solver_name << "," << i << "," << j << ","
         << Sij.real() << "," << Sij.imag() << ","
         << Sji.real() << "," << Sji.imag() << ","
         << asym.real() << "," << asym.imag() << "," << std::abs(asym) << "," << rel(Sij, Sji, asym) << ","
         << asym_lumped.real() << "," << asym_lumped.imag() << "," << std::abs(asym_lumped) << "," << rel(Sij_lumped, Sji_lumped, asym_lumped) << ","
         << green_float.real() << "," << green_float.imag() << "," << std::abs(green_float) << ","
         << green_free.real() << "," << green_free.imag() << "," << std::abs(green_free) << ","
         << green_wall.real() << "," << green_wall.imag() << "," << std::abs(green_wall) << ","
         << green_all.real() << "," << green_all.imag() << "," << std::abs(green_all) << ","
         << std::abs(asym - green_float) << "," << std::abs(green_all - green_sum) << "\n";
    }
  }
}

void append_radiation_effective_ntd_face_contrib(const std::filesystem::path& path,
                                                 double omega,
                                                 const char* solver_name,
                                                 const std::unordered_map<int, bem_frequency_domain::LinearSolution>& solutions,
                                                 const FaceSets& face_sets,
                                                 const Tddd& ref_point,
                                                 double moment_sign) {
  if (solutions.size() < 2)
    return;

  std::vector<int> keys;
  keys.reserve(solutions.size());
  for (const auto& [dof, _] : solutions)
    keys.push_back(dof);
  std::sort(keys.begin(), keys.end());

  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header)
    fs << "omega,solver,dof_i,dof_j,face_index,cx,cy,cz,nx,ny,nz,area,"
          "green_re,green_im,green_abs,"
          "Sij_face_re,Sij_face_im,Sji_face_re,Sji_face_im,"
          "asym_face_re,asym_face_im,asym_face_abs,closure_face_abs\n";

  for (std::size_t a = 0; a < keys.size(); ++a) {
    for (std::size_t b = a + 1; b < keys.size(); ++b) {
      const int i = keys[a];
      const int j = keys[b];
      const auto& sol_i = solutions.at(i);
      const auto& sol_j = solutions.at(j);

      for (auto* f : face_sets.float_surface) {
        const auto iti = sol_i.face_field.find(f);
        const auto itj = sol_j.face_field.find(f);
        if (iti == sol_i.face_field.end() || itj == sol_j.face_field.end())
          continue;

        auto [p0, p1, p2] = f->getPoints();
        const auto& x0 = p0->X;
        const auto& x1 = p1->X;
        const auto& x2 = p2->X;
        const auto& fi = iti->second;
        const auto& fj = itj->second;
        Complex green{0.0, 0.0};
        Complex Sij{0.0, 0.0};
        Complex Sji{0.0, 0.0};

        for (const auto& l : bary) {
          const double l0 = l[0], l1 = l[1], l2 = l[2];
          const Tddd xq = l0 * x0 + l1 * x1 + l2 * x2;
          const Complex phi_i = l0 * fi.phi[0] + l1 * fi.phi[1] + l2 * fi.phi[2];
          const Complex phi_j = l0 * fj.phi[0] + l1 * fj.phi[1] + l2 * fj.phi[2];
          const Complex phin_i = l0 * fi.phin[0] + l1 * fi.phin[1] + l2 * fi.phin[2];
          const Complex phin_j = l0 * fj.phin[0] + l1 * fj.phin[1] + l2 * fj.phin[2];
          const double mi = rigid_mode_normal_component(i, xq, f->normal, ref_point, moment_sign);
          const double mj = rigid_mode_normal_component(j, xq, f->normal, ref_point, moment_sign);
          green += w * (phi_j * phin_i - phi_i * phin_j) * f->area;
          Sij += w * phi_j * mi * f->area;
          Sji += w * phi_i * mj * f->area;
        }

        const Complex asym = Sij - Sji;
        fs << std::scientific << std::setprecision(12)
           << omega << "," << solver_name << "," << i << "," << j << ","
           << f->index << ","
           << f->centroid[0] << "," << f->centroid[1] << "," << f->centroid[2] << ","
           << f->normal[0] << "," << f->normal[1] << "," << f->normal[2] << ","
           << f->area << ","
           << green.real() << "," << green.imag() << "," << std::abs(green) << ","
           << Sij.real() << "," << Sij.imag() << ","
           << Sji.real() << "," << Sji.imag() << ","
           << asym.real() << "," << asym.imag() << "," << std::abs(asym) << ","
           << std::abs(asym - green) << "\n";
      }
    }
  }
}

int radiation_surface_kind(const BEM_DOF_Base* node, const networkFace* f, const FaceSets& face_sets);
const char* radiation_surface_name(int k);

void append_radiation_algebraic_reciprocity(const std::filesystem::path& path,
                                            double omega,
                                            const char* solver_name,
                                            const std::unordered_map<int, std::vector<Complex>>& rhs_by_dof,
                                            const std::unordered_map<int, std::vector<Complex>>& u_by_dof,
                                            const std::unordered_map<int, bem_frequency_domain::LinearSolution>& solutions,
                                            const FaceSets& face_sets,
                                            const Tddd& ref_point,
                                            double moment_sign) {
  if (rhs_by_dof.size() < 2 || u_by_dof.size() < 2)
    return;

  std::vector<int> keys;
  keys.reserve(rhs_by_dof.size());
  for (const auto& [dof, rhs] : rhs_by_dof) {
    const auto uit = u_by_dof.find(dof);
    if (uit == u_by_dof.end() || uit->second.size() != rhs.size())
      continue;
    keys.push_back(dof);
  }
  std::sort(keys.begin(), keys.end());
  if (keys.size() < 2)
    return;

  auto dot = [](const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex out{0.0, 0.0};
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i)
      out += a[i] * b[i];
    return out;
  };
  auto hdot = [](const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex out{0.0, 0.0};
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i)
      out += std::conj(a[i]) * b[i];
    return out;
  };
  auto rel = [](Complex a, Complex b, Complex diff) {
    return std::abs(diff) / std::max({std::abs(a), std::abs(b), 1e-300});
  };

  std::unordered_map<int, std::array<Complex, 6>> s_cols;
  for (int dof : keys) {
    const auto sit = solutions.find(dof);
    if (sit != solutions.end())
      s_cols.emplace(dof, integrate_reciprocity_from_solution(sit->second, face_sets.float_surface, ref_point, moment_sign));
  }

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,solver,dof_i,dof_j,n,"
          "bi_dot_uj_re,bi_dot_uj_im,bj_dot_ui_re,bj_dot_ui_im,dot_asym_re,dot_asym_im,dot_asym_abs,dot_asym_rel,"
          "bi_hdot_uj_re,bi_hdot_uj_im,bj_hdot_ui_re,bj_hdot_ui_im,hdot_asym_re,hdot_asym_im,hdot_asym_abs,hdot_asym_rel,"
          "surface_Sij_re,surface_Sij_im,surface_Sji_re,surface_Sji_im,surface_asym_abs,surface_asym_rel\n";
  }

  for (std::size_t a = 0; a < keys.size(); ++a) {
    for (std::size_t b = a + 1; b < keys.size(); ++b) {
      const int i = keys[a];
      const int j = keys[b];
      const auto& rhs_i = rhs_by_dof.at(i);
      const auto& rhs_j = rhs_by_dof.at(j);
      const auto& u_i = u_by_dof.at(i);
      const auto& u_j = u_by_dof.at(j);
      const Complex bi_uj = dot(rhs_i, u_j);
      const Complex bj_ui = dot(rhs_j, u_i);
      const Complex dot_asym = bi_uj - bj_ui;
      const Complex hbi_uj = hdot(rhs_i, u_j);
      const Complex hbj_ui = hdot(rhs_j, u_i);
      const Complex hdot_asym = hbi_uj - hbj_ui;

      Complex Sij{0.0, 0.0};
      Complex Sji{0.0, 0.0};
      double surface_asym_abs = 0.0;
      double surface_asym_rel = 0.0;
      const auto si = s_cols.find(j);
      const auto sj = s_cols.find(i);
      if (si != s_cols.end() && sj != s_cols.end()) {
        Sij = si->second[static_cast<std::size_t>(i)];
        Sji = sj->second[static_cast<std::size_t>(j)];
        const Complex surface_asym = Sij - Sji;
        surface_asym_abs = std::abs(surface_asym);
        surface_asym_rel = rel(Sij, Sji, surface_asym);
      }

      fs << std::scientific << std::setprecision(12)
         << omega << "," << solver_name << "," << i << "," << j << "," << rhs_i.size() << ","
         << bi_uj.real() << "," << bi_uj.imag() << ","
         << bj_ui.real() << "," << bj_ui.imag() << ","
         << dot_asym.real() << "," << dot_asym.imag() << ","
         << std::abs(dot_asym) << "," << rel(bi_uj, bj_ui, dot_asym) << ","
         << hbi_uj.real() << "," << hbi_uj.imag() << ","
         << hbj_ui.real() << "," << hbj_ui.imag() << ","
         << hdot_asym.real() << "," << hdot_asym.imag() << ","
         << std::abs(hdot_asym) << "," << rel(hbi_uj, hbj_ui, hdot_asym) << ","
         << Sij.real() << "," << Sij.imag() << ","
         << Sji.real() << "," << Sji.imag() << ","
         << surface_asym_abs << "," << surface_asym_rel << "\n";
    }
  }
}

void append_radiation_reciprocity_pair_summary(const std::filesystem::path& path,
                                               double omega,
                                               const char* solver_name,
                                               const std::unordered_map<int, std::vector<Complex>>& rhs_by_dof,
                                               const std::unordered_map<int, std::vector<Complex>>& u_by_dof,
                                               const std::unordered_map<int, bem_frequency_domain::LinearSolution>& solutions,
                                               const std::unordered_set<networkFace*>& all_faces,
                                               const FaceSets& face_sets,
                                               const Tddd& ref_point,
                                               double moment_sign) {
  if (rhs_by_dof.size() < 2 || u_by_dof.size() < 2 || solutions.size() < 2)
    return;

  std::vector<int> keys;
  keys.reserve(rhs_by_dof.size());
  for (const auto& [dof, rhs] : rhs_by_dof) {
    const auto uit = u_by_dof.find(dof);
    const auto sit = solutions.find(dof);
    if (uit == u_by_dof.end() || sit == solutions.end() || uit->second.size() != rhs.size())
      continue;
    keys.push_back(dof);
  }
  std::sort(keys.begin(), keys.end());
  if (keys.size() < 2)
    return;

  auto dot = [](const std::vector<Complex>& a, const std::vector<Complex>& b) {
    Complex out{0.0, 0.0};
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i)
      out += a[i] * b[i];
    return out;
  };
  auto rel = [](Complex a, Complex b, Complex diff) {
    return std::abs(diff) / std::max({std::abs(a), std::abs(b), 1e-300});
  };
  auto classify = [](double surface_abs,
                     double algebraic_abs,
                     double algebraic_closure_abs,
                     double green_float_abs,
                     double green_other_abs,
                     double green_closure_abs) -> const char* {
    const double scale = std::max({surface_abs, algebraic_abs, 1e-300});
    if (surface_abs <= 1e-8 * scale && algebraic_abs <= 1e-8 * scale && green_float_abs <= 1e-8 * scale)
      return "symmetric_within_roundoff";
    if (surface_abs > 0.0 && green_closure_abs <= 0.1 * surface_abs && green_other_abs <= 1e-6 * surface_abs)
      return "surface_asym_is_float_green_identity_residual";
    if (surface_abs > 0.0 && green_other_abs > 1e-6 * surface_abs)
      return "other_boundary_green_contribution";
    if (algebraic_abs >= 0.5 * surface_abs && algebraic_closure_abs <= 0.5 * scale)
      return "linear_system_asymmetry_tracks_surface";
    if (green_closure_abs >= 0.5 * std::max({surface_abs, green_float_abs, 1e-300}))
      return "surface_green_closure_gap";
    if (surface_abs >= 10.0 * std::max(algebraic_abs, 1e-300))
      return "surface_projection_or_mode_dominant";
    if (algebraic_closure_abs >= 0.5 * scale)
      return "surface_vs_algebraic_closure_gap";
    return "mixed";
  };

  std::unordered_set<networkFace*> wall_faces;
  wall_faces.reserve(all_faces.size());
  for (auto* f : all_faces) {
    if (!face_sets.float_surface.contains(f) && !face_sets.free_surface.contains(f))
      wall_faces.emplace(f);
  }

  std::unordered_map<int, std::array<Complex, 6>> s_cols;
  for (int dof : keys)
    s_cols.emplace(dof, integrate_reciprocity_from_solution(solutions.at(dof), face_sets.float_surface, ref_point, moment_sign));

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,solver,dof_i,dof_j,n,"
          "surface_asym_re,surface_asym_im,surface_asym_abs,surface_asym_rel,"
          "algebraic_asym_re,algebraic_asym_im,algebraic_asym_abs,algebraic_asym_rel,"
          "algebraic_closure_re,algebraic_closure_im,algebraic_closure_abs,"
          "green_float_re,green_float_im,green_float_abs,"
          "green_free_re,green_free_im,green_free_abs,"
          "green_wall_re,green_wall_im,green_wall_abs,"
          "green_all_re,green_all_im,green_all_abs,"
          "green_float_closure_re,green_float_closure_im,green_float_closure_abs,classification,"
          "Sij_re,Sij_im,Sji_re,Sji_im,bi_dot_uj_re,bi_dot_uj_im,bj_dot_ui_re,bj_dot_ui_im\n";
  }

  for (std::size_t a = 0; a < keys.size(); ++a) {
    for (std::size_t b = a + 1; b < keys.size(); ++b) {
      const int i = keys[a];
      const int j = keys[b];
      const auto& rhs_i = rhs_by_dof.at(i);
      const auto& rhs_j = rhs_by_dof.at(j);
      const auto& u_i = u_by_dof.at(i);
      const auto& u_j = u_by_dof.at(j);

      const Complex bi_uj = dot(rhs_i, u_j);
      const Complex bj_ui = dot(rhs_j, u_i);
      const Complex algebraic_asym = bi_uj - bj_ui;
      const Complex Sij = s_cols.at(j)[static_cast<std::size_t>(i)];
      const Complex Sji = s_cols.at(i)[static_cast<std::size_t>(j)];
      const Complex surface_asym = Sij - Sji;
      const Complex algebraic_closure = surface_asym - algebraic_asym;
      const Complex green_float = integrate_green_identity_pair(solutions.at(i), solutions.at(j), face_sets.float_surface);
      const Complex green_free = integrate_green_identity_pair(solutions.at(i), solutions.at(j), face_sets.free_surface);
      const Complex green_wall = integrate_green_identity_pair(solutions.at(i), solutions.at(j), wall_faces);
      const Complex green_all = integrate_green_identity_pair(solutions.at(i), solutions.at(j), all_faces);
      const Complex green_float_closure = surface_asym - green_float;
      const double green_other_abs = std::abs(green_free) + std::abs(green_wall);
      const char* klass = classify(std::abs(surface_asym),
                                   std::abs(algebraic_asym),
                                   std::abs(algebraic_closure),
                                   std::abs(green_float),
                                   green_other_abs,
                                   std::abs(green_float_closure));

      fs << std::scientific << std::setprecision(12)
         << omega << "," << solver_name << "," << i << "," << j << "," << rhs_i.size() << ","
         << surface_asym.real() << "," << surface_asym.imag() << ","
         << std::abs(surface_asym) << "," << rel(Sij, Sji, surface_asym) << ","
         << algebraic_asym.real() << "," << algebraic_asym.imag() << ","
         << std::abs(algebraic_asym) << "," << rel(bi_uj, bj_ui, algebraic_asym) << ","
         << algebraic_closure.real() << "," << algebraic_closure.imag() << "," << std::abs(algebraic_closure) << ","
         << green_float.real() << "," << green_float.imag() << "," << std::abs(green_float) << ","
         << green_free.real() << "," << green_free.imag() << "," << std::abs(green_free) << ","
         << green_wall.real() << "," << green_wall.imag() << "," << std::abs(green_wall) << ","
         << green_all.real() << "," << green_all.imag() << "," << std::abs(green_all) << ","
         << green_float_closure.real() << "," << green_float_closure.imag() << "," << std::abs(green_float_closure) << ","
         << klass << ","
         << Sij.real() << "," << Sij.imag() << ","
         << Sji.real() << "," << Sji.imag() << ","
         << bi_uj.real() << "," << bi_uj.imag() << ","
         << bj_ui.real() << "," << bj_ui.imag() << "\n";
    }
  }
}

void append_radiation_adjoint_residual_summary(const std::filesystem::path& path,
                                               double omega,
                                               int dof_col,
                                               const char* solver_name,
                                               const bem_frequency_domain::Solution& sol,
                                               const FaceSets& face_sets) {
  if (sol.adjoint_residual.empty() || sol.adjoint_residual.size() != sol.id_by_index.size())
    return;

  struct Stat {
    std::uint64_t samples = 0;
    long double sum_abs = 0.0L;
    long double l2 = 0.0L;
    double max_abs = 0.0;
    std::size_t worst_index = 0;
    Complex worst_value{0.0, 0.0};
  };

  constexpr int surface_count = 4;
  constexpr int key_count = 2; // face-keyed or smooth/default key
  std::array<Stat, surface_count * key_count> stats{};
  long double total_l2 = 0.0L;

  auto stat_index = [](int surface, int key_kind) { return surface * key_count + key_kind; };
  auto key_name = [](int key_kind) { return key_kind == 0 ? "face_key" : "default_key"; };

  for (std::size_t i = 0; i < sol.adjoint_residual.size(); ++i) {
    const auto* node = std::get<0>(sol.id_by_index[i]);
    const auto* f = std::get<1>(sol.id_by_index[i]);
    const int surface = radiation_surface_kind(node, f, face_sets);
    const int key_kind = f ? 0 : 1;
    const Complex value = sol.adjoint_residual[i];
    const double av = std::abs(value);
    total_l2 += static_cast<long double>(av) * av;
    auto& st = stats[static_cast<std::size_t>(stat_index(surface, key_kind))];
    ++st.samples;
    st.sum_abs += av;
    st.l2 += static_cast<long double>(av) * av;
    if (av > st.max_abs) {
      st.max_abs = av;
      st.worst_index = i;
      st.worst_value = value;
    }
  }

  const double total_l2_norm = std::sqrt(static_cast<double>(total_l2));
  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,component_surface,key_kind,samples,sum_abs,l2_norm,relative_l2,max_abs,"
          "worst_index,worst_re,worst_im,total_l2_norm\n";
  }

  for (int surface = 0; surface < surface_count; ++surface) {
    for (int key_kind = 0; key_kind < key_count; ++key_kind) {
      const auto& st = stats[static_cast<std::size_t>(stat_index(surface, key_kind))];
      if (st.samples == 0)
        continue;
      const double l2_norm = std::sqrt(static_cast<double>(st.l2));
      fs << std::scientific << std::setprecision(12)
         << omega << "," << dof_col << "," << solver_name << ","
         << radiation_surface_name(surface) << "," << key_name(key_kind) << ","
         << st.samples << "," << static_cast<double>(st.sum_abs) << ","
         << l2_norm << ","
         << (total_l2_norm > 0.0 ? l2_norm / total_l2_norm : 0.0) << ","
         << st.max_abs << "," << st.worst_index << ","
         << st.worst_value.real() << "," << st.worst_value.imag() << ","
         << total_l2_norm << "\n";
    }
  }
}

void append_radiation_bie_residual_summary(const std::filesystem::path& path,
                                           double omega,
                                           int dof_col,
                                           const char* solver_name,
                                           const bem_frequency_domain::Solution& sol,
                                           const FaceSets& face_sets,
                                           std::function<bool(const BEM_DOF_Base&, const networkFace*)> use_bie_row_for_interface) {
  if (sol.bie_residual.empty() || sol.bie_residual.size() != sol.id_by_index.size())
    return;

  struct Stat {
    std::uint64_t samples = 0;
    long double sum_abs = 0.0L;
    long double l2 = 0.0L;
    double max_abs = 0.0;
    std::size_t worst_index = 0;
    Complex worst_value{0.0, 0.0};
  };

  constexpr int surface_count = 4;
  constexpr int key_count = 2;      // face-keyed or smooth/default key
  constexpr int row_kind_count = 2; // BIE row or interface constraint row
  std::array<Stat, surface_count * key_count * row_kind_count> stats{};
  long double total_l2 = 0.0L;

  auto stat_index = [](int surface, int key_kind, int row_kind) {
    return (surface * key_count + key_kind) * row_kind_count + row_kind;
  };
  auto key_name = [](int key_kind) { return key_kind == 0 ? "face_key" : "default_key"; };
  const bool have_solve_row_kind = sol.interface_constraint_row.size() == sol.bie_residual.size();
  auto row_kind = [&](std::size_t i, const BEM_DOF_Base* node, const networkFace* f) {
    if (have_solve_row_kind)
      return sol.interface_constraint_row[i] ? 1 : 0;
    const bool forced_bie = node && use_bie_row_for_interface && use_bie_row_for_interface(*node, f);
    // solve_linear_bvp restores boundary flags before this post-solve summary is written,
    // so use the stable face-keyed interface layout rather than the restored Neumann flag.
    return (node && node->BCInterface && f && !forced_bie) ? 1 : 0;
  };
  auto row_kind_name = [](int k) { return k == 1 ? "interface_constraint" : "bie"; };

  for (std::size_t i = 0; i < sol.bie_residual.size(); ++i) {
    const auto* node = std::get<0>(sol.id_by_index[i]);
    const auto* f = std::get<1>(sol.id_by_index[i]);
    const int surface = radiation_surface_kind(node, f, face_sets);
    const int key_kind = f ? 0 : 1;
    const int row = row_kind(i, node, f);
    const Complex value = sol.bie_residual[i];
    const double av = std::abs(value);
    total_l2 += static_cast<long double>(av) * av;
    auto& st = stats[static_cast<std::size_t>(stat_index(surface, key_kind, row))];
    ++st.samples;
    st.sum_abs += av;
    st.l2 += static_cast<long double>(av) * av;
    if (av > st.max_abs) {
      st.max_abs = av;
      st.worst_index = i;
      st.worst_value = value;
    }
  }

  const double total_l2_norm = std::sqrt(static_cast<double>(total_l2));
  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,component_surface,key_kind,row_kind,samples,sum_abs,l2_norm,relative_l2,max_abs,"
          "worst_index,worst_re,worst_im,total_l2_norm,reported_bie_residual_l2\n";
  }

  for (int surface = 0; surface < surface_count; ++surface) {
    for (int key_kind = 0; key_kind < key_count; ++key_kind) {
      for (int row = 0; row < row_kind_count; ++row) {
        const auto& st = stats[static_cast<std::size_t>(stat_index(surface, key_kind, row))];
        if (st.samples == 0)
          continue;
        const double l2_norm = std::sqrt(static_cast<double>(st.l2));
        fs << std::scientific << std::setprecision(12)
           << omega << "," << dof_col << "," << solver_name << ","
           << radiation_surface_name(surface) << "," << key_name(key_kind) << "," << row_kind_name(row) << ","
           << st.samples << "," << static_cast<double>(st.sum_abs) << ","
           << l2_norm << ","
           << (total_l2_norm > 0.0 ? l2_norm / total_l2_norm : 0.0) << ","
           << st.max_abs << "," << st.worst_index << ","
           << st.worst_value.real() << "," << st.worst_value.imag() << ","
           << total_l2_norm << "," << sol.bie_residual_l2 << "\n";
      }
    }
  }
}

void append_radiation_adjoint_column_detail(const std::filesystem::path& path,
                                            double omega,
                                            int dof_col,
                                            const char* solver_name,
                                            const std::vector<Complex>& A_col_major,
                                            const std::vector<Complex>& rhs,
                                            const std::vector<Complex>& u,
                                            const std::vector<bem_frequency_domain::Id>& id_by_index,
                                            const std::unordered_set<networkFace*>& robin_faces,
                                            const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                            const FaceSets& face_sets,
                                            std::size_t top_k = 8) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || A_col_major.size() != n * n || rhs.size() != n || u.size() != n)
    return;

  auto aidx = [&](std::size_t row, std::size_t col) { return row + col * n; };
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto col_kind_name = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (is_robin_id(node, f))
      return "robin_phin";
    if (isDirichletBieDofKey(node, f))
      return "dirichlet_phin";
    return "neumann_phi";
  };
  auto row_kind = [](const BEM_DOF_Base* node, const networkFace* f) {
    return (node && node->BCInterface && isNeumannBieDofKey(node, f)) ? 1 : 0;
  };
  auto row_kind_name = [](int k) { return k == 1 ? "interface_constraint" : "bie"; };
  auto key_kind_name = [](const networkFace* f) { return f ? "face_key" : "default_key"; };
  auto active_neumann_count = [](const BEM_DOF_Base* node) {
    int active = 0;
    if (!node)
      return active;
    for (const auto& [df, d] : node->dofs) {
      if (d.index >= 0 && isNeumannBieDofKey(node, df))
        ++active;
    }
    return active;
  };
  auto surface_counts = [&](const BEM_DOF_Base* node) {
    std::array<std::size_t, 4> out{0, 0, 0, 0}; // float, free, wall, robin
    if (!node)
      return out;
    for (auto* adj : node->getBoundaryFaces()) {
      if (face_sets.float_surface.contains(adj))
        ++out[0];
      else if (face_sets.free_surface.contains(adj))
        ++out[1];
      else
        ++out[2];
      if (robin_faces.contains(adj))
        ++out[3];
    }
    return out;
  };

  struct Component {
    std::size_t index = 0;
    Complex atu{0.0, 0.0};
    Complex residual{0.0, 0.0};
  };
  std::vector<Component> components;
  components.reserve(n);
  for (std::size_t col = 0; col < n; ++col) {
    Complex atu{0.0, 0.0};
    for (std::size_t row = 0; row < n; ++row)
      atu += A_col_major[aidx(row, col)] * u[row];
    components.push_back(Component{col, atu, atu - rhs[col]});
  }
  std::sort(components.begin(), components.end(), [](const Component& a, const Component& b) {
    return std::abs(a.residual) > std::abs(b.residual);
  });
  if (components.size() > top_k)
    components.resize(top_k);

  struct Group {
    std::uint64_t entries = 0;
    Complex sum{0.0, 0.0};
    long double sum_abs = 0.0L;
    double max_abs = 0.0;
    std::size_t worst_row = 0;
    Complex worst_value{0.0, 0.0};
  };
  constexpr int surface_count = 4;
  constexpr int row_kind_count = 2;
  constexpr int key_kind_count = 2;
  auto group_index = [](int surface, int rk, int key_kind) {
    return (surface * row_kind_count + rk) * key_kind_count + key_kind;
  };

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,rank,component_index,component_surface,component_key_kind,component_unknown,"
          "x,y,z,is_bc_interface,active_neumann_dofs,float_faces,free_faces,wall_faces,robin_adjacent_faces,"
          "rhs_re,rhs_im,atu_re,atu_im,res_re,res_im,res_abs,"
          "row_surface,row_kind,row_key_kind,entries,contrib_re,contrib_im,contrib_abs,sum_abs,max_abs,"
          "worst_row,worst_row_surface,worst_row_kind,worst_row_key_kind,worst_re,worst_im\n";
  }

  for (std::size_t rank = 0; rank < components.size(); ++rank) {
    const auto& comp = components[rank];
    const auto* col_node = std::get<0>(id_by_index[comp.index]);
    const auto* col_face = std::get<1>(id_by_index[comp.index]);
    const auto col_x = col_node ? col_node->getPosition() : Tddd{0.0, 0.0, 0.0};
    const auto counts = surface_counts(col_node);
    const int col_surface = radiation_surface_kind(col_node, col_face, face_sets);
    const int col_key_kind = col_face ? 0 : 1;

    std::array<Group, surface_count * row_kind_count * key_kind_count> groups{};
    for (std::size_t row = 0; row < n; ++row) {
      const auto* row_node = std::get<0>(id_by_index[row]);
      const auto* row_face = std::get<1>(id_by_index[row]);
      const int rs = radiation_surface_kind(row_node, row_face, face_sets);
      const int rk = row_kind(row_node, row_face);
      const int kk = row_face ? 0 : 1;
      const Complex value = A_col_major[aidx(row, comp.index)] * u[row];
      const double av = std::abs(value);
      auto& g = groups[static_cast<std::size_t>(group_index(rs, rk, kk))];
      ++g.entries;
      g.sum += value;
      g.sum_abs += av;
      if (av > g.max_abs) {
        g.max_abs = av;
        g.worst_row = row;
        g.worst_value = value;
      }
    }

    for (int rs = 0; rs < surface_count; ++rs) {
      for (int rk = 0; rk < row_kind_count; ++rk) {
        for (int kk = 0; kk < key_kind_count; ++kk) {
          const auto& g = groups[static_cast<std::size_t>(group_index(rs, rk, kk))];
          if (g.entries == 0 || (std::abs(g.sum) == 0.0 && g.sum_abs == 0.0))
            continue;
          const auto* worst_node = std::get<0>(id_by_index[g.worst_row]);
          const auto* worst_face = std::get<1>(id_by_index[g.worst_row]);
          fs << std::scientific << std::setprecision(12)
             << omega << "," << dof_col << "," << solver_name << "," << rank << ","
             << comp.index << ","
             << radiation_surface_name(col_surface) << "," << key_kind_name(col_face) << ","
             << col_kind_name(col_node, col_face) << ","
             << col_x[0] << "," << col_x[1] << "," << col_x[2] << ","
             << (col_node && col_node->BCInterface ? 1 : 0) << ","
             << active_neumann_count(col_node) << ","
             << counts[0] << "," << counts[1] << "," << counts[2] << "," << counts[3] << ","
             << rhs[comp.index].real() << "," << rhs[comp.index].imag() << ","
             << comp.atu.real() << "," << comp.atu.imag() << ","
             << comp.residual.real() << "," << comp.residual.imag() << "," << std::abs(comp.residual) << ","
             << radiation_surface_name(rs) << "," << row_kind_name(rk) << "," << (kk == 0 ? "face_key" : "default_key") << ","
             << g.entries << "," << g.sum.real() << "," << g.sum.imag() << "," << std::abs(g.sum) << ","
             << static_cast<double>(g.sum_abs) << "," << g.max_abs << ","
             << g.worst_row << ","
             << radiation_surface_name(radiation_surface_kind(worst_node, worst_face, face_sets)) << ","
             << row_kind_name(row_kind(worst_node, worst_face)) << "," << key_kind_name(worst_face) << ","
             << g.worst_value.real() << "," << g.worst_value.imag() << "\n";
        }
      }
    }
  }
}

void write_radiation_time_domain_geometry_compare(const std::filesystem::path& path,
                                                  const std::unordered_set<networkFace*>& float_faces,
                                                  const Tddd& ref_point,
                                                  double moment_sign) {
  std::ofstream fs(path);
  fs << "face_index,cx,cy,cz,"
        "area_face,area_time_domain,area_rel_diff,"
        "nx_face,ny_face,nz_face,nx_time_domain,ny_time_domain,nz_time_domain,"
        "normal_dot,normal_angle_rad,normal_diff_abs";
  for (int dof = 0; dof < 6; ++dof)
    fs << ",mode" << dof << "_face,mode" << dof << "_time_domain,mode" << dof << "_abs_diff";
  fs << "\n";

  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};

  for (auto* f : float_faces) {
    auto [p0, p1, p2] = f->getPoints();
    const T3Tddd X012{p0->X, p1->X, p2->X};
    const auto intpX = interpolationTriangleLinear0101(X012);
    const Tddd n_td = TriangleNormal(X012);

    double area_td = 0.0;
    for (const auto& [x0, x1, w0w1] : __GWGW10__Tuple)
      area_td += intpX.J(x0, x1) * w0w1;

    std::array<double, 6> mode_face{};
    std::array<double, 6> mode_td{};
    mode_face.fill(0.0);
    mode_td.fill(0.0);

    for (const auto& l : bary) {
      const double l0 = l[0], l1 = l[1], l2 = l[2];
      const Tddd xq = l0 * p0->X + l1 * p1->X + l2 * p2->X;
      for (int dof = 0; dof < 6; ++dof)
        mode_face[static_cast<std::size_t>(dof)] += w * rigid_mode_normal_component(dof, xq, f->normal, ref_point, moment_sign) * f->area;
    }

    for (const auto& [x0, x1, w0w1] : __GWGW10__Tuple) {
      const Tddd xq = intpX(x0, x1);
      const double jac_w = intpX.J(x0, x1) * w0w1;
      for (int dof = 0; dof < 6; ++dof)
        mode_td[static_cast<std::size_t>(dof)] += rigid_mode_normal_component(dof, xq, n_td, ref_point, moment_sign) * jac_w;
    }

    const double area_rel = std::abs(f->area - area_td) / std::max({std::abs(f->area), std::abs(area_td), 1e-300});
    const double normal_dot = std::clamp(Dot(f->normal, n_td), -1.0, 1.0);
    const double normal_angle = std::acos(normal_dot);
    const double normal_diff = Norm(f->normal - n_td);

    fs << std::scientific << std::setprecision(12)
       << f->index << ","
       << f->centroid[0] << "," << f->centroid[1] << "," << f->centroid[2] << ","
       << f->area << "," << area_td << "," << area_rel << ","
       << f->normal[0] << "," << f->normal[1] << "," << f->normal[2] << ","
       << n_td[0] << "," << n_td[1] << "," << n_td[2] << ","
       << normal_dot << "," << normal_angle << "," << normal_diff;
    for (int dof = 0; dof < 6; ++dof) {
      const double a = mode_face[static_cast<std::size_t>(dof)];
      const double b = mode_td[static_cast<std::size_t>(dof)];
      fs << "," << a << "," << b << "," << std::abs(a - b);
    }
    fs << "\n";
  }
}

std::vector<double> lumped_area_for_ids(const std::vector<bem_frequency_domain::Id>& id_by_index,
                                        const std::unordered_set<networkFace*>& faces);

bool faces_share_point(const networkFace* a, const networkFace* b);
int radiation_surface_kind(const BEM_DOF_Base* node, const networkFace* f, const FaceSets& face_sets);
const char* radiation_surface_name(int k);

void append_bie_operator_symmetry(const std::filesystem::path& path,
                                  double omega,
                                  const char* solver_name,
                                  const BEM_BVP& bvp,
                                  const std::vector<bem_frequency_domain::Id>& id_by_index,
                                  const std::unordered_set<networkFace*>& faces) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || bvp.IGIGn.size() != n)
    return;

  std::unordered_map<bem_frequency_domain::Id, std::size_t, bem_frequency_domain::IdHash, bem_frequency_domain::IdEq> index;
  index.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    index.emplace(id_by_index[i], i);

  std::vector<double> lumped_area(n, 0.0);
  auto add_lumped_area = [&](BEM_DOF_Base* node, networkFace* f, double value) {
    bem_frequency_domain::Id key{node, f};
    auto it = index.find(key);
    if (it == index.end()) {
      key = {node, nullptr};
      it = index.find(key);
    }
    if (it != index.end())
      lumped_area[it->second] += value;
  };

  for (auto* f : faces) {
    if (!f)
      continue;
    auto [p0, p1, p2] = f->getPoints();
    const double a3 = f->area / 3.0;
    add_lumped_area(p0, f, a3);
    add_lumped_area(p1, f, a3);
    add_lumped_area(p2, f, a3);
  }

  struct Stat {
    std::uint64_t pairs = 0;
    double max_abs = 0.0;
    double max_rel = 0.0;
    std::size_t wi = 0;
    std::size_t wj = 0;
    double vij = 0.0;
    double vji = 0.0;
    double ai = 0.0;
    double aj = 0.0;
  };

  auto update = [](Stat& st, std::size_t i, std::size_t j, double a, double b, double wi, double wj) {
    ++st.pairs;
    const double diff = std::abs(a - b);
    const double denom = std::max({std::abs(a), std::abs(b), 1e-300});
    const double rel = diff / denom;
    if (rel > st.max_rel || (rel == st.max_rel && diff > st.max_abs)) {
      st.max_abs = diff;
      st.max_rel = rel;
      st.wi = i;
      st.wj = j;
      st.vij = a;
      st.vji = b;
      st.ai = wi;
      st.aj = wj;
    }
  };

  Stat ig, ign, ig_weighted, ign_weighted;
  for (std::size_t i = 0; i < n; ++i) {
    if (bvp.IGIGn[i].size() != n)
      return;
    for (std::size_t j = i + 1; j < n; ++j) {
      const double wi = lumped_area[i];
      const double wj = lumped_area[j];
      const double gij = bvp.IGIGn[i][j][0];
      const double gji = bvp.IGIGn[j][i][0];
      const double hij = bvp.IGIGn[i][j][1];
      const double hji = bvp.IGIGn[j][i][1];
      update(ig, i, j, gij, gji, wi, wj);
      update(ign, i, j, hij, hji, wi, wj);
      if (wi > 0.0 && wj > 0.0) {
        update(ig_weighted, i, j, wi * gij, wj * gji, wi, wj);
        update(ign_weighted, i, j, wi * hij, wj * hji, wi, wj);
      }
    }
  }

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header)
    fs << "omega,solver,metric,pairs,max_abs,max_rel,worst_i,worst_j,value_ij,value_ji,"
          "lumped_area_i,lumped_area_j\n";

  auto write = [&](const char* metric, const Stat& st) {
    fs << std::scientific << std::setprecision(12)
       << omega << "," << solver_name << "," << metric << ","
       << st.pairs << "," << st.max_abs << "," << st.max_rel << ","
       << st.wi << "," << st.wj << "," << st.vij << "," << st.vji << ","
       << st.ai << "," << st.aj << "\n";
  };
  write("IG", ig);
  write("IG_lumped_area_weighted", ig_weighted);
  write("IGn_source_normal_self_comparison_not_adjoint", ign);
  write("IGn_source_normal_self_lumped_not_adjoint", ign_weighted);
}

void append_bie_operator_bins(const std::filesystem::path& path,
                              double omega,
                              const char* solver_name,
                              const BEM_BVP& bvp,
                              const std::vector<bem_frequency_domain::Id>& id_by_index,
                              const std::unordered_set<networkFace*>& faces,
                              const FaceSets& face_sets) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || bvp.IGIGn.size() != n)
    return;
  const auto area = lumped_area_for_ids(id_by_index, faces);

  auto surface_kind = [&](const BEM_DOF_Base* node, const networkFace* f) -> int {
    if (f) {
      if (face_sets.float_surface.contains(const_cast<networkFace*>(f)))
        return 0;
      if (face_sets.free_surface.contains(const_cast<networkFace*>(f)))
        return 1;
      return 2;
    }
    if (!node)
      return 3;
    bool has_float = false;
    bool has_free = false;
    bool has_wall = false;
    for (auto* adj : node->getBoundaryFaces()) {
      if (face_sets.float_surface.contains(adj))
        has_float = true;
      else if (face_sets.free_surface.contains(adj))
        has_free = true;
      else
        has_wall = true;
    }
    const int count = static_cast<int>(has_float) + static_cast<int>(has_free) + static_cast<int>(has_wall);
    if (count != 1)
      return 3;
    if (has_float)
      return 0;
    if (has_free)
      return 1;
    return 2;
  };
  auto surface_name = [](int k) -> const char* {
    switch (k) {
    case 0:
      return "float";
    case 1:
      return "free";
    case 2:
      return "wall";
    default:
      return "mixed";
    }
  };

  auto distance_bin = [&](std::size_t i, std::size_t j, double distance) -> int {
    const auto* ni = std::get<0>(id_by_index[i]);
    const auto* nj = std::get<0>(id_by_index[j]);
    const auto* fi = std::get<1>(id_by_index[i]);
    const auto* fj = std::get<1>(id_by_index[j]);
    if (ni && ni == nj)
      return 0;
    if (fi && fj && faces_share_point(fi, fj))
      return 1;
    const double h = std::max(std::sqrt(std::max(area[i], 0.0)), std::sqrt(std::max(area[j], 0.0)));
    if (h <= 0.0)
      return 4;
    return (distance <= 5.0 * h) ? 2 : 3;
  };
  auto bin_name = [](int k) -> const char* {
    switch (k) {
    case 0:
      return "same_node";
    case 1:
      return "adjacent_face";
    case 2:
      return "near";
    case 3:
      return "far";
    default:
      return "zero_area";
    }
  };

  struct Stat {
    std::uint64_t pairs = 0;
    long double rel2 = 0.0L;
    double max_abs = 0.0;
    double max_rel = 0.0;
    std::size_t wi = 0;
    std::size_t wj = 0;
    double vij = 0.0;
    double vji = 0.0;
    double ai = 0.0;
    double aj = 0.0;
    double distance = 0.0;
  };
  constexpr int surface_count = 4;
  constexpr int bin_count = 5;
  std::array<std::array<Stat, bin_count>, surface_count * surface_count> ig{};
  std::array<std::array<Stat, bin_count>, surface_count * surface_count> ign{};

  auto update = [](Stat& st, std::size_t i, std::size_t j, double a, double b, double ai, double aj, double distance) {
    ++st.pairs;
    const double diff = std::abs(a - b);
    const double denom = std::max({std::abs(a), std::abs(b), 1e-300});
    const double rel = diff / denom;
    st.rel2 += static_cast<long double>(rel) * rel;
    if (rel > st.max_rel || (rel == st.max_rel && diff > st.max_abs)) {
      st.max_abs = diff;
      st.max_rel = rel;
      st.wi = i;
      st.wj = j;
      st.vij = a;
      st.vji = b;
      st.ai = ai;
      st.aj = aj;
      st.distance = distance;
    }
  };

  for (std::size_t i = 0; i < n; ++i) {
    if (area[i] <= 0.0 || bvp.IGIGn[i].size() != n)
      continue;
    const auto* ni = std::get<0>(id_by_index[i]);
    const auto* fi = std::get<1>(id_by_index[i]);
    const int si = surface_kind(ni, fi);
    const Tddd xi = ni ? ni->getPosition() : Tddd{0.0, 0.0, 0.0};
    for (std::size_t j = i + 1; j < n; ++j) {
      if (area[j] <= 0.0 || bvp.IGIGn[j].size() != n)
        continue;
      const auto* nj = std::get<0>(id_by_index[j]);
      const auto* fj = std::get<1>(id_by_index[j]);
      const int sj = surface_kind(nj, fj);
      const Tddd xj = nj ? nj->getPosition() : Tddd{0.0, 0.0, 0.0};
      const double distance = Norm(xj - xi);
      const int db = distance_bin(i, j, distance);
      const int key = si * surface_count + sj;
      update(ig[static_cast<std::size_t>(key)][static_cast<std::size_t>(db)],
             i,
             j,
             area[i] * bvp.IGIGn[i][j][0],
             area[j] * bvp.IGIGn[j][i][0],
             area[i],
             area[j],
             distance);
      update(ign[static_cast<std::size_t>(key)][static_cast<std::size_t>(db)],
             i,
             j,
             area[i] * bvp.IGIGn[i][j][1],
             area[j] * bvp.IGIGn[j][i][1],
             area[i],
             area[j],
             distance);
    }
  }

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header)
    fs << "omega,solver,metric,row_surface,col_surface,distance_bin,pairs,max_abs,max_rel,rms_rel,"
          "worst_i,worst_j,value_ij,value_ji,lumped_area_i,lumped_area_j,distance\n";

  auto write = [&](const char* metric, int si, int sj, int db, const Stat& st) {
    if (st.pairs == 0)
      return;
    const double rms_rel = std::sqrt(static_cast<double>(st.rel2 / static_cast<long double>(st.pairs)));
    fs << std::scientific << std::setprecision(12)
       << omega << "," << solver_name << "," << metric << ","
       << surface_name(si) << "," << surface_name(sj) << "," << bin_name(db) << ","
       << st.pairs << "," << st.max_abs << "," << st.max_rel << "," << rms_rel << ","
       << st.wi << "," << st.wj << "," << st.vij << "," << st.vji << ","
       << st.ai << "," << st.aj << "," << st.distance << "\n";
  };
  for (int si = 0; si < surface_count; ++si) {
    for (int sj = 0; sj < surface_count; ++sj) {
      const int key = si * surface_count + sj;
      for (int db = 0; db < bin_count; ++db) {
        write("IG_lumped_area_weighted", si, sj, db, ig[static_cast<std::size_t>(key)][static_cast<std::size_t>(db)]);
        write("IGn_lumped_area_weighted", si, sj, db, ign[static_cast<std::size_t>(key)][static_cast<std::size_t>(db)]);
      }
    }
  }
}

int radiation_surface_kind(const BEM_DOF_Base* node, const networkFace* f, const FaceSets& face_sets) {
  if (f) {
    if (face_sets.float_surface.contains(const_cast<networkFace*>(f)))
      return 0;
    if (face_sets.free_surface.contains(const_cast<networkFace*>(f)))
      return 1;
    return 2;
  }
  if (!node)
    return 3;
  bool has_float = false;
  bool has_free = false;
  bool has_wall = false;
  for (auto* adj : node->getBoundaryFaces()) {
    if (face_sets.float_surface.contains(adj))
      has_float = true;
    else if (face_sets.free_surface.contains(adj))
      has_free = true;
    else
      has_wall = true;
  }
  const int count = static_cast<int>(has_float) + static_cast<int>(has_free) + static_cast<int>(has_wall);
  if (count != 1)
    return 3;
  if (has_float)
    return 0;
  if (has_free)
    return 1;
  return 2;
}

const char* radiation_surface_name(int k) {
  switch (k) {
  case 0:
    return "float";
  case 1:
    return "free";
  case 2:
    return "wall";
  default:
    return "mixed";
  }
}

void append_radiation_linear_system_blocks(const std::filesystem::path& path,
                                           double omega,
                                           int dof_col,
                                           const char* solver_name,
                                           const std::vector<Complex>& A_col_major,
                                           const std::vector<Complex>& b,
                                           const std::vector<bem_frequency_domain::Id>& id_by_index,
                                           const std::unordered_set<networkFace*>& robin_faces,
                                           const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                           const FaceSets& face_sets,
                                           bool robin_unknown_phi,
                                           double kappa_scale) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || A_col_major.size() != n * n || b.size() != n)
    return;

  struct Stat {
    std::uint64_t entries = 0;
    long double sum_abs = 0.0L;
    long double frob2 = 0.0L;
    double max_abs = 0.0;
    std::size_t worst_row = 0;
    std::size_t worst_col = 0;
    Complex worst_value{0.0, 0.0};
  };
  struct RhsStat {
    std::uint64_t rows = 0;
    long double sum_abs = 0.0L;
    long double frob2 = 0.0L;
    double max_abs = 0.0;
  };

  constexpr int surface_count = 4;
  constexpr int row_kind_count = 2; // BIE row, interface constraint row
  constexpr int col_kind_count = 4; // neumann_phi, dirichlet_phin, robin_phin, robin_phi
  std::array<Stat, surface_count * surface_count * row_kind_count * col_kind_count> stats{};
  std::array<RhsStat, surface_count * row_kind_count> rhs_stats{};

  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto col_kind = [&](const BEM_DOF_Base* node, const networkFace* f) {
    const bool robin = is_robin_id(node, f);
    if (robin)
      return robin_unknown_phi ? 3 : 2;
    if (isDirichletBieDofKey(node, f))
      return 1;
    return 0;
  };
  auto col_kind_name = [](int k) {
    switch (k) {
    case 0:
      return "neumann_phi";
    case 1:
      return "dirichlet_phin";
    case 2:
      return "robin_phin";
    default:
      return "robin_phi";
    }
  };
  auto row_kind = [](const BEM_DOF_Base* node, const networkFace* f) {
    return (node && node->BCInterface && isNeumannBieDofKey(node, f)) ? 1 : 0;
  };
  auto row_kind_name = [](int k) { return k == 1 ? "interface_constraint" : "bie"; };
  auto stat_index = [](int rs, int cs, int rk, int ck) {
    return (((rs * surface_count + cs) * row_kind_count + rk) * col_kind_count + ck);
  };
  auto rhs_index = [](int rs, int rk) { return rs * row_kind_count + rk; };
  auto mid = [&](std::size_t row, std::size_t col) { return row + col * n; };

  for (std::size_t i = 0; i < n; ++i) {
    const auto* ni = std::get<0>(id_by_index[i]);
    const auto* fi = std::get<1>(id_by_index[i]);
    const int rs = radiation_surface_kind(ni, fi, face_sets);
    const int rk = row_kind(ni, fi);
    auto& rhs = rhs_stats[static_cast<std::size_t>(rhs_index(rs, rk))];
    const double rb = std::abs(b[i]);
    ++rhs.rows;
    rhs.sum_abs += rb;
    rhs.frob2 += static_cast<long double>(rb) * rb;
    rhs.max_abs = std::max(rhs.max_abs, rb);
    for (std::size_t j = 0; j < n; ++j) {
      const auto* nj = std::get<0>(id_by_index[j]);
      const auto* fj = std::get<1>(id_by_index[j]);
      const int cs = radiation_surface_kind(nj, fj, face_sets);
      const int ck = col_kind(nj, fj);
      auto& st = stats[static_cast<std::size_t>(stat_index(rs, cs, rk, ck))];
      const Complex v = A_col_major[mid(i, j)];
      const double av = std::abs(v);
      ++st.entries;
      st.sum_abs += av;
      st.frob2 += static_cast<long double>(av) * av;
      if (av > st.max_abs) {
        st.max_abs = av;
        st.worst_row = i;
        st.worst_col = j;
        st.worst_value = v;
      }
    }
  }

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,row_surface,col_surface,row_kind,col_unknown,"
          "entries,sum_abs,frob_norm,max_abs,worst_row,worst_col,worst_re,worst_im,"
          "rhs_rows,sum_rhs_abs,rhs_frob_norm,rhs_max_abs,robin_unknown_phi,kappa_scale\n";
  }
  for (int rs = 0; rs < surface_count; ++rs) {
    for (int cs = 0; cs < surface_count; ++cs) {
      for (int rk = 0; rk < row_kind_count; ++rk) {
        const auto& rhs = rhs_stats[static_cast<std::size_t>(rhs_index(rs, rk))];
        for (int ck = 0; ck < col_kind_count; ++ck) {
          const auto& st = stats[static_cast<std::size_t>(stat_index(rs, cs, rk, ck))];
          if (st.entries == 0)
            continue;
          fs << std::scientific << std::setprecision(12)
             << omega << "," << dof_col << "," << solver_name << ","
             << radiation_surface_name(rs) << "," << radiation_surface_name(cs) << ","
             << row_kind_name(rk) << "," << col_kind_name(ck) << ","
             << st.entries << ","
             << static_cast<double>(st.sum_abs) << ","
             << std::sqrt(static_cast<double>(st.frob2)) << ","
             << st.max_abs << "," << st.worst_row << "," << st.worst_col << ","
             << st.worst_value.real() << "," << st.worst_value.imag() << ","
             << rhs.rows << ","
             << static_cast<double>(rhs.sum_abs) << ","
             << std::sqrt(static_cast<double>(rhs.frob2)) << ","
             << rhs.max_abs << ","
             << (robin_unknown_phi ? 1 : 0) << "," << kappa_scale << "\n";
        }
      }
    }
  }
}

void append_radiation_linear_system_id_map(const std::filesystem::path& path,
                                           double omega,
                                           int dof_col,
                                           const char* solver_name,
                                           const std::vector<bem_frequency_domain::Id>& id_by_index,
                                           const std::unordered_set<networkFace*>& robin_faces,
                                           const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                           const FaceSets& face_sets,
                                           bool robin_unknown_phi,
                                           double kappa_scale) {
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto col_kind = [&](const BEM_DOF_Base* node, const networkFace* f) {
    const bool robin = is_robin_id(node, f);
    if (robin)
      return robin_unknown_phi ? 3 : 2;
    if (isDirichletBieDofKey(node, f))
      return 1;
    return 0;
  };
  auto col_kind_name = [](int k) {
    switch (k) {
    case 0:
      return "neumann_phi";
    case 1:
      return "dirichlet_phin";
    case 2:
      return "robin_phin";
    default:
      return "robin_phi";
    }
  };
  auto row_kind = [](const BEM_DOF_Base* node, const networkFace* f) {
    return (node && node->BCInterface && isNeumannBieDofKey(node, f)) ? 1 : 0;
  };
  auto row_kind_name = [](int k) { return k == 1 ? "interface_constraint" : "bie"; };

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,index,node_ptr,face_ptr,x,y,z,face_cx,face_cy,face_cz,"
          "face_nx,face_ny,face_nz,face_area,surface,row_kind,col_unknown,"
          "is_bc_interface,is_neumann_key,is_dirichlet_key,is_robin_id,"
          "boundary_faces,float_faces,free_faces,wall_faces,robin_adjacent_faces,"
          "neumann_adjacent_faces,dirichlet_adjacent_faces,active_neumann_dofs,active_dirichlet_index,"
          "robin_unknown_phi,kappa_scale\n";
  }

  for (std::size_t i = 0; i < id_by_index.size(); ++i) {
    const auto* node = std::get<0>(id_by_index[i]);
    const auto* f = std::get<1>(id_by_index[i]);
    if (!node)
      continue;
    std::size_t n_float = 0;
    std::size_t n_free = 0;
    std::size_t n_wall = 0;
    std::size_t n_robin = 0;
    std::size_t n_neumann = 0;
    std::size_t n_dirichlet = 0;
    for (auto* adj : node->getBoundaryFaces()) {
      if (face_sets.float_surface.contains(adj))
        ++n_float;
      else if (face_sets.free_surface.contains(adj))
        ++n_free;
      else
        ++n_wall;
      if (robin_faces.contains(adj))
        ++n_robin;
      if (isNeumannBoundaryState(node, adj))
        ++n_neumann;
      if (isDirichletBoundaryState(node, adj))
        ++n_dirichlet;
    }
    int active_neumann = 0;
    for (const auto& [df, d] : node->dofs) {
      if (d.index >= 0 && isNeumannBieDofKey(node, df))
        ++active_neumann;
    }
    const auto* active_dirichlet = node->findActiveBieDof(nullptr);
    const auto x = node->getPosition();
    const Tddd fc = f ? f->center : Tddd{0.0, 0.0, 0.0};
    const Tddd fn = f ? f->normal : Tddd{0.0, 0.0, 0.0};
    const double fa = f ? f->area : 0.0;

    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << "," << i << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << reinterpret_cast<std::uintptr_t>(f) << ","
       << x[0] << "," << x[1] << "," << x[2] << ","
       << fc[0] << "," << fc[1] << "," << fc[2] << ","
       << fn[0] << "," << fn[1] << "," << fn[2] << "," << fa << ","
       << radiation_surface_name(radiation_surface_kind(node, f, face_sets)) << ","
       << row_kind_name(row_kind(node, f)) << ","
       << col_kind_name(col_kind(node, f)) << ","
       << (node->BCInterface ? 1 : 0) << ","
       << (isNeumannBieDofKey(node, f) ? 1 : 0) << ","
       << (isDirichletBieDofKey(node, f) ? 1 : 0) << ","
       << (is_robin_id(node, f) ? 1 : 0) << ","
       << node->getBoundaryFaces().size() << ","
       << n_float << "," << n_free << "," << n_wall << "," << n_robin << ","
       << n_neumann << "," << n_dirichlet << ","
       << active_neumann << ","
       << (active_dirichlet ? active_dirichlet->index : -1) << ","
       << (robin_unknown_phi ? 1 : 0) << "," << kappa_scale << "\n";
  }
}

void append_radiation_interface_constraints(const std::filesystem::path& path,
                                            double omega,
                                            int dof_col,
                                            const char* solver_name,
                                            const std::vector<Complex>& A_col_major,
                                            const std::vector<Complex>& b,
                                            const std::vector<bem_frequency_domain::Id>& id_by_index,
                                            const std::unordered_set<networkFace*>& robin_faces,
                                            const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                            const FaceSets& face_sets,
                                            double gravity,
                                            bool robin_unknown_phi,
                                            double kappa_scale) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || A_col_major.size() != n * n || b.size() != n)
    return;

  auto aidx = [&](std::size_t row, std::size_t col) { return row + col * n; };
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto surface_counts = [&](const BEM_DOF_Base* node) {
    std::array<std::size_t, 4> out{0, 0, 0, 0}; // float, free, wall, robin
    if (!node)
      return out;
    for (auto* adj : node->getBoundaryFaces()) {
      if (face_sets.float_surface.contains(adj))
        ++out[0];
      else if (face_sets.free_surface.contains(adj))
        ++out[1];
      else
        ++out[2];
      if (robin_faces.contains(adj))
        ++out[3];
    }
    return out;
  };
  auto active_neumann_count = [](const BEM_DOF_Base* node) {
    int active = 0;
    if (!node)
      return active;
    for (const auto& [df, d] : node->dofs) {
      if (d.index >= 0 && isNeumannBieDofKey(node, df))
        ++active;
    }
    return active;
  };

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,row_index,node_ptr,row_face_ptr,x,y,z,row_surface,"
          "active_dirichlet_index,active_dirichlet_surface,active_neumann_dofs,"
          "float_faces,free_faces,wall_faces,robin_adjacent_faces,is_robin_node,is_robin_row,"
          "self_re,self_im,self_abs,dir_re,dir_im,dir_abs,rhs_re,rhs_im,rhs_abs,"
          "kappa_re,kappa_im,kappa_abs,expected_dir_re,expected_dir_im,expected_mismatch_abs,"
          "dir_over_self_abs,row_nonzero_count,row_max_abs,row_sum_abs,robin_unknown_phi,kappa_scale\n";
  }

  for (std::size_t i = 0; i < n; ++i) {
    const auto* node = std::get<0>(id_by_index[i]);
    const auto* f_row = std::get<1>(id_by_index[i]);
    if (!node || !(node->BCInterface && isNeumannBieDofKey(node, f_row)))
      continue;
    const auto* d_dir = node->findActiveBieDof(nullptr);
    if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
      continue;
    const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);
    const auto* f_dir = std::get<1>(id_by_index[j_dir]);
    const bool robin_node = [&] {
      const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
      return it != node_is_robin.end() && it->second;
    }();
    const bool robin_row = is_robin_id(node, f_row);
    Complex kappa{0.0, 0.0};
    if (robin_node) {
      const Complex s(0.0, -omega);
      kappa = kappa_scale * (-(s * s) / gravity);
    }
    const Complex self = A_col_major[aidx(i, i)];
    const Complex dir = A_col_major[aidx(i, j_dir)];
    const Complex expected_dir = (std::abs(kappa) > 1e-15) ? -self / kappa : Complex{0.0, 0.0};
    const auto counts = surface_counts(node);
    const auto x = node->getPosition();
    std::size_t row_nonzero = 0;
    double row_max = 0.0;
    long double row_sum = 0.0L;
    for (std::size_t j = 0; j < n; ++j) {
      const double av = std::abs(A_col_major[aidx(i, j)]);
      if (av > 0.0) {
        ++row_nonzero;
        row_sum += av;
        row_max = std::max(row_max, av);
      }
    }
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << "," << i << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << reinterpret_cast<std::uintptr_t>(f_row) << ","
       << x[0] << "," << x[1] << "," << x[2] << ","
       << radiation_surface_name(radiation_surface_kind(node, f_row, face_sets)) << ","
       << j_dir << ","
       << radiation_surface_name(radiation_surface_kind(node, f_dir, face_sets)) << ","
       << active_neumann_count(node) << ","
       << counts[0] << "," << counts[1] << "," << counts[2] << "," << counts[3] << ","
       << (robin_node ? 1 : 0) << "," << (robin_row ? 1 : 0) << ","
       << self.real() << "," << self.imag() << "," << std::abs(self) << ","
       << dir.real() << "," << dir.imag() << "," << std::abs(dir) << ","
       << b[i].real() << "," << b[i].imag() << "," << std::abs(b[i]) << ","
       << kappa.real() << "," << kappa.imag() << "," << std::abs(kappa) << ","
       << expected_dir.real() << "," << expected_dir.imag() << ","
       << std::abs(dir - expected_dir) << ","
       << (std::abs(self) > 0.0 ? std::abs(dir) / std::abs(self) : 0.0) << ","
       << row_nonzero << "," << row_max << "," << static_cast<double>(row_sum) << ","
       << (robin_unknown_phi ? 1 : 0) << "," << kappa_scale << "\n";
  }
}

void append_radiation_interface_adjoint_pairs(const std::filesystem::path& path,
                                              double omega,
                                              int dof_col,
                                              const char* solver_name,
                                              const std::vector<Complex>& A_col_major,
                                              const std::vector<Complex>& rhs,
                                              const std::vector<Complex>& u,
                                              const std::vector<bem_frequency_domain::Id>& id_by_index,
                                              const std::unordered_set<networkFace*>& robin_faces,
                                              const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                              const FaceSets& face_sets,
                                              double gravity,
                                              bool robin_unknown_phi,
                                              double kappa_scale) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || A_col_major.size() != n * n || rhs.size() != n || u.size() != n)
    return;

  auto aidx = [&](std::size_t row, std::size_t col) { return row + col * n; };
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto col_kind_name = [&](const BEM_DOF_Base* node, const networkFace* f) {
    const bool robin = is_robin_id(node, f);
    if (robin)
      return robin_unknown_phi ? "robin_phi" : "robin_phin";
    if (isDirichletBieDofKey(node, f))
      return "dirichlet_phin";
    return "neumann_phi";
  };
  auto active_neumann_count = [](const BEM_DOF_Base* node) {
    int active = 0;
    if (!node)
      return active;
    for (const auto& [df, d] : node->dofs) {
      if (d.index >= 0 && isNeumannBieDofKey(node, df))
        ++active;
    }
    return active;
  };

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,row_index,dir_index,node_ptr,row_face_ptr,x,y,z,"
          "row_surface,dir_surface,row_unknown,dir_unknown,active_neumann_dofs,"
          "self_re,self_im,dir_re,dir_im,kappa_re,kappa_im,"
          "u_row_re,u_row_im,u_dir_re,u_dir_im,rhs_row_re,rhs_row_im,rhs_dir_re,rhs_dir_im,"
          "row_res_re,row_res_im,row_res_abs,"
          "dir_adjoint_res_re,dir_adjoint_res_im,dir_adjoint_res_abs,"
          "row_to_dir_adjoint_re,row_to_dir_adjoint_im,row_to_dir_adjoint_abs,"
          "row_to_self_adjoint_re,row_to_self_adjoint_im,row_to_self_adjoint_abs,"
          "robin_unknown_phi,kappa_scale\n";
  }

  for (std::size_t i = 0; i < n; ++i) {
    const auto* node = std::get<0>(id_by_index[i]);
    const auto* f_row = std::get<1>(id_by_index[i]);
    if (!node || !(node->BCInterface && isNeumannBieDofKey(node, f_row)))
      continue;
    const auto* d_dir = node->findActiveBieDof(nullptr);
    if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
      continue;
    const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);
    const auto* f_dir = std::get<1>(id_by_index[j_dir]);

    Complex row_res{0.0, 0.0};
    Complex dir_adjoint_res{0.0, 0.0};
    for (std::size_t j = 0; j < n; ++j)
      row_res += A_col_major[aidx(i, j)] * u[j];
    row_res -= rhs[i];
    for (std::size_t r = 0; r < n; ++r)
      dir_adjoint_res += A_col_major[aidx(r, j_dir)] * u[r];
    dir_adjoint_res -= rhs[j_dir];

    const bool robin_node = [&] {
      const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
      return it != node_is_robin.end() && it->second;
    }();
    Complex kappa{0.0, 0.0};
    if (robin_node) {
      const Complex s(0.0, -omega);
      kappa = kappa_scale * (-(s * s) / gravity);
    }

    const Complex self = A_col_major[aidx(i, i)];
    const Complex dir = A_col_major[aidx(i, j_dir)];
    const Complex row_to_dir = dir * u[i];
    const Complex row_to_self = self * u[i];
    const auto x = node->getPosition();
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << ","
       << i << "," << j_dir << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << reinterpret_cast<std::uintptr_t>(f_row) << ","
       << x[0] << "," << x[1] << "," << x[2] << ","
       << radiation_surface_name(radiation_surface_kind(node, f_row, face_sets)) << ","
       << radiation_surface_name(radiation_surface_kind(node, f_dir, face_sets)) << ","
       << col_kind_name(node, f_row) << "," << col_kind_name(node, f_dir) << ","
       << active_neumann_count(node) << ","
       << self.real() << "," << self.imag() << ","
       << dir.real() << "," << dir.imag() << ","
       << kappa.real() << "," << kappa.imag() << ","
       << u[i].real() << "," << u[i].imag() << ","
       << u[j_dir].real() << "," << u[j_dir].imag() << ","
       << rhs[i].real() << "," << rhs[i].imag() << ","
       << rhs[j_dir].real() << "," << rhs[j_dir].imag() << ","
       << row_res.real() << "," << row_res.imag() << "," << std::abs(row_res) << ","
       << dir_adjoint_res.real() << "," << dir_adjoint_res.imag() << "," << std::abs(dir_adjoint_res) << ","
       << row_to_dir.real() << "," << row_to_dir.imag() << "," << std::abs(row_to_dir) << ","
       << row_to_self.real() << "," << row_to_self.imag() << "," << std::abs(row_to_self) << ","
       << (robin_unknown_phi ? 1 : 0) << "," << kappa_scale << "\n";
  }
}

void append_radiation_interface_constraint_balance(const std::filesystem::path& path,
                                                   double omega,
                                                   int dof_col,
                                                   const char* solver_name,
                                                   const std::vector<Complex>& A_col_major,
                                                   const std::vector<Complex>& rhs,
                                                   const std::vector<Complex>& u,
                                                   const std::vector<bem_frequency_domain::Id>& id_by_index,
                                                   const std::unordered_set<networkFace*>& robin_faces,
                                                   const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                                   const FaceSets& face_sets,
                                                   bool robin_unknown_phi,
                                                   const std::string& interface_row_scale_mode,
                                                   bool interface_robin_adjoint_probe) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || A_col_major.size() != n * n || rhs.size() != n || u.size() != n)
    return;

  auto aidx = [&](std::size_t row, std::size_t col) { return row + col * n; };
  auto key_kind_name = [](const networkFace* f) { return f ? "face_key" : "default_key"; };
  auto row_kind = [](const BEM_DOF_Base* node, const networkFace* f) {
    return (node && node->BCInterface && isNeumannBieDofKey(node, f)) ? 1 : 0;
  };
  auto row_kind_name = [](int k) { return k == 1 ? "interface_constraint" : "bie"; };
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto unknown_name = [&](const BEM_DOF_Base* node, const networkFace* f) {
    const bool robin = is_robin_id(node, f);
    if (robin)
      return robin_unknown_phi ? "robin_phi" : "robin_phin";
    if (isDirichletBieDofKey(node, f))
      return "dirichlet_phin";
    return "neumann_phi";
  };
  auto column_adjoint_residual = [&](std::size_t col) {
    Complex out{0.0, 0.0};
    for (std::size_t row = 0; row < n; ++row)
      out += A_col_major[aidx(row, col)] * u[row];
    return out - rhs[col];
  };

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,row_index,dir_index,node_ptr,row_face_ptr,x,y,z,"
          "row_surface,row_kind,row_key_kind,row_unknown,dir_surface,dir_key_kind,dir_unknown,"
          "row_norm,row_sum_abs,row_rhs_re,row_rhs_im,row_rhs_abs,row_solution_re,row_solution_im,"
          "row_res_re,row_res_im,row_res_abs,row_adjoint_res_re,row_adjoint_res_im,row_adjoint_res_abs,"
          "dir_adjoint_res_re,dir_adjoint_res_im,dir_adjoint_res_abs,"
          "max_contrib_col,max_contrib_surface,max_contrib_key_kind,max_contrib_unknown,"
          "max_contrib_re,max_contrib_im,max_contrib_abs,"
          "interface_row_scale_mode,interface_robin_adjoint_probe\n";
  }

  for (std::size_t i = 0; i < n; ++i) {
    const auto* node = std::get<0>(id_by_index[i]);
    const auto* f_row = std::get<1>(id_by_index[i]);
    if (!node || row_kind(node, f_row) != 1)
      continue;
    const auto* d_dir = node->findActiveBieDof(nullptr);
    if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
      continue;
    const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);
    const auto* f_dir = std::get<1>(id_by_index[j_dir]);

    long double row_l2 = 0.0L;
    long double row_sum = 0.0L;
    Complex row_res{0.0, 0.0};
    std::size_t max_col = 0;
    Complex max_contrib{0.0, 0.0};
    double max_contrib_abs = -1.0;
    for (std::size_t j = 0; j < n; ++j) {
      const Complex aij = A_col_major[aidx(i, j)];
      const double av = std::abs(aij);
      row_l2 += static_cast<long double>(av) * av;
      row_sum += av;
      const Complex contrib = aij * u[j];
      row_res += contrib;
      const double cav = std::abs(contrib);
      if (cav > max_contrib_abs) {
        max_contrib_abs = cav;
        max_contrib = contrib;
        max_col = j;
      }
    }
    row_res -= rhs[i];
    const Complex row_adjoint = column_adjoint_residual(i);
    const Complex dir_adjoint = column_adjoint_residual(j_dir);
    const auto x = node->getPosition();
    const auto* max_node = std::get<0>(id_by_index[max_col]);
    const auto* max_face = std::get<1>(id_by_index[max_col]);

    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << ","
       << i << "," << j_dir << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << reinterpret_cast<std::uintptr_t>(f_row) << ","
       << x[0] << "," << x[1] << "," << x[2] << ","
       << radiation_surface_name(radiation_surface_kind(node, f_row, face_sets)) << ","
       << row_kind_name(row_kind(node, f_row)) << "," << key_kind_name(f_row) << ","
       << unknown_name(node, f_row) << ","
       << radiation_surface_name(radiation_surface_kind(node, f_dir, face_sets)) << ","
       << key_kind_name(f_dir) << "," << unknown_name(node, f_dir) << ","
       << std::sqrt(static_cast<double>(row_l2)) << "," << static_cast<double>(row_sum) << ","
       << rhs[i].real() << "," << rhs[i].imag() << "," << std::abs(rhs[i]) << ","
       << u[i].real() << "," << u[i].imag() << ","
       << row_res.real() << "," << row_res.imag() << "," << std::abs(row_res) << ","
       << row_adjoint.real() << "," << row_adjoint.imag() << "," << std::abs(row_adjoint) << ","
       << dir_adjoint.real() << "," << dir_adjoint.imag() << "," << std::abs(dir_adjoint) << ","
       << max_col << ","
       << radiation_surface_name(radiation_surface_kind(max_node, max_face, face_sets)) << ","
       << key_kind_name(max_face) << "," << unknown_name(max_node, max_face) << ","
       << max_contrib.real() << "," << max_contrib.imag() << "," << std::abs(max_contrib) << ","
       << (interface_row_scale_mode.empty() ? "default" : interface_row_scale_mode) << ","
       << (interface_robin_adjoint_probe ? 1 : 0) << "\n";
  }
}

void append_radiation_robin_interface_adjoint_pairs(const std::filesystem::path& path,
                                                    double omega,
                                                    int dof_col,
                                                    const char* solver_name,
                                                    const std::vector<Complex>& A_col_major,
                                                    const std::vector<Complex>& rhs,
                                                    const std::vector<Complex>& u,
                                                    const std::vector<bem_frequency_domain::Id>& id_by_index,
                                                    const std::unordered_set<networkFace*>& robin_faces,
                                                    const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                                    const FaceSets& face_sets,
                                                    bool robin_unknown_phi,
                                                    const std::string& interface_row_scale_mode,
                                                    bool interface_robin_adjoint_probe) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || A_col_major.size() != n * n || rhs.size() != n || u.size() != n)
    return;

  auto aidx = [&](std::size_t row, std::size_t col) { return row + col * n; };
  auto key_kind_name = [](const networkFace* f) { return f ? "face_key" : "default_key"; };
  auto row_kind = [](const BEM_DOF_Base* node, const networkFace* f) {
    return (node && node->BCInterface && isNeumannBieDofKey(node, f)) ? 1 : 0;
  };
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto unknown_name = [&](const BEM_DOF_Base* node, const networkFace* f) {
    const bool robin = is_robin_id(node, f);
    if (robin)
      return robin_unknown_phi ? "robin_phi" : "robin_phin";
    if (isDirichletBieDofKey(node, f))
      return "dirichlet_phin";
    return "neumann_phi";
  };
  auto rel = [](Complex a, Complex b, Complex diff) {
    return std::abs(diff) / std::max({std::abs(a), std::abs(b), 1e-300});
  };

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,row_index,dir_index,node_ptr,row_face_ptr,x,y,z,"
          "row_surface,row_kind,row_key_kind,row_unknown,dir_surface,dir_key_kind,dir_unknown,"
          "A_row_dir_re,A_row_dir_im,A_dir_row_re,A_dir_row_im,"
          "pair_diff_re,pair_diff_im,pair_diff_abs,pair_diff_rel,"
          "A_row_row_re,A_row_row_im,A_dir_dir_re,A_dir_dir_im,"
          "rhs_row_re,rhs_row_im,rhs_dir_re,rhs_dir_im,u_row_re,u_row_im,u_dir_re,u_dir_im,"
          "interface_row_scale_mode,interface_robin_adjoint_probe\n";
  }

  for (std::size_t i = 0; i < n; ++i) {
    const auto* node = std::get<0>(id_by_index[i]);
    const auto* f_row = std::get<1>(id_by_index[i]);
    if (!node || row_kind(node, f_row) != 1 || !is_robin_id(node, nullptr))
      continue;
    const auto* d_dir = node->findActiveBieDof(nullptr);
    if (!d_dir || d_dir->index < 0 || static_cast<std::size_t>(d_dir->index) >= n)
      continue;
    const std::size_t j_dir = static_cast<std::size_t>(d_dir->index);
    const auto* f_dir = std::get<1>(id_by_index[j_dir]);
    const Complex a_row_dir = A_col_major[aidx(i, j_dir)];
    const Complex a_dir_row = A_col_major[aidx(j_dir, i)];
    const Complex diff = a_row_dir - a_dir_row;
    const auto x = node->getPosition();

    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << ","
       << i << "," << j_dir << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << reinterpret_cast<std::uintptr_t>(f_row) << ","
       << x[0] << "," << x[1] << "," << x[2] << ","
       << radiation_surface_name(radiation_surface_kind(node, f_row, face_sets)) << ","
       << "interface_constraint," << key_kind_name(f_row) << "," << unknown_name(node, f_row) << ","
       << radiation_surface_name(radiation_surface_kind(node, f_dir, face_sets)) << ","
       << key_kind_name(f_dir) << "," << unknown_name(node, f_dir) << ","
       << a_row_dir.real() << "," << a_row_dir.imag() << ","
       << a_dir_row.real() << "," << a_dir_row.imag() << ","
       << diff.real() << "," << diff.imag() << ","
       << std::abs(diff) << "," << rel(a_row_dir, a_dir_row, diff) << ","
       << A_col_major[aidx(i, i)].real() << "," << A_col_major[aidx(i, i)].imag() << ","
       << A_col_major[aidx(j_dir, j_dir)].real() << "," << A_col_major[aidx(j_dir, j_dir)].imag() << ","
       << rhs[i].real() << "," << rhs[i].imag() << ","
       << rhs[j_dir].real() << "," << rhs[j_dir].imag() << ","
       << u[i].real() << "," << u[i].imag() << ","
       << u[j_dir].real() << "," << u[j_dir].imag() << ","
       << (interface_row_scale_mode.empty() ? "default" : interface_row_scale_mode) << ","
       << (interface_robin_adjoint_probe ? 1 : 0) << "\n";
  }
}

void append_radiation_condensed_system_reciprocity(const std::filesystem::path& path,
                                                   double omega,
                                                   int dof_col,
                                                   const char* solver_name,
                                                   const std::vector<Complex>& A_col_major,
                                                   const std::vector<Complex>& rhs,
                                                   const std::vector<Complex>& u,
                                                   const std::vector<bem_frequency_domain::Id>& retained_ids,
                                                   const FaceSets& face_sets,
                                                   const std::unordered_set<networkFace*>& robin_faces,
                                                   const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                                   bool robin_unknown_phi,
                                                   const std::string& interface_row_scale_mode,
                                                   bool interface_robin_adjoint_probe,
                                                   std::size_t original_size,
                                                   std::size_t eliminated_count) {
  const std::size_t n = retained_ids.size();
  if (n == 0 || A_col_major.size() != n * n || rhs.size() != n || u.size() != n)
    return;

  auto aidx = [&](std::size_t row, std::size_t col) { return row + col * n; };
  auto key_kind_name = [](const networkFace* f) { return f ? "face_key" : "default_key"; };
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto unknown_name = [&](const BEM_DOF_Base* node, const networkFace* f) {
    const bool robin = is_robin_id(node, f);
    if (robin)
      return robin_unknown_phi ? "robin_phi" : "robin_phin";
    if (isDirichletBieDofKey(node, f))
      return "dirichlet_phin";
    return "neumann_phi";
  };
  auto rel = [](Complex a, Complex b, Complex diff) {
    return std::abs(diff) / std::max({std::abs(a), std::abs(b), 1e-300});
  };

  double max_abs = 0.0;
  double max_rel = 0.0;
  std::size_t worst_i = 0;
  std::size_t worst_j = 0;
  Complex worst_aij{0.0, 0.0};
  Complex worst_aji{0.0, 0.0};
  Complex worst_diff{0.0, 0.0};
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const Complex aij = A_col_major[aidx(i, j)];
      const Complex aji = A_col_major[aidx(j, i)];
      const Complex diff = aij - aji;
      const double abs_diff = std::abs(diff);
      const double rel_diff = rel(aij, aji, diff);
      if (rel_diff > max_rel || (rel_diff == max_rel && abs_diff > max_abs)) {
        max_abs = abs_diff;
        max_rel = rel_diff;
        worst_i = i;
        worst_j = j;
        worst_aij = aij;
        worst_aji = aji;
        worst_diff = diff;
      }
    }
  }

  long double rhs_l2 = 0.0L;
  long double u_l2 = 0.0L;
  long double residual_l2 = 0.0L;
  for (std::size_t i = 0; i < n; ++i) {
    rhs_l2 += static_cast<long double>(std::abs(rhs[i])) * std::abs(rhs[i]);
    u_l2 += static_cast<long double>(std::abs(u[i])) * std::abs(u[i]);
    Complex residual = -rhs[i];
    for (std::size_t j = 0; j < n; ++j)
      residual += A_col_major[aidx(i, j)] * u[j];
    residual_l2 += static_cast<long double>(std::abs(residual)) * std::abs(residual);
  }

  const auto* row_node = std::get<0>(retained_ids[worst_i]);
  const auto* row_face = std::get<1>(retained_ids[worst_i]);
  const auto* col_node = std::get<0>(retained_ids[worst_j]);
  const auto* col_face = std::get<1>(retained_ids[worst_j]);

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,n_retained,n_original,n_eliminated,"
          "matrix_asym_abs,matrix_asym_rel,worst_row,worst_col,"
          "A_row_col_re,A_row_col_im,A_col_row_re,A_col_row_im,diff_re,diff_im,"
          "worst_row_surface,worst_row_key_kind,worst_row_unknown,"
          "worst_col_surface,worst_col_key_kind,worst_col_unknown,"
          "rhs_l2,u_l2,residual_l2,interface_row_scale_mode,interface_robin_adjoint_probe\n";
  }
  fs << std::scientific << std::setprecision(12)
     << omega << "," << dof_col << "," << solver_name << ","
     << n << "," << original_size << "," << eliminated_count << ","
     << max_abs << "," << max_rel << "," << worst_i << "," << worst_j << ","
     << worst_aij.real() << "," << worst_aij.imag() << ","
     << worst_aji.real() << "," << worst_aji.imag() << ","
     << worst_diff.real() << "," << worst_diff.imag() << ","
     << radiation_surface_name(radiation_surface_kind(row_node, row_face, face_sets)) << ","
     << key_kind_name(row_face) << "," << unknown_name(row_node, row_face) << ","
     << radiation_surface_name(radiation_surface_kind(col_node, col_face, face_sets)) << ","
     << key_kind_name(col_face) << "," << unknown_name(col_node, col_face) << ","
     << std::sqrt(static_cast<double>(rhs_l2)) << ","
     << std::sqrt(static_cast<double>(u_l2)) << ","
     << std::sqrt(static_cast<double>(residual_l2)) << ","
     << (interface_row_scale_mode.empty() ? "default" : interface_row_scale_mode) << ","
     << (interface_robin_adjoint_probe ? 1 : 0) << "\n";
}

void append_radiation_condensed_interface_map(const std::filesystem::path& path,
                                              double omega,
                                              int dof_col,
                                              const char* solver_name,
                                              const std::vector<bem_frequency_domain::Id>& original_ids,
                                              const std::vector<bem_frequency_domain::Id>& retained_ids,
                                              const std::vector<int>& retained_index,
                                              const std::vector<bem_frequency_domain::CondensedInterfaceExpr>& eliminated,
                                              const std::unordered_set<networkFace*>& robin_faces,
                                              const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin,
                                              const FaceSets& face_sets,
                                              bool robin_unknown_phi,
                                              const std::string& interface_row_scale_mode,
                                              bool interface_robin_adjoint_probe) {
  const std::size_t n = original_ids.size();
  if (retained_index.size() != n || eliminated.size() != n)
    return;

  auto key_kind_name = [](const networkFace* f) { return f ? "face_key" : "default_key"; };
  auto is_robin_id = [&](const BEM_DOF_Base* node, const networkFace* f) {
    if (f)
      return robin_faces.contains(const_cast<networkFace*>(f));
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto is_robin_node = [&](const BEM_DOF_Base* node) {
    const auto it = node_is_robin.find(const_cast<BEM_DOF_Base*>(node));
    return it != node_is_robin.end() && it->second;
  };
  auto unknown_name = [&](const BEM_DOF_Base* node, const networkFace* f) {
    const bool robin = is_robin_id(node, f);
    if (robin)
      return robin_unknown_phi ? "robin_phi" : "robin_phin";
    if (isDirichletBieDofKey(node, f))
      return "dirichlet_phin";
    return "neumann_phi";
  };

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header) {
    fs << "omega,dof,solver,original_index,retained_index,is_eliminated,dir_index,dir_retained_index,"
          "node_ptr,face_ptr,x,y,z,surface,key_kind,unknown,is_robin_node,is_robin_id,"
          "dir_node_ptr,dir_face_ptr,dir_surface,dir_key_kind,dir_unknown,"
          "scale_re,scale_im,constant_re,constant_im,interface_row_scale_mode,interface_robin_adjoint_probe\n";
  }

  for (std::size_t i = 0; i < n; ++i) {
    const auto* node = std::get<0>(original_ids[i]);
    const auto* face = std::get<1>(original_ids[i]);
    const auto& expr = eliminated[i];
    const int dir_index = expr.eliminated ? static_cast<int>(expr.dir_index) : -1;
    const int dir_retained = (expr.eliminated && expr.dir_index < retained_index.size()) ? retained_index[expr.dir_index] : -1;
    const auto* dir_node = (expr.eliminated && expr.dir_index < original_ids.size()) ? std::get<0>(original_ids[expr.dir_index]) : nullptr;
    const auto* dir_face = (expr.eliminated && expr.dir_index < original_ids.size()) ? std::get<1>(original_ids[expr.dir_index]) : nullptr;
    const auto x = node ? node->getPosition() : Tddd{0.0, 0.0, 0.0};

    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << ","
       << i << "," << retained_index[i] << "," << (expr.eliminated ? 1 : 0) << ","
       << dir_index << "," << dir_retained << ","
       << reinterpret_cast<std::uintptr_t>(node) << ","
       << reinterpret_cast<std::uintptr_t>(face) << ","
       << x[0] << "," << x[1] << "," << x[2] << ","
       << radiation_surface_name(radiation_surface_kind(node, face, face_sets)) << ","
       << key_kind_name(face) << "," << unknown_name(node, face) << ","
       << (is_robin_node(node) ? 1 : 0) << "," << (is_robin_id(node, face) ? 1 : 0) << ","
       << reinterpret_cast<std::uintptr_t>(dir_node) << ","
       << reinterpret_cast<std::uintptr_t>(dir_face) << ","
       << radiation_surface_name(radiation_surface_kind(dir_node, dir_face, face_sets)) << ","
       << key_kind_name(dir_face) << "," << unknown_name(dir_node, dir_face) << ","
       << expr.scale.real() << "," << expr.scale.imag() << ","
       << expr.constant.real() << "," << expr.constant.imag() << ","
       << (interface_row_scale_mode.empty() ? "default" : interface_row_scale_mode) << ","
       << (interface_robin_adjoint_probe ? 1 : 0) << "\n";
  }
}

std::vector<double> lumped_area_for_ids(const std::vector<bem_frequency_domain::Id>& id_by_index,
                                        const std::unordered_set<networkFace*>& faces) {
  const std::size_t n = id_by_index.size();
  std::unordered_map<bem_frequency_domain::Id, std::size_t, bem_frequency_domain::IdHash, bem_frequency_domain::IdEq> index;
  index.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    index.emplace(id_by_index[i], i);

  std::vector<double> lumped_area(n, 0.0);
  auto add_lumped_area = [&](BEM_DOF_Base* node, networkFace* f, double value) {
    bem_frequency_domain::Id key{node, f};
    auto it = index.find(key);
    if (it == index.end()) {
      key = {node, nullptr};
      it = index.find(key);
    }
    if (it != index.end())
      lumped_area[it->second] += value;
  };

  for (auto* f : faces) {
    if (!f)
      continue;
    auto [p0, p1, p2] = f->getPoints();
    const double a3 = f->area / 3.0;
    add_lumped_area(p0, f, a3);
    add_lumped_area(p1, f, a3);
    add_lumped_area(p2, f, a3);
  }
  return lumped_area;
}

void symmetrize_area_weighted_bie_operator(BEM_BVP& bvp,
                                           const std::vector<bem_frequency_domain::Id>& id_by_index,
                                           const std::unordered_set<networkFace*>& faces) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || bvp.IGIGn.size() != n)
    return;
  const auto area = lumped_area_for_ids(id_by_index, faces);

  for (std::size_t i = 0; i < n; ++i) {
    if (area[i] <= 0.0 || bvp.IGIGn[i].size() != n)
      continue;
    for (std::size_t j = i + 1; j < n; ++j) {
      if (area[j] <= 0.0 || bvp.IGIGn[j].size() != n)
        continue;
      for (int k = 0; k < 2; ++k) {
        const double kij = area[i] * bvp.IGIGn[i][j][static_cast<std::size_t>(k)];
        const double kji = area[j] * bvp.IGIGn[j][i][static_cast<std::size_t>(k)];
        const double ksym = 0.5 * (kij + kji);
        bvp.IGIGn[i][j][static_cast<std::size_t>(k)] = ksym / area[i];
        bvp.IGIGn[j][i][static_cast<std::size_t>(k)] = ksym / area[j];
      }
    }
  }
}

bool faces_share_point(const networkFace* a, const networkFace* b) {
  if (!a || !b)
    return false;
  auto [a0, a1, a2] = a->getPoints();
  auto [b0, b1, b2] = b->getPoints();
  const std::array<const networkPoint*, 3> ap{a0, a1, a2};
  const std::array<const networkPoint*, 3> bp{b0, b1, b2};
  for (const auto* pa : ap)
    for (const auto* pb : bp)
      if (pa == pb)
        return true;
  return false;
}

double integrate_panel_single_layer_at_point(const Tddd& x, const networkFace* f) {
  if (!f)
    return 0.0;

  // Six-point Dunavant rule on a linear triangle.  The weights below sum to
  // 1/2 on the reference triangle, so multiply by 2*area.
  constexpr std::array<std::array<double, 4>, 6> q = {{
      {0.816847572980459, 0.091576213509771, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.054975871827661},
      {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.111690794839006},
  }};

  auto [p0, p1, p2] = f->getPoints();
  double out = 0.0;
  for (const auto& e : q) {
    const Tddd y = e[0] * p0->X + e[1] * p1->X + e[2] * p2->X;
    const double r = Norm(y - x);
    if (r > 1e-14)
      out += e[3] / r;
  }
  return 2.0 * f->area * out;
}

double integrate_panel_pair_single_layer_quad6(const networkFace* fi, const networkFace* fj) {
  if (!fi || !fj)
    return 0.0;

  // Tensor-product Dunavant rule for the constant-panel double integral
  // int_{fi} int_{fj} 1/|x-y| dS_y dS_x. This is a diagnostic integrator;
  // shared/self panels still need a singular Duffy treatment for accuracy.
  constexpr std::array<std::array<double, 4>, 6> q = {{
      {0.816847572980459, 0.091576213509771, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.054975871827661},
      {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.111690794839006},
  }};

  auto [i0, i1, i2] = fi->getPoints();
  auto [j0, j1, j2] = fj->getPoints();
  double out = 0.0;
  for (const auto& qi : q) {
    const Tddd x = qi[0] * i0->X + qi[1] * i1->X + qi[2] * i2->X;
    for (const auto& qj : q) {
      const Tddd y = qj[0] * j0->X + qj[1] * j1->X + qj[2] * j2->X;
      const double r = Norm(y - x);
      if (r > 1e-14)
        out += qi[3] * qj[3] / r;
    }
  }
  return 4.0 * fi->area * fj->area * out;
}

double integrate_panel_double_layer_source_normal_at_point(const Tddd& x, const networkFace* f) {
  if (!f)
    return 0.0;

  constexpr std::array<std::array<double, 4>, 6> q = {{
      {0.816847572980459, 0.091576213509771, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.054975871827661},
      {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.111690794839006},
  }};

  auto [p0, p1, p2] = f->getPoints();
  double out = 0.0;
  for (const auto& e : q) {
    const Tddd y = e[0] * p0->X + e[1] * p1->X + e[2] * p2->X;
    const Tddd rvec = x - y;
    const double r = Norm(rvec);
    if (r > 1e-14)
      out += e[3] * Dot(f->normal, rvec) / (r * r * r);
  }
  return 2.0 * f->area * out;
}

double integrate_panel_double_layer_target_normal_at_point(const Tddd& x, const Tddd& target_normal, const networkFace* f) {
  if (!f)
    return 0.0;

  constexpr std::array<std::array<double, 4>, 6> q = {{
      {0.816847572980459, 0.091576213509771, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.054975871827661},
      {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.054975871827661},
      {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.111690794839006},
      {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.111690794839006},
  }};

  auto [p0, p1, p2] = f->getPoints();
  double out = 0.0;
  for (const auto& e : q) {
    const Tddd y = e[0] * p0->X + e[1] * p1->X + e[2] * p2->X;
    const Tddd rvec = y - x;
    const double r = Norm(rvec);
    if (r > 1e-14)
      out += e[3] * Dot(target_normal, rvec) / (r * r * r);
  }
  return 2.0 * f->area * out;
}

void write_panel_operator_symmetry(const std::filesystem::path& path,
                                   double omega,
                                   const char* solver_name,
                                   const std::unordered_set<networkFace*>& faces) {
  std::vector<networkFace*> panels;
  panels.reserve(faces.size());
  for (auto* f : faces) {
    if (f)
      panels.push_back(f);
  }
  std::sort(panels.begin(), panels.end(), [](const networkFace* a, const networkFace* b) {
    return a->index < b->index;
  });

  struct Stat {
    std::uint64_t pairs = 0;
    double max_abs = 0.0;
    double max_rel = 0.0;
    int worst_i = -1;
    int worst_j = -1;
    double value_ij = 0.0;
    double value_ji = 0.0;
    double area_i = 0.0;
    double area_j = 0.0;
  };

  auto update = [](Stat& st, const networkFace* fi, const networkFace* fj, double vij, double vji) {
    ++st.pairs;
    const double diff = std::abs(vij - vji);
    const double rel = diff / std::max({std::abs(vij), std::abs(vji), 1e-300});
    if (rel > st.max_rel || (rel == st.max_rel && diff > st.max_abs)) {
      st.max_abs = diff;
      st.max_rel = rel;
      st.worst_i = fi ? fi->index : -1;
      st.worst_j = fj ? fj->index : -1;
      st.value_ij = vij;
      st.value_ji = vji;
      st.area_i = fi ? fi->area : 0.0;
      st.area_j = fj ? fj->area : 0.0;
    }
  };

  Stat quad_all;
  Stat quad_disjoint;
  Stat quad_far_disjoint;
  Stat centroid_all;
  Stat centroid_disjoint;
  Stat centroid_far_disjoint;

  for (std::size_t i = 0; i < panels.size(); ++i) {
    auto* fi = panels[i];
    for (std::size_t j = i + 1; j < panels.size(); ++j) {
      auto* fj = panels[j];
      const double gij = integrate_panel_single_layer_at_point(fi->center, fj);
      const double gji = integrate_panel_single_layer_at_point(fj->center, fi);
      const double weighted_ij = fi->area * gij;
      const double weighted_ji = fj->area * gji;
      update(quad_all, fi, fj, weighted_ij, weighted_ji);

      const double rc = Norm(fj->center - fi->center);
      const double centroid_ij = rc > 1e-14 ? fi->area * fj->area / rc : 0.0;
      const double centroid_ji = rc > 1e-14 ? fj->area * fi->area / rc : 0.0;
      update(centroid_all, fi, fj, centroid_ij, centroid_ji);

      if (!faces_share_point(fi, fj)) {
        update(quad_disjoint, fi, fj, weighted_ij, weighted_ji);
        update(centroid_disjoint, fi, fj, centroid_ij, centroid_ji);
        const double distance = Norm(fj->center - fi->center);
        const double panel_length = std::max(std::sqrt(std::max(fi->area, 0.0)), std::sqrt(std::max(fj->area, 0.0)));
        if (distance > 5.0 * panel_length) {
          update(quad_far_disjoint, fi, fj, weighted_ij, weighted_ji);
          update(centroid_far_disjoint, fi, fj, centroid_ij, centroid_ji);
        }
      }
    }
  }

  std::ofstream fs(path);
  fs << "omega,solver,surface,panel_count,metric,pairs,max_abs,max_rel,"
        "worst_face_i,worst_face_j,value_ij,value_ji,area_i,area_j\n";
  auto write = [&](const char* metric, const Stat& st) {
    fs << std::scientific << std::setprecision(12)
       << omega << "," << solver_name << ",float,"
       << panels.size() << "," << metric << ","
       << st.pairs << "," << st.max_abs << "," << st.max_rel << ","
       << st.worst_i << "," << st.worst_j << ","
       << st.value_ij << "," << st.value_ji << ","
       << st.area_i << "," << st.area_j << "\n";
  };
  write("single_layer_quad6_area_weighted_all", quad_all);
  write("single_layer_quad6_area_weighted_disjoint", quad_disjoint);
  write("single_layer_quad6_area_weighted_far_disjoint", quad_far_disjoint);
  write("single_layer_centroid_area_weighted_all", centroid_all);
  write("single_layer_centroid_area_weighted_disjoint", centroid_disjoint);
  write("single_layer_centroid_area_weighted_far_disjoint", centroid_far_disjoint);
}

void write_panel_pair_integral_check(const std::filesystem::path& path,
                                     double omega,
                                     const char* solver_name,
                                     const std::unordered_set<networkFace*>& faces,
                                     const Tddd& ref_point,
                                     double moment_sign) {
  std::vector<networkFace*> panels;
  panels.reserve(faces.size());
  for (auto* f : faces) {
    if (f)
      panels.push_back(f);
  }
  std::sort(panels.begin(), panels.end(), [](const networkFace* a, const networkFace* b) {
    return a->index < b->index;
  });

  struct Stat {
    std::uint64_t pairs = 0;
    double max_abs = 0.0;
    double max_rel = 0.0;
    int worst_i = -1;
    int worst_j = -1;
    double value_ij = 0.0;
    double value_ji = 0.0;
  };
  auto update = [](Stat& st, const networkFace* fi, const networkFace* fj, double vij, double vji) {
    ++st.pairs;
    const double diff = std::abs(vij - vji);
    const double rel = diff / std::max({std::abs(vij), std::abs(vji), 1e-300});
    if (rel > st.max_rel || (rel == st.max_rel && diff > st.max_abs)) {
      st.max_abs = diff;
      st.max_rel = rel;
      st.worst_i = fi ? fi->index : -1;
      st.worst_j = fj ? fj->index : -1;
      st.value_ij = vij;
      st.value_ji = vji;
    }
  };
  auto mode_at_center = [&](const networkFace* f) {
    std::array<double, 6> m{};
    for (int d = 0; d < 6; ++d)
      m[static_cast<std::size_t>(d)] = rigid_mode_normal_component(d, f->center, f->normal, ref_point, moment_sign);
    return m;
  };
  auto add_modal = [](Matrix6& M, const std::array<double, 6>& mi, const std::array<double, 6>& mj, double k) {
    for (int r = 0; r < 6; ++r) {
      for (int c = 0; c < 6; ++c)
        M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] += mi[static_cast<std::size_t>(r)] * k * mj[static_cast<std::size_t>(c)];
    }
  };

  Stat all;
  Stat disjoint;
  Stat far;
  Matrix6 modal_all = make_zero_matrix6();
  Matrix6 modal_disjoint = make_zero_matrix6();
  Matrix6 modal_far = make_zero_matrix6();
  std::uint64_t modal_pairs_all = 0;
  std::uint64_t modal_pairs_disjoint = 0;
  std::uint64_t modal_pairs_far = 0;

  std::vector<std::array<double, 6>> modes;
  modes.reserve(panels.size());
  for (auto* f : panels)
    modes.push_back(mode_at_center(f));

  for (std::size_t i = 0; i < panels.size(); ++i) {
    auto* fi = panels[i];
    for (std::size_t j = i + 1; j < panels.size(); ++j) {
      auto* fj = panels[j];
      const double kij = integrate_panel_pair_single_layer_quad6(fi, fj);
      const double kji = integrate_panel_pair_single_layer_quad6(fj, fi);
      update(all, fi, fj, kij, kji);
      add_modal(modal_all, modes[i], modes[j], kij);
      add_modal(modal_all, modes[j], modes[i], kij);
      modal_pairs_all += 2;

      const bool disjoint_pair = !faces_share_point(fi, fj);
      const double distance = Norm(fj->center - fi->center);
      const double panel_length = std::max(std::sqrt(std::max(fi->area, 0.0)), std::sqrt(std::max(fj->area, 0.0)));
      const bool far_pair = disjoint_pair && distance > 5.0 * panel_length;
      if (disjoint_pair) {
        update(disjoint, fi, fj, kij, kji);
        add_modal(modal_disjoint, modes[i], modes[j], kij);
        add_modal(modal_disjoint, modes[j], modes[i], kij);
        modal_pairs_disjoint += 2;
      }
      if (far_pair) {
        update(far, fi, fj, kij, kji);
        add_modal(modal_far, modes[i], modes[j], kij);
        add_modal(modal_far, modes[j], modes[i], kij);
        modal_pairs_far += 2;
      }
    }
  }

  auto modal_stat = [](const Matrix6& M) {
    Stat st;
    for (int i = 0; i < 6; ++i) {
      for (int j = i + 1; j < 6; ++j) {
        const double a = M[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        const double b = M[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
        ++st.pairs;
        const double diff = std::abs(a - b);
        const double rel = diff / std::max({std::abs(a), std::abs(b), 1e-300});
        if (rel > st.max_rel || (rel == st.max_rel && diff > st.max_abs)) {
          st.max_abs = diff;
          st.max_rel = rel;
          st.worst_i = i;
          st.worst_j = j;
          st.value_ij = a;
          st.value_ji = b;
        }
      }
    }
    return st;
  };

  std::ofstream fs(path);
  fs << "omega,solver,surface,metric,panel_count,pairs,max_abs,max_rel,worst_i,worst_j,value_ij,value_ji,note\n";
  auto write = [&](const char* metric, std::uint64_t pairs, const Stat& st, const char* note) {
    fs << std::scientific << std::setprecision(12)
       << omega << "," << solver_name << ",float," << metric << ","
       << panels.size() << "," << pairs << ","
       << st.max_abs << "," << st.max_rel << ","
       << st.worst_i << "," << st.worst_j << ","
       << st.value_ij << "," << st.value_ji << ","
       << note << "\n";
  };
  write("constant_panel_double_integral_quad6_all", all.pairs, all, "tensor_product_no_duffy");
  write("constant_panel_double_integral_quad6_disjoint", disjoint.pairs, disjoint, "tensor_product_no_duffy");
  write("constant_panel_double_integral_quad6_far_disjoint", far.pairs, far, "tensor_product_no_duffy");
  write("constant_panel_double_integral_quad6_modal_all", modal_pairs_all, modal_stat(modal_all), "center_modes_tensor_product_no_duffy");
  write("constant_panel_double_integral_quad6_modal_disjoint", modal_pairs_disjoint, modal_stat(modal_disjoint), "center_modes_tensor_product_no_duffy");
  write("constant_panel_double_integral_quad6_modal_far_disjoint", modal_pairs_far, modal_stat(modal_far), "center_modes_tensor_product_no_duffy");
}

void write_panel_double_layer_adjoint(const std::filesystem::path& path,
                                      double omega,
                                      const char* solver_name,
                                      const std::unordered_set<networkFace*>& faces) {
  std::vector<networkFace*> panels;
  panels.reserve(faces.size());
  for (auto* f : faces) {
    if (f)
      panels.push_back(f);
  }
  std::sort(panels.begin(), panels.end(), [](const networkFace* a, const networkFace* b) {
    return a->index < b->index;
  });

  struct Stat {
    std::uint64_t pairs = 0;
    double max_abs = 0.0;
    double max_rel = 0.0;
    int worst_i = -1;
    int worst_j = -1;
    double value_ij = 0.0;
    double value_ji_or_adjoint = 0.0;
    double area_i = 0.0;
    double area_j = 0.0;
  };

  auto update = [](Stat& st, const networkFace* fi, const networkFace* fj, double a, double b) {
    ++st.pairs;
    const double diff = std::abs(a - b);
    const double rel = diff / std::max({std::abs(a), std::abs(b), 1e-300});
    if (rel > st.max_rel || (rel == st.max_rel && diff > st.max_abs)) {
      st.max_abs = diff;
      st.max_rel = rel;
      st.worst_i = fi ? fi->index : -1;
      st.worst_j = fj ? fj->index : -1;
      st.value_ij = a;
      st.value_ji_or_adjoint = b;
      st.area_i = fi ? fi->area : 0.0;
      st.area_j = fj ? fj->area : 0.0;
    }
  };

  Stat quad_source_self;
  Stat quad_source_adjoint;
  Stat quad_source_adjoint_disjoint;
  Stat quad_source_adjoint_far;
  Stat centroid_source_self;
  Stat centroid_source_adjoint;

  for (std::size_t i = 0; i < panels.size(); ++i) {
    auto* fi = panels[i];
    for (std::size_t j = i + 1; j < panels.size(); ++j) {
      auto* fj = panels[j];
      const double dij = integrate_panel_double_layer_source_normal_at_point(fi->center, fj);
      const double dji = integrate_panel_double_layer_source_normal_at_point(fj->center, fi);
      const double weighted_ij = fi->area * dij;
      const double weighted_ji = fj->area * dji;
      update(quad_source_self, fi, fj, weighted_ij, weighted_ji);

      const double adjoint_ji = integrate_panel_double_layer_target_normal_at_point(fj->center, fj->normal, fi);
      const double weighted_adjoint_ji = fj->area * adjoint_ji;
      update(quad_source_adjoint, fi, fj, weighted_ij, weighted_adjoint_ji);

      const double rc = Norm(fj->center - fi->center);
      if (rc > 1e-14) {
        const Tddd rij = fi->center - fj->center;
        const Tddd rji = fj->center - fi->center;
        const double centroid_ij = fi->area * fj->area * Dot(fj->normal, rij) / (rc * rc * rc);
        const double centroid_ji = fj->area * fi->area * Dot(fi->normal, rji) / (rc * rc * rc);
        const double centroid_adjoint_ji = fj->area * fi->area * Dot(fj->normal, rij) / (rc * rc * rc);
        update(centroid_source_self, fi, fj, centroid_ij, centroid_ji);
        update(centroid_source_adjoint, fi, fj, centroid_ij, centroid_adjoint_ji);
      }

      if (!faces_share_point(fi, fj)) {
        update(quad_source_adjoint_disjoint, fi, fj, weighted_ij, weighted_adjoint_ji);
        const double panel_length = std::max(std::sqrt(std::max(fi->area, 0.0)), std::sqrt(std::max(fj->area, 0.0)));
        if (rc > 5.0 * panel_length)
          update(quad_source_adjoint_far, fi, fj, weighted_ij, weighted_adjoint_ji);
      }
    }
  }

  std::ofstream fs(path);
  fs << "omega,solver,surface,panel_count,metric,pairs,max_abs,max_rel,"
        "worst_face_i,worst_face_j,value_ij,value_ji_or_adjoint,area_i,area_j\n";
  auto write = [&](const char* metric, const Stat& st) {
    fs << std::scientific << std::setprecision(12)
       << omega << "," << solver_name << ",float,"
       << panels.size() << "," << metric << ","
       << st.pairs << "," << st.max_abs << "," << st.max_rel << ","
       << st.worst_i << "," << st.worst_j << ","
       << st.value_ij << "," << st.value_ji_or_adjoint << ","
       << st.area_i << "," << st.area_j << "\n";
  };
  write("double_layer_source_normal_self_sym_quad6_all", quad_source_self);
  write("double_layer_source_normal_vs_adjoint_quad6_all", quad_source_adjoint);
  write("double_layer_source_normal_vs_adjoint_quad6_disjoint", quad_source_adjoint_disjoint);
  write("double_layer_source_normal_vs_adjoint_quad6_far_disjoint", quad_source_adjoint_far);
  write("double_layer_source_normal_self_sym_centroid_all", centroid_source_self);
  write("double_layer_source_normal_vs_adjoint_centroid_all", centroid_source_adjoint);
}

void write_panel_mode_projection(const std::filesystem::path& path,
                                 double omega,
                                 const char* solver_name,
                                 const std::unordered_set<networkFace*>& faces,
                                 const Tddd& ref_point,
                                 double moment_sign) {
  std::vector<networkFace*> panels;
  panels.reserve(faces.size());
  for (auto* f : faces) {
    if (f)
      panels.push_back(f);
  }
  std::sort(panels.begin(), panels.end(), [](const networkFace* a, const networkFace* b) {
    return a->index < b->index;
  });

  struct PanelMode {
    networkFace* f = nullptr;
    std::array<double, 6> m{};
  };
  std::vector<PanelMode> modes;
  modes.reserve(panels.size());
  for (auto* f : panels) {
    PanelMode pm;
    pm.f = f;
    for (int d = 0; d < 6; ++d)
      pm.m[static_cast<std::size_t>(d)] = rigid_mode_normal_component(d, f->center, f->normal, ref_point, moment_sign);
    modes.push_back(pm);
  }

  struct Projection {
    const char* metric = "";
    std::uint64_t pairs = 0;
    Matrix6 values = make_zero_matrix6();
  };
  Projection all{"single_layer_quad6_panel_mode_all"};
  Projection disjoint{"single_layer_quad6_panel_mode_disjoint"};
  Projection far{"single_layer_quad6_panel_mode_far_disjoint"};
  Projection centroid_all{"single_layer_centroid_panel_mode_all"};
  Projection centroid_disjoint{"single_layer_centroid_panel_mode_disjoint"};
  Projection centroid_far{"single_layer_centroid_panel_mode_far_disjoint"};
  Projection sym_all{"single_layer_quad6_symmetrized_panel_mode_all"};
  Projection sym_disjoint{"single_layer_quad6_symmetrized_panel_mode_disjoint"};
  Projection sym_far{"single_layer_quad6_symmetrized_panel_mode_far_disjoint"};

  auto add_pair = [](Projection& proj, const PanelMode& row, const PanelMode& col, double g_row_col) {
    ++proj.pairs;
    for (int i = 0; i < 6; ++i) {
      const double wi = row.f->area * row.m[static_cast<std::size_t>(i)];
      for (int j = 0; j < 6; ++j)
        proj.values[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
            wi * g_row_col * col.m[static_cast<std::size_t>(j)];
    }
  };
  auto add_kernel_pair = [](Projection& proj, const PanelMode& row, const PanelMode& col, double weighted_kernel) {
    ++proj.pairs;
    for (int i = 0; i < 6; ++i) {
      const double mi = row.m[static_cast<std::size_t>(i)];
      for (int j = 0; j < 6; ++j)
        proj.values[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
            mi * weighted_kernel * col.m[static_cast<std::size_t>(j)];
    }
  };

  for (std::size_t i = 0; i < modes.size(); ++i) {
    for (std::size_t j = 0; j < modes.size(); ++j) {
      const double gij = integrate_panel_single_layer_at_point(modes[i].f->center, modes[j].f);
      add_pair(all, modes[i], modes[j], gij);
      const double distance = Norm(modes[j].f->center - modes[i].f->center);
      const double centroid_gij = distance > 1e-14 ? modes[j].f->area / distance : 0.0;
      add_pair(centroid_all, modes[i], modes[j], centroid_gij);
      if (!faces_share_point(modes[i].f, modes[j].f)) {
        add_pair(disjoint, modes[i], modes[j], gij);
        add_pair(centroid_disjoint, modes[i], modes[j], centroid_gij);
        const double panel_length = std::max(std::sqrt(std::max(modes[i].f->area, 0.0)), std::sqrt(std::max(modes[j].f->area, 0.0)));
        if (distance > 5.0 * panel_length)
          add_pair(far, modes[i], modes[j], gij);
        if (distance > 5.0 * panel_length)
          add_pair(centroid_far, modes[i], modes[j], centroid_gij);
      }
    }
  }
  for (std::size_t i = 0; i < modes.size(); ++i) {
    for (std::size_t j = i; j < modes.size(); ++j) {
      const double gij = integrate_panel_single_layer_at_point(modes[i].f->center, modes[j].f);
      const double gji = integrate_panel_single_layer_at_point(modes[j].f->center, modes[i].f);
      const double kij = modes[i].f->area * gij;
      const double kji = modes[j].f->area * gji;
      const double ksym = 0.5 * (kij + kji);
      const bool disjoint_pair = !faces_share_point(modes[i].f, modes[j].f);
      const double distance = Norm(modes[j].f->center - modes[i].f->center);
      const double panel_length = std::max(std::sqrt(std::max(modes[i].f->area, 0.0)), std::sqrt(std::max(modes[j].f->area, 0.0)));
      const bool far_pair = disjoint_pair && distance > 5.0 * panel_length;

      add_kernel_pair(sym_all, modes[i], modes[j], ksym);
      if (i != j)
        add_kernel_pair(sym_all, modes[j], modes[i], ksym);
      if (disjoint_pair) {
        add_kernel_pair(sym_disjoint, modes[i], modes[j], ksym);
        if (i != j)
          add_kernel_pair(sym_disjoint, modes[j], modes[i], ksym);
      }
      if (far_pair) {
        add_kernel_pair(sym_far, modes[i], modes[j], ksym);
        if (i != j)
          add_kernel_pair(sym_far, modes[j], modes[i], ksym);
      }
    }
  }

  std::ofstream fs(path);
  fs << "omega,solver,surface,metric,panel_count,pairs,dof_row,dof_col,value\n";
  auto write = [&](const Projection& proj) {
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        fs << std::scientific << std::setprecision(12)
           << omega << "," << solver_name << ",float," << proj.metric << ","
           << panels.size() << "," << proj.pairs << ","
           << i << "," << j << ","
           << proj.values[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] << "\n";
      }
    }
  };
  write(all);
  write(disjoint);
  write(far);
  write(centroid_all);
  write(centroid_disjoint);
  write(centroid_far);
  write(sym_all);
  write(sym_disjoint);
  write(sym_far);
}

double matrix6_symmetry_max_rel(const Matrix6& m) {
  double max_rel = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int j = i + 1; j < 6; ++j) {
      const double a = m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      const double b = m[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
      const double rel = std::abs(a - b) / std::max({std::abs(a), std::abs(b), 1e-300});
      max_rel = std::max(max_rel, rel);
    }
  }
  return max_rel;
}

double matrix6_max_abs(const Matrix6& m) {
  double out = 0.0;
  for (const auto& row : m)
    for (double v : row)
      out = std::max(out, std::abs(v));
  return out;
}

bool invert_matrix6(Matrix6 a, Matrix6& inv) {
  inv = make_zero_matrix6();
  for (int i = 0; i < 6; ++i)
    inv[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1.0;

  const double scale = std::max(matrix6_max_abs(a), 1.0);
  for (int col = 0; col < 6; ++col) {
    int pivot = col;
    double pivot_abs = std::abs(a[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)]);
    for (int r = col + 1; r < 6; ++r) {
      const double v = std::abs(a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)]);
      if (v > pivot_abs) {
        pivot_abs = v;
        pivot = r;
      }
    }
    if (pivot_abs <= 1e-13 * scale)
      return false;
    if (pivot != col) {
      std::swap(a[static_cast<std::size_t>(pivot)], a[static_cast<std::size_t>(col)]);
      std::swap(inv[static_cast<std::size_t>(pivot)], inv[static_cast<std::size_t>(col)]);
    }
    const double diag = a[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)];
    for (int j = 0; j < 6; ++j) {
      a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)] /= diag;
      inv[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)] /= diag;
    }
    for (int r = 0; r < 6; ++r) {
      if (r == col)
        continue;
      const double factor = a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)];
      if (factor == 0.0)
        continue;
      for (int j = 0; j < 6; ++j) {
        a[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)] -=
            factor * a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)];
        inv[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)] -=
            factor * inv[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)];
      }
    }
  }
  return true;
}

double matrix6_inverse_residual_max_abs(const Matrix6& a, const Matrix6& inv) {
  double out = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double v = (i == j) ? -1.0 : 0.0;
      for (int k = 0; k < 6; ++k)
        v += a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] *
             inv[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
      out = std::max(out, std::abs(v));
    }
  }
  return out;
}

void write_panel_modal_solve(const std::filesystem::path& path,
                             double omega,
                             const char* solver_name,
                             const std::unordered_set<networkFace*>& faces,
                             const Tddd& ref_point,
                             double moment_sign) {
  std::vector<networkFace*> panels;
  panels.reserve(faces.size());
  for (auto* f : faces) {
    if (f)
      panels.push_back(f);
  }
  std::sort(panels.begin(), panels.end(), [](const networkFace* a, const networkFace* b) {
    return a->index < b->index;
  });
  if (panels.empty())
    return;

  struct PanelMode {
    networkFace* f = nullptr;
    std::array<double, 6> m{};
  };
  std::vector<PanelMode> modes;
  modes.reserve(panels.size());
  for (auto* f : panels) {
    PanelMode pm;
    pm.f = f;
    for (int d = 0; d < 6; ++d)
      pm.m[static_cast<std::size_t>(d)] = rigid_mode_normal_component(d, f->center, f->normal, ref_point, moment_sign);
    modes.push_back(pm);
  }

  struct Projection {
    const char* metric = "";
    std::uint64_t pairs = 0;
    Matrix6 op = make_zero_matrix6();
  };
  Projection quad{"single_layer_quad6_modal_solve"};
  Projection centroid{"single_layer_centroid_modal_solve"};
  Projection sym_quad{"single_layer_quad6_symmetrized_modal_solve"};

  auto add_weighted_kernel = [](Projection& proj, const PanelMode& row, const PanelMode& col, double weighted_kernel) {
    ++proj.pairs;
    for (int i = 0; i < 6; ++i) {
      const double mi = row.m[static_cast<std::size_t>(i)];
      for (int j = 0; j < 6; ++j)
        proj.op[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
            mi * weighted_kernel * col.m[static_cast<std::size_t>(j)];
    }
  };

  for (std::size_t i = 0; i < modes.size(); ++i) {
    for (std::size_t j = 0; j < modes.size(); ++j) {
      const double gij = integrate_panel_single_layer_at_point(modes[i].f->center, modes[j].f);
      add_weighted_kernel(quad, modes[i], modes[j], modes[i].f->area * gij);

      const double distance = Norm(modes[j].f->center - modes[i].f->center);
      const double centroid_kernel = distance > 1e-14 ? modes[i].f->area * modes[j].f->area / distance : 0.0;
      add_weighted_kernel(centroid, modes[i], modes[j], centroid_kernel);
    }
  }
  for (std::size_t i = 0; i < modes.size(); ++i) {
    for (std::size_t j = i; j < modes.size(); ++j) {
      const double gij = integrate_panel_single_layer_at_point(modes[i].f->center, modes[j].f);
      const double gji = integrate_panel_single_layer_at_point(modes[j].f->center, modes[i].f);
      const double kij = modes[i].f->area * gij;
      const double kji = modes[j].f->area * gji;
      const double ksym = 0.5 * (kij + kji);
      add_weighted_kernel(sym_quad, modes[i], modes[j], ksym);
      if (i != j)
        add_weighted_kernel(sym_quad, modes[j], modes[i], ksym);
    }
  }

  std::ofstream fs(path);
  fs << "omega,solver,surface,metric,panel_count,pairs,regularization,dof_row,dof_col,"
        "operator_value,inverse_value,operator_sym_rel,inverse_sym_rel,inverse_residual_max_abs,invertible\n";

  auto write_projection = [&](Projection proj) {
    Matrix6 solve_op = proj.op;
    const double scale = std::max(matrix6_max_abs(solve_op), 1.0);
    const double regularization = 1e-10 * scale;
    for (int i = 0; i < 6; ++i)
      solve_op[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] += regularization;

    Matrix6 inv = make_zero_matrix6();
    const bool invertible = invert_matrix6(solve_op, inv);
    const double op_sym = matrix6_symmetry_max_rel(solve_op);
    const double inv_sym = invertible ? matrix6_symmetry_max_rel(inv) : 0.0;
    const double residual = invertible ? matrix6_inverse_residual_max_abs(solve_op, inv) : std::numeric_limits<double>::infinity();

    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        fs << std::scientific << std::setprecision(12)
           << omega << "," << solver_name << ",float," << proj.metric << ","
           << panels.size() << "," << proj.pairs << "," << regularization << ","
           << i << "," << j << ","
           << solve_op[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] << ","
           << inv[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] << ","
           << op_sym << "," << inv_sym << "," << residual << ","
           << (invertible ? 1 : 0) << "\n";
      }
    }
  };
  write_projection(quad);
  write_projection(centroid);
  write_projection(sym_quad);
}

void append_bie_row_multiplicity(const std::filesystem::path& path,
                                 double omega,
                                 const char* solver_name,
                                 const BEM_BVP& bvp,
                                 const std::vector<bem_frequency_domain::Id>& id_by_index) {
  const std::size_t n = id_by_index.size();
  if (n == 0 || bvp.IGIGn.size() != n)
    return;

  std::unordered_map<const BEM_DOF_Base*, std::vector<std::size_t>> by_node;
  by_node.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto* node = std::get<0>(id_by_index[i]);
    if (node)
      by_node[node].push_back(i);
  }

  struct PairStat {
    std::uint64_t groups = 0;
    std::uint64_t row_pairs = 0;
    double max_offdiag_ig_abs = 0.0;
    double max_offdiag_ign_abs = 0.0;
    double max_mutual_ig_abs = 0.0;
    double max_mutual_ign_abs = 0.0;
    std::size_t worst_offdiag_i = 0;
    std::size_t worst_offdiag_j = 0;
    std::size_t worst_offdiag_k = 0;
    std::size_t worst_mutual_i = 0;
    std::size_t worst_mutual_j = 0;
  } st;

  for (const auto& [node, ids] : by_node) {
    (void)node;
    if (ids.size() < 2)
      continue;
    ++st.groups;
    for (std::size_t a = 0; a < ids.size(); ++a) {
      const std::size_t i = ids[a];
      if (bvp.IGIGn[i].size() != n)
        return;
      for (std::size_t b = a + 1; b < ids.size(); ++b) {
        const std::size_t j = ids[b];
        if (bvp.IGIGn[j].size() != n)
          return;
        ++st.row_pairs;
        const double mutual_ig = std::max(std::abs(bvp.IGIGn[i][j][0] - bvp.IGIGn[j][i][0]),
                                          std::abs(bvp.IGIGn[i][i][0] - bvp.IGIGn[j][j][0]));
        const double mutual_ign = std::max(std::abs(bvp.IGIGn[i][j][1] - bvp.IGIGn[j][i][1]),
                                           std::abs(bvp.IGIGn[i][i][1] - bvp.IGIGn[j][j][1]));
        if (mutual_ig > st.max_mutual_ig_abs || mutual_ign > st.max_mutual_ign_abs) {
          st.max_mutual_ig_abs = std::max(st.max_mutual_ig_abs, mutual_ig);
          st.max_mutual_ign_abs = std::max(st.max_mutual_ign_abs, mutual_ign);
          st.worst_mutual_i = i;
          st.worst_mutual_j = j;
        }
        for (std::size_t k = 0; k < n; ++k) {
          if (k == i || k == j)
            continue;
          const double dig = std::abs(bvp.IGIGn[i][k][0] - bvp.IGIGn[j][k][0]);
          const double dign = std::abs(bvp.IGIGn[i][k][1] - bvp.IGIGn[j][k][1]);
          if (dig > st.max_offdiag_ig_abs || dign > st.max_offdiag_ign_abs) {
            st.max_offdiag_ig_abs = std::max(st.max_offdiag_ig_abs, dig);
            st.max_offdiag_ign_abs = std::max(st.max_offdiag_ign_abs, dign);
            st.worst_offdiag_i = i;
            st.worst_offdiag_j = j;
            st.worst_offdiag_k = k;
          }
        }
      }
    }
  }

  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header)
    fs << "omega,solver,multiple_id_groups,row_pairs,max_offdiag_ig_abs,max_offdiag_ign_abs,"
          "max_mutual_ig_abs,max_mutual_ign_abs,worst_offdiag_i,worst_offdiag_j,worst_offdiag_k,"
          "worst_mutual_i,worst_mutual_j\n";
  fs << std::scientific << std::setprecision(12)
     << omega << "," << solver_name << ","
     << st.groups << "," << st.row_pairs << ","
     << st.max_offdiag_ig_abs << "," << st.max_offdiag_ign_abs << ","
     << st.max_mutual_ig_abs << "," << st.max_mutual_ign_abs << ","
     << st.worst_offdiag_i << "," << st.worst_offdiag_j << "," << st.worst_offdiag_k << ","
     << st.worst_mutual_i << "," << st.worst_mutual_j << "\n";
}

void append_phi_jump_diagnostic(const std::filesystem::path& path,
                                double omega,
                                int dof,
                                const char* solver_name,
                                const bem_frequency_domain::Solution& sol,
                                const FaceSets& face_sets) {
  std::unordered_map<const BEM_DOF_Base*, std::vector<std::pair<networkFace*, Complex>>> by_node;
  by_node.reserve(sol.n);
  for (std::size_t i = 0; i < sol.n; ++i) {
    auto [node, f] = sol.id_by_index[i];
    if (!node || !f)
      continue;
    if (!face_sets.float_surface.contains(f))
      continue;
    by_node[node].push_back({f, sol.phi[i]});
  }

  std::uint64_t groups = 0;
  std::uint64_t pairs = 0;
  double max_abs = 0.0;
  double max_rel = 0.0;
  const BEM_DOF_Base* worst_node = nullptr;
  networkFace* worst_fi = nullptr;
  networkFace* worst_fj = nullptr;
  Complex worst_phi_i{0.0, 0.0};
  Complex worst_phi_j{0.0, 0.0};

  for (const auto& [node, vals] : by_node) {
    if (vals.size() < 2)
      continue;
    ++groups;
    for (std::size_t a = 0; a < vals.size(); ++a) {
      for (std::size_t b = a + 1; b < vals.size(); ++b) {
        ++pairs;
        const Complex diff = vals[a].second - vals[b].second;
        const double abs_diff = std::abs(diff);
        const double rel_diff = abs_diff / std::max({std::abs(vals[a].second), std::abs(vals[b].second), 1e-300});
        if (rel_diff > max_rel || (rel_diff == max_rel && abs_diff > max_abs)) {
          max_abs = abs_diff;
          max_rel = rel_diff;
          worst_node = node;
          worst_fi = vals[a].first;
          worst_fj = vals[b].first;
          worst_phi_i = vals[a].second;
          worst_phi_j = vals[b].second;
        }
      }
    }
  }

  const Tddd x = worst_node ? worst_node->getPosition() : Tddd{0.0, 0.0, 0.0};
  const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream fs(path, std::ios::app);
  if (write_header)
    fs << "omega,dof,solver,float_multi_phi_groups,pairs,max_abs,max_rel,"
          "worst_x,worst_y,worst_z,worst_face_i,worst_face_j,phi_i_re,phi_i_im,phi_j_re,phi_j_im\n";
  fs << std::scientific << std::setprecision(12)
     << omega << "," << dof << "," << solver_name << ","
     << groups << "," << pairs << "," << max_abs << "," << max_rel << ","
     << x[0] << "," << x[1] << "," << x[2] << ","
     << (worst_fi ? worst_fi->index : -1) << "," << (worst_fj ? worst_fj->index : -1) << ","
     << worst_phi_i.real() << "," << worst_phi_i.imag() << ","
     << worst_phi_j.real() << "," << worst_phi_j.imag() << "\n";
}

void write_radiation_mode_consistency(const std::filesystem::path& path,
                                      const std::unordered_set<networkFace*>& float_faces,
                                      const Tddd& ref_point,
                                      double rotation_sign,
                                      double moment_sign) {
  std::ofstream fs(path);
  fs << "dof,q_l2,m_l2,q_minus_m_l2,q_plus_m_l2,max_abs_q_minus_m,max_abs_q_plus_m,"
        "ref_x,ref_y,ref_z,rotation_sign,moment_sign\n";

  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};

  for (int dof = 0; dof < 6; ++dof) {
    long double q2 = 0.0L;
    long double m2 = 0.0L;
    long double diff2 = 0.0L;
    long double sum2 = 0.0L;
    double max_diff = 0.0;
    double max_sum = 0.0;

    for (auto* f : float_faces) {
      auto [p0, p1, p2] = f->getPoints();
      const auto& x0 = p0->X;
      const auto& x1 = p1->X;
      const auto& x2 = p2->X;
      const auto& n = f->normal;
      for (const auto& l : bary) {
        const double l0 = l[0], l1 = l[1], l2 = l[2];
        const Tddd xq = l0 * x0 + l1 * x1 + l2 * x2;
        const double q = Dot(velocity_unit_dof(dof, xq, ref_point, rotation_sign), n);
        const double m = rigid_mode_normal_component(dof, xq, n, ref_point, moment_sign);
        const double weight = w * f->area;
        q2 += static_cast<long double>(weight) * q * q;
        m2 += static_cast<long double>(weight) * m * m;
        diff2 += static_cast<long double>(weight) * (q - m) * (q - m);
        sum2 += static_cast<long double>(weight) * (q + m) * (q + m);
        max_diff = std::max(max_diff, std::abs(q - m));
        max_sum = std::max(max_sum, std::abs(q + m));
      }
    }

    fs << std::scientific << std::setprecision(12)
       << dof << ","
       << std::sqrt(static_cast<double>(q2)) << ","
       << std::sqrt(static_cast<double>(m2)) << ","
       << std::sqrt(static_cast<double>(diff2)) << ","
       << std::sqrt(static_cast<double>(sum2)) << ","
       << max_diff << "," << max_sum << ","
       << ref_point[0] << "," << ref_point[1] << "," << ref_point[2] << ","
       << rotation_sign << "," << moment_sign << "\n";
  }
}

void write_radiation_surface_summary(const std::filesystem::path& path,
                                     const Network& water,
                                     const FaceSets& face_sets,
                                     const Tddd& ref_point,
                                     double moment_sign) {
  std::ofstream fs(path);
  fs << "surface,face_count,area,nx_area,ny_area,nz_area,mx_area,my_area,mz_area,"
        "ref_x,ref_y,ref_z,moment_sign\n";

  auto write_one = [&](const std::string& name, const std::function<bool(networkFace*)>& pred) {
    std::size_t count = 0;
    double area = 0.0;
    Tddd nint{0.0, 0.0, 0.0};
    Tddd mint{0.0, 0.0, 0.0};
    for (auto* f : water.getBoundaryFaces()) {
      if (!pred(f))
        continue;
      ++count;
      area += f->area;
      nint += f->area * f->normal;
      mint += moment_sign * f->area * Cross(f->center - ref_point, f->normal);
    }
    fs << std::scientific << std::setprecision(12)
       << name << "," << count << "," << area << ","
       << nint[0] << "," << nint[1] << "," << nint[2] << ","
       << mint[0] << "," << mint[1] << "," << mint[2] << ","
       << ref_point[0] << "," << ref_point[1] << "," << ref_point[2] << ","
       << moment_sign << "\n";
  };

  write_one("float", [&](networkFace* f) { return face_sets.float_surface.contains(f); });
  write_one("free", [&](networkFace* f) { return face_sets.free_surface.contains(f); });
  write_one("wall", [&](networkFace* f) {
    return !face_sets.float_surface.contains(f) && !face_sets.free_surface.contains(f);
  });
  write_one("all", [&](networkFace*) { return true; });
}

void write_radiation_bie_dof_diagnostics(const std::filesystem::path& summary_path,
                                         const std::filesystem::path& interface_path,
                                         const std::vector<Network*>& nets,
                                         const FaceSets& face_sets,
                                         const bem_frequency_domain::BoundaryData& bc) {
  bem_frequency_domain::BoundaryStateGuard guard(nets);
  std::unordered_set<networkFace*> robin_faces;
  bem_frequency_domain::apply_face_bc_to_mesh(nets, bc, robin_faces);
  const std::size_t n = setNodeFaceIndices(nets, [&](const BEM_DOF_Base* node) {
    return node && bc.force_neumann_multiple && bc.force_neumann_multiple(*node);
  });

  auto has_robin_face = [&](const BEM_DOF_Base* node) {
    for (auto* f : node->getBoundaryFaces()) {
      if (robin_faces.contains(f))
        return true;
    }
    return false;
  };
  auto face_kind = [&](networkFace* f) -> std::string {
    if (!f)
      return "none";
    if (face_sets.float_surface.contains(f))
      return "float";
    if (face_sets.free_surface.contains(f))
      return "free";
    return "wall";
  };

  std::size_t point_ids = 0;
  std::size_t line_ids = 0;
  std::size_t active_neumann = 0;
  std::size_t active_dirichlet = 0;
  std::size_t active_robin_dirichlet = 0;
  std::size_t neumann_float = 0;
  std::size_t neumann_free = 0;
  std::size_t neumann_wall = 0;
  std::size_t interface_entities = 0;
  std::size_t interface_neumann_rows = 0;
  std::size_t interface_robin_constraints = 0;
  std::size_t interface_true_dirichlet_constraints = 0;
  std::size_t interface_missing_dirichlet = 0;

  auto scan_node = [&](const BEM_DOF_Base* node, const char* entity_type) {
    if (!node)
      return;
    const bool is_interface = node->BCInterface;
    if (is_interface)
      ++interface_entities;
    for (const auto& [f, d] : node->dofs) {
      if (d.index < 0)
        continue;
      if (std::string(entity_type) == "point")
        ++point_ids;
      else
        ++line_ids;
      if (isDirichletBieDofKey(node, f)) {
        ++active_dirichlet;
        if (has_robin_face(node))
          ++active_robin_dirichlet;
      }
      if (isNeumannBieDofKey(node, f)) {
        ++active_neumann;
        const auto kind = face_kind(f);
        if (kind == "float")
          ++neumann_float;
        else if (kind == "free")
          ++neumann_free;
        else if (kind == "wall")
          ++neumann_wall;
        if (is_interface) {
          ++interface_neumann_rows;
          const auto* d_dir = node->findActiveBieDof(nullptr);
          if (!d_dir) {
            ++interface_missing_dirichlet;
          } else if (has_robin_face(node)) {
            ++interface_robin_constraints;
          } else {
            ++interface_true_dirichlet_constraints;
          }
        }
      }
    }
  };

  std::ofstream ifs(interface_path);
  ifs << "entity_type,entity_ptr,x,y,z,boundary_faces,neumann_faces,dirichlet_faces,robin_faces,"
         "active_dirichlet_index,active_neumann_dofs\n";
  auto write_interface = [&](const BEM_DOF_Base* node, const char* entity_type) {
    if (!node || !node->BCInterface)
      return;
    std::size_t n_neu = 0;
    std::size_t n_dir = 0;
    std::size_t n_rob = 0;
    for (auto* f : node->getBoundaryFaces()) {
      if (isNeumannBoundaryState(node, f))
        ++n_neu;
      if (isDirichletBoundaryState(node, f))
        ++n_dir;
      if (robin_faces.contains(f))
        ++n_rob;
    }
    int active_neu = 0;
    for (const auto& [f, d] : node->dofs) {
      if (d.index >= 0 && isNeumannBieDofKey(node, f))
        ++active_neu;
    }
    const auto* d_dir = node->findActiveBieDof(nullptr);
    const auto& x = node->getPosition();
    ifs << entity_type << ","
        << reinterpret_cast<std::uintptr_t>(node) << ","
        << std::scientific << std::setprecision(12)
        << x[0] << "," << x[1] << "," << x[2] << ","
        << node->getBoundaryFaces().size() << ","
        << n_neu << "," << n_dir << "," << n_rob << ","
        << (d_dir ? d_dir->index : -1) << ","
        << active_neu << "\n";
  };

  for (auto* net : nets) {
    if (!net)
      continue;
    for (auto* p : net->getBoundaryPoints()) {
      scan_node(p, "point");
      write_interface(p, "point");
    }
    for (auto* l : net->getBoundaryLines()) {
      scan_node(l, "line");
      write_interface(l, "line");
    }
  }

  std::ofstream sfs(summary_path);
  sfs << "key,value\n";
  sfs << "active_ids," << n << "\n";
  sfs << "point_ids," << point_ids << "\n";
  sfs << "line_ids," << line_ids << "\n";
  sfs << "active_neumann," << active_neumann << "\n";
  sfs << "active_dirichlet," << active_dirichlet << "\n";
  sfs << "active_robin_dirichlet," << active_robin_dirichlet << "\n";
  sfs << "neumann_float," << neumann_float << "\n";
  sfs << "neumann_free," << neumann_free << "\n";
  sfs << "neumann_wall," << neumann_wall << "\n";
  sfs << "robin_faces," << robin_faces.size() << "\n";
  sfs << "interface_entities," << interface_entities << "\n";
  sfs << "interface_neumann_rows," << interface_neumann_rows << "\n";
  sfs << "interface_robin_constraints," << interface_robin_constraints << "\n";
  sfs << "interface_true_dirichlet_constraints," << interface_true_dirichlet_constraints << "\n";
  sfs << "interface_missing_dirichlet," << interface_missing_dirichlet << "\n";

  guard.restore();
}

void append_radiation_bc_residual(std::ostream& fs,
                                  double omega,
                                  int dof_col,
                                  const char* solver_name,
                                  const bem_frequency_domain::LinearSolution& sol,
                                  const FaceSets& face_sets,
                                  const Tddd& ref_point,
                                  double rotation_sign,
                                  const std::function<Complex(const BEM_DOF_Base&)>& free_surface_kappa) {
  struct Accum {
    std::size_t samples = 0;
    long double res2 = 0.0L;
    long double expected2 = 0.0L;
    long double phin2 = 0.0L;
    double max_abs_res = 0.0;
  };

  // Groups:
  // 0 float_all, 1 float_interface_points, 2 float_noninterface_points,
  // 3 free_all, 4 free_interface_points, 5 wall_all, 6 wall_interface_points.
  std::array<Accum, 7> acc{};
  auto surface_index = [&](const networkFace* f) -> int {
    if (face_sets.float_surface.contains(const_cast<networkFace*>(f)))
      return 0;
    if (face_sets.free_surface.contains(const_cast<networkFace*>(f)))
      return 1;
    return 2;
  };
  auto surface_name = [](int i) -> const char* {
    switch (i) {
    case 0:
      return "float";
    case 1:
      return "float_interface";
    case 2:
      return "float_noninterface";
    case 3:
      return "free";
    case 4:
      return "free_interface";
    case 5:
      return "wall";
    case 6:
      return "wall_interface";
    default:
      return "unknown";
    }
  };
  auto add_sample = [&](int group, const Complex& res, const Complex& expected, const Complex& phin) {
    auto& a = acc[static_cast<std::size_t>(group)];
    ++a.samples;
    a.res2 += std::norm(res);
    a.expected2 += std::norm(expected);
    a.phin2 += std::norm(phin);
    a.max_abs_res = std::max(a.max_abs_res, std::abs(res));
  };
  auto is_geometric_interface_point = [&](const networkPoint* p) {
    if (!p)
      return false;
    bool touches_float = false;
    bool touches_free = false;
    bool touches_wall = false;
    for (auto* adj : p->getBoundaryFaces()) {
      if (face_sets.float_surface.contains(adj))
        touches_float = true;
      else if (face_sets.free_surface.contains(adj))
        touches_free = true;
      else
        touches_wall = true;
    }
    return (touches_float && touches_free) || (touches_free && touches_wall);
  };

  for (const auto& [f, field] : sol.face_field) {
    if (!f)
      continue;
    auto [p0, p1, p2] = f->getPoints();
    const std::array<networkPoint*, 3> points{p0, p1, p2};
    const int si = surface_index(f);
    for (std::size_t k = 0; k < 3; ++k) {
      Complex expected{0.0, 0.0};
      if (si == 0) {
        expected = Complex{Dot(velocity_unit_dof(dof_col, points[k]->X, ref_point, rotation_sign), f->normal), 0.0};
      } else if (si == 1) {
        expected = free_surface_kappa(*points[k]) * field.phi[k];
      }
      const Complex res = field.phin[k] - expected;
      if (si == 0) {
        add_sample(0, res, expected, field.phin[k]);
        add_sample(is_geometric_interface_point(points[k]) ? 1 : 2, res, expected, field.phin[k]);
      } else if (si == 1) {
        add_sample(3, res, expected, field.phin[k]);
        if (is_geometric_interface_point(points[k]))
          add_sample(4, res, expected, field.phin[k]);
      } else {
        add_sample(5, res, expected, field.phin[k]);
        if (is_geometric_interface_point(points[k]))
          add_sample(6, res, expected, field.phin[k]);
      }
    }
  }

  for (int i = 0; i < 7; ++i) {
    const auto& a = acc[static_cast<std::size_t>(i)];
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << "," << surface_name(i) << ","
       << a.samples << ","
       << std::sqrt(static_cast<double>(a.res2)) << ","
       << a.max_abs_res << ","
       << std::sqrt(static_cast<double>(a.expected2)) << ","
       << std::sqrt(static_cast<double>(a.phin2)) << ","
       << ref_point[0] << "," << ref_point[1] << "," << ref_point[2] << ","
       << rotation_sign << "\n";
  }
}

struct FreeSurfaceAnimationPoint {
  Tddd x{};
  Complex phi_total{0.0, 0.0};
  Complex eta_total_hat{0.0, 0.0};
  double mu = 0.0;
  std::vector<std::tuple<std::string, Complex, Complex>> components;
};

Complex eta_hat_from_free_surface_phi(const Complex& phi, double mu, double omega, double gravity) {
  return -(Complex{mu, -omega} / gravity) * phi;
}

void write_free_surface_eta_vtu(const std::filesystem::path& path,
                                const std::vector<FreeSurfaceAnimationPoint>& points,
                                const std::vector<std::array<int, 3>>& cells,
                                double omega,
                                double time,
                                double displacement_scale) {
  const double phase = omega * time;
  const Complex osc{std::cos(phase), -std::sin(phase)};

  std::ofstream fs(path);
  if (!fs)
    throw std::runtime_error("write_free_surface_eta_vtu: cannot open " + path.string());

  fs << std::scientific << std::setprecision(12);
  fs << "<?xml version=\"1.0\"?>\n";
  fs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  fs << "  <UnstructuredGrid>\n";
  fs << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";

  fs << "      <Points>\n";
  fs << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto& p : points) {
    const double eta_t = (p.eta_total_hat * osc).real();
    fs << "          " << p.x[0] << " " << p.x[1] << " " << (p.x[2] + displacement_scale * eta_t) << "\n";
  }
  fs << "        </DataArray>\n";
  fs << "      </Points>\n";

  fs << "      <PointData Scalars=\"eta_total\">\n";
  auto write_scalar = [&](const char* name, const std::function<double(const FreeSurfaceAnimationPoint&)>& value) {
    fs << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
    for (const auto& p : points)
      fs << "          " << value(p) << "\n";
    fs << "        </DataArray>\n";
  };
  write_scalar("eta", [&](const FreeSurfaceAnimationPoint& p) { return (p.eta_total_hat * osc).real(); });
  write_scalar("eta_total", [&](const FreeSurfaceAnimationPoint& p) { return (p.eta_total_hat * osc).real(); });
  write_scalar("eta_total_hat_re", [&](const FreeSurfaceAnimationPoint& p) { return p.eta_total_hat.real(); });
  write_scalar("eta_total_hat_im", [&](const FreeSurfaceAnimationPoint& p) { return p.eta_total_hat.imag(); });
  write_scalar("eta_total_abs", [&](const FreeSurfaceAnimationPoint& p) { return std::abs(p.eta_total_hat); });
  write_scalar("eta_total_phase", [&](const FreeSurfaceAnimationPoint& p) { return std::atan2(p.eta_total_hat.imag(), p.eta_total_hat.real()); });
  write_scalar("phi_total_re", [&](const FreeSurfaceAnimationPoint& p) { return p.phi_total.real(); });
  write_scalar("phi_total_im", [&](const FreeSurfaceAnimationPoint& p) { return p.phi_total.imag(); });
  write_scalar("phi_total_abs", [&](const FreeSurfaceAnimationPoint& p) { return std::abs(p.phi_total); });
  if (!points.empty()) {
    for (std::size_t ci = 0; ci < points.front().components.size(); ++ci) {
      const std::string name = std::get<0>(points.front().components[ci]);
      const std::string eta_name = "eta_" + name;
      const std::string eta_re_name = "eta_" + name + "_hat_re";
      const std::string eta_im_name = "eta_" + name + "_hat_im";
      const std::string eta_abs_name = "eta_" + name + "_abs";
      const std::string eta_phase_name = "eta_" + name + "_phase";
      const std::string phi_re_name = "phi_" + name + "_re";
      const std::string phi_im_name = "phi_" + name + "_im";
      const std::string phi_abs_name = "phi_" + name + "_abs";
      write_scalar(eta_name.c_str(), [&, ci](const FreeSurfaceAnimationPoint& p) {
        return (std::get<2>(p.components[ci]) * osc).real();
      });
      write_scalar(eta_re_name.c_str(), [ci](const FreeSurfaceAnimationPoint& p) {
        return std::get<2>(p.components[ci]).real();
      });
      write_scalar(eta_im_name.c_str(), [ci](const FreeSurfaceAnimationPoint& p) {
        return std::get<2>(p.components[ci]).imag();
      });
      write_scalar(eta_abs_name.c_str(), [ci](const FreeSurfaceAnimationPoint& p) {
        return std::abs(std::get<2>(p.components[ci]));
      });
      write_scalar(eta_phase_name.c_str(), [ci](const FreeSurfaceAnimationPoint& p) {
        const auto& eta = std::get<2>(p.components[ci]);
        return std::atan2(eta.imag(), eta.real());
      });
      write_scalar(phi_re_name.c_str(), [ci](const FreeSurfaceAnimationPoint& p) {
        return std::get<1>(p.components[ci]).real();
      });
      write_scalar(phi_im_name.c_str(), [ci](const FreeSurfaceAnimationPoint& p) {
        return std::get<1>(p.components[ci]).imag();
      });
      write_scalar(phi_abs_name.c_str(), [ci](const FreeSurfaceAnimationPoint& p) {
        return std::abs(std::get<1>(p.components[ci]));
      });
    }
  }
  write_scalar("mu", [&](const FreeSurfaceAnimationPoint& p) { return p.mu; });
  fs << "      </PointData>\n";

  fs << "      <Cells>\n";
  fs << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto& c : cells)
    fs << "          " << c[0] << " " << c[1] << " " << c[2] << "\n";
  fs << "        </DataArray>\n";
  fs << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < cells.size(); ++i)
    fs << "          " << (3 * (i + 1)) << "\n";
  fs << "        </DataArray>\n";
  fs << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < cells.size(); ++i)
    fs << "          5\n";
  fs << "        </DataArray>\n";
  fs << "      </Cells>\n";
  fs << "    </Piece>\n";
  fs << "  </UnstructuredGrid>\n";
  fs << "</VTKFile>\n";
}

void write_free_surface_eta_animation_components(
    const std::filesystem::path& root,
    const std::string& label,
    double omega,
    const bem_frequency_domain::LinearSolution& geometry_sol,
    const std::vector<std::pair<std::string, const bem_frequency_domain::LinearSolution*>>& component_solutions,
    const std::function<double(const Tddd&)>& mu_at,
    double gravity,
    int frames,
    double displacement_scale) {
  if (frames <= 0)
    return;
  if (!(omega > 0.0) || !(gravity > 0.0))
    throw std::runtime_error("write_free_surface_eta_animation_components: omega and gravity must be positive");
  if (!std::isfinite(displacement_scale))
    throw std::runtime_error("write_free_surface_eta_animation_components: displacement scale must be finite");

  std::filesystem::create_directories(root);

  std::vector<std::pair<const networkFace*, const bem_frequency_domain::FaceField*>> faces;
  faces.reserve(geometry_sol.face_field.size());
  for (const auto& [f, field] : geometry_sol.face_field) {
    if (f)
      faces.push_back({f, &field});
  }
  std::sort(faces.begin(), faces.end(), [](const auto& a, const auto& b) {
    return reinterpret_cast<std::uintptr_t>(a.first) < reinterpret_cast<std::uintptr_t>(b.first);
  });

  std::vector<FreeSurfaceAnimationPoint> points;
  std::vector<std::array<int, 3>> cells;
  points.reserve(3 * faces.size());
  cells.reserve(faces.size());

  std::vector<std::pair<std::string, const bem_frequency_domain::LinearSolution*>> valid_components;
  valid_components.reserve(component_solutions.size());
  for (const auto& comp : component_solutions) {
    if (!comp.first.empty() && comp.second)
      valid_components.push_back(comp);
  }

  for (const auto& [f, field] : faces) {
    auto [p0, p1, p2] = f->getPoints();
    const std::array<networkPoint*, 3> face_points{p0, p1, p2};
    std::array<int, 3> cell{};
    for (int k = 0; k < 3; ++k) {
      const auto* p = face_points[static_cast<std::size_t>(k)];
      const Tddd x = p ? p->X : Tddd{0.0, 0.0, 0.0};
      const double mu = mu_at ? mu_at(x) : 0.0;
      const Complex phi_total = field->phi[static_cast<std::size_t>(k)];
      const Complex eta_total_hat = eta_hat_from_free_surface_phi(phi_total, mu, omega, gravity);
      std::vector<std::tuple<std::string, Complex, Complex>> components;
      components.reserve(valid_components.size());
      for (const auto& [name, comp_sol] : valid_components) {
        Complex phi_component{0.0, 0.0};
        const auto it = comp_sol->face_field.find(const_cast<networkFace*>(f));
        if (it != comp_sol->face_field.end())
          phi_component = it->second.phi[static_cast<std::size_t>(k)];
        components.push_back({name, phi_component, eta_hat_from_free_surface_phi(phi_component, mu, omega, gravity)});
      }
      cell[static_cast<std::size_t>(k)] = static_cast<int>(points.size());
      points.push_back(FreeSurfaceAnimationPoint{x, phi_total, eta_total_hat, mu, std::move(components)});
    }
    cells.push_back(cell);
  }

  const double period = 2.0 * M_PI / omega;
  const std::string prefix = label + "_eta_total";
  std::vector<std::string> frame_files;
  frame_files.reserve(static_cast<std::size_t>(frames));

  for (int frame = 0; frame < frames; ++frame) {
    const double time = period * static_cast<double>(frame) / static_cast<double>(frames);
    std::ostringstream name;
    name << prefix << "_" << std::setw(4) << std::setfill('0') << frame << ".vtu";
    frame_files.push_back(name.str());
    write_free_surface_eta_vtu(root / frame_files.back(), points, cells, omega, time, displacement_scale);
  }

  const auto pvd_path = root / (prefix + ".pvd");
  std::ofstream pvd(pvd_path);
  if (!pvd)
    throw std::runtime_error("write_free_surface_eta_animation_components: cannot open " + pvd_path.string());
  pvd << std::scientific << std::setprecision(12);
  pvd << "<?xml version=\"1.0\"?>\n";
  pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  pvd << "  <Collection>\n";
  for (int frame = 0; frame < frames; ++frame) {
    const double time = period * static_cast<double>(frame) / static_cast<double>(frames);
    pvd << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\""
        << frame_files[static_cast<std::size_t>(frame)] << "\"/>\n";
  }
  pvd << "  </Collection>\n";
  pvd << "</VTKFile>\n";

  const auto meta_path = root / (prefix + "_metadata.csv");
  std::ofstream meta(meta_path);
  meta << "key,value\n"
       << "label," << label << "\n"
       << "omega," << std::scientific << std::setprecision(12) << omega << "\n"
       << "period," << period << "\n"
       << "frames," << frames << "\n"
       << "points," << points.size() << "\n"
       << "cells," << cells.size() << "\n"
       << "displacement_scale," << displacement_scale << "\n"
       << "geometry,eta_total\n"
       << "eta_relation,eta_hat=-(mu-i*omega)*phi/g\n";
  meta << "components";
  for (const auto& [name, _] : valid_components)
    meta << "," << name;
  meta << "\n";

  std::cout << "wrote free-surface eta animation: " << pvd_path << std::endl;
}

CVector6 solve_rao_displacement(double omega,
                                const std::vector<int>& dofs,
                                const std::array<double, 6>& mass_diag,
                                const Matrix6& restoring,
                                const Matrix6& added_mass,
                                const Matrix6& damping,
                                const CVector6& excitation,
                                CVector6& residual_out) {
  if (dofs.empty())
    throw std::runtime_error("solve_rao_displacement: no active dofs");

  const Complex I(0.0, 1.0);
  const int nd = static_cast<int>(dofs.size());
  std::vector<Complex> A_col_major(static_cast<std::size_t>(nd * nd), Complex{0.0, 0.0});
  std::vector<Complex> rhs(static_cast<std::size_t>(nd), Complex{0.0, 0.0});

  for (int r = 0; r < nd; ++r) {
    const int row = dofs[static_cast<std::size_t>(r)];
    rhs[static_cast<std::size_t>(r)] = excitation[static_cast<std::size_t>(row)];
    for (int c = 0; c < nd; ++c) {
      const int col = dofs[static_cast<std::size_t>(c)];
      const double m = (row == col) ? mass_diag[static_cast<std::size_t>(row)] : 0.0;
      const Complex d = Complex{restoring[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]
                                    - omega * omega * (m + added_mass[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]),
                                0.0}
                        - I * omega * damping[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
      A_col_major[static_cast<std::size_t>(r + c * nd)] = d;
    }
  }

  bem_frequency_domain::lapack_zlu lu(nd, A_col_major);
  lu.solve_in_place(rhs);

  CVector6 xi{};
  xi.fill(Complex{0.0, 0.0});
  for (int i = 0; i < nd; ++i)
    xi[static_cast<std::size_t>(dofs[static_cast<std::size_t>(i)])] = rhs[static_cast<std::size_t>(i)];

  residual_out.fill(Complex{0.0, 0.0});
  for (int r = 0; r < nd; ++r) {
    const int row = dofs[static_cast<std::size_t>(r)];
    Complex lhs{0.0, 0.0};
    for (int c = 0; c < nd; ++c) {
      const int col = dofs[static_cast<std::size_t>(c)];
      const double m = (row == col) ? mass_diag[static_cast<std::size_t>(row)] : 0.0;
      const Complex d = Complex{restoring[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]
                                    - omega * omega * (m + added_mass[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]),
                                0.0}
                        - I * omega * damping[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
      lhs += d * xi[static_cast<std::size_t>(col)];
    }
    residual_out[static_cast<std::size_t>(row)] = lhs - excitation[static_cast<std::size_t>(row)];
  }

  return xi;
}

} // namespace

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cerr << "usage: ./main_freq_domain <input_dir> [--omega w1 w2 ...] [--dofs i j ...]\n"
                 "  [--qtf] [--qtf-wave] [--qtf-unit-wave] [--qtf-newman] [--qtf-check] [--qtf-full] [--rao] [--wave-excitation-group-probe]\n"
                 "  [--radiation-complex-lu|--radiation-complex-gmres|--radiation-complex-fmm-gmres|--freq-complex-fmm-gmres] [--complex-gmres-tol tol] [--complex-gmres-max-iter n] [--complex-gmres-restart n]\n"
                 "  [--complex-gmres-preconditioner diagonal|fmm_near_diagonal|row_norm|none]\n"
                 "  [--freq-sponge mu_max r_start length [order]]\n"
                 "  [--freq-wall-sponge mu_max width [order]]\n"
                 "  [--freq-absorber-sponge-by-wavelength mu_max n_lambda [order]]\n"
                 "  [--freq-absorber-wall-max-sponge absorber_mu n_lambda absorber_order wall_mu wall_width wall_order]\n"
	                 "  [--freq-sponge-profile power|smoothstep|smoothstep_power|floor_power|floor_smoothstep_power]\n"
	                 "  [--freq-sponge-floor eps]\n"
	                 "  [--fmm-max-level n] [--fmm-bucket-max-points n] [--fmm-mac-theta theta] [--nearfield-mode mode]\n"
	                 "  [--radiation-fmm-no-coordinate-scaling|--freq-fmm-no-coordinate-scaling]\n"
	                 "  [--lu-dunavant-quadrature far_points near_points] (0 keeps legacy GW for that range)\n"
                 "  [--radiation-debug] [--radiation-ref-point x y z]\n"
		                 "  [--radiation-rotation-sign +/-1] [--radiation-moment-sign +/-1]\n"
		                 "  [--radiation-postprocess-probe]\n"
		                 "  [--radiation-pitch-probe]\n"
			                 "  [--radiation-coupling-probe]\n"
				                 "  [--radiation-fmm-matvec-compare]\n"
				                 "  [--wave-heading-deg heading]\n"
				                 "  [--freq-open-boundary none|sommerfeld] [--freq-open-boundary-kappa-scale s]\n"
				                 "  [--freq-open-boundary-wave-number real|complex_sponge] [--freq-open-boundary-mu-scale s]\n"
				                 "  [--freq-open-boundary-curvature none|cylindrical] [--freq-open-boundary-radial-projection]\n"
				                 "  [--freq-output-free-surface-animation frames] [--freq-free-surface-animation-scale s]\n"
			                 "  [--radiation-face-neumann] [--radiation-symmetrized-lu]\n"
                 "  [--radiation-interface-bie-wall-free] [--radiation-interface-bie-float-free]\n"
                 "  [--radiation-condense-interface]\n"
                 "  [--radiation-interface-row-scale unit|area|lumped_area]\n"
                 "  [--radiation-interface-robin-adjoint-probe]\n"
                 "  [--freq-bc-contact-check] [--freq-bc-contact-check-only]\n"
                 "  [--export-wetted-body-mesh]\n";
    return 2;
  }

  SimulationSettings setting(argv[1], SimulationSettings::DomainMode::Frequency);
  use_linear_element = setting.bem.element.linear;
  use_pseudo_quadratic_element = setting.bem.element.pseudo_quadratic;
  use_quadratic_linear_hybrid = setting.bem.element.quadratic_linear_hybrid;
  use_true_quadratic_element = setting.bem.element.true_quadratic || use_quadratic_linear_hybrid;

  solver_type = "GMRES";
  preconditioner_type = setting.bem.solver.preconditioner_type;
  ilu_neighborhood_type = setting.bem.solver.ilu_neighborhood_type;
  ilu_kring_num = setting.bem.solver.ilu_kring_num;
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
  enable_pressure_detachment = setting.bem.solver.enable_pressure_detachment;
  detachment_consecutive_steps = setting.bem.solver.detachment_consecutive_steps;
	  nearfield_mode = setting.bem.solver.nearfield_mode;
	  fmm_max_level = setting.bem.solver.fmm_max_level;
	  fmm_bucket_max_points = setting.bem.solver.fmm_bucket_max_points;
	  milu_omega = setting.bem.solver.milu_omega;
	  g_p2m_quadrature_points = setting.bem.solver.p2m_quadrature_points;
	  g_mac_theta = setting.bem.solver.mac_theta;
	  fmm_max_level = parse_option_int(argc, argv, "--fmm-max-level", fmm_max_level);
	  fmm_bucket_max_points = parse_option_int(argc, argv, "--fmm-bucket-max-points", fmm_bucket_max_points);
	  g_mac_theta = parse_option_double(argc, argv, "--fmm-mac-theta", g_mac_theta);
	  nearfield_mode = parse_option_string(argc, argv, "--nearfield-mode", nearfield_mode);
	  if (fmm_max_level < 0)
	    throw std::runtime_error("--fmm-max-level must be non-negative");
	  if (fmm_bucket_max_points <= 0)
	    throw std::runtime_error("--fmm-bucket-max-points must be positive");
	  if (!(std::isfinite(g_mac_theta) && g_mac_theta > 0.0))
	    throw std::runtime_error("--fmm-mac-theta must be positive and finite");

  const auto omegas = parse_omegas(argc, argv, setting.frequency.omegas);
  std::vector<int> dofs = parse_dofs(argc, argv, setting.frequency.dofs);
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--lu-dunavant-quadrature") {
      if (i + 2 >= argc || std::string(argv[i + 1]).starts_with("--") || std::string(argv[i + 2]).starts_with("--"))
        throw std::runtime_error("--lu-dunavant-quadrature requires far_points near_points");
      g_lu_far_dunavant_points = std::stoi(argv[++i]);
      g_lu_near_dunavant_points = std::stoi(argv[++i]);
      if (g_lu_far_dunavant_points < 0 || g_lu_near_dunavant_points < 0)
        throw std::runtime_error("--lu-dunavant-quadrature values must be non-negative; use 0 to keep legacy GW");
      if (g_lu_far_dunavant_points > 0)
        (void)getDunavantP2MRule(g_lu_far_dunavant_points);
      if (g_lu_near_dunavant_points > 0)
        (void)getDunavantP2MRule(g_lu_near_dunavant_points);
    }
  }
  const bool dofs_explicit = has_flag(argc, argv, "--dofs") || !setting.frequency.dofs.empty();
  const bool qtf_radiation = has_flag(argc, argv, "--qtf");
  const bool qtf_wave = has_flag(argc, argv, "--qtf-wave");
  const bool qtf_unit_wave = has_flag(argc, argv, "--qtf-unit-wave");
  const bool qtf_newman = has_flag(argc, argv, "--qtf-newman");
  const bool qtf_check = has_flag(argc, argv, "--qtf-check");
  const bool qtf_full = has_flag(argc, argv, "--qtf-full") || qtf_newman || qtf_check;
  const bool qtf_symmetry = !qtf_full;
  const bool wave_excitation_group_probe = has_flag(argc, argv, "--wave-excitation-group-probe");
	  const bool rao = has_flag(argc, argv, "--rao");
	  const bool radiation_complex_lu = has_flag(argc, argv, "--radiation-complex-lu");
	  const bool radiation_complex_gmres = has_flag(argc, argv, "--radiation-complex-gmres");
	  const bool freq_complex_fmm_gmres = has_flag(argc, argv, "--freq-complex-fmm-gmres") ||
	                                      has_flag(argc, argv, "--radiation-complex-fmm-gmres");
	  const bool radiation_complex_fmm_gmres = freq_complex_fmm_gmres;
	  const bool wave_complex_fmm_gmres = freq_complex_fmm_gmres;
	  const bool fmm_true_quad_straight_geometry_probe = has_flag(argc, argv, "--fmm-true-quad-straight-geometry-probe");
	  const bool radiation_complex_solver = radiation_complex_lu || radiation_complex_gmres || radiation_complex_fmm_gmres;
		  const bool radiation_debug = has_flag(argc, argv, "--radiation-debug");
		  const bool radiation_postprocess_probe = has_flag(argc, argv, "--radiation-postprocess-probe");
			  const bool radiation_pitch_probe = has_flag(argc, argv, "--radiation-pitch-probe");
			  const bool radiation_coupling_probe = has_flag(argc, argv, "--radiation-coupling-probe");
			  const bool radiation_fmm_matvec_compare = has_flag(argc, argv, "--radiation-fmm-matvec-compare");
			  const bool freq_fmm_no_coordinate_scaling = has_flag(argc, argv, "--freq-fmm-no-coordinate-scaling") ||
			                                             has_flag(argc, argv, "--radiation-fmm-no-coordinate-scaling");
			  const bool radiation_fmm_no_coordinate_scaling = freq_fmm_no_coordinate_scaling;
			  const bool wave_fmm_no_coordinate_scaling = freq_fmm_no_coordinate_scaling;
	  const bool radiation_face_neumann = has_flag(argc, argv, "--radiation-face-neumann");
	  const bool radiation_symmetrized_lu = has_flag(argc, argv, "--radiation-symmetrized-lu");
	  const bool radiation_free_surface_dirichlet_zero = has_flag(argc, argv, "--radiation-free-surface-dirichlet-zero");
	  const bool radiation_robin_phi_unknown_legacy = has_flag(argc, argv, "--radiation-robin-phi-unknown");
	  const std::string radiation_robin_unknown =
	      parse_option_string(argc,
	                          argv,
	                          "--radiation-robin-unknown",
	                          radiation_robin_phi_unknown_legacy ? "phi" : "phin");
	  if (radiation_robin_phi_unknown_legacy && radiation_robin_unknown != "phi")
	    throw std::runtime_error("--radiation-robin-phi-unknown conflicts with --radiation-robin-unknown " +
	                             radiation_robin_unknown);
	  if (radiation_robin_unknown == "mixed")
	    throw std::runtime_error("--radiation-robin-unknown mixed is not implemented safely yet; use phin or phi");
	  if (radiation_robin_unknown != "phin" && radiation_robin_unknown != "phi")
	    throw std::runtime_error("--radiation-robin-unknown must be one of phin, phi");
	  const auto radiation_robin_unknown_policy =
	      (radiation_robin_unknown == "phi") ? bem_frequency_domain::RobinUnknownPolicy::Phi
	                                         : bem_frequency_domain::RobinUnknownPolicy::Phin;
	  const bool radiation_robin_phi_unknown =
	      radiation_robin_unknown_policy == bem_frequency_domain::RobinUnknownPolicy::Phi;
	  const bool radiation_interface_bie_wall_free = has_flag(argc, argv, "--radiation-interface-bie-wall-free");
	  const bool radiation_interface_bie_float_free = has_flag(argc, argv, "--radiation-interface-bie-float-free");
	  const bool radiation_condense_interface = has_flag(argc, argv, "--radiation-condense-interface");
		  const bool radiation_interface_robin_adjoint_probe = has_flag(argc, argv, "--radiation-interface-robin-adjoint-probe");
		  const bool freq_bc_contact_check_only = has_flag(argc, argv, "--freq-bc-contact-check-only");
		  const bool freq_bc_contact_check = has_flag(argc, argv, "--freq-bc-contact-check") || freq_bc_contact_check_only;
	  const bool export_wetted_body_mesh = has_flag(argc, argv, "--export-wetted-body-mesh");
	  const std::string radiation_interface_row_scale = parse_option_string(argc, argv, "--radiation-interface-row-scale", "default");
	  const double radiation_rotation_sign = parse_option_double(argc, argv, "--radiation-rotation-sign", 1.0);
	  const double radiation_moment_sign = parse_option_double(argc, argv, "--radiation-moment-sign", 1.0);
	  const double radiation_robin_kappa_scale = parse_option_double(argc, argv, "--radiation-robin-kappa-scale", 1.0);
	  const double complex_gmres_tol = parse_option_double(argc, argv, "--complex-gmres-tol", solver_tol);
	  const int complex_gmres_max_iter = parse_option_int(argc, argv, "--complex-gmres-max-iter", solver_max_iter);
	  const int complex_gmres_restart = parse_option_int(argc, argv, "--complex-gmres-restart", solver_restart);
			  const std::string complex_gmres_preconditioner = parse_option_string(argc, argv, "--complex-gmres-preconditioner", "diagonal");
			  const std::optional<double> wave_heading_deg = parse_optional_double(argc, argv, "--wave-heading-deg");
			  if (has_flag(argc, argv, "--freq-wave-relax-incident") || has_flag(argc, argv, "--freq-relax-incident"))
			    throw std::runtime_error("--freq-wave-relax-incident was removed; use --freq-open-boundary sommerfeld for the outgoing-boundary probe");
			  const std::string freq_open_boundary_mode =
			      has_flag(argc, argv, "--freq-open-boundary-sommerfeld")
			          ? "sommerfeld"
			          : parse_option_string(argc, argv, "--freq-open-boundary", "none");
			  const bool freq_open_boundary_sommerfeld = (freq_open_boundary_mode == "sommerfeld");
			  const double freq_open_boundary_kappa_scale = parse_option_double(argc, argv, "--freq-open-boundary-kappa-scale", 1.0);
				  const std::string freq_open_boundary_wave_number =
				      parse_option_string(argc, argv, "--freq-open-boundary-wave-number", "real");
				  const double freq_open_boundary_mu_scale = parse_option_double(argc, argv, "--freq-open-boundary-mu-scale", 1.0);
				  const std::string freq_open_boundary_curvature =
				      parse_option_string(argc, argv, "--freq-open-boundary-curvature", "none");
				  const bool freq_open_boundary_radial_projection =
				      has_flag(argc, argv, "--freq-open-boundary-radial-projection");
				  const int freq_free_surface_animation_frames =
				      parse_option_int(argc, argv, "--freq-output-free-surface-animation", 0);
				  const double freq_free_surface_animation_scale =
				      parse_option_double(argc, argv, "--freq-free-surface-animation-scale", 1.0);
				  require_sign_option(radiation_rotation_sign, "--radiation-rotation-sign");
	  require_sign_option(radiation_moment_sign, "--radiation-moment-sign");
	  if (!std::isfinite(radiation_robin_kappa_scale) || radiation_robin_kappa_scale <= 0.0)
	    throw std::runtime_error("--radiation-robin-kappa-scale must be a positive finite value");
	  if (!(radiation_interface_row_scale == "default" ||
	        radiation_interface_row_scale == "unit" ||
	        radiation_interface_row_scale == "area" ||
	        radiation_interface_row_scale == "lumped_area"))
	    throw std::runtime_error("--radiation-interface-row-scale must be one of default, unit, area, lumped_area");
	  const bool compute_wave = qtf_wave || rao;
	  const int complex_solver_count = (radiation_complex_lu ? 1 : 0) + (radiation_complex_gmres ? 1 : 0) + (radiation_complex_fmm_gmres ? 1 : 0);
	  if (complex_solver_count > 1)
	    throw std::runtime_error("frequency_domain_main: choose only one of --radiation-complex-lu, --radiation-complex-gmres, and --radiation-complex-fmm-gmres");
	  if (radiation_symmetrized_lu && !radiation_complex_lu)
	    throw std::runtime_error("frequency_domain_main: --radiation-symmetrized-lu requires --radiation-complex-lu");
	  if ((radiation_complex_gmres || radiation_complex_fmm_gmres) && (complex_gmres_tol <= 0.0 || !std::isfinite(complex_gmres_tol)))
	    throw std::runtime_error("--complex-gmres-tol must be positive and finite");
	  if ((radiation_complex_gmres || radiation_complex_fmm_gmres) && (complex_gmres_max_iter <= 0 || complex_gmres_restart <= 0))
	    throw std::runtime_error("--complex-gmres-max-iter and --complex-gmres-restart must be positive");
		  if ((radiation_complex_gmres || radiation_complex_fmm_gmres) &&
		      complex_gmres_preconditioner != "diagonal" &&
	      complex_gmres_preconditioner != "fmm_near_diagonal" &&
	      complex_gmres_preconditioner != "row_norm" &&
	      complex_gmres_preconditioner != "none")
	    throw std::runtime_error("--complex-gmres-preconditioner must be one of: diagonal, fmm_near_diagonal, row_norm, none");
		  if (radiation_complex_gmres && complex_gmres_preconditioner == "fmm_near_diagonal")
		    throw std::runtime_error("--complex-gmres-preconditioner fmm_near_diagonal requires --freq-complex-fmm-gmres");
		  if (!(freq_open_boundary_mode == "none" || freq_open_boundary_mode == "sommerfeld"))
		    throw std::runtime_error("--freq-open-boundary must be one of: none, sommerfeld");
		  if (!(std::isfinite(freq_open_boundary_kappa_scale) && freq_open_boundary_kappa_scale > 0.0))
		    throw std::runtime_error("--freq-open-boundary-kappa-scale must be positive and finite");
			  if (!(freq_open_boundary_wave_number == "real" || freq_open_boundary_wave_number == "complex_sponge"))
			    throw std::runtime_error("--freq-open-boundary-wave-number must be one of: real, complex_sponge");
			  if (!(std::isfinite(freq_open_boundary_mu_scale) && freq_open_boundary_mu_scale >= 0.0))
			    throw std::runtime_error("--freq-open-boundary-mu-scale must be non-negative and finite");
			  if (!(freq_open_boundary_curvature == "none" || freq_open_boundary_curvature == "cylindrical"))
			    throw std::runtime_error("--freq-open-boundary-curvature must be one of: none, cylindrical");
			  if (freq_free_surface_animation_frames < 0)
			    throw std::runtime_error("--freq-output-free-surface-animation must be non-negative");
			  if (!std::isfinite(freq_free_surface_animation_scale))
			    throw std::runtime_error("--freq-free-surface-animation-scale must be finite");
			  if (freq_open_boundary_sommerfeld && !radiation_complex_solver)
			    throw std::runtime_error("--freq-open-boundary sommerfeld requires a complex radiation solver; use --radiation-complex-lu for the current standard path");
			  if (radiation_fmm_matvec_compare && !radiation_complex_fmm_gmres)
			    throw std::runtime_error("--radiation-fmm-matvec-compare requires --radiation-complex-fmm-gmres");
		  if (fmm_true_quad_straight_geometry_probe)
	    setenv("BEM_FMM_TRUE_QUAD_STRAIGHT_GEOMETRY_PROBE", "1", 1);
	  else
	    unsetenv("BEM_FMM_TRUE_QUAD_STRAIGHT_GEOMETRY_PROBE");

  FrequencySpongeRuntime sponge{
      setting.frequency.sponge.enabled,
      FrequencySpongeRuntime::Mode::Radial,
      setting.frequency.sponge.r_start,
      setting.frequency.sponge.length,
      setting.frequency.sponge.mu_max,
      setting.frequency.sponge.order,
      3.0,
      "power",
      0.0,
      0.0,
      0.0,
      2};
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--freq-sponge-off") {
      sponge.enabled = false;
    } else if (arg == "--freq-sponge") {
      if (i + 3 >= argc)
        throw std::runtime_error("--freq-sponge requires mu_max r_start length [order]");
      sponge.enabled = true;
      sponge.mode = FrequencySpongeRuntime::Mode::Radial;
      sponge.mu_max = std::stod(argv[++i]);
      sponge.r_start = std::stod(argv[++i]);
      sponge.length = std::stod(argv[++i]);
      if (i + 1 < argc && !std::string(argv[i + 1]).starts_with("--"))
        sponge.order = std::stoi(argv[++i]);
    } else if (arg == "--freq-wall-sponge") {
      if (i + 2 >= argc)
        throw std::runtime_error("--freq-wall-sponge requires mu_max width [order]");
      sponge.enabled = true;
      sponge.mode = FrequencySpongeRuntime::Mode::WallDistance;
      sponge.mu_max = std::stod(argv[++i]);
      sponge.r_start = 0.0;
      sponge.length = std::stod(argv[++i]);
      if (i + 1 < argc && !std::string(argv[i + 1]).starts_with("--"))
        sponge.order = std::stoi(argv[++i]);
    } else if (arg == "--freq-absorber-sponge-by-wavelength") {
      if (i + 2 >= argc)
        throw std::runtime_error("--freq-absorber-sponge-by-wavelength requires mu_max n_lambda [order]");
      sponge.enabled = true;
      sponge.mode = FrequencySpongeRuntime::Mode::AbsorberWavelength;
      sponge.mu_max = std::stod(argv[++i]);
      sponge.r_start = 0.0;
      sponge.length = 0.0;
      sponge.n_wavelengths = std::stod(argv[++i]);
      if (i + 1 < argc && !std::string(argv[i + 1]).starts_with("--"))
        sponge.order = std::stoi(argv[++i]);
    } else if (arg == "--freq-absorber-wall-max-sponge") {
      if (i + 6 >= argc)
        throw std::runtime_error("--freq-absorber-wall-max-sponge requires absorber_mu n_lambda absorber_order wall_mu wall_width wall_order");
      sponge.enabled = true;
      sponge.mode = FrequencySpongeRuntime::Mode::AbsorberWallMax;
      sponge.mu_max = std::stod(argv[++i]);
      sponge.r_start = 0.0;
      sponge.length = 0.0;
      sponge.n_wavelengths = std::stod(argv[++i]);
      sponge.order = std::stoi(argv[++i]);
      sponge.wall_mu_max = std::stod(argv[++i]);
      sponge.wall_width = std::stod(argv[++i]);
      sponge.wall_order = std::stoi(argv[++i]);
    } else if (arg == "--freq-sponge-profile") {
      if (i + 1 >= argc || std::string(argv[i + 1]).starts_with("--"))
        throw std::runtime_error("--freq-sponge-profile requires a profile name");
      sponge.profile = argv[++i];
    } else if (arg == "--freq-sponge-floor") {
      if (i + 1 >= argc || std::string(argv[i + 1]).starts_with("--"))
        throw std::runtime_error("--freq-sponge-floor requires eps");
      sponge.floor_fraction = std::stod(argv[++i]);
    }
  }
  if (sponge.order < 1)
    throw std::runtime_error("frequency_domain_main: sponge order must be >= 1");
  if (sponge.enabled && !(sponge.mu_max >= 0.0))
    throw std::runtime_error("frequency_domain_main: enabled sponge requires mu_max >= 0");
  if (sponge.enabled && (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWavelength ||
                         sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax)) {
    if (!(sponge.n_wavelengths > 0.0) || !std::isfinite(sponge.n_wavelengths))
      throw std::runtime_error("frequency_domain_main: absorber wavelength sponge requires n_lambda > 0");
    if (setting.AbsorberObject.empty())
      throw std::runtime_error("frequency_domain_main: absorber wavelength sponge requires at least one Absorber object in settings");
    if (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax) {
      if (!(sponge.wall_mu_max >= 0.0) || !std::isfinite(sponge.wall_mu_max))
        throw std::runtime_error("frequency_domain_main: --freq-absorber-wall-max-sponge requires wall_mu >= 0");
      if (!(sponge.wall_width > 0.0) || !std::isfinite(sponge.wall_width))
        throw std::runtime_error("frequency_domain_main: --freq-absorber-wall-max-sponge requires wall_width > 0");
      if (sponge.wall_order < 1)
        throw std::runtime_error("frequency_domain_main: --freq-absorber-wall-max-sponge requires wall_order >= 1");
    }
  } else if (sponge.enabled && !(sponge.length > 0.0)) {
    throw std::runtime_error("frequency_domain_main: enabled radial/wall sponge requires length > 0");
  }
  if (!(sponge.profile == "power" ||
        sponge.profile == "smoothstep" ||
        sponge.profile == "smoothstep_power" ||
        sponge.profile == "floor_power" ||
        sponge.profile == "floor_smoothstep_power"))
    throw std::runtime_error("--freq-sponge-profile must be one of: power, smoothstep, smoothstep_power, floor_power, floor_smoothstep_power");
  if (!std::isfinite(sponge.floor_fraction) || sponge.floor_fraction < 0.0 || sponge.floor_fraction > 1.0)
    throw std::runtime_error("--freq-sponge-floor must be in [0, 1]");

  if (setting.FluidObject.empty())
    throw std::runtime_error("frequency_domain_main: no FluidObject found");

  Network* water = setting.FluidObject.front();
  if (!water)
    throw std::runtime_error("frequency_domain_main: water network is null");

  Network* float_body = nullptr;
  for (auto* rb : setting.RigidBodyObject) {
    if (!rb)
      continue;
    if (rb->getName() == "float") {
      float_body = rb;
      break;
    }
  }
  if (!float_body) {
    for (auto* rb : setting.RigidBodyObject) {
      if (!rb)
        continue;
      if (std::ranges::any_of(rb->isFixed, [](bool v) { return !v; })) {
        float_body = rb;
        break;
      }
    }
  }
  if (!float_body)
    throw std::runtime_error("frequency_domain_main: could not find a movable rigid body (float)");

  if (rao) {
    if (!dofs_explicit)
      dofs = nonfixed_dofs(*float_body);
    if (dofs.empty())
      throw std::runtime_error("frequency_domain_main: --rao requested but the selected float has no non-fixed DOFs");
    for (int dof : dofs) {
      if (float_body->isFixed[static_cast<std::size_t>(dof)])
        throw std::runtime_error("frequency_domain_main: --rao selected DOF " + std::to_string(dof) + " is fixed");
    }
    if (setting.frequency.rao_restoring.empty())
      throw std::runtime_error("frequency_domain_main: --rao requires freq_rao_restoring with 36 row-major values in settings.json");
  }

	  const double rho = setting.common.water_density;
	  const double g = setting.common.gravity;
	  const Matrix6 rao_restoring = rao ? parse_row_major_matrix6(setting.frequency.rao_restoring, "freq_rao_restoring") : make_zero_matrix6();
	  const std::array<double, 6> rao_mass_diag = rigid_body_mass_diagonal(*float_body);
	  const Tddd radiation_ref_point = parse_option_tddd(argc, argv, "--radiation-ref-point", float_body->COM);

  std::filesystem::path outdir = setting.common.output_directory / "frequency_domain";
  std::filesystem::create_directories(outdir);

	  std::unordered_map<double, Matrix6> added_mass_by_omega;
	  std::unordered_map<double, Matrix6> damping_by_omega;
	  std::unordered_map<double, Matrix6> added_mass_raw_by_omega;
	  std::unordered_map<double, Matrix6> damping_raw_by_omega;
  std::unordered_map<int, std::vector<bem_frequency_domain::LinearSolution>> qtf_solutions;
  std::map<std::pair<double, int>, bem_frequency_domain::LinearSolution> radiation_free_surface_solutions;

  // Guard & override boundary-state flags (keep this binary self-contained).
  bem_frequency_domain::BoundaryStateGuard guard(setting.FluidObject);

  water->setGeometricPropertiesForce();
  float_body->setGeometricPropertiesForce();

	  const auto face_sets = classify_faces_deepcwind(*water, *float_body);
	  if (freq_bc_contact_check || radiation_debug || export_wetted_body_mesh) {
	    const auto bc_check = write_frequency_contact_bc_check(outdir,
	                                                          setting,
	                                                          water,
	                                                          float_body,
	                                                          face_sets,
	                                                          radiation_robin_unknown_policy,
	                                                          radiation_face_neumann);
	    std::cout << "[freq_bc_contact_check] time_active_ids=" << bc_check.time_active_ids
	              << " freq_active_ids=" << bc_check.freq_active_ids
	              << " node_face_mismatches=" << bc_check.node_face_mismatches
	              << " face_mismatches=" << bc_check.face_mismatches
	              << " BCInterface_no_contact={p:" << bc_check.pre_bvp_stats.bcinterface_points_no_contact
	              << "/" << bc_check.pre_bvp_stats.bcinterface_points
	              << ",l:" << bc_check.pre_bvp_stats.bcinterface_lines_no_contact
	              << "/" << bc_check.pre_bvp_stats.bcinterface_lines << "}"
	              << std::endl;
	    if (export_wetted_body_mesh) {
	      const auto mesh_summary = export_wetted_body_mesh_from_frequency_bc(outdir, setting, water, float_body, face_sets);
	      std::cout << "[export_wetted_body_mesh] vertices=" << mesh_summary.vertices
	                << " faces=" << mesh_summary.faces
	                << " area=" << mesh_summary.area
	                << " faces_without_float_contact=" << mesh_summary.faces_without_float_contact
	                << " outdir=" << outdir << std::endl;
	    }
		    if (freq_bc_contact_check_only)
		      return 0;
		  }
	  std::ofstream ofs(outdir / "radiation_coeffs.csv");
	  ofs << "omega,dof_row,dof_col,A,B\n";
		  std::ofstream raw_coeffs_fs(outdir / "radiation_coeffs_raw.csv");
		  raw_coeffs_fs << "omega,dof_row,dof_col,A,B,solver,sponge_enabled,sponge_mu_max,sponge_r_start,sponge_length,sponge_order\n";
		  std::ofstream coeff_symmetry_fs(outdir / "radiation_coeffs_symmetry.csv");
		  coeff_symmetry_fs << "omega,dof_i,dof_j,metric,value_ij,value_ji,value_sym,antisym_abs,antisym_rel,has_ij,has_ji,solver\n";
		  std::ofstream zfs(outdir / "radiation_impedance.csv");
		  zfs << "omega,dof_row,dof_col,Z_re,Z_im,A,B,solver,sponge_enabled,sponge_mu_max,sponge_r_start,sponge_length,sponge_order\n";
		  std::ofstream reciprocity_fs(outdir / "radiation_reciprocity.csv");
		  reciprocity_fs << "omega,dof_row,dof_col,S_re,S_im,solver,sponge_enabled,sponge_mu_max,sponge_r_start,sponge_length,sponge_order,"
		                    "ref_x,ref_y,ref_z,rotation_sign,moment_sign\n";
		  std::ofstream bc_residual_fs;
		  if (radiation_debug) {
		    bc_residual_fs.open(outdir / "radiation_bc_residual.csv");
		    bc_residual_fs << "omega,dof,solver,surface,samples,res_l2,max_abs_res,expected_l2,phin_l2,"
		                      "ref_x,ref_y,ref_z,rotation_sign\n";
		  }
		  std::unordered_set<networkFace*> all_boundary_faces;
	  for (auto* f : water->getBoundaryFaces())
	    all_boundary_faces.emplace(f);
	  auto should_use_bie_interface = [&](const BEM_DOF_Base& node, const networkFace*) {
	    if (!(radiation_interface_bie_wall_free || radiation_interface_bie_float_free))
	      return false;
	    bool touches_float = false;
	    bool touches_free = false;
	    bool touches_wall = false;
	    for (auto* adj : node.getBoundaryFaces()) {
	      if (face_sets.float_surface.contains(adj))
	        touches_float = true;
	      else if (face_sets.free_surface.contains(adj))
	        touches_free = true;
	      else
	        touches_wall = true;
	    }
	    if (radiation_interface_bie_wall_free && !touches_float && touches_free && touches_wall)
	      return true;
	    if (radiation_interface_bie_float_free && touches_float && touches_free)
	      return true;
	    return false;
	  };
	  if (radiation_debug) {
	    bem_frequency_domain::BoundaryData debug_bc;
	    debug_bc.face_bc = [&](const networkFace& f) -> bem_frequency_domain::FaceBC {
	      if (face_sets.free_surface.contains(const_cast<networkFace*>(&f)))
	        return radiation_free_surface_dirichlet_zero ? bem_frequency_domain::FaceBC::Dirichlet
	                                                    : bem_frequency_domain::FaceBC::Robin;
	      return bem_frequency_domain::FaceBC::Neumann;
	    };
	    debug_bc.neumann_phin = [&](const BEM_DOF_Base&, const networkFace*) -> bem_frequency_domain::Complex {
	      return bem_frequency_domain::Complex{0.0, 0.0};
	    };
	    debug_bc.dirichlet_phi = [&](const BEM_DOF_Base&) -> bem_frequency_domain::Complex {
	      return bem_frequency_domain::Complex{0.0, 0.0};
	    };
	    debug_bc.robin_unknown_policy = radiation_robin_unknown_policy;
	    debug_bc.use_bie_row_for_interface = should_use_bie_interface;
	    if (radiation_face_neumann) {
	      debug_bc.force_neumann_multiple = [&](const BEM_DOF_Base& node) -> bool {
	        return std::ranges::any_of(node.getBoundaryFaces(), [&](const networkFace* f) {
	          return face_sets.float_surface.contains(const_cast<networkFace*>(f));
	        });
	      };
	    }
	    write_radiation_mode_consistency(outdir / "radiation_mode_consistency.csv",
	                                     face_sets.float_surface,
	                                     radiation_ref_point,
	                                     radiation_rotation_sign,
	                                     radiation_moment_sign);
	    write_radiation_surface_summary(outdir / "radiation_surface_summary.csv",
	                                    *water,
	                                    face_sets,
	                                    radiation_ref_point,
	                                    radiation_moment_sign);
	    write_radiation_time_domain_geometry_compare(outdir / "radiation_time_domain_geometry_compare.csv",
	                                                 face_sets.float_surface,
	                                                 radiation_ref_point,
	                                                 radiation_moment_sign);
	    write_panel_operator_symmetry(outdir / "radiation_panel_operator_symmetry.csv",
	                                  omegas.empty() ? 0.0 : omegas.front(),
	                                  "constant_panel",
	                                  face_sets.float_surface);
	    write_panel_pair_integral_check(outdir / "radiation_panel_pair_integral_check.csv",
	                                    omegas.empty() ? 0.0 : omegas.front(),
	                                    "constant_panel",
	                                    face_sets.float_surface,
	                                    radiation_ref_point,
	                                    radiation_moment_sign);
	    write_panel_double_layer_adjoint(outdir / "radiation_panel_double_layer_adjoint.csv",
	                                     omegas.empty() ? 0.0 : omegas.front(),
	                                     "constant_panel",
	                                     face_sets.float_surface);
	    write_panel_mode_projection(outdir / "radiation_panel_mode_projection.csv",
	                                omegas.empty() ? 0.0 : omegas.front(),
	                                "constant_panel",
	                                face_sets.float_surface,
	                                radiation_ref_point,
	                                radiation_moment_sign);
	    write_panel_modal_solve(outdir / "radiation_panel_modal_solve.csv",
	                            omegas.empty() ? 0.0 : omegas.front(),
	                            "constant_panel",
	                            face_sets.float_surface,
	                            radiation_ref_point,
	                            radiation_moment_sign);
	    write_radiation_bie_dof_diagnostics(outdir / "radiation_bie_dof_summary.csv",
	                                        outdir / "radiation_bie_interface.csv",
	                                        setting.FluidObject,
	                                        face_sets,
	                                        debug_bc);
	  }
  const double water_depth = std::get<1>(water->bounds[2]) - std::get<0>(water->bounds[2]);
  auto absorber_distance = [&](const Tddd& x) -> std::pair<bool, double> {
    bool inside = false;
    double min_distance = std::numeric_limits<double>::infinity();
    for (auto* absorber : setting.AbsorberObject) {
      if (!absorber)
        continue;
      if (!absorber->InsideQ(x))
        continue;
      inside = true;
      auto [near_f, x_nearest] = absorber->Nearest(x);
      if (near_f)
        min_distance = std::min(min_distance, Norm(x - x_nearest));
    }
    if (!inside)
      return {false, 0.0};
    if (!std::isfinite(min_distance))
      min_distance = 0.0;
    return {true, min_distance};
  };
		  auto sponge_mu_at_position = [&](const Tddd& x, double omega) -> double {
	    if (!sponge.enabled)
	      return 0.0;
    double s = 0.0;
    auto wall_mu_at = [&]() -> double {
      const double x_min = std::get<0>(water->bounds[0]);
      const double x_max = std::get<1>(water->bounds[0]);
      const double y_min = std::get<0>(water->bounds[1]);
      const double y_max = std::get<1>(water->bounds[1]);
      const double d_wall = std::min({x[0] - x_min, x_max - x[0], x[1] - y_min, y_max - x[1]});
      const double width = (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax) ? sponge.wall_width : sponge.length;
	      const double mu_max = (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax) ? sponge.wall_mu_max : sponge.mu_max;
	      const int order = (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax) ? sponge.wall_order : sponge.order;
	      const double sw = std::clamp((width - d_wall) / width, 0.0, 1.0);
	      return mu_max * sponge_profile_value_with_order(sw, sponge, order);
    };
    auto absorber_mu_at = [&]() -> double {
      auto [inside, sd] = absorber_distance(x);
      if (!inside)
        return 0.0;
      const double k = finite_depth_wave_number(omega, water_depth, g);
      const double wavelength = 2.0 * M_PI / k;
      const double width = sponge.n_wavelengths * wavelength;
      const double sa = std::clamp(sd / width, 0.0, 1.0);
      return sponge.mu_max * sponge_profile_value(sa, sponge);
    };
    if (sponge.mode == FrequencySpongeRuntime::Mode::WallDistance) {
      return wall_mu_at();
    } else if (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWavelength) {
      return absorber_mu_at();
    } else if (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax) {
      return std::max(absorber_mu_at(), wall_mu_at());
    } else {
      const double dx = x[0] - float_body->COM[0];
      const double dy = x[1] - float_body->COM[1];
      const double r = std::sqrt(dx * dx + dy * dy);
      s = std::clamp((r - sponge.r_start) / sponge.length, 0.0, 1.0);
    }
	    return sponge.mu_max * sponge_profile_value(s, sponge);
	  };
		  auto make_sponge_mu = [&](double omega) -> std::function<double(const BEM_DOF_Base&)> {
		    if (!sponge.enabled)
		      return {};
	    return [&, omega](const BEM_DOF_Base& node) -> double {
	      return sponge_mu_at_position(node.getPosition(), omega);
	    };
	  };
	  auto touches_open_boundary = [&](const BEM_DOF_Base& node) -> bool {
	    for (auto* adj : node.getBoundaryFaces()) {
	      if (adj && face_sets.outer_wall.contains(adj))
	        return true;
	    }
	    return false;
	  };
	  auto touches_free_surface = [&](const BEM_DOF_Base& node) -> bool {
	    for (auto* adj : node.getBoundaryFaces()) {
	      if (adj && face_sets.free_surface.contains(adj))
	        return true;
	    }
	    return false;
	  };
	  auto open_boundary_applies_to = [&](const BEM_DOF_Base& node, const networkFace* f) -> bool {
	    if (!freq_open_boundary_sommerfeld)
	      return false;
	    if (f)
	      return face_sets.outer_wall.contains(const_cast<networkFace*>(f));
	    // For a face-keyless ID on a free-surface/wall corner, keep the
	    // free-surface Robin coefficient.  Face-keyed wall rows still get the
	    // outgoing coefficient.
	    return touches_open_boundary(node) && !touches_free_surface(node);
	  };
	  auto complex_wave_number_from_sponge_mu = [&](double omega, double mu) -> bem_frequency_domain::Complex {
	    using bem_frequency_domain::Complex;
	    if (!(omega > 0.0))
	      return Complex{0.0, 0.0};
	    if (!(mu > 0.0))
	      return Complex{finite_depth_wave_number(omega, water_depth, g), 0.0};
	    const Complex s(mu, -omega);
	    const Complex rhs = -(s * s) / g;
	    Complex k{finite_depth_wave_number(omega, water_depth, g), 0.0};
	    if (std::abs(k) < 1e-14)
	      k = std::sqrt(rhs);
	    for (int iter = 0; iter < 40; ++iter) {
	      const Complex kh = k * water_depth;
	      const Complex th = std::tanh(kh);
	      const Complex sech2 = Complex{1.0, 0.0} - th * th;
	      const Complex f = k * th - rhs;
	      const Complex df = th + k * water_depth * sech2;
	      if (std::abs(df) < 1e-20)
	        break;
	      const Complex dk = f / df;
	      k -= dk;
	      if (std::abs(dk) <= 1e-12 * std::max(1.0, std::abs(k)))
	        break;
	    }
	    // f(k)=k tanh(kh) is even.  Keep the outgoing branch with positive real
	    // part; for the damped sponge branch this also gives positive Im(k).
	    if (k.real() < 0.0)
	      k = -k;
	    return k;
	  };
		  auto open_boundary_wave_number_at_position = [&](const Tddd& x, double omega) -> bem_frequency_domain::Complex {
		    if (freq_open_boundary_wave_number == "real")
		      return bem_frequency_domain::Complex{finite_depth_wave_number(omega, water_depth, g), 0.0};
		    const double mu = freq_open_boundary_mu_scale * sponge_mu_at_position(x, omega);
		    return complex_wave_number_from_sponge_mu(omega, mu);
		  };
		  auto open_boundary_radius_at_position = [&](const Tddd& x) -> double {
		    const double dx = x[0] - float_body->COM[0];
		    const double dy = x[1] - float_body->COM[1];
		    return std::sqrt(dx * dx + dy * dy);
		  };
		  auto open_boundary_normal_at = [&](const BEM_DOF_Base& node, const networkFace* f) -> Tddd {
		    if (f)
		      return f->normal;
		    Tddd normal = {0.0, 0.0, 0.0};
		    for (auto* adj : node.getBoundaryFaces()) {
		      if (adj && face_sets.outer_wall.contains(adj))
		        normal += adj->normal;
		    }
		    const double nm = Norm(normal);
		    if (nm > 1e-15)
		      return normal / nm;
		    return {0.0, 0.0, 0.0};
		  };
		  auto open_boundary_radial_projection_at = [&](const BEM_DOF_Base& node, const networkFace* f) -> double {
		    if (!freq_open_boundary_radial_projection)
		      return 1.0;
		    const auto x = node.getPosition();
		    const double r = open_boundary_radius_at_position(x);
		    if (!(r > 1e-12))
		      return 1.0;
		    const Tddd er = {(x[0] - float_body->COM[0]) / r,
		                     (x[1] - float_body->COM[1]) / r,
		                     0.0};
		    const Tddd n = open_boundary_normal_at(node, f);
		    const double ndotr = Dot(n, er);
		    return std::isfinite(ndotr) ? ndotr : 1.0;
		  };
		  auto open_boundary_cylindrical_decay_at_position = [&](const Tddd& x) -> double {
		    if (freq_open_boundary_curvature != "cylindrical")
		      return 0.0;
		    const double r = open_boundary_radius_at_position(x);
		    return (r > 1e-12) ? 0.5 / r : 0.0;
		  };
		  auto open_boundary_kappa_at = [&](const BEM_DOF_Base& node, const networkFace* f, double omega) -> bem_frequency_domain::Complex {
		    const auto x = node.getPosition();
		    const auto k_open = open_boundary_wave_number_at_position(x, omega);
		    bem_frequency_domain::Complex kappa =
		        bem_frequency_domain::Complex{0.0, 1.0} * k_open -
		        bem_frequency_domain::Complex{open_boundary_cylindrical_decay_at_position(x), 0.0};
		    kappa *= open_boundary_radial_projection_at(node, f);
		    return freq_open_boundary_kappa_scale * kappa;
		  };
		  auto configure_open_boundary_robin =
		      [&](bem_frequency_domain::BoundaryData& bc, const bem_frequency_domain::LinearFSBC& fsbc, double omega) {
		        if (!freq_open_boundary_sommerfeld)
		          return;
		        bc.robin_kappa = [&, omega](const BEM_DOF_Base& node, const networkFace* f) -> bem_frequency_domain::Complex {
		          if (open_boundary_applies_to(node, f)) {
		            return open_boundary_kappa_at(node, f, omega);
		          }
		          return fsbc.kappa_at(node);
		        };
		      };
	  {
    const auto profile_path = outdir / "frequency_sponge_profile.csv";
    std::ofstream spfs(profile_path);
		    spfs << "x,y,z,r,d_wall,absorber_inside,absorber_distance,profile_omega,profile_wavelength,profile_width,profile_s,profile_value,mu\n";
    const double x_min = std::get<0>(water->bounds[0]);
    const double x_max = std::get<1>(water->bounds[0]);
    const double y_min = std::get<0>(water->bounds[1]);
    const double y_max = std::get<1>(water->bounds[1]);
    const double profile_omega = omegas.empty() ? 0.0 : omegas.front();
    double profile_wavelength = 0.0;
    double profile_width = sponge.length;
    if (sponge.enabled && (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWavelength ||
                           sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax) &&
        profile_omega > 0.0) {
      profile_wavelength = 2.0 * M_PI / finite_depth_wave_number(profile_omega, water_depth, g);
      profile_width = sponge.n_wavelengths * profile_wavelength;
    }
    for (auto* f : face_sets.free_surface) {
      const auto c = f->center;
      const double dx = c[0] - float_body->COM[0];
      const double dy = c[1] - float_body->COM[1];
      const double d_wall = std::min({c[0] - x_min, x_max - c[0], c[1] - y_min, y_max - c[1]});
      auto [absorber_inside, absorber_sd] = absorber_distance(c);
      double profile_s = 0.0;
      if (sponge.enabled && sponge.mode == FrequencySpongeRuntime::Mode::WallDistance) {
        profile_s = std::clamp((sponge.length - d_wall) / sponge.length, 0.0, 1.0);
      } else if (sponge.enabled &&
                 (sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWavelength ||
                  sponge.mode == FrequencySpongeRuntime::Mode::AbsorberWallMax) &&
                 absorber_inside) {
        profile_s = profile_width > 0.0 ? std::clamp(absorber_sd / profile_width, 0.0, 1.0) : 0.0;
      } else if (sponge.enabled && sponge.length > 0.0) {
        profile_s = std::clamp((std::sqrt(dx * dx + dy * dy) - sponge.r_start) / sponge.length, 0.0, 1.0);
      }
      const double profile_value = sponge.enabled ? sponge_profile_value(profile_s, sponge) : 0.0;
      spfs << std::scientific << std::setprecision(12)
	           << c[0] << "," << c[1] << "," << c[2] << ","
	           << std::sqrt(dx * dx + dy * dy) << "," << d_wall << ","
		           << (absorber_inside ? 1 : 0) << "," << absorber_sd << ","
		           << profile_omega << "," << profile_wavelength << "," << profile_width << ","
		           << profile_s << "," << profile_value << ","
		           << sponge_mu_at_position(c, profile_omega) << "\n";
	    }
	    std::cout << "wrote: " << profile_path << std::endl;
	  }
	  {
	    const auto open_path = outdir / "frequency_open_boundary_summary.csv";
	    std::ofstream obs(open_path);
	    double outer_wall_area = 0.0;
	    for (auto* f : face_sets.outer_wall)
	      outer_wall_area += f ? f->area : 0.0;
		    obs << "mode,wave_number_mode,omega,k_real,k_imag,scale,mu_scale,"
		        << "wall_mu_min,wall_mu_mean,wall_mu_max,outer_wall_faces,outer_wall_area,corner_policy,"
		        << "curvature_mode,radial_projection,radial_factor_min,radial_factor_mean,radial_factor_max,"
		        << "curvature_decay_min,curvature_decay_mean,curvature_decay_max,kappa_re_mean,kappa_im_mean\n";
		    for (double omega : omegas) {
		      double mu_min = std::numeric_limits<double>::infinity();
		      double mu_max = 0.0;
		      double mu_sum = 0.0;
		      double radial_min = std::numeric_limits<double>::infinity();
		      double radial_max = -std::numeric_limits<double>::infinity();
		      double radial_sum = 0.0;
		      double decay_min = std::numeric_limits<double>::infinity();
		      double decay_max = 0.0;
		      double decay_sum = 0.0;
		      std::size_t mu_count = 0;
		      bem_frequency_domain::Complex k_sum{0.0, 0.0};
		      bem_frequency_domain::Complex kappa_sum{0.0, 0.0};
		      for (auto* f : face_sets.outer_wall) {
		        if (!f)
		          continue;
		        const double mu = freq_open_boundary_mu_scale * sponge_mu_at_position(f->center, omega);
		        const auto k_eff = open_boundary_wave_number_at_position(f->center, omega);
		        const double r = open_boundary_radius_at_position(f->center);
		        double radial_factor = 1.0;
		        if (freq_open_boundary_radial_projection && r > 1e-12) {
		          const Tddd er = {(f->center[0] - float_body->COM[0]) / r,
		                           (f->center[1] - float_body->COM[1]) / r,
		                           0.0};
		          radial_factor = Dot(f->normal, er);
		        }
		        const double decay = open_boundary_cylindrical_decay_at_position(f->center);
		        const bem_frequency_domain::Complex kappa =
		            freq_open_boundary_kappa_scale * radial_factor *
		            (bem_frequency_domain::Complex{0.0, 1.0} * k_eff -
		             bem_frequency_domain::Complex{decay, 0.0});
		        mu_min = std::min(mu_min, mu);
		        mu_max = std::max(mu_max, mu);
		        mu_sum += mu;
		        radial_min = std::min(radial_min, radial_factor);
		        radial_max = std::max(radial_max, radial_factor);
		        radial_sum += radial_factor;
		        decay_min = std::min(decay_min, decay);
		        decay_max = std::max(decay_max, decay);
		        decay_sum += decay;
		        k_sum += k_eff;
		        kappa_sum += kappa;
		        ++mu_count;
		      }
		      if (mu_count == 0)
		        mu_min = 0.0;
		      const double mu_mean = mu_count ? mu_sum / static_cast<double>(mu_count) : 0.0;
		      if (mu_count == 0) {
		        radial_min = 0.0;
		        radial_max = 0.0;
		        decay_min = 0.0;
		      }
		      const double radial_mean = mu_count ? radial_sum / static_cast<double>(mu_count) : 0.0;
		      const double decay_mean = mu_count ? decay_sum / static_cast<double>(mu_count) : 0.0;
		      const auto k_mean = mu_count ? k_sum / static_cast<double>(mu_count)
		                                   : bem_frequency_domain::Complex{0.0, 0.0};
		      const auto kappa_mean = mu_count ? kappa_sum / static_cast<double>(mu_count)
		                                       : bem_frequency_domain::Complex{0.0, 0.0};
		      obs << freq_open_boundary_mode << "," << freq_open_boundary_wave_number << ","
		          << std::scientific << std::setprecision(12)
		          << omega << "," << k_mean.real() << "," << k_mean.imag() << ","
		          << freq_open_boundary_kappa_scale << "," << freq_open_boundary_mu_scale << ","
		          << mu_min << "," << mu_mean << "," << mu_max << ","
		          << face_sets.outer_wall.size() << "," << outer_wall_area
		          << ",face_keyed_wall_robin_free_wall_corner_prefers_free_surface,"
		          << freq_open_boundary_curvature << "," << (freq_open_boundary_radial_projection ? 1 : 0) << ","
		          << radial_min << "," << radial_mean << "," << radial_max << ","
		          << decay_min << "," << decay_mean << "," << decay_max << ","
		          << kappa_mean.real() << "," << kappa_mean.imag() << "\n";
	    }
	    std::cout << "wrote: " << open_path << std::endl;
	  }
	  const bool lu_dunavant_enabled = g_lu_far_dunavant_points > 0 || g_lu_near_dunavant_points > 0;
  const char* radiation_solver_name = radiation_complex_lu ? "complex_lu" : (radiation_complex_gmres ? "complex_gmres_dense" : (radiation_complex_fmm_gmres ? "complex_fmm_gmres" : "real_gmres"));
  const char* wave_solver_name = wave_complex_fmm_gmres ? "complex_fmm_gmres" : (radiation_complex_gmres ? "complex_gmres_dense" : "complex_lu");
  {
    const auto settings_path = outdir / "frequency_solver_settings.csv";
    std::ofstream fs(settings_path);
    const char* element_order = use_true_quadratic_element ? "true_quadratic" : (use_pseudo_quadratic_element ? "pseudo_quadratic" : "linear");
    const char* postprocess_field_order = use_true_quadratic_element ? "true_quadratic_6node" : "linear_3node";
    fs << "key,value\n"
       << "element_order," << element_order << "\n"
       << "postprocess_field_order," << postprocess_field_order << "\n"
		       << "radiation_solver," << radiation_solver_name << "\n"
		       << "wave_solver," << (compute_wave ? wave_solver_name : "none") << "\n"
			       << "wave_formulation,scattered_potential\n"
			       << "incident_wave_added_after_solve," << (compute_wave ? 1 : 0) << "\n"
			       << "sponge_applies_to,scattered_and_radiation_potential\n"
			       << "open_boundary_mode," << freq_open_boundary_mode << "\n"
			       << "open_boundary_kappa_scale," << freq_open_boundary_kappa_scale << "\n"
				       << "open_boundary_wave_number," << freq_open_boundary_wave_number << "\n"
				       << "open_boundary_mu_scale," << freq_open_boundary_mu_scale << "\n"
				       << "open_boundary_curvature," << freq_open_boundary_curvature << "\n"
				       << "open_boundary_radial_projection," << (freq_open_boundary_radial_projection ? 1 : 0) << "\n"
				       << "open_boundary_condition,"
				       << (freq_open_boundary_sommerfeld
				               ? "phin=scale*radial_factor*(i*k-curvature_decay)*phi_on_outer_wall"
				               : "none")
				       << "\n"
			       << "free_surface_animation_frames," << freq_free_surface_animation_frames << "\n"
			       << "free_surface_animation_scale," << freq_free_surface_animation_scale << "\n"
		       << "complex_gmres_tol," << complex_gmres_tol << "\n"
       << "complex_gmres_max_iter," << complex_gmres_max_iter << "\n"
       << "complex_gmres_restart," << complex_gmres_restart << "\n"
       << "complex_gmres_preconditioner," << complex_gmres_preconditioner << "\n"
	       << "complex_fmm_matrix_free," << ((radiation_complex_fmm_gmres || (compute_wave && wave_complex_fmm_gmres)) ? 1 : 0) << "\n"
	       << "fmm_setup_source," << ((radiation_complex_fmm_gmres || (compute_wave && wave_complex_fmm_gmres)) ? "time_domain_shared" : "none") << "\n"
	       << "fmm_true_quad_straight_geometry_probe," << (fmm_true_quad_straight_geometry_probe ? 1 : 0) << "\n"
	       << "radiation_fmm_matvec_compare," << (radiation_fmm_matvec_compare ? 1 : 0) << "\n"
	       << "radiation_fmm_coordinate_scaling," << (radiation_fmm_no_coordinate_scaling ? 0 : 1) << "\n"
	       << "wave_fmm_coordinate_scaling," << (wave_fmm_no_coordinate_scaling ? 0 : 1) << "\n"
	       << "fmm_max_level," << fmm_max_level << "\n"
	       << "fmm_bucket_max_points," << fmm_bucket_max_points << "\n"
	       << "fmm_mac_theta," << g_mac_theta << "\n"
	       << "nearfield_mode," << nearfield_mode << "\n"
	       << "lu_quadrature_mode," << (lu_dunavant_enabled ? "diagnostic_dunavant" : "default_gw") << "\n"
       << "lu_far_dunavant_points," << g_lu_far_dunavant_points << "\n"
       << "lu_near_dunavant_points," << g_lu_near_dunavant_points << "\n"
       << "lu_far_default,gw1xgw1\n"
       << "lu_near_default,gw5xgw5\n"
       << "lu_singular_default,duffy_gw5xgw5\n"
       << "sponge_enabled," << (sponge.enabled ? 1 : 0) << "\n"
       << "sponge_mode," << frequency_sponge_mode_name(sponge.mode) << "\n"
       << "sponge_mu_max," << sponge.mu_max << "\n"
       << "sponge_r_start," << sponge.r_start << "\n"
       << "sponge_length," << sponge.length << "\n"
       << "sponge_n_wavelengths," << sponge.n_wavelengths << "\n"
       << "sponge_order," << sponge.order << "\n"
       << "sponge_profile," << sponge.profile << "\n"
       << "sponge_floor_fraction," << sponge.floor_fraction << "\n"
       << "sponge_wall_mu_max," << sponge.wall_mu_max << "\n"
       << "sponge_wall_width," << sponge.wall_width << "\n"
       << "sponge_wall_order," << sponge.wall_order << "\n"
       << "sponge_water_depth," << water_depth << "\n"
	       << "sponge_absorber_count," << setting.AbsorberObject.size() << "\n"
		       << "radiation_robin_unknown_policy,"
		       << bem_frequency_domain::robin_unknown_policy_name(radiation_robin_unknown_policy) << "\n"
		       << "wave_robin_unknown_policy,"
		       << bem_frequency_domain::robin_unknown_policy_name(radiation_robin_unknown_policy) << "\n"
		       << "radiation_robin_phi_unknown_legacy," << (radiation_robin_phi_unknown_legacy ? 1 : 0) << "\n"
	       << "radiation_interface_row_scale," << radiation_interface_row_scale << "\n"
	       << "radiation_interface_robin_adjoint_probe," << (radiation_interface_robin_adjoint_probe ? 1 : 0) << "\n"
	       << "radiation_postprocess_probe," << (radiation_postprocess_probe ? 1 : 0) << "\n"
	       << "radiation_pitch_probe," << (radiation_pitch_probe ? 1 : 0) << "\n"
	       << "radiation_coupling_probe," << (radiation_coupling_probe ? 1 : 0) << "\n"
	       << "wave_excitation_group_probe," << (wave_excitation_group_probe ? 1 : 0) << "\n";
    std::cout << "wrote: " << settings_path << std::endl;
  }
  std::cout << "frequency sponge: enabled=" << sponge.enabled
            << " mode=" << frequency_sponge_mode_name(sponge.mode)
            << " mu_max=" << sponge.mu_max
            << " r_start=" << sponge.r_start
            << " length=" << sponge.length
            << " n_wavelengths=" << sponge.n_wavelengths
            << " order=" << sponge.order
            << " profile=" << sponge.profile
            << " floor_fraction=" << sponge.floor_fraction
            << " wall_mu_max=" << sponge.wall_mu_max
            << " wall_width=" << sponge.wall_width
            << " wall_order=" << sponge.wall_order
            << " water_depth=" << water_depth
            << " absorber_count=" << setting.AbsorberObject.size()
            << " radiation_solver=" << (radiation_complex_lu ? "complex_lu" : (radiation_complex_gmres ? "complex_gmres_dense" : (radiation_complex_fmm_gmres ? "complex_fmm_gmres" : "real_gmres")))
            << " wave_solver=" << (compute_wave ? wave_solver_name : "none")
            << " complex_gmres_tol=" << complex_gmres_tol
            << " complex_gmres_max_iter=" << complex_gmres_max_iter
            << " complex_gmres_restart=" << complex_gmres_restart
	            << " complex_gmres_preconditioner=" << complex_gmres_preconditioner
	            << " fmm_max_level=" << fmm_max_level
	            << " fmm_bucket_max_points=" << fmm_bucket_max_points
	            << " fmm_mac_theta=" << g_mac_theta
	            << " nearfield_mode=" << nearfield_mode
	            << " lu_quadrature_mode=" << (lu_dunavant_enabled ? "diagnostic_dunavant" : "default_gw")
            << " lu_far_dunavant_points=" << g_lu_far_dunavant_points
            << " lu_near_dunavant_points=" << g_lu_near_dunavant_points
            << " radiation_debug=" << radiation_debug
            << " radiation_postprocess_probe=" << radiation_postprocess_probe
            << " radiation_pitch_probe=" << radiation_pitch_probe
            << " radiation_coupling_probe=" << radiation_coupling_probe
            << " radiation_ref_point=(" << radiation_ref_point[0] << "," << radiation_ref_point[1] << "," << radiation_ref_point[2] << ")"
            << " radiation_rotation_sign=" << radiation_rotation_sign
            << " radiation_moment_sign=" << radiation_moment_sign
            << " radiation_face_neumann=" << radiation_face_neumann
            << " radiation_symmetrized_lu=" << radiation_symmetrized_lu
            << " radiation_free_surface_dirichlet_zero=" << radiation_free_surface_dirichlet_zero
            << " radiation_robin_unknown_policy=" << bem_frequency_domain::robin_unknown_policy_name(radiation_robin_unknown_policy)
            << " radiation_robin_phi_unknown=" << radiation_robin_phi_unknown
            << " radiation_robin_phi_unknown_legacy=" << radiation_robin_phi_unknown_legacy
            << " radiation_robin_kappa_scale=" << radiation_robin_kappa_scale
            << " radiation_interface_bie_wall_free=" << radiation_interface_bie_wall_free
            << " radiation_interface_bie_float_free=" << radiation_interface_bie_float_free
            << " radiation_condense_interface=" << radiation_condense_interface
            << " radiation_interface_row_scale=" << radiation_interface_row_scale
	            << " radiation_interface_robin_adjoint_probe=" << radiation_interface_robin_adjoint_probe
	            << " fmm_true_quad_straight_geometry_probe=" << fmm_true_quad_straight_geometry_probe
	            << " radiation_fmm_matvec_compare=" << radiation_fmm_matvec_compare
			            << " radiation_fmm_coordinate_scaling=" << (radiation_fmm_no_coordinate_scaling ? 0 : 1)
			            << " wave_fmm_coordinate_scaling=" << (wave_fmm_no_coordinate_scaling ? 0 : 1)
				            << " open_boundary_mode=" << freq_open_boundary_mode
				            << " open_boundary_wave_number=" << freq_open_boundary_wave_number
				            << " open_boundary_mu_scale=" << freq_open_boundary_mu_scale
				            << " open_boundary_curvature=" << freq_open_boundary_curvature
				            << " open_boundary_radial_projection=" << (freq_open_boundary_radial_projection ? 1 : 0)
				            << std::endl;
	  auto write_reciprocity_column = [&](double omega, int dof_col, const std::array<Complex, 6>& S_col, const char* solver_name) {
    for (int dof_row = 0; dof_row < 6; ++dof_row) {
      reciprocity_fs << std::scientific << std::setprecision(12)
                     << omega << "," << dof_row << "," << dof_col << ","
                     << S_col[static_cast<std::size_t>(dof_row)].real() << ","
	                     << S_col[static_cast<std::size_t>(dof_row)].imag() << ","
	                     << solver_name << "," << (sponge.enabled ? 1 : 0) << ","
	                     << sponge.mu_max << "," << sponge.r_start << ","
	                     << sponge.length << "," << sponge.order << ","
	                     << radiation_ref_point[0] << "," << radiation_ref_point[1] << "," << radiation_ref_point[2] << ","
	                     << radiation_rotation_sign << "," << radiation_moment_sign << "\n";
    }
  };
  auto append_linear_solver_residual = [&](const char* problem,
                                           double omega,
                                           int dof_col,
                                           const char* solver_name,
                                           const bem_frequency_domain::Solution& sol) {
    const auto path = outdir / "linear_solver_residual.csv";
    const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
    std::ofstream fs(path, std::ios::app);
    if (!exists) {
      fs << "problem,omega,dof_col,solver,effective_preconditioner,iterations,"
            "absolute_residual,relative_residual,converged\n";
    }
    fs << problem << "," << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << solver_name << ","
       << sol.linear_solver_effective_preconditioner << ","
       << sol.linear_solver_iterations << ","
       << sol.linear_solver_residual_l2 << ","
       << sol.linear_solver_relative_residual_l2 << ","
       << (sol.linear_solver_converged ? 1 : 0) << "\n";
  };
  const char* complex_lu_solver_name = radiation_complex_fmm_gmres
                                           ? "complex_fmm_gmres"
                                           : (radiation_complex_gmres
                                                  ? (radiation_condense_interface ? "complex_gmres_dense_condensed_interface" : "complex_gmres_dense")
                                                  : (radiation_condense_interface ? "complex_lu_condensed_interface" : "complex_lu"));
  auto append_fmm_nearfield_summary = [&](double omega, int dof_col, const bem_frequency_domain::Solution& sol) {
    if (sol.fmm_setup_source.empty())
      return;
    const auto path = outdir / "fmm_nearfield_summary.csv";
    const bool exists = std::filesystem::exists(path) && std::filesystem::file_size(path) > 0;
    std::ofstream fs(path, std::ios::app);
    if (!exists) {
      fs << "omega,dof_col,solver,setup_source,coordinate_scaling,coordinate_scale_factor,morton_reindex,"
            "reused_sources,reused_static_fmm,targets,sources,total_near_terms,mean_nnz_per_target,"
            "max_nnz_per_target,vertex_targets,vertex_total_near_terms,vertex_mean_nnz_per_target,vertex_max_nnz_per_target,"
            "midpoint_targets,midpoint_total_near_terms,midpoint_mean_nnz_per_target,midpoint_max_nnz_per_target,"
            "max_source_offset,mean_source_offset,p95_source_offset,near_all_to_all_warning,vertex_all_to_all_warning,midpoint_all_to_all_warning\n";
    }
    const bool all_to_all_warning =
        sol.n > 0 &&
        sol.fmm_mean_near_terms_per_target >= 0.95 * static_cast<double>(sol.n);
    const bool vertex_all_to_all_warning =
        sol.n > 0 &&
        sol.fmm_vertex_targets > 0 &&
        sol.fmm_vertex_mean_near_terms_per_target >= 0.95 * static_cast<double>(sol.n);
    const bool midpoint_all_to_all_warning =
        sol.n > 0 &&
        sol.fmm_midpoint_targets > 0 &&
        sol.fmm_midpoint_mean_near_terms_per_target >= 0.95 * static_cast<double>(sol.n);
    fs << std::scientific << std::setprecision(12)
       << omega << "," << dof_col << "," << complex_lu_solver_name << ","
       << sol.fmm_setup_source << ","
       << (sol.fmm_coordinate_scaling ? 1 : 0) << ","
       << sol.fmm_coordinate_scale_factor << ","
       << (sol.fmm_morton_reindex ? 1 : 0) << ","
       << (sol.fmm_reused_sources ? 1 : 0) << ","
       << (sol.fmm_reused_static ? 1 : 0) << ","
       << sol.fmm_targets << ","
       << sol.fmm_sources << ","
       << sol.fmm_total_near_terms << ","
       << sol.fmm_mean_near_terms_per_target << ","
       << sol.fmm_max_near_terms_per_target << ","
       << sol.fmm_vertex_targets << ","
       << sol.fmm_vertex_total_near_terms << ","
       << sol.fmm_vertex_mean_near_terms_per_target << ","
       << sol.fmm_vertex_max_near_terms_per_target << ","
       << sol.fmm_midpoint_targets << ","
       << sol.fmm_midpoint_total_near_terms << ","
       << sol.fmm_midpoint_mean_near_terms_per_target << ","
       << sol.fmm_midpoint_max_near_terms_per_target << ","
       << sol.fmm_max_source_offset << ","
       << sol.fmm_mean_source_offset << ","
       << sol.fmm_p95_source_offset << ","
       << (all_to_all_warning ? 1 : 0) << ","
       << (vertex_all_to_all_warning ? 1 : 0) << ","
       << (midpoint_all_to_all_warning ? 1 : 0) << "\n";
    if (use_true_quadratic_element) {
      const auto geom_path = outdir / "true_quadratic_geometry_summary.csv";
      const bool geom_exists = std::filesystem::exists(geom_path) && std::filesystem::file_size(geom_path) > 0;
      std::ofstream gfs(geom_path, std::ios::app);
      if (!geom_exists) {
        gfs << "omega,dof_col,solver,straight_geometry_probe,lines,uninitialized_xmid,nonfinite_xmid,"
               "max_endpoint_midpoint_diff,mean_endpoint_midpoint_diff,scaled_mismatch_suspect,"
               "max_source_offset,mean_source_offset,p95_source_offset,top_source_offsets\n";
      }
      std::size_t lines = 0, uninitialized = 0, nonfinite = 0, scaled_mismatch = 0;
      double max_diff = 0.0, sum_diff = 0.0;
      for (auto* l : water->getBoundaryLines()) {
        ++lines;
        auto [pa, pb] = l->getPoints();
        const Tddd linear_mid = 0.5 * (pa->X + pb->X);
        const bool xmid_zero = (l->X_mid[0] == 0.0 && l->X_mid[1] == 0.0 && l->X_mid[2] == 0.0);
        if (xmid_zero && (linear_mid[0] != 0.0 || linear_mid[1] != 0.0 || linear_mid[2] != 0.0))
          ++uninitialized;
        if (!std::isfinite(l->X_mid[0]) || !std::isfinite(l->X_mid[1]) || !std::isfinite(l->X_mid[2]))
          ++nonfinite;
        const double diff = Norm(l->X_mid - linear_mid);
        max_diff = std::max(max_diff, diff);
        sum_diff += diff;
        const double xmid_norm = Norm(l->X_mid);
        const double linear_norm = Norm(linear_mid);
        if (linear_norm > 1e-12 && (xmid_norm / linear_norm < 0.01 || xmid_norm / linear_norm > 100.0))
          ++scaled_mismatch;
      }
      std::ostringstream top_offsets;
      for (std::size_t i = 0; i < sol.fmm_top_source_offsets.size(); ++i) {
        if (i)
          top_offsets << ";";
        top_offsets << sol.fmm_top_source_offsets[i].first << ":" << sol.fmm_top_source_offsets[i].second;
      }
      gfs << std::scientific << std::setprecision(12)
          << omega << "," << dof_col << "," << complex_lu_solver_name << ","
          << (fmm_true_quad_straight_geometry_probe ? 1 : 0) << ","
          << lines << "," << uninitialized << "," << nonfinite << ","
          << max_diff << "," << (lines ? sum_diff / static_cast<double>(lines) : 0.0) << ","
          << scaled_mismatch << ","
          << sol.fmm_max_source_offset << ","
          << sol.fmm_mean_source_offset << ","
          << sol.fmm_p95_source_offset << ","
          << top_offsets.str() << "\n";
    }
  };

  // Set all boundary node-face states as Neumann (unknown: phi).  Radiation
  // uses the real GMRES path and injects the free-surface Robin relation in the
  // matvec postprocess below.
  for (auto* f : water->getBoundaryFaces()) {
    if (use_true_quadratic_element) {
      f->isTrueQuadraticElement = true;
      f->isPseudoQuadraticElement = false;
      f->isLinearElement = false;
    } else if (use_pseudo_quadratic_element) {
      f->isTrueQuadraticElement = false;
      f->isPseudoQuadraticElement = true;
      f->isLinearElement = false;
    } else {
      f->isTrueQuadraticElement = false;
      f->isPseudoQuadraticElement = false;
      f->isLinearElement = true;
    }
  }

	  bem_frequency_domain::BoundaryData real_bc;
	  real_bc.face_bc = [&](const networkFace& f) -> bem_frequency_domain::FaceBC {
	    if (face_sets.free_surface.contains(const_cast<networkFace*>(&f)))
	      return radiation_free_surface_dirichlet_zero ? bem_frequency_domain::FaceBC::Dirichlet
	                                                  : bem_frequency_domain::FaceBC::Robin;
	    return bem_frequency_domain::FaceBC::Neumann;
	  };
	  real_bc.robin_unknown_policy = radiation_robin_unknown_policy;
	  if (radiation_face_neumann) {
	    real_bc.force_neumann_multiple = [&](const BEM_DOF_Base& node) -> bool {
	      return std::ranges::any_of(node.getBoundaryFaces(), [&](const networkFace* f) {
	        return face_sets.float_surface.contains(const_cast<networkFace*>(f));
	      });
	    };
	  }
	  std::unordered_set<networkFace*> real_robin_faces;
	  bem_frequency_domain::apply_face_bc_to_mesh(setting.FluidObject, real_bc, real_robin_faces);

	  auto force_radiation_neumann_multiple = [&](const BEM_DOF_Base* node) -> bool {
	    return node && real_bc.force_neumann_multiple && real_bc.force_neumann_multiple(*node);
	  };
	  const std::size_t n = setNodeFaceIndices(setting.FluidObject, force_radiation_neumann_multiple);

  // Build per-index Robin mask. solveGMRES may reindex unknowns internally, so
  // the mask is rebuilt through BEM_BVP::after_reindex_unknowns before solve.
  std::vector<unsigned char> is_robin(static_cast<std::size_t>(n), 0);
  const double z_free = std::get<1>(water->bounds[2]);
  auto rebuild_robin_mask = [&]() {
    std::fill(is_robin.begin(), is_robin.end(), 0);
    auto mark_robin_dofs = [&](const auto* node) {
    for (const auto& [f, d] : node->dofs) {
      const int i = d.index;
      if (i < 0 || static_cast<std::size_t>(i) >= n)
        continue;
	      bool robin = false;
	      if (f) {
	        robin = real_robin_faces.contains(f);
	      } else {
	        robin = std::ranges::any_of(node->getBoundaryFaces(), [&](const networkFace* adj) {
	          return real_robin_faces.contains(const_cast<networkFace*>(adj));
	        });
      }
      is_robin[static_cast<std::size_t>(i)] = robin ? 1 : 0;
    }
    };
    for (auto* p : water->getBoundaryPoints()) {
      mark_robin_dofs(p);
    }
    for (auto* l : water->getBoundaryLines())
      mark_robin_dofs(l);
    // Robin mask for midpoint DOFs (true_quadratic).
    if (use_true_quadratic_element) {
      for (auto* l : water->getBoundaryLines()) {
        int i = l->midpoint_index;
        if (i < 0 || static_cast<std::size_t>(i) >= n)
          continue;
        // Edge is on free surface if all adjacent faces are free surface.
        bool robin = true;
        for (auto* f : l->Faces) {
          if (!face_sets.free_surface.contains(f)) {
            robin = false;
            break;
          }
        }
        is_robin[static_cast<std::size_t>(i)] = robin ? 1 : 0;
      }
    }
  };
  rebuild_robin_mask();

  {
    int n_robin = 0, n_total = 0;
    for (std::size_t i = 0; i < n; ++i) {
      ++n_total;
      if (is_robin[i])
        ++n_robin;
    }
    std::cout << "unknowns n=" << n << " (robin=" << n_robin << " non-robin=" << (n_total - n_robin)
              << "), omegas=" << omegas.size() << ", dofs=" << dofs.size() << std::endl;
  }

  for (double omega : omegas) {
    if (!(omega > 0.0))
      continue;
    const double kappa = (g == 0.0) ? 0.0 : (omega * omega / g);
    auto& added_mass = added_mass_by_omega[omega];
	    auto& damping = damping_by_omega[omega];
	    auto& added_mass_raw = added_mass_raw_by_omega[omega];
	    auto& damping_raw = damping_raw_by_omega[omega];
	    added_mass = make_zero_matrix6();
	    damping = make_zero_matrix6();
	    added_mass_raw = make_zero_matrix6();
	    damping_raw = make_zero_matrix6();
	    Matrix6Mask has_radiation_coeff = make_false_matrix6_mask();
		    std::unordered_map<int, bem_frequency_domain::LinearSolution> debug_radiation_all_solutions;
		    std::unordered_map<int, std::vector<Complex>> debug_linear_rhs_by_dof;
		    std::unordered_map<int, std::vector<Complex>> debug_linear_u_by_dof;
		    std::unordered_map<int, std::vector<Complex>> debug_condensed_rhs_by_dof;
		    std::unordered_map<int, std::vector<Complex>> debug_condensed_u_by_dof;
		    bool wrote_bie_operator_symmetry = false;
		    bool wrote_linear_system_blocks = false;

	    std::cout << "\n=== omega=" << omega << " kappa=" << kappa << " ===" << std::endl;

    // Compute selected columns of the radiation impedance.
    for (int dof_col : dofs) {
      if (radiation_complex_solver) {
        bem_frequency_domain::LinearFSBC fsbc;
        fsbc.omega = omega;
        fsbc.gravity = g;
        fsbc.kappa_scale = radiation_robin_kappa_scale;
        fsbc.sponge_mu = make_sponge_mu(omega);

        bem_frequency_domain::BoundaryData bc;
	        bc.face_bc = [&](const networkFace& f) -> bem_frequency_domain::FaceBC {
	          if (face_sets.free_surface.contains(const_cast<networkFace*>(&f)))
	            return radiation_free_surface_dirichlet_zero ? bem_frequency_domain::FaceBC::Dirichlet
	                                                        : bem_frequency_domain::FaceBC::Robin;
	          if (freq_open_boundary_sommerfeld && face_sets.outer_wall.contains(const_cast<networkFace*>(&f)))
	            return bem_frequency_domain::FaceBC::Robin;
	          return bem_frequency_domain::FaceBC::Neumann;
	        };
	        configure_open_boundary_robin(bc, fsbc, omega);
	        bc.neumann_phin = [&](const BEM_DOF_Base& node, const networkFace* f) -> bem_frequency_domain::Complex {
          bool on_float = false;
          Tddd normal = {0.0, 0.0, 0.0};
          if (f) {
            on_float = face_sets.float_surface.contains(const_cast<networkFace*>(f));
            normal = f->normal;
          } else {
            for (auto* adj : node.getBoundaryFaces()) {
              if (face_sets.float_surface.contains(adj)) {
                on_float = true;
                normal += adj->normal;
              }
            }
            const double nm = Norm(normal);
            if (nm > 1e-15)
              normal /= nm;
          }
          if (!on_float)
            return bem_frequency_domain::Complex{0.0, 0.0};
	          const Tddd v = velocity_unit_dof(dof_col, node.getPosition(), radiation_ref_point, radiation_rotation_sign);
	          return bem_frequency_domain::Complex{Dot(v, normal), 0.0};
        };
        bc.dirichlet_phi = [&](const BEM_DOF_Base&) -> bem_frequency_domain::Complex {
          return bem_frequency_domain::Complex{0.0, 0.0};
        };
        bc.robin_unknown_policy = radiation_robin_unknown_policy;
        bc.use_bie_row_for_interface = should_use_bie_interface;
        bc.eliminate_interface_constraints = radiation_condense_interface;
        bc.interface_row_scale_mode = radiation_interface_row_scale;
        bc.interface_robin_adjoint_probe = radiation_interface_robin_adjoint_probe;
        bc.compute_adjoint_residual = radiation_debug;
        bc.linear_solver = radiation_complex_fmm_gmres ? "fmm_gmres" : (radiation_complex_gmres ? "gmres" : "lu");
	        bc.gmres_tol = complex_gmres_tol;
	        bc.gmres_max_iter = complex_gmres_max_iter;
	        bc.gmres_restart = complex_gmres_restart;
	        bc.gmres_preconditioner = complex_gmres_preconditioner;
	        bc.fmm_coordinate_scaling = !radiation_fmm_no_coordinate_scaling;
	        if (radiation_fmm_matvec_compare) {
	          const int debug_dof_col = dof_col;
	          bc.debug_fmm_matvec_compare =
	              [&, debug_dof_col](const std::vector<Complex>& dense_A,
	                                  const std::vector<Complex>& dense_rhs,
	                                  const std::function<std::vector<Complex>(const std::vector<Complex>&)>& fmm_matvec,
	                                  const std::vector<Complex>& fmm_rhs,
	                                  const std::vector<Complex>& u,
	                                  const std::vector<bem_frequency_domain::Id>& id_by_index,
	                                  const std::unordered_set<networkFace*>& robin_faces,
	                                  const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin) {
	                append_radiation_fmm_matvec_compare(outdir,
	                                                    omega,
	                                                    debug_dof_col,
	                                                    complex_lu_solver_name,
	                                                    dense_A,
	                                                    dense_rhs,
	                                                    fmm_matvec,
	                                                    fmm_rhs,
	                                                    u,
	                                                    id_by_index,
	                                                    robin_faces,
	                                                    node_is_robin,
	                                                    face_sets);
	              };
	        }
	        if (radiation_face_neumann) {
	          bc.force_neumann_multiple = [&](const BEM_DOF_Base& node) -> bool {
            return std::ranges::any_of(node.getBoundaryFaces(), [&](const networkFace* f) {
              return face_sets.float_surface.contains(const_cast<networkFace*>(f));
            });
          };
        }
        if (radiation_debug) {
          bc.debug_bie_operator = [&](const BEM_BVP& bvp, const std::vector<bem_frequency_domain::Id>& id_by_index) {
            if (wrote_bie_operator_symmetry)
              return;
            append_bie_operator_symmetry(outdir / "radiation_bie_operator_symmetry.csv",
                                         omega,
                                         complex_lu_solver_name,
                                         bvp,
                                         id_by_index,
                                         all_boundary_faces);
            append_bie_row_multiplicity(outdir / "radiation_bie_row_multiplicity.csv",
                                        omega,
                                        complex_lu_solver_name,
                                        bvp,
                                        id_by_index);
            append_bie_operator_bins(outdir / "radiation_bie_operator_bins.csv",
                                     omega,
                                     complex_lu_solver_name,
                                     bvp,
                                     id_by_index,
                                     all_boundary_faces,
                                     face_sets);
            wrote_bie_operator_symmetry = true;
          };
          const int debug_dof_col = dof_col;
          bc.debug_linear_system = [&, debug_dof_col](const std::vector<Complex>& A,
                                                      const std::vector<Complex>& b,
                                                      const std::vector<bem_frequency_domain::Id>& id_by_index,
                                                      const std::unordered_set<networkFace*>& robin_faces,
                                                      const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin) {
            debug_linear_rhs_by_dof[debug_dof_col] = b;
	            if (wrote_linear_system_blocks)
	              return;
            append_radiation_linear_system_blocks(outdir / "radiation_linear_system_blocks.csv",
                                                  omega,
                                                  debug_dof_col,
                                                  complex_lu_solver_name,
                                                  A,
                                                  b,
                                                  id_by_index,
                                                  robin_faces,
                                                  node_is_robin,
                                                  face_sets,
                                                  radiation_robin_phi_unknown,
                                                  radiation_robin_kappa_scale);
            append_radiation_linear_system_id_map(outdir / "radiation_linear_system_id_map.csv",
                                                  omega,
                                                  debug_dof_col,
                                                  complex_lu_solver_name,
                                                  id_by_index,
                                                  robin_faces,
                                                  node_is_robin,
                                                  face_sets,
                                                  radiation_robin_phi_unknown,
                                                  radiation_robin_kappa_scale);
            append_radiation_interface_constraints(outdir / "radiation_interface_constraints.csv",
                                                   omega,
                                                   debug_dof_col,
                                                   complex_lu_solver_name,
                                                   A,
                                                   b,
                                                   id_by_index,
                                                   robin_faces,
                                                   node_is_robin,
                                                   face_sets,
                                                   g,
                                                   radiation_robin_phi_unknown,
                                                   radiation_robin_kappa_scale);
            wrote_linear_system_blocks = true;
          };
          bc.debug_solved_linear_system = [&, debug_dof_col](const std::vector<Complex>& A,
                                                             const std::vector<Complex>& rhs,
                                                             const std::vector<Complex>& u,
                                                             const std::vector<bem_frequency_domain::Id>& id_by_index,
                                                             const std::unordered_set<networkFace*>& robin_faces,
                                                             const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin) {
            append_radiation_adjoint_column_detail(outdir / "radiation_adjoint_column_detail.csv",
                                                   omega,
                                                   debug_dof_col,
                                                   complex_lu_solver_name,
                                                   A,
                                                   rhs,
                                                   u,
                                                   id_by_index,
                                                   robin_faces,
                                                   node_is_robin,
                                                   face_sets);
            append_radiation_interface_adjoint_pairs(outdir / "radiation_interface_adjoint_pairs.csv",
                                                     omega,
                                                     debug_dof_col,
                                                     complex_lu_solver_name,
                                                     A,
                                                     rhs,
                                                     u,
                                                     id_by_index,
                                                     robin_faces,
                                                     node_is_robin,
                                                     face_sets,
                                                     g,
                                                     radiation_robin_phi_unknown,
                                                     radiation_robin_kappa_scale);
            append_radiation_interface_constraint_balance(outdir / "radiation_interface_constraint_balance.csv",
                                                          omega,
                                                          debug_dof_col,
                                                          complex_lu_solver_name,
                                                          A,
                                                          rhs,
                                                          u,
                                                          id_by_index,
                                                          robin_faces,
                                                          node_is_robin,
                                                          face_sets,
                                                          radiation_robin_phi_unknown,
                                                          radiation_interface_row_scale,
                                                          radiation_interface_robin_adjoint_probe);
            append_radiation_robin_interface_adjoint_pairs(outdir / "radiation_robin_interface_adjoint_pairs.csv",
                                                           omega,
                                                           debug_dof_col,
                                                           complex_lu_solver_name,
                                                           A,
                                                           rhs,
                                                           u,
                                                           id_by_index,
                                                           robin_faces,
                                                           node_is_robin,
                                                           face_sets,
                                                           radiation_robin_phi_unknown,
	                                                           radiation_interface_row_scale,
	                                                           radiation_interface_robin_adjoint_probe);
	          };
	          bc.debug_condensed_solved_linear_system =
	              [&, debug_dof_col](const std::vector<Complex>& A_condensed,
	                                  const std::vector<Complex>& rhs_condensed,
	                                  const std::vector<Complex>& u_condensed,
	                                  const std::vector<bem_frequency_domain::Id>& original_ids,
	                                  const std::vector<bem_frequency_domain::Id>& retained_ids,
	                                  const std::vector<int>& retained_index,
	                                  const std::vector<bem_frequency_domain::CondensedInterfaceExpr>& eliminated,
	                                  const std::unordered_set<networkFace*>& robin_faces,
	                                  const std::unordered_map<BEM_DOF_Base*, bool>& node_is_robin) {
	            debug_condensed_rhs_by_dof[debug_dof_col] = rhs_condensed;
	            debug_condensed_u_by_dof[debug_dof_col] = u_condensed;
	            std::size_t eliminated_count = 0;
	            for (const auto& expr : eliminated)
	              if (expr.eliminated)
	                ++eliminated_count;
	            append_radiation_condensed_system_reciprocity(outdir / "radiation_condensed_system_reciprocity.csv",
	                                                          omega,
	                                                          debug_dof_col,
	                                                          complex_lu_solver_name,
	                                                          A_condensed,
	                                                          rhs_condensed,
	                                                          u_condensed,
	                                                          retained_ids,
	                                                          face_sets,
	                                                          robin_faces,
	                                                          node_is_robin,
	                                                          radiation_robin_phi_unknown,
	                                                          radiation_interface_row_scale,
	                                                          radiation_interface_robin_adjoint_probe,
	                                                          original_ids.size(),
	                                                          eliminated_count);
	            append_radiation_condensed_interface_map(outdir / "radiation_condensed_interface_map.csv",
	                                                     omega,
	                                                     debug_dof_col,
	                                                     complex_lu_solver_name,
	                                                     original_ids,
	                                                     retained_ids,
	                                                     retained_index,
	                                                     eliminated,
	                                                     robin_faces,
	                                                     node_is_robin,
	                                                     face_sets,
	                                                     radiation_robin_phi_unknown,
	                                                     radiation_interface_row_scale,
	                                                     radiation_interface_robin_adjoint_probe);
	          };
	        }

        const auto rad_raw = bem_frequency_domain::solve_linear_bvp({water}, fsbc, bc);
        if (radiation_complex_fmm_gmres)
          append_fmm_nearfield_summary(omega, dof_col, rad_raw);
        append_linear_solver_residual("radiation", omega, dof_col, complex_lu_solver_name, rad_raw);
        if (radiation_debug)
          debug_linear_u_by_dof[dof_col] = rad_raw.u;
        const auto rad_sol = bem_frequency_domain::capture_linear_solution(omega, rad_raw, face_sets.float_surface);
        if (freq_free_surface_animation_frames > 0) {
          const auto free_sol = bem_frequency_domain::capture_linear_solution(omega, rad_raw, face_sets.free_surface);
          radiation_free_surface_solutions[{omega, dof_col}] = free_sol;
          const std::string label = "radiation_dof" + std::to_string(dof_col) + "_omega_" + safe_number_label(omega);
          write_free_surface_eta_animation_components(
              outdir / "free_surface_eta_animation",
              label,
              omega,
              free_sol,
              {{"radiation_dof" + std::to_string(dof_col), &free_sol}},
              [&](const Tddd& x) { return sponge_mu_at_position(x, omega); },
              g,
              freq_free_surface_animation_frames,
              freq_free_surface_animation_scale);
        }
        if (radiation_debug) {
          append_radiation_adjoint_residual_summary(outdir / "radiation_adjoint_residual_summary.csv",
                                                    omega,
                                                    dof_col,
                                                    complex_lu_solver_name,
                                                    rad_raw,
                                                    face_sets);
          append_radiation_bie_residual_summary(outdir / "radiation_bie_residual_summary.csv",
                                                omega,
                                                dof_col,
                                                complex_lu_solver_name,
                                                rad_raw,
                                                face_sets,
                                                should_use_bie_interface);
          append_phi_jump_diagnostic(outdir / "radiation_phi_jump.csv",
                                     omega,
                                     dof_col,
                                     complex_lu_solver_name,
                                     rad_raw,
                                     face_sets);
          const auto rad_sol_all = bem_frequency_domain::capture_linear_solution(omega, rad_raw, all_boundary_faces);
          debug_radiation_all_solutions.emplace(dof_col, rad_sol_all);
	          append_radiation_bc_residual(bc_residual_fs,
	                                       omega,
                                       dof_col,
                                       complex_lu_solver_name,
                                       rad_sol_all,
                                       face_sets,
                                       radiation_ref_point,
                                       radiation_rotation_sign,
                                       [&](const BEM_DOF_Base& node) { return fsbc.kappa_at(node); });
        }
        if (qtf_radiation)
          qtf_solutions[dof_col].push_back(rad_sol);
        if (radiation_postprocess_probe) {
          append_radiation_postprocess_probe(outdir,
                                             omega,
                                             dof_col,
                                             complex_lu_solver_name,
                                             rad_sol,
                                             face_sets.float_surface,
                                             float_body->COM,
                                             rho,
                                             radiation_moment_sign);
        }
	        if (radiation_pitch_probe) {
	          append_radiation_pitch_probe(outdir,
	                                       omega,
	                                       dof_col,
	                                       complex_lu_solver_name,
                                       rad_raw,
                                       rad_sol,
                                       face_sets.float_surface,
                                       float_body->COM,
                                       radiation_ref_point,
                                       rho,
	                                       radiation_rotation_sign,
	                                       radiation_moment_sign);
	        }
	        if (radiation_coupling_probe) {
	          append_radiation_coupling_probe(outdir,
	                                          omega,
	                                          dof_col,
	                                          complex_lu_solver_name,
	                                          rad_sol,
	                                          face_sets.float_surface,
	                                          radiation_ref_point,
	                                          float_body->COM,
	                                          rho,
	                                          radiation_moment_sign);
	        }
		        const auto S_col = integrate_reciprocity_from_solution(rad_sol, face_sets.float_surface, radiation_ref_point, radiation_moment_sign);
	        write_reciprocity_column(omega, dof_col, S_col, complex_lu_solver_name);
	        if (radiation_debug) {
	          const auto S_lumped_col = integrate_reciprocity_from_solution_lumped(rad_sol, face_sets.float_surface, radiation_ref_point, radiation_moment_sign);
	          write_reciprocity_column(omega, dof_col, S_lumped_col, "complex_lu_lumped");
	        }
	        if (radiation_symmetrized_lu) {
	          auto bc_sym = bc;
	          bc_sym.debug_bie_operator = {};
	          bc_sym.transform_bie_operator = [&](BEM_BVP& bvp, const std::vector<bem_frequency_domain::Id>& id_by_index) {
	            symmetrize_area_weighted_bie_operator(bvp, id_by_index, all_boundary_faces);
	          };
	          const auto sym_raw = bem_frequency_domain::solve_linear_bvp({water}, fsbc, bc_sym);
	          const auto sym_sol = bem_frequency_domain::capture_linear_solution(omega, sym_raw, face_sets.float_surface);
	          const auto S_sym_col = integrate_reciprocity_from_solution(sym_sol, face_sets.float_surface, radiation_ref_point, radiation_moment_sign);
	          write_reciprocity_column(omega, dof_col, S_sym_col, "complex_lu_symop");
	          if (radiation_debug) {
	            const auto sym_sol_all = bem_frequency_domain::capture_linear_solution(omega, sym_raw, all_boundary_faces);
	            append_radiation_bc_residual(bc_residual_fs,
	                                         omega,
	                                         dof_col,
	                                         "complex_lu_symop",
	                                         sym_sol_all,
	                                         face_sets,
	                                         radiation_ref_point,
	                                         radiation_rotation_sign,
	                                         [&](const BEM_DOF_Base& node) { return fsbc.kappa_at(node); });
	            const auto sym_residual_path = outdir / "radiation_symmetrized_lu_residual.csv";
	            const bool write_header = !std::filesystem::exists(sym_residual_path) || std::filesystem::file_size(sym_residual_path) == 0;
	            std::ofstream sfs(sym_residual_path, std::ios::app);
	            if (write_header)
	              sfs << "omega,dof,solver,bie_residual_l2\n";
	            sfs << std::scientific << std::setprecision(12)
	                << omega << "," << dof_col << ",complex_lu_symop," << sym_raw.bie_residual_l2 << "\n";
	          }
	          std::cout << "dof " << dof_col << ": complex_lu_symop residual=" << sym_raw.bie_residual_l2 << std::endl;
	        }
	        const auto Z_col = bem_frequency_domain::integrate_linear_pressure_force(rad_sol, face_sets.float_surface, radiation_ref_point, rho, radiation_moment_sign);
        for (int dof_row = 0; dof_row < 6; ++dof_row) {
          const double A = Z_col[dof_row].imag() / omega;
          const double B = -Z_col[dof_row].real();
          added_mass_raw[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] = A;
          damping_raw[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] = B;
          has_radiation_coeff[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] = true;
          raw_coeffs_fs << std::scientific << std::setprecision(12)
                        << omega << "," << dof_row << "," << dof_col << "," << A << "," << B << ","
                        << complex_lu_solver_name << "," << (sponge.enabled ? 1 : 0) << "," << sponge.mu_max << ","
                        << sponge.r_start << "," << sponge.length << "," << sponge.order << "\n";
          zfs << std::scientific << std::setprecision(12)
              << omega << "," << dof_row << "," << dof_col << ","
              << Z_col[dof_row].real() << "," << Z_col[dof_row].imag() << ","
              << A << "," << B << "," << complex_lu_solver_name << ","
              << (sponge.enabled ? 1 : 0) << "," << sponge.mu_max << ","
              << sponge.r_start << "," << sponge.length << "," << sponge.order << "\n";
        }
        std::cout << "dof " << dof_col << ": " << complex_lu_solver_name
                  << " bie_residual=" << rad_raw.bie_residual_l2
                  << " linear_residual=" << rad_raw.linear_solver_residual_l2
                  << " iter=" << rad_raw.linear_solver_iterations
                  << " converged=" << (rad_raw.linear_solver_converged ? 1 : 0) << std::endl;
        continue;
      }

      // Set boundary values.
      auto clear_node_values = [](auto* node) {
        std::get<0>(node->phiphin) = 0.0;
        std::get<1>(node->phiphin) = 0.0;
        std::get<0>(node->phiphin_t) = 0.0;
        std::get<1>(node->phiphin_t) = 0.0;
        for (auto& [f, d] : node->dofs) {
          (void)f;
          d.phi = 0.0;
          d.phin = 0.0;
          d.phi_t = 0.0;
          d.phin_t = 0.0;
        }
      };
      for (auto* p : water->getBoundaryPoints())
        clear_node_values(p);
      for (auto* l : water->getBoundaryLines())
        clear_node_values(l);

		      const Tddd com = radiation_ref_point;
	      std::unordered_map<bem_frequency_domain::Id, double, bem_frequency_domain::IdHash, bem_frequency_domain::IdEq> prescribed_neumann_phin;
	      prescribed_neumann_phin.reserve(n);
	      auto remember_prescribed_phin = [&](BEM_DOF_Base* node, networkFace* f, double phin) {
	        prescribed_neumann_phin[std::make_tuple(node, f)] = phin;
	      };
	      for (auto* p : water->getBoundaryPoints()) {
	        for (auto& [f, d] : p->dofs) {
	          if (d.index < 0)
	            continue;
          // Unknown phi (initial guess).
          d.phi = 0.0;
          d.phi_t = 0.0;

          // Known phin.
          double phin = 0.0;
          bool on_float = false;
          if (f) {
            on_float = face_sets.float_surface.contains(f);
          } else {
            // smooth patch point: treat as float if it is below free surface and within float bbox
            const auto bbox = compute_bbox(*float_body);
            on_float = (p->X[2] < z_free - 1e-8) && (bbox.min[0] - 1e-6 <= p->X[0] && p->X[0] <= bbox.max[0] + 1e-6) && (bbox.min[1] - 1e-6 <= p->X[1] && p->X[1] <= bbox.max[1] + 1e-6) && (bbox.min[2] - 1e-6 <= p->X[2] && p->X[2] <= bbox.max[2] + 1e-6);
          }

          if (on_float) {
	            const Tddd v = velocity_unit_dof(dof_col, p->X, com, radiation_rotation_sign);
            const Tddd nrm = f ? f->normal : p->getNormalNeumann_BEM();
            phin = Dot(v, nrm);
          } else {
            // tank walls/bottom: impermeable (phin=0), free surface Robin handled in postprocess.
            phin = 0.0;
          }

	          d.phin = phin;
	          remember_prescribed_phin(p, f, phin);
	          d.phin_t = 0.0;
	          if (f == nullptr) {
	            std::get<0>(p->phiphin) = d.phi;
            std::get<1>(p->phiphin) = d.phin;
            std::get<0>(p->phiphin_t) = d.phi_t;
            std::get<1>(p->phiphin_t) = d.phin_t;
          }
        }
      }

      // Initialize midpoint boundary values for true_quadratic.
      if (use_true_quadratic_element) {
        for (auto* l : water->getBoundaryLines()) {
          if (l->midpoint_index < 0)
            continue;
          auto [pA, pB] = l->getPoints();
          const Tddd midX = 0.5 * (pA->X + pB->X);
          for (auto& [f, d] : l->dofs) {
            if (d.index < 0)
              continue;
            bool on_float = false;
            Tddd nrm = {0.0, 0.0, 0.0};
            if (f) {
              on_float = face_sets.float_surface.contains(f);
              nrm = f->normal;
            } else {
              on_float = true;
              for (auto* adj : l->Faces) {
                if (!face_sets.float_surface.contains(adj)) {
                  on_float = false;
                  break;
                }
                nrm += adj->normal;
              }
              const double nm = Norm(nrm);
              if (nm > 1e-15)
                nrm /= nm;
            }
	            d.phi = 0.0;
	            d.phi_t = 0.0;
		            d.phin = on_float ? Dot(velocity_unit_dof(dof_col, midX, com, radiation_rotation_sign), nrm) : 0.0;
	            remember_prescribed_phin(l, f, d.phin);
	            d.phin_t = 0.0;
	            if (f == nullptr) {
	              std::get<0>(l->phiphin) = d.phi;
              std::get<1>(l->phiphin) = d.phin;
              std::get<0>(l->phiphin_t) = d.phi_t;
              std::get<1>(l->phiphin_t) = d.phin_t;
            }
          }
        }
      }

      BEM_BVP bvp(setting.FluidObject);
      bvp.matrix_size = static_cast<int>(n);
      bvp.after_reindex_unknowns = rebuild_robin_mask;

      double time_setup = 0.0, time_solve = 0.0;
      TimeWatch watch;

      auto postprocess = [&]() {
        const double kappa_solve =
            (bvp.use_coordinate_scaling_ && bvp.coordinate_scale_factor_ > 1e-10)
                ? bvp.coordinate_scale_factor_ * kappa
                : kappa;
        for (auto& [isDirichlet, i, phi, phin] : bvp.cache_DorN_phi_phin) {
          if (i < 0 || static_cast<std::size_t>(i) >= is_robin.size())
            continue;
          if (!is_robin[static_cast<std::size_t>(i)])
            continue;
          if (isDirichlet) {
            phi = phin / kappa_solve;
            bvp.cache_phi_val_D_by_index[i] = phi;
          } else {
            phin = kappa_solve * phi;
            bvp.cache_phin_val_D_by_index[i] = phin;
          }
        }
      };

      bvp.solveGMRES(watch, time_setup, time_solve, postprocess);

	      // storePhiPhin writes the solved unknowns back to the mesh, but the
	      // mesh-side known Neumann values can be stale afterward. Restore them
	      // before diagnostics/integration read phi/phin from the mesh.
	      auto restore_post_solve_phin = [&](auto* node) {
	        for (auto& [f, d] : node->dofs) {
	          const int i = d.index;
	          if (i < 0 || static_cast<std::size_t>(i) >= is_robin.size())
	            continue;
	          if (is_robin[static_cast<std::size_t>(i)]) {
	            if (isDirichletBieDofKey(node, f))
	              d.phi = d.phin / kappa;
	            else
	              d.phin = kappa * d.phi;
	          } else if (isNeumannBieDofKey(node, f)) {
	            const auto it = prescribed_neumann_phin.find(std::make_tuple(static_cast<BEM_DOF_Base*>(node), f));
	            if (it != prescribed_neumann_phin.end())
	              d.phin = it->second;
	          }
	        }
	        if (const auto* d0 = node->findActiveBieDof(nullptr)) {
	          std::get<0>(node->phiphin) = d0->phi;
	          std::get<1>(node->phiphin) = d0->phin;
	        }
	      };
	      for (auto* p : water->getBoundaryPoints())
	        restore_post_solve_phin(p);
	      for (auto* l : water->getBoundaryLines())
	        restore_post_solve_phin(l);

      if (qtf_radiation) {
        qtf_solutions[dof_col].push_back(bem_frequency_domain::capture_linear_solution(omega, face_sets.float_surface));
      }

	      if (radiation_debug) {
	        const auto rad_sol_all = bem_frequency_domain::capture_linear_solution(omega, all_boundary_faces);
	        debug_radiation_all_solutions.emplace(dof_col, rad_sol_all);
	        append_radiation_bc_residual(bc_residual_fs,
	                                     omega,
                                     dof_col,
                                     "real_gmres",
                                     rad_sol_all,
                                     face_sets,
                                     com,
                                     radiation_rotation_sign,
                                     [&](const BEM_DOF_Base&) { return Complex{kappa, 0.0}; });
      }

	      const auto S_col = integrate_reciprocity_from_mesh(face_sets.float_surface, com, radiation_moment_sign);
	      write_reciprocity_column(omega, dof_col, S_col, "real_gmres");
	      const auto Z_col = integrate_complex_pressure_force(*water, face_sets.float_surface, com, rho, omega, radiation_moment_sign);
      // Our convention: Z = i*omega*A - B  (e^{-i omega t}).
      for (int dof_row = 0; dof_row < 6; ++dof_row) {
        const double A = Z_col[dof_row].imag() / omega;
        const double B = -Z_col[dof_row].real();
        added_mass_raw[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] = A;
        damping_raw[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] = B;
        has_radiation_coeff[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] = true;
        raw_coeffs_fs << std::scientific << std::setprecision(12)
                      << omega << "," << dof_row << "," << dof_col << "," << A << "," << B << ",real_gmres,"
                      << (sponge.enabled ? 1 : 0) << "," << sponge.mu_max << ","
                      << sponge.r_start << "," << sponge.length << "," << sponge.order << "\n";
        zfs << std::scientific << std::setprecision(12)
            << omega << "," << dof_row << "," << dof_col << ","
            << Z_col[dof_row].real() << "," << Z_col[dof_row].imag() << ","
            << A << "," << B << ",real_gmres,"
            << (sponge.enabled ? 1 : 0) << "," << sponge.mu_max << ","
            << sponge.r_start << "," << sponge.length << "," << sponge.order << "\n";
      }

	      std::cout << "dof " << dof_col << ": setup=" << time_setup << " solve=" << time_solve << " residual=" << bvp.last_gmres_residual_norm << " iter=" << bvp.last_gmres_total_iter << std::endl;
	    }
	    if (radiation_debug) {
	      append_radiation_green_identity(outdir / "radiation_green_identity.csv",
	                                      omega,
	                                      radiation_complex_solver ? complex_lu_solver_name : "real_gmres",
	                                      debug_radiation_all_solutions,
	                                      all_boundary_faces,
	                                      face_sets);
	      append_radiation_effective_ntd_map(outdir / "radiation_effective_ntd_map.csv",
	                                         omega,
	                                         radiation_complex_solver ? complex_lu_solver_name : "real_gmres",
	                                         debug_radiation_all_solutions,
	                                         all_boundary_faces,
	                                         face_sets,
	                                         radiation_ref_point,
	                                         radiation_moment_sign);
	      append_radiation_effective_ntd_face_contrib(outdir / "radiation_effective_ntd_face_contrib.csv",
	                                                  omega,
	                                                  radiation_complex_solver ? complex_lu_solver_name : "real_gmres",
	                                                  debug_radiation_all_solutions,
	                                                  face_sets,
	                                                  radiation_ref_point,
	                                                  radiation_moment_sign);
	      if (radiation_complex_solver) {
	        append_radiation_algebraic_reciprocity(outdir / "radiation_algebraic_reciprocity.csv",
	                                               omega,
	                                               complex_lu_solver_name,
	                                               debug_linear_rhs_by_dof,
	                                               debug_linear_u_by_dof,
	                                               debug_radiation_all_solutions,
	                                               face_sets,
	                                               radiation_ref_point,
	                                               radiation_moment_sign);
		        append_radiation_reciprocity_pair_summary(outdir / "radiation_reciprocity_pair_summary.csv",
		                                                  omega,
		                                                  complex_lu_solver_name,
		                                                  debug_linear_rhs_by_dof,
		                                                  debug_linear_u_by_dof,
	                                                  debug_radiation_all_solutions,
	                                                  all_boundary_faces,
		                                                  face_sets,
		                                                  radiation_ref_point,
		                                                  radiation_moment_sign);
		        append_radiation_algebraic_reciprocity(outdir / "radiation_condensed_algebraic_reciprocity.csv",
		                                               omega,
		                                               "complex_lu_condensed_effective",
		                                               debug_condensed_rhs_by_dof,
		                                               debug_condensed_u_by_dof,
		                                               debug_radiation_all_solutions,
		                                               face_sets,
		                                               radiation_ref_point,
		                                               radiation_moment_sign);
		      }
		    }

	    added_mass = symmetrize_radiation_matrix6(added_mass_raw, has_radiation_coeff);
	    damping = symmetrize_radiation_matrix6(damping_raw, has_radiation_coeff);

	    struct RadiationCoeffAsymmetry {
	      int i = -1;
	      int j = -1;
	      double abs_value = 0.0;
	      double rel_value = 0.0;
	    };
	    RadiationCoeffAsymmetry worst_A;
	    RadiationCoeffAsymmetry worst_B;
	    const char* active_radiation_solver = radiation_complex_solver ? complex_lu_solver_name : "real_gmres";
	    auto write_coeff_symmetry_metric = [&](const char* metric, const Matrix6& raw, const Matrix6& sym, RadiationCoeffAsymmetry& worst) {
	      for (int i = 0; i < 6; ++i) {
	        for (int j = i; j < 6; ++j) {
	          const bool has_ij = has_radiation_coeff[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
	          const bool has_ji = has_radiation_coeff[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
	          if (!has_ij && !has_ji)
	            continue;
	          const double value_ij = has_ij ? raw[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] : 0.0;
	          const double value_ji = has_ji ? raw[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] : 0.0;
	          const double sym_value = sym[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
	          const double abs_asym = (has_ij && has_ji) ? std::abs(value_ij - value_ji) : 0.0;
	          const double rel_asym = (has_ij && has_ji) ? relative_pair_error(value_ij, value_ji) : 0.0;
	          if (i != j && has_ij && has_ji && abs_asym > worst.abs_value)
	            worst = {i, j, abs_asym, rel_asym};
	          coeff_symmetry_fs << std::scientific << std::setprecision(12)
	                            << omega << "," << i << "," << j << "," << metric << ","
	                            << value_ij << "," << value_ji << "," << sym_value << ","
	                            << abs_asym << "," << rel_asym << ","
	                            << (has_ij ? 1 : 0) << "," << (has_ji ? 1 : 0) << ","
	                            << active_radiation_solver << "\n";
	        }
	      }
	    };
	    write_coeff_symmetry_metric("A", added_mass_raw, added_mass, worst_A);
	    write_coeff_symmetry_metric("B", damping_raw, damping, worst_B);

	    for (int dof_col : dofs) {
	      if (dof_col < 0 || dof_col >= 6)
	        continue;
	      for (int dof_row = 0; dof_row < 6; ++dof_row) {
	        ofs << std::scientific << std::setprecision(12)
	            << omega << "," << dof_row << "," << dof_col << ","
	            << added_mass[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] << ","
	            << damping[static_cast<std::size_t>(dof_row)][static_cast<std::size_t>(dof_col)] << "\n";
	      }
	    }
	    if (worst_A.i >= 0 || worst_B.i >= 0) {
	      std::cout << "omega " << omega << ": radiation_coeffs.csv uses symmetrized A/B";
	      if (worst_A.i >= 0)
	        std::cout << "; worst raw A asym pair " << worst_A.i << worst_A.j
	                  << " abs=" << worst_A.abs_value << " rel=" << worst_A.rel_value;
	      if (worst_B.i >= 0)
	        std::cout << "; worst raw B asym pair " << worst_B.i << worst_B.j
	                  << " abs=" << worst_B.abs_value << " rel=" << worst_B.rel_value;
	      std::cout << std::endl;
	    }
	  }

  if (qtf_radiation) {
    for (const auto& [dof, solutions] : qtf_solutions) {
      if (solutions.empty())
        continue;
      const auto result = bem_frequency_domain::compute_qtf_result(solutions, face_sets.float_surface, float_body->COM, rho, qtf_symmetry);
      const std::string tag = "dof" + std::to_string(dof);
      const auto minus_path = outdir / ("qtf_minus_" + tag + ".csv");
      const auto plus_path = outdir / ("qtf_plus_" + tag + ".csv");
      bem_frequency_domain::write_qtf_csv(minus_path, result, false, qtf_symmetry);
      bem_frequency_domain::write_qtf_csv(plus_path, result, true, qtf_symmetry);
      if (qtf_newman) {
        const auto newman_path = outdir / ("qtf_minus_newman_" + tag + ".csv");
        bem_frequency_domain::write_qtf_newman_csv(newman_path, result, qtf_symmetry);
        const auto newman_report = bem_frequency_domain::qtf_newman_report(result);
        const auto newman_report_path = outdir / ("qtf_newman_" + tag + ".txt");
        std::ofstream nfs(newman_report_path);
        nfs << std::scientific << std::setprecision(6) << "max_abs=" << newman_report.max_abs << "\n"
            << "max_rel=" << newman_report.max_rel << "\n";
        std::cout << "wrote: " << newman_path << "\n"
                  << "       " << newman_report_path << std::endl;
      }
      if (qtf_check) {
        const auto report = bem_frequency_domain::qtf_symmetry_report(result);
        const auto report_path = outdir / ("qtf_symmetry_" + tag + ".txt");
        std::ofstream rfs(report_path);
        rfs << std::scientific << std::setprecision(6) << "max_abs_qminus=" << report.max_abs_qminus << "\n"
            << "max_abs_qplus=" << report.max_abs_qplus << "\n"
            << "max_sym_abs_qminus=" << report.max_sym_abs_qminus << "\n"
            << "max_sym_abs_qplus=" << report.max_sym_abs_qplus << "\n"
            << "max_sym_rel_qminus=" << report.max_sym_rel_qminus << "\n"
            << "max_sym_rel_qplus=" << report.max_sym_rel_qplus << "\n";
        std::cout << "wrote: " << report_path << std::endl;
      }
      std::cout << "wrote: " << minus_path << "\n"
                << "       " << plus_path << std::endl;
    }
  }

  if (compute_wave) {
    auto wave_base = water->water_wave_theory;
    if (qtf_unit_wave)
      wave_base.A = 1.0;
    if (wave_heading_deg)
      wave_base.theta = (*wave_heading_deg) / 180.0 * M_PI;
    if (!(wave_base.h > 0.0)) {
      const double z_min = std::get<0>(water->bounds[2]);
      wave_base.h = std::max(1e-6, z_free - z_min);
      wave_base.bottom_z = z_min;
    }
    if (!(wave_base.A > 0.0))
      throw std::runtime_error("frequency_domain_main: wave_theory amplitude is zero (use --qtf-unit-wave to force A=1)");

    const auto waterline_segments = bem_frequency_domain::collect_waterline_segments(*water, face_sets.float_surface, face_sets.free_surface);

	    std::vector<double> wave_omegas;
	    std::vector<bem_frequency_domain::LinearSolution> wave_incident_solutions;
	    std::vector<bem_frequency_domain::LinearSolution> wave_diffraction_solutions;
	    std::vector<bem_frequency_domain::LinearSolution> wave_solutions;
	    std::vector<bem_frequency_domain::Solution> wave_scat_solutions;
	    wave_omegas.reserve(omegas.size());
	    wave_incident_solutions.reserve(omegas.size());
	    wave_diffraction_solutions.reserve(omegas.size());
	    wave_solutions.reserve(omegas.size());
	    wave_scat_solutions.reserve(omegas.size());

    for (double omega : omegas) {
      if (!(omega > 0.0))
        continue;
      auto wave = bem_frequency_domain::make_wave_for_omega(wave_base, omega);
      wave.A = wave_base.A;
      wave.theta = wave_base.theta;
      wave.phase_shift = wave_base.phase_shift;
      wave.bottom_z = wave_base.bottom_z;

      const auto inc = bem_frequency_domain::build_incident_solution(wave, face_sets.float_surface);

		      bem_frequency_domain::LinearFSBC fsbc;
		      fsbc.omega = omega;
		      fsbc.gravity = g;
		      fsbc.sponge_mu = make_sponge_mu(omega);
			      auto average_normal_at = [&](const BEM_DOF_Base& node, const auto& accept_face, const Tddd& fallback) -> Tddd {
			        Tddd normal = {0.0, 0.0, 0.0};
			        for (auto* adj : node.getBoundaryFaces()) {
		          if (adj && accept_face(adj))
		            normal += adj->normal;
		        }
		        const double nm = Norm(normal);
		        if (nm > 1e-15)
			          return normal / nm;
			        return fallback;
			      };

		      bem_frequency_domain::BoundaryData bc;
	      bc.face_bc = [&](const networkFace& f) -> bem_frequency_domain::FaceBC {
	        if (face_sets.free_surface.contains(const_cast<networkFace*>(&f)))
	          return bem_frequency_domain::FaceBC::Robin;
	        if (freq_open_boundary_sommerfeld && face_sets.outer_wall.contains(const_cast<networkFace*>(&f)))
	          return bem_frequency_domain::FaceBC::Robin;
	        return bem_frequency_domain::FaceBC::Neumann;
	      };
	      configure_open_boundary_robin(bc, fsbc, omega);
			      bc.neumann_phin = [&](const BEM_DOF_Base& node, const networkFace* f) -> bem_frequency_domain::Complex {
			        bool on_float = false;
		        if (f) {
		          on_float = face_sets.float_surface.contains(const_cast<networkFace*>(f));
		        } else {
		          for (auto* adj : node.getBoundaryFaces()) {
		            if (face_sets.float_surface.contains(adj)) {
		              on_float = true;
		            }
		          }
		        }
			        if (!on_float)
			          return bem_frequency_domain::Complex{0.0, 0.0};
		        const Tddd normal = f ? f->normal
		                              : average_normal_at(node,
		                                                  [&](const networkFace* face) {
		                                                    return face && face_sets.float_surface.contains(const_cast<networkFace*>(face));
		                                                  },
		                                                  {0.0, 0.0, 0.0});
		        return -bem_frequency_domain::incident_phin_hat(wave, node.getPosition(), normal);
	      };
	      bc.dirichlet_phi = [&](const BEM_DOF_Base&) -> bem_frequency_domain::Complex { return bem_frequency_domain::Complex{0.0, 0.0}; };
	      bc.robin_unknown_policy = radiation_robin_unknown_policy;
	      bc.robin_unknown_phi = radiation_robin_phi_unknown;
	      bc.linear_solver = wave_complex_fmm_gmres ? "fmm_gmres" : (radiation_complex_gmres ? "gmres" : "lu");
      bc.gmres_tol = complex_gmres_tol;
      bc.gmres_max_iter = complex_gmres_max_iter;
      bc.gmres_restart = complex_gmres_restart;
      bc.gmres_preconditioner = complex_gmres_preconditioner;
      bc.fmm_coordinate_scaling = !wave_fmm_no_coordinate_scaling;

	      const auto scat = bem_frequency_domain::solve_linear_bvp({water}, fsbc, bc);
	      if (wave_complex_fmm_gmres)
	        append_fmm_nearfield_summary(omega, -1, scat);
	      append_linear_solver_residual("wave", omega, -1, wave_solver_name, scat);
	      const auto scat_sol = bem_frequency_domain::capture_linear_solution(omega, scat, face_sets.float_surface);

		      const auto total = bem_frequency_domain::combine_linear_solutions(
		          omega, {{&inc, bem_frequency_domain::Complex{1.0, 0.0}}, {&scat_sol, bem_frequency_domain::Complex{1.0, 0.0}}});

	      if (freq_free_surface_animation_frames > 0) {
	        const auto inc_free = bem_frequency_domain::build_incident_solution(wave, face_sets.free_surface);
	        const auto scat_free = bem_frequency_domain::capture_linear_solution(omega, scat, face_sets.free_surface);
	        const auto total_free = bem_frequency_domain::combine_linear_solutions(
	            omega, {{&inc_free, bem_frequency_domain::Complex{1.0, 0.0}}, {&scat_free, bem_frequency_domain::Complex{1.0, 0.0}}});
	        const std::string omega_label = safe_number_label(omega);
	        std::vector<std::pair<std::string, const bem_frequency_domain::LinearSolution*>> components{
	            {"incident", &inc_free},
	            {"scattering", &scat_free},
	        };
	        for (int dof_col : dofs) {
	          const auto it = radiation_free_surface_solutions.find({omega, dof_col});
	          if (it != radiation_free_surface_solutions.end())
	            components.push_back({"radiation_dof" + std::to_string(dof_col), &it->second});
	        }
	        write_free_surface_eta_animation_components(
	            outdir / "free_surface_eta_animation",
	            "wave_complete_omega_" + omega_label,
	            omega,
	            total_free,
	            components,
	            [&](const Tddd& x) { return sponge_mu_at_position(x, omega); },
	            g,
	            freq_free_surface_animation_frames,
	            freq_free_surface_animation_scale);
	      }

	      wave_omegas.push_back(omega);
	      wave_incident_solutions.push_back(inc);
	      wave_diffraction_solutions.push_back(scat_sol);
	      wave_solutions.push_back(total);
	      wave_scat_solutions.push_back(scat);
	    }

	    if (!wave_solutions.empty()) {
	      std::vector<std::array<bem_frequency_domain::Complex, 6>> wave_excitation;
	      std::vector<std::array<bem_frequency_domain::Complex, 6>> wave_excitation_incident;
	      std::vector<std::array<bem_frequency_domain::Complex, 6>> wave_excitation_diffraction;
	      wave_excitation.reserve(wave_solutions.size());
	      wave_excitation_incident.reserve(wave_incident_solutions.size());
	      wave_excitation_diffraction.reserve(wave_diffraction_solutions.size());
	      const std::size_t n_wave_force = std::min({wave_solutions.size(),
	                                                 wave_incident_solutions.size(),
	                                                 wave_diffraction_solutions.size()});
	      for (std::size_t i = 0; i < n_wave_force; ++i) {
	        wave_excitation.push_back(bem_frequency_domain::integrate_linear_pressure_force(wave_solutions[i], face_sets.float_surface, float_body->COM, rho));
	        wave_excitation_incident.push_back(bem_frequency_domain::integrate_linear_pressure_force(wave_incident_solutions[i], face_sets.float_surface, float_body->COM, rho));
	        wave_excitation_diffraction.push_back(bem_frequency_domain::integrate_linear_pressure_force(wave_diffraction_solutions[i], face_sets.float_surface, float_body->COM, rho));
	      }
	      {
	        const auto excitation_path = outdir / "wave_excitation.csv";
	        std::ofstream wfs(excitation_path);
        wfs << "omega,fx_re,fx_im,fy_re,fy_im,fz_re,fz_im,mx_re,mx_im,my_re,my_im,mz_re,mz_im\n";
        const std::size_t n = std::min(wave_omegas.size(), wave_excitation.size());
        for (std::size_t i = 0; i < n; ++i) {
          const auto& f = wave_excitation[i];
          wfs << std::scientific << std::setprecision(12) << wave_omegas[i];
          for (int k = 0; k < 6; ++k) {
            wfs << "," << f[k].real() << "," << f[k].imag();
          }
          wfs << "\n";
	        }
	        std::cout << "wrote: " << excitation_path << std::endl;
	      }
	      {
	        const auto components_path = outdir / "wave_excitation_components.csv";
	        std::ofstream cfs(components_path);
	        cfs << "omega,component,dof,value_re,value_im,amp,phase,balance_abs,balance_rel\n";
	        const std::size_t n = std::min({wave_omegas.size(),
	                                         wave_excitation.size(),
	                                         wave_excitation_incident.size(),
	                                         wave_excitation_diffraction.size()});
	        auto write_component = [&](double omega,
	                                   const char* component,
	                                   int dof,
	                                   const bem_frequency_domain::Complex& z,
	                                   double balance_abs,
	                                   double balance_rel) {
	          cfs << std::scientific << std::setprecision(12)
	              << omega << "," << component << "," << dof << ","
	              << z.real() << "," << z.imag() << ","
	              << std::abs(z) << "," << std::atan2(z.imag(), z.real()) << ","
	              << balance_abs << "," << balance_rel << "\n";
	        };
	        for (std::size_t i = 0; i < n; ++i) {
	          for (int dof = 0; dof < 6; ++dof) {
	            const auto total = wave_excitation[i][static_cast<std::size_t>(dof)];
	            const auto incident = wave_excitation_incident[i][static_cast<std::size_t>(dof)];
	            const auto diffraction = wave_excitation_diffraction[i][static_cast<std::size_t>(dof)];
	            const auto balance = total - incident - diffraction;
	            const double scale = std::max({std::abs(total), std::abs(incident) + std::abs(diffraction), 1e-300});
	            const double balance_abs = std::abs(balance);
	            const double balance_rel = balance_abs / scale;
	            write_component(wave_omegas[i], "incident", dof, incident, 0.0, 0.0);
	            write_component(wave_omegas[i], "diffraction", dof, diffraction, 0.0, 0.0);
	            write_component(wave_omegas[i], "total", dof, total, 0.0, 0.0);
	            write_component(wave_omegas[i], "total_minus_incident_minus_diffraction", dof, balance, balance_abs, balance_rel);
	          }
	        }
	        std::cout << "wrote: " << components_path << std::endl;
	      }
	      {
	        const auto components_path = outdir / "wave_excitation_components_standard.csv";
	        std::ofstream cfs(components_path);
	        cfs << "omega,component,symbol,dof,value_re,value_im,amp,phase,legacy_component,balance_abs,balance_rel\n";
	        const std::size_t n = std::min({wave_omegas.size(),
	                                         wave_excitation.size(),
	                                         wave_excitation_incident.size(),
	                                         wave_excitation_diffraction.size()});
	        auto write_component = [&](double omega,
	                                   const char* component,
	                                   const char* symbol,
	                                   int dof,
	                                   const bem_frequency_domain::Complex& z,
	                                   const char* legacy_component,
	                                   double balance_abs,
	                                   double balance_rel) {
	          cfs << std::scientific << std::setprecision(12)
	              << omega << "," << component << "," << symbol << "," << dof << ","
	              << z.real() << "," << z.imag() << ","
	              << std::abs(z) << "," << std::atan2(z.imag(), z.real()) << ","
	              << legacy_component << "," << balance_abs << "," << balance_rel << "\n";
	        };
	        for (std::size_t i = 0; i < n; ++i) {
	          for (int dof = 0; dof < 6; ++dof) {
	            const auto excitation = wave_excitation[i][static_cast<std::size_t>(dof)];
	            const auto froude_krylov = wave_excitation_incident[i][static_cast<std::size_t>(dof)];
	            const auto scattering = wave_excitation_diffraction[i][static_cast<std::size_t>(dof)];
	            const auto closure = excitation - froude_krylov - scattering;
	            const double scale = std::max({std::abs(excitation), std::abs(froude_krylov) + std::abs(scattering), 1e-300});
	            const double balance_abs = std::abs(closure);
	            const double balance_rel = balance_abs / scale;
	            write_component(wave_omegas[i], "froude_krylov", "phi_I", dof, froude_krylov, "incident", 0.0, 0.0);
	            write_component(wave_omegas[i], "scattering", "phi_S", dof, scattering, "diffraction", 0.0, 0.0);
	            write_component(wave_omegas[i], "excitation", "phi_D=phi_I+phi_S", dof, excitation, "total", 0.0, 0.0);
	            write_component(wave_omegas[i], "closure", "phi_D-phi_I-phi_S", dof, closure, "total_minus_incident_minus_diffraction", balance_abs, balance_rel);
	          }
	        }
	        std::cout << "wrote: " << components_path << std::endl;
	      }
	      if (wave_excitation_group_probe) {
	        const std::size_t n = std::min({wave_omegas.size(),
	                                         wave_incident_solutions.size(),
	                                         wave_diffraction_solutions.size(),
	                                         wave_solutions.size()});
	        for (std::size_t i = 0; i < n; ++i) {
	          append_wave_excitation_group_probe(outdir,
	                                             wave_omegas[i],
	                                             wave_solver_name,
	                                             wave_incident_solutions[i],
	                                             wave_diffraction_solutions[i],
	                                             wave_solutions[i],
	                                             face_sets.float_surface,
	                                             float_body->COM,
	                                             rho);
	        }
	        std::cout << "wrote: " << (outdir / "wave_excitation_face_group_contrib.csv") << "\n"
	                  << "       " << (outdir / "wave_excitation_face_contrib.csv") << std::endl;
	      }

	      if (rao) {
        const auto rao_path = outdir / "rao.csv";
        const auto residual_path = outdir / "rao_residual.csv";
        const auto rao_raw_path = outdir / "rao_raw.csv";
        const auto residual_raw_path = outdir / "rao_residual_raw.csv";
        const auto rao_delta_path = outdir / "rao_sym_delta.csv";
        std::ofstream rfs(rao_path);
        std::ofstream rrfs(residual_path);
        std::ofstream raw_rfs(rao_raw_path);
        std::ofstream raw_rrfs(residual_raw_path);
        std::ofstream delta_rfs(rao_delta_path);
        rfs << "omega,dof,xi_re,xi_im,amp,phase\n";
        rrfs << "omega,dof,res_re,res_im,abs\n";
        raw_rfs << "omega,dof,xi_re,xi_im,amp,phase\n";
        raw_rrfs << "omega,dof,res_re,res_im,abs\n";
        delta_rfs << "omega,dof,xi_sym_re,xi_sym_im,xi_raw_re,xi_raw_im,delta_re,delta_im,delta_abs,delta_rel\n";

        const std::size_t n_out = std::min(wave_omegas.size(), wave_excitation.size());
        for (std::size_t i = 0; i < n_out; ++i) {
          const double omega = wave_omegas[i];
          const auto a_it = added_mass_by_omega.find(omega);
          const auto b_it = damping_by_omega.find(omega);
          if (a_it == added_mass_by_omega.end() || b_it == damping_by_omega.end())
            throw std::runtime_error("frequency_domain_main: missing radiation coefficients for RAO omega=" + std::to_string(omega));
          const auto a_raw_it = added_mass_raw_by_omega.find(omega);
          const auto b_raw_it = damping_raw_by_omega.find(omega);
          if (a_raw_it == added_mass_raw_by_omega.end() || b_raw_it == damping_raw_by_omega.end())
            throw std::runtime_error("frequency_domain_main: missing raw radiation coefficients for RAO omega=" + std::to_string(omega));

          CVector6 residual{};
          const CVector6 xi = solve_rao_displacement(omega, dofs, rao_mass_diag, rao_restoring, a_it->second, b_it->second, wave_excitation[i], residual);
          CVector6 raw_residual{};
          const CVector6 xi_raw = solve_rao_displacement(omega, dofs, rao_mass_diag, rao_restoring, a_raw_it->second, b_raw_it->second, wave_excitation[i], raw_residual);
          for (int dof : dofs) {
            const auto z = xi[static_cast<std::size_t>(dof)];
            rfs << std::scientific << std::setprecision(12) << omega << "," << dof << "," << z.real() << "," << z.imag() << "," << std::abs(z) << "," << std::atan2(z.imag(), z.real()) << "\n";
            const auto r = residual[static_cast<std::size_t>(dof)];
            rrfs << std::scientific << std::setprecision(12) << omega << "," << dof << "," << r.real() << "," << r.imag() << "," << std::abs(r) << "\n";
            const auto zr = xi_raw[static_cast<std::size_t>(dof)];
            raw_rfs << std::scientific << std::setprecision(12) << omega << "," << dof << "," << zr.real() << "," << zr.imag() << "," << std::abs(zr) << "," << std::atan2(zr.imag(), zr.real()) << "\n";
            const auto rr = raw_residual[static_cast<std::size_t>(dof)];
            raw_rrfs << std::scientific << std::setprecision(12) << omega << "," << dof << "," << rr.real() << "," << rr.imag() << "," << std::abs(rr) << "\n";
            const auto delta = z - zr;
            const double delta_rel = std::abs(delta) / std::max({std::abs(z), std::abs(zr), 1e-300});
            delta_rfs << std::scientific << std::setprecision(12)
                      << omega << "," << dof << "," << z.real() << "," << z.imag() << ","
                      << zr.real() << "," << zr.imag() << "," << delta.real() << "," << delta.imag() << ","
                      << std::abs(delta) << "," << delta_rel << "\n";
          }
        }
        std::cout << "wrote: " << rao_path << "\n"
                  << "       " << residual_path << "\n"
                  << "       " << rao_raw_path << "\n"
                  << "       " << residual_raw_path << "\n"
                  << "       " << rao_delta_path << std::endl;
      }

      if (qtf_wave) {
        const auto result = bem_frequency_domain::compute_qtf_result(wave_solutions, face_sets.float_surface, float_body->COM, rho, qtf_symmetry);
        auto result_total = result;
        if (!waterline_segments.empty()) {
          const auto wl_result = bem_frequency_domain::compute_waterline_qtf(wave_omegas, wave_scat_solutions, wave_base, waterline_segments, float_body->COM, rho, g, qtf_symmetry);
          for (std::size_t i = 0; i < result_total.qminus.size(); ++i) {
            for (std::size_t j = 0; j < result_total.qminus[i].size(); ++j) {
              for (std::size_t k = 0; k < 6; ++k) {
                result_total.qminus[i][j][k] += wl_result.qminus[i][j][k];
                result_total.qplus[i][j][k] += wl_result.qplus[i][j][k];
              }
            }
          }
          if (qtf_check) {
            const auto wl_minus_path = outdir / "qtf_minus_wave_wl.csv";
            const auto wl_plus_path = outdir / "qtf_plus_wave_wl.csv";
            bem_frequency_domain::write_qtf_csv(wl_minus_path, wl_result, false, qtf_symmetry);
            bem_frequency_domain::write_qtf_csv(wl_plus_path, wl_result, true, qtf_symmetry);
            std::cout << "wrote: " << wl_minus_path << "\n"
                      << "       " << wl_plus_path << std::endl;
          }
        }
        const auto minus_path = outdir / "qtf_minus_wave.csv";
        const auto plus_path = outdir / "qtf_plus_wave.csv";
        bem_frequency_domain::write_qtf_csv(minus_path, result_total, false, qtf_symmetry);
        bem_frequency_domain::write_qtf_csv(plus_path, result_total, true, qtf_symmetry);
        if (qtf_newman) {
          const auto newman_path = outdir / "qtf_minus_newman_wave.csv";
          bem_frequency_domain::write_qtf_newman_csv(newman_path, result_total, qtf_symmetry);
          const auto newman_report = bem_frequency_domain::qtf_newman_report(result_total);
          const auto newman_report_path = outdir / "qtf_newman_wave.txt";
          std::ofstream nfs(newman_report_path);
          nfs << std::scientific << std::setprecision(6) << "max_abs=" << newman_report.max_abs << "\n"
              << "max_rel=" << newman_report.max_rel << "\n";
          std::cout << "wrote: " << newman_path << "\n"
                    << "       " << newman_report_path << std::endl;
        }
        if (qtf_check) {
          const auto report = bem_frequency_domain::qtf_symmetry_report(result_total);
          const auto report_path = outdir / "qtf_symmetry_wave.txt";
          std::ofstream rfs(report_path);
          rfs << std::scientific << std::setprecision(6) << "max_abs_qminus=" << report.max_abs_qminus << "\n"
              << "max_abs_qplus=" << report.max_abs_qplus << "\n"
              << "max_sym_abs_qminus=" << report.max_sym_abs_qminus << "\n"
              << "max_sym_abs_qplus=" << report.max_sym_abs_qplus << "\n"
              << "max_sym_rel_qminus=" << report.max_sym_rel_qminus << "\n"
              << "max_sym_rel_qplus=" << report.max_sym_rel_qplus << "\n";
          std::cout << "wrote: " << report_path << std::endl;
        }
        std::cout << "wrote: " << minus_path << "\n"
                  << "       " << plus_path << std::endl;
      }
    }
  }

  guard.restore();
  std::cout << "wrote: " << (outdir / "radiation_coeffs.csv") << "\n"
            << "       " << (outdir / "radiation_coeffs_raw.csv") << "\n"
            << "       " << (outdir / "radiation_coeffs_symmetry.csv") << std::endl;
  return 0;
}
