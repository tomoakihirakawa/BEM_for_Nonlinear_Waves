#include "pch.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Global variables expected by legacy headers.
bool use_linear_element = false;
bool use_pseudo_quadratic_element = false;
bool use_true_quadratic_element = false;
enum class NodeRelocationMethod { none, ALE, interpolation };
enum class NodeRelocationSurface { linear, pseudo_quadratic, true_quadratic };
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
std::string nearfield_mode = "scalar";
int g_p2m_quadrature_points = 6;
double g_mac_theta = 0.25;

int fmm_max_level = 7;
int fmm_bucket_max_points = 50;

int time_step = 0;
double simulation_time = 0.0;

#define BEM
#include "Network.hpp"

#include "BEM_freqency_domain.hpp"
#include "BEM_inputfile_reader.hpp"
#include "BEM_qtf.hpp"
#include "BEM_setBoundaryTypes.hpp"
#include "BEM_solveBVP.hpp"

namespace {

using Complex = std::complex<double>;

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

struct FaceSets {
  std::unordered_set<networkFace*> free_surface;  // Robin
  std::unordered_set<networkFace*> float_surface; // body boundary (radiation forcing)
};

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
    if (on_outer_x || on_outer_y || on_bottom) {
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

Tddd velocity_unit_dof(int dof, const Tddd& x, const Tddd& com) {
  if (dof < 0 || dof > 5)
    return {0.0, 0.0, 0.0};
  if (dof <= 2) {
    Tddd v{0.0, 0.0, 0.0};
    v[dof] = 1.0;
    return v;
  }
  Tddd omega{0.0, 0.0, 0.0};
  omega[dof - 3] = 1.0;
  return Cross(omega, x - com);
}

double get_phi_on_face(const networkPoint* p, const networkFace* f) {
  if (!p)
    return 0.0;
  if (f && p->phiOnFace.contains(const_cast<networkFace*>(f)))
    return p->phiOnFace.at(const_cast<networkFace*>(f));
  if (p->phiOnFace.contains(nullptr))
    return p->phiOnFace.at(nullptr);
  return 0.0;
}

std::array<Complex, 6> integrate_complex_pressure_force(const Network& water, const std::unordered_set<networkFace*>& float_faces, const Tddd& float_com, double rho, double omega) {
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

      const Tddd rq = xq - float_com;
      const std::array<Complex, 3> tq = {
          rq[1] * fq[2] - rq[2] * fq[1],
          rq[2] * fq[0] - rq[0] * fq[2],
          rq[0] * fq[1] - rq[1] * fq[0],
      };
      out[3] += w * tq[0] * f->area;
      out[4] += w * tq[1] * f->area;
      out[5] += w * tq[2] * f->area;
    }
  }

  (void)water;
  return out;
}

} // namespace

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cerr << "usage: ./main_freq_domain <input_dir> [--omega w1 w2 ...] [--dofs i j ...]\n"
                 "  [--qtf] [--qtf-wave] [--qtf-unit-wave] [--qtf-newman] [--qtf-check] [--qtf-full]\n";
    return 2;
  }

  SimulationSettings setting(argv[1], SimulationSettings::DomainMode::Frequency);
  use_pseudo_quadratic_element = setting.bem.element.pseudo_quadratic;
  use_true_quadratic_element = setting.bem.element.true_quadratic;

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
  nearfield_mode = setting.bem.solver.nearfield_mode;
  fmm_max_level = setting.bem.solver.fmm_max_level;
  fmm_bucket_max_points = setting.bem.solver.fmm_bucket_max_points;
  milu_omega = setting.bem.solver.milu_omega;
  g_p2m_quadrature_points = setting.bem.solver.p2m_quadrature_points;
  g_mac_theta = setting.bem.solver.mac_theta;

  const auto omegas = parse_omegas(argc, argv, setting.frequency.omegas);
  const auto dofs = parse_dofs(argc, argv, setting.frequency.dofs);
  const bool qtf_radiation = has_flag(argc, argv, "--qtf");
  const bool qtf_wave = has_flag(argc, argv, "--qtf-wave");
  const bool qtf_unit_wave = has_flag(argc, argv, "--qtf-unit-wave");
  const bool qtf_newman = has_flag(argc, argv, "--qtf-newman");
  const bool qtf_check = has_flag(argc, argv, "--qtf-check");
  const bool qtf_full = has_flag(argc, argv, "--qtf-full") || qtf_newman || qtf_check;
  const bool qtf_symmetry = !qtf_full;

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

  const double rho = setting.common.water_density;
  const double g = setting.common.gravity;

  std::filesystem::path outdir = setting.common.output_directory / "frequency_domain";
  std::filesystem::create_directories(outdir);

  std::ofstream ofs(outdir / "radiation_coeffs.csv");
  ofs << "omega,dof_row,dof_col,A,B\n";

  std::unordered_map<int, std::vector<bem_frequency_domain::LinearSolution>> qtf_solutions;

  // Guard & override boundary-state flags (keep this binary self-contained).
  bem_frequency_domain::BoundaryStateGuard guard(setting.FluidObject);

  water->setGeometricPropertiesForce();
  float_body->setGeometricPropertiesForce();

  const auto face_sets = classify_faces_deepcwind(*water, *float_body);
  std::unordered_set<const networkPoint*> float_points;
  float_points.reserve(face_sets.float_surface.size() * 3);
  for (auto* f : face_sets.float_surface) {
    auto [p0, p1, p2] = f->getPoints();
    float_points.emplace(p0);
    float_points.emplace(p1);
    float_points.emplace(p2);
  }

  // Set all boundary faces as Neumann (unknown: phi).
  for (auto* f : water->getBoundaryFaces()) {
    f->Dirichlet = false;
    f->Neumann = true;
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
  for (auto* l : water->getLines()) {
    l->Dirichlet = false;
    l->Neumann = true;
    l->CORNER = false;
  }
  for (auto* p : water->getPoints()) {
    p->Dirichlet = false;
    p->Neumann = true;
    p->CORNER = false;
    setIsMultipleNode(p);
  }

  const std::size_t n = setNodeFaceIndices(setting.FluidObject);

  // Build per-index Robin mask.
  std::vector<unsigned char> is_robin(static_cast<std::size_t>(n), 0);
  const double z_free = std::get<1>(water->bounds[2]);
  for (auto* p : water->getBoundaryPoints()) {
    for (const auto& [f, i] : p->f2Index) {
      if (i < 0 || static_cast<std::size_t>(i) >= n)
        continue;
      bool robin = false;
      if (f) {
        robin = face_sets.free_surface.contains(f);
      } else {
        robin = is_close(p->X[2], z_free, 1e-6);
      }
      is_robin[static_cast<std::size_t>(i)] = robin ? 1 : 0;
    }
  }
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

    std::cout << "\n=== omega=" << omega << " kappa=" << kappa << " ===" << std::endl;

    // Compute selected columns of the radiation impedance.
    for (int dof_col : dofs) {
      // Set boundary values.
      for (auto* p : water->getPoints()) {
        p->phiOnFace.clear();
        p->phinOnFace.clear();
        p->phitOnFace.clear();
        p->phintOnFace.clear();
      }

      const Tddd com = float_body->COM;
      for (auto* p : water->getBoundaryPoints()) {
        for (const auto& [f, i] : p->f2Index) {
          // Unknown phi (initial guess).
          p->phiOnFace[f] = 0.0;
          p->phitOnFace[f] = 0.0;

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
            const Tddd v = velocity_unit_dof(dof_col, p->X, com);
            const Tddd nrm = f ? f->normal : p->getNormalNeumann_BEM();
            phin = Dot(v, nrm);
          } else {
            // tank walls/bottom: impermeable (phin=0), free surface Robin handled in postprocess.
            phin = 0.0;
          }

          p->phinOnFace[f] = phin;
          p->phintOnFace[f] = 0.0;
        }
      }

      // Initialize midpoint boundary values for true_quadratic.
      if (use_true_quadratic_element) {
        for (auto* l : water->getBoundaryLines()) {
          if (l->midpoint_index < 0)
            continue;
          auto [pA, pB] = l->getPoints();
          const Tddd midX = 0.5 * (pA->X + pB->X);
          // Check if edge is on float surface.
          bool on_float = true;
          for (auto* f : l->Faces) {
            if (!face_sets.float_surface.contains(f)) {
              on_float = false;
              break;
            }
          }
          if (on_float) {
            const Tddd v = velocity_unit_dof(dof_col, midX, com);
            // Average normal from adjacent faces.
            Tddd nrm = {0.0, 0.0, 0.0};
            for (auto* f : l->Faces)
              nrm += f->normal;
            double nm = Norm(nrm);
            if (nm > 1e-15)
              nrm /= nm;
            l->phiphin[1] = Dot(v, nrm);
          } else {
            l->phiphin[1] = 0.0;
          }
          l->phiphin[0] = 0.0; // initial guess for unknown phi
        }
      }

      BEM_BVP bvp(setting.FluidObject);
      bvp.matrix_size = static_cast<int>(n);

      double time_setup = 0.0, time_solve = 0.0;
      TimeWatch watch;

      auto postprocess = [&]() {
        for (auto& [isDirichlet, i, phi, phin] : bvp.cache_DorN_phi_phin) {
          if (isDirichlet)
            continue;
          if (i < 0 || static_cast<std::size_t>(i) >= is_robin.size())
            continue;
          if (!is_robin[static_cast<std::size_t>(i)])
            continue;
          phin = kappa * phi;
          bvp.cache_phin_val_D_by_index[i] = phin;
        }
      };

      bvp.solveGMRES(watch, time_setup, time_solve, postprocess);

      // Re-evaluate Robin BC at midpoints using the solved phi_mid.
      // solveSystemGMRES restores pre-solve known BCs, but for Robin edges
      // the correct phin is kappa * phi_solved, not the initial value.
      if (use_true_quadratic_element) {
        for (auto* l : water->getBoundaryLines()) {
          int i = l->midpoint_index;
          if (i < 0 || static_cast<std::size_t>(i) >= is_robin.size())
            continue;
          if (is_robin[static_cast<std::size_t>(i)])
            l->phiphin[1] = kappa * l->phiphin[0];
        }
      }

      if (qtf_radiation) {
        qtf_solutions[dof_col].push_back(bem_frequency_domain::capture_linear_solution(omega, face_sets.float_surface));
      }

      const auto Z_col = integrate_complex_pressure_force(*water, face_sets.float_surface, com, rho, omega);
      // Our convention: Z = i*omega*A - B  (e^{-i omega t}).
      for (int dof_row = 0; dof_row < 6; ++dof_row) {
        const double A = Z_col[dof_row].imag() / omega;
        const double B = -Z_col[dof_row].real();
        ofs << omega << "," << dof_row << "," << dof_col << "," << A << "," << B << "\n";
      }

      std::cout << "dof " << dof_col << ": setup=" << time_setup << " solve=" << time_solve << " residual=" << bvp.last_gmres_residual_norm << " iter=" << bvp.last_gmres_total_iter << std::endl;
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

  if (qtf_wave) {
    auto wave_base = water->water_wave_theory;
    if (qtf_unit_wave)
      wave_base.A = 1.0;
    if (!(wave_base.h > 0.0)) {
      const double z_min = std::get<0>(water->bounds[2]);
      wave_base.h = std::max(1e-6, z_free - z_min);
      wave_base.bottom_z = z_min;
    }
    if (!(wave_base.A > 0.0))
      throw std::runtime_error("frequency_domain_main: wave_theory amplitude is zero (use --qtf-unit-wave to force A=1)");

    const auto waterline_segments = bem_frequency_domain::collect_waterline_segments(*water, face_sets.float_surface, face_sets.free_surface);

    std::vector<double> wave_omegas;
    std::vector<bem_frequency_domain::LinearSolution> wave_solutions;
    std::vector<bem_frequency_domain::Solution> wave_scat_solutions;
    wave_omegas.reserve(omegas.size());
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

      bem_frequency_domain::BoundaryData bc;
      bc.face_bc = [&](const networkFace& f) -> bem_frequency_domain::FaceBC {
        if (face_sets.free_surface.contains(const_cast<networkFace*>(&f)))
          return bem_frequency_domain::FaceBC::Robin;
        return bem_frequency_domain::FaceBC::Neumann;
      };
      bc.neumann_phin = [&](const networkPoint& p, const networkFace* f) -> bem_frequency_domain::Complex {
        bool on_float = false;
        if (f) {
          on_float = face_sets.float_surface.contains(const_cast<networkFace*>(f));
        } else {
          on_float = float_points.contains(&p);
        }
        if (!on_float)
          return bem_frequency_domain::Complex{0.0, 0.0};
        const Tddd normal = f ? f->normal : p.getNormalNeumann_BEM();
        return -bem_frequency_domain::incident_phin_hat(wave, p.X, normal);
      };
      bc.dirichlet_phi = [&](const networkPoint&) -> bem_frequency_domain::Complex { return bem_frequency_domain::Complex{0.0, 0.0}; };

      const auto scat = bem_frequency_domain::solve_linear_bvp({water}, fsbc, bc);
      const auto scat_sol = bem_frequency_domain::capture_linear_solution(omega, scat, face_sets.float_surface);

      const auto total = bem_frequency_domain::combine_linear_solutions(omega, {{&inc, bem_frequency_domain::Complex{1.0, 0.0}}, {&scat_sol, bem_frequency_domain::Complex{1.0, 0.0}}});

      wave_omegas.push_back(omega);
      wave_solutions.push_back(total);
      wave_scat_solutions.push_back(scat);
    }

    if (!wave_solutions.empty()) {
      std::vector<std::array<bem_frequency_domain::Complex, 6>> wave_excitation;
      wave_excitation.reserve(wave_solutions.size());
      for (const auto& sol : wave_solutions) {
        wave_excitation.push_back(bem_frequency_domain::integrate_linear_pressure_force(sol, face_sets.float_surface, float_body->COM, rho));
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

  guard.restore();
  std::cout << "wrote: " << (outdir / "radiation_coeffs.csv") << std::endl;
  return 0;
}
