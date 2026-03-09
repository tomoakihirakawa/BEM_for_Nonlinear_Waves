#include "BEM_freqency_domain.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>

// -----------------------------------------------------------------------------
// Required global knobs (declared as extern in existing headers).
// Keep defaults so this test does not affect the legacy code paths.
// -----------------------------------------------------------------------------
std::string solver_type = "LU";
std::string coupling_type = "NONE";
std::string preconditioner_type = "NONE";
std::string ilu_neighborhood_type = "BUCKETS";
int ilu_kring_num = 1;
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
double coupling_tol = 1e-10;
double simulation_time = 0.0;
std::vector<double> coupling_params;
bool use_pseudo_quadratic_element = false;

int main() {
  try {
    const std::filesystem::path tank_obj = std::filesystem::path(__FILE__).parent_path() / "../../obj/Ruehl2016/tank.obj";
    Network tank(tank_obj.string(), "tank");
    tank.setGeometricPropertiesForce();

    // Bounds
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();
    for (const auto *p : tank.getPoints()) {
      min_x = std::min(min_x, std::get<0>(p->X));
      max_x = std::max(max_x, std::get<0>(p->X));
      min_z = std::min(min_z, std::get<2>(p->X));
      max_z = std::max(max_z, std::get<2>(p->X));
    }
    const double Lx = std::max(1e-12, max_x - min_x);
    const double Lz = std::max(1e-12, max_z - min_z);
    const double tol_x = 1e-8 * Lx + 1e-12;
    const double tol_z = 1e-8 * Lz + 1e-12;

    using namespace bem_frequency_domain;

    LinearFSBC fsbc;
    fsbc.omega = 2.0;      // rad/s
    fsbc.gravity = 9.81;   // m/s^2
    fsbc.sponge_mu = {};   // no sponge

    BoundaryData bc;
    bc.face_bc = [&](const networkFace &f) -> FaceBC {
      const double zc = std::get<2>(f.centroid);
      if (std::abs(zc - max_z) < tol_z)
        return FaceBC::Robin; // free surface (Robin)
      return FaceBC::Neumann; // walls/bottom
    };

    // Neumann BC: drive one side as a simple wavemaker-like flux.
    bc.neumann_phin = [&](const networkPoint &p, const networkFace *f) -> Complex {
      const double x = (f != nullptr) ? std::get<0>(f->centroid) : std::get<0>(p.X);
      if (std::abs(x - min_x) < tol_x)
        return Complex{1.0, 0.0};
      return Complex{0.0, 0.0};
    };

    // No true Dirichlet faces in this test.
    bc.dirichlet_phi = [&](const networkPoint &) -> Complex { return Complex{0.0, 0.0}; };

    std::cout << "Solving linear frequency-domain BVP on: " << tank_obj << "\n"
              << "  omega=" << fsbc.omega << "  g=" << fsbc.gravity << "\n"
              << "  bounds: x=[" << min_x << "," << max_x << "], z=[" << min_z << "," << max_z << "]" << std::endl;

    auto sol = solve_linear_bvp({&tank}, fsbc, bc);

    double max_abs_phi = 0.0, max_abs_phin = 0.0;
    for (std::size_t i = 0; i < sol.n; ++i) {
      max_abs_phi = std::max(max_abs_phi, std::abs(sol.phi[i]));
      max_abs_phin = std::max(max_abs_phin, std::abs(sol.phin[i]));
    }

    std::cout << "n=" << sol.n << "\n"
              << "max|phi|=" << max_abs_phi << "\n"
              << "max|phin|=" << max_abs_phin << "\n"
              << "BIE residual L2=" << sol.bie_residual_l2 << std::endl;

    // Minimal sanity checks
    if (!std::isfinite(max_abs_phi) || !std::isfinite(max_abs_phin) || !std::isfinite(sol.bie_residual_l2)) {
      std::cerr << "Non-finite result detected." << std::endl;
      return 2;
    }
    if (max_abs_phi == 0.0 && max_abs_phin == 0.0) {
      std::cerr << "Trivial (all-zero) solution; expected non-zero response." << std::endl;
      return 3;
    }

    std::cout << "OK" << std::endl;
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
}
