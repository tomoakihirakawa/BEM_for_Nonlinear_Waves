/*

cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=main_time_domain.cpp

*/
#include "pch.hpp"
#include <cctype>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <execinfo.h>
#include <string>
#include <unistd.h>

/*DOC_EXTRACT 0_0_BEM

# BEM-MEL

<img src="./sample_Goring1979.gif">

*/

// #define _debugging_

bool use_linear_element = false;
bool use_pseudo_quadratic_element = false;
bool use_true_quadratic_element = false;
enum class NodeRelocationMethod { none, ALE, interpolation };
enum class NodeRelocationSurface { linear, pseudo_quadratic, true_quadratic };
NodeRelocationMethod node_relocation_method = NodeRelocationMethod::none;
NodeRelocationSurface node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
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
double g_mac_theta = 0.25;

int time_step;
double simulation_time = 0;

#define BEM
#define simulation
//
#include "Network.hpp"

JSONoutput jsonout;

// pvd cpg_pvd("./vtu/bem.pvd");

#include "BEM.hpp"
#include "BEM_inputfile_reader.hpp"
// 追加
#include "OutputCommon.hpp"
#include "OutputJSON.hpp"
#include "OutputParaview.hpp"

#include "BEM_internal_flow.hpp"
#include "VPM.hpp"
#include "BEM_checkpoint.hpp"

#include "BEM_remesh_main.hpp"

namespace {
void crashBacktraceHandler(int sig) {
  const char* sig_name = strsignal(sig);
  if (!sig_name)
    sig_name = "UNKNOWN";

  ::write(STDERR_FILENO, "\n=== crash ===\n", 13);
  ::write(STDERR_FILENO, sig_name, std::strlen(sig_name));
  ::write(STDERR_FILENO, "\n", 1);

  void* frames[128];
  const int n = ::backtrace(frames, 128);
  ::backtrace_symbols_fd(frames, n, STDERR_FILENO);
  ::write(STDERR_FILENO, "=== end ===\n", 12);

  std::_Exit(128 + sig);
}

void installCrashBacktraceIfRequested() {
  const char* env = std::getenv("BEM_BACKTRACE");
  if (!env || std::strcmp(env, "1") != 0)
    return;

  std::signal(SIGSEGV, crashBacktraceHandler);
  std::signal(SIGABRT, crashBacktraceHandler);
  std::signal(SIGBUS, crashBacktraceHandler);
}

void enableRealFieldFmmDefaultForTimeDomain() {
  // Time-domain boundary data/solutions are typically real.
  // Prefer the real-field (m<0 conjugate symmetry) optimization by default for time-domain runs,
  // while still allowing users to override via `BEM_FMM_REALFIELD_M_CONJ=0/1`.
  if (std::getenv("BEM_FMM_REALFIELD_M_CONJ") != nullptr)
    return;
#if defined(_WIN32)
  _putenv_s("BEM_FMM_REALFIELD_M_CONJ", "1");
#else
  ::setenv("BEM_FMM_REALFIELD_M_CONJ", "1", 0 /* overwrite */);
#endif
}

struct PenetrationReport {
  std::size_t point_count = 0;
  std::size_t face_count = 0;
  double max_distance = 0.0;
  double worst_tolerance = 0.0;
  double worst_local_length = 0.0;
  std::string worst_body;
};

double representativeLength(const networkPoint* p) {
  if (!p)
    return 0.0;
  double sum = 0.0;
  std::size_t count = 0;
  for (const auto* l : p->getLines()) {
    if (!l)
      continue;
    const double len = l->length();
    if (!(len > 0.0))
      continue;
    sum += len;
    ++count;
  }
  if (count > 0)
    return sum / static_cast<double>(count);
  if (p->contact_range > 0.0)
    return p->contact_range;
  return 0.0;
}

double representativeLength(const networkLine* l) {
  if (!l)
    return 0.0;
  const double len = l->length();
  if (len > 0.0)
    return len;
  const auto [p0, p1] = l->getPoints();
  return std::max(representativeLength(p0), representativeLength(p1));
}

double penetrationTolerance(const double h_local) {
  constexpr double c_rel = 0.25;
  constexpr double c_abs_eps = 1e-6;
  return (h_local > 0.0) ? c_rel * h_local + c_abs_eps : c_abs_eps;
}

std::vector<Network*> penetrationCandidateSolids(const std::vector<Network*>& fluids,
                                                 const std::vector<Network*>& solids) {
  static std::unordered_map<const Network*, bool> enclosure_cache;
  static std::unordered_set<const Network*> logged;

  std::vector<Network*> out;
  out.reserve(solids.size());

  for (auto* solid : solids) {
    if (!solid)
      continue;

    auto [it, inserted] = enclosure_cache.emplace(solid, false);
    if (inserted) {
      std::size_t inside_count = 0;
      std::size_t total_count = 0;
      for (const auto* fluid : fluids) {
        if (!fluid)
          continue;
        for (const auto* p : fluid->getBoundaryPoints()) {
          ++total_count;
          if (solid->InsideQ(p->X))
            ++inside_count;
        }
      }
      const double inside_fraction =
          (total_count > 0) ? static_cast<double>(inside_count) / static_cast<double>(total_count) : 0.0;
      it->second = inside_fraction > 0.2;

      if (!logged.contains(solid)) {
        std::cout << Yellow << "[structure_contact] body " << solid->getName()
                  << ": inside_fraction=" << inside_fraction
                  << (it->second ? " -> enclosure, excluded from penetration rejection"
                                 : " -> obstacle, monitored for penetration")
                  << colorReset << std::endl;
        logged.insert(solid);
      }
    }

    if (!it->second)
      out.push_back(solid);
  }

  return out;
}

[[nodiscard]] PenetrationReport detectStructurePenetration(const std::vector<Network*>& fluids,
                                                           const std::vector<Network*>& solids) {
  PenetrationReport report;
  for (const auto* water : fluids) {
    for (const auto* p : water->getBoundaryPoints()) {
      const double h_local = representativeLength(p);
      const double tolerance = penetrationTolerance(h_local);
      for (const auto* solid : solids) {
        if (!solid->InsideQ(p->X))
          continue;
        auto [near_f, near_x] = solid->Nearest(p->X);
        const double dist = near_f ? Norm(p->X - near_x) : 0.0;
        if (dist <= tolerance)
          continue;
        ++report.point_count;
        if (dist > report.max_distance) {
          report.max_distance = dist;
          report.worst_tolerance = tolerance;
          report.worst_local_length = h_local;
          report.worst_body = solid->getName();
        }
        break;
      }
    }
    if (use_true_quadratic_element) {
      // Midpoint DOFs exist only on true-quadratic elements.
      for (const auto* l : water->getBoundaryLines()) {
        const bool has_true_quad = std::ranges::any_of(l->getBoundaryFaces(), [](const auto* f) { return f->isTrueQuadraticElement; });
        if (!has_true_quad)
          continue;
        const double h_local = representativeLength(l);
        const double tolerance = penetrationTolerance(h_local);
        for (const auto* solid : solids) {
          if (!solid->InsideQ(l->X_mid))
            continue;
          auto [near_f, near_x] = solid->Nearest(l->X_mid);
          const double dist = near_f ? Norm(l->X_mid - near_x) : 0.0;
          if (dist <= tolerance)
            continue;
          ++report.face_count;
          if (dist > report.max_distance) {
            report.max_distance = dist;
            report.worst_tolerance = tolerance;
            report.worst_local_length = h_local;
            report.worst_body = solid->getName();
          }
          break;
        }
      }
    }
  }
  return report;
}

void throwIfStructurePenetrated(const std::vector<Network*>& fluids,
                                const std::vector<Network*>& solids,
                                const int time_step,
                                const std::string& stage) {
  const auto active_solids = penetrationCandidateSolids(fluids, solids);
  if (active_solids.empty())
    return;
  const auto report = detectStructurePenetration(fluids, active_solids);
  if (report.point_count == 0 && report.face_count == 0)
    return;
  throw step_failure("structure penetration detected at " + stage +
                     " on time_step " + std::to_string(time_step) +
                     " (points=" + std::to_string(report.point_count) +
                     ", mids=" + std::to_string(report.face_count) +
                     ", max_distance=" + std::to_string(report.max_distance) +
                     (report.worst_body.empty() ? "" : ", body=" + report.worst_body) +
                     ", tolerance=" + std::to_string(report.worst_tolerance) +
                     ", h_local=" + std::to_string(report.worst_local_length) + ")");
}
} // namespace

// ---------------------------------------------------------------------------
// Export kernel scatter data: G and dG/dn vs distance for a representative
// water surface node. Called once at time_step == 0 for error analysis.
// ---------------------------------------------------------------------------
// Dunavant quadrature rules for triangles (barycentric coordinates)
// Format: {lambda0, lambda1, lambda2, weight}, weights sum to 1.0
// Integral: area * sum(w_i * f(lambda_i * X))
static constexpr std::array<std::array<double, 4>, 1> __dunavant_1__ = {{{1. / 3., 1. / 3., 1. / 3., 1.0}}};
static constexpr std::array<std::array<double, 4>, 3> __dunavant_3__ = {{{2. / 3., 1. / 6., 1. / 6., 1. / 3.},
                                                                         {1. / 6., 2. / 3., 1. / 6., 1. / 3.},
                                                                         {1. / 6., 1. / 6., 2. / 3., 1. / 3.}}};
// Dunavant / Hammer-Marlowe-Stroud 4-point rule, degree 3 (n=4)
// Note: centroid weight is negative, which is mathematically valid
static constexpr std::array<std::array<double, 4>, 4> __dunavant_4__ = {{{1. / 3., 1. / 3., 1. / 3., -27. / 48.},
                                                                         {3. / 5., 1. / 5., 1. / 5., 25. / 48.},
                                                                         {1. / 5., 3. / 5., 1. / 5., 25. / 48.},
                                                                         {1. / 5., 1. / 5., 3. / 5., 25. / 48.}}};
// Dunavant 6-point rule, degree 4 (n=5), all positive weights
static constexpr std::array<std::array<double, 4>, 6> __dunavant_6__ = {{{0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322},
                                                                         {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.109951743655322},
                                                                         {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.109951743655322},
                                                                         {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011},
                                                                         {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.223381589678011},
                                                                         {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.223381589678011}}};
static constexpr std::array<std::array<double, 4>, 7> __dunavant_7__ = {{{1. / 3., 1. / 3., 1. / 3., 0.225},
                                                                         {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506},
                                                                         {0.470142064105115, 0.059715871789770, 0.470142064105115, 0.132394152788506},
                                                                         {0.470142064105115, 0.470142064105115, 0.059715871789770, 0.132394152788506},
                                                                         {0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827},
                                                                         {0.101286507323456, 0.797426985353087, 0.101286507323456, 0.125939180544827},
                                                                         {0.101286507323456, 0.101286507323456, 0.797426985353087, 0.125939180544827}}};
// D12: 12-point rule (degree 6), weights sum to 1.0
static constexpr std::array<std::array<double, 4>, 12> __dunavant_12__ = {{{0.501426509658179, 0.249286745170910, 0.249286745170910, 0.116786275726379},
                                                                           {0.249286745170910, 0.501426509658179, 0.249286745170910, 0.116786275726379},
                                                                           {0.249286745170910, 0.249286745170910, 0.501426509658179, 0.116786275726379},
                                                                           {0.873821971016996, 0.063089014491502, 0.063089014491502, 0.050844906370207},
                                                                           {0.063089014491502, 0.873821971016996, 0.063089014491502, 0.050844906370207},
                                                                           {0.063089014491502, 0.063089014491502, 0.873821971016996, 0.050844906370207},
                                                                           {0.053145049844817, 0.310352451033784, 0.636502499121399, 0.082851075618374},
                                                                           {0.053145049844817, 0.636502499121399, 0.310352451033784, 0.082851075618374},
                                                                           {0.310352451033784, 0.053145049844817, 0.636502499121399, 0.082851075618374},
                                                                           {0.310352451033784, 0.636502499121399, 0.053145049844817, 0.082851075618374},
                                                                           {0.636502499121399, 0.053145049844817, 0.310352451033784, 0.082851075618374},
                                                                           {0.636502499121399, 0.310352451033784, 0.053145049844817, 0.082851075618374}}};
// D13: 13-point rule (degree 7), weights sum to 1.0 (negative center weight)
static constexpr std::array<std::array<double, 4>, 13> __dunavant_13__ = {{{1. / 3., 1. / 3., 1. / 3., -0.149570044467682},
                                                                           {0.479308067841920, 0.260345966079040, 0.260345966079040, 0.175615257433208},
                                                                           {0.260345966079040, 0.479308067841920, 0.260345966079040, 0.175615257433208},
                                                                           {0.260345966079040, 0.260345966079040, 0.479308067841920, 0.175615257433208},
                                                                           {0.869739794195568, 0.065130102902216, 0.065130102902216, 0.053347235608838},
                                                                           {0.065130102902216, 0.869739794195568, 0.065130102902216, 0.053347235608838},
                                                                           {0.065130102902216, 0.065130102902216, 0.869739794195568, 0.053347235608838},
                                                                           {0.048690315425316, 0.312865496004874, 0.638444188569810, 0.077113760890257},
                                                                           {0.048690315425316, 0.638444188569810, 0.312865496004874, 0.077113760890257},
                                                                           {0.312865496004874, 0.048690315425316, 0.638444188569810, 0.077113760890257},
                                                                           {0.312865496004874, 0.638444188569810, 0.048690315425316, 0.077113760890257},
                                                                           {0.638444188569810, 0.048690315425316, 0.312865496004874, 0.077113760890257},
                                                                           {0.638444188569810, 0.312865496004874, 0.048690315425316, 0.077113760890257}}};
// D16: 16-point rule (degree 8), weights sum to 1.0
static constexpr std::array<std::array<double, 4>, 16> __dunavant_16__ = {{{1. / 3., 1. / 3., 1. / 3., 0.144315607677787},
                                                                           {0.081414823414554, 0.459292588292723, 0.459292588292723, 0.095091634267285},
                                                                           {0.459292588292723, 0.081414823414554, 0.459292588292723, 0.095091634267285},
                                                                           {0.459292588292723, 0.459292588292723, 0.081414823414554, 0.095091634267285},
                                                                           {0.658861384496480, 0.170569307751760, 0.170569307751760, 0.103217370534718},
                                                                           {0.170569307751760, 0.658861384496480, 0.170569307751760, 0.103217370534718},
                                                                           {0.170569307751760, 0.170569307751760, 0.658861384496480, 0.103217370534718},
                                                                           {0.898905543365938, 0.050547228317031, 0.050547228317031, 0.032458497623198},
                                                                           {0.050547228317031, 0.898905543365938, 0.050547228317031, 0.032458497623198},
                                                                           {0.050547228317031, 0.050547228317031, 0.898905543365938, 0.032458497623198},
                                                                           {0.008394777409958, 0.263112829634638, 0.728492392955404, 0.027230314174435},
                                                                           {0.008394777409958, 0.728492392955404, 0.263112829634638, 0.027230314174435},
                                                                           {0.263112829634638, 0.008394777409958, 0.728492392955404, 0.027230314174435},
                                                                           {0.263112829634638, 0.728492392955404, 0.008394777409958, 0.027230314174435},
                                                                           {0.728492392955404, 0.008394777409958, 0.263112829634638, 0.027230314174435},
                                                                           {0.728492392955404, 0.263112829634638, 0.008394777409958, 0.027230314174435}}};
// D19: 19-point rule (degree 9), weights sum to 1.0
static constexpr std::array<std::array<double, 4>, 19> __dunavant_19__ = {{{1. / 3., 1. / 3., 1. / 3., 0.097135796282799},
                                                                           {0.020634961602525, 0.489682519198738, 0.489682519198738, 0.031334700227139},
                                                                           {0.489682519198738, 0.020634961602525, 0.489682519198738, 0.031334700227139},
                                                                           {0.489682519198738, 0.489682519198738, 0.020634961602525, 0.031334700227139},
                                                                           {0.125820817014127, 0.437089591492937, 0.437089591492937, 0.077827541004774},
                                                                           {0.437089591492937, 0.125820817014127, 0.437089591492937, 0.077827541004774},
                                                                           {0.437089591492937, 0.437089591492937, 0.125820817014127, 0.077827541004774},
                                                                           {0.623592928761935, 0.188203535619033, 0.188203535619033, 0.079647738927210},
                                                                           {0.188203535619033, 0.623592928761935, 0.188203535619033, 0.079647738927210},
                                                                           {0.188203535619033, 0.188203535619033, 0.623592928761935, 0.079647738927210},
                                                                           {0.910540973211095, 0.044729513394453, 0.044729513394453, 0.025577675658698},
                                                                           {0.044729513394453, 0.910540973211095, 0.044729513394453, 0.025577675658698},
                                                                           {0.044729513394453, 0.044729513394453, 0.910540973211095, 0.025577675658698},
                                                                           {0.036838412054736, 0.221962989160766, 0.741198598784498, 0.043283539377289},
                                                                           {0.036838412054736, 0.741198598784498, 0.221962989160766, 0.043283539377289},
                                                                           {0.221962989160766, 0.036838412054736, 0.741198598784498, 0.043283539377289},
                                                                           {0.221962989160766, 0.741198598784498, 0.036838412054736, 0.043283539377289},
                                                                           {0.741198598784498, 0.036838412054736, 0.221962989160766, 0.043283539377289},
                                                                           {0.741198598784498, 0.221962989160766, 0.036838412054736, 0.043283539377289}}};
// D25: 25-point rule (degree 10), weights sum to 1.0
static constexpr std::array<std::array<double, 4>, 25> __dunavant_25__ = {{{1. / 3., 1. / 3., 1. / 3., 0.090817990382754},
                                                                           {0.028844733232685, 0.485577633383657, 0.485577633383657, 0.036725957756467},
                                                                           {0.485577633383657, 0.028844733232685, 0.485577633383657, 0.036725957756467},
                                                                           {0.485577633383657, 0.485577633383657, 0.028844733232685, 0.036725957756467},
                                                                           {0.781036849029926, 0.109481575485037, 0.109481575485037, 0.045321059435528},
                                                                           {0.109481575485037, 0.781036849029926, 0.109481575485037, 0.045321059435528},
                                                                           {0.109481575485037, 0.109481575485037, 0.781036849029926, 0.045321059435528},
                                                                           {0.141707219414880, 0.307939838764121, 0.550352941820999, 0.072757916845420},
                                                                           {0.141707219414880, 0.550352941820999, 0.307939838764121, 0.072757916845420},
                                                                           {0.307939838764121, 0.141707219414880, 0.550352941820999, 0.072757916845420},
                                                                           {0.307939838764121, 0.550352941820999, 0.141707219414880, 0.072757916845420},
                                                                           {0.550352941820999, 0.141707219414880, 0.307939838764121, 0.072757916845420},
                                                                           {0.550352941820999, 0.307939838764121, 0.141707219414880, 0.072757916845420},
                                                                           {0.025003534762686, 0.246672560639903, 0.728323904597411, 0.028327242531057},
                                                                           {0.025003534762686, 0.728323904597411, 0.246672560639903, 0.028327242531057},
                                                                           {0.246672560639903, 0.025003534762686, 0.728323904597411, 0.028327242531057},
                                                                           {0.246672560639903, 0.728323904597411, 0.025003534762686, 0.028327242531057},
                                                                           {0.728323904597411, 0.025003534762686, 0.246672560639903, 0.028327242531057},
                                                                           {0.728323904597411, 0.246672560639903, 0.025003534762686, 0.028327242531057},
                                                                           {0.009540815400299, 0.066803251012200, 0.923655933587500, 0.009421666963733},
                                                                           {0.009540815400299, 0.923655933587500, 0.066803251012200, 0.009421666963733},
                                                                           {0.066803251012200, 0.009540815400299, 0.923655933587500, 0.009421666963733},
                                                                           {0.066803251012200, 0.923655933587500, 0.009540815400299, 0.009421666963733},
                                                                           {0.923655933587500, 0.009540815400299, 0.066803251012200, 0.009421666963733},
                                                                           {0.923655933587500, 0.066803251012200, 0.009540815400299, 0.009421666963733}}};

void exportKernelScatterData(const std::vector<Network*>& FluidObject,
                             const std::filesystem::path& output_directory) {
  if (FluidObject.empty())
    return;

  auto* water = FluidObject[0];
  auto surfaces = water->getBoundaryFaces();
  auto points = water->getBoundaryPoints();

  if (surfaces.empty() || points.empty())
    return;

  // Characteristic element size: h = sqrt(mean face area)
  double total_area = 0;
  for (const auto& f : surfaces)
    total_area += f->area;
  const double h = std::sqrt(total_area / surfaces.size());

  // Find target point closest to the centroid of the water surface
  Tddd centroid = {0., 0., 0.};
  for (const auto& p : points)
    centroid = centroid + p->X;
  centroid = centroid / static_cast<double>(points.size());

  networkPoint* target = nullptr;
  double min_dist = 1e20;
  for (const auto& p : points) {
    double d = Norm(p->X - centroid);
    if (d < min_dist) {
      min_dist = d;
      target = p;
    }
  }
  if (!target)
    return;

  std::cout << "[exportKernelScatterData] target: (" << target->X[0] << ", " << target->X[1] << ", " << target->X[2]
            << "), h = " << h << ", N_faces = " << surfaces.size() << std::endl;

  std::ofstream ofs(output_directory / "kernel_scatter.csv");
  ofs << std::setprecision(15);
  ofs << "face_index,type,is_singular,r,r_over_h,"
      << "intG_gw1,intGn_gw1,intG_gw3,intGn_gw3,intG_gw5,intGn_gw5,intG_gw11,intGn_gw11,"
      << "intG_d1,intGn_d1,intG_d3,intGn_d3,intG_d4,intGn_d4,intG_d6,intGn_d6,intG_d7,intGn_d7,"
      << "intG_d12,intGn_d12,intG_d13,intGn_d13,intG_d16,intGn_d16,intG_d19,intGn_d19,intG_d25,intGn_d25,"
      << "area,area_over_h2\n";

  constexpr double inv4pi = 1.0 / (4.0 * M_PI);

  int face_idx = 0;
  for (const auto& f : surfaces) {
    auto [p0, p1, p2] = f->getPoints();
    const Tddd X0 = p0->X, X1 = p1->X, X2 = p2->X;

    const Tddd cross = Cross(X1 - X0, X2 - X0);
    const double J_det = Norm(cross);

    const double r = Norm((X0 + X1 + X2) / 3. - target->X);

    // Detect singular face (target is a vertex of this face)
    const bool is_singular = (p0 == target || p1 == target || p2 == target);

    // --- Tensor product GW quadrature (square-to-triangle mapping) ---
    // intG = sum ww*(1-t0) * J_det * G,  intGn = -sum ww*(1-t0) * dot(R, cross) * Gn_factor
    double intG_gw1 = 0., intGn_gw1 = 0.;
    for (const auto& [t0, t1, ww] : __array_GW1xGW1__) {
      const double weight = ww * (1. - t0);
      const auto N = ModTriShape<3>(t0, t1);
      const Tddd X = N[0] * X0 + N[1] * X1 + N[2] * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_gw1 += weight * J_det * inv4pi / nr;
      intGn_gw1 -= weight * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_gw3 = 0., intGn_gw3 = 0.;
    for (const auto& [t0, t1, ww] : __array_GW3xGW3__) {
      const double weight = ww * (1. - t0);
      const auto N = ModTriShape<3>(t0, t1);
      const Tddd X = N[0] * X0 + N[1] * X1 + N[2] * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_gw3 += weight * J_det * inv4pi / nr;
      intGn_gw3 -= weight * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_gw5 = 0., intGn_gw5 = 0.;
    for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
      const double weight = ww * (1. - t0);
      const auto N = ModTriShape<3>(t0, t1);
      const Tddd X = N[0] * X0 + N[1] * X1 + N[2] * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_gw5 += weight * J_det * inv4pi / nr;
      intGn_gw5 -= weight * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    // --- Dunavant quadrature (direct barycentric coordinates) ---
    // intG = (J_det/2) * sum w * G,  intGn = -(1/2) * sum w * dot(R, cross) * Gn_factor
    const double half_J = J_det * 0.5;

    double intG_d1 = 0., intGn_d1 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_1__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d1 += half_J * w * inv4pi / nr;
      intGn_d1 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d3 = 0., intGn_d3 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_3__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d3 += half_J * w * inv4pi / nr;
      intGn_d3 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d4 = 0., intGn_d4 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_4__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d4 += half_J * w * inv4pi / nr;
      intGn_d4 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d6 = 0., intGn_d6 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_6__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d6 += half_J * w * inv4pi / nr;
      intGn_d6 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d7 = 0., intGn_d7 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_7__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d7 += half_J * w * inv4pi / nr;
      intGn_d7 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d12 = 0., intGn_d12 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_12__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d12 += half_J * w * inv4pi / nr;
      intGn_d12 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d13 = 0., intGn_d13 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_13__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d13 += half_J * w * inv4pi / nr;
      intGn_d13 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d16 = 0., intGn_d16 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_16__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d16 += half_J * w * inv4pi / nr;
      intGn_d16 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d19 = 0., intGn_d19 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_19__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d19 += half_J * w * inv4pi / nr;
      intGn_d19 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    double intG_d25 = 0., intGn_d25 = 0.;
    for (const auto& [l0, l1, l2, w] : __dunavant_25__) {
      const Tddd X = l0 * X0 + l1 * X1 + l2 * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_d25 += half_J * w * inv4pi / nr;
      intGn_d25 -= 0.5 * w * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    // --- Reference: GW11x11 (121 points) ---
    double intG_gw11 = 0., intGn_gw11 = 0.;
    for (const auto& [t0, t1, ww] : __array_GW11xGW11__) {
      const double weight = ww * (1. - t0);
      const auto N = ModTriShape<3>(t0, t1);
      const Tddd X = N[0] * X0 + N[1] * X1 + N[2] * X2;
      const Tddd R = X - target->X;
      const double nr = Norm(R);
      if (nr < 1e-20)
        continue;
      intG_gw11 += weight * J_det * inv4pi / nr;
      intGn_gw11 -= weight * Dot(R, cross) * inv4pi / (nr * nr * nr);
    }

    const char* type = f->Dirichlet ? "D" : "N";
    ofs << face_idx << "," << type << "," << (is_singular ? 1 : 0) << ","
        << r << "," << r / h << ","
        << intG_gw1 << "," << intGn_gw1 << ","
        << intG_gw3 << "," << intGn_gw3 << ","
        << intG_gw5 << "," << intGn_gw5 << ","
        << intG_gw11 << "," << intGn_gw11 << ","
        << intG_d1 << "," << intGn_d1 << ","
        << intG_d3 << "," << intGn_d3 << ","
        << intG_d4 << "," << intGn_d4 << ","
        << intG_d6 << "," << intGn_d6 << ","
        << intG_d7 << "," << intGn_d7 << ","
        << intG_d12 << "," << intGn_d12 << ","
        << intG_d13 << "," << intGn_d13 << ","
        << intG_d16 << "," << intGn_d16 << ","
        << intG_d19 << "," << intGn_d19 << ","
        << intG_d25 << "," << intGn_d25 << ","
        << f->area << "," << f->area / (h * h) << "\n";
    face_idx++;
  }

  ofs.close();
  std::cout << "[exportKernelScatterData] wrote " << face_idx << " faces to "
            << (output_directory / "kernel_scatter.csv").string() << std::endl;
}

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

  /* --------------------------------------------------------------------------*/
  /*                           Set up logging to file
  /* --------------------------------------------------------------------------*/
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
  const std::filesystem::path output_directory = setting.common.output_directory;
  use_linear_element = setting.bem.element.linear;
  use_pseudo_quadratic_element = setting.bem.element.pseudo_quadratic;
  use_true_quadratic_element = setting.bem.element.true_quadratic;
  // Node relocation: method
  {
    using M = SimulationSettings::TimeDomainSettings::NodeRelocationSettings::Method;
    switch (setting.time.node_relocation.method) {
      case M::ALE:            node_relocation_method = NodeRelocationMethod::ALE; break;
      case M::interpolation:  node_relocation_method = NodeRelocationMethod::interpolation; break;
      default:                node_relocation_method = NodeRelocationMethod::none; break;
    }
  }
  // Node relocation: surface precision (auto-resolve from element type unless explicitly set)
  if (setting.time.node_relocation.surface_explicitly_set) {
    using S = SimulationSettings::TimeDomainSettings::NodeRelocationSettings::Surface;
    switch (setting.time.node_relocation.surface) {
      case S::linear:           node_relocation_surface = NodeRelocationSurface::linear; break;
      case S::true_quadratic:   node_relocation_surface = NodeRelocationSurface::true_quadratic; break;
      default:                  node_relocation_surface = NodeRelocationSurface::pseudo_quadratic; break;
    }
  } else {
    if (use_true_quadratic_element)
      node_relocation_surface = NodeRelocationSurface::true_quadratic;
    else
      node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
  }
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

  std::map<std::string, outputInfo> NetOutputInfo = setting.NetOutputInfo;
  const bool use_VPM = setting.vpm.enabled;
  const bool shell_visualization = setting.remeshing.shell_visualization;
  const bool front_advancing_debug = setting.remeshing.front_advancing_debug;
  const auto subsurface_altitude_reject = setting.remeshing.subsurface_altitude_reject;
  const bool initial_mesh_pre_relax = setting.remeshing.initial_mesh_pre_relax.enabled;
  const int initial_mesh_pre_relax_loop = setting.remeshing.initial_mesh_pre_relax.loop;
  const double initial_mesh_pre_relax_coef = setting.remeshing.initial_mesh_pre_relax.coef;
  const std::size_t VPM_wall_min_absorb_receivers = setting.vpm.wall_min_absorb_receivers;
  const double VPM_wall_min_absorb_total_weight = setting.vpm.wall_min_absorb_total_weight;
  const double VPM_sigma_factor = setting.vpm.sigma_factor;
  const std::string VPM_stretching_scheme = setting.vpm.stretching_scheme;
  const std::string VPM_PSE_correction = setting.vpm.PSE_correction;
  const double min_edge_length = setting.remeshing.min_edge_length;
  const int checkpoint_interval = setting.checkpoint.interval;
  const std::string restart_from_checkpoint = setting.checkpoint.restart_from;
  const int checkpoint_max_keep = setting.checkpoint.max_keep;
  std::vector<Network*> FluidObject = setting.FluidObject;
  std::vector<Network*> RigidBodyObject = setting.RigidBodyObject;
  std::vector<Network*> SoftBodyObject = setting.SoftBodyObject;
  std::vector<Network*> AbsorberObject = setting.AbsorberObject;
  std::vector<JSON> MeasurementJSONs = setting.MeasurementJSONs;
  /* --------------------------------------------------------------------------*/

#ifdef USE_TETGEN
  if (tetrahedralize)
    for (auto& network : FluidObject)
      network->tetrahedralize();
#endif

  /* ----------------- Create output directory and copy files -----------------*/

  std::filesystem::create_directories(output_directory);
  std::filesystem::copy_file(setting.settings_file_path, output_directory / "settings.json", std::filesystem::copy_options::overwrite_existing);
  std::filesystem::copy_file("./main_time_domain.cpp", output_directory / "main_time_domain.cpp", std::filesystem::copy_options::overwrite_existing);
  std::filesystem::copy_file("./main.cpp", output_directory / "main.cpp", std::filesystem::copy_options::overwrite_existing);
  //
  std::regex pattern("^BEM.*\\.hpp$");
  for (auto& entry : std::filesystem::directory_iterator("."))
    if (std::regex_match(entry.path().filename().string(), pattern))
      std::filesystem::copy_file(entry.path(), output_directory / entry.path().filename(), std::filesystem::copy_options::overwrite_existing);

  /* --------------------------------------------------------------------------*/

  // auto water = FluidObject[0];
  PVDWriter cornerPointsPVD(output_directory / "cornerPointsPVD.pvd");
  PVDWriter DirichletSurfacePVD(output_directory / "DirichletSurface.pvd");
  PVDWriter vpm_pvd(output_directory / "vpm.pvd");
  Print("setting done");

  /* --------------------------------------------------------------------------*/
  /* --------------------------------------------------------------------------*/
  /* --------------------------------------------------------------------------*/

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
    //  * ------------------------------------------------------ */
    //  *                         メインループ                   */
    //  * ------------------------------------------------------ */
    Buckets<networkFace*> Buckets_Fluid_Faces;
    Buckets<networkPoint*> Buckets_Fluid_Points;

    VortexMethod vpm;
    if (VPM_stretching_scheme == "transpose" || VPM_stretching_scheme == "Transpose" || VPM_stretching_scheme == "TRANSPOSE") {
      vpm.setStretchingScheme(VortexMethod::StretchingScheme::Transpose);
    } else if (VPM_stretching_scheme == "standard" || VPM_stretching_scheme == "Standard" || VPM_stretching_scheme == "STANDARD") {
      vpm.setStretchingScheme(VortexMethod::StretchingScheme::Standard);
    } else {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Unknown VPM_stretching_scheme: " + VPM_stretching_scheme);
    }
    if (VPM_PSE_correction == "none" || VPM_PSE_correction == "None" || VPM_PSE_correction == "NONE" || VPM_PSE_correction == "0") {
      vpm.setPSECorrectionMode(VortexMethod::PSECorrectionMode::None);
    } else if (VPM_PSE_correction == "grad" || VPM_PSE_correction == "Grad" || VPM_PSE_correction == "GRAD" || VPM_PSE_correction == "gradient" || VPM_PSE_correction == "Gradient" || VPM_PSE_correction == "GRADIENT" || VPM_PSE_correction == "1") {
      vpm.setPSECorrectionMode(VortexMethod::PSECorrectionMode::Gradient);
    } else if (VPM_PSE_correction == "curvature" || VPM_PSE_correction == "Curvature" || VPM_PSE_correction == "CURVATURE" || VPM_PSE_correction == "laplacian" || VPM_PSE_correction == "Laplacian" || VPM_PSE_correction == "LAPLACIAN" || VPM_PSE_correction == "2") {
      vpm.setPSECorrectionMode(VortexMethod::PSECorrectionMode::Curvature);
    } else {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Unknown VPM_PSE_correction: " + VPM_PSE_correction);
    }
    double time_setup = 0.0, time_solve = 0.0, unknownsizse = 0.0;

    // --- Checkpoint restart ---
    int start_time_step = 0;
    double last_successful_dt = max_dt;
    if (!restart_from_checkpoint.empty()) {
      auto [rstep, rtime, rdt] =
          BEM_Checkpoint::readCheckpoint(restart_from_checkpoint, FluidObject, RigidBodyObject, SoftBodyObject,
                                         vpm, use_true_quadratic_element);
      start_time_step = rstep + 1; // resume from the NEXT step
      simulation_time = rtime;
      last_successful_dt = rdt;
      std::cout << Green << "[restart] Loaded checkpoint: step=" << rstep
                << ", time=" << rtime << ", dt=" << rdt << colorReset << std::endl;

      // Output initial state from checkpoint so it appears in ParaView/JSON output
      {
        OutputContext ctx{.dt = rdt, .time_step = rstep, .simulation_time = rtime,
                         .output_directory = output_directory,
                         .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};
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
      int RK_order = 4;

      double spacing = 0.;
      for (auto& water : FluidObject)
        spacing += Mean(extLength(water->getLines())) * 10;
      spacing /= FluidObject.size();
      std::vector<Network*> AllObjects = Join(FluidObject, RigidBodyObject, SoftBodyObject);

      // --------------------------------------------------------------------------
      // Initial mesh pre-relaxation (relocation-style smoothing), before the first remesh/collapse.
      // IMPORTANT: boundary types (Dirichlet/Neumann/CORNER) must be decided at this point because
      // the mesh smoothing (calculateVecToSurface) applies different constraints depending on them.
      // --------------------------------------------------------------------------
      if (time_step == 0 && start_time_step == 0 && initial_mesh_pre_relax && initial_mesh_pre_relax_loop > 0 && initial_mesh_pre_relax_coef > 0.0) {
        TimeWatch pre_watch;
        constexpr int pre_relax_cycles = 1; // flip -> smoothing -> flip -> ... (at t=0 only)
        std::cout << Green << "[initial_mesh_pre_relax] begin (loop=" << initial_mesh_pre_relax_loop << ", coef=" << initial_mesh_pre_relax_coef << ", cycles=" << pre_relax_cycles << ")" << colorReset << std::endl;

        auto write_pre_relax_mesh = [&](Network* net, const std::string& tag) {
          const std::string suffix = tag.empty() ? "" : ("_" + tag);

          // VTU
          if (auto it = NetOutputInfo.find(net->getName()); it != NetOutputInfo.end()) {
            const auto& info = it->second;
            std::filesystem::path filename = info.vtu_file_name + std::string("pre_relax") + suffix + ".vtu";
            mk_vtu(output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, /*dt=*/0.0));
            if (info.PVD) {
              info.PVD->push(filename, simulation_time);
              info.PVD->output();
            }
          } else {
            std::filesystem::path filename = net->getName() + std::string("_pre_relax") + suffix + ".vtu";
            mk_vtu(output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, /*dt=*/0.0));
          }

          // OBJ (debug)
          {
            std::filesystem::path filename = net->getName() + std::string("_pre_relax") + suffix + ".obj";
            std::ofstream ofs(output_directory / filename);
            if (ofs)
              createOBJ(ofs, *net);
          }
        };

        // Prepare spatial acceleration for contact queries used by setBoundaryTypes / getNearestContactFace.
        _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);

        const std::vector<Network*> contact_objects = Join(RigidBodyObject, SoftBodyObject);

        for (auto& water : FluidObject) {
          // Use a dummy RK state (Heun/RK not needed): we just want RK_with_Ubuff/RK_without_Ubuff helpers.
          // Keep steps=1 so getX() corresponds to a full-step update (not dt/2).
          constexpr double relax_dt = 1.0;
          const double rad = std::acos(-1.0) / 180.0;
          const Tdd flip_limitD = {20.0 * rad /*target n diff*/, 20.0 * rad /*acceptable normal change*/};
          const Tdd flip_limitN = {5.0 * rad /*target n diff*/, 5.0 * rad /*acceptable normal change*/}; // Stricter for Neumann to prevent large deformations

          auto init_relax_rk = [&]() {
            for (const auto& p : water->getPoints()) {
              p->u_node.fill(0.0);
              p->RK_X.initialize(relax_dt, simulation_time, ToX(p), 1);
            }
          };

          // Alternate: flip -> smoothing -> flip -> ... to improve mesh quality before the first collapse.
          setBoundaryTypes(water, contact_objects);
          if (use_true_quadratic_element)
            computeAllCornerMidpointOffsets(water);
          // Pre-relaxation always uses pseudo_quadratic surface regardless of element type.
          // At this point, X_mid values are at straight-edge midpoints (no curvature info),
          // so true_quadratic surface would produce flat surfaces. pseudo_quadratic uses
          // neighbor reconstruction for better curvature approximation.
          auto saved_surface = node_relocation_surface;
          if (node_relocation_surface == NodeRelocationSurface::true_quadratic)
            node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
          for (int c = 0; c < pre_relax_cycles; ++c) {
            init_relax_rk();
            water->setGeometricPropertiesForce();
            flipIfBatched(*water, flip_limitD, flip_limitN, "pre-relax");
            calculateVecToSurface(*water, initial_mesh_pre_relax_loop, initial_mesh_pre_relax_coef);
            for (const auto& p : water->getPoints())
              p->setXSingle(RK_with_Ubuff(p));
            std::string cycle_tag = "cycle" + std::to_string(c + 1);
          }
          node_relocation_surface = saved_surface;
          water->setGeometricPropertiesForce();
          write_pre_relax_mesh(water, "after_relax");
        }

        std::cout << Green << "[initial_mesh_pre_relax] done" << Blue << "\nElapsed time: " << Red << pre_watch() << colorReset << " s\n";
      }

      if (time_step == 0 && start_time_step == 0 && shell_visualization) {
        OutputContext ctx{.dt = 0.0, .time_step = 0, .simulation_time = simulation_time, .output_directory = output_directory, .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};
        OutputParaView::write_shell_step(ctx, NetOutputInfo, FluidObject);
      }
      if (time_step == 0 && start_time_step == 0 && front_advancing_debug) {
        OutputContext ctx{.dt = 0.0, .time_step = 0, .simulation_time = simulation_time, .output_directory = output_directory, .cpu_clock_start = cpu_clock_start, .wall_clock_start = wall_clock_start};
        OutputParaView::write_aft_candidates_step(ctx, NetOutputInfo, FluidObject);
      }

      // Variables that must survive across retry attempts
      std::vector<double> convergence;
      std::vector<double> ElapsedTimeSetup, ElapsedTimeSolve;
      double ElapsedTimeNodeRelocation = 0, ElapsedTimeTotal = 0;
      std::unordered_map<std::string, std::unordered_map<int, int>> repeated_reject_faces;

      auto parse_int_after = [](const std::string& text, const std::string& key) -> std::optional<int> {
        const auto pos = text.find(key);
        if (pos == std::string::npos)
          return std::nullopt;
        std::size_t begin = pos + key.size();
        std::size_t end = begin;
        while (end < text.size() && std::isdigit(static_cast<unsigned char>(text[end])))
          ++end;
        if (end == begin)
          return std::nullopt;
        return std::stoi(text.substr(begin, end - begin));
      };

      auto parse_water_name_from_reason = [](const std::string& reason) -> std::optional<std::string> {
        const auto on_pos = reason.find(" on ");
        if (on_pos == std::string::npos)
          return std::nullopt;
        const auto at_pos = reason.find(" at time_step", on_pos + 4);
        if (at_pos == std::string::npos || at_pos <= on_pos + 4)
          return std::nullopt;
        return reason.substr(on_pos + 4, at_pos - (on_pos + 4));
      };

      auto remember_rejected_face = [&](const step_failure& e) {
        const bool is_face_collapse_failure =
            e.reason.find("tiny face area ratio") != std::string::npos ||
            e.reason.find("subsurface face altitude ratio") != std::string::npos;
        if (!is_face_collapse_failure)
          return;
        const auto face_index = parse_int_after(e.reason, "face_index=");
        const auto water_name = parse_water_name_from_reason(e.reason);
        if (!face_index || !water_name)
          return;
        auto& count = repeated_reject_faces[*water_name][*face_index];
        ++count;
        std::cout << Yellow << "[step_reject] remember pathological face: water=" << *water_name
                  << " face=" << *face_index << " count=" << count
                  << colorReset << std::endl;
      };

      auto collapse_repeatedly_rejected_faces = [&](Network& water, const int step_retry) {
        if (step_retry <= 0)
          return;
        auto it = repeated_reject_faces.find(water.getName());
        if (it == repeated_reject_faces.end())
          return;
        bool changed = false;
        for (const auto& [face_index, count] : it->second) {
          if (count < 2)
            continue;
          if (collapseFaceByIndexIfPossible(water, face_index)) {
            changed = true;
            std::cout << Yellow << "[step_reject] collapsed repeatedly rejected face on retry: water="
                      << water.getName() << " face=" << face_index
                      << " count=" << count << " retry=" << step_retry
                      << colorReset << std::endl;
          }
        }
        if (changed) {
          water.setGeometricPropertiesForce();
          water.checkConnectivity();
        }
      };

      auto sync_body_states_from_rk = [&]() {
        for (auto* net : RigidBodyObject) {
          if (net->RK_COM.steps > 0)
            net->COM = net->RK_COM.getX();
          if (net->RK_Q.steps > 0)
            net->Q = Normalize(net->RK_Q.getX());
          if (net->RK_Velocity.steps > 0)
            net->velocity = net->RK_Velocity.getX();
          net->setGeometricPropertiesForce();
        }

        for (auto* net : SoftBodyObject) {
          if (net->RK_COM.steps > 0)
            net->COM = net->RK_COM.getX();
          if (net->RK_Q.steps > 0)
            net->Q = Normalize(net->RK_Q.getX());
          if (net->RK_Velocity.steps > 0)
            net->velocity = net->RK_Velocity.getX();
          net->setGeometricPropertiesForce();
        }
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
        sync_body_states_from_rk();
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
        if (use_VPM)
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
        sync_body_states_from_rk();
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

        if (use_VPM) {
          std::string filename = "vpm_" + tag + ".vtp";
          std::filesystem::path path = output_directory / filename;
          std::ofstream ofs(path);
          if (ofs) {
            const auto& particles = vpm.getParticles();
            const size_t num_particles = particles.size();
            ofs << "<?xml version=\"1.0\"?>\n";
            ofs << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
            ofs << "  <PolyData>\n";
            ofs << "    <Piece NumberOfPoints=\"" << num_particles << "\" NumberOfVerts=\"" << num_particles << "\">\n";
            ofs << "      <Points>\n";
            ofs << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (const auto& p : particles)
              ofs << "          " << p.x[0] << " " << p.x[1] << " " << p.x[2] << "\n";
            ofs << "        </DataArray>\n";
            ofs << "      </Points>\n";
            ofs << "      <Verts>\n";
            ofs << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
            for (size_t i = 0; i < num_particles; ++i)
              ofs << "          " << i << "\n";
            ofs << "        </DataArray>\n";
            ofs << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
            for (size_t i = 0; i < num_particles; ++i)
              ofs << "          " << i + 1 << "\n";
            ofs << "        </DataArray>\n";
            ofs << "      </Verts>\n";
            ofs << "      <PointData>\n";
            ofs << "        <DataArray type=\"Float64\" Name=\"Vorticity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (const auto& p : particles)
              ofs << "          " << p.alpha[0] << " " << p.alpha[1] << " " << p.alpha[2] << "\n";
            ofs << "        </DataArray>\n";
            ofs << "      </PointData>\n";
            ofs << "    </Piece>\n";
            ofs << "  </PolyData>\n";
            ofs << "</VTKFile>\n";
          }
        }

        const auto checkpoint_path = output_directory / ("checkpoint_" + tag + ".bin");
        BEM_Checkpoint::writeCheckpointToPath(checkpoint_path, time_step, snapshot_time, dt,
                                             FluidObject, RigidBodyObject, SoftBodyObject, vpm,
                                             use_true_quadratic_element);
      };

      // === Step rejection: save state before remesh, retry with smaller dt on failure ===
      constexpr int max_step_retries = 5;
      constexpr double dt_reduction_factor = 0.5;
      double dt_override = 0.0; // non-zero means use this dt instead of CFL
      auto step_checkpoint_path = output_directory / "_step_retry.bin";

      // Save state before this step (for rollback on failure)
      {
        sync_body_states_from_rk();
        // Write checkpoint with current step's state using the standard function
        // but to a temporary filename via the step number 999999
        BEM_Checkpoint::writeCheckpoint(output_directory, 999999,
                                        simulation_time, last_successful_dt,
                                        FluidObject, RigidBodyObject, SoftBodyObject, vpm,
                                        use_true_quadratic_element);
        auto temp_path = output_directory / "checkpoint_999999.bin";
        if (std::filesystem::exists(temp_path))
          std::filesystem::rename(temp_path, step_checkpoint_path);
      }

      int current_rk_step = 0;
      for (int step_retry = 0; step_retry <= max_step_retries; ++step_retry) {
        if (step_retry > 0) {
          // Rollback: reload mesh state from step checkpoint
          std::cout << Red << "[step_reject] retry " << step_retry << "/" << max_step_retries
                    << " with dt_override=" << dt_override << colorReset << std::endl;
          auto [rstep, rtime, rdt] =
              BEM_Checkpoint::readCheckpoint(step_checkpoint_path, FluidObject, RigidBodyObject, SoftBodyObject,
                                             vpm, use_true_quadratic_element);
          simulation_time = rtime;
          // Re-initialize BVP after mesh reconstruction (unique_ptr prevents copy/move)
          BVP.WATERS = FluidObject;
          BVP.ilu_preconditioner.reset();
          BVP.ilut_preconditioner.reset();
          BVP.schwarz_preconditioner.reset();
          BVP.ilu_topology_hash_cached = 0; // force rebuild
          BVP.output_directory = output_directory;
          // Recompute AllObjects after mesh reconstruction
          AllObjects = Join(FluidObject, RigidBodyObject, SoftBodyObject);
        }

        try {
      _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);
      for (auto& water : FluidObject) {
        water->setContactFaces(Join(RigidBodyObject, SoftBodyObject));
      }
      for (auto& water : FluidObject) {
        remesh_for_main_loop(*water, time_step, min_edge_length, tetrahedralize, surface_flip, setting.remeshing.collision,
                             setting.remeshing.surface_split, setting.remeshing.surface_collapse);
        collapse_repeatedly_rejected_faces(*water, step_retry);
        monitorTinyFaceRelativeToLocalMean(*water, time_step, std::nullopt);
        throwIfSubsurfaceFaceAltitudeTooSmall(*water, time_step, std::nullopt, subsurface_altitude_reject);
        throwIfTinyFaceRelativeToLocalMean(*water, time_step, std::nullopt);
      }
      throwIfStructurePenetrated(FluidObject, Join(RigidBodyObject, SoftBodyObject), time_step, "post-remesh");

      // --------------------------------------------------------------------------
      // Initial condition: apply known analytical solution at t=0
      // Move free surface nodes to η(x) and set φ from theory.
      // --------------------------------------------------------------------------
      if (time_step == 0 && start_time_step == 0) {
        bool has_any_ic = std::ranges::any_of(FluidObject, [](const Network* w) { return w->ic_eta && w->ic_phi; });
        if (has_any_ic) {
          TimeWatch ic_watch;

          // Ensure boundary types are set (needed for Dirichlet/Neumann classification).
          // Must always run here because remesh_for_main_loop may have invalidated
          // boundary types even if initial_mesh_pre_relax already called setBoundaryTypes.
          _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);
          for (auto& water : FluidObject) {
            water->setContactFaces(Join(RigidBodyObject, SoftBodyObject));
            setBoundaryTypes(water, Join(RigidBodyObject, SoftBodyObject));
            if (use_true_quadratic_element)
              computeAllCornerMidpointOffsets(water);
          }

          for (auto& water : FluidObject) {
            if (!water->ic_eta || !water->ic_phi)
              continue;

            std::cout << Green << "[initial_condition] applying to " << water->getName() << colorReset << std::endl;

            // Step 1: Move Dirichlet (free surface) nodes to z = η(x, y, t=0)
            for (const auto& p : water->getPoints()) {
              if (p->Dirichlet || p->CORNER) {
                auto X = ToX(p);
                double eta = water->ic_eta(X, 0.0);
                p->setXSingle({X[0], X[1], eta});
              }
            }
            water->setGeometricPropertiesForce();

            // Step 2: Mesh smoothing to improve mesh quality after deformation
            {
              // Use pseudo_quadratic surface for IC smoothing (same reason as pre-relaxation:
              // X_mid has no curvature info at this point).
              auto saved_surface = node_relocation_surface;
              if (node_relocation_surface == NodeRelocationSurface::true_quadratic)
                node_relocation_surface = NodeRelocationSurface::pseudo_quadratic;
              constexpr double relax_dt = 1.0;
              for (const auto& p : water->getPoints()) {
                p->u_node.fill(0.0);
                p->RK_X.initialize(relax_dt, 0.0, ToX(p), 1);
              }
              calculateVecToSurface(*water, 30, 0.05);
              for (const auto& p : water->getPoints()) {
                p->setXSingle(RK_with_Ubuff(p));
              }
              water->setGeometricPropertiesForce();
              node_relocation_surface = saved_surface;
            }

            // Step 3: Set φ from theory at the final (post-smoothing) node positions
            {
              int n_dirichlet = 0, n_corner = 0, n_neumann = 0, n_other = 0;
              double phi_min = 1e30, phi_max = -1e30;
              for (const auto& p : water->getPoints()) {
                if (p->Dirichlet)
                  ++n_dirichlet;
                else if (p->CORNER)
                  ++n_corner;
                else if (p->Neumann)
                  ++n_neumann;
                else
                  ++n_other;

                if (p->Dirichlet || p->CORNER) {
                  auto X = ToX(p);
                  double phi_val = water->ic_phi(X, 0.0);
                  std::get<0>(p->phiphin) = phi_val;
                  p->phi_Dirichlet = phi_val;
                  phi_min = std::min(phi_min, phi_val);
                  phi_max = std::max(phi_max, phi_val);
                }
              }
              std::cout << Green << "[initial_condition] node classification: "
                        << "Dirichlet=" << n_dirichlet << " CORNER=" << n_corner
                        << " Neumann=" << n_neumann << " other=" << n_other
                        << "\n  phi applied to " << (n_dirichlet + n_corner) << " / "
                        << (n_dirichlet + n_corner + n_neumann + n_other) << " nodes"
                        << "\n  phi range: [" << phi_min << ", " << phi_max << "]"
                        << colorReset << std::endl;
            }

            // Step 4: Set phi_mid and X_mid from theory (true quadratic edge midpoints)
            if (use_true_quadratic_element) {
              int n_mid_set = 0;
              double phi_mid_min = 1e30, phi_mid_max = -1e30;
              double max_dev_from_avg = 0;
              for (auto* l : water->getBoundaryLines()) {
                auto [pA, pB] = l->getPoints();
                // Update X_mid to midpoint of deformed mesh
                l->X_mid = 0.5 * (pA->X + pB->X);
                if (l->CORNER)
                  l->X_mid += l->corner_midpoint_offset;
                l->Xtarget = l->X_mid;

                if (l->Dirichlet || l->CORNER) {
                  double phi_val = water->ic_phi(l->X_mid, 0.0);
                  double phi_avg = 0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
                  l->phiphin[0] = phi_val;
                  phi_mid_min = std::min(phi_mid_min, phi_val);
                  phi_mid_max = std::max(phi_mid_max, phi_val);
                  max_dev_from_avg = std::max(max_dev_from_avg, std::abs(phi_val - phi_avg));
                  ++n_mid_set;
                }
              }
              std::cout << Green << "[initial_condition] phi_mid set from theory: n=" << n_mid_set
                        << " range=[" << phi_mid_min << ", " << phi_mid_max << "]"
                        << " max|phi_mid - avg(phiA,phiB)|=" << max_dev_from_avg
                        << colorReset << std::endl;
            }

            std::cout << Green << "[initial_condition] done for " << water->getName()
                      << Blue << "\nElapsed time: " << Red << ic_watch() << colorReset << " s\n";
          }
        }
      }

      //
      CoordinateBounds bounds_org(FluidObject[0]->bounds);
      for (auto& water : FluidObject)
        bounds_org += water->bounds;

      CoordinateBounds bounds = bounds_org.scaledBounds(1.5);

      Buckets_Fluid_Faces.initialize(bounds, bounds.getScale() / 10.);
      Buckets_Fluid_Points.initialize(bounds, bounds.getScale() / 10.);

      // # ------------------------------------------------------ */
      // #                       刻み時間の決定                   */
      // # ------------------------------------------------------ */
      for (auto& water : FluidObject) {
        show_info(*water);

        auto dt_cfl = dt_CFL(*water, max_dt, .5);
        if (dt > dt_cfl)
          dt = dt_cfl;
        if (time_step <= 2)
          dt = dt / 10.;
      }

      if (use_VPM) {
        const auto vpm_dt = vpm.suggestTimeStep(dt /*max_dt (do not increase dt)*/);
        if (vpm_dt.dt < dt) {
          Print("[dt] limited by VPM:", "dt=", vpm_dt.dt, "dt_strain=", vpm_dt.dt_strain, "dt_diffusion=", vpm_dt.dt_diffusion, "dt_move=", vpm_dt.dt_move, "max_strain_rate=", vpm_dt.max_strain_rate, "max_u=", vpm_dt.max_u, "min_sigma=", vpm_dt.min_sigma, "particles=", vpm_dt.particle_count);
        }
        dt = std::min(dt, vpm_dt.dt);
      }
      if (dt < 1E-13)
        dt = 1E-13;
      const double dt_growth_limit = 1.25;
      if (last_successful_dt > 0.0) {
        const double dt_capped = last_successful_dt * dt_growth_limit;
        if (dt > dt_capped) {
          std::cout << Yellow << "[dt] growth limited: " << dt << " -> " << dt_capped
                    << " (last_successful_dt=" << last_successful_dt << ")" << colorReset << std::endl;
          dt = dt_capped;
        }
      }
      // Apply dt override from step rejection (smaller dt for retry)
      if (dt_override > 0 && dt_override < dt) {
        std::cout << Yellow << "[step_reject] dt reduced: " << dt << " -> " << dt_override << colorReset << std::endl;
        dt = dt_override;
      }
      Print("===========================================================================");
      Print("       dt :", Red, std::setprecision(10), dt, colorReset);
      Print("time_step :", Red, time_step, colorReset);
      Print("real time :", Red, simulation_time, colorReset);
      Print("---------------------------------------------------------------------------");

      // ---------------------------------------------------------------------------
      // Helper: robust inside test for closed surfaces.
      // `Network::InsideQ()` depends on the orientation of boundary faces (winding number sign).
      // Use |windingNumber| so that globally inverted meshes are still handled correctly.
      // ---------------------------------------------------------------------------
      auto inside_closed_surface = [&](const Network* net, const Tddd& X) -> bool {
        if (!net)
          return false;
        if (!net->CoordinateBounds::InsideQ(X))
          return false;
        return std::abs(net->windingNumber(X)) >= 0.5;
      };

      auto push_out_of_body = [&](const Network* rb, VortexParticle& p, double delta) {
        if (!rb)
          return;
        auto [f, nearest_pos] = rb->Nearest(p.x);
        if (!f)
          return;
        Tddd n = f->normal;
        const double nn = Norm(n);
        if (!(nn > 0.0))
          return;
        n = n / nn;

        // project to the plane of the nearest face, then offset to the outside side.
        const double d = Dot(p.x - nearest_pos, n);
        const Tddd x_plane = p.x - d * n;

        // choose the side which is outside (robust against normal orientation).
        Tddd x_try = x_plane + delta * n;
        if (inside_closed_surface(rb, x_try))
          x_try = x_plane - delta * n;

        // if still inside (e.g., highly curved / concave), try a few larger offsets.
        if (inside_closed_surface(rb, x_try)) {
          double dd = delta;
          for (int k = 0; k < 4 && inside_closed_surface(rb, x_try); ++k) {
            dd *= 2.0;
            x_try = x_plane + dd * n;
            if (inside_closed_surface(rb, x_try))
              x_try = x_plane - dd * n;
          }
        }

        p.x = x_try;
      };

      //% ---------------------------------------------------------------------------
      //% VPM operator splitting (Strang): diffusion half-step at the beginning.
      //% ---------------------------------------------------------------------------
      if (use_VPM) {
        vpm.applyDiffusionEuler(0.5 * dt);

        // Wall boundary handling for particles (push-out / reflection).
        vpm.processParticles([&](VortexParticle& p) {
          for (const auto& rb : RigidBodyObject) {
            if (!inside_closed_surface(rb, p.x))
              continue;
            const double delta = std::max(0.3 * p.sigma, 1e-12);
            push_out_of_body(rb, p, delta);
          }
        });
      }

      // b@ ------------------------------------------------------------- */
      // b@           初期値問題を解く（時間微分方程式を数値積分する）    */
      // b@ ------------------------------------------------------------- */

      for (auto& water : FluidObject) {
        for (const auto& p : water->getPoints()) {
          p->RK_phi.initialize(dt, simulation_time, std::get<0>(p->phiphin), RK_order);
          p->RK_X.initialize(dt, simulation_time, getPosition(p), RK_order);
        }
      }

      // Initialize midpoint RK4 (true quadratic)
      // Note: initialize ALL boundary lines, not just Dirichlet ones,
      // because setBoundaryTypes() (called later inside the RK loop)
      // may change which lines are Dirichlet.
      if (use_true_quadratic_element) {
        for (auto& water : FluidObject) {
          for (auto* l : water->getBoundaryLines()) {
            // Keep evolved midpoint positions across time steps.
            // Only initialize from endpoints for newly created lines.
            if (l->RK_X.steps == 0) {
              auto [pA, pB] = l->getPoints();
              l->X_mid = 0.5 * (pA->X + pB->X);
              if (l->CORNER)
                l->X_mid += l->corner_midpoint_offset;
            }
            l->Xtarget = l->X_mid;
            l->RK_phi.initialize(dt, simulation_time, l->phiphin[0], RK_order);
            l->RK_X.initialize(dt, simulation_time, l->X_mid, RK_order);
          }
        }
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
      if (use_VPM)
        vpm.initializeRK(dt, simulation_time, RK_order);

      // b@ -------------------------------------------------------------- */

      auto rebuild_fluid_global_bounds = [&]() {
        CoordinateBounds bounds_org(FluidObject[0]->bounds);
        for (auto& water : FluidObject)
          bounds_org += water->bounds;
        return CoordinateBounds(bounds_org.scaledBounds(1.5));
      };

      auto rebuild_fluid_global_point_bucket = [&]() {
        const auto bounds = rebuild_fluid_global_bounds();
        Buckets_Fluid_Points.initialize(bounds, bounds.getScale() / 10.);
        for (auto& water : FluidObject)
          for (const auto& p : water->getPoints())
            Buckets_Fluid_Points.add(ToX(p), p);
      };

      auto rebuild_fluid_global_buckets = [&]() {
        const auto bounds = rebuild_fluid_global_bounds();
        Buckets_Fluid_Faces.initialize(bounds, bounds.getScale() / 10.);
        Buckets_Fluid_Points.initialize(bounds, bounds.getScale() / 10.);
        for (auto& water : FluidObject) {
          for (const auto& p : water->getPoints())
            Buckets_Fluid_Points.add(ToX(p), p);
          for (const auto& f : water->getBoundaryFaces())
            if (f)
              Buckets_Fluid_Faces.add(f->getXtuple(), f);
        }
      };

      rebuild_fluid_global_point_bucket();

      // b@ -------------------------------------------------------------- */

      int RK_step = 0;

      TimeWatch watch;

      auto log_corner_connected_neumann_lines_after_boundary_types = [&](const Network* water, const char* phase) {
        if (!water)
          return;
        auto local_line_length = [](const networkLine* const l) {
          static std::unordered_set<networkPoint*> adjacent_points;
          adjacent_points.clear();
          for (const auto& f : l->getBoundaryFaces()) {
            auto points = f->getPoints();
            adjacent_points.insert(points.begin(), points.end());
          }
          static std::unordered_set<networkLine*> adjacent_lines;
          adjacent_lines.clear();
          for (const auto& p : adjacent_points)
            for (const auto& L : p->getBoundaryLines())
              if (L != l)
                adjacent_lines.insert(L);
          if (adjacent_lines.empty())
            return 0.0;
          auto ret = 0.;
          for (const auto& L : adjacent_lines)
            ret += L->length();
          return ret / adjacent_lines.size();
        };
        auto is_corner_connected_neumann_line = [](const networkLine* l) {
          if (!l || !l->Neumann || l->CORNER)
            return false;
          auto [p0, p1] = l->getPoints();
          return p0 && p1 && (p0->CORNER || p1->CORNER);
        };
        std::vector<networkLine*> candidates;
        for (auto* l : water->getBoundaryLines())
          if (is_corner_connected_neumann_line(l))
            candidates.emplace_back(l);
        std::ranges::sort(candidates, [](const auto* a, const auto* b) { return a->length() < b->length(); });
        if (candidates.size() > 5)
          candidates.resize(5);

        for (auto* l : candidates) {
          auto [p0, p1] = l->getPoints();
          const auto faces = l->getBoundaryFaces();
          const double len = l->length();
          const double local_mean_len = local_line_length(l);
          double alt0 = -1.0, alt1 = -1.0, alt_threshold = -1.0;
          double alt_ratio0 = -1.0, alt_ratio1 = -1.0;
          double aspect_ratio0 = -1.0, aspect_ratio1 = -1.0;
          double min_angle_deg0 = -1.0, min_angle_deg1 = -1.0;
          double area0 = -1.0, area1 = -1.0;
          double mean_area0 = -1.0, mean_area1 = -1.0;
          int common_points = -1;
          int opp0_lines = -1, opp1_lines = -1;
          bool face0_dirichlet = false, face0_neumann = false;
          bool face1_dirichlet = false, face1_neumann = false;
          if (faces.size() == 2 && faces[0] && faces[1] && p0 && p1) {
            auto* f0 = faces[0];
            auto* f1 = faces[1];
            face0_dirichlet = f0->Dirichlet;
            face0_neumann = f0->Neumann;
            face1_dirichlet = f1->Dirichlet;
            face1_neumann = f1->Neumann;
            auto face_shape_metrics = [](const networkFace* f) {
              std::array<double, 2> out{-1.0, -1.0};
              if (!f)
                return out;
              const auto pts = f->getPoints();
              const Tddd a = pts[0]->X, b = pts[1]->X, c = pts[2]->X;
              const double l01 = Norm(a - b);
              const double l12 = Norm(b - c);
              const double l20 = Norm(c - a);
              const double max_edge = std::max({l01, l12, l20});
              const double area = boundaryFaceArea(const_cast<networkFace*>(f));
              const double altitude = (max_edge > 1e-20) ? 2.0 * area / max_edge : 0.0;
              out[0] = (altitude > 1e-20) ? max_edge / altitude : 1E+100;
              auto clamp_cos = [](double x) { return std::max(-1.0, std::min(1.0, x)); };
              const double ang0 = std::acos(clamp_cos(Dot((b - a) / l01, (c - a) / l20))) * 180.0 / M_PI;
              const double ang1 = std::acos(clamp_cos(Dot((a - b) / l01, (c - b) / l12))) * 180.0 / M_PI;
              const double ang2 = 180.0 - ang0 - ang1;
              out[1] = std::min({ang0, ang1, ang2});
              return out;
            };
            auto [a, this0, b, l1, p2, l2] = f0->getPointsAndLines(l);
            auto [q0, this1, q1, e1, q2, e2] = f1->getPointsAndLines(l);
            if (this0 == l && this1 == l && a == q1 && b == q0) {
              area0 = boundaryFaceArea(f0);
              area1 = boundaryFaceArea(f1);
              const double max_edge0 = std::max({Norm(f0->getPoints()[0]->X - f0->getPoints()[1]->X),
                                                 Norm(f0->getPoints()[1]->X - f0->getPoints()[2]->X),
                                                 Norm(f0->getPoints()[2]->X - f0->getPoints()[0]->X)});
              const double max_edge1 = std::max({Norm(f1->getPoints()[0]->X - f1->getPoints()[1]->X),
                                                 Norm(f1->getPoints()[1]->X - f1->getPoints()[2]->X),
                                                 Norm(f1->getPoints()[2]->X - f1->getPoints()[0]->X)});
              alt0 = (max_edge0 > 1e-20) ? 2.0 * area0 / max_edge0 : 0.0;
              alt1 = (max_edge1 > 1e-20) ? 2.0 * area1 / max_edge1 : 0.0;
              alt_threshold = (local_mean_len > 0.0) ? 0.1 * local_mean_len : -1.0;
              alt_ratio0 = (local_mean_len > 0.0) ? alt0 / local_mean_len : -1.0;
              alt_ratio1 = (local_mean_len > 0.0) ? alt1 / local_mean_len : -1.0;
              auto shape0 = face_shape_metrics(f0);
              auto shape1 = face_shape_metrics(f1);
              aspect_ratio0 = shape0[0];
              aspect_ratio1 = shape1[0];
              min_angle_deg0 = shape0[1];
              min_angle_deg1 = shape1[1];
              mean_area0 = localMeanFaceArea(f0);
              mean_area1 = localMeanFaceArea(f1);
              common_points = static_cast<int>(Intersection(p0->getNeighborPointsOnSurfaces(),
                                                            p1->getNeighborPointsOnSurfaces()).size());
              opp0_lines = p2 ? static_cast<int>(p2->getLines().size()) : -1;
              opp1_lines = q2 ? static_cast<int>(q2->getLines().size()) : -1;
            }
          }
          std::cout << Magenta << "[corner_neumann_debug] " << phase
                    << " p0=" << (p0 ? p0->index : -1)
                    << " p1=" << (p1 ? p1->index : -1)
                    << " len=" << len
                    << " local_mean_len=" << local_mean_len
                    << " faces=" << faces.size()
                    << " alt0=" << alt0
                    << " alt1=" << alt1
                    << " alt_threshold=" << alt_threshold
                    << " alt_ratio0=" << alt_ratio0
                    << " alt_ratio1=" << alt_ratio1
                    << " aspect_ratio0=" << aspect_ratio0
                    << " aspect_ratio1=" << aspect_ratio1
                    << " min_angle_deg0=" << min_angle_deg0
                    << " min_angle_deg1=" << min_angle_deg1
                    << " area0=" << area0
                    << " area1=" << area1
                    << " local_mean_area0=" << mean_area0
                    << " local_mean_area1=" << mean_area1
                    << " common_points=" << common_points
                    << " line_flags={D:" << l->Dirichlet << ",N:" << l->Neumann << ",C:" << l->CORNER << "}"
                    << " face0_flags={D:" << face0_dirichlet << ",N:" << face0_neumann << "}"
                    << " face1_flags={D:" << face1_dirichlet << ",N:" << face1_neumann << "}"
                    << " endpoint0_flags={D:" << (p0 ? p0->Dirichlet : false)
                    << ",N:" << (p0 ? p0->Neumann : false)
                    << ",C:" << (p0 ? p0->CORNER : false) << "}"
                    << " endpoint1_flags={D:" << (p1 ? p1->Dirichlet : false)
                    << ",N:" << (p1 ? p1->Neumann : false)
                    << ",C:" << (p1 ? p1->CORNER : false) << "}"
                    << " endpoint_corner={" << (p0 ? p0->CORNER : false) << "," << (p1 ? p1->CORNER : false) << "}"
                    << " endpoint_lines={" << (p0 ? p0->getLines().size() : 0) << "," << (p1 ? p1->getLines().size() : 0) << "}"
                    << " opposite_lines={" << opp0_lines << "," << opp1_lines << "}"
                    << " x0=" << (p0 ? p0->X : Tddd{0., 0., 0.})
                    << " x1=" << (p1 ? p1->X : Tddd{0., 0., 0.})
                    << colorReset << std::endl;
        }
      };

      /* -------------------------------------------------------------------------- */

      do {
        auto RK_time = (*(Buckets_Fluid_Points.data1D).begin())->RK_X.gett(); //%各ルンゲクッタの時刻を使う
        current_rk_step = ++RK_step;
        std::cout << "RK_step = " << RK_step << "/" << RK_order << ", RK_time = " << RK_time << ", simulation_time = " << simulation_time << std::endl;

        if (RK_step == 1) {
          _Pragma("omp parallel for") for (const auto& net : AllObjects) net->makeBuckets(net->getScale() / 10.);

          for (auto& water : FluidObject)
            water->setContactFaces(Join(RigidBodyObject, SoftBodyObject));

          /*
        RKで節点が大きく移動する可能性がるので，setBoundaryTypesをここに移動した．
        RKの間，flipはされない．
        setBoundaryTypesは，接触面を更新するためにここに移動した．
        ただ，今は大きく動かないとしてはじめだけにしている．
        */
          std::cout << Green << "makeBuckets" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
          //! 体積を保存するようにリメッシュする必要があるだろう．
          for (auto& water : FluidObject) {
            setBoundaryTypes(water, Join(RigidBodyObject, SoftBodyObject));
            water->setMinDepthFromCORNER();
            if (use_true_quadratic_element)
              computeAllCornerMidpointOffsets(water);
            const auto collapse_result = collapseCornerConnectedNeumannLinesAfterBoundaryTypes(*water, time_step);
            if (collapse_result.changed) {
              water->setContactFaces(Join(RigidBodyObject, SoftBodyObject));
              setBoundaryTypes(water, Join(RigidBodyObject, SoftBodyObject));
              water->setMinDepthFromCORNER();
              if (use_true_quadratic_element)
                computeAllCornerMidpointOffsets(water);
              std::cout << Magenta << "[corner_neumann_post_type] recategorized after "
                        << collapse_result.collapsed << " collapse(s)" << colorReset << std::endl;
            } else
              log_corner_connected_neumann_lines_after_boundary_types(water, "post-setBoundaryTypes");
          }
          rebuild_fluid_global_buckets();
          std::cout << Green << "setBoundaryTypes" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

          // Export kernel scatter data at the first time step (after boundary types are set)
          if (time_step == 0)
            exportKernelScatterData(FluidObject, output_directory);
        }

        // BIEで解く際も利用するので，吸収体に吸収される点を記録しておく
        for (const auto& p : Buckets_Fluid_Points.data1D) {
          p->absorbedBy = nullptr;
          for (const auto& net : AbsorberObject)
            if (net->InsideQ(p->X))
              p->absorbedBy = net;
        }

        /* --------------------------------------------------------------------------*/
        /*                       境界値問題の設定と解法                              */
        /* --------------------------------------------------------------------------*/

        setBodyVelocity(Join(RigidBodyObject, SoftBodyObject)); // まず，流体ではなく，剛体・柔体オブジェクトの速度を節点に設定する．
        std::cout << Green << "setBodyVelocity" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

        setPhiPhinOnFace(FluidObject);
        std::cout << Green << "setPhiPhinOnFace" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

        //% --------------------------------------------------------------------------*/
        //%              VPM Coupling: 2. 渦による誘導速度を境界条件に反映            */
        //% --------------------------------------------------------------------------*/
        // BEMのノイマン境界条件: dphi/dn = (u_body - u_omega) . n
        // setPhiPhinOnFace ですでに u_body . n が設定されているため、u_omega . n を引く

        for (auto water : FluidObject) {
          for (const auto& p : water->getPoints())
            p->u_omega_VPM.fill(0.0);
          for (auto* l : water->getBoundaryLines())
            l->u_omega_VPM.fill(0.0);
        }

        if (use_VPM) {
          for (auto water : FluidObject) {
            _Pragma("omp parallel for") for (auto p : ToVector(water->getPoints())) {
              p->u_omega_VPM = vpm.computeVelocity(p->X);
              if (!p->Neumann)
                continue;
              for (auto& [face, phin] : p->phinOnFace) {
                if (face)
                  phin -= Dot(p->u_omega_VPM, face->normal);
                else
                  phin -= Dot(p->u_omega_VPM, p->getNormalNeumann_BEM());
              }
              if (p->phinOnFace.count(nullptr))
                std::get<1>(p->phiphin) = p->phinOnFace.at(nullptr);
            }
            auto boundary_lines = water->getBoundaryLines();
            _Pragma("omp parallel for") for (size_t il = 0; il < boundary_lines.size(); ++il) {
              auto* l = boundary_lines[il];
              bool has_true_quad = std::ranges::any_of(l->getBoundaryFaces(), [](const auto* f) { return f->isTrueQuadraticElement; });
              if (has_true_quad)
                l->u_omega_VPM = vpm.computeVelocity(l->X_mid);
              if (!l->Neumann)
                continue;
              for (auto& [face, phin] : l->phinOnFace) {
                if (face)
                  phin -= Dot(l->u_omega_VPM, face->normal);
                else
                  phin -= Dot(l->u_omega_VPM, l->getNormal());
              }
              if (l->phinOnFace.count(nullptr))
                std::get<1>(l->phiphin) = l->phinOnFace.at(nullptr);
            }
          }
        }

        auto [time_setup_, time_solve_, unknownsize_] = BVP.solve();
        if (!BVP.last_gmres_converged) {
          const double gmres_retry_threshold = std::max(1e-3, solver_tol * 1e6);
          std::cerr << Yellow << "[GMRES] time_step " << time_step
                    << " RK_step " << RK_step
                    << ": residual=" << BVP.last_gmres_residual_norm
                    << " iter=" << BVP.last_gmres_total_iter
                    << colorReset << std::endl;
          BVP.logLastGmresHotspots(std::cerr, 8);
          if (BVP.last_gmres_residual_norm > gmres_retry_threshold) {
            throw step_failure("GMRES not converged at RK_step " + std::to_string(RK_step)
                               + ", residual=" + std::to_string(BVP.last_gmres_residual_norm));
          }
          std::cerr << Yellow << "[GMRES] RK_step " << RK_step
                    << ": residual " << BVP.last_gmres_residual_norm
                    << " above solver_tol=" << solver_tol
                    << " but below retry threshold=" << gmres_retry_threshold
                    << ", continuing" << colorReset << std::endl;
        }
        time_setup = time_setup_;
        time_solve = time_solve_;
        ElapsedTimeSetup.push_back(time_setup);
        ElapsedTimeSolve.push_back(time_solve);
        unknownsizse = unknownsize_;
        std::cout << Green << "BVP.solve -> {Φ,Φn}が決まる" << Blue << "\nElapsed time: " << Red << watch() << " s\n";

        for (auto water : FluidObject)
          set_u_potential_BEM(*water);
        std::cout << Green << "set_u_potential_BEM" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

        for (auto water : FluidObject)
          set_u_total(*water);
        std::cout << Green << "set_u_total" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

        const bool do_node_relocation = (node_relocation_method != NodeRelocationMethod::none)
                          && (RK_step == RK_order && time_step % node_relocation_period == 0);
        for (auto water : FluidObject)
          setNodeVelocity(*water, do_node_relocation ? 30 : 0, 0.05);
        std::cout << Green << "setNodeVelocity" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";

        //% --------------------------------------------------------------------------*/
        //%              VPM Coupling: 3. Time Integration (Advection+Stretching)     */
        //% --------------------------------------------------------------------------*/
        // Advance particles with the same RK stages as the main loop.
        // Diffusion is NOT included here (operator splitting).
        if (use_VPM) {
          vpm.set_u_potential_BEM([&](const std::array<double, 3>& x) { return getBEMVelocityAt_cached(x, FluidObject); });
          vpm.calcVelocityAndStretching(); // sets u_total and d_alpha_dt (stretching only)
          vpm.pushRK();                    // stage update for x and alpha

          // No-through safety: particles may enter rigid bodies during advection/stretching.
          // Push them out immediately so the next stage does not evaluate velocities inside solids.
          vpm.processParticles([&](VortexParticle& p) {
            for (const auto& rb : RigidBodyObject) {
              if (!inside_closed_surface(rb, p.x))
                continue;
              const double delta = std::max(0.3 * p.sigma, 1e-12);
              push_out_of_body(rb, p, delta);
            }
          });
        }

        /* --------------------------------------------------------------------------*/
        /*                                加速度の計算                               */
        /* --------------------------------------------------------------------------*/

        _Pragma("omp parallel for") for (const auto& net : Join(RigidBodyObject, SoftBodyObject)) { // \label{BEM:impose_velocity}
          //  重心位置と姿勢の時間発展
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

        /* --------------------------------------------------------------------------*/

        std::vector<Network*> movableObjects;
        for (auto& net : Join(RigidBodyObject, SoftBodyObject))
          if (std::ranges::any_of(net->isFixed, [](const auto& v) { return v == false; })) {
            movableObjects.push_back(net);
          }

        if (coupling_type == "NONE" || movableObjects.empty()) {
          convergence.clear();
          std::cout << Green << "Skipping BVP.solveForPhiPhin_t (no movable bodies or coupling=NONE)." << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
        } else {
          convergence = BVP.solveForPhiPhin_t(movableObjects);
          std::cout << Green << "BVP.solveForPhiPhin_t-> {Φt,Φtn}とnet->accelerationが決まる" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
        }

        /*DOC_EXTRACT 0_4_0_FLOATING_BODY_SIMULATION

        ### 浮体の重心位置・姿勢・速度の更新

        浮体の重心位置は，重心に関する運動方程式を解くことで求める．
        姿勢は，角運動量に関する運動方程式などを使って，各加速度を求める．姿勢はクオータニオンを使って表現する．

        */

        _Pragma("omp parallel for") for (const auto& net : Join(RigidBodyObject, SoftBodyObject)) { // \label{BEM:impose_velocity}
          //  重心位置と姿勢の時間発展
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

        // b$ --------------------------------------------------------------------------
        // b$                                波の吸収（ダンピング領域）
        // b$ --------------------------------------------------------------------------

        /*DOC_EXTRACT 0_4_1_UPDATE_POSITION

        ### 流体の$\phi$時間発展，$\phi_n$の時間発展はない

        ### 波の吸収（ダンピング領域）

        $$
        \begin{aligned}
        \gamma &= 1 - 2 \frac{\text{horizontal distance from the center of the
        absorber}}{\text{width of the absorber}} \\
        \phi_{\rm ref} &= \frac{\sum \phi \cdot \text{area}}{\sum \text{area}}
        \end{aligned}
        $$
        */

        // signed distanceの計算
        _Pragma("omp parallel for") for (const auto& p : ToVector(Buckets_Fluid_Points.data1D)) {
          if (p->absorbedBy != nullptr) {
            double min_distance = 1E+20;
            auto [f, X_nearest] = p->absorbedBy->Nearest(p->X);
            p->signed_distance_vector = p->X - X_nearest;
            p->signed_distance = Norm(p->signed_distance_vector);
          } else
            p->signed_distance = 0;
        }

        // 中点の signed_distance 計算（true quadratic）
        if (use_true_quadratic_element) {
          for (auto water : FluidObject) {
            for (auto* l : water->getBoundaryLines()) {
              if (!(l->Dirichlet || l->CORNER))
                continue;
              auto [pA, pB] = l->getPoints();
              if (pA->absorbedBy && pB->absorbedBy) {
                auto [f, X_near] = pA->absorbedBy->Nearest(l->X_mid);
                l->signed_distance = Norm(l->X_mid - X_near);
                l->absorbedBy = pA->absorbedBy;
              } else if (pA->absorbedBy || pB->absorbedBy) {
                auto absorber = pA->absorbedBy ? pA->absorbedBy : pB->absorbedBy;
                auto [f, X_near] = absorber->Nearest(l->X_mid);
                l->signed_distance = Norm(l->X_mid - X_near);
                l->absorbedBy = absorber;
              } else {
                l->signed_distance = 0.;
                l->absorbedBy = nullptr;
              }
            }
          }
        }

        // phiの平均
        double mean_phi = 0., total_area = 0;
        for (const auto& f : Buckets_Fluid_Faces.data1D) {
          auto [p0, p1, p2] = f->getPoints();
          mean_phi += (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3 * f->area;
          total_area += f->area;
        }
        if (!(total_area > 0.) || !std::isfinite(total_area))
          mean_phi = 0.;
        else
          mean_phi /= total_area;

        const bool debug_nan_guard = []() {
          if (const char* env = std::getenv("BEM_RELOCATION_NAN_DEBUG"))
            return std::string(env) != "0";
          return false;
        }();
        std::size_t nan_guard_point_count = 0;
        std::size_t nan_guard_phi_count = 0;

        auto pushPosition = [](RungeKutta<Tddd>& RK, const Tddd& u_mesh, const Tddd& fallback, std::size_t& nan_count) -> Tddd {
          RK.push(u_mesh);
          auto X = RK.getX();
          if (!isFinite(X)) {
            RK.repush(Tddd{0., 0., 0.});
            X = RK.getX();
            if (!isFinite(X))
              X = fallback;
            ++nan_count;
          }
          return X;
        };

        auto pushPhi = [](RungeKutta<double>& RK, double dphi_dt, std::size_t& nan_count) -> double {
          if (!std::isfinite(dphi_dt)) {
            dphi_dt = 0.;
            ++nan_count;
          }
          RK.push(dphi_dt);
          double phi = RK.getX();
          if (!std::isfinite(phi)) {
            RK.repush(0.);
            phi = RK.getX();
            if (!std::isfinite(phi))
              phi = 0.;
            ++nan_count;
          }
          return phi;
        };

        // 吸収帯の gamma, eta z-補正, ref_phi を計算する共通関数
        // u_mesh[2] を直接補正し、{gamma, ref_phi} を返す
        auto computeAbsorption = [&mean_phi](
                                     const Network* absorber, double signed_distance,
                                     RungeKutta<Tddd>& RK, Tddd& u_mesh, bool do_eta_phi,
                                     std::size_t& nan_pos_count, std::size_t& nan_phi_count) -> std::pair<double, double> {
          double gamma = absorber->absorb_gamma(signed_distance);
          if (!std::isfinite(gamma)) {
            gamma = 0.;
            ++nan_phi_count;
          }
          double ref_phi = 0.;
          if (do_eta_phi) {
            ref_phi = mean_phi;
            auto nextX = RK.getX(u_mesh);
            const double nextT = RK.getTimeAtNextStep();
            if (isFinite(nextX)) {
              const double eta = absorber->absorb_eta(nextX, nextT);
              const double dt_rk = RK.getdt();
              const double to_eta = eta - nextX[2];
              if (std::isfinite(dt_rk) && dt_rk > 0. && std::isfinite(to_eta)) {
                const double u2_new = u_mesh[2] + gamma * to_eta / dt_rk;
                if (std::isfinite(u2_new))
                  u_mesh[2] = u2_new;
                else
                  ++nan_pos_count;
              } else
                ++nan_pos_count;
              auto nextX_abs = RK.getX(u_mesh);
              const double phi_abs = isFinite(nextX_abs) ? absorber->absorb_phi(nextX_abs, nextT) : std::numeric_limits<double>::quiet_NaN();
              if (std::isfinite(phi_abs))
                ref_phi = phi_abs + mean_phi;
              else
                ++nan_phi_count;
            } else {
              ++nan_pos_count;
            }
          }
          return {gamma, ref_phi};
        };

        for (const auto& p : Buckets_Fluid_Points.data1D) {
          p->U_absorbed.fill(0.);
          double gamma = 0, ref_phi = 0;
          if (p->absorbedBy != nullptr) {
            bool has_dirichlet = std::ranges::any_of(p->getBoundaryFaces(), [](auto f) { return f->Dirichlet; });
            auto uz_before = p->u_node[2];
            std::tie(gamma, ref_phi) = computeAbsorption(p->absorbedBy, p->signed_distance,
                                                         p->RK_X, p->u_node, has_dirichlet,
                                                         nan_guard_point_count, nan_guard_phi_count);
            p->U_absorbed[2] = p->u_node[2] - uz_before;
          }

          if (!isFinite(p->u_node)) {
            p->u_node = isFinite(p->u_total) ? p->u_total : Tddd{0., 0., 0.};
            ++nan_guard_point_count;
          }
          if (!isFinite(p->RK_X.getX(p->u_node))) {
            p->u_node = {0., 0., 0.};
            ++nan_guard_point_count;
          }

          if (!p->Neumann) {
            const Tddd& U_dphi = (node_relocation_method == NodeRelocationMethod::interpolation) ? p->u_total : p->u_node;
            double dphi_dt = p->DphiDt_damped({gamma, ref_phi}, U_dphi, 0.);
            std::get<0>(p->phiphin) = p->phi_Dirichlet = pushPhi(p->RK_phi, dphi_dt, nan_guard_phi_count);
          }
          p->setXSingle(pushPosition(p->RK_X, p->u_node, getPosition(p), nan_guard_point_count));
          p->phi_tmp = 0;
        }
        if (debug_nan_guard && (nan_guard_point_count > 0 || nan_guard_phi_count > 0)) {
          std::cout << Magenta << "[ALE:nan-guard:point] corrected_u=" << nan_guard_point_count
                    << " corrected_phi=" << nan_guard_phi_count << colorReset << std::endl;
        }

        if (use_true_quadratic_element)
          for (auto water : FluidObject)
            for (auto* l : water->getBoundaryLines()) {
              l->U_absorbed.fill(0.);
              double gamma = 0, ref_phi = 0;
              if (l->absorbedBy != nullptr) {
                bool has_dirichlet = std::ranges::any_of(l->getBoundaryFaces(), [](auto f) { return f->Dirichlet; });
                auto uz_before = l->u_node[2];
                std::tie(gamma, ref_phi) = computeAbsorption(l->absorbedBy, l->signed_distance,
                                                             l->RK_X, l->u_node, has_dirichlet,
                                                             nan_guard_point_count, nan_guard_phi_count);
                l->U_absorbed[2] = l->u_node[2] - uz_before;
              }

              if (!isFinite(l->u_node)) {
                l->u_node = isFinite(l->u_total) ? l->u_total : Tddd{0., 0., 0.};
                ++nan_guard_point_count;
              }
              if (!isFinite(l->RK_X.getX(l->u_node))) {
                l->u_node = {0., 0., 0.};
                ++nan_guard_point_count;
              }

              if (!l->Neumann) {
                const Tddd& U_dphi = (node_relocation_method == NodeRelocationMethod::interpolation) ? l->u_total : l->u_node;
                double dphi_dt = l->DphiDt_damped({gamma, ref_phi}, U_dphi, 0.);
                std::get<0>(l->phiphin) = l->phi_Dirichlet = pushPhi(l->RK_phi, dphi_dt, nan_guard_phi_count);
              }
              l->setXSingle(pushPosition(l->RK_X, l->u_node, getPosition(l), nan_guard_point_count));
              l->phi_tmp = 0;
            }
        if (debug_nan_guard && (nan_guard_point_count > 0 || nan_guard_phi_count > 0)) {
          std::cout << Magenta << "[ALE:nan-guard:point] corrected_u=" << nan_guard_point_count
                    << " corrected_phi=" << nan_guard_phi_count << colorReset << std::endl;
        }

        // // 中点の RK push（true quadratic）
        // if (use_true_quadratic_element) {
        //   for (auto water : FluidObject) {
        //     std::size_t nan_guard_mid_count = 0;
        //     std::size_t nan_guard_mid_phi_count = 0;

        //     for (auto* l : water->getBoundaryLines()) {
        //       if (l->Dirichlet || l->CORNER) {
        //         double gamma_mid = 0., ref_phi_mid = 0.;
        //         if (l->absorbedBy)
        //           std::tie(gamma_mid, ref_phi_mid) = computeAbsorption(l->absorbedBy, l->signed_distance,
        //                                                                l->RK_X, l->u_node, true,
        //                                                                nan_guard_mid_count, nan_guard_mid_phi_count);
        //         double dphi_mid_dt = -gamma_mid * (l->phiphin[0] - ref_phi_mid) + DphiDt_at_midpoint(l);
        //         l->phiphin[0] = pushPhi(l->RK_phi, dphi_mid_dt, nan_guard_mid_phi_count);
        //       }
        //       // Position evolution (all boundary edge midpoints)
        //       auto [pA_line, pB_line] = l->getPoints();
        //       Tddd mid_fallback = 0.5 * (pA_line->X + pB_line->X);
        //       if (l->CORNER)
        //         mid_fallback += l->corner_midpoint_offset;
        //       l->X_mid = pushPosition(l->RK_X, l->u_node, mid_fallback, nan_guard_mid_count);
        //       l->Xtarget = l->X_mid;
        //     }
        //     if (debug_nan_guard && (nan_guard_mid_count > 0 || nan_guard_mid_phi_count > 0)) {
        //       std::cout << Magenta << "[ALE:nan-guard:mid] water=" << water->getName()
        //                 << " corrected_u=" << nan_guard_mid_count
        //                 << " corrected_phi=" << nan_guard_mid_phi_count << colorReset << std::endl;
        //     }
        //   }
        // }

        // b$

        for (auto water : FluidObject)
          std::cout << Green << "name:" << water->getName() << ": setBounds" << colorReset << std::endl;

        for (auto net : AllObjects)
          net->setGeometricPropertiesForce();
        for (auto water : FluidObject)
          monitorTinyFaceRelativeToLocalMean(*water, time_step, RK_step);
        for (auto water : FluidObject)
          throwIfSubsurfaceFaceAltitudeTooSmall(*water, time_step, RK_step, subsurface_altitude_reject);
        for (auto water : FluidObject)
          throwIfTinyFaceRelativeToLocalMean(*water, time_step, RK_step);

        for (auto water : FluidObject) {
          auto name = water->getName();
          std::ofstream ofs(output_directory / (name + std::to_string(RK_step) + ".obj"));
          createOBJ(ofs, *water);
          ofs.close();
        }

        std::cout << Blue << "Elapsed time: " << Red << watch() << colorReset << " s\n";
      } while (!((*(Buckets_Fluid_Points.data1D).begin())->RK_X.finished));

      sync_body_states_from_rk();
      for (auto net : AllObjects)
        net->setGeometricPropertiesForce();
      throwIfStructurePenetrated(FluidObject, Join(RigidBodyObject, SoftBodyObject), time_step, "post-RK");

      // POST-RK4 interpolation: use cached (face, t0, t1) from calculateVecToSurface
      // to interpolate phi from the DphiDt(u_total) distribution to smoothed positions
      if (node_relocation_method == NodeRelocationMethod::interpolation && time_step % node_relocation_period == 0) {
        for (auto water : FluidObject) {
          // 1. Save DphiDt(u_total) phi distribution before overwriting
          std::unordered_map<networkPoint*, double> lagrangian_phi;
          for (auto* p : water->getPoints())
            if (!p->Neumann || p->CORNER)
              lagrangian_phi[p] = std::get<0>(p->phiphin);

          std::unordered_map<networkLine*, double> lagrangian_phi_mid;
          if (use_true_quadratic_element)
            for (auto* l : water->getBoundaryLines())
              if (l->Dirichlet || l->CORNER)
                lagrangian_phi_mid[l] = l->phiphin[0];

          // 2. Points: interpolate phi using cached (face, t0, t1)
          for (auto* p : water->getPoints()) {
            if (p->Neumann && !p->CORNER) continue;
            if (!p->relocation_face) continue;

            auto* f = p->relocation_face;
            auto [t0, t1] = p->relocation_param;
            double phi_new;

            if (node_relocation_surface == NodeRelocationSurface::true_quadratic) {
              auto [p0, l0, p1, l1, p2, l2] = f->PLPLPL;
              auto N = TriShape<6>(t0, t1);
              phi_new = N[0] * lagrangian_phi[p0] + N[1] * lagrangian_phi[p1] + N[2] * lagrangian_phi[p2]
                      + N[3] * lagrangian_phi_mid[l0] + N[4] * lagrangian_phi_mid[l1] + N[5] * lagrangian_phi_mid[l2];
            } else if (node_relocation_surface == NodeRelocationSurface::pseudo_quadratic && f->dodecaPoints[0]) {
              phi_new = f->dodecaPoints[0]->interpolate(t0, t1,
                  [&](networkPoint* q) -> double {
                    auto it = lagrangian_phi.find(q);
                    return (it != lagrangian_phi.end()) ? it->second : 0.0;
                  });
            } else {
              auto [q0, q1, q2] = f->getPoints();
              double wc = 1.0 - t0 - t1;
              phi_new = t0 * lagrangian_phi[q0] + t1 * lagrangian_phi[q1] + wc * lagrangian_phi[q2];
            }
            std::get<0>(p->phiphin) = p->phi_Dirichlet = phi_new;
          }

          // 3. Lines: true_quadratic interpolation using cached (face, t0, t1)
          if (use_true_quadratic_element) {
            for (auto* l : water->getBoundaryLines()) {
              if (!(l->Dirichlet || l->CORNER)) continue;
              if (!l->relocation_face) continue;

              auto* f = l->relocation_face;
              auto [t0, t1] = l->relocation_param;
              auto [p0, l0, p1, l1, p2, l2] = f->PLPLPL;
              auto N = TriShape<6>(t0, t1);
              l->phiphin[0] = N[0] * lagrangian_phi[p0] + N[1] * lagrangian_phi[p1] + N[2] * lagrangian_phi[p2]
                            + N[3] * lagrangian_phi_mid[l0] + N[4] * lagrangian_phi_mid[l1] + N[5] * lagrangian_phi_mid[l2];
            }
          }
        }
        std::cout << Green << "Interpolation relocation completed" << Blue << "\nElapsed time: " << Red << watch() << colorReset << " s\n";
      }

      std::cout << Red << "Total elapsed time: " << Red << watch()[1] << colorReset << " s\n";
      ElapsedTimeTotal = watch()[1];

      /* ---------------------------------- 所要時間 ---------------------------------- */
      std::cout << "============================================================================" << std::endl;
      std::cout << "ElapsedTimeSetup=" << ElapsedTimeSetup << std::endl;
      std::cout << "ElapsedTimeSolve=" << ElapsedTimeSolve << std::endl;
      std::cout << "ElapsedTimeNodeRelocation=" << ElapsedTimeNodeRelocation << std::endl;
      std::cout << "ElapsedTimeTotal=" << ElapsedTimeTotal << std::endl;
      std::cout << "============================================================================" << std::endl;
      /* --------------------------------------------------------------------------*/

      std::cout << Green << "simulation_timeを取得" << colorReset << std::endl;
      simulation_time = (*(Buckets_Fluid_Points.data1D).begin())->RK_X.gett();

      //@ --------------------------------------------------------------------------*/
      //@                                VPM Step                                   */
      //@ --------------------------------------------------------------------------*/
      // VPM time-step (Strang splitting):
      // 1) diffusion half-step (done at the beginning of the time-step)
      // 2) advection + stretching (advanced in the RK loop; no diffusion)
      // 3) boundary vorticity injection (absorb -> shed)
      // 4) diffusion half-step (new + old particles)
      if (use_VPM) {
        // Safety: ensure wall/body velocities are evaluated at the final RK stage (t_{n+1})
        // before computing slip velocity and injecting wall vorticity.
        setBodyVelocity(Join(RigidBodyObject, SoftBodyObject));

        // (3) Shedding / absorption at the wall
        for (auto water : FluidObject) {
          // Per-body `shed_vortices=false` means *no* vorticity is generated at all.
          // If `shed_vortices=true`, we still avoid creating new particles when near-wall receivers exist.
          vpm.injectWallVorticityFluxPSE(water, dt, VPM_wall_min_absorb_receivers, VPM_wall_min_absorb_total_weight, VPM_sigma_factor, true /*allow_shed*/);
        }

        // (4) Diffusion (half-step) including the newly injected particles
        vpm.applyDiffusionEuler(0.5 * dt);

        // --------------------------------------------------------------------------
        // 粒子のメンテナンス1: 安全弁（壁面内に入った粒子を反射させる）
        // --------------------------------------------------------------------------
        vpm.processParticles([&](VortexParticle& p) {
          for (const auto& rb : RigidBodyObject) {
            if (inside_closed_surface(rb, p.x)) {
              // 壁面内部に入り込んだ場合、最近傍点を使って反射させる
              // x_new = x - 2 * ( (x - x_surface) . n ) * n
              const double delta = std::max(0.3 * p.sigma, 1e-12);
              push_out_of_body(rb, p, delta);
            }
          }
        });

        // --------------------------------------------------------------------------
        // 粒子のメンテナンス2: 領域外や物体内部(反射しきれなかったもの)の粒子を削除
        // --------------------------------------------------------------------------
        vpm.removeParticles([&](const VortexParticle& p) {
          // 1. 剛体内部への貫入判定
          for (const auto& rb : RigidBodyObject) {
            if (inside_closed_surface(rb, p.x))
              return true;
          }
          // 2. 流体領域外判定 (Waterの外部に出たら削除)
          bool inside_any_water = false;
          for (const auto& water : FluidObject) {
            if (inside_closed_surface(water, p.x)) {
              inside_any_water = true;
              break;
            }
          }
          return !inside_any_water; // どの水領域にも入っていなければ削除
        });

        // Note: Diffusion is handled by operator splitting:
        // - dt/2 at the beginning of the time-step
        // - dt/2 here (after shedding), including newly injected particles
      }

      //! 速度ポテンシャルの平均を0にする
      {
        double mean_phi = 0., total_area = 0;
        for (const auto& f : Buckets_Fluid_Faces.data1D) {
          auto [p0, p1, p2] = f->getPoints();
          mean_phi += (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3 * f->area;
          total_area += f->area;
        }
        mean_phi /= total_area;
        // mean_phi /= Buckets_Fluid_Points.data1D.size();
        for (const auto& p : Buckets_Fluid_Points.data1D) {
          p->phi_Dirichlet -= mean_phi;
          // p->phi_Neumann -= mean_phi;
          std::get<0>(p->phiphin) -= mean_phi;
        }
        // Midpoint phi_mid must be shifted by the same amount (true quadratic elements).
        // All midpoints (Dirichlet, Neumann, CORNER) must be shifted so that
        // the midpoint phi stays consistent with the shifted vertex phi.
        if (use_true_quadratic_element) {
          for (auto* water : FluidObject)
            for (auto* l : water->getBoundaryLines())
              l->phiphin[0] -= mean_phi;
        }
      }

      //! クォータニオンの正規化
      for (const auto& net : RigidBodyObject)
        net->Q.normalize();

      // Fixed-body excitation-force runs: we skip BVP.solveForPhiPhin_t, so pressure is not updated.
      // For force output, approximate phi_t with a first-order finite difference at the output time-step.
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

      /* ------------------------------------------------------ */

      // Step succeeded — break out of retry loop
      last_successful_dt = dt;
      break;

        } catch (step_failure& e) {
          std::cerr << Red << "[step_reject] time_step " << time_step << ": " << e.reason << colorReset << std::endl;
          remember_rejected_face(e);
          if (step_retry >= max_step_retries) {
            std::cerr << Yellow << "[step_reject] writing failure snapshot before abort" << colorReset << std::endl;
            write_failure_snapshot(current_rk_step, step_retry);
            std::cerr << Red << "[step_reject] max retries exceeded, giving up" << colorReset << std::endl;
            throw; // propagate as fatal
          }
          // Reduce dt for next retry
          if (dt_override <= 0)
            dt_override = std::max(1E-13, std::min(dt, last_successful_dt) * dt_reduction_factor);
          else
            dt_override = std::max(1E-13, dt_override * dt_reduction_factor);
          continue; // retry
        }
      } // end retry loop

      // Clean up temporary checkpoint
      std::filesystem::remove(step_checkpoint_path);

      write_step_outputs(false);
      /* --------------------------------------------------------------------------*/
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  };
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
