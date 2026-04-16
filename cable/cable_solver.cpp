/**
 * @file cable_solver.cpp
 * @brief Cable solver — JSON I/O wrapper around LumpedCableSystem.
 *
 * Usage:
 *   ./cable_solver input.json output_dir/
 *
 * Build:
 *   mkdir -p build && cd build
 *   cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=cable/cable_solver.cpp -DOUTPUT_NAME=cable_solver ../..
 *   make -j$(sysctl -n hw.logicalcpu)
 *
 * Input JSON formats accepted (auto-detected):
 *
 *   1. SETTINGS MODE
 *      Has `input_files` key listing per-cable JSON files (BEM pattern).
 *      Shared scalars (gravity, mode, etc.) are defaults for each cable.
 *      Each cable file is loaded and solved independently.
 *
 *   2. PER-CABLE FORMAT (new)
 *      Has `end_a_position` key. Self-contained single-cable JSON.
 *      Optional `end_a_body`/`end_b_body` for BEM linkage labels.
 *      Optional `end_a_motion`/`end_b_motion` for dynamic mode (Phase B).
 *
 *   3. MULTI-LINE BEM-COMPATIBLE FORMAT
 *      One or more `mooring_<name>` / `cable_<name>` keys with 13-element
 *      flat arrays. Same schema as the BEM input reader.
 *
 *   4. SINGLE-LINE LEGACY FORMAT
 *      Flat keys: point_a, point_b, cable_length, n_segments, etc.
 */

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "LumpedCable.hpp"
#include "WindField.hpp"
#include "basic.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "basic_vectors.hpp"
#include "my_vtk.hpp"
#include "vtkWriter.hpp"

namespace fs = std::filesystem;

/* -------------------------------------------------------------------------- */
/*                          Topological point order                            */
/* -------------------------------------------------------------------------- */

// Walk the cable from firstPoint to lastPoint in topological order.
static std::vector<networkPoint*> getOrderedPoints(const LumpedCable* mooring) {
   std::vector<networkPoint*> ordered;
   auto current = mooring->firstPoint;
   networkPoint* prev = nullptr;
   while (current) {
      ordered.push_back(current);
      if (current == mooring->lastPoint) break;
      networkPoint* next = nullptr;
      for (auto l : current->getLines()) {
         auto other = (*l)(current);
         if (other != prev) {
            next = other;
            break;
         }
      }
      prev = current;
      current = next;
   }
   return ordered;
}

/* -------------------------------------------------------------------------- */
/*                      ParaView (.vtp) writer for a cable                    */
/* -------------------------------------------------------------------------- */

// Zero-pad an integer to `width` digits (e.g. 42 → "0042").
static std::string zfill_int(int value, int width = 5) {
   std::ostringstream ss;
   ss << std::setw(width) << std::setfill('0') << value;
   return ss.str();
}

// Write a single-timestep .vtp file containing:
//   - Points: ordered node positions along the cable
//   - Lines:  (n-1) segments, one cell per segment
//   - CellData "tension" [N]       (per-segment)
//   - CellData "strain"  [-]       (per-segment; tension-only clamped to >=0)
//   - PointData "velocity_mag" [m/s] (per-node)
//
// `points` is the ordered list returned by getOrderedPoints().
// `seg_tension` and `seg_strain` must have size `points.size() - 1`.
static void writeCableVTP(const fs::path& path,
                          const std::vector<networkPoint*>& points,
                          const std::vector<double>& seg_tension,
                          const std::vector<double>& seg_strain) {
   vtkPolygonWriter<networkPoint*> w;
   w.reserve(points.size());

   for (auto* p : points)
      w.add(p);

   // One polyline per segment so per-segment CellData maps 1:1 with cells.
   for (size_t i = 0; i + 1 < points.size(); ++i)
      w.addLine(std::vector<networkPoint*>{points[i], points[i + 1]});

   w.addCellData("tension", seg_tension);
   w.addCellData("strain", seg_strain);

   std::unordered_map<networkPoint*, double> vel_mag;
   vel_mag.reserve(points.size());
   for (auto* p : points) {
      double vm = std::sqrt(p->velocity[0] * p->velocity[0]
                            + p->velocity[1] * p->velocity[1]
                            + p->velocity[2] * p->velocity[2]);
      vel_mag[p] = vm;
   }
   w.addPointData("velocity_mag", vel_mag);

   std::ofstream ofs(path);
   w.write(ofs);
   ofs.close();
}

/* -------------------------------------------------------------------------- */
/*                         JSON result writer (single line)                    */
/* -------------------------------------------------------------------------- */

static void writeSingleResultJSON(const fs::path& filepath,
                                  const LumpedCable* mooring,
                                  double computation_time_ms,
                                  bool converged) {
   auto points = getOrderedPoints(mooring);
   int n_nodes = points.size();
   double natural_length = mooring->natural_length();

   std::ofstream ofs(filepath);
   if (!ofs.is_open()) {
      std::cerr << "ERROR: cannot open " << filepath << std::endl;
      return;
   }

   ofs << std::setprecision(10);
   ofs << "{\n";
   ofs << "  \"n_nodes\": " << n_nodes << ",\n";

   // positions
   ofs << "  \"positions\": [";
   for (int i = 0; i < n_nodes; ++i) {
      auto X = points[i]->X;
      ofs << "[" << std::get<0>(X) << "," << std::get<1>(X) << "," << std::get<2>(X) << "]";
      if (i < n_nodes - 1) ofs << ",";
   }
   ofs << "],\n";

   // Compute tensions directly from node positions (not via RK_X_sub)
   double stiffness_val = 0;
   for (auto l : mooring->getLines()) { stiffness_val = l->stiffness; break; }

   std::vector<double> seg_tension(n_nodes - 1, 0.);
   for (int i = 0; i < n_nodes - 1; ++i) {
      Tddd dv = points[i + 1]->X - points[i]->X;
      double dist = Norm(dv);
      double strain = (dist - natural_length) / natural_length;
      if (strain > 0.)
         seg_tension[i] = stiffness_val * strain;
   }

   // Per-node tension: average of adjacent segments
   ofs << "  \"tensions\": [";
   for (int i = 0; i < n_nodes; ++i) {
      double t = 0;
      int count = 0;
      if (i > 0) { t += seg_tension[i - 1]; count++; }
      if (i < n_nodes - 1) { t += seg_tension[i]; count++; }
      if (count > 0) t /= count;
      ofs << t;
      if (i < n_nodes - 1) ofs << ",";
   }
   ofs << "],\n";

   double top_tension = (n_nodes > 1) ? seg_tension.back() : 0;
   double bottom_tension = (n_nodes > 1) ? seg_tension.front() : 0;
   double max_tension = *std::max_element(seg_tension.begin(), seg_tension.end());

   ofs << "  \"top_tension\": " << top_tension << ",\n";
   ofs << "  \"bottom_tension\": " << bottom_tension << ",\n";
   ofs << "  \"max_tension\": " << max_tension << ",\n";
   ofs << "  \"converged\": " << (converged ? "true" : "false") << ",\n";
   ofs << "  \"computation_time_ms\": " << computation_time_ms << "\n";
   ofs << "}\n";
   ofs.close();

   std::cout << "Result written to " << filepath << std::endl;
}

/* -------------------------------------------------------------------------- */
/*                         JSON result writer (multi line)                    */
/* -------------------------------------------------------------------------- */

static void writeMultiResultJSON(const fs::path& filepath,
                                 const LumpedCableSystem& system,
                                 double computation_time_ms,
                                 bool converged) {
   std::ofstream ofs(filepath);
   if (!ofs.is_open()) {
      std::cerr << "ERROR: cannot open " << filepath << std::endl;
      return;
   }

   ofs << std::setprecision(10);
   ofs << "{\n";
   ofs << "  \"n_cables\": " << system.size() << ",\n";
   ofs << "  \"converged\": " << (converged ? "true" : "false") << ",\n";
   ofs << "  \"computation_time_ms\": " << computation_time_ms << ",\n";
   ofs << "  \"cables\": {\n";

   const auto& cables = system.cables();
   for (size_t ci = 0; ci < cables.size(); ++ci) {
      auto* mooring = cables[ci];
      auto points = getOrderedPoints(mooring);
      int n_nodes = points.size();
      double natural_length = mooring->natural_length();

      double stiffness_val = 0;
      for (auto l : mooring->getLines()) { stiffness_val = l->stiffness; break; }

      std::vector<double> seg_tension(n_nodes - 1, 0.);
      for (int i = 0; i < n_nodes - 1; ++i) {
         Tddd dv = points[i + 1]->X - points[i]->X;
         double dist = Norm(dv);
         double strain = (dist - natural_length) / natural_length;
         if (strain > 0.)
            seg_tension[i] = stiffness_val * strain;
      }
      double top_tension = (n_nodes > 1) ? seg_tension.back() : 0;
      double bottom_tension = (n_nodes > 1) ? seg_tension.front() : 0;
      double max_tension = *std::max_element(seg_tension.begin(), seg_tension.end());

      ofs << "    \"" << system.nameOf(ci) << "\": {\n";
      ofs << "      \"n_nodes\": " << n_nodes << ",\n";

      ofs << "      \"positions\": [";
      for (int i = 0; i < n_nodes; ++i) {
         auto X = points[i]->X;
         ofs << "[" << std::get<0>(X) << "," << std::get<1>(X) << "," << std::get<2>(X) << "]";
         if (i < n_nodes - 1) ofs << ",";
      }
      ofs << "],\n";

      ofs << "      \"tensions\": [";
      for (int i = 0; i < n_nodes; ++i) {
         double t = 0;
         int count = 0;
         if (i > 0) { t += seg_tension[i - 1]; count++; }
         if (i < n_nodes - 1) { t += seg_tension[i]; count++; }
         if (count > 0) t /= count;
         ofs << t;
         if (i < n_nodes - 1) ofs << ",";
      }
      ofs << "],\n";

      ofs << "      \"top_tension\": " << top_tension << ",\n";
      ofs << "      \"bottom_tension\": " << bottom_tension << ",\n";
      ofs << "      \"max_tension\": " << max_tension << "\n";
      ofs << "    }" << (ci + 1 < cables.size() ? "," : "") << "\n";
   }

   ofs << "  }\n";
   ofs << "}\n";
   ofs.close();

   std::cout << "Multi-cable result written to " << filepath << std::endl;
}

/* -------------------------------------------------------------------------- */
/*                         PrescribedMotion                                    */
/* -------------------------------------------------------------------------- */

struct PrescribedMotion {
   enum Kind { Fixed, Sinusoidal, Cantilever };
   Kind kind = Fixed;
   int dof = 2;           // 0=surge, 1=sway, 2=heave

   // Mode components. For backward compat, modes[0] is always the primary.
   // Sinusoidal uses modes[0] only (unless sinusoidal_multi injection added
   // extra components). Cantilever may have modes[0]=mode1, modes[1]=mode2.
   std::vector<double> amplitudes  = {0.0};
   std::vector<double> freqs_hz    = {0.0};
   std::vector<double> phases_rad  = {0.0};

   // Cantilever-only geometry (shared across all modes).
   Tddd fix_point = {0., 0., 0.};
   Tddd axis_dir = {0., 0., 1.};
   Tddd bend_dir = {1., 0., 0.};
   double length = 1.0;

   // Backward-compat scalar accessors (return modes[0]).
   double amplitude()  const { return amplitudes.empty()  ? 0. : amplitudes[0]; }
   double freq_hz_0()  const { return freqs_hz.empty()    ? 0. : freqs_hz[0]; }
   double phase_rad_0() const { return phases_rad.empty() ? 0. : phases_rad[0]; }

   Tddd velocity(double t) const {
      if (kind == Fixed) return {0., 0., 0.};
      double v = 0.;
      for (size_t n = 0; n < amplitudes.size(); ++n) {
         double w = 2. * M_PI * freqs_hz[n];
         v += w * amplitudes[n] * std::cos(w * t + phases_rad[n]);
      }
      Tddd result = {0., 0., 0.};
      if (dof == 0) std::get<0>(result) = v;
      else if (dof == 1) std::get<1>(result) = v;
      else std::get<2>(result) = v;
      return result;
   }

   Tddd displacement(double t) const {
      if (kind == Fixed) return {0., 0., 0.};
      double d = 0.;
      for (size_t n = 0; n < amplitudes.size(); ++n) {
         double w = 2. * M_PI * freqs_hz[n];
         d += amplitudes[n] * std::sin(w * t + phases_rad[n]);
      }
      Tddd result = {0., 0., 0.};
      if (dof == 0) std::get<0>(result) = d;
      else if (dof == 1) std::get<1>(result) = d;
      else std::get<2>(result) = d;
      return result;
   }
};

static int parseDof(const std::string& s) {
   if (s == "surge") return 0;
   if (s == "sway") return 1;
   return 2;  // heave
}

// Euler-Bernoulli clamped-free beam mode shape ψ_n(ξ), ξ ∈ [0, 1].
// Characteristic equation: cos(βL) cosh(βL) + 1 = 0.
// mode ∈ {1, 2, 3}. Normalized so ψ_n(1) = 1 (tip = tip amplitude).
static double cantileverMode(int mode, double xi) {
   static constexpr double BETAS[]  = {1.87510407, 4.69409113, 7.85475744};
   static constexpr double SIGMAS[] = {0.73409551, 1.01846732, 0.99922450};
   if (mode < 1 || mode > 3) return 0.0;
   if (xi <= 0.0) return 0.0;
   const double beta = BETAS[mode - 1];
   const double sigma = SIGMAS[mode - 1];
   auto psi_raw = [&](double s) {
      double bs = beta * s;
      return (std::cosh(bs) - std::cos(bs)) - sigma * (std::sinh(bs) - std::sin(bs));
   };
   // Pre-computed tip values ψ_raw(1) for normalization.
   static const double TIP[3] = {
      (std::cosh(BETAS[0]) - std::cos(BETAS[0])) - SIGMAS[0] * (std::sinh(BETAS[0]) - std::sin(BETAS[0])),
      (std::cosh(BETAS[1]) - std::cos(BETAS[1])) - SIGMAS[1] * (std::sinh(BETAS[1]) - std::sin(BETAS[1])),
      (std::cosh(BETAS[2]) - std::cos(BETAS[2])) - SIGMAS[2] * (std::sinh(BETAS[2]) - std::sin(BETAS[2])),
   };
   if (xi >= 1.0) return 1.0;
   return psi_raw(xi) / TIP[mode - 1];
}

static double cantileverFirstMode(double xi) { return cantileverMode(1, xi); }

// Map bend_dir to nearest cardinal DOF index (0=x, 1=y, 2=z).
static int cardinalDof(const Tddd& v) {
   double ax = std::abs(std::get<0>(v));
   double ay = std::abs(std::get<1>(v));
   double az = std::abs(std::get<2>(v));
   if (ax >= ay && ax >= az) return 0;
   if (ay >= az) return 1;
   return 2;
}

// Compute the local sinusoidal amplitude(s) for a cantilever body at a given
// attachment point. Returns one value per mode. Mode n uses cantileverMode(n, ξ).
static std::vector<double> cantileverLocalAmplitudes(
    const PrescribedMotion& m, const Tddd& attach_pos) {
   if (m.kind != PrescribedMotion::Cantilever) return {};
   Tddd r = attach_pos - m.fix_point;
   double s = std::get<0>(r) * std::get<0>(m.axis_dir)
            + std::get<1>(r) * std::get<1>(m.axis_dir)
            + std::get<2>(r) * std::get<2>(m.axis_dir);
   double xi = (m.length > 0.) ? (s / m.length) : 0.;
   std::vector<double> out;
   out.reserve(m.amplitudes.size());
   for (size_t n = 0; n < m.amplitudes.size(); ++n)
      out.push_back(m.amplitudes[n] * cantileverMode(static_cast<int>(n + 1), xi));
   return out;
}

// Backward-compat: returns mode-1 amplitude only.
static double cantileverLocalAmplitude(const PrescribedMotion& m,
                                       const Tddd& attach_pos) {
   auto amps = cantileverLocalAmplitudes(m, attach_pos);
   return amps.empty() ? 0. : amps[0];
}

// Parse BEM-style `velocity` field from a body JSON.
// Supported forms:
//   "velocity": "fixed"                                   (string)
//   "velocity": "floating"                                (treated as fixed for now)
//   "velocity": ["sinusoidal"|"sin"|"cos", start_time,
//                amplitude, period, axis_x, axis_y, axis_z]
// Returns PrescribedMotion with kind=Sinusoidal if parseable; Fixed otherwise.
// `motion_type_out` receives a human-readable label for logging.
static PrescribedMotion parseBodyVelocity(const JSON& body_json,
                                          std::string& motion_type_out) {
   PrescribedMotion pm;
   motion_type_out = "fixed";
   if (!body_json.find("velocity")) return pm;
   auto v = body_json.at("velocity");
   if (v.empty()) return pm;

   std::string mtype = v[0];
   if (mtype == "fixed" || mtype == "floating") {
      motion_type_out = mtype;
      return pm;  // Fixed motion
   }
   if (mtype == "sinusoidal" || mtype == "sin" || mtype == "cos") {
      if (v.size() < 7) {
         std::cerr << "WARNING: velocity array for sinusoidal needs 7 elements "
                   << "[type, start_time, amplitude, period, ax, ay, az], got "
                   << v.size() << std::endl;
         return pm;
      }
      double start_time = std::stod(v[1]);
      double amplitude = std::stod(v[2]);
      double period = std::stod(v[3]);
      double ax = std::stod(v[4]);
      double ay = std::stod(v[5]);
      double az = std::stod(v[6]);
      double absx = std::abs(ax), absy = std::abs(ay), absz = std::abs(az);
      int dof = (absx >= absy && absx >= absz) ? 0 : (absy >= absz ? 1 : 2);
      pm.kind = PrescribedMotion::Sinusoidal;
      pm.dof = dof;
      double fhz = (period > 0.) ? 1. / period : 0.;
      pm.amplitudes  = {amplitude};
      pm.freqs_hz    = {fhz};
      pm.phases_rad  = {-2. * M_PI * fhz * start_time};
      motion_type_out = "sinusoidal";
      return pm;
   }
   if (mtype == "cantilever") {
      // 14-element: mode 1 only (backward compat)
      //   ["cantilever", start_time, tip_amp, period,
      //    fix(3), axis(3), bend(3), length]
      // 17-element: mode 1 + mode 2
      //   [...14 as above..., tip_amp2, period2, phase2_rad]
      if (v.size() < 14) {
         std::cerr << "WARNING: velocity array for cantilever needs 14 or 17 elements, got "
                   << v.size() << std::endl;
         return pm;
      }
      double start_time = std::stod(v[1]);
      double tip_amp = std::stod(v[2]);
      double period = std::stod(v[3]);
      Tddd fix = {std::stod(v[4]), std::stod(v[5]), std::stod(v[6])};
      Tddd axis = {std::stod(v[7]), std::stod(v[8]), std::stod(v[9])};
      Tddd bend = {std::stod(v[10]), std::stod(v[11]), std::stod(v[12])};
      double beam_length = std::stod(v[13]);
      auto normalize = [](Tddd& u) {
         double n = std::sqrt(std::get<0>(u)*std::get<0>(u)
                            + std::get<1>(u)*std::get<1>(u)
                            + std::get<2>(u)*std::get<2>(u));
         if (n > 0.) { std::get<0>(u)/=n; std::get<1>(u)/=n; std::get<2>(u)/=n; }
      };
      normalize(axis);
      normalize(bend);
      double freq1 = (period > 0.) ? 1. / period : 0.;
      double phase1 = -2. * M_PI * freq1 * start_time;

      pm.kind = PrescribedMotion::Cantilever;
      pm.dof = cardinalDof(bend);
      pm.amplitudes  = {tip_amp};
      pm.freqs_hz    = {freq1};
      pm.phases_rad  = {phase1};
      pm.fix_point = fix;
      pm.axis_dir = axis;
      pm.bend_dir = bend;
      pm.length = beam_length;

      if (v.size() >= 17) {
         double tip_amp2 = std::stod(v[14]);
         double period2  = std::stod(v[15]);
         double phase2   = std::stod(v[16]);
         double freq2 = (period2 > 0.) ? 1. / period2 : 0.;
         pm.amplitudes.push_back(tip_amp2);
         pm.freqs_hz.push_back(freq2);
         pm.phases_rad.push_back(phase2 - 2. * M_PI * freq2 * start_time);
      }
      motion_type_out = (pm.amplitudes.size() > 1) ? "cantilever(2mode)" : "cantilever";
      return pm;
   }
   std::cerr << "WARNING: unsupported velocity type '" << mtype
             << "' in body file; treating as fixed" << std::endl;
   return pm;
}

/* -------------------------------------------------------------------------- */
/*             Per-cable solver (new per-cable JSON format)                    */
/* -------------------------------------------------------------------------- */

// Helper: create getDouble/getInt/getArray3/getString lambdas from a JSON object
// with optional fallback to a parent JSON (for settings-mode defaults).
struct JSONReader {
   const JSON& input;
   const JSON* defaults;  // nullable — settings.json shared scalars

   double getDouble(const std::string& key, double def) const {
      if (input.find(key)) return std::stod(input.at(key)[0]);
      if (defaults && defaults->find(key)) return std::stod(defaults->at(key)[0]);
      return def;
   }
   int getInt(const std::string& key, int def) const {
      if (input.find(key)) return std::stoi(input.at(key)[0]);
      if (defaults && defaults->find(key)) return std::stoi(defaults->at(key)[0]);
      return def;
   }
   Tddd getArray3(const std::string& key, Tddd def) const {
      if (input.find(key)) {
         auto v = input.at(key);
         if (v.size() >= 3) return {std::stod(v[0]), std::stod(v[1]), std::stod(v[2])};
      }
      return def;
   }
   std::string getString(const std::string& key, const std::string& def) const {
      if (input.find(key)) return input.at(key)[0];
      if (defaults && defaults->find(key)) return defaults->at(key)[0];
      return def;
   }
   std::vector<double> getDoubleArray(const std::string& key) const {
      std::vector<double> out;
      const JSON* src = nullptr;
      if (input.find(key)) src = &input;
      else if (defaults && defaults->find(key)) src = defaults;
      if (!src) return out;
      for (const auto& s : src->at(key))
         out.push_back(std::stod(s));
      return out;
   }
};

// Resolve fluid preset and wind configuration from a JSONReader.
// Returns (rho, Cd, wind_fn) tuple. Wind is nullptr for "none" or missing.
//
// Settings.json flat-key schema:
//   "fluid":                       "water" | "air" (default "water")
//   "fluid_density":               [kg/m^3]  (optional, overrides preset)
//   "drag_Cd":                     [-]       (optional, overrides preset)
//   "wind_type":                   "none" | "uniform" | "AR1" (default "none")
//   "wind_U_mean":                 [Ux, Uy, Uz] [m/s]
//   "wind_turbulence_intensity":   sigma_u / |U_mean|  (AR1 only)
//   "wind_integral_time_scale":    T_L [s]             (AR1 only)
//   "wind_seed":                   RNG seed            (optional)
struct FluidConfig {
   double rho;
   double Cd;
   std::function<Tddd(const Tddd&, double)> wind_fn;
};

static FluidConfig resolveFluidConfig(const JSONReader& r) {
   FluidConfig cfg;
   std::string preset = r.getString("fluid", "water");
   if (preset == "air") {
      cfg.rho = 1.225;
      cfg.Cd = 1.2;
   } else {  // "water" or unknown -> water defaults (legacy behavior)
      cfg.rho = _WATER_DENSITY_;
      cfg.Cd = 2.5;
   }
   cfg.rho = r.getDouble("fluid_density", cfg.rho);
   cfg.Cd = r.getDouble("drag_Cd", cfg.Cd);

   std::string wtype = r.getString("wind_type", "none");
   if (wtype == "uniform") {
      Tddd U = r.getArray3("wind_U_mean", {0., 0., 0.});
      cfg.wind_fn = WindField::makeUniform(U);
   } else if (wtype == "AR1") {
      Tddd U = r.getArray3("wind_U_mean", {0., 0., 0.});
      double TI = r.getDouble("wind_turbulence_intensity", 0.15);
      double TL = r.getDouble("wind_integral_time_scale", 5.0);
      unsigned seed = static_cast<unsigned>(
          r.getInt("wind_seed", static_cast<int>(std::time(nullptr))));
      cfg.wind_fn = WindField::makeAR1(U, TI, TL, seed);
   }
   return cfg;
}

static void applyFluidConfig(LumpedCableSystem& sys, const FluidConfig& cfg) {
   for (auto* c : sys.cables()) {
      c->FluidDensity = cfg.rho;
      c->DragForceCoefficient = cfg.Cd;
      c->wind_field = cfg.wind_fn;
   }
   if (cfg.wind_fn)
      std::cout << "[fluid] rho=" << cfg.rho << " Cd=" << cfg.Cd << " wind=ON" << std::endl;
   else
      std::cout << "[fluid] rho=" << cfg.rho << " Cd=" << cfg.Cd << " wind=OFF" << std::endl;
}

static int solvePerCable(const fs::path& input_path,
                         const fs::path& output_dir,
                         const JSON* settings_defaults,
                         PVDWriter* master_pvd = nullptr,
                         int master_part = 0) {
   JSON input(input_path.string());
   JSONReader r{input, settings_defaults};

   std::string name = r.getString("name", input_path.stem().string());
   Tddd pos_a = r.getArray3("end_a_position", {0., 0., 0.});
   Tddd pos_b = r.getArray3("end_b_position", {0., 0., 0.});
   double cable_length = r.getDouble("cable_length", 100.);
   int n_points = r.getInt("n_points", 41);
   double line_density = r.getDouble("line_density", 348.5);
   double EA = r.getDouble("EA", 14e8);
   double damping = r.getDouble("damping", 0.5);
   double diameter = r.getDouble("diameter", 0.132);
   double gravity = r.getDouble("gravity", 9.81);
   std::string mode = r.getString("mode", "equilibrium");
   int max_steps = r.getInt("max_equilibrium_steps", 500000);
   double tol = r.getDouble("equilibrium_tol", 0.01);
   int snapshot_interval = r.getInt("snapshot_interval", 10000);

   // --- Initial condition: "length" (default) or "tension" ---
   std::string initial_condition = r.getString("initial_condition", "length");
   double target_tension_top = r.getDouble("tension_top", 0.);
   double target_tension_bot = r.getDouble("tension_bottom", 0.);
   // backward compat: single "tension" key sets both
   double target_tension_single = r.getDouble("tension", 0.);
   if (target_tension_top == 0. && target_tension_bot == 0. && target_tension_single > 0.) {
      target_tension_top = target_tension_single;
      target_tension_bot = target_tension_single;
   }
   double target_tension_avg = (target_tension_top + target_tension_bot) / 2.;

   _GRAVITY_ = gravity;
   _GRAVITY3_ = {0., 0., -gravity};

   double chord_length = Norm(pos_b - pos_a);
   double natural_length = cable_length;  // default: cable_length IS the natural length

   if (initial_condition == "tension" && target_tension_avg > 0.) {
      // Iterative solver: find natural_length such that the equilibrium
      // average tension matches target. Uses secant method.
      std::cout << "=== Cable Solver (per-cable, tension-driven) ===" << std::endl;
      std::cout << "Name: " << name << std::endl;
      std::cout << "End A: " << pos_a << "  End B: " << pos_b << std::endl;
      std::cout << "Chord length: " << chord_length << " m" << std::endl;
      std::cout << "Target tension_top: " << target_tension_top / 1e3 << " kN"
                << "  tension_bottom: " << target_tension_bot / 1e3 << " kN"
                << "  avg: " << target_tension_avg / 1e3 << " kN" << std::endl;

      // Initial guess from linear approximation
      double L0_a = chord_length / (1. + target_tension_avg / EA);
      double L0_b = L0_a * 0.999;  // slightly different for secant

      // Returns {T_top, T_bot} for a given natural length L0.
      auto solveAndGetTensions = [&](double L0) -> std::pair<double, double> {
         LumpedCableSystem sys;
         sys.addCable(name,
                      CableAttachment::worldFixed(pos_a),
                      CableAttachment::worldFixed(pos_b),
                      L0, n_points,
                      CableProperties{line_density, EA, damping, diameter});
         sys.solveEquilibrium(tol, max_steps, 0);
         auto* cable = sys.cables().front();
         auto points = getOrderedPoints(cable);
         int n = points.size();
         double nat_len = cable->natural_length();
         double stiffness_val = 0;
         for (auto l : cable->getLines()) { stiffness_val = l->stiffness; break; }
         if (n < 2) return {0., 0.};
         Tddd dv_top = points[n - 1]->X - points[n - 2]->X;
         double strain_top = (Norm(dv_top) - nat_len) / nat_len;
         double T_top = (strain_top > 0.) ? stiffness_val * strain_top : 0.;
         Tddd dv_bot = points[1]->X - points[0]->X;
         double strain_bot = (Norm(dv_bot) - nat_len) / nat_len;
         double T_bot = (strain_bot > 0.) ? stiffness_val * strain_bot : 0.;
         return {T_top, T_bot};
      };

      auto computeRMS = [&](double T_top, double T_bot) -> double {
         double e_top = T_top - target_tension_top;
         double e_bot = T_bot - target_tension_bot;
         return std::sqrt((e_top * e_top + e_bot * e_bot) / 2.);
      };

      // Secant method on the average tension (which controls the overall
      // tension level). The top/bottom split is determined by physics.
      // RMS of (top, bottom) errors is evaluated as quality metric.
      auto [Ttop_a, Tbot_a] = solveAndGetTensions(L0_a);
      auto [Ttop_b, Tbot_b] = solveAndGetTensions(L0_b);
      double Tavg_a = (Ttop_a + Tbot_a) / 2.;
      double Tavg_b = (Ttop_b + Tbot_b) / 2.;

      std::cout << "  iter 0: L0=" << L0_a
                << " T_top=" << Ttop_a / 1e3 << " T_bot=" << Tbot_a / 1e3
                << " RMS=" << computeRMS(Ttop_a, Tbot_a) / 1e3 << " kN" << std::endl;
      std::cout << "  iter 1: L0=" << L0_b
                << " T_top=" << Ttop_b / 1e3 << " T_bot=" << Tbot_b / 1e3
                << " RMS=" << computeRMS(Ttop_b, Tbot_b) / 1e3 << " kN" << std::endl;

      double tension_tol = target_tension_avg * 1e-4;
      int max_iter = 20;
      for (int it = 2; it < max_iter; ++it) {
         double r_a = Tavg_a - target_tension_avg;
         double r_b = Tavg_b - target_tension_avg;
         if (std::abs(r_b) < tension_tol) { L0_a = L0_b; break; }
         if (std::abs(r_a - r_b) < 1e-30) break;
         double L0_new = L0_b - r_b * (L0_b - L0_a) / (r_b - r_a);
         L0_new = std::max(L0_new, chord_length * 0.5);
         L0_new = std::min(L0_new, chord_length * 1.5);
         L0_a = L0_b; Tavg_a = Tavg_b;
         L0_b = L0_new;
         auto [tt, tb] = solveAndGetTensions(L0_b);
         Ttop_b = tt; Tbot_b = tb;
         Tavg_b = (Ttop_b + Tbot_b) / 2.;
         double rms = computeRMS(Ttop_b, Tbot_b);
         std::cout << "  iter " << it << ": L0=" << L0_b
                   << " T_top=" << Ttop_b / 1e3 << " T_bot=" << Tbot_b / 1e3
                   << " RMS=" << rms / 1e3 << " kN"
                   << " (" << rms / target_tension_avg * 100 << "%)" << std::endl;
         if (std::abs(Tavg_b - target_tension_avg) < tension_tol) break;
      }

      natural_length = L0_b;
      double final_rms = computeRMS(Ttop_b, Tbot_b);
      std::cout << "Converged natural length: " << natural_length << " m" << std::endl;
      std::cout << "Final: T_top=" << Ttop_b / 1e3 << " kN (target "
                << target_tension_top / 1e3 << ")  T_bot=" << Tbot_b / 1e3
                << " kN (target " << target_tension_bot / 1e3 << ")" << std::endl;
      std::cout << "RMS error: " << final_rms / 1e3 << " kN ("
                << final_rms / target_tension_avg * 100 << "%)" << std::endl;
      std::cout << "====================" << std::endl;

   } else {
      std::cout << "=== Cable Solver (per-cable) ===" << std::endl;
      std::cout << "Name: " << name << std::endl;
      std::cout << "End A: " << pos_a << "  End B: " << pos_b << std::endl;
      std::cout << "Natural length: " << natural_length << " m" << std::endl;
      std::cout << "n_points: " << n_points << std::endl;
      std::cout << "EA: " << EA << " N  density: " << line_density << " kg/m" << std::endl;
      std::cout << "Gravity: " << gravity << " m/s^2  Mode: " << mode << std::endl;
      std::cout << "====================" << std::endl;
   }

   if (mode != "equilibrium" && mode != "dynamic") {
      std::cerr << "Unknown mode: " << mode << std::endl;
      return 1;
   }

   // --- Parse endpoint motion ---
   auto parseMotion = [&](const std::string& prefix) -> PrescribedMotion {
      PrescribedMotion m;
      std::string mt = r.getString(prefix + "_motion", "fixed");
      if (mt == "sinusoidal") {
         m.kind = PrescribedMotion::Sinusoidal;
         m.dof = parseDof(r.getString(prefix + "_motion_dof", "heave"));
         m.amplitudes  = {r.getDouble(prefix + "_motion_amplitude", 0.)};
         m.freqs_hz    = {r.getDouble(prefix + "_motion_frequency", 0.)};
         m.phases_rad  = {r.getDouble(prefix + "_motion_phase", 0.)};
      } else if (mt == "sinusoidal_multi") {
         m.kind = PrescribedMotion::Sinusoidal;
         m.dof = parseDof(r.getString(prefix + "_motion_dof", "heave"));
         m.amplitudes  = r.getDoubleArray(prefix + "_motion_amplitudes");
         m.freqs_hz    = r.getDoubleArray(prefix + "_motion_frequencies");
         m.phases_rad  = r.getDoubleArray(prefix + "_motion_phases");
      }
      return m;
   };

   PrescribedMotion motion_a = parseMotion("end_a");
   PrescribedMotion motion_b = parseMotion("end_b");

   // --- Build cable system with appropriate attachments ---
   // For dynamic endpoints, create empty Network objects as bodies.
   std::unique_ptr<Network> body_a, body_b;

   auto makeEndpoint = [](const Tddd& pos, const PrescribedMotion& motion,
                          std::unique_ptr<Network>& body_out,
                          const std::string& body_name) -> CableAttachment {
      if (motion.kind == PrescribedMotion::Fixed) {
         return CableAttachment::worldFixed(pos);
      }
      body_out = std::make_unique<Network>("file_name_is_not_given", body_name);
      body_out->isRigidBody = true;
      body_out->isSoftBody = false;
      body_out->isFloatingBody = false;
      body_out->isAbsorber = false;
      body_out->COM = body_out->ICOM = pos;
      // Initialize RK accumulators so currentWorldPosition() works
      // before the dynamic loop starts (needed during solveEquilibrium).
      body_out->RK_COM.initialize(1./*dummy dt*/, 0., pos, 4);
      T4d identity_q = {1., 0., 0., 0.};
      body_out->RK_Q.initialize(1., 0., identity_q, 4);
      return CableAttachment::onBody(body_out.get(), pos);
   };

   CableAttachment att_a = makeEndpoint(pos_a, motion_a, body_a, "end_a");
   CableAttachment att_b = makeEndpoint(pos_b, motion_b, body_b, "end_b");

   LumpedCableSystem system;
   system.addCable(name, att_a, att_b,
                   natural_length, n_points,
                   CableProperties{line_density, EA, damping, diameter});

   // Apply fluid/wind settings (water+no-wind by default, preserves legacy behavior).
   // During solveEquilibrium() the LumpedCable RAII guards override these
   // temporarily (water density, Cd=1000, no wind), so equilibrium converges
   // identically regardless of target fluid.
   applyFluidConfig(system, resolveFluidConfig(r));

   auto t_start = std::chrono::high_resolution_clock::now();

   if (mode == "equilibrium") {
      // ======================== EQUILIBRIUM MODE ========================

      auto snapshot_cb = [&](int iter, double max_vel,
                             const std::map<std::string, std::vector<std::array<double, 3>>>& positions) {
         auto it = positions.begin();
         if (it == positions.end()) return;
         const auto& pts = it->second;
         std::cout << "SNAPSHOT {\"iter\":" << iter
                   << ",\"norm_v\":" << max_vel
                   << ",\"positions\":[";
         for (size_t k = 0; k < pts.size(); ++k) {
            std::cout << "[" << pts[k][0] << "," << pts[k][1] << "," << pts[k][2] << "]";
            if (k + 1 < pts.size()) std::cout << ",";
         }
         std::cout << "]}" << std::endl;
      };

      std::cout << "Finding equilibrium..." << std::endl;
      bool converged = system.solveEquilibrium(tol, max_steps, snapshot_interval, snapshot_cb);

      if (converged)
         std::cout << "Converged." << std::endl;
      else
         std::cout << "WARNING: did not converge after " << max_steps << " steps" << std::endl;

      auto t_end_eq = std::chrono::high_resolution_clock::now();
      double elapsed_ms = std::chrono::duration<double, std::milli>(t_end_eq - t_start).count();
      std::cout << "Computation time: " << elapsed_ms << " ms" << std::endl;

      fs::path result_path = output_dir / (name + "_result.json");
      writeSingleResultJSON(result_path, system.cables().front(), elapsed_ms, converged);

   } else {
      // ======================== DYNAMIC MODE ========================

      double dt = r.getDouble("dt", 0.01);
      double t_end_time = r.getDouble("t_end", 1.0);
      double output_interval = r.getDouble("output_interval", dt);
      int rk_order = 4;

      std::cout << "Dynamic mode: dt=" << dt << " t_end=" << t_end_time
                << " output_interval=" << output_interval << std::endl;

      // First solve static equilibrium at initial positions. We suppress the
      // wind field here so the initial condition is the gravity-only
      // catenary (fast pseudo-relaxation via Cd=1000 RAII in solveEquilibrium).
      // The dynamic loop below re-applies the user's wind_field.
      std::cout << "Initial equilibrium (wind suppressed for catenary IC)..." << std::endl;
      {
         std::vector<std::unique_ptr<LumpedCable::WindFieldGuard>> tmp_wind_guards;
         tmp_wind_guards.reserve(system.cables().size());
         for (auto* c : system.cables())
            tmp_wind_guards.emplace_back(
                std::make_unique<LumpedCable::WindFieldGuard>(c, nullptr));
         system.solveEquilibrium(tol, max_steps, 0 /* no snapshots */);
      }
      std::cout << "Initial equilibrium done." << std::endl;

      // Collect bodies that need RK driving
      std::vector<std::pair<Network*, PrescribedMotion*>> driven_bodies;
      if (body_a && motion_a.kind != PrescribedMotion::Fixed)
         driven_bodies.push_back({body_a.get(), &motion_a});
      if (body_b && motion_b.kind != PrescribedMotion::Fixed)
         driven_bodies.push_back({body_b.get(), &motion_b});

      // Compute CFL dt for the cable (same as solveEquilibrium).
      // The user-specified dt is the output/body time step. Internally
      // we sub-step at the cable's CFL-limited dt.
      auto* cable = system.cables().front();
      double wave_speed = std::sqrt(EA / line_density);
      double dt_cfl = cable->natural_length() / wave_speed;
      int substeps_per_dt = std::max(1, (int)std::ceil(dt / dt_cfl));
      double dt_sub = dt / substeps_per_dt;

      std::cout << "CFL: dt_cfl=" << dt_cfl << " substeps_per_dt=" << substeps_per_dt
                << " dt_sub=" << dt_sub << std::endl;

      // Time series storage
      std::vector<double> time_series;
      std::vector<double> top_tension_series;
      std::vector<double> bot_tension_series;
      std::vector<double> max_tension_series;

      double next_output_time = 0.;
      int step_count = 0;

      // ParaView (.vtp + .pvd) output. Enabled by default; can be disabled
      // via `output_paraview: false` in the input JSON. Per-cable .pvd lives
      // at <output_dir>/<name>.pvd; .vtp files go under
      // <output_dir>/paraview/<name>_<frame>.vtp.
      //
      // `paraview_interval` controls VTP output cadence in simulation
      // seconds. It is independent of `output_interval` (which governs the
      // JSON tension time series). Defaults to 10 × output_interval so a
      // typical run produces ~10 VTP frames per cable — enough for a
      // recognizable animation, small enough to avoid filesystem bloat.
      bool output_paraview = (r.getString("output_paraview", "true") != "false");
      double paraview_interval = r.getDouble("paraview_interval", 10.0 * output_interval);
      if (paraview_interval <= 0.)
         paraview_interval = 10.0 * output_interval;
      double next_paraview_time = 0.;
      std::optional<PVDWriter> per_cable_pvd;
      if (output_paraview) {
         fs::create_directories(output_dir / "paraview");
         per_cable_pvd.emplace((output_dir / (name + ".pvd")).string());
      }
      int snap_index = 0;
      // PVD flush cadence: rewrite the .pvd every N frames for crash
      // resilience, plus once at simulation end. Per-frame flushes are
      // O(N²) cumulative I/O and needless.
      constexpr int PVD_FLUSH_EVERY = 25;

      // Helper to compute current tensions and per-segment arrays. Returns
      // (top, bottom, max) and fills out `seg_t` and `seg_eps` with the
      // per-segment tension [N] and strain [-] (tension-only clamped so
      // negative strains map to zero tension; strain itself is kept signed
      // for diagnostics).
      auto computeTensions = [&](std::vector<double>& seg_t,
                                 std::vector<double>& seg_eps)
          -> std::tuple<double, double, double> {
         auto points = getOrderedPoints(cable);
         int n = points.size();
         double nat_len = cable->natural_length();
         double stiffness_val = 0;
         for (auto l : cable->getLines()) { stiffness_val = l->stiffness; break; }
         seg_t.assign(n > 1 ? n - 1 : 0, 0.);
         seg_eps.assign(n > 1 ? n - 1 : 0, 0.);
         for (int i = 0; i < n - 1; ++i) {
            Tddd dv = points[i + 1]->X - points[i]->X;
            double dist = Norm(dv);
            double strain = (dist - nat_len) / nat_len;
            seg_eps[i] = strain;
            if (strain > 0.) seg_t[i] = stiffness_val * strain;
         }
         double top = (n > 1) ? seg_t.back() : 0;
         double bot = (n > 1) ? seg_t.front() : 0;
         double mx = (n > 1) ? *std::max_element(seg_t.begin(), seg_t.end()) : 0;
         return {top, bot, mx};
      };

      // Write one ParaView frame: .vtp + PVD push (per-cable and master).
      // The PVD files are flushed every PVD_FLUSH_EVERY frames; a final
      // flush happens after the time loop (see below).
      auto writeParaviewFrame = [&](double t_frame,
                                    const std::vector<networkPoint*>& points,
                                    const std::vector<double>& seg_t,
                                    const std::vector<double>& seg_eps) {
         if (!output_paraview) return;
         fs::path vtp_rel = fs::path("paraview") / (name + "_" + zfill_int(snap_index) + ".vtp");
         writeCableVTP(output_dir / vtp_rel, points, seg_t, seg_eps);
         per_cable_pvd->push(vtp_rel.string(), t_frame);
         if (master_pvd) {
            // master_pvd is shared when solvePerCable is called concurrently
            // from the settings-mode dynamic loop; serialize the append/flush.
            _Pragma("omp critical(master_pvd)")
            {
               master_pvd->push(vtp_rel.string(), t_frame, master_part);
            }
         }
         if (snap_index % PVD_FLUSH_EVERY == 0) {
            per_cable_pvd->output();
            if (master_pvd) {
               _Pragma("omp critical(master_pvd)")
               {
                  master_pvd->output();
               }
            }
         }
         ++snap_index;
      };

      // Record initial state
      std::vector<double> seg_t_cur, seg_eps_cur;
      {
         auto [top, bot, mx] = computeTensions(seg_t_cur, seg_eps_cur);
         time_series.push_back(0.);
         top_tension_series.push_back(top);
         bot_tension_series.push_back(bot);
         max_tension_series.push_back(mx);
         next_output_time = output_interval;
      }

      // Emit a t=0 snapshot of the initial equilibrium state so playback
      // starts cleanly from the static catenary.
      if (snapshot_interval > 0) {
         auto points = getOrderedPoints(cable);
         std::cout << "SNAPSHOT {\"cable\":\"" << name << "\""
                   << ",\"t\":0,\"iter\":-1"
                   << ",\"positions\":[";
         for (size_t k = 0; k < points.size(); ++k) {
            auto X = points[k]->X;
            std::cout << "[" << std::get<0>(X) << "," << std::get<1>(X)
                      << "," << std::get<2>(X) << "]";
            if (k + 1 < points.size()) std::cout << ",";
         }
         std::cout << "]}" << std::endl;
      }
      if (output_paraview) {
         auto points = getOrderedPoints(cable);
         writeParaviewFrame(0., points, seg_t_cur, seg_eps_cur);
         next_paraview_time = paraview_interval;
      }

      // Dynamic time loop with CFL sub-stepping.
      // Outer loop: user dt (for output and body motion).
      // Inner loop: CFL-limited sub-steps for cable stability.
      for (double t = 0.; t < t_end_time - dt * 0.5; t += dt, step_count++) {
         for (int sub = 0; sub < substeps_per_dt; ++sub) {
            double t_sub = t + sub * dt_sub;

            // Update driven body positions to analytical pose at t_sub + dt_sub
            for (auto& [body, motion] : driven_bodies) {
               body->COM = body->ICOM + motion->displacement(t_sub + dt_sub);
            }

            // Single step of cable dynamics using step() directly
            // (advanceRKStage internally calls cable->step which does a full
            // internal RK4 at the given dt_sub, respecting CFL stability).
            auto setBCfixed = [&](networkPoint* p) {
               if (p == cable->firstPoint) {
                  auto target = att_a.kind == CableAttachment::WorldFixed
                     ? att_a.world_position : body_a->COM;
                  Tddd v = (target - p->X) / dt_sub;
                  p->acceleration.fill(0);
                  std::get<0>(p->velocity) = std::get<0>(v);
                  std::get<1>(p->velocity) = std::get<1>(v);
                  std::get<2>(p->velocity) = std::get<2>(v);
               } else if (p == cable->lastPoint) {
                  auto target = att_b.kind == CableAttachment::WorldFixed
                     ? att_b.world_position : body_b->COM;
                  Tddd v = (target - p->X) / dt_sub;
                  p->acceleration.fill(0);
                  std::get<0>(p->velocity) = std::get<0>(v);
                  std::get<1>(p->velocity) = std::get<1>(v);
                  std::get<2>(p->velocity) = std::get<2>(v);
               }
            };

            cable->step(t_sub, dt_sub, setBCfixed);

            // Apply RK result
            cable->applyMooringSimulationResult();

            // Snap endpoints to exact positions
            if (att_a.kind == CableAttachment::WorldFixed)
               cable->firstPoint->setX(att_a.world_position);
            else if (body_a)
               cable->firstPoint->setX(body_a->COM);

            if (att_b.kind == CableAttachment::WorldFixed)
               cable->lastPoint->setX(att_b.world_position);
            else if (body_b)
               cable->lastPoint->setX(body_b->COM);
         }

         // Output at user dt intervals
         double t_now = t + dt;
         bool did_output_step = false;
         if (t_now >= next_output_time - dt * 0.01) {
            auto [top, bot, mx] = computeTensions(seg_t_cur, seg_eps_cur);
            time_series.push_back(t_now);
            top_tension_series.push_back(top);
            bot_tension_series.push_back(bot);
            max_tension_series.push_back(mx);
            next_output_time += output_interval;
            did_output_step = true;
         }

         // ParaView frame at its own (coarser) cadence. Independent of the
         // JSON time series so a 10-frame animation can coexist with a
         // 100-sample tension curve.
         if (output_paraview && t_now >= next_paraview_time - dt * 0.01) {
            if (!did_output_step)
               (void)computeTensions(seg_t_cur, seg_eps_cur);
            auto points = getOrderedPoints(cable);
            writeParaviewFrame(t_now, points, seg_t_cur, seg_eps_cur);
            next_paraview_time += paraview_interval;
         }

         // SNAPSHOT (stdout, for GUI memory playback)
         if (snapshot_interval > 0 && step_count % snapshot_interval == 0) {
            auto points = getOrderedPoints(cable);
            std::cout << "SNAPSHOT {\"cable\":\"" << name << "\""
                      << ",\"t\":" << t_now << ",\"iter\":" << step_count
                      << ",\"positions\":[";
            for (size_t k = 0; k < points.size(); ++k) {
               auto X = points[k]->X;
               std::cout << "[" << std::get<0>(X) << "," << std::get<1>(X)
                         << "," << std::get<2>(X) << "]";
               if (k + 1 < points.size()) std::cout << ",";
            }
            std::cout << "]}" << std::endl;
         }
         (void)did_output_step;  // reserved for future use
      }

      auto t_end_dyn = std::chrono::high_resolution_clock::now();
      double elapsed_ms = std::chrono::duration<double, std::milli>(t_end_dyn - t_start).count();
      std::cout << "Dynamic simulation complete. Steps: " << step_count
                << "  Elapsed: " << elapsed_ms << " ms" << std::endl;

      // Final flush of the per-cable PVD so the last batch of frames is
      // visible in ParaView even if PVD_FLUSH_EVERY didn't hit at the end.
      if (output_paraview && per_cable_pvd) {
         per_cable_pvd->output();
         std::cout << "ParaView: " << snap_index << " frame(s) written "
                   << "(interval=" << paraview_interval << "s)" << std::endl;
      }

      // Write dynamic result JSON
      {
         fs::path result_path = output_dir / (name + "_result.json");
         std::ofstream ofs(result_path);
         ofs << std::setprecision(10);
         ofs << "{\n";
         ofs << "  \"name\": \"" << name << "\",\n";
         ofs << "  \"mode\": \"dynamic\",\n";
         ofs << "  \"dt\": " << dt << ",\n";
         ofs << "  \"t_end\": " << t_end_time << ",\n";
         ofs << "  \"n_output_steps\": " << time_series.size() << ",\n";

         ofs << "  \"time\": [";
         for (size_t i = 0; i < time_series.size(); ++i)
            ofs << time_series[i] << (i + 1 < time_series.size() ? "," : "");
         ofs << "],\n";

         ofs << "  \"top_tension\": [";
         for (size_t i = 0; i < top_tension_series.size(); ++i)
            ofs << top_tension_series[i] << (i + 1 < top_tension_series.size() ? "," : "");
         ofs << "],\n";

         ofs << "  \"bottom_tension\": [";
         for (size_t i = 0; i < bot_tension_series.size(); ++i)
            ofs << bot_tension_series[i] << (i + 1 < bot_tension_series.size() ? "," : "");
         ofs << "],\n";

         ofs << "  \"max_tension\": [";
         for (size_t i = 0; i < max_tension_series.size(); ++i)
            ofs << max_tension_series[i] << (i + 1 < max_tension_series.size() ? "," : "");
         ofs << "],\n";

         // Final cable shape
         auto points = getOrderedPoints(system.cables().front());
         ofs << "  \"positions_final\": [";
         for (size_t i = 0; i < points.size(); ++i) {
            auto X = points[i]->X;
            ofs << "[" << std::get<0>(X) << "," << std::get<1>(X) << "," << std::get<2>(X) << "]";
            if (i + 1 < points.size()) ofs << ",";
         }
         ofs << "],\n";

         ofs << "  \"computation_time_ms\": " << elapsed_ms << "\n";
         ofs << "}\n";
         ofs.close();
         std::cout << "Dynamic result written to " << result_path << std::endl;
      }
   }

   return 0;
}

/* -------------------------------------------------------------------------- */
/*                                  main                                      */
/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[]) {
   if (argc < 3) {
      std::cerr << "Usage: " << argv[0] << " input.json output_dir/" << std::endl;
      return 1;
   }

   fs::path input_path(argv[1]);
   fs::path output_dir(argv[2]);

   if (!fs::exists(input_path)) {
      std::cerr << "ERROR: input file not found: " << input_path << std::endl;
      return 1;
   }
   fs::create_directories(output_dir);

   /* ---------------------- read input JSON ---------------------- */

   JSON input(input_path.string());

   /* --- FORMAT 1: Settings mode (has `input_files` key) --- */
   if (input.find("input_files")) {
      fs::path input_dir = input_path.parent_path();
      auto cable_files = input.at("input_files");

      JSONReader sr{input, nullptr};
      double gravity = sr.getDouble("gravity", 9.81);
      std::string mode = sr.getString("mode", "equilibrium");
      int max_steps = sr.getInt("max_equilibrium_steps", 500000);
      double tol = sr.getDouble("equilibrium_tol", 0.01);
      int snapshot_interval = sr.getInt("snapshot_interval", 10000);

      _GRAVITY_ = gravity;
      _GRAVITY3_ = {0., 0., -gravity};

      // ---- Parse bodies section ----
      struct BodyDef {
         std::string name;
         std::unique_ptr<Network> network;
         PrescribedMotion motion;
      };
      std::map<std::string, std::shared_ptr<BodyDef>> bodies;

      // Pre-scan input_files: discriminate between cable files and body files.
      // A body file is a JSON with `type` field containing "RigidBody" (BEM
      // convention). Each body file becomes a BodyDef in the `bodies` map.
      // All other files are treated as cables.
      std::vector<std::string> actual_cable_files;
      for (const auto& fname : cable_files) {
         fs::path p = input_dir / fname;
         if (!fs::exists(p)) {
            // keep the missing file in cable list so the per-cable loop reports it
            actual_cable_files.push_back(fname);
            continue;
         }
         JSON j(p.string());
         JSONReader jr{j, nullptr};
         std::string type = jr.getString("type", "");
         if (type.find("RigidBody") != std::string::npos) {
            std::string bname = jr.getString("name", fs::path(fname).stem().string());
            std::string motion_type;
            PrescribedMotion pm = parseBodyVelocity(j, motion_type);
            auto bd = std::make_shared<BodyDef>();
            bd->name = bname;
            bd->motion = pm;
            bodies[bname] = bd;
            std::cout << "Body: " << bname << "  motion=" << motion_type
                      << "  (from " << fname << ")" << std::endl;
         } else {
            actual_cable_files.push_back(fname);
         }
      }
      cable_files = actual_cable_files;

      // ---- If no bodies defined or equilibrium: solve each cable independently ----
      if (bodies.empty() || mode == "equilibrium") {
         std::cout << "=== Cable Solver (settings mode, " << cable_files.size()
                   << " cable files, equilibrium) ===" << std::endl;
         int failures = 0;
         // Parallel across independent cable files. Each solvePerCable creates
         // its own LumpedCableSystem (private RK state, private output file)
         // so they do not share any per-cable state. The only shared writes
         // are the global _GRAVITY_ / _GRAVITY3_ scalars, which all cables
         // set to the same value from settings.json — a benign race.
         // stdout snapshot lines may interleave across threads; each line is
         // self-contained JSON with a cable name so downstream consumers can
         // demux. For clean logs set OMP_NUM_THREADS=1.
         _Pragma("omp parallel for reduction(+:failures) schedule(dynamic)")
         for (size_t idx = 0; idx < cable_files.size(); ++idx) {
            const auto& fname = cable_files[idx];
            fs::path cable_path = input_dir / fname;
            if (!fs::exists(cable_path)) {
               std::cerr << "ERROR: cable file not found: " << cable_path << std::endl;
               failures++;
               continue;
            }
            std::cout << "\n--- " << fname << " ---" << std::endl;
            int ret = solvePerCable(cable_path, output_dir, &input);
            if (ret != 0) failures++;
         }
         std::cout << "\n=== Settings mode complete: " << cable_files.size()
                   << " cables, " << failures << " failures ===" << std::endl;
         return failures > 0 ? 1 : 0;
      }

      // ---- Dynamic mode with bodies: BEM-pattern shared time loop ----
      std::cout << "=== Cable Solver (settings mode, " << cable_files.size()
                << " cables, dynamic with " << bodies.size() << " bodies) ===" << std::endl;

      double dt = sr.getDouble("dt", 0.01);
      double t_end_time = sr.getDouble("t_end", 1.0);
      double output_interval_val = sr.getDouble("output_interval", dt);

      // First pass: solve equilibrium for each cable independently to get
      // initial state, then collect them into a shared system for dynamic.
      // (Tension initial condition is handled inside solvePerCable.)

      // For now, run equilibrium per cable, then dynamic per cable with body motion.
      // Full multi-cable shared time loop is a future enhancement.
      // Master PVD aggregates all cables into one ParaView file (one `part`
      // per cable). This is the single file to open in ParaView.
      PVDWriter master_pvd((output_dir / "cables.pvd").string());
      int failures = 0;
      // Parallelize across cables. master_pvd is shared — access is guarded
      // by `omp critical(master_pvd)` inside solvePerCable (writeParaviewFrame).
      // `cable_part` is derived from loop index so no counter race. stdout may
      // interleave across threads; each line carries a cable name for demux.
      _Pragma("omp parallel for reduction(+:failures) schedule(dynamic)")
      for (size_t idx = 0; idx < cable_files.size(); ++idx) {
         const auto& fname = cable_files[idx];
         int cable_part = static_cast<int>(idx);
         fs::path cable_path = input_dir / fname;
         if (!fs::exists(cable_path)) {
            std::cerr << "ERROR: cable file not found: " << cable_path << std::endl;
            failures++;
            continue;
         }

         // Read the cable JSON to find body references and attach positions.
         JSON cable_json(cable_path.string());
         JSONReader cr{cable_json, &input};
         std::string end_a_body_name = cr.getString("end_a_body", "");
         std::string end_b_body_name = cr.getString("end_b_body", "");
         Tddd pos_a = cr.getArray3("end_a_position", {0., 0., 0.});
         Tddd pos_b = cr.getArray3("end_b_position", {0., 0., 0.});

         // If a body is referenced, inject its motion into the cable JSON
         // as end_X_motion keys so solvePerCable picks them up.
         // This bridges the body-based design to the per-cable solver.
         JSON augmented = cable_json;
         auto setKey = [&](const std::string& k, const std::string& v) {
            augmented.map_S_S[k] = std::vector<std::string>{v};
         };
         auto setKeyArray = [&](const std::string& k, const std::vector<double>& vals) {
            std::vector<std::string> sv;
            sv.reserve(vals.size());
            for (double v : vals) sv.push_back(std::to_string(v));
            augmented.map_S_S[k] = sv;
         };

         auto injectMotion = [&](const std::string& prefix,
                                 const std::string& body_name,
                                 const Tddd& attach_pos) {
            if (body_name.empty()) return;
            auto it = bodies.find(body_name);
            if (it == bodies.end()) return;
            auto& pm = it->second->motion;
            std::string dof_str = (pm.dof == 0) ? "surge" : (pm.dof == 1) ? "sway" : "heave";

            if (pm.kind == PrescribedMotion::Cantilever) {
               auto local_amps = cantileverLocalAmplitudes(pm, attach_pos);
               if (local_amps.size() <= 1) {
                  // Single-mode: backward compat sinusoidal
                  setKey(prefix + "_motion", "sinusoidal");
                  setKey(prefix + "_motion_dof", dof_str);
                  setKey(prefix + "_motion_amplitude", std::to_string(local_amps.empty() ? 0. : local_amps[0]));
                  setKey(prefix + "_motion_frequency", std::to_string(pm.freq_hz_0()));
                  setKey(prefix + "_motion_phase", std::to_string(pm.phase_rad_0()));
               } else {
                  // Multi-mode: emit sinusoidal_multi with array keys
                  setKey(prefix + "_motion", "sinusoidal_multi");
                  setKey(prefix + "_motion_dof", dof_str);
                  setKeyArray(prefix + "_motion_amplitudes", local_amps);
                  setKeyArray(prefix + "_motion_frequencies", pm.freqs_hz);
                  setKeyArray(prefix + "_motion_phases", pm.phases_rad);
               }
            } else if (pm.kind == PrescribedMotion::Sinusoidal) {
               setKey(prefix + "_motion", "sinusoidal");
               setKey(prefix + "_motion_dof", dof_str);
               setKey(prefix + "_motion_amplitude", std::to_string(pm.amplitude()));
               setKey(prefix + "_motion_frequency", std::to_string(pm.freq_hz_0()));
               setKey(prefix + "_motion_phase", std::to_string(pm.phase_rad_0()));
            }
         };
         injectMotion("end_a", end_a_body_name, pos_a);
         injectMotion("end_b", end_b_body_name, pos_b);

         // Also inject dynamic mode keys
         setKey("mode", "dynamic");
         setKey("dt", std::to_string(dt));
         setKey("t_end", std::to_string(t_end_time));
         setKey("output_interval", std::to_string(output_interval_val));

         // Write augmented JSON to temp file and solve
         fs::path tmp_cable = output_dir / (fs::path(fname).stem().string() + "_input.json");
         {
            std::ofstream ofs(tmp_cable);
            ofs << "{\n";
            bool first = true;
            for (auto& [k, v] : augmented()) {
               if (!first) ofs << ",\n";
               first = false;
               if (v.size() == 1) {
                  // Try to write as number if possible, otherwise as string
                  bool is_num = true;
                  try { std::stod(v[0]); } catch (...) { is_num = false; }
                  if (is_num && k != "name" && k != "mode" && k != "initial_condition"
                      && k.find("_body") == std::string::npos
                      && k.find("_motion") == std::string::npos
                      && k.find("_dof") == std::string::npos) {
                     ofs << "  \"" << k << "\": " << v[0];
                  } else {
                     ofs << "  \"" << k << "\": \"" << v[0] << "\"";
                  }
               } else {
                  ofs << "  \"" << k << "\": [";
                  for (size_t i = 0; i < v.size(); ++i) {
                     bool is_num = true;
                     try { std::stod(v[i]); } catch (...) { is_num = false; }
                     if (is_num) ofs << v[i]; else ofs << "\"" << v[i] << "\"";
                     if (i + 1 < v.size()) ofs << ", ";
                  }
                  ofs << "]";
               }
            }
            ofs << "\n}\n";
         }

         std::cout << "\n--- " << fname << " (body-driven dynamic) ---" << std::endl;
         int ret = solvePerCable(tmp_cable, output_dir, &input, &master_pvd, cable_part);
         if (ret != 0) failures++;
      }

      // Final flush of the master PVD (solvePerCable already flushes per
      // frame, but write once more so the file reflects the final state).
      master_pvd.output();
      std::cout << "Master ParaView file: " << (output_dir / "cables.pvd") << std::endl;

      std::cout << "\n=== Settings mode complete: " << cable_files.size()
                << " cables, " << failures << " failures ===" << std::endl;
      return failures > 0 ? 1 : 0;
   }

   /* --- FORMAT 2: Per-cable (has `end_a_position` key) --- */
   if (input.find("end_a_position")) {
      return solvePerCable(input_path, output_dir, nullptr);
   }

   /* --- FORMAT 3 & 4: Legacy multi-line / single-line (unchanged) --- */

   auto getDouble = [&](const std::string& key, double default_val) -> double {
      if (input.find(key))
         return std::stod(input.at(key)[0]);
      return default_val;
   };

   auto getInt = [&](const std::string& key, int default_val) -> int {
      if (input.find(key))
         return std::stoi(input.at(key)[0]);
      return default_val;
   };

   auto getArray3 = [&](const std::string& key, Tddd default_val) -> Tddd {
      if (input.find(key)) {
         auto v = input.at(key);
         if (v.size() >= 3)
            return {std::stod(v[0]), std::stod(v[1]), std::stod(v[2])};
      }
      return default_val;
   };

   auto getString = [&](const std::string& key, const std::string& default_val) -> std::string {
      if (input.find(key))
         return input.at(key)[0];
      return default_val;
   };

   // Top-level scalars (shared across all cables in multi-line mode)
   double gravity = getDouble("gravity", 9.81);
   std::string mode = getString("mode", "equilibrium");
   int max_steps = getInt("max_equilibrium_steps", 500000);
   double tol = getDouble("equilibrium_tol", 0.01);
   int snapshot_interval = getInt("snapshot_interval", 10000);

   _GRAVITY_ = gravity;
   _GRAVITY3_ = {0., 0., -gravity};

   // Detect input format: scan for any `mooring_*` or `cable_*` keys.
   // The 13-element flat array form is the multi-line / BEM-compatible mode.
   std::vector<std::pair<std::string, std::vector<std::string>>> multi_keys;
   input.for_each([&](const std::string& key, const std::vector<std::string>& value) {
      // BEM uses substring `"mooring"` (`key.contains("mooring")`); we accept
      // both prefixes here for symmetry with the project naming.
      bool is_mooring = key.find("mooring") != std::string::npos
                        && key.find("mooring_") != std::string::npos;
      bool is_cable = key.find("cable_") == 0;  // strict prefix
      if ((is_mooring || is_cable) && value.size() == 13) {
         multi_keys.emplace_back(key, value);
      }
   });

   const bool is_multi_mode = !multi_keys.empty();

   /* ---------------------- build cable system ---------------------- */

   LumpedCableSystem system;

   if (is_multi_mode) {
      std::cout << "=== Cable Solver (multi-line, " << multi_keys.size() << " cables) ===" << std::endl;
      std::cout << "Mode: " << mode << std::endl;
      std::cout << "Gravity: " << gravity << " m/s^2" << std::endl;
      std::cout << "====================" << std::endl;

      for (const auto& [key, value] : multi_keys) {
         const std::string name = value[0];
         const Tddd point_a = {std::stod(value[1]), std::stod(value[2]), std::stod(value[3])};
         const Tddd point_b = {std::stod(value[4]), std::stod(value[5]), std::stod(value[6])};
         const double cable_length = std::stod(value[7]);
         const int n_points = std::stoi(value[8]);
         const double line_density = std::stod(value[9]);
         const double EA = std::stod(value[10]);
         const double damping = std::stod(value[11]);
         const double diameter = std::stod(value[12]);

         std::cout << "  cable [" << name << "] " << point_a << " → " << point_b
                   << "  L=" << cable_length << "  n_points=" << n_points
                   << "  EA=" << EA << "  density=" << line_density << std::endl;

         system.addCable(name,
                         CableAttachment::worldFixed(point_a),
                         CableAttachment::worldFixed(point_b),
                         cable_length,
                         n_points,
                         CableProperties{line_density, EA, damping, diameter});
      }
   } else {
      // Single-line legacy format
      Tddd point_a = getArray3("point_a", {500., 0., -58.});
      Tddd point_b = getArray3("point_b", {0., 0., 0.});
      double cable_length = getDouble("cable_length", 522.);
      int n_segments = getInt("n_segments", 40);
      double line_density = getDouble("line_density", 348.5);
      double EA = getDouble("EA", 14e8);
      double damping = getDouble("damping", 0.5);
      double diameter = getDouble("diameter", 0.132);
      int n_points = n_segments + 1;

      std::cout << "=== Cable Solver ===" << std::endl;
      std::cout << "Point A: " << point_a << std::endl;
      std::cout << "Point B: " << point_b << std::endl;
      std::cout << "Cable length: " << cable_length << " m" << std::endl;
      std::cout << "Segments: " << n_segments << std::endl;
      std::cout << "Line density: " << line_density << " kg/m" << std::endl;
      std::cout << "EA: " << EA << " N" << std::endl;
      std::cout << "Damping: " << damping << std::endl;
      std::cout << "Diameter: " << diameter << " m" << std::endl;
      std::cout << "Gravity: " << gravity << " m/s^2" << std::endl;
      std::cout << "Mode: " << mode << std::endl;
      std::cout << "====================" << std::endl;

      system.addCable("cable",
                      CableAttachment::worldFixed(point_a),
                      CableAttachment::worldFixed(point_b),
                      cable_length,
                      n_points,
                      CableProperties{line_density, EA, damping, diameter});
   }

   // Apply fluid/wind settings to all cables built above. Missing keys =>
   // water density, Cd=2.5, no wind (legacy behavior).
   {
      JSONReader fluid_r{input, nullptr};
      applyFluidConfig(system, resolveFluidConfig(fluid_r));
   }

   if (mode != "equilibrium") {
      std::cerr << "Unknown mode: " << mode << std::endl;
      return 1;
   }

   /* ---------------------- solve ---------------------- */

   auto t_start = std::chrono::high_resolution_clock::now();

   // SNAPSHOT callback streamed to stdout. Single-line mode keeps the
   // legacy schema; multi-line mode adds a `cable` field per snapshot so
   // pycable can demux which line is being updated.
   auto snapshot_cb = [&](int iter, double max_vel,
                          const std::map<std::string, std::vector<std::array<double, 3>>>& positions) {
      if (is_multi_mode) {
         for (const auto& [name, pts] : positions) {
            std::cout << "SNAPSHOT {\"iter\":" << iter
                      << ",\"cable\":\"" << name << "\""
                      << ",\"norm_v\":" << max_vel
                      << ",\"positions\":[";
            for (size_t k = 0; k < pts.size(); ++k) {
               std::cout << "[" << pts[k][0] << "," << pts[k][1] << "," << pts[k][2] << "]";
               if (k + 1 < pts.size()) std::cout << ",";
            }
            std::cout << "]}" << std::endl;
         }
      } else {
         // Single-line: pick the only cable
         auto it = positions.begin();
         if (it == positions.end()) return;
         const auto& pts = it->second;
         std::cout << "SNAPSHOT {\"iter\":" << iter
                   << ",\"norm_v\":" << max_vel
                   << ",\"positions\":[";
         for (size_t k = 0; k < pts.size(); ++k) {
            std::cout << "[" << pts[k][0] << "," << pts[k][1] << "," << pts[k][2] << "]";
            if (k + 1 < pts.size()) std::cout << ",";
         }
         std::cout << "]}" << std::endl;
      }
   };

   std::cout << "Finding equilibrium..." << std::endl;
   bool converged = system.solveEquilibrium(tol, max_steps, snapshot_interval, snapshot_cb);

   if (converged)
      std::cout << "Converged." << std::endl;
   else
      std::cout << "WARNING: did not converge after " << max_steps << " steps" << std::endl;

   auto t_end = std::chrono::high_resolution_clock::now();
   double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
   std::cout << "Computation time: " << elapsed_ms << " ms" << std::endl;

   /* ---------------------- write result ---------------------- */

   if (is_multi_mode)
      writeMultiResultJSON(output_dir / "result.json", system, elapsed_ms, converged);
   else
      writeSingleResultJSON(output_dir / "result.json", system.cables().front(), elapsed_ms, converged);

   return 0;
}
