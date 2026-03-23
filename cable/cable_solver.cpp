/**
 * @file cable_solver.cpp
 * @brief Cable playground solver — JSON I/O wrapper around MooringLine.
 *
 * Usage:
 *   ./cable_solver input.json output_dir/
 *
 * Build:
 *   mkdir -p build && cd build
 *   cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=cable/cable_solver.cpp -DOUTPUT_NAME=cable_solver ../..
 *   make -j$(sysctl -n hw.logicalcpu)
 */

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "MooringLine.hpp"
#include "basic.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "basic_vectors.hpp"

namespace fs = std::filesystem;

/* -------------------------------------------------------------------------- */
/*                         JSON result writer                                 */
/* -------------------------------------------------------------------------- */

// Walk the mooring line from firstPoint to lastPoint in topological order.
std::vector<networkPoint*> getOrderedPoints(const MooringLine* mooring) {
   std::vector<networkPoint*> ordered;
   auto current = mooring->firstPoint;
   networkPoint* prev = nullptr;
   while (current) {
      ordered.push_back(current);
      if (current == mooring->lastPoint) break;
      // Find the next point by following the line that doesn't go back to prev
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

void writeResultJSON(const fs::path& filepath,
                     const MooringLine* mooring,
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
   // Get EA from any line element
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

   // Per-node tension: average of adjacent segments (for visualization)
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

   // Summary
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

   // Read parameters
   Tddd point_a = getArray3("point_a", {500., 0., -58.});
   Tddd point_b = getArray3("point_b", {0., 0., 0.});
   double cable_length = getDouble("cable_length", 522.);
   int n_segments = getInt("n_segments", 40);
   double line_density = getDouble("line_density", 348.5);
   double EA = getDouble("EA", 14e8);
   double damping = getDouble("damping", 0.5);
   double diameter = getDouble("diameter", 0.132);
   double gravity = getDouble("gravity", 9.81);
   std::string mode = getString("mode", "equilibrium");

   // Apply gravity
   _GRAVITY_ = gravity;
   _GRAVITY3_ = {0., 0., -gravity};

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

   /* ---------------------- create mooring line ---------------------- */

   auto mooring = new MooringLine(point_a, point_b, cable_length, n_points);
   mooring->setDensityStiffnessDampingDiameter(line_density, EA, damping, diameter);

   for (auto p : mooring->getPoints()) {
      p->velocity.fill(0);
      p->acceleration.fill(0);
   }

   /* ---------------------- boundary condition ---------------------- */

   auto boundary_condition = [&](networkPoint* p) {
      if (p == mooring->firstPoint) {
         p->acceleration.fill(0);
         p->velocity.fill(0);
      }
      if (p == mooring->lastPoint) {
         p->acceleration.fill(0);
         p->velocity.fill(0);
      }
   };

   /* ---------------------- solve ---------------------- */

   auto t_start = std::chrono::high_resolution_clock::now();
   bool converged = true;

   if (mode == "equilibrium") {
      std::cout << "Finding equilibrium..." << std::endl;

      // Direct RK4 equilibrium solver — bypasses simulate()'s sub-stepping
      // ramp which wastes ~50% of computation restarting from dt_cfl=1e-8
      // every call. Here we use a fixed CFL-based dt with no ramp.
      double saved_Cd = mooring->DragForceCoefficient;
      double eq_drag = 1000.;
      mooring->DragForceCoefficient = eq_drag;

      auto ordered = getOrderedPoints(mooring);
      int n_nodes = ordered.size();
      double natural_len = mooring->natural_length();
      double wave_speed = std::sqrt(EA / line_density);
      double dt_cfl = 1.0 * natural_len / wave_speed;  // CFL=1.0 (RK4 allows this)

      int max_steps = getInt("max_equilibrium_steps", 500000);
      double tol = getDouble("equilibrium_tol", 0.01);
      int snapshot_interval = getInt("snapshot_interval", 10000);

      std::cout << "dt_cfl=" << dt_cfl << " wave_speed=" << wave_speed
                << " max_steps=" << max_steps << std::endl;

      for (int step = 0; step < max_steps; ++step) {
         // Initialize RK4 for this step
         for (auto& p : ordered) {
            p->RK_velocity_sub.initialize(dt_cfl, 0, p->velocityTranslational(), 4);
            p->RK_X_sub.initialize(dt_cfl, 0, p->X, 4);
         }

         // 4 RK4 stages
         while (true) {
            for (auto& p : ordered) {
               auto a = (p->getTension() + p->getDragForce(eq_drag) + p->getGravitationalForce()) / p->mass;
               std::get<0>(p->acceleration) = std::get<0>(a);
               std::get<1>(p->acceleration) = std::get<1>(a);
               std::get<2>(p->acceleration) = std::get<2>(a);
               boundary_condition(p);
            }
            for (auto& p : ordered) {
               p->RK_X_sub.push(p->RK_velocity_sub.get_x());
               p->RK_velocity_sub.push(p->accelTranslational());
            }
            if (ordered[0]->RK_X_sub.finished) break;
         }

         // Apply RK4 result to actual positions/velocities
         double max_vel = 0;
         for (auto& p : ordered) {
            p->setX(p->RK_X_sub.get_x());
            auto v = p->RK_velocity_sub.get_x();
            std::get<0>(p->velocity) = std::get<0>(v);
            std::get<1>(p->velocity) = std::get<1>(v);
            std::get<2>(p->velocity) = std::get<2>(v);
            double vnorm = Norm(v);
            if (vnorm > max_vel) max_vel = vnorm;
         }

         // Snapshot for GUI
         if (step % snapshot_interval == 0 || (step > 1000 && max_vel < tol)) {
            std::cout << "SNAPSHOT {\"iter\":" << step
                      << ",\"norm_v\":" << max_vel
                      << ",\"positions\":[";
            for (int k = 0; k < n_nodes; ++k) {
               auto X = ordered[k]->X;
               std::cout << "[" << std::get<0>(X) << "," << std::get<1>(X) << "," << std::get<2>(X) << "]";
               if (k < n_nodes - 1) std::cout << ",";
            }
            std::cout << "]}" << std::endl;
         }

         if (step > 1000 && max_vel < tol) {
            std::cout << "Converged at step " << step << " (max_vel=" << max_vel << ")" << std::endl;
            break;
         }
         if (step == max_steps - 1) {
            std::cout << "WARNING: did not converge after " << max_steps
                      << " steps (max_vel=" << max_vel << ")" << std::endl;
            converged = false;
         }
      }
      mooring->DragForceCoefficient = saved_Cd;
      std::cout << "Equilibrium done." << std::endl;
   } else {
      std::cerr << "Unknown mode: " << mode << std::endl;
      return 1;
   }

   auto t_end = std::chrono::high_resolution_clock::now();
   double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

   std::cout << "Computation time: " << elapsed_ms << " ms" << std::endl;

   /* ---------------------- write result ---------------------- */

   writeResultJSON(output_dir / "result.json", mooring, elapsed_ms, converged);

   delete mooring;
   return 0;
}
