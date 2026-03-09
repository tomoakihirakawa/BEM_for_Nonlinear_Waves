#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "Network.hpp"
#include "VPM.hpp"

/*
 * BEM Checkpoint file format (version 3)
 *
 * [Header]
 *   magic:           char[8]   "BEMCKPT\0"
 *   version:         uint32
 *   time_step:       int32
 *   simulation_time: double
 *   dt:              double
 *
 * [Fluid Objects]
 *   num_fluids: uint32
 *   per fluid: {
 *     name_length: uint32,  name: char[name_length]
 *     num_points: uint32
 *     per point: { X: double[3], phiphin: double[2], phiphin_t: double[2], phi_Dirichlet: double }
 *     num_faces: uint32
 *     per face:  { indices: uint32[3] }
 *     has_quadratic: uint8
 *     if has_quadratic:
 *       num_lines: uint32
 *       per line: { endpoint_indices: uint32[2], X_mid: double[3], phiphin: double[2], phiphin_t: double[2], corner_offset: double[3] }
 *   }
 *
 * [Rigid Bodies]
 *   num_rigid: uint32
 *   per body: { name_length: uint32, name: char[], COM: double[3], Q: double[4], velocity: double[6], acceleration: double[6],
 *               num_points: uint32, positions: double[3*num_points] }
 *
 * [Soft Bodies]
 *   num_soft: uint32
 *   per body: { name_length: uint32, name: char[], COM: double[3], Q: double[4], velocity: double[6], num_points: uint32, positions: double[3*num_points] }
 *
 * [VPM Particles]
 *   num_particles: uint32
 *   per particle: { x: double[3], alpha: double[3], sigma: double, volume: double }
 */

namespace BEM_Checkpoint {

static constexpr char MAGIC[8] = {'B', 'E', 'M', 'C', 'K', 'P', 'T', '\0'};
static constexpr uint32_t VERSION = 3;

// ============================================================================
// Low-level I/O helpers
// ============================================================================

inline void write_raw(FILE* fp, const void* data, size_t bytes) {
   if (std::fwrite(data, 1, bytes, fp) != bytes)
      throw std::runtime_error("BEM_Checkpoint: write error");
}

inline void read_raw(FILE* fp, void* data, size_t bytes) {
   if (std::fread(data, 1, bytes, fp) != bytes)
      throw std::runtime_error("BEM_Checkpoint: read error (unexpected EOF)");
}

template <typename T>
void write_val(FILE* fp, const T& v) { write_raw(fp, &v, sizeof(T)); }

template <typename T>
T read_val(FILE* fp) {
   T v;
   read_raw(fp, &v, sizeof(T));
   return v;
}

inline void write_string(FILE* fp, const std::string& s) {
   uint32_t len = static_cast<uint32_t>(s.size());
   write_val(fp, len);
   write_raw(fp, s.data(), len);
}

inline std::string read_string(FILE* fp) {
   uint32_t len = read_val<uint32_t>(fp);
   std::string s(len, '\0');
   read_raw(fp, s.data(), len);
   return s;
}

// ============================================================================
// writeCheckpoint
// ============================================================================

inline void writeCheckpointToPath(const std::filesystem::path& path,
                                  int time_step,
                                  double simulation_time,
                                  double dt,
                                  const std::vector<Network*>& FluidObject,
                                  const std::vector<Network*>& RigidBodyObject,
                                  const std::vector<Network*>& SoftBodyObject,
                                  const VortexMethod& vpm,
                                  bool use_true_quadratic_element) {
   FILE* fp = std::fopen(path.c_str(), "wb");
   if (!fp)
      throw std::runtime_error("BEM_Checkpoint: cannot open " + path.string());

   // --- Header ---
   write_raw(fp, MAGIC, 8);
   write_val(fp, VERSION);
   write_val(fp, static_cast<int32_t>(time_step));
   write_val(fp, simulation_time);
   write_val(fp, dt);

   // --- Fluid Objects ---
   write_val(fp, static_cast<uint32_t>(FluidObject.size()));
   for (const auto* water : FluidObject) {
      write_string(fp, water->getName());

      // Build stable point ordering: unordered_set → vector + index map
      const auto& pts_set = water->getPoints();
      std::vector<networkPoint*> pts(pts_set.begin(), pts_set.end());
      std::map<networkPoint*, uint32_t> pt_idx;
      for (uint32_t i = 0; i < pts.size(); ++i)
         pt_idx[pts[i]] = i;

      // Points
      write_val(fp, static_cast<uint32_t>(pts.size()));
      for (const auto* p : pts) {
         Tddd X = p->getXtuple();
         write_raw(fp, &X, sizeof(Tddd));
         write_raw(fp, &p->phiphin, sizeof(Tdd));
         write_raw(fp, &p->phiphin_t, sizeof(Tdd));
         write_val(fp, p->phi_Dirichlet);
      }

      // Faces
      auto faces = water->getBoundaryFaces();
      write_val(fp, static_cast<uint32_t>(faces.size()));
      for (const auto* f : faces) {
         auto [p0, p1, p2] = f->getPoints();
         uint32_t idx[3] = {pt_idx.at(p0), pt_idx.at(p1), pt_idx.at(p2)};
         write_raw(fp, idx, sizeof(idx));
      }

      // Quadratic line data
      uint8_t has_quad = use_true_quadratic_element ? 1 : 0;
      write_val(fp, has_quad);
      if (use_true_quadratic_element) {
         auto lines = water->getBoundaryLines();
         write_val(fp, static_cast<uint32_t>(lines.size()));
         for (const auto* l : lines) {
            auto [pA, pB] = l->getPoints();
            uint32_t ep[2] = {pt_idx.at(pA), pt_idx.at(pB)};
            write_raw(fp, ep, sizeof(ep));
            write_raw(fp, &l->X_mid, sizeof(Tddd));
            write_raw(fp, &l->phiphin, sizeof(Tdd));
            write_raw(fp, &l->phiphin_t, sizeof(Tdd));
            write_raw(fp, &l->corner_midpoint_offset, sizeof(Tddd));
         }
      }
   }

   // --- Rigid Bodies ---
   write_val(fp, static_cast<uint32_t>(RigidBodyObject.size()));
   for (const auto* net : RigidBodyObject) {
      write_string(fp, net->getName());
      write_raw(fp, &net->COM, sizeof(Tddd));
      T4d q = net->Q();
      write_raw(fp, &q, sizeof(T4d));
      write_raw(fp, &net->velocity, sizeof(T6d));
      write_raw(fp, &net->acceleration, sizeof(T6d));
      const auto& pts_set = net->getPoints();
      std::vector<networkPoint*> pts(pts_set.begin(), pts_set.end());
      std::ranges::sort(pts, [](const networkPoint* a, const networkPoint* b) {
         return std::tie(a->initialX[0], a->initialX[1], a->initialX[2]) <
                std::tie(b->initialX[0], b->initialX[1], b->initialX[2]);
      });
      write_val(fp, static_cast<uint32_t>(pts.size()));
      for (const auto* p : pts) {
         Tddd X = p->getXtuple();
         write_raw(fp, &X, sizeof(Tddd));
      }
   }

   // --- Soft Bodies ---
   write_val(fp, static_cast<uint32_t>(SoftBodyObject.size()));
   for (const auto* net : SoftBodyObject) {
      write_string(fp, net->getName());
      write_raw(fp, &net->COM, sizeof(Tddd));
      T4d q = net->Q();
      write_raw(fp, &q, sizeof(T4d));
      write_raw(fp, &net->velocity, sizeof(T6d));
      const auto& pts_set = net->getPoints();
      std::vector<networkPoint*> pts(pts_set.begin(), pts_set.end());
      write_val(fp, static_cast<uint32_t>(pts.size()));
      for (const auto* p : pts) {
         Tddd X = p->getXtuple();
         write_raw(fp, &X, sizeof(Tddd));
      }
   }

   // --- VPM Particles ---
   const auto& particles = vpm.getParticles();
   write_val(fp, static_cast<uint32_t>(particles.size()));
   for (const auto& p : particles) {
      write_raw(fp, &p.x, sizeof(std::array<double, 3>));
      write_raw(fp, &p.alpha, sizeof(std::array<double, 3>));
      write_val(fp, p.sigma);
      write_val(fp, p.volume);
   }

   std::fclose(fp);
   std::cout << Green << "[checkpoint] Saved: " << path.filename().string()
             << " (t=" << simulation_time << ", step=" << time_step << ")"
             << colorReset << std::endl;
}

inline void writeCheckpoint(const std::filesystem::path& dir,
                            int time_step,
                            double simulation_time,
                            double dt,
                            const std::vector<Network*>& FluidObject,
                            const std::vector<Network*>& RigidBodyObject,
                            const std::vector<Network*>& SoftBodyObject,
                            const VortexMethod& vpm,
                            bool use_true_quadratic_element) {
   char fname[64];
   std::snprintf(fname, sizeof(fname), "checkpoint_%06d.bin", time_step);
   writeCheckpointToPath(dir / fname, time_step, simulation_time, dt,
                         FluidObject, RigidBodyObject, SoftBodyObject, vpm,
                         use_true_quadratic_element);
}

// ============================================================================
// readCheckpointHeader  (lightweight: reads only the header)
// ============================================================================

struct CheckpointHeader {
   uint32_t version;
   int32_t time_step;
   double simulation_time;
   double dt;
};

inline bool readCheckpointHeader(const std::filesystem::path& filepath, CheckpointHeader& hdr) {
   FILE* fp = std::fopen(filepath.c_str(), "rb");
   if (!fp) return false;

   char magic[8];
   if (std::fread(magic, 1, 8, fp) != 8 || std::memcmp(magic, MAGIC, 8) != 0) {
      std::fclose(fp);
      return false;
   }
   hdr.version = read_val<uint32_t>(fp);
   hdr.time_step = read_val<int32_t>(fp);
   hdr.simulation_time = read_val<double>(fp);
   hdr.dt = read_val<double>(fp);
   std::fclose(fp);
   return true;
}

// ============================================================================
// readCheckpoint
// ============================================================================

inline std::tuple<int, double, double>
readCheckpoint(const std::filesystem::path& filepath,
               std::vector<Network*>& FluidObject,
               const std::vector<Network*>& RigidBodyObject,
               const std::vector<Network*>& SoftBodyObject,
               VortexMethod& vpm,
               bool use_true_quadratic_element) {
   FILE* fp = std::fopen(filepath.c_str(), "rb");
   if (!fp)
      throw std::runtime_error("BEM_Checkpoint: cannot open " + filepath.string());

   // --- Header ---
   char magic[8];
   read_raw(fp, magic, 8);
   if (std::memcmp(magic, MAGIC, 8) != 0)
      throw std::runtime_error("BEM_Checkpoint: invalid magic in " + filepath.string());

   uint32_t version = read_val<uint32_t>(fp);
   if (version < 1 || version > VERSION)
      throw std::runtime_error("BEM_Checkpoint: unsupported version " + std::to_string(version));

   int32_t ckpt_step = read_val<int32_t>(fp);
   double ckpt_time = read_val<double>(fp);
   double ckpt_dt = read_val<double>(fp);

   // --- Fluid Objects ---
   uint32_t num_fluids = read_val<uint32_t>(fp);
   if (num_fluids != FluidObject.size())
      throw std::runtime_error("BEM_Checkpoint: fluid object count mismatch (file=" +
                               std::to_string(num_fluids) + ", sim=" +
                               std::to_string(FluidObject.size()) + ")");

   for (size_t fi = 0; fi < num_fluids; ++fi) {
      std::string ckpt_name = read_string(fp);

      // Match by index (assume same order as settings)
      Network* old_net = FluidObject[fi];
      if (old_net->getName() != ckpt_name)
         std::cerr << Yellow << "[checkpoint] WARNING: name mismatch: expected '"
                   << old_net->getName() << "', got '" << ckpt_name << "'" << colorReset << std::endl;

      // Read point data
      uint32_t num_points = read_val<uint32_t>(fp);
      struct PointRecord {
         Tddd X;
         Tdd phiphin;
         Tdd phiphin_t;
         double phi_Dirichlet;
      };
      std::vector<PointRecord> point_records(num_points);
      for (uint32_t i = 0; i < num_points; ++i) {
         auto& r = point_records[i];
         read_raw(fp, &r.X, sizeof(Tddd));
         read_raw(fp, &r.phiphin, sizeof(Tdd));
         read_raw(fp, &r.phiphin_t, sizeof(Tdd));
         r.phi_Dirichlet = read_val<double>(fp);
      }

      // Read face indices
      uint32_t num_faces = read_val<uint32_t>(fp);
      std::vector<std::vector<int>> face_indices(num_faces);
      for (uint32_t i = 0; i < num_faces; ++i) {
         uint32_t idx[3];
         read_raw(fp, idx, sizeof(idx));
         face_indices[i] = {static_cast<int>(idx[0]), static_cast<int>(idx[1]), static_cast<int>(idx[2])};
      }

      // Read quadratic data
      uint8_t has_quad = read_val<uint8_t>(fp);
      struct LineRecord {
         uint32_t ep[2];
         Tddd X_mid;
         Tdd phiphin;
         Tdd phiphin_t;
         Tddd corner_offset;
      };
      std::vector<LineRecord> line_records;
      if (has_quad) {
         uint32_t num_lines = read_val<uint32_t>(fp);
         line_records.resize(num_lines);
         for (uint32_t i = 0; i < num_lines; ++i) {
            auto& r = line_records[i];
            read_raw(fp, r.ep, sizeof(r.ep));
            read_raw(fp, &r.X_mid, sizeof(Tddd));
            if (version >= 3) {
               read_raw(fp, &r.phiphin, sizeof(Tdd));
               read_raw(fp, &r.phiphin_t, sizeof(Tdd));
            } else {
               r.phiphin = {read_val<double>(fp), 0.0};
               r.phiphin_t = {0.0, 0.0};
            }
            read_raw(fp, &r.corner_offset, sizeof(Tddd));
         }
      }

      // --- Reconstruct Network ---
      // Preserve settings from old network
      JSON saved_inputJSON = old_net->inputJSON;
      std::string saved_name = old_net->getName();
      auto saved_velocity_name_start = old_net->velocity_name_start;
      auto saved_isFixed = old_net->isFixed;
      // ic_phi/ic_eta are std::function lambdas that capture shared_ptr,
      // so they survive the old network deletion.
      auto saved_ic_phi = old_net->ic_phi;
      auto saved_ic_eta = old_net->ic_eta;
      auto saved_absorb_velocity = old_net->absorb_velocity;
      auto saved_absorb_gradPhi_t = old_net->absorb_gradPhi_t;
      auto saved_absorb_phi = old_net->absorb_phi;
      auto saved_absorb_eta = old_net->absorb_eta;
      auto saved_absorb_gamma = old_net->absorb_gamma;
      auto saved_water_wave_theory = old_net->water_wave_theory;

      // Create new Network
      Network* new_net = new Network(); // empty constructor

      // Set vertices
      std::vector<Tddd> vertices(num_points);
      for (uint32_t i = 0; i < num_points; ++i)
         vertices[i] = point_records[i].X;
      V_netPp points = new_net->setPoints(vertices);

      // Set faces
      new_net->setFaces(face_indices, points);

      // Set physical quantities on each point
      // points[i] corresponds to vertex i — but setPoints may merge duplicates.
      // Since checkpoint saves unique points, there should be no merging.
      for (uint32_t i = 0; i < num_points; ++i) {
         networkPoint* p = points[i];
         if (!p) continue;
         p->phiphin = point_records[i].phiphin;
         p->phiphin_t = point_records[i].phiphin_t;
         p->phi_Dirichlet = point_records[i].phi_Dirichlet;
      }

      // Restore quadratic line data
      if (has_quad && use_true_quadratic_element) {
         // Build point-pair → line map
         std::map<std::pair<networkPoint*, networkPoint*>, networkLine*> edge_map;
         for (auto* l : new_net->getBoundaryLines()) {
            auto [pA, pB] = l->getPoints();
            // Store both orderings
            edge_map[{pA, pB}] = l;
            edge_map[{pB, pA}] = l;
         }
         for (const auto& lr : line_records) {
            networkPoint* pA = points[lr.ep[0]];
            networkPoint* pB = points[lr.ep[1]];
            auto it = edge_map.find({pA, pB});
            if (it != edge_map.end()) {
               networkLine* l = it->second;
               l->X_mid = lr.X_mid;
               l->phiphin = lr.phiphin;
               l->phiphin_t = lr.phiphin_t;
               l->corner_midpoint_offset = lr.corner_offset;
            }
         }
      } else if (!has_quad && use_true_quadratic_element) {
         std::cerr << Yellow
                   << "[checkpoint] WARNING: loading linear checkpoint into true quadratic mode for '"
                   << saved_name
                   << "'. Reconstructing midpoint state from endpoint values."
                   << colorReset << std::endl;
         for (auto* l : new_net->getBoundaryLines()) {
            auto [pA, pB] = l->getPoints();
            if (!pA || !pB)
               continue;
            const Tddd X_mid = 0.5 * (pA->X + pB->X);
            l->setXSingle(X_mid);
            l->corner_midpoint_offset = {0., 0., 0.};
            l->phiphin = 0.5 * (pA->phiphin + pB->phiphin);
            l->phiphin_t = 0.5 * (pA->phiphin_t + pB->phiphin_t);
         }
      }

      // Restore metadata
      new_net->setName(saved_name);
      new_net->inputJSON = saved_inputJSON;
      new_net->velocity_name_start = saved_velocity_name_start;
      new_net->isFixed = saved_isFixed;
      new_net->isFluid = true;
      new_net->isRigidBody = false;
      new_net->isSoftBody = false;
      new_net->isAbsorber = false;
      new_net->ic_phi = saved_ic_phi;
      new_net->ic_eta = saved_ic_eta;
      new_net->absorb_velocity = saved_absorb_velocity;
      new_net->absorb_gradPhi_t = saved_absorb_gradPhi_t;
      new_net->absorb_phi = saved_absorb_phi;
      new_net->absorb_eta = saved_absorb_eta;
      new_net->absorb_gamma = saved_absorb_gamma;
      new_net->water_wave_theory = saved_water_wave_theory;

      // Replace in FluidObject array
      delete old_net;
      FluidObject[fi] = new_net;
   }

   // --- Rigid Bodies ---
   uint32_t num_rigid = read_val<uint32_t>(fp);
   if (num_rigid != RigidBodyObject.size())
      throw std::runtime_error("BEM_Checkpoint: rigid body count mismatch");

   for (size_t ri = 0; ri < num_rigid; ++ri) {
      std::string ckpt_name = read_string(fp);
      Network* net = RigidBodyObject[ri];

      Tddd com;
      T4d q;
      T6d vel, acc;
      read_raw(fp, &com, sizeof(Tddd));
      read_raw(fp, &q, sizeof(T4d));
      read_raw(fp, &vel, sizeof(T6d));
      read_raw(fp, &acc, sizeof(T6d));
      std::vector<Tddd> point_positions;
      if (version >= 2) {
         uint32_t num_pts = read_val<uint32_t>(fp);
         point_positions.resize(num_pts);
         for (uint32_t i = 0; i < num_pts; ++i)
            read_raw(fp, &point_positions[i], sizeof(Tddd));
      }

      // Restore rigid-body state without mutating initialX.
      // initialX is the reference geometry used by rigidTransformation(),
      // so overwriting it during rollback causes later pose commits to jump.
      net->COM = com;
      net->Q = Quaternion(q);
      net->velocity = vel;
      net->acceleration = acc;
      if (version >= 2) {
         const auto& pts_set = net->getPoints();
         std::vector<networkPoint*> pts(pts_set.begin(), pts_set.end());
         std::ranges::sort(pts, [](const networkPoint* a, const networkPoint* b) {
            return std::tie(a->initialX[0], a->initialX[1], a->initialX[2]) <
                   std::tie(b->initialX[0], b->initialX[1], b->initialX[2]);
         });
         if (point_positions.size() != pts.size())
            throw std::runtime_error("BEM_Checkpoint: rigid body point count mismatch for " + ckpt_name);
         for (size_t i = 0; i < pts.size(); ++i)
            pts[i]->setXSingle(point_positions[i]);
      } else {
         for (const auto* p : net->getPoints())
            const_cast<networkPoint*>(p)->setXSingle(net->rigidTransformation(p->initialX));
      }
      net->setGeometricPropertiesForce();
   }

   // --- Soft Bodies ---
   uint32_t num_soft = read_val<uint32_t>(fp);
   if (num_soft != SoftBodyObject.size())
      throw std::runtime_error("BEM_Checkpoint: soft body count mismatch");

   for (size_t si = 0; si < num_soft; ++si) {
      std::string ckpt_name = read_string(fp);
      Network* net = SoftBodyObject[si];

      Tddd com;
      T4d q;
      T6d vel;
      read_raw(fp, &com, sizeof(Tddd));
      read_raw(fp, &q, sizeof(T4d));
      read_raw(fp, &vel, sizeof(T6d));

      net->COM = com;
      net->Q = Quaternion(q);
      net->velocity = vel;

      uint32_t num_pts = read_val<uint32_t>(fp);
      // Restore point positions in iteration order
      // (same order as was saved: unordered_set → vector)
      const auto& pts_set = net->getPoints();
      std::vector<networkPoint*> pts(pts_set.begin(), pts_set.end());
      if (num_pts != pts.size())
         throw std::runtime_error("BEM_Checkpoint: soft body point count mismatch for " + ckpt_name);
      for (uint32_t i = 0; i < num_pts; ++i) {
         Tddd X;
         read_raw(fp, &X, sizeof(Tddd));
         pts[i]->setXSingle(X);
      }
      net->setGeometricPropertiesForce();
   }

   // --- VPM Particles ---
   uint32_t num_particles = read_val<uint32_t>(fp);
   auto& plist = vpm.getParticles();
   plist.clear();
   for (uint32_t i = 0; i < num_particles; ++i) {
      std::array<double, 3> x, alpha;
      read_raw(fp, &x, sizeof(x));
      read_raw(fp, &alpha, sizeof(alpha));
      double sigma = read_val<double>(fp);
      double volume = read_val<double>(fp);
      vpm.addParticle(x, alpha, sigma, volume);
   }

   std::fclose(fp);
   std::cout << Green << "[checkpoint] Loaded: " << filepath.filename().string()
             << " (t=" << ckpt_time << ", step=" << ckpt_step << ")"
             << colorReset << std::endl;

   return {static_cast<int>(ckpt_step), ckpt_time, ckpt_dt};
}

// ============================================================================
// pruneCheckpoints — keep only the most recent `max_keep` checkpoint files
// ============================================================================

inline void pruneCheckpoints(const std::filesystem::path& dir, int max_keep) {
   if (max_keep <= 0) return; // 0 = keep all

   std::regex pattern(R"(checkpoint_(\d+)\.bin)");
   std::vector<std::pair<int, std::filesystem::path>> checkpoints;

   for (const auto& entry : std::filesystem::directory_iterator(dir)) {
      if (!entry.is_regular_file()) continue;
      std::string fname = entry.path().filename().string();
      std::smatch m;
      if (std::regex_match(fname, m, pattern)) {
         int step = std::stoi(m[1].str());
         checkpoints.emplace_back(step, entry.path());
      }
   }

   if (static_cast<int>(checkpoints.size()) <= max_keep) return;

   // Sort by step descending
   std::sort(checkpoints.begin(), checkpoints.end(),
             [](const auto& a, const auto& b) { return a.first > b.first; });

   // Remove older ones
   for (size_t i = max_keep; i < checkpoints.size(); ++i) {
      std::filesystem::remove(checkpoints[i].second);
      std::cout << "[checkpoint] Pruned: " << checkpoints[i].second.filename().string() << std::endl;
   }
}

} // namespace BEM_Checkpoint
