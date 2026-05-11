#pragma once

#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "Network.hpp"
#include "my_vtk.hpp"

namespace BEMMeshPipeline::RemeshDebug {

struct CandidatePatchRecord {
  T3Tddd tri;
  int attempt_id;
  int op;
  int reason;
  int success;
  int reject_code;
};

struct RemeshedPatchRecord {
  T3Tddd tri;
  int attempt_id;
  int op;
  int reason;
};

struct TriggerEdgeRecord {
  Tddd p0, p1;
  int attempt_id;
  int op;
  int reason;
  int success;
  int reject_code;
  int bc_type; // 0=Neumann, 1=Dirichlet, 2=Interface
};

struct SubsurfaceRejectFaceRecord {
  T3Tddd tri;
  int trigger_bc;      // 0=Neumann, 1=Dirichlet, 2=Interface
  int worst_side;      // 0=Neumann side, 1=Dirichlet side, 2=invalid
  int worst_line_flag; // 0=Neumann line, 1=Dirichlet line, 2=Interface line, 3=invalid
  double altitude_rel;
  double angle_deg;
};

inline int lineFlagIndex(const networkLine* l) {
  if (!l) return 3;
  if (l->BCInterface) return 2;
  if (l->Dirichlet) return 1;
  return 0;
}

inline T3Tddd snapshotFaceTri(const networkFace* f) {
  if (!f)
    return T3Tddd{{Tddd{0., 0., 0.}, Tddd{0., 0., 0.}, Tddd{0., 0., 0.}}};
  auto [p0, p1, p2] = f->getPoints();
  return {p0 ? p0->X : Tddd{0., 0., 0.},
          p1 ? p1->X : Tddd{0., 0., 0.},
          p2 ? p2->X : Tddd{0., 0., 0.}};
}

class Recorder {
 public:
  Recorder(bool enabled_in,
           const Network& water_in,
           int time_step_in,
           std::string output_directory_in,
           double simulation_time_in,
           PVDWriter* candidate_patches_pvd_in,
           PVDWriter* remeshed_patches_pvd_in,
           PVDWriter* trigger_edges_pvd_in)
      : enabled_(enabled_in),
        water_(water_in),
        time_step_(time_step_in),
        output_directory_(std::move(output_directory_in)),
        simulation_time_(simulation_time_in),
        candidate_patches_pvd_(candidate_patches_pvd_in),
        remeshed_patches_pvd_(remeshed_patches_pvd_in),
        trigger_edges_pvd_(trigger_edges_pvd_in) {}

  bool enabled() const { return enabled_; }

  std::vector<TriggerEdgeRecord>& triggerEdgeRecords() { return trigger_edge_records_; }
  const std::vector<TriggerEdgeRecord>& triggerEdgeRecords() const { return trigger_edge_records_; }

  void collectCandidatePatch(const Network& patch, int op, int reason, int aid,
                             int success, int reject_code) {
    if (!enabled_) return;
    for (auto* f : patch.getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      candidate_patch_records_.push_back({{p0->X, p1->X, p2->X}, aid, op, reason, success, reject_code});
    }
  }

  void collectRemeshedPatch(const Network& patch, int op, int reason, int aid) {
    if (!enabled_) return;
    for (auto* f : patch.getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      remeshed_patch_records_.push_back({{p0->X, p1->X, p2->X}, aid, op, reason});
    }
  }

  std::vector<T3Tddd> snapshotPatchTris(const Network& patch) const {
    std::vector<T3Tddd> tris;
    if (!enabled_) return tris;
    for (auto* f : patch.getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      tris.push_back({p0->X, p1->X, p2->X});
    }
    return tris;
  }

  void commitRemeshedSnapshot(const std::vector<T3Tddd>& tris, int op, int reason, int aid) {
    if (!enabled_) return;
    for (const auto& t : tris)
      remeshed_patch_records_.push_back({t, aid, op, reason});
  }

  void collectTriggerEdge(const networkLine* l, int op, int reason, int aid,
                          int success, int reject_code) {
    if (!enabled_ || !l) return;
    auto [p0, p1] = l->getPoints();
    if (!p0 || !p1) return;
    const int bc = l->BCInterface ? 2 : (l->Dirichlet ? 1 : 0);
    trigger_edge_records_.push_back({p0->X, p1->X, aid, op, reason, success, reject_code, bc});
  }

  void recordSubsurfaceRejectFace(const T3Tddd& tri, int trigger_bc,
                                  int worst_side, int worst_line_flag,
                                  double altitude_rel, double angle_deg) {
    if (!enabled_) return;
    subsurface_reject_face_records_.push_back({
        tri,
        trigger_bc,
        worst_side,
        worst_line_flag,
        altitude_rel,
        angle_deg});
  }

  void markAttemptSuccess(int aid, int reject_success_code) {
    if (!enabled_) return;
    for (auto& r : trigger_edge_records_)
      if (r.attempt_id == aid) { r.success = 1; r.reject_code = reject_success_code; }
    for (auto& r : candidate_patch_records_)
      if (r.attempt_id == aid) { r.success = 1; r.reject_code = reject_success_code; }
  }

  void markAttemptReplaceFailed(int aid, int replace_failed_code) {
    if (!enabled_) return;
    for (auto& r : trigger_edge_records_)
      if (r.attempt_id == aid) r.reject_code = replace_failed_code;
    for (auto& r : candidate_patch_records_)
      if (r.attempt_id == aid) r.reject_code = replace_failed_code;
  }

  void flushAll() {
    if (!enabled_) return;
    flushCandidatePatchVTU();
    flushRemeshedPatchVTU();
    flushTriggerEdgesVTU();
    flushSubsurfaceRejectFacesVTU();
  }

 private:
  template <class Records>
  void writeTrianglePatchVTU(const std::string& tag,
                             const Records& records,
                             bool has_success,
                             PVDWriter* pvd) {
    if (!enabled_ || records.empty())
      return;
    const std::string vtu_name = water_.getName() + "_" + tag + "_" + std::to_string(time_step_) + ".vtu";
    const std::string filename = output_directory_.empty()
                                     ? "/tmp/" + vtu_name
                                     : output_directory_ + "/" + vtu_name;
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp)
      return;
    const int n_cells = static_cast<int>(records.size());
    const int n_points = n_cells * 3;
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", n_cells, n_points);
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' format='ascii'>\n");
    for (const auto& r : records) {
      const auto& [a, b, c] = r.tri;
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(a), (float)std::get<1>(a), (float)std::get<2>(a));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(b), (float)std::get<1>(b), (float)std::get<2>(b));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(c), (float)std::get<1>(c), (float)std::get<2>(c));
    }
    fprintf(fp, "</DataArray>\n</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d %d %d\n", i * 3, i * 3 + 1, i * 3 + 2);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d\n", (i + 1) * 3);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "5\n");
    fprintf(fp, "</DataArray>\n</Cells>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type='Int32' Name='attempt_id' format='ascii'>\n");
    for (const auto& r : records) fprintf(fp, "%d\n", r.attempt_id);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='op' format='ascii'>\n");
    for (const auto& r : records) fprintf(fp, "%d\n", r.op);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='reason' format='ascii'>\n");
    for (const auto& r : records) fprintf(fp, "%d\n", r.reason);
    fprintf(fp, "</DataArray>\n");
    if (has_success) {
      fprintf(fp, "<DataArray type='Int32' Name='success' format='ascii'>\n");
      for (const auto& r : records) {
        if constexpr (requires { r.success; })
          fprintf(fp, "%d\n", r.success);
      }
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "<DataArray type='Int32' Name='reject_code' format='ascii'>\n");
      for (const auto& r : records) {
        if constexpr (requires { r.reject_code; })
          fprintf(fp, "%d\n", r.reject_code);
      }
      fprintf(fp, "</DataArray>\n");
    }
    fprintf(fp, "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(fp);
    std::cout << "[patch debug] wrote " << filename << " (" << n_cells << " faces)" << std::endl;
    if (pvd) {
      pvd->push(vtu_name, simulation_time_);
      pvd->output();
    }
  }

  void flushCandidatePatchVTU() {
    writeTrianglePatchVTU("candidate_patches", candidate_patch_records_, true, candidate_patches_pvd_);
  }

  void flushRemeshedPatchVTU() {
    writeTrianglePatchVTU("remeshed_patches", remeshed_patch_records_, false, remeshed_patches_pvd_);
  }

  void flushTriggerEdgesVTU() {
    if (!enabled_ || trigger_edge_records_.empty())
      return;
    const std::string vtu_name = water_.getName() + "_trigger_edges_" + std::to_string(time_step_) + ".vtu";
    const std::string filename = output_directory_.empty()
                                     ? "/tmp/" + vtu_name
                                     : output_directory_ + "/" + vtu_name;
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp)
      return;
    const int n_cells = static_cast<int>(trigger_edge_records_.size());
    const int n_points = n_cells * 2;
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", n_cells, n_points);
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' format='ascii'>\n");
    for (const auto& r : trigger_edge_records_) {
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(r.p0), (float)std::get<1>(r.p0), (float)std::get<2>(r.p0));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(r.p1), (float)std::get<1>(r.p1), (float)std::get<2>(r.p1));
    }
    fprintf(fp, "</DataArray>\n</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d %d\n", i * 2, i * 2 + 1);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d\n", (i + 1) * 2);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "3\n");
    fprintf(fp, "</DataArray>\n</Cells>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type='Int32' Name='attempt_id' format='ascii'>\n");
    for (const auto& r : trigger_edge_records_) fprintf(fp, "%d\n", r.attempt_id);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='op' format='ascii'>\n");
    for (const auto& r : trigger_edge_records_) fprintf(fp, "%d\n", r.op);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='reason' format='ascii'>\n");
    for (const auto& r : trigger_edge_records_) fprintf(fp, "%d\n", r.reason);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='success' format='ascii'>\n");
    for (const auto& r : trigger_edge_records_) fprintf(fp, "%d\n", r.success);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='reject_code' format='ascii'>\n");
    for (const auto& r : trigger_edge_records_) fprintf(fp, "%d\n", r.reject_code);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='bc_type' format='ascii'>\n");
    for (const auto& r : trigger_edge_records_) fprintf(fp, "%d\n", r.bc_type);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(fp);
    std::cout << "[patch debug] wrote " << filename << " (" << n_cells << " edges)" << std::endl;
    if (trigger_edges_pvd_) {
      trigger_edges_pvd_->push(vtu_name, simulation_time_);
      trigger_edges_pvd_->output();
    }
  }

  void flushSubsurfaceRejectFacesVTU() {
    if (!enabled_ || subsurface_reject_face_records_.empty())
      return;
    const std::string vtu_name = water_.getName() + "_subsurface_reject_faces_" + std::to_string(time_step_) + ".vtu";
    const std::string filename = output_directory_.empty()
                                     ? "/tmp/" + vtu_name
                                     : output_directory_ + "/" + vtu_name;
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp)
      return;
    const int n_cells = static_cast<int>(subsurface_reject_face_records_.size());
    const int n_points = n_cells * 3;
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", n_cells, n_points);
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' format='ascii'>\n");
    for (const auto& r : subsurface_reject_face_records_) {
      const auto& [a, b, c] = r.tri;
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(a), (float)std::get<1>(a), (float)std::get<2>(a));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(b), (float)std::get<1>(b), (float)std::get<2>(b));
      fprintf(fp, "%f %f %f\n", (float)std::get<0>(c), (float)std::get<1>(c), (float)std::get<2>(c));
    }
    fprintf(fp, "</DataArray>\n</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d %d %d\n", i * 3, i * 3 + 1, i * 3 + 2);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "%d\n", (i + 1) * 3);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (int i = 0; i < n_cells; ++i) fprintf(fp, "5\n");
    fprintf(fp, "</DataArray>\n</Cells>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type='Int32' Name='trigger_bc' format='ascii'>\n");
    for (const auto& r : subsurface_reject_face_records_) fprintf(fp, "%d\n", r.trigger_bc);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='worst_face_side' format='ascii'>\n");
    for (const auto& r : subsurface_reject_face_records_) fprintf(fp, "%d\n", r.worst_side);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='worst_line_flag' format='ascii'>\n");
    for (const auto& r : subsurface_reject_face_records_) fprintf(fp, "%d\n", r.worst_line_flag);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Float32' Name='altitude_rel' format='ascii'>\n");
    for (const auto& r : subsurface_reject_face_records_) fprintf(fp, "%g\n", (float)r.altitude_rel);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Float32' Name='edge_angle_deg' format='ascii'>\n");
    for (const auto& r : subsurface_reject_face_records_) fprintf(fp, "%g\n", (float)r.angle_deg);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(fp);
    std::cout << "[patch debug] wrote " << filename << " (" << n_cells << " subsurface reject faces)" << std::endl;
  }

  bool enabled_ = false;
  const Network& water_;
  int time_step_ = 0;
  std::string output_directory_;
  double simulation_time_ = 0.0;
  PVDWriter* candidate_patches_pvd_ = nullptr;
  PVDWriter* remeshed_patches_pvd_ = nullptr;
  PVDWriter* trigger_edges_pvd_ = nullptr;
  std::vector<CandidatePatchRecord> candidate_patch_records_;
  std::vector<RemeshedPatchRecord> remeshed_patch_records_;
  std::vector<TriggerEdgeRecord> trigger_edge_records_;
  std::vector<SubsurfaceRejectFaceRecord> subsurface_reject_face_records_;
};

} // namespace BEMMeshPipeline::RemeshDebug
