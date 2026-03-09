#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <set>

#include "BEM.hpp"
#include "OutputCommon.hpp"
#include "vtkWriter.hpp"
#include "VPM.hpp"

namespace OutputParaView {

// Element-type-aware VTU output: linear faces → VTK_TRIANGLE (type 5),
// pseudo-quadratic/true-quadratic faces → VTK_QUADRATIC_TRIANGLE (type 22, 6 nodes).
inline void mk_vtu_quadratic(const std::string &filename, const V_netFp &Faces, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  try {
#if defined(debug_mk_vtu)
    std::cout << Magenta << filename << colorReset;
    std::cout << "  Faces.size() : " << std::to_string(Faces.size()) << colorReset << " ";
#endif
    struct PointEntry {
      Tddd position;
      networkPoint *ptr;      // non-null for vertex nodes
      networkPoint *mid_pA;   // midpoint endpoint A
      networkPoint *mid_pB;   // midpoint endpoint B
      networkLine *mid_line;  // midpoint edge
      bool is_true_quad_mid;
    };

    std::vector<PointEntry> all_points;
    std::vector<int> cell_sizes;
    std::vector<uint8_t> cell_types;

    for (const auto &f : Faces) {
      auto [p0, p1, p2] = f->getPoints();
      all_points.push_back({p0->getXtuple(), p0, nullptr, nullptr, nullptr, false});
      all_points.push_back({p1->getXtuple(), p1, nullptr, nullptr, nullptr, false});
      all_points.push_back({p2->getXtuple(), p2, nullptr, nullptr, nullptr, false});

      if (f->isPseudoQuadraticElement || f->isTrueQuadraticElement) {
        auto [p0_, l0, p1_, l1, p2_, l2] = f->PLPLPL;
        Tddd mid0, mid1, mid2;

        if (f->isTrueQuadraticElement) {
          mid0 = l0->X_mid;
          mid1 = l1->X_mid;
          mid2 = l2->X_mid;
        } else {
          // Pseudo-quadratic: use DodecaPoints
          // dodecaPoints[0] origin=p0: p0 at (1,0), p1 at (0,1), p2 at (0,0)
          if (f->dodecaPoints[0]) {
            auto &dp = f->dodecaPoints[0];
            mid0 = dp->X(0.5, 0.5) + dp->corner_offset(0.5, 0.5); // l0(p0-p1)
            mid1 = dp->X(0.0, 0.5) + dp->corner_offset(0.0, 0.5); // l1(p1-p2)
            mid2 = dp->X(0.5, 0.0) + dp->corner_offset(0.5, 0.0); // l2(p2-p0)
          } else {
            mid0 = 0.5 * (p0->getXtuple() + p1->getXtuple());
            mid1 = 0.5 * (p1->getXtuple() + p2->getXtuple());
            mid2 = 0.5 * (p2->getXtuple() + p0->getXtuple());
          }
        }

        auto [pa0, pb0] = l0->getPoints();
        auto [pa1, pb1] = l1->getPoints();
        auto [pa2, pb2] = l2->getPoints();
        all_points.push_back({mid0, nullptr, pa0, pb0, l0, f->isTrueQuadraticElement});
        all_points.push_back({mid1, nullptr, pa1, pb1, l1, f->isTrueQuadraticElement});
        all_points.push_back({mid2, nullptr, pa2, pb2, l2, f->isTrueQuadraticElement});

        cell_sizes.push_back(6);
        cell_types.push_back(22); // VTK_QUADRATIC_TRIANGLE
      } else {
        cell_sizes.push_back(3);
        cell_types.push_back(5); // VTK_TRIANGLE
      }
    }

    const int num_points = (int)all_points.size();
    const int num_cells = (int)cell_sizes.size();

    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened");

    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", num_cells, num_points);

    // Points
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
    for (const auto &pe : all_points)
      fprintf(fp, "%s %s %s ",
              NumtoString(std::get<0>(pe.position)).c_str(),
              NumtoString(std::get<1>(pe.position)).c_str(),
              NumtoString(std::get<2>(pe.position)).c_str());
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "</Points>\n");

    // PointData
    if (!VV_name_comp_mapPVd.empty()) {
      fprintf(fp, "<PointData>\n");
      for (const auto &V_var : VV_name_comp_mapPVd) {
        std::string Name = std::get<std::string>(V_var[0]);

        if (V_var.size() > 1 && std::holds_alternative<uomap_P_d>(V_var[1])) {
          const auto &smap = std::get<uomap_P_d>(V_var[1]);
          fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='%s' format='ascii'>\n", Name.c_str());
          for (const auto &pe : all_points) {
            if (pe.ptr) {
              auto it = smap.find(pe.ptr);
              if (it != smap.end() && isFinite(it->second))
                fprintf(fp, "%s ", NumtoString(it->second).c_str());
              else
                fprintf(fp, "NaN ");
            } else {
              double val = 0.;
              bool found = false;
              if (pe.mid_line) {
                if (Name == "direction_info_count") { val = pe.mid_line->debug_direction_info_count; found = true; }
                else if (Name == "contact_faces_count") { val = pe.mid_line->debug_contact_faces_count; found = true; }
                else if (Name == "body_vertices_count") { val = pe.mid_line->debug_body_vertices_count; found = true; }
                else if (Name == "isInContact_pass_count") { val = pe.mid_line->debug_isInContact_pass_count; found = true; }
                else if (pe.is_true_quad_mid) {
                  if (Name == "φ") { val = pe.mid_line->phiphin[0]; found = true; }
                  else if (Name == "φn") {
                    if (pe.mid_line->phinOnFace.count(nullptr))
                      val = pe.mid_line->phinOnFace.at(nullptr);
                    else {
                      // CORNER midpoint: area-weighted average for visualization only
                      double wa = 0., wp = 0.;
                      for (auto* f : pe.mid_line->getBoundaryFaces())
                        if (f && pe.mid_line->phinOnFace.count(f)) { wp += pe.mid_line->phinOnFace.at(f) * f->area; wa += f->area; }
                      val = (wa > 0.) ? wp / wa : pe.mid_line->phiphin[1];
                    }
                    found = true;
                  }
                  else if (Name == "φt") { val = pe.mid_line->phiphin_t[0]; found = true; }
                  else if (Name == "φnt") { val = pe.mid_line->phiphin_t[1]; found = true; }
                  else if (Name == "diag") { val = pe.mid_line->diag_coeff_BEM; found = true; }
                }
              }
              if (!found) {
                auto itA = smap.find(pe.mid_pA);
                auto itB = smap.find(pe.mid_pB);
                if (itA != smap.end() && itB != smap.end() &&
                    isFinite(itA->second) && isFinite(itB->second)) {
                  val = 0.5 * (itA->second + itB->second);
                  found = true;
                }
              }
              if (found && isFinite(val))
                fprintf(fp, "%s ", NumtoString(val).c_str());
              else
                fprintf(fp, "NaN ");
            }
          }
          fprintf(fp, "\n</DataArray>\n");
        } else if (V_var.size() > 1 && std::holds_alternative<uomap_P_Tddd>(V_var[1])) {
          const auto &vmap = std::get<uomap_P_Tddd>(V_var[1]);
          fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='%s' format='ascii'>\n", Name.c_str());
          for (const auto &pe : all_points) {
            if (pe.ptr) {
              auto it = vmap.find(pe.ptr);
              if (it != vmap.end()) {
                const auto &v = it->second;
                fprintf(fp, "%s %s %s ",
                        isFinite(std::get<0>(v)) ? NumtoString(std::get<0>(v)).c_str() : "NaN",
                        isFinite(std::get<1>(v)) ? NumtoString(std::get<1>(v)).c_str() : "NaN",
                        isFinite(std::get<2>(v)) ? NumtoString(std::get<2>(v)).c_str() : "NaN");
              } else {
                fprintf(fp, "NaN NaN NaN ");
              }
            } else {
              auto itA = vmap.find(pe.mid_pA);
              auto itB = vmap.find(pe.mid_pB);
              if (itA != vmap.end() && itB != vmap.end()) {
                auto avg = 0.5 * (itA->second + itB->second);
                fprintf(fp, "%s %s %s ",
                        isFinite(std::get<0>(avg)) ? NumtoString(std::get<0>(avg)).c_str() : "NaN",
                        isFinite(std::get<1>(avg)) ? NumtoString(std::get<1>(avg)).c_str() : "NaN",
                        isFinite(std::get<2>(avg)) ? NumtoString(std::get<2>(avg)).c_str() : "NaN");
              } else {
                fprintf(fp, "NaN NaN NaN ");
              }
            }
          }
          fprintf(fp, "\n</DataArray>\n");
        }
      }
      fprintf(fp, "</PointData>\n");
    }

    // Cells
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (int i = 0; i < num_points; ++i)
      fprintf(fp, "%d ", i);
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    { int sum = 0; for (const auto &sz : cell_sizes) fprintf(fp, "%d ", sum += sz); }
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (const auto &t : cell_types)
      fprintf(fp, "%d ", t);
    fprintf(fp, "\n</DataArray>\n");
    fprintf(fp, "</Cells>\n");

    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
#if defined(debug_mk_vtu)
    std::cout << Red << "|" << colorReset << std::endl;
#endif
    fclose(fp);
  } catch (std::exception &e) {
    std::cerr << e.what() << colorReset << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  }
}

enum class ShellPointLayer : int {
  Original = 0,
  Water = 1,
  Air = 2,
};

enum class ShellFaceKind : int {
  seed_interface = 0,
  water_outer = 1,
  air_outer = 2,
  internal = 3,
};

struct ShellVisualizationData {
  std::unique_ptr<Network> net;
  std::unordered_map<networkFace *, ShellFaceKind> face_kind;
  std::string warning;
  explicit operator bool() const noexcept { return static_cast<bool>(net); }
};

inline bool shell_precheck_normals(Network *net, std::string &warning) {
  for (auto *line : net->getBoundaryLines()) {
    if (!line || line->CORNER)
      continue;
    auto faces = line->getBoundaryFaces();
    if (faces.size() != 2 || !faces[0] || !faces[1])
      continue;
    const double dot = Dot(faces[0]->normal, faces[1]->normal);
    if (dot < -0.1) {
      warning = "skip shell visualization for " + net->getName() + ": reversed boundary-face normals detected on non-CORNER edge";
      return false;
    }
  }
  return true;
}

inline ShellVisualizationData build_shell_visualization(Network *seed_net) {
  ShellVisualizationData out;
  if (!seed_net)
    return out;

  auto seed_faces = seed_net->getBoundaryFaces();
  if (seed_faces.empty()) {
    out.warning = "skip shell visualization for " + seed_net->getName() + ": no boundary faces";
    return out;
  }

  if (!shell_precheck_normals(seed_net, out.warning))
    return out;

  std::unordered_map<networkPoint *, Tddd> normal_sum;
  std::unordered_map<networkPoint *, double> edge_sum;
  std::unordered_map<networkPoint *, int> edge_count;
  std::unordered_map<networkPoint *, double> min_altitude;
  normal_sum.reserve(seed_faces.size() * 3);
  edge_sum.reserve(seed_faces.size() * 3);
  edge_count.reserve(seed_faces.size() * 3);
  min_altitude.reserve(seed_faces.size() * 3);

  for (auto *face : seed_faces) {
    if (!face)
      continue;
    const auto [p0, p1, p2] = face->getPoints();
    if (!p0 || !p1 || !p2)
      continue;
    const double area = TriangleArea(p0->X, p1->X, p2->X);
    const Tddd normal = face->normal;
    if (!(area > 0.) || !isFinite(area) || !isFinite(normal)) {
      out.warning = "skip shell visualization for " + seed_net->getName() + ": invalid boundary-face geometry";
      return out;
    }

    normal_sum[p0] += area * normal;
    normal_sum[p1] += area * normal;
    normal_sum[p2] += area * normal;

    const double l01 = Norm(p1->X - p0->X);
    const double l12 = Norm(p2->X - p1->X);
    const double l20 = Norm(p0->X - p2->X);
    edge_sum[p0] += l01 + l20;
    edge_sum[p1] += l01 + l12;
    edge_sum[p2] += l12 + l20;
    edge_count[p0] += 2;
    edge_count[p1] += 2;
    edge_count[p2] += 2;

    const double eps = 1E-14;
    const double h0 = 2. * area / std::max(l12, eps);
    const double h1 = 2. * area / std::max(l20, eps);
    const double h2 = 2. * area / std::max(l01, eps);
    auto update_min_alt = [&](networkPoint *p, double h) {
      auto it = min_altitude.find(p);
      if (it == min_altitude.end())
        min_altitude.emplace(p, h);
      else
        it->second = std::min(it->second, h);
    };
    update_min_alt(p0, h0);
    update_min_alt(p1, h1);
    update_min_alt(p2, h2);
  }

  std::unordered_map<networkPoint *, Tddd> vertex_normal;
  std::unordered_map<networkPoint *, double> thickness;
  vertex_normal.reserve(normal_sum.size());
  thickness.reserve(normal_sum.size());
  for (const auto &[point, sum] : normal_sum) {
    const double nrm = Norm(sum);
    const auto it_count = edge_count.find(point);
    const auto it_edge = edge_sum.find(point);
    const auto it_alt = min_altitude.find(point);
    if (!(nrm > 0.) || !isFinite(nrm) || it_count == edge_count.end() || it_count->second <= 0 ||
        it_edge == edge_sum.end() || !(it_edge->second > 0.) ||
        it_alt == min_altitude.end() || !(it_alt->second > 0.)) {
      out.warning = "skip shell visualization for " + seed_net->getName() + ": invalid vertex-normal or thickness data";
      return out;
    }
    vertex_normal[point] = sum / nrm;
    const double mean_edge = it_edge->second / static_cast<double>(it_count->second);
    const double h = std::min(0.2 * mean_edge, 0.3 * it_alt->second);
    if (!(h > 0.) || !isFinite(h)) {
      out.warning = "skip shell visualization for " + seed_net->getName() + ": shell thickness collapsed";
      return out;
    }
    thickness[point] = h;
  }

  struct FaceShellPoints {
    std::array<int, 3> original;
    std::array<int, 3> water;
    std::array<int, 3> air;
  };

  std::vector<Tddd> coords;
  std::vector<FaceShellPoints> shell_points_per_face;
  coords.reserve(seed_faces.size() * 9);
  shell_points_per_face.reserve(seed_faces.size());

  for (auto *face : seed_faces) {
    const auto points = face->getPoints();
    FaceShellPoints ids;
    for (int i = 0; i < 3; ++i) {
      auto *p = points[i];
      const auto &n = vertex_normal.at(p);
      const double h = thickness.at(p);
      ids.original[i] = static_cast<int>(coords.size());
      coords.push_back(p->X);
      ids.water[i] = static_cast<int>(coords.size());
      coords.push_back(p->X - h * n);
      ids.air[i] = static_cast<int>(coords.size());
      coords.push_back(p->X + h * n);
    }
    shell_points_per_face.emplace_back(ids);
  }

  out.net = std::make_unique<Network>("file_name_is_not_given", seed_net->getName() + "_shell");
  auto shell_points = out.net->setPoints(coords);
  std::unordered_map<networkPoint *, ShellPointLayer> point_layer;
  point_layer.reserve(shell_points.size());

  auto shell_point = [&](int idx) -> networkPoint * {
    if (idx < 0 || idx >= static_cast<int>(shell_points.size()))
      return nullptr;
    return shell_points[idx];
  };

  for (const auto &ids : shell_points_per_face) {
    for (int i = 0; i < 3; ++i) {
      point_layer[shell_point(ids.original[i])] = ShellPointLayer::Original;
      point_layer[shell_point(ids.water[i])] = ShellPointLayer::Water;
      point_layer[shell_point(ids.air[i])] = ShellPointLayer::Air;
    }
  }

  auto append_prism = [&](const std::array<networkPoint *, 3> &base,
                          const std::array<networkPoint *, 3> &offset,
                          TetraState state) -> bool {
    const std::array<T_4P, 3> prism_tets = {{
        T_4P{base[0], base[1], base[2], offset[0]},
        T_4P{base[1], offset[1], offset[2], offset[0]},
        T_4P{base[1], base[2], offset[2], offset[0]},
    }};
    for (const auto &tet_points : prism_tets) {
      auto [created, tet] = genTetra(out.net.get(), tet_points);
      (void)created;
      if (!tet)
        return false;
      tet->tetra_state = state;
    }
    return true;
  };

  for (const auto &ids : shell_points_per_face) {
    const std::array<networkPoint *, 3> original = {shell_point(ids.original[0]), shell_point(ids.original[1]), shell_point(ids.original[2])};
    const std::array<networkPoint *, 3> water = {shell_point(ids.water[0]), shell_point(ids.water[1]), shell_point(ids.water[2])};
    const std::array<networkPoint *, 3> air = {shell_point(ids.air[0]), shell_point(ids.air[1]), shell_point(ids.air[2])};
    if (!append_prism(original, water, TetraState::Water) ||
        !append_prism(original, air, TetraState::Air)) {
      out.warning = "skip shell visualization for " + seed_net->getName() + ": failed to create shell tetrahedra";
      out.net.reset();
      return out;
    }
  }

  for (auto *face : out.net->getFaces()) {
    if (!face)
      continue;
    std::array<int, 3> counts = {0, 0, 0};
    for (auto *p : face->getPoints()) {
      auto it = point_layer.find(p);
      if (it != point_layer.end())
        counts[static_cast<int>(it->second)]++;
    }
    auto kind = ShellFaceKind::internal;
    if (counts[static_cast<int>(ShellPointLayer::Original)] == 3) {
      kind = ShellFaceKind::seed_interface;
      face->is_active_face = true;
    } else if (counts[static_cast<int>(ShellPointLayer::Water)] == 3) {
      kind = ShellFaceKind::water_outer;
      face->is_active_face = false;
    } else if (counts[static_cast<int>(ShellPointLayer::Air)] == 3) {
      kind = ShellFaceKind::air_outer;
      face->is_active_face = false;
    } else {
      face->is_active_face = false;
    }
    out.face_kind[face] = kind;
  }

  return out;
}

inline void write_shell_faces_vtu(const std::filesystem::path &path,
                                  const std::vector<networkFace *> &faces,
                                  const std::unordered_map<networkFace *, ShellFaceKind> &face_kind_map) {
  std::unordered_map<networkPoint *, int> point_index;
  std::vector<networkPoint *> ordered_points;
  ordered_points.reserve(faces.size() * 3);
  auto add_point = [&](networkPoint *p) {
    if (!p || point_index.contains(p))
      return;
    point_index.emplace(p, static_cast<int>(ordered_points.size()));
    ordered_points.emplace_back(p);
  };
  for (auto *face : faces)
    for (auto *p : face->getPoints())
      add_point(p);

  std::ofstream ofs(path);
  ofs << "<?xml version=\"1.0\"?>\n";
  ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  ofs << "<UnstructuredGrid>\n";
  ofs << "<Piece NumberOfPoints=\"" << ordered_points.size() << "\" NumberOfCells=\"" << faces.size() << "\">\n";
  ofs << "<Points>\n";
  ofs << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Position\" format=\"ascii\">\n";
  for (auto *p : ordered_points)
    ofs << NumtoString(p->X[0]) << " " << NumtoString(p->X[1]) << " " << NumtoString(p->X[2]) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "</Points>\n";
  ofs << "<CellData>\n";
  ofs << "<DataArray type=\"Int32\" Name=\"is_active_face\" format=\"ascii\">\n";
  for (auto *face : faces)
    ofs << (face && face->is_active_face ? 1 : 0) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray type=\"Int32\" Name=\"face_kind\" format=\"ascii\">\n";
  for (auto *face : faces) {
    const auto it = face_kind_map.find(face);
    const auto kind = (it == face_kind_map.end()) ? ShellFaceKind::internal : it->second;
    ofs << static_cast<int>(kind) << " ";
  }
  ofs << "\n</DataArray>\n";
  ofs << "</CellData>\n";
  ofs << "<Cells>\n";
  ofs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (auto *face : faces) {
    for (auto *p : face->getPoints())
      ofs << point_index.at(p) << " ";
  }
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < faces.size(); ++i)
    ofs << 3 * (i + 1) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < faces.size(); ++i)
    ofs << "5 ";
  ofs << "\n</DataArray>\n";
  ofs << "</Cells>\n";
  ofs << "</Piece>\n";
  ofs << "</UnstructuredGrid>\n";
  ofs << "</VTKFile>\n";
}

inline void write_shell_tets_vtu(const std::filesystem::path &path,
                                 const std::vector<networkTetra *> &tets) {
  std::unordered_map<networkPoint *, int> point_index;
  std::vector<networkPoint *> ordered_points;
  ordered_points.reserve(tets.size() * 4);
  auto add_point = [&](networkPoint *p) {
    if (!p || point_index.contains(p))
      return;
    point_index.emplace(p, static_cast<int>(ordered_points.size()));
    ordered_points.emplace_back(p);
  };
  for (auto *tet : tets)
    for (auto *p : tet->Points)
      add_point(p);

  std::ofstream ofs(path);
  ofs << "<?xml version=\"1.0\"?>\n";
  ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  ofs << "<UnstructuredGrid>\n";
  ofs << "<Piece NumberOfPoints=\"" << ordered_points.size() << "\" NumberOfCells=\"" << tets.size() << "\">\n";
  ofs << "<Points>\n";
  ofs << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Position\" format=\"ascii\">\n";
  for (auto *p : ordered_points)
    ofs << NumtoString(p->X[0]) << " " << NumtoString(p->X[1]) << " " << NumtoString(p->X[2]) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "</Points>\n";
  ofs << "<CellData>\n";
  ofs << "<DataArray type=\"Int32\" Name=\"tetra_state\" format=\"ascii\">\n";
  for (auto *tet : tets)
    ofs << static_cast<int>(tet ? tet->tetra_state : TetraState::Air) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "</CellData>\n";
  ofs << "<Cells>\n";
  ofs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (auto *tet : tets) {
    for (auto *p : tet->Points)
      ofs << point_index.at(p) << " ";
  }
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < tets.size(); ++i)
    ofs << 4 * (i + 1) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < tets.size(); ++i)
    ofs << "10 ";
  ofs << "\n</DataArray>\n";
  ofs << "</Cells>\n";
  ofs << "</Piece>\n";
  ofs << "</UnstructuredGrid>\n";
  ofs << "</VTKFile>\n";
}

inline std::vector<networkTetra *> filter_shell_tets_by_state(const std::vector<networkTetra *> &tets,
                                                              const TetraState state) {
  std::vector<networkTetra *> filtered;
  filtered.reserve(tets.size());
  for (auto *tet : tets)
    if (tet && tet->tetra_state == state)
      filtered.emplace_back(tet);
  return filtered;
}

inline void write_shell_step(const OutputContext &ctx,
                             const std::map<std::string, outputInfo> &NetOutputInfo,
                             const std::vector<Network *> &FluidObject) {
  for (auto *net : FluidObject) {
    auto shell = build_shell_visualization(net);
    if (!shell) {
      if (!shell.warning.empty())
        std::cout << Yellow << shell.warning << colorReset << std::endl;
      continue;
    }

    const auto face_key = net->getName() + "_shell_faces";
    const auto tet_key = net->getName() + "_shell_tets";
    const auto inner_tet_key = net->getName() + "_shell_inner_tets";
    const auto outer_tet_key = net->getName() + "_shell_outer_tets";
    const auto face_it = NetOutputInfo.find(face_key);
    const auto tet_it = NetOutputInfo.find(tet_key);
    const auto inner_tet_it = NetOutputInfo.find(inner_tet_key);
    const auto outer_tet_it = NetOutputInfo.find(outer_tet_key);
    const std::filesystem::path face_filename = ((face_it != NetOutputInfo.end()) ? face_it->second.vtu_file_name : face_key + "_") + std::to_string(ctx.time_step) + ".vtu";
    const std::filesystem::path tet_filename = ((tet_it != NetOutputInfo.end()) ? tet_it->second.vtu_file_name : tet_key + "_") + std::to_string(ctx.time_step) + ".vtu";
    const std::filesystem::path inner_tet_filename = ((inner_tet_it != NetOutputInfo.end()) ? inner_tet_it->second.vtu_file_name : inner_tet_key + "_") + std::to_string(ctx.time_step) + ".vtu";
    const std::filesystem::path outer_tet_filename = ((outer_tet_it != NetOutputInfo.end()) ? outer_tet_it->second.vtu_file_name : outer_tet_key + "_") + std::to_string(ctx.time_step) + ".vtu";

    std::vector<networkFace *> faces;
    std::vector<networkTetra *> tets;
    faces.reserve(shell.net->getFaces().size());
    tets.reserve(shell.net->getTetras().size());
    for (auto *face : shell.net->getFaces())
      faces.emplace_back(face);
    for (auto *tet : shell.net->getTetras())
      tets.emplace_back(tet);
    const auto inner_tets = filter_shell_tets_by_state(tets, TetraState::Water);
    const auto outer_tets = filter_shell_tets_by_state(tets, TetraState::Air);

    write_shell_faces_vtu(ctx.output_directory / face_filename, faces, shell.face_kind);
    write_shell_tets_vtu(ctx.output_directory / tet_filename, tets);
    write_shell_tets_vtu(ctx.output_directory / inner_tet_filename, inner_tets);
    write_shell_tets_vtu(ctx.output_directory / outer_tet_filename, outer_tets);

    if (face_it != NetOutputInfo.end() && face_it->second.PVD) {
      face_it->second.PVD->push(face_filename, ctx.simulation_time);
      face_it->second.PVD->output();
    }
    if (tet_it != NetOutputInfo.end() && tet_it->second.PVD) {
      tet_it->second.PVD->push(tet_filename, ctx.simulation_time);
      tet_it->second.PVD->output();
    }
    if (inner_tet_it != NetOutputInfo.end() && inner_tet_it->second.PVD) {
      inner_tet_it->second.PVD->push(inner_tet_filename, ctx.simulation_time);
      inner_tet_it->second.PVD->output();
    }
    if (outer_tet_it != NetOutputInfo.end() && outer_tet_it->second.PVD) {
      outer_tet_it->second.PVD->push(outer_tet_filename, ctx.simulation_time);
      outer_tet_it->second.PVD->output();
    }
  }
}

struct AFTCandidateTet {
  T4Tddd vertices;
  double edge_angle_deg = 0.;
  double volume_rel = 0.;
  double min_altitude_rel = 0.;
  int source_boundary_line_id = -1;
};

struct AFTCandidateStats {
  std::size_t total_edges = 0;
  std::size_t angle_pass = 0;
  std::size_t volume_reject = 0;
  std::size_t altitude_reject = 0;
  std::size_t intersection_reject = 0;
  std::size_t accepted = 0;
};

struct AFTCandidateVisualizationData {
  std::vector<AFTCandidateTet> accepted;
  AFTCandidateStats stats;
};

inline std::array<networkPoint *, 4> aft_points(networkLine *line, networkFace *f0, networkFace *f1) {
  auto [p0, p1] = line->getPoints();
  auto opposite = [](networkFace *face, networkPoint *a, networkPoint *b) -> networkPoint * {
    if (!face)
      return nullptr;
    auto [x0, x1, x2] = face->getPoints();
    for (auto *p : {x0, x1, x2})
      if (p != a && p != b)
        return p;
    return nullptr;
  };
  return {p0, p1, opposite(f0, p0, p1), opposite(f1, p0, p1)};
}

inline std::array<std::uintptr_t, 4> aft_candidate_key(const std::array<networkPoint *, 4> &points) {
  std::array<std::uintptr_t, 4> key = {
      reinterpret_cast<std::uintptr_t>(points[0]),
      reinterpret_cast<std::uintptr_t>(points[1]),
      reinterpret_cast<std::uintptr_t>(points[2]),
      reinterpret_cast<std::uintptr_t>(points[3])};
  std::sort(key.begin(), key.end());
  return key;
}

inline double aft_signed_volume(const T4Tddd &tet) {
  const auto &[a, b, c, d] = tet;
  return Dot(a - d, Cross(b - d, c - d)) / 6.;
}

inline double aft_max_edge_length(const T4Tddd &tet) {
  const auto &[a, b, c, d] = tet;
  return std::max({Norm(a - b), Norm(a - c), Norm(a - d), Norm(b - c), Norm(b - d), Norm(c - d)});
}

inline double aft_min_altitude(const T4Tddd &tet, const double abs_volume) {
  const auto &[a, b, c, d] = tet;
  const double eps = 1E-20;
  const std::array<double, 4> areas = {
      TriangleArea(b, c, d),
      TriangleArea(a, c, d),
      TriangleArea(a, b, d),
      TriangleArea(a, b, c)};
  double min_h = 1E+100;
  for (const auto area : areas) {
    if (!(area > eps) || !std::isfinite(area))
      return 0.;
    min_h = std::min(min_h, 3. * abs_volume / area);
  }
  return min_h;
}

inline bool aft_face_shares_any_candidate_vertex(networkFace *face, const std::array<networkPoint *, 4> &points) {
  if (!face)
    return false;
  auto [f0, f1, f2] = face->getPoints();
  for (auto *fp : {f0, f1, f2})
    for (auto *p : points)
      if (fp == p)
        return true;
  return false;
}

inline bool aft_new_face_intersects_boundary(Network *net,
                                             const T3Tddd &new_face,
                                             const CoordinateBounds &query_bounds,
                                             const double query_range,
                                             const std::array<networkPoint *, 4> &candidate_points,
                                             networkFace *seed_f0,
                                             networkFace *seed_f1) {
  std::unordered_set<networkFace *> nearby_faces;
  auto gather = [&](const Tddd &x) {
    net->BucketFaces.apply(x, query_range, [&](networkFace *face) {
      if (face)
        nearby_faces.emplace(face);
    });
  };
  gather(std::get<0>(new_face));
  gather(std::get<1>(new_face));
  gather(std::get<2>(new_face));
  gather(query_bounds.getCenter());

  for (auto *face : nearby_faces) {
    if (!face || face == seed_f0 || face == seed_f1)
      continue;
    if (aft_face_shares_any_candidate_vertex(face, candidate_points))
      continue;
    if (!IntersectQ(query_bounds.bounds, face->bounds))
      continue;
    auto [p0, p1, p2] = face->getPoints();
    try {
      if (IntersectionTriangles(new_face, T3Tddd{{p0->X, p1->X, p2->X}}).isIntersecting)
        return true;
    } catch (...) {
      return true;
    }
  }
  return false;
}

inline AFTCandidateVisualizationData build_aft_candidate_visualization(Network *seed_net) {
  AFTCandidateVisualizationData out;
  if (!seed_net)
    return out;

  const auto boundary_lines = seed_net->getBoundaryLines();
  if (boundary_lines.empty())
    return out;

  const double bucket_spacing = std::max(seed_net->getScale() / 20., 1E-6);
  seed_net->makeBucketFaces(bucket_spacing);

  std::set<std::array<std::uintptr_t, 4>> seen_candidates;
  int source_boundary_line_id = 0;
  for (auto *line : boundary_lines) {
    if (!line)
      continue;
    auto boundary_faces = line->getBoundaryFaces();
    if (boundary_faces.size() != 2 || !boundary_faces[0] || !boundary_faces[1])
      continue;

    ++out.stats.total_edges;
    const int current_line_id = source_boundary_line_id++;
    auto candidate_points = aft_points(line, boundary_faces[0], boundary_faces[1]);
    if (std::ranges::any_of(candidate_points, [](auto *p) { return p == nullptr; }))
      continue;

    const auto key = aft_candidate_key(candidate_points);
    if (!seen_candidates.emplace(key).second)
      continue;

    if (std::set<networkPoint *>{candidate_points.begin(), candidate_points.end()}.size() != 4)
      continue;

    const double dot = std::clamp(Dot(boundary_faces[0]->normal, boundary_faces[1]->normal), -1., 1.);
    const double alpha_deg = std::acos(dot) * 180. / M_PI;
    const double theta = 180. - alpha_deg;
    if (!(1. < theta && theta < 120.))
      continue;
    ++out.stats.angle_pass;

    T4Tddd tet = {candidate_points[0]->X, candidate_points[1]->X, candidate_points[2]->X, candidate_points[3]->X};
    double signed_volume = aft_signed_volume(tet);
    if (signed_volume < 0.) {
      std::swap(candidate_points[2], candidate_points[3]);
      tet = {candidate_points[0]->X, candidate_points[1]->X, candidate_points[2]->X, candidate_points[3]->X};
      signed_volume = -signed_volume;
    }
    const double L = aft_max_edge_length(tet);
    if (!(L > 0.) || !std::isfinite(L)) {
      ++out.stats.volume_reject;
      continue;
    }

    const double volume_rel = signed_volume / std::pow(L, 3);
    if (!(volume_rel > 1E-6) || !std::isfinite(volume_rel)) {
      ++out.stats.volume_reject;
      continue;
    }

    const double min_altitude = aft_min_altitude(tet, signed_volume);
    const double min_altitude_rel = min_altitude / L;
    if (!(min_altitude_rel >= 0.05) || !std::isfinite(min_altitude_rel)) {
      ++out.stats.altitude_reject;
      continue;
    }

    const T3Tddd new_face0 = {candidate_points[0]->X, candidate_points[2]->X, candidate_points[3]->X};
    const T3Tddd new_face1 = {candidate_points[1]->X, candidate_points[2]->X, candidate_points[3]->X};
    const CoordinateBounds bounds0(new_face0);
    const CoordinateBounds bounds1(new_face1);
    const double query_range = std::max({L, bounds0.getScale(), bounds1.getScale(), 1E-9});
    if (aft_new_face_intersects_boundary(seed_net, new_face0, bounds0, query_range, candidate_points, boundary_faces[0], boundary_faces[1]) ||
        aft_new_face_intersects_boundary(seed_net, new_face1, bounds1, query_range, candidate_points, boundary_faces[0], boundary_faces[1])) {
      ++out.stats.intersection_reject;
      continue;
    }

    out.accepted.push_back({tet, theta, volume_rel, min_altitude_rel, current_line_id});
  }

  out.stats.accepted = out.accepted.size();
  return out;
}

inline void write_aft_candidates_vtu(const std::filesystem::path &filename,
                                     const std::vector<AFTCandidateTet> &candidates) {
  std::ofstream ofs(filename);
  if (!ofs)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename.string() + " can not be opened");

  const auto num_cells = static_cast<int>(candidates.size());
  const auto num_points = 4 * num_cells;

  ofs << "<?xml version='1.0' encoding='UTF-8'?>\n";
  ofs << "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n";
  ofs << "<UnstructuredGrid>\n";
  ofs << "<Piece NumberOfCells='" << num_cells << "' NumberOfPoints='" << num_points << "'>\n";
  ofs << "<Points>\n";
  ofs << "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n";
  for (const auto &candidate : candidates)
    for (const auto &x : candidate.vertices)
      ofs << NumtoString(x[0]) << " " << NumtoString(x[1]) << " " << NumtoString(x[2]) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "</Points>\n";

  ofs << "<Cells>\n";
  ofs << "<DataArray type='Int32' Name='connectivity' format='ascii'>\n";
  for (int i = 0; i < num_points; ++i)
    ofs << i << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray type='Int32' Name='offsets' format='ascii'>\n";
  for (int i = 1; i <= num_cells; ++i)
    ofs << 4 * i << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray type='UInt8' Name='types' format='ascii'>\n";
  for (int i = 0; i < num_cells; ++i)
    ofs << "10 ";
  ofs << "\n</DataArray>\n";
  ofs << "</Cells>\n";

  ofs << "<CellData>\n";
  ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='edge_angle_deg' format='ascii'>\n";
  for (const auto &candidate : candidates)
    ofs << NumtoString(candidate.edge_angle_deg) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='volume_rel' format='ascii'>\n";
  for (const auto &candidate : candidates)
    ofs << NumtoString(candidate.volume_rel) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='min_altitude_rel' format='ascii'>\n";
  for (const auto &candidate : candidates)
    ofs << NumtoString(candidate.min_altitude_rel) << " ";
  ofs << "\n</DataArray>\n";
  ofs << "<DataArray NumberOfComponents='1' type='Int32' Name='source_boundary_line_id' format='ascii'>\n";
  for (const auto &candidate : candidates)
    ofs << candidate.source_boundary_line_id << " ";
  ofs << "\n</DataArray>\n";
  ofs << "</CellData>\n";

  ofs << "</Piece>\n";
  ofs << "</UnstructuredGrid>\n";
  ofs << "</VTKFile>\n";
}

inline void write_aft_candidates_step(const OutputContext &ctx,
                                      const std::map<std::string, outputInfo> &NetOutputInfo,
                                      const std::vector<Network *> &FluidObject) {
  for (auto *net : FluidObject) {
    if (!net)
      continue;
    auto aft = build_aft_candidate_visualization(net);
    const auto key = net->getName() + "_aft_candidates";
    const auto it = NetOutputInfo.find(key);
    const std::filesystem::path filename = ((it != NetOutputInfo.end()) ? it->second.vtu_file_name : key + "_") + std::to_string(ctx.time_step) + ".vtu";
    write_aft_candidates_vtu(ctx.output_directory / filename, aft.accepted);
    std::cout << Cyan << "[aft_debug] " << net->getName()
              << " total_edges=" << aft.stats.total_edges
              << " angle_pass=" << aft.stats.angle_pass
              << " volume_reject=" << aft.stats.volume_reject
              << " altitude_reject=" << aft.stats.altitude_reject
              << " intersection_reject=" << aft.stats.intersection_reject
              << " accepted=" << aft.stats.accepted
              << colorReset << std::endl;
    if (it != NetOutputInfo.end() && it->second.PVD) {
      it->second.PVD->push(filename, ctx.simulation_time);
      it->second.PVD->output();
    }
  }
}

void write_step(const OutputContext &ctx, const std::map<std::string, outputInfo> &NetOutputInfo, const std::vector<Network *> &FluidObject, const std::vector<Network *> &RigidBodyObject, const std::vector<Network *> &SoftBodyObject, const std::unordered_set<networkFace *> & /*allFaces*/) {
  // Fluid meshes
  for (auto *net : FluidObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu_quadratic(ctx.output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, ctx.dt));
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }
  }

  // Fluid tetrahedral meshes
  for (auto *net : FluidObject) {
    auto it = NetOutputInfo.find(net->getName() + "_tetra");
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    std::ofstream ofs(ctx.output_directory / filename);
    vtkUnstructuredGridWrite(ofs, net->getTetras());
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }
  }

  // Rigid meshes
  for (auto *net : RigidBodyObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu_quadratic(ctx.output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, ctx.dt));
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }

    // Mooring lines (VTP)
    for (auto *mooring : net->mooringLines) {
      auto itm = NetOutputInfo.find(mooring->getName());
      if (itm == NetOutputInfo.end())
        continue;
      const auto &minfo = itm->second;

      std::filesystem::path filename_m = minfo.vtu_file_name + std::to_string(ctx.time_step) + ".vtp";
      std::ofstream ofs(ctx.output_directory / filename_m);
      vtkPolygonWrite(ofs, mooring->getLines());
      ofs.close();
      if (minfo.PVD) {
        minfo.PVD->push(filename_m, ctx.simulation_time);
        minfo.PVD->output();
      }
    }
  }

  // Soft meshes
  for (auto *net : SoftBodyObject) {
    auto it = NetOutputInfo.find(net->getName());
    if (it == NetOutputInfo.end())
      continue;
    const auto &info = it->second;

    std::filesystem::path filename = info.vtu_file_name + std::to_string(ctx.time_step) + ".vtu";
    mk_vtu_quadratic(ctx.output_directory / filename, net->getBoundaryFaces(), dataForOutput(net, ctx.dt));
    if (info.PVD) {
      info.PVD->push(filename, ctx.simulation_time);
      info.PVD->output();
    }
  }
}

void write_vpm(const OutputContext &ctx, const VortexMethod &vpm, PVDWriter &pvd) {
  std::string filename = "vpm_" + std::to_string(ctx.time_step) + ".vtp";
  std::filesystem::path path = ctx.output_directory / filename;

  std::ofstream ofs(path);
  if (!ofs) {
    std::cerr << "Error: Cannot open file for writing: " << path << std::endl;
    return;
  }

  const auto &particles = vpm.getParticles();
  const size_t num_particles = particles.size();

  ofs << "<?xml version=\"1.0\"?>\n";
  ofs << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  ofs << "  <PolyData>\n";
  ofs << "    <Piece NumberOfPoints=\"" << num_particles << "\" NumberOfVerts=\"" << num_particles << "\">\n";

  // Points coordinates
  ofs << "      <Points>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.x[0] << " " << p.x[1] << " " << p.x[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "      </Points>\n";

  // Vertices (cells)
  ofs << "      <Verts>\n";
  ofs << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (size_t i = 0; i < num_particles; ++i) {
    ofs << "          " << i << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 0; i < num_particles; ++i) {
    ofs << "          " << i + 1 << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "      </Verts>\n";

  // PointData
  ofs << "      <PointData>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"Vorticity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.alpha[0] << " " << p.alpha[1] << " " << p.alpha[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"u_total\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.u_total[0] << " " << p.u_total[1] << " " << p.u_total[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"u_potential_BEM\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.u_potential_BEM[0] << " " << p.u_potential_BEM[1] << " " << p.u_potential_BEM[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"u_omega_VPM\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.u_omega_VPM[0] << " " << p.u_omega_VPM[1] << " " << p.u_omega_VPM[2] << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "        <DataArray type=\"Float64\" Name=\"Sigma\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for (const auto &p : particles) {
    ofs << "          " << p.sigma << "\n";
  }
  ofs << "        </DataArray>\n";
  ofs << "      </PointData>\n";

  ofs << "    </Piece>\n";
  ofs << "  </PolyData>\n";
  ofs << "</VTKFile>\n";

  ofs.close();

  pvd.push(filename, ctx.simulation_time);
  pvd.output();
}

} // namespace OutputParaView
