#ifndef my_vtk_H
#define my_vtk_H

#include "Network.hpp"

using map_P_Vd = std::map<netP *, V_d>;
using map_P_Vi = std::map<netP *, V_i>;

/*
VTK output entry points are currently split by responsibility.

- `vtkUnstructuredGridWriter` / `vtkPolygonWriter` in vtkWriter.hpp
  Generic typed writers. Prefer these when building new output code.
- `OutputParaView::mk_vtu_quadratic` in bem/OutputParaview.hpp
  Canonical BEM boundary-surface writer. It is quadratic-element aware and supports CellData.
- `mk_vtu(...)` overloads in this file
  Legacy convenience/debug helpers for quick linear VTU dumps of points, lines, and faces.

New feature work should usually go to the first two layers; the overloads below are kept for
existing call sites and lightweight debugging.
*/

/*DataArray_detail
networkPointをvtuとして出力するための関数．
引数variant内の`map`に，networkPointの値が欠けていたとしても，`NaN`を代わりに返す．
DataArray_detail */
/*DataArray_code*/
using V_str = std::vector<std::string>;

inline void mk_pvd(const std::string &filename, const std::map<std::string, double> &map_s_d) {
  Print("Creating " + filename + " ..." + colorReset);
  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  // check if file is opened
  if (fp == NULL)
    throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened\nMake sure you have the output directory!"));
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n");
  fprintf(fp, "  <Collection>\n");
  for (const auto &[file, time] : map_s_d) {
    fprintf(fp, "    <DataSet file=\"%s\" group=\"\" part=\"0\" timestep=\"%.15lf\"/>", file.c_str(), time);
    fprintf(fp, "\n");
  }
  fprintf(fp, "  </Collection>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
  Print("Done." + colorReset);
}

struct vtu_set {
  V_str vtu_names;
  V_d times;
  int index;
  vtu_set(const int index_IN) : index(index_IN) {};
  void push(const std::string &vtu_name, const double time) {
    this->vtu_names.emplace_back(vtu_name);
    this->times.emplace_back(time);
  };
};

class pvd {
public:
  std::map<std::string /*setname*/, vtu_set *> vtu_sets;
  V_str vtu_names;
  V_d times;
  std::string pvd_name;

  pvd(const std::string &pvd_name_IN = "") : vtu_names({}), times({}), pvd_name(pvd_name_IN), vtu_sets({}) {};

  void push(const std::string &vtu_name, const double time) {
    this->vtu_names.emplace_back(vtu_name);
    this->times.emplace_back(time);
  };

  void push(const std::string &name /*mapのkey,pvdの名前*/, const std::string &vtu_name, const double time) {
    if (this->vtu_sets.find(name) == vtu_sets.end()) {
      this->vtu_sets[name] = new vtu_set(this->vtu_sets.size());
      push(name, vtu_name, time);
    } else {
      (*(this->vtu_sets[name])).push(vtu_name, time);
    }
  };

  void output(const std::string &pvd_name_IN = "") const {
    std::string filename;
    if (!pvd_name_IN.empty())
      filename = pvd_name_IN;
    else if (!this->pvd_name.empty())
      filename = this->pvd_name;
    else
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "file name is not given"));

    Print("Creating " + filename + " ..." + colorReset);
    FILE *fp;
    fp = fopen(filename.c_str(), "wb");
    // check if file is opened
    if (fp == NULL)
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened\nMake sure you have the output directory!"));
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n");
    fprintf(fp, "  <Collection>\n");
    for (auto i = 0; i < times.size(); i++) {
      fprintf(fp, "    <DataSet file=\"%s\" group=\"\" part=\"0\" timestep=\"%.15lf\" volume=\"0.01\"/>", this->vtu_names[i].c_str(), this->times[i]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "  </Collection>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    Print("Done." + colorReset);
  };

  void output_(const std::string &pvd_name_IN = "") const {
    std::string filename;
    if (!pvd_name_IN.empty())
      filename = pvd_name_IN;
    else if (!this->pvd_name.empty())
      filename = this->pvd_name;
    else
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "file name is not given"));

    Print("Creating " + filename + " ..." + colorReset);
    FILE *fp;
    fp = fopen(filename.c_str(), "wb");
    // check if file is opened
    if (fp == NULL)
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened\nMake sure you have the output directory!"));
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n");
    fprintf(fp, "  <Collection>\n");
    for (const auto &[name, p_vtu] : this->vtu_sets) {
      for (auto i = 0; i < (*p_vtu).times.size(); i++) {
        fprintf(fp, "    <DataSet group=\"\" part=\"%d\" name=\"%s\" timestep=\"%.15lf\" file=\"%s\"/>", (*p_vtu).index, name.c_str(), (*p_vtu).times[i], (*p_vtu).vtu_names[i].c_str());
        fprintf(fp, "\n");
      }
    }
    fprintf(fp, "  </Collection>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    Print("Done." + colorReset);
  };
};

inline std::string NumtoString(const auto &x) {
  std::stringstream ss;
  ss << std::setprecision(8) << (float)(Between(x, {-1E-13, 1E-13}) ? 0. : x);
  return ss.str();
};

/* ------------------------------------------------------ */
#define debug_PVDWriter
struct PVDWriter {
  std::string pvdFileName;
  std::vector<std::tuple<std::string, double>> V_vtuFileName_time;
  //
  PVDWriter(const std::string &pvdFileNameIN) : pvdFileName(pvdFileNameIN) {};
  void push(const std::string &vtuFileName, const double time) { this->V_vtuFileName_time.emplace_back(std::make_tuple(vtuFileName, time)); };
  //
  void output(std::string filename = "") const {
    if (filename.empty())
      filename = this->pvdFileName;
#if defined(debug_PVDWriter)
    std::cout << Magenta << filename << std::endl;
#endif
    Print("Creating " + filename + " ..." + colorReset);
    FILE *fp;
    fp = fopen(filename.c_str(), "wb");
    // check if file is opened
    if (fp == NULL)
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened\nMake sure you have the output directory!"));
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n");
    fprintf(fp, "  <Collection>\n");
    for (const auto &name_time : V_vtuFileName_time)
      fprintf(fp, "    <DataSet file=\"%s\" group=\"\" part=\"0\" timestep=\"%.15lf\"/>\n", std::get<0>(name_time).c_str(), std::get<1>(name_time));
    fprintf(fp, "  </Collection>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    Print("Done." + colorReset);
  };
};
/* ------------------------------------------------------ */
template <class T> void writeDataArray(FILE *fp, const std::vector<T> &Points, const std::string &Name, const int NumberOfComponents, std::unordered_map<T, V_d> &map_p_v) {
  std::string str = "<DataArray NumberOfComponents='" + std::to_string(NumberOfComponents) + "' type='Float32' Name='" + Name + "' format='ascii'>\n";
  fprintf(fp, "%s", str.c_str());
  for (const auto &p : Points) {
    if (map_p_v.find(p) != map_p_v.end()) {
      auto vec = map_p_v[p];
      for (auto i = 0; i < NumberOfComponents; i++) {
        if (i < vec.size()) {
          if (vec[i] > 1E+15 || vec[i] < -1E+15)
            fprintf(fp, "NaN "); // 値が大きすぎるときは，NaNにする
          else
            fprintf(fp, "%s ", NumtoString(vec[i]).c_str());
        } else
          fprintf(fp, "NaN ");
      }
    } else
      for (auto i = 0; i < NumberOfComponents; i++)
        fprintf(fp, "NaN ");
    // fprintf(fp, "\n");
  }
  fprintf(fp, "\n</DataArray>\n");
};
// map Tddd
template <class T> void writeDataArray(FILE *fp, const std::vector<T> &Points, const std::string &Name, const std::unordered_map<T, Tddd> &map_p_v) {
  try {
    int NumberOfComponents = 3;
    std::string str = "<DataArray NumberOfComponents='" + std::to_string(NumberOfComponents) + "' type='Float32' Name='" + Name + "' format='ascii'>\n";
    fprintf(fp, "%s", str.c_str());
    for (const auto &p : Points) {
      if (map_p_v.find(p) != map_p_v.end()) {
        auto vec = map_p_v.at(p);
        if (!isFinite(std::get<0>(vec)))
          fprintf(fp, "NaN "); // 値が大きすぎるときは，NaNにする
        else
          fprintf(fp, "%s ", NumtoString(std::get<0>(vec)).c_str());
        if (!isFinite(std::get<1>(vec)))
          fprintf(fp, "NaN "); // 値が大きすぎるときは，NaNにする
        else
          fprintf(fp, "%s ", NumtoString(std::get<1>(vec)).c_str());
        if (!isFinite(std::get<2>(vec)))
          fprintf(fp, "NaN "); // 値が大きすぎるときは，NaNにする
        else
          fprintf(fp, "%s ", NumtoString(std::get<2>(vec)).c_str());
      } else
        for (auto i = 0; i < NumberOfComponents; i++)
          fprintf(fp, "NaN ");
      // fprintf(fp, "\n");
    }
    fprintf(fp, "\n</DataArray>\n");
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};
// map double
template <class T> void writeDataArray(FILE *fp, const std::vector<T> &Points, const std::string &Name, const std::unordered_map<T, double> &map_p_v) {
  try {
    int NumberOfComponents = 1;
    std::string str = "<DataArray NumberOfComponents='" + std::to_string(NumberOfComponents) + "' type='Float32' Name='" + Name + "' format='ascii'>\n";
    fprintf(fp, "%s", str.c_str());
    for (const auto &p : Points) {
      if (map_p_v.find(p) != map_p_v.end()) {
        auto value = map_p_v.at(p);
        if (!isFinite(value))
          fprintf(fp, "NaN "); // 値が大きすぎるときは，NaNにする
        else
          fprintf(fp, "%s ", NumtoString(value).c_str());
      } else
        for (auto i = 0; i < NumberOfComponents; i++)
          fprintf(fp, "NaN ");
      // fprintf(fp, "\n");
    }
    fprintf(fp, "\n</DataArray>\n");
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};
/* ------------------------------------------------------ */
// using uomap_P_Vd = std::unordered_map<networkPoint *, std::vector<double>>;
using uomap_P_Tddd = std::unordered_map<networkPoint *, Tddd>;
using uomap_P_d = std::unordered_map<networkPoint *, double>;
using uomap_F_d = std::unordered_map<networkFace *, double>;
// DOF-keyed maps: BEM_DOF_Base* accepts both networkPoint* and networkLine*
using uomap_DOF_Tddd = std::unordered_map<BEM_DOF_Base *, Tddd>;
using uomap_DOF_d = std::unordered_map<BEM_DOF_Base *, double>;
using VarForOutput = std::variant<std::string, int, uomap_P_Tddd, uomap_P_d, uomap_F_d, uomap_DOF_Tddd, uomap_DOF_d>;
using V_VarForOutput = std::vector<VarForOutput>;
using VV_VarForOutput = std::vector<V_VarForOutput>;
/* ------------------------------------------------------ */

template <class T> void DataArray(FILE *fp, const std::vector<T> &Points, const VV_VarForOutput &VV_name_comp_mapPVd) {
  try {
    for (auto &V_name_comp_mapPVd : VV_name_comp_mapPVd) {
      // V_name_comp_mapPVdの先頭は必ずstringでなければならない．
      // それ以降は，多様．
      std::string Name = std::get<std::string>(V_name_comp_mapPVd[0]);
      if (V_name_comp_mapPVd.size() > 1 && std::holds_alternative<std::unordered_map<T, double>>(V_name_comp_mapPVd[1])) {
        writeDataArray(fp, Points, Name, std::get<std::unordered_map<T, double>>(V_name_comp_mapPVd[1]));
      } else if (V_name_comp_mapPVd.size() > 0 && std::holds_alternative<std::unordered_map<T, Tddd>>(V_name_comp_mapPVd[1])) {
        writeDataArray(fp, Points, Name, std::get<std::unordered_map<T, Tddd>>(V_name_comp_mapPVd[1]));
      } else if (V_name_comp_mapPVd.size() > 1 && std::holds_alternative<uomap_F_d>(V_name_comp_mapPVd[1])) {
        // Generic point-based VTU writers do not support CellData. Face-based
        // diagnostics are emitted only by element-aware writers such as
        // OutputParaView::mk_vtu_quadratic().
        continue;
      } else {
        std::stringstream ss;
        ss << Name << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
    };
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  }
};
//-----------------------------------------
#define debug_mk_vtu
inline void mk_vtu2(const std::string &filename, const VV_d &VV_points) {
#if defined(debug_mk_vtu)
  std::cout << Magenta << filename << colorReset;
  std::cout << "  VV_points.size() : " << std::to_string(VV_points.size()) << colorReset << " ";
#endif
  auto Points = Flatten(VV_points);
  V_i PointsSize;
#if defined(debug_mk_vtu)
  std::cout << red << "|";
#endif
  for (const auto &V_points : VV_points) {
    PointsSize.push_back(V_points.size());
  }
#if defined(debug_mk_vtu)
  std::cout << red << "|";
#endif
  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened\nMake sure you have the output directory!");
  } else {
#if defined(debug_mk_vtu)
    std::cout << red << "|";
#endif
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    {
#if defined(debug_mk_vtu)
      std::cout << red << "|";
#endif
      fprintf(fp, "<UnstructuredGrid>\n");
      {
#if defined(debug_mk_vtu)
        std::cout << red << "|";
#endif
        fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", (int)VV_points.size(), (int)Points.size());
        {
#if defined(debug_mk_vtu)
          std::cout << red << "|";
#endif
          fprintf(fp, "<Points>\n");
          { // x,y,z: 3 components
            fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
            V_d xyz(3, 0.);
            // for (const auto &xyz : VV_points)
            // 	fprintf(fp, "%lf %lf %lf\n", xyz[0], xyz[1], xyz[2]);
            for (const auto &xyz : VV_points)
              fprintf(fp, "%s %s %s ", NumtoString(xyz[0]).c_str(), NumtoString(xyz[1]).c_str(), NumtoString(xyz[2]).c_str());
            fprintf(fp, "\n</DataArray>\n");
          }
          fprintf(fp, "</Points>\n");

// std::cout << red << "|";
// if (!VV_name_comp_mapPVd.empty())
// {
// 	fprintf(fp, "<PointData>\n");
// 	DataArray(fp, Points, VV_name_comp_mapPVd);
// 	fprintf(fp, "</PointData>\n");
// }
#if defined(debug_mk_vtu)
          std::cout << red << "|";
#endif
          fprintf(fp, "<Cells>\n");
          {
            fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
            int i(0);
            for (const auto &V_points : VV_points) {
              for (auto ii = 0; ii < V_points.size(); ++ii)
                fprintf(fp, "%d ", i++);
              fprintf(fp, "\n");
            }

            fprintf(fp, "</DataArray>\n");
            fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
            {
              int sum(0);
              for (const auto &i : PointsSize)
                fprintf(fp, "%d\n", sum = sum + i);
            }

            fprintf(fp, "</DataArray>\n");
            fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
            for (const auto &f : VV_points)
              fprintf(fp, "7 ");
            fprintf(fp, "\n</DataArray>\n");
          }
          fprintf(fp, "</Cells>\n");
        }
        fprintf(fp, "</Piece>\n");
      }
      fprintf(fp, "</UnstructuredGrid>\n");
    }
    fprintf(fp, "</VTKFile>\n");
#if defined(debug_mk_vtu)
    std::cout << Red << "|" << colorReset << std::endl;
#endif
  }
  fclose(fp);
}
/* ------------------------------------------------------ */
inline void mk_vtu(const std::string &filename, const std::vector<T4Tddd> &VV_points) {
  const int N = 4;
  try {
#if defined(debug_mk_vtu)
    std::cout << Magenta << filename << colorReset;
    std::cout << "  VV_points.size() : " << std::to_string(VV_points.size()) << colorReset << " ";
#endif
    V_netPp Points;
#if defined(debug_mk_vtu)
    std::cout << red << "|";
#endif
    FILE *fp;
    fp = fopen(filename.c_str(), "wb");
    if (fp == NULL)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened\nMake sure you have the output directory!");
    else {
#if defined(debug_mk_vtu)
      std::cout << red << "|";
#endif
      fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
      fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
      {
#if defined(debug_mk_vtu)
        std::cout << red << "|";
#endif
        fprintf(fp, "<UnstructuredGrid>\n");
        {
#if defined(debug_mk_vtu)
          std::cout << red << "|";
#endif
          fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", (int)VV_points.size(), N * (int)VV_points.size());
          {
#if defined(debug_mk_vtu)
            std::cout << red << "|";
#endif
            fprintf(fp, "<Points>\n");
            { // x,y,z: 3 components
              fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
              for (const auto &T8Tddds : VV_points) {
                auto [X0, X1, X2, X3] = T8Tddds;
                fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f ", std::get<0>(X0), std::get<1>(X0), std::get<2>(X0), std::get<0>(X1), std::get<1>(X1), std::get<2>(X1), std::get<0>(X2), std::get<1>(X2), std::get<2>(X2), std::get<0>(X3), std::get<1>(X3), std::get<2>(X3));
              }
              fprintf(fp, "\n</DataArray>\n");
            }
            fprintf(fp, "</Points>\n");

#if defined(debug_mk_vtu)
            std::cout << red << "|";
#endif
            fprintf(fp, "<Cells>\n");
            {
              fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
              for (auto i = 0; i < VV_points.size(); ++i)
                fprintf(fp, "%d %d %d %d ", N * i, N * i + 1, N * i + 2, N * i + 3);
              fprintf(fp, "\n</DataArray>\n");
              fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
              for (auto i = 0; i < VV_points.size(); ++i)
                fprintf(fp, "%d ", N * (i + 1));
              fprintf(fp, "\n</DataArray>\n");
              fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
              for (const auto &f : VV_points)
                fprintf(fp, "9 ");
              fprintf(fp, "\n</DataArray>\n");
            }
            fprintf(fp, "</Cells>\n");
          }
          fprintf(fp, "</Piece>\n");
        }
        fprintf(fp, "</UnstructuredGrid>\n");
      }
      fprintf(fp, "</VTKFile>\n");
#if defined(debug_mk_vtu)
      std::cout << Red << "|" << colorReset << std::endl;
#endif
    }
    fclose(fp);
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};
/* ------------------------------------------------------ */
inline void mk_vtu(const std::string &filename, const VV_netPp &VV_points, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  try {
#if defined(debug_mk_vtu)
    std::cout << Magenta << filename << colorReset;
    std::cout << "  VV_points.size() : " << std::to_string(VV_points.size()) << colorReset << " ";
#endif
    V_netPp Points;
    V_i PointsSize;
#if defined(debug_mk_vtu)
    std::cout << red << "|";
#endif
    for (const auto &V_points : VV_points) {
      PointsSize.push_back(V_points.size());
      for (const auto &p : V_points)
        Points.emplace_back(p);
    }
#if defined(debug_mk_vtu)
    std::cout << red << "|";
#endif
    FILE *fp;
    fp = fopen(filename.c_str(), "wb");
    if (fp == NULL) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, filename + " can not be opened\nMake sure you have the output directory!");
    } else {
#if defined(debug_mk_vtu)
      std::cout << red << "|";
#endif
      fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
      fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
      {
#if defined(debug_mk_vtu)
        std::cout << red << "|";
#endif
        fprintf(fp, "<UnstructuredGrid>\n");
        {
#if defined(debug_mk_vtu)
          std::cout << red << "|";
#endif
          fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", (int)VV_points.size(), (int)Points.size());
          {
#if defined(debug_mk_vtu)
            std::cout << red << "|";
#endif
            fprintf(fp, "<Points>\n");
            { // x,y,z: 3 components
              fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
              Tddd xyz;
              for (const auto &V_points : VV_points)
                for (const auto &p : V_points) {
                  xyz = p->getXtuple();
                  // fprintf(fp, "%f %f %f ", std::get<0>(xyz), std::get<1>(xyz), std::get<2>(xyz));
                  fprintf(fp, "%s %s %s ", NumtoString(std::get<0>(xyz)).c_str(), NumtoString(std::get<1>(xyz)).c_str(), NumtoString(std::get<2>(xyz)).c_str());
                }

              fprintf(fp, "\n</DataArray>\n");
            }
            fprintf(fp, "</Points>\n");

#if defined(debug_mk_vtu)
            std::cout << red << "|";
#endif
            if (!VV_name_comp_mapPVd.empty()) {
              fprintf(fp, "<PointData>\n");
              DataArray(fp, Points, VV_name_comp_mapPVd);
              fprintf(fp, "</PointData>\n");
            }
#if defined(debug_mk_vtu)
            std::cout << red << "|";
#endif
            fprintf(fp, "<Cells>\n");
            {
              fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
              int i(0);
              for (const auto &V_points : VV_points) {
                for (const auto &p : V_points)
                  fprintf(fp, "%d ", i++);
                // fprintf(fp, "\n");
              }

              fprintf(fp, "\n</DataArray>\n");
              fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
              {
                int sum(0);
                for (const auto &i : PointsSize)
                  fprintf(fp, "%d ", sum = sum + i);
              }

              fprintf(fp, "\n</DataArray>\n");
              fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
              for (auto ii = 0; ii < VV_points.size(); ++ii)
                fprintf(fp, "7 ");
              fprintf(fp, "\n</DataArray>\n");
            }
            fprintf(fp, "</Cells>\n");
          }
          fprintf(fp, "</Piece>\n");
        }
        fprintf(fp, "</UnstructuredGrid>\n");
      }
      fprintf(fp, "</VTKFile>\n");
#if defined(debug_mk_vtu)
      std::cout << Red << "|" << colorReset << std::endl;
#endif
    }
    fclose(fp);
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
}
/* ------------------------------------------------------ */
inline void mk_vtu(const std::string &filename, const std::vector<std::unordered_set<networkPoint *>> &V_points, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  VV_netPp ps;
  for (const auto &uo_ps : V_points) {
    ps.push_back({});
    for (const auto &p : uo_ps)
      ps.rbegin()->emplace_back(p);
  }
  mk_vtu(filename, ps, VV_name_comp_mapPVd);
}
/* ------------------------------------------------------ */
inline void mk_vtu(const std::string &filename, const V_netFp &Faces, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  VV_netPp fs;
  for (const auto &f : Faces)
    fs.emplace_back(ToVector(f->getPoints()));
  mk_vtu(filename, fs, VV_name_comp_mapPVd);
}
inline void mk_vtu(const std::string &filename, const std::unordered_set<networkFace *> &Faces, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  VV_netPp VVps(Faces.size());
  int i = 0;
  for (const auto &f : Faces)
    VVps[i++] = ToVector(f->getPoints());
  mk_vtu(filename, VVps, VV_name_comp_mapPVd);
}
/* ------------------------------------------------------ */
inline void mk_vtu(const std::string &filename, const V_netLp &Lines, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  VV_netPp vvp;
  for (const auto &l : Lines)
    vvp.emplace_back(ToVector(l->getPoints()));
  mk_vtu(filename, vvp, VV_name_comp_mapPVd);
}
/*DataArray_code*/

inline void mk_vtu(const std::string &filename, const VVV_d &vvv_IN, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  Network net;

  VV_netPp vvp({});
  for (const auto &vv : vvv_IN) {
    V_netPp ps({});
    for (const auto &v : vv)
      if (isFinite(v))
        ps.emplace_back(new networkPoint(&net, {v[0], v[1], v[2]}));
    vvp.emplace_back(ps);
  }

  mk_vtu(filename, vvp, VV_name_comp_mapPVd);

  for (const auto &vv : vvp)
    for (const auto &v : vv)
      delete v;
};
inline void mk_vtu(const std::string &filename, const std::vector<std::vector<Tddd>> &vvv_IN, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  Network net;

  VV_netPp vvp({});
  for (const auto &vv : vvv_IN) {
    V_netPp ps({});
    for (const auto &v : vv)
      if (isFinite(v))
        ps.emplace_back(new networkPoint(&net, {std::get<0>(v), std::get<1>(v), std::get<2>(v)}));
    vvp.emplace_back(ps);
  }

  mk_vtu(filename, vvp, VV_name_comp_mapPVd);

  for (const auto &vv : vvp)
    for (const auto &v : vv)
      delete v;
};

inline void mk_vtu(const std::string &filename, const std::vector<Tddd> &vv, const VV_VarForOutput &VV_name_comp_mapPVd = {}) {
  Network net;

  VV_netPp vvp({});
  V_netPp ps({});
  for (const auto &v : vv)
    if (isFinite(v))
      ps.emplace_back(new networkPoint(&net, {std::get<0>(v), std::get<1>(v), std::get<2>(v)}));
  vvp.emplace_back(ps);

  mk_vtu(filename, vvp, VV_name_comp_mapPVd);

  for (const auto &vv : vvp)
    for (const auto &v : vv)
      delete v;
};
/*DataArray_code*/

#endif
