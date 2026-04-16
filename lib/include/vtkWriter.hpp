#ifndef vtkWriter_H
#define vtkWriter_H
#include "XML.hpp"
#include "basic_arithmetic_vector_operations.hpp"

inline Tddd ToX(const std::shared_ptr<Tddd> &X) { return *X; };

/*

vtuとvtpの違いは，PolysかCellsか
そして，vtkUnstructuredGridWriterには
```xml
<DataArray type="UInt8" Name="types" format="ascii">
10
</DataArray>
```
をCellsタグの中に追加するだけ．四面体（Tetrahedron）の場合はtypeが10、三角形（Triangle）の場合はtypeが5．

*/

template <typename T> struct vtkUnstructuredGridWriter : XMLElement {
  std::shared_ptr<XMLElement> UnstructuredGrid, Piece, Points, PointData, CellData, Cells, Lines;
  using VertexId = std::tuple<T, int>;
  std::unordered_map<VertexId, std::tuple<int, Tddd>> vertices;
  void reserve(const int &n) { this->vertices.reserve(n); };
  std::vector<std::vector<VertexId>> connectivity_lines;
  std::vector<std::array<VertexId, 3>> connectivity3;
  std::vector<std::array<VertexId, 4>> connectivity4;
  void clear() { vertices.clear(); };
  void add(const T &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); };
  void add(const T &v, const T &u) {
    std::get<1>(this->vertices[{v, 0}]) = ToX(v);
    std::get<1>(this->vertices[{u, 0}]) = ToX(u);
  };
  void add(const T &v, const T &u, const T &w) {
    std::get<1>(this->vertices[{v, 0}]) = ToX(v);
    std::get<1>(this->vertices[{u, 0}]) = ToX(u);
    std::get<1>(this->vertices[{w, 0}]) = ToX(w);
  };
  void add(const T &v, const T &u, const T &w, const T &x) {
    std::get<1>(this->vertices[{v, 0}]) = ToX(v);
    std::get<1>(this->vertices[{u, 0}]) = ToX(u);
    std::get<1>(this->vertices[{w, 0}]) = ToX(w);
    std::get<1>(this->vertices[{x, 0}]) = ToX(x);
  };
  void add(const std::vector<std::vector<T>> &verts) {
    for (const auto &V : verts)
      for (const auto &v : V)
        std::get<1>(this->vertices[{v, 0}]) = ToX(v);
  };
  template <std::size_t N> void add(const std::array<T, N> &verts) {
    std::ranges::for_each(verts, [&](const auto &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); });
  };

  void add(const std::vector<T> &verts) {
    std::ranges::for_each(verts, [&](const auto &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); });
  };
  void add(const std::unordered_set<T> &verts) {
    std::ranges::for_each(verts, [&](const auto &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); });
  };
  void addLine(const std::vector<T> &conns) {
    std::vector<VertexId> tmp(conns.size());
    for (auto i = 0; i < conns.size(); ++i)
      tmp[i] = {conns[i], 0};
    this->connectivity_lines.emplace_back(tmp);
  };
  // void addLine(const T &p0, const T &p1) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}}); };
  // void addLine(const T &p0, const T &p1, const T &p2) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}, {p2, 0}}); };
  // void addLine(const T &p0, const T &p1, const T &p2, const T &p3) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}, {p2, 0}, {p3, 0}}); };
  // Simplifying addLine methods
  template <typename... Args> void addLine(Args... args) {
    std::vector<VertexId> line = {{args, 0}...};
    this->connectivity_lines.emplace_back(std::move(line));
  }
  void addLine(const std::vector<std::vector<T>> &conns) { this->addLine(conns); };

  template <std::size_t N> void addGrid(const std::array<T, N> &conns) {
    static_assert(N == 3 || N == 4 || N == 8, "addGrid supports only 3, 4 and 8 element arrays");
    if constexpr (N == 3) {
      this->connectivity3.push_back({{{std::get<0>(conns), 0}, {std::get<1>(conns), 0}, {std::get<2>(conns), 0}}});
    } else if constexpr (N == 4) {
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<1>(conns), 0}, {std::get<2>(conns), 0}, {std::get<3>(conns), 0}}});
    } else if constexpr (N == 8) {
      /*
      cube case
      0----------1
      |\         |\
      | \        | \
      |  \       |  \
      |   4------|---5
      |   |      |   |
      3---|------2   |
       \  |       \  |
        \ |        \ |
         7----------6
      */
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<1>(conns), 0}, {std::get<2>(conns), 0}, {std::get<3>(conns), 0}}}); // 0-1-2-3
      this->connectivity4.push_back({{{std::get<4>(conns), 0}, {std::get<7>(conns), 0}, {std::get<6>(conns), 0}, {std::get<5>(conns), 0}}}); // 4-7-6-5
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<3>(conns), 0}, {std::get<7>(conns), 0}, {std::get<4>(conns), 0}}}); // 0-3-7-4
      this->connectivity4.push_back({{{std::get<1>(conns), 0}, {std::get<5>(conns), 0}, {std::get<6>(conns), 0}, {std::get<2>(conns), 0}}}); // 1-5-6-2
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<4>(conns), 0}, {std::get<5>(conns), 0}, {std::get<1>(conns), 0}}}); // 0-4-5-1
      this->connectivity4.push_back({{{std::get<3>(conns), 0}, {std::get<2>(conns), 0}, {std::get<6>(conns), 0}, {std::get<7>(conns), 0}}}); // 3-2-6-7
    }
  }

  template <std::size_t N> void addGrid(const std::vector<std::array<T, N>> &conns) {
    static_assert(N == 3 || N == 4 || N == 8, "addGrid supports only 3, 4 and 8 element arrays");
    for (const auto &c : conns)
      this->addGrid(c);
  }

  vtkUnstructuredGridWriter()
      : XMLElement("VTKFile", {{"type", "UnstructuredGrid"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}}), //<- constant setting
        UnstructuredGrid(new XMLElement("UnstructuredGrid")), Piece(new XMLElement("Piece", {{"NumberOfPoints", "4"}, {"NumberOfVerts", "0"}, {"NumberOfLines", "0"}, {"NumberOfStrips", "0"}, {"NumberOfCells", "4"}})), PointData(new XMLElement("PointData")), CellData(new XMLElement("CellData", {{"Scalars", "cell_scalars"}, {"Normals", "cell_normals"}})), Points(new XMLElement("Points")), Lines(new XMLElement("Lines")), Cells(new XMLElement("Cells")) {
    this->XMLElement::add(UnstructuredGrid);
    UnstructuredGrid->add(Piece);
    Piece->add(Points);
    Piece->add(PointData);
    Piece->add(CellData);
    Piece->add(Cells);
    Piece->add(Lines);
    /* -------------------------------------------------------------------------- */
    {
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"NumberOfComponents", "3"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        int i = 0;
        for (auto &[_, INT_X] : this->vertices) {
          // Incrementing index while assigning
          std::get<0>(INT_X) = i++;
          // Use std::transform for cleaner syntax
          std::transform(std::get<1>(INT_X).begin(), std::get<1>(INT_X).end(), std::ostream_iterator<float>(ofs, " "), [](const auto &x) {
            if (Between(x, {-1E-14, 1E-14}))
              return static_cast<float>(0);
            else
              return static_cast<float>(x);
          });
        }
      };
      this->Points->add(DataArray);
    }
    /* ------------------------------------------------------ */
    /*                        Polygons                        */
    /* ------------------------------------------------------ */
    { // connectivity
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        for (const auto &[ID0, ID1, ID2] : this->connectivity3)
          ofs << std::get<0>(vertices[ID0]) << " " << std::get<0>(vertices[ID1]) << " " << std::get<0>(vertices[ID2]) << " ";
        for (const auto &[ID0, ID1, ID2, ID3] : this->connectivity4)
          ofs << std::get<0>(vertices[ID0]) << " " << std::get<0>(vertices[ID1]) << " " << std::get<0>(vertices[ID2]) << " " << std::get<0>(vertices[ID3]) << " ";
      };
      this->Cells->add(DataArray);
    }
    { // offsets
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        int j = 0;
        for (auto i = 0; i < this->connectivity3.size(); ++i)
          ofs << (j += 3) << " ";
        for (auto i = 0; i < this->connectivity4.size(); ++i)
          ofs << (j += 4) << " ";
      };
      this->Cells->add(DataArray);
    }
    { // offsets
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "types"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        int j = 0;
        for (auto i = 0; i < this->connectivity3.size(); ++i)
          ofs << 5 << " ";
        for (auto i = 0; i < this->connectivity4.size(); ++i)
          ofs << 10 << " ";
      };
      this->Cells->add(DataArray);
    }
    /* ------------------------------------------------------ */
    /*                          Lines                         */
    /* ------------------------------------------------------ */
    { // connectivity
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        for (const auto &verts : this->connectivity_lines)
          for (const auto &ID : verts)
            ofs << std::get<0>(vertices[ID]) << " ";
      };
      this->Lines->add(DataArray);
    }
    { // offsets
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        int j = 0;
        for (auto i = 0; i < this->connectivity_lines.size(); ++i)
          ofs << (j += this->connectivity_lines[i].size()) << " ";
      };
      this->Lines->add(DataArray);
    }
  };
  /* ------------------------------------------------------ */
  std::unordered_map<std::string, std::unordered_map<T, Tddd>> data_Tddd;
  std::unordered_map<std::string, std::unordered_map<T, double>> data_double;
  // std::unordered_map<std::string, std::unordered_map<T, int>> data_int;
  void addPointData(const std::string &name, const std::unordered_map<T, Tddd> &V) {
    auto &data = data_Tddd[name];
    data.insert(V.cbegin(), V.cend());
    /* ------------------------------------------------------ */
    std::shared_ptr<XMLElement> DataArray = nullptr;
    for (const auto &elem : this->PointData->elements_shared)
      if (elem->attributes["Name"] == name) {
        DataArray = elem;
        // break;
      }
    /* ------------------------------------------------------ */
    if (DataArray == nullptr) {
      // std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", name}, {"format", "ascii"}, {"NumberOfComponents", "3"}}));
      DataArray = std::make_shared<XMLElement>("DataArray", std::map<std::string, std::string>{{"type", "Float32"}, {"Name", name}, {"format", "ascii"}, {"NumberOfComponents", "3"}});
      this->PointData->add(DataArray);
    }
    DataArray->writer = [&](std::stringstream &ofs) {
      auto it = data.begin();
      for (const auto &[ID, value] : this->vertices) {
        it = data.find(std::get<0>(ID));
        if (it != data.end()) {
          std::ranges::for_each(it->second, [&](const auto &x) {
            if (isFinite(x))
              ofs << std::setprecision(7) << (Between(x, {-1E-14, 1E-14}) ? 0 : (float)x) << " ";
            else
              ofs << std::setprecision(7) << "NaN ";
          });
        } else
          ofs << "NaN NaN NaN ";
      }
    };
  };
  void addPointData(const std::string &name, const std::unordered_map<T, double> &V) {
    auto &data = data_double[name];
    data.insert(V.cbegin(), V.cend());
    /* ------------------------------------------------------ */
    std::shared_ptr<XMLElement> DataArray = nullptr;
    for (const auto &elem : this->PointData->elements_shared)
      if (elem->attributes["Name"] == name) {
        DataArray = elem;
        // break;
      }
    /* ------------------------------------------------------ */
    if (DataArray == nullptr) {
      // std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", name}, {"format", "ascii"}}));
      DataArray = std::make_shared<XMLElement>("DataArray", std::map<std::string, std::string>{{"type", "Float32"}, {"Name", name}, {"format", "ascii"}});
      this->PointData->add(DataArray);
    }
    DataArray->writer = [&](std::stringstream &ofs) {
      auto it = data.begin();
      for (const auto &[ID, value] : this->vertices) {
        it = data.find(std::get<0>(ID));
        if (it != data.end()) {
          if (isFinite(it->second))
            ofs << std::setprecision(6) << (Between(it->second, {-1E-14, 1E-14}) ? 0 : (float)it->second) << " ";
          else
            ofs << std::setprecision(6) << "NaN ";
        } else
          ofs << "NaN ";
      }
    };
  };

  /* -------------------------------- addWithPosition -------------------------- */
  void addWithPosition(const T &v, const Tddd &position, int idx = 0) {
    std::get<1>(this->vertices[{v, idx}]) = position;
  }

  /* --------------------------------- addCellData --------------------------- */
  std::unordered_map<std::string, std::vector<double>> celldata_double;

  void addCellData(const std::string &name, const std::vector<double> &V) {
    celldata_double[name] = V;
    std::shared_ptr<XMLElement> DA = nullptr;
    for (const auto &elem : this->CellData->elements_shared)
      if (elem->attributes["Name"] == name)
        DA = elem;
    if (!DA) {
      DA = std::make_shared<XMLElement>("DataArray", std::map<std::string, std::string>{{"type", "Float32"}, {"Name", name}, {"format", "ascii"}});
      this->CellData->add(DA);
    }
    auto &data = celldata_double[name];
    DA->writer = [&](std::stringstream &ss) {
      for (const auto &val : data) {
        if (isFinite(val))
          ss << std::setprecision(6) << (Between(val, {-1E-14, 1E-14}) ? 0 : (float)val) << " ";
        else
          ss << "NaN ";
      }
    };
  }

  /* ---------------------------------- write --------------------------------- */

  template <typename V> V &write(V &ofs) const {
    ofs << "<?xml version=\"1.0\"?>\n";
    Piece->attributes["NumberOfPoints"] = std::to_string(this->vertices.size());
    Piece->attributes["NumberOfCells"] = std::to_string(this->connectivity3.size() + this->connectivity4.size());
    Piece->attributes["NumberOfLines"] = std::to_string(this->connectivity_lines.size());
    if (this->elements.empty() && this->elements_shared.empty()) {
      if (writer) {
        std::stringstream ss;
        writer(ss);
        ofs << this->getStartTag() << "\n" << ss.str() << "\n" << this->getEndTag() << "\n";
      } else
        ofs << this->getStartTag() << "\n" << this->getEndTag() << "\n";
    } else {
      ofs << this->getStartTag() << "\n";
      for (const auto &elem : this->elements)
        elem->write(ofs);
      for (const auto &elem : this->elements_shared)
        elem->write(ofs);
      ofs << this->getEndTag() << "\n";
    }
    return ofs;
  };
};

template <typename T> void vtkUnstructuredGridWrite(std::ofstream &ofs, const std::vector<std::array<T, 3>> &triangles, const std::vector<std::pair<std::string, std::unordered_map<T, double>>> &pointScalars = {}) {
  vtkUnstructuredGridWriter<T> vtu;
  vtu.reserve(triangles.size() * 3);
  for (const auto &tri : triangles) {
    vtu.add(tri);     // register points
    vtu.addGrid(tri); // push triangle connectivity (size=3)
  }
  for (const auto &[name, map] : pointScalars)
    vtu.addPointData(name, map);
  vtu.write(ofs);
}

template <typename T> void vtkUnstructuredGridWrite(std::ofstream &ofs, const std::vector<std::array<T, 4>> &tetras, const std::vector<std::pair<std::string, std::unordered_map<T, double>>> &pointScalars = {}) {
  vtkUnstructuredGridWriter<T> vtu;
  vtu.reserve(tetras.size() * 4);
  for (const auto &tet : tetras) {
    vtu.add(tet);
    vtu.addGrid(tet); // size=4 -> tetra
  }
  for (const auto &[name, map] : pointScalars)
    vtu.addPointData(name, map);
  vtu.write(ofs);
}

template <typename T> void vtkUnstructuredGridWrite(std::ofstream &ofs, const std::vector<T> &pointsOnly, const std::vector<std::pair<std::string, std::unordered_map<T, double>>> &pointScalars = {}) {
  vtkUnstructuredGridWriter<T> vtu;
  vtu.add(pointsOnly);
  for (const auto &[name, map] : pointScalars)
    vtu.addPointData(name, map);
  vtu.write(ofs);
}

// Overload using tuple (name, map) list with Tddd data (vector fields)
template <typename T> void vtkUnstructuredGridWrite(std::ofstream &ofs, const std::vector<std::array<T, 4>> &tetras, const std::vector<std::pair<std::string, std::unordered_map<T, Tddd>>> &pointVectors, const std::vector<std::pair<std::string, std::unordered_map<T, double>>> &pointScalars) {
  vtkUnstructuredGridWriter<T> vtu;
  for (const auto &tet : tetras) {
    vtu.add(tet);
    vtu.addGrid(tet);
  }
  for (const auto &[name, vmap] : pointVectors)
    vtu.addPointData(name, vmap);
  for (const auto &[name, smap] : pointScalars)
    vtu.addPointData(name, smap);
  vtu.write(ofs);
}

//@ -------------------------------------------------------------------------- */
//@ -------------------------------------------------------------------------- */
//@ -------------------------------------------------------------------------- */

template <typename T> struct vtkPolygonWriter : XMLElement {
  /*
  verticesはインデックスと座標ベクトルを含んでいる．
  その組み合わせをidとしている．インデックスは大抵は０．
  ポリゴンやラインを加える際は，このidの列をaddする．
  vtpでは座標ベクトルのインデックスを与えることで，ポリゴンを作るが，idが与えられれば，
  VertexId={インデックス，座標}なので，インデックスもわかるので，自動でインデックスを取得し出力してくれる．
  具体的な使い方については，octreeを見てください．
  */
  /* ------------------------------------------------------ */
  std::shared_ptr<XMLElement> PolyData, Piece, Points, PointData, CellData, Polys, Lines;
  using VertexId = std::tuple<T, int>;
  std::unordered_map<VertexId, std::tuple<int, Tddd>> vertices;
  void reserve(const int &n) { this->vertices.reserve(n); };
  std::vector<std::vector<VertexId>> connectivity_lines;
  std::vector<std::array<VertexId, 3>> connectivity3;
  std::vector<std::array<VertexId, 4>> connectivity4;
  void clear() { vertices.clear(); };
  void add(const T &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); };
  void add(const T &v, const T &u) {
    std::get<1>(this->vertices[{v, 0}]) = ToX(v);
    std::get<1>(this->vertices[{u, 0}]) = ToX(u);
  };
  void add(const T &v, const T &u, const T &w) {
    std::get<1>(this->vertices[{v, 0}]) = ToX(v);
    std::get<1>(this->vertices[{u, 0}]) = ToX(u);
    std::get<1>(this->vertices[{w, 0}]) = ToX(w);
  };
  void add(const T &v, const T &u, const T &w, const T &x) {
    std::get<1>(this->vertices[{v, 0}]) = ToX(v);
    std::get<1>(this->vertices[{u, 0}]) = ToX(u);
    std::get<1>(this->vertices[{w, 0}]) = ToX(w);
    std::get<1>(this->vertices[{x, 0}]) = ToX(x);
  };
  void add(const std::vector<std::vector<T>> &verts) {
    for (const auto &V : verts)
      for (const auto &v : V)
        std::get<1>(this->vertices[{v, 0}]) = ToX(v);
  };
  template <std::size_t N> void add(const std::array<T, N> &verts) {
    std::ranges::for_each(verts, [&](const auto &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); });
  };

  void add(const std::vector<T> &verts) {
    std::ranges::for_each(verts, [&](const auto &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); });
  };
  void add(const std::unordered_set<T> &verts) {
    std::ranges::for_each(verts, [&](const auto &v) { std::get<1>(this->vertices[{v, 0}]) = ToX(v); });
  };
  void addLine(const std::vector<T> &conns) {
    std::vector<VertexId> tmp(conns.size());
    for (auto i = 0; i < conns.size(); ++i)
      tmp[i] = {conns[i], 0};
    this->connectivity_lines.emplace_back(tmp);
  };
  // void addLine(const T &p0, const T &p1) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}}); };
  // void addLine(const T &p0, const T &p1, const T &p2) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}, {p2, 0}}); };
  // void addLine(const T &p0, const T &p1, const T &p2, const T &p3) { this->connectivity_lines.push_back({{p0, 0}, {p1, 0}, {p2, 0}, {p3, 0}}); };
  // Simplifying addLine methods
  template <typename... Args> void addLine(Args... args) {
    std::vector<VertexId> line = {{args, 0}...};
    this->connectivity_lines.emplace_back(std::move(line));
  }
  void addLine(const std::vector<std::vector<T>> &conns) { this->addLine(conns); };

  template <std::size_t N> void addPolygon(const std::array<T, N> &conns) {
    if constexpr (N == 3) {
      this->connectivity3.push_back({{{std::get<0>(conns), 0}, {std::get<1>(conns), 0}, {std::get<2>(conns), 0}}});
    } else if constexpr (N == 4) {
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<1>(conns), 0}, {std::get<2>(conns), 0}, {std::get<3>(conns), 0}}});
    } else if constexpr (N == 8) {
      /*
      cube case
      0----------1
      |\         |\
      | \        | \
      |  \       |  \
      |   4------|---5
      |   |      |   |
      3---|------2   |
       \  |       \  |
        \ |        \ |
         7----------6
      */
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<1>(conns), 0}, {std::get<2>(conns), 0}, {std::get<3>(conns), 0}}}); // 0-1-2-3
      this->connectivity4.push_back({{{std::get<4>(conns), 0}, {std::get<7>(conns), 0}, {std::get<6>(conns), 0}, {std::get<5>(conns), 0}}}); // 4-7-6-5
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<3>(conns), 0}, {std::get<7>(conns), 0}, {std::get<4>(conns), 0}}}); // 0-3-7-4
      this->connectivity4.push_back({{{std::get<1>(conns), 0}, {std::get<5>(conns), 0}, {std::get<6>(conns), 0}, {std::get<2>(conns), 0}}}); // 1-5-6-2
      this->connectivity4.push_back({{{std::get<0>(conns), 0}, {std::get<4>(conns), 0}, {std::get<5>(conns), 0}, {std::get<1>(conns), 0}}}); // 0-4-5-1
      this->connectivity4.push_back({{{std::get<3>(conns), 0}, {std::get<2>(conns), 0}, {std::get<6>(conns), 0}, {std::get<7>(conns), 0}}}); // 3-2-6-7
    } else {
      std::vector<VertexId> tmp(conns.size());
      for (auto i = 0; i < conns.size(); ++i)
        tmp[i] = {conns[i], 0};
      this->connectivity_lines.emplace_back(tmp);
      // static_assert(N == 3 || N == 4 || N == 8, "addPolygon supports only 3, 4 and 8 element arrays");
    }
  }

  template <std::size_t N> void addPolygon(const std::vector<std::array<T, N>> &conns) {
    // static_assert(N == 3 || N == 4 || N == 8, "addPolygon supports only 3, 4 and 8 element arrays");
    for (const auto &c : conns)
      this->addPolygon(c);
  }

  vtkPolygonWriter()
      : XMLElement("VTKFile", {{"type", "PolyData"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}}), //<- constant setting
        PolyData(new XMLElement("PolyData")), Piece(new XMLElement("Piece", {{"NumberOfPoints", "4"}, {"NumberOfVerts", "0"}, {"NumberOfLines", "0"}, {"NumberOfStrips", "0"}, {"NumberOfPolys", "4"}})), PointData(new XMLElement("PointData")), CellData(new XMLElement("CellData", {{"Scalars", "cell_scalars"}, {"Normals", "cell_normals"}})), Points(new XMLElement("Points")), Lines(new XMLElement("Lines")), Polys(new XMLElement("Polys")) {
    this->XMLElement::add(PolyData);
    PolyData->add(Piece);
    Piece->add(Points);
    Piece->add(PointData);
    Piece->add(CellData);
    Piece->add(Polys);
    Piece->add(Lines);
    /* -------------------------------------------------------------------------- */
    {
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"NumberOfComponents", "3"}, {"format", "ascii"}}));
      // DataArray->writer = [&](std::stringstream &ofs) {
      //    int i = 0;
      //    for (auto &[_, INT_X] : this->vertices) {
      //       // 書き出すともに番号を保存
      //       std::get<0>(INT_X) = i++;
      //       std::ranges::for_each(std::get<1>(INT_X), [&](auto &x) { ofs << std::setprecision(7) << (float)x << " "; });
      //    }
      // };
      // this->Points->add(DataArray);
      //
      // Using make_shared for efficiency and simplicity
      // auto DataArray = std::make_shared<XMLElement>("DataArray", {{"type", "Float32"}, {"NumberOfComponents", "3"}, {"format", "ascii"}});

      DataArray->writer = [&](std::stringstream &ofs) {
        int i = 0;
        for (auto &[_, INT_X] : this->vertices) {
          // Incrementing index while assigning
          std::get<0>(INT_X) = i++;

          // Use std::transform for cleaner syntax
          std::transform(std::get<1>(INT_X).begin(), std::get<1>(INT_X).end(), std::ostream_iterator<float>(ofs, " "), [](const auto &x) {
            if (Between(x, {-1E-14, 1E-14}))
              return static_cast<float>(0);
            else
              return static_cast<float>(x);
          });
        }
      };

      this->Points->add(DataArray);
    }
    /* ------------------------------------------------------ */
    /*                        Polygons                        */
    /* ------------------------------------------------------ */
    { // connectivity
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        for (const auto &[ID0, ID1, ID2] : this->connectivity3)
          ofs << std::get<0>(vertices[ID0]) << " " << std::get<0>(vertices[ID1]) << " " << std::get<0>(vertices[ID2]) << " ";
        for (const auto &[ID0, ID1, ID2, ID3] : this->connectivity4)
          ofs << std::get<0>(vertices[ID0]) << " " << std::get<0>(vertices[ID1]) << " " << std::get<0>(vertices[ID2]) << " " << std::get<0>(vertices[ID3]) << " ";
      };
      this->Polys->add(DataArray);
    }
    { // offsets
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        int j = 0;
        for (auto i = 0; i < this->connectivity3.size(); ++i)
          ofs << (j += 3) << " ";
        for (auto i = 0; i < this->connectivity4.size(); ++i)
          ofs << (j += 4) << " ";
      };
      this->Polys->add(DataArray);
    }
    /* ------------------------------------------------------ */
    /*                          Lines                         */
    /* ------------------------------------------------------ */
    { // connectivity
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "connectivity"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        for (const auto &verts : this->connectivity_lines)
          for (const auto &ID : verts)
            ofs << std::get<0>(vertices[ID]) << " ";
      };
      this->Lines->add(DataArray);
    }
    { // offsets
      std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Int32"}, {"Name", "offsets"}, {"format", "ascii"}}));
      DataArray->writer = [&](std::stringstream &ofs) {
        int j = 0;
        for (auto i = 0; i < this->connectivity_lines.size(); ++i)
          ofs << (j += this->connectivity_lines[i].size()) << " ";
      };
      this->Lines->add(DataArray);
    }
  };
  /* ------------------------------------------------------ */
  std::unordered_map<std::string, std::unordered_map<T, Tddd>> data_Tddd;
  std::unordered_map<std::string, std::unordered_map<T, double>> data_double;
  // std::unordered_map<std::string, std::unordered_map<T, int>> data_int;
  void addPointData(const std::string &name, const std::unordered_map<T, Tddd> &V) {
    auto &data = data_Tddd[name];
    data.insert(V.cbegin(), V.cend());
    /* ------------------------------------------------------ */
    std::shared_ptr<XMLElement> DataArray = nullptr;
    for (const auto &elem : this->PointData->elements_shared)
      if (elem->attributes["Name"] == name) {
        DataArray = elem;
        // break;
      }
    /* ------------------------------------------------------ */
    if (DataArray == nullptr) {
      // std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", name}, {"format", "ascii"}, {"NumberOfComponents", "3"}}));
      DataArray = std::make_shared<XMLElement>("DataArray", std::map<std::string, std::string>{{"type", "Float32"}, {"Name", name}, {"format", "ascii"}, {"NumberOfComponents", "3"}});
      this->PointData->add(DataArray);
    }
    DataArray->writer = [&](std::stringstream &ofs) {
      auto it = data.begin();
      for (const auto &[ID, value] : this->vertices) {
        it = data.find(std::get<0>(ID));
        if (it != data.end()) {
          std::ranges::for_each(it->second, [&](const auto &x) {
            if (isFinite(x))
              ofs << std::setprecision(6) << (Between(x, {-1E-14, 1E-14}) ? 0 : (float)x) << " ";
            else
              ofs << std::setprecision(6) << "NaN ";
          });
        } else
          ofs << "NaN NaN NaN ";
      }
    };
  };
  void addPointData(const std::string &name, const std::unordered_map<T, double> &V) {
    auto &data = data_double[name];
    data.insert(V.cbegin(), V.cend());
    /* ------------------------------------------------------ */
    std::shared_ptr<XMLElement> DataArray = nullptr;
    for (const auto &elem : this->PointData->elements_shared)
      if (elem->attributes["Name"] == name) {
        DataArray = elem;
        // break;
      }
    /* ------------------------------------------------------ */
    if (DataArray == nullptr) {
      // std::shared_ptr<XMLElement> DataArray(new XMLElement("DataArray", {{"type", "Float32"}, {"Name", name}, {"format", "ascii"}}));
      DataArray = std::make_shared<XMLElement>("DataArray", std::map<std::string, std::string>{{"type", "Float32"}, {"Name", name}, {"format", "ascii"}});
      this->PointData->add(DataArray);
    }
    DataArray->writer = [&](std::stringstream &ofs) {
      auto it = data.begin();
      for (const auto &[ID, value] : this->vertices) {
        it = data.find(std::get<0>(ID));
        if (it != data.end()) {
          if (isFinite(it->second))
            ofs << std::setprecision(6) << (Between(it->second, {-1E-14, 1E-14}) ? 0 : (float)it->second) << " ";
          else
            ofs << std::setprecision(6) << "NaN ";
        } else
          ofs << "NaN ";
      }
    };
  };

  /* --------------------------------- addCellData --------------------------- */
  // Writes one scalar per cell. The ordering follows the writer of Lines/Polys
  // connectivity: Polys first (Polys N=3, then N=4), then Lines. For a
  // pure-Lines use-case (cable segments), V must match connectivity_lines.
  std::unordered_map<std::string, std::vector<double>> celldata_double;
  void addCellData(const std::string &name, const std::vector<double> &V) {
    celldata_double[name] = V;
    std::shared_ptr<XMLElement> DA = nullptr;
    for (const auto &elem : this->CellData->elements_shared)
      if (elem->attributes["Name"] == name)
        DA = elem;
    if (!DA) {
      DA = std::make_shared<XMLElement>("DataArray", std::map<std::string, std::string>{{"type", "Float32"}, {"Name", name}, {"format", "ascii"}});
      this->CellData->add(DA);
    }
    auto &data = celldata_double[name];
    DA->writer = [&](std::stringstream &ss) {
      for (const auto &val : data) {
        if (isFinite(val))
          ss << std::setprecision(6) << (Between(val, {-1E-14, 1E-14}) ? 0 : (float)val) << " ";
        else
          ss << "NaN ";
      }
    };
  }

  /* ---------------------------------- write --------------------------------- */

  template <typename V> V &write(V &ofs) const {
    ofs << "<?xml version=\"1.0\"?>\n";
    Piece->attributes["NumberOfPoints"] = std::to_string(this->vertices.size());
    Piece->attributes["NumberOfPolys"] = std::to_string(this->connectivity3.size() + this->connectivity4.size());
    Piece->attributes["NumberOfLines"] = std::to_string(this->connectivity_lines.size());
    if (this->elements.empty() && this->elements_shared.empty()) {
      if (writer) {
        std::stringstream ss;
        writer(ss);
        ofs << this->getStartTag() << "\n" << ss.str() << "\n" << this->getEndTag() << "\n";
      } else
        ofs << this->getStartTag() << "\n" << this->getEndTag() << "\n";
    } else {
      ofs << this->getStartTag() << "\n";
      for (const auto &elem : this->elements)
        elem->write(ofs);
      for (const auto &elem : this->elements_shared)
        elem->write(ofs);
      ofs << this->getEndTag() << "\n";
    }
    return ofs;
  };
};

/* -------------------------------------------------------------------------- */

/* ------------------------------------------------------ */
/*                         T8Tddd                         */
/* ------------------------------------------------------ */

/*
cube case
0----------1
|\         |\
| \        | \
|  \       |  \
|   4------|---5
|   |      |   |
3---|------2   |
 \  |       \  |
  \ |        \ |
   7----------6
*/

template <size_t N> void vtkPolygonWrite(std::ofstream &ofs, const std::vector<std::array<Tddd, N>> &V) {
  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
  vtp.reserve(V.size());
  for (const auto &X : V) {
    std::array<std::shared_ptr<Tddd>, N> sharedPtrArray;
    for (size_t i = 0; i < N; ++i) {
      vtp.add(sharedPtrArray[i] = std::make_shared<Tddd>(X[i]));
    }
    vtp.addPolygon(sharedPtrArray);
  }
  vtp.write(ofs);
}

template <size_t N> void vtkPolygonWrite(std::ofstream &ofs, const std::array<Tddd, N> &V) {
  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
  vtp.reserve(V.size());
  std::array<std::shared_ptr<Tddd>, N> sharedPtrArray;
  for (size_t i = 0; i < N; ++i) {
    vtp.add(sharedPtrArray[i] = std::make_shared<Tddd>(V[i]));
  }
  vtp.addPolygon(sharedPtrArray);
  vtp.write(ofs);
}

// /* ------------------------------------------------------ */
// /*                         T4Tddd                         */
// /* ------------------------------------------------------ */
// void vtkPolygonWrite(std::ofstream &ofs, const std::vector<T4Tddd> &V) {
//    vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
//    for (const auto &[X0, X1, X2, X3] : V) {
//       std::shared_ptr<Tddd> x0(new Tddd(X0));
//       vtp.add(x0);
//       std::shared_ptr<Tddd> x1(new Tddd(X1));
//       vtp.add(x1);
//       std::shared_ptr<Tddd> x2(new Tddd(X2));
//       vtp.add(x2);
//       std::shared_ptr<Tddd> x3(new Tddd(X3));
//       vtp.add(x3);
//       vtp.addPolygon(std::array<std::shared_ptr<Tddd>, 4>{x0, x1, x2, x3});
//    }
//    vtp.write(ofs);
// };

// /* ------------------------------------------------------ */
// /*                         T3Tddd                         */
// /* ------------------------------------------------------ */
// void vtkPolygonWrite(std::ofstream &ofs, const std::vector<T3Tddd> &V) {
//    vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
//    for (const auto &[X0, X1, X2] : V) {
//       std::shared_ptr<Tddd> x0(new Tddd(X0));
//       vtp.add(x0);
//       std::shared_ptr<Tddd> x1(new Tddd(X1));
//       vtp.add(x1);
//       std::shared_ptr<Tddd> x2(new Tddd(X2));
//       vtp.add(x2);
//       vtp.addPolygon(std::array<std::shared_ptr<Tddd>, 3>{x0, x1, x2});
//    }
//    vtp.write(ofs);
// };

// /* ------------------------------------------------------ */
// /*                          Tddd                          */
// /* ------------------------------------------------------ */

inline void vtkPolygonWrite(std::ofstream &ofs, const std::vector<Tddd> &V) {
  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
  vtp.reserve(V.size());
  for (const auto &X : V) {
    std::shared_ptr<Tddd> x(new Tddd(X));
    vtp.add(x);
  }
  vtp.write(ofs);
};

// void vtkPolygonWrite(std::ofstream &ofs, const auto &V) {
//    vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
//    vtp.reserve(V.size());
//    for (const auto &X : V) {
//       std::shared_ptr<Tddd> x(new Tddd(ToX(X)));
//       vtp.add(x);
//    }
//    vtp.write(ofs);
// };

inline void vtkPolygonWrite(std::ofstream &ofs, const auto &V) {
  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
  vtp.reserve(V.size());
  for (const auto &PorF : V) {
    auto X = ToX(PorF);
    if constexpr (std::is_same_v<decltype(X), Tddd>) {
      auto x = std::make_shared<Tddd>(X);
      vtp.add(x);
    } else if constexpr (std::is_same_v<decltype(X), T2Tddd>) {
      auto [X0, X1] = X;
      auto x0 = std::make_shared<Tddd>(X0);
      auto x1 = std::make_shared<Tddd>(X1);
      vtp.add(x0);
      vtp.add(x1);
      vtp.addPolygon(std::array<std::shared_ptr<Tddd>, 2>{x0, x1});
    } else if constexpr (std::is_same_v<decltype(X), T3Tddd>) {
      auto [X0, X1, X2] = X;
      auto x0 = std::make_shared<Tddd>(X0);
      auto x1 = std::make_shared<Tddd>(X1);
      auto x2 = std::make_shared<Tddd>(X2);
      vtp.add(x0);
      vtp.add(x1);
      vtp.add(x2);
      vtp.addPolygon(std::array<std::shared_ptr<Tddd>, 3>{x0, x1, x2});
    }
  }
  vtp.write(ofs);
}

// template <typename U>
// void vtkPolygonWrite(std::ofstream &ofs, const auto &V) {
//    vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
//    vtp.reserve(V.size());
//    for (const auto &PorF : V) {
//       auto X = ToX(PorF);
//       if constexpr (std::is_same_v<decltype(X), Tddd>) {
//          auto x = std::make_shared<Tddd>(X);
//          vtp.add(x);
//       } else if constexpr (std::is_same_v<decltype(X), T2Tddd>) {
//          auto [X0, X1] = X;
//          auto x0 = std::make_shared<Tddd>(X0);
//          auto x1 = std::make_shared<Tddd>(X1);
//          vtp.add(x0);
//          vtp.add(x1);
//          vtp.addPolygon(std::array<std::shared_ptr<Tddd>, 2>{x0, x1});
//       } else if constexpr (std::is_same_v<decltype(X), T3Tddd>) {
//          auto [X0, X1, X2] = X;
//          auto x0 = std::make_shared<Tddd>(X0);
//          auto x1 = std::make_shared<Tddd>(X1);
//          auto x2 = std::make_shared<Tddd>(X2);
//          vtp.add(x0);
//          vtp.add(x1);
//          vtp.add(x2);
//          vtp.addPolygon(std::array<std::shared_ptr<Tddd>, 3>{x0, x1, x2});
//       }
//    }
//    vtp.write(ofs);
// }

/* ------------------------------------------------------ */
/*                         T2Tddd                         */
/* ------------------------------------------------------ */
inline void vtkPolygonWrite(std::ofstream &ofs, const std::vector<T2Tddd> &V) {
  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
  vtp.reserve(V.size());
  for (const auto &[X0, X1] : V) {
    std::shared_ptr<Tddd> x0(new Tddd(X0));
    vtp.add(x0);
    std::shared_ptr<Tddd> x1(new Tddd(X1));
    vtp.add(x1);
    vtp.addLine({x0, x1});
  }
  vtp.write(ofs);
};
// void vtkPolygonWrite(const std::string &name, const std::vector<T2Tddd> &V) {
//    std::ofstream ofs(name);
//    vtkPolygonWrite(ofs, V);
//    ofs.close();
// };
template <typename T, typename U> void vtkPolygonWrite(std::ofstream &ofs, const std::vector<std::tuple<T, T>> &V, const std::unordered_map<T, U> &data_double) {
  vtkPolygonWriter<T> vtp;
  vtp.reserve(V.size());
  for (const auto &[X0, X1] : V) {
    vtp.add(X0, X1);
    vtp.addLine(std::vector<T>{X0, X1});
  }
  vtp.addPointData("data", data_double);
  vtp.write(ofs);
};
/* ------------------------------------------------------ */
inline void vtkPolygonWrite(std::ofstream &ofs, const Tddd &X) {
  vtkPolygonWriter<std::shared_ptr<Tddd>> vtp;
  std::shared_ptr<Tddd> x(new Tddd(X));
  vtp.add(x);
  vtp.write(ofs);
};

// void vtkPolygonWrite(const std::string &name, const Tddd &V) {
//    std::ofstream ofs(name);
//    vtkPolygonWrite(ofs, V);
//    ofs.close();
// };

#include <concepts>
#include <type_traits>

template <class C>
concept PolygonPointContainer = requires(const C &c) {
  typename C::value_type;
  { c.begin() } -> std::input_iterator;
  { c.end() } -> std::input_iterator;
  { c.size() } -> std::convertible_to<std::size_t>;
};

// data マップ付き版
template <PolygonPointContainer C, class Map>
  requires std::same_as<typename C::value_type, typename Map::key_type>
void vtkPolygonWrite(std::ofstream &ofs, const C &V, const Map &data_map) {
  using T = typename C::value_type;
  vtkPolygonWriter<T> vtp;
  vtp.reserve(V.size());
  for (const auto &x : V)
    vtp.add(x);
  vtp.addPointData("data", data_map);
  vtp.write(ofs);
}

// template <typename T, typename U>
// void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<T> &V, const std::unordered_map<T, U> &data_double) {
//    vtkPolygonWriter<T> vtp;
//    vtp.reserve(V.size());
//    for (const auto &X : V)
//       vtp.add(X);
//    vtp.addPointData("data", data_double);
//    vtp.write(ofs);
// };

template <typename T, typename U> void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<T> &V, const std::vector<std::tuple<std::string, std::unordered_map<T, U>>> &name_uo_data = {}) {
  vtkPolygonWriter<T> vtp;
  for (const auto &X : V)
    vtp.add(X);

  for (const auto &[name, data] : name_uo_data)
    vtp.addPointData(name, data);
  vtp.write(ofs);
};

template <typename T, typename U> void vtkPolygonWrite(std::ofstream &ofs, const std::unordered_set<T> &V, const std::vector<std::pair<std::string, std::unordered_map<T, U>>> &name_uo_data = {}) {
  vtkPolygonWriter<T> vtp;
  for (const auto &X : V)
    vtp.add(X);

  for (const auto &[name, data] : name_uo_data)
    vtp.addPointData(name, data);
  vtp.write(ofs);
};

#endif