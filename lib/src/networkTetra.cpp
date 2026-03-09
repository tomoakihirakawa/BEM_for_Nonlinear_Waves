// #pragma once
#include "Network.hpp"
#include "pch.hpp"

networkTetra::networkTetra(Network *network_IN, const T_4P &points, const T_6L &lines, const T_4F &faces) : Tetrahedron(ToX(points)), network(network_IN), Points(points), Lines(lines), Faces(faces) {
  for (const auto &f : faces)
    if (f == nullptr)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "A face of a tetrahedron cannot be nullptr.");

  // 1. 面に四面体を格納
  for (const auto &f : faces) {
    if (f->Tetras[0] == nullptr)
      f->Tetras[0] = this;
    else if (f->Tetras[1] == nullptr)
      f->Tetras[1] = this;
    else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "A face cannot belong to more than two tetrahedra.");
  }

  // 2. 線分に四面体を格納
  for (auto &l : this->Lines) {
    l->Tetras.emplace_back(this);
  }

  // 3. 頂点に四面体を格納
  for (auto &p : this->Points) {
    p->Tetras.emplace_back(this);
  }

  // 4. ネットワークに四面体を格納
  this->network->add(this);
};

/* -------------------------------------------------------------------------- */

networkTetra::~networkTetra() {

  // 1. FacesからのTetraの削除
  for (auto &f : this->Faces) {
    if (f) {
      if (f->Tetras[0] == this) {
        f->Tetras[0] = nullptr;
        // std::cout << "(face:" << f << ")->Tetras:" << f->Tetras << std::endl;
      } else if (f->Tetras[1] == this) {
        f->Tetras[1] = nullptr;
        // std::cout << "(face:" << f << ")->Tetras:" << f->Tetras << std::endl;
      }
    }
    f = nullptr;
  }

  // 2. LinesからのTetraの削除
  for (auto &l : this->Lines) {
    erase_element(l->Tetras, this); // l->Tetrasがvector/unordered_setかに応じて動作
    l = nullptr;
  }

  // 3. PointsからのTetraの削除
  for (auto &p : this->Points) {
    erase_element(p->Tetras, this); // p->Tetrasがvector/unordered_setかに応じて動作
    p = nullptr;
  }

  // 4. ネットワークからのTetra削除
  network->erase_element(this); // unordered_set対応版erase_element
  network = nullptr;

  //   std::cout << "networkTetra::~networkTetra() done\n";
}
