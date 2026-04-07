#pragma once

#include "BEM_midpoint_debug.hpp"
#include "BEM_node_face_state.hpp"
#include "BEM_pressure_detachment.hpp"
#include <cassert>
#include <unordered_set>

extern int time_step;
extern bool enable_pressure_detachment;
extern int detachment_consecutive_steps;

enum class LineBoundaryType : std::uint8_t {
  Undefined = 0,
  Dirichlet = 1,
  Neumann = 2,
  Corner = 3,
};

inline LineBoundaryType getEdgeNodeBoundaryTypeFromFaceStates(const networkLine* l) {
  if (!l)
    return LineBoundaryType::Undefined;

  const auto bfs = l->getBoundaryFaces();
  if (bfs.empty())
    return LineBoundaryType::Undefined;

  const bool all_neumann = std::ranges::all_of(bfs, [&](const auto& f) {
    return isNeumannBoundaryState(l, f);
  });
  if (all_neumann)
    return LineBoundaryType::Neumann;

  const bool all_dirichlet = std::ranges::all_of(bfs, [&](const auto& f) {
    return isDirichletBoundaryState(l, f);
  });
  if (all_dirichlet)
    return LineBoundaryType::Dirichlet;

  return LineBoundaryType::Corner;
}

inline NodeFaceBoundaryType getVertexNodeBoundaryTypeFromFaceStates(const networkPoint* p) {
  if (!p)
    return NodeFaceBoundaryType::Undefined;

  const auto bfs = p->getBoundaryFaces();
  if (bfs.empty())
    return NodeFaceBoundaryType::Undefined;

  const bool all_neumann = std::ranges::all_of(bfs, [&](const auto& f) {
    return isNeumannBoundaryState(p, f);
  });
  if (all_neumann)
    return NodeFaceBoundaryType::Neumann;

  const bool all_dirichlet = std::ranges::all_of(bfs, [&](const auto& f) {
    return isDirichletBoundaryState(p, f);
  });
  if (all_dirichlet)
    return NodeFaceBoundaryType::Dirichlet;

  return NodeFaceBoundaryType::Undefined;
}

inline bool hasAnyNeumannBoundaryState(const auto* entity) {
  const auto bfs = entity->getBoundaryFaces();
  return std::ranges::any_of(bfs, [&](const auto* f) {
    return getNodeFaceBoundaryType(entity, f) == NodeFaceBoundaryType::Neumann;
  });
}

inline bool hasAnyDirichletBoundaryState(const auto* entity) {
  const auto bfs = entity->getBoundaryFaces();
  return std::ranges::any_of(bfs, [&](const auto* f) {
    return getNodeFaceBoundaryType(entity, f) == NodeFaceBoundaryType::Dirichlet;
  });
}

template <class Entity>
inline void prunePressureDetachedFaceStates(Entity* entity, const std::unordered_set<networkFace*>& alive_faces) {
  if (!entity)
    return;
  // Prune detached_by_pressure in dofs for faces that are dead or lost contact
  for (auto it = entity->dofs.begin(); it != entity->dofs.end();) {
    auto* f = it->first;
    auto& d = it->second;
    if (!d.detached_by_pressure) {
      ++it;
      continue;
    }
    const bool face_alive = (f != nullptr) && alive_faces.count(f);
    const bool still_geometrically_in_contact = face_alive && d.nearestContactFace() != nullptr;
    if (!face_alive || !still_geometrically_in_contact)
      d.detached_by_pressure = false;
    ++it;
  }
}

// Pressure-detach flag is persistent while geometric contact remains.
// Here we only prune stale node-face states whose boundary face disappeared
// or whose geometric contact has already been lost.
inline void prunePressureDetachedTraceStates(Network* water) {
  const auto faces = water->getBoundaryFaces();
  const std::unordered_set<networkFace*> alive_faces(faces.begin(), faces.end());
  for (auto* p : water->getBoundaryPoints())
    prunePressureDetachedFaceStates(p, alive_faces);
  for (auto* l : water->getBoundaryLines())
    prunePressureDetachedFaceStates(l, alive_faces);
}

// face D/N is no longer stored as primary state. The only face-level quantity
// updated here is penetratedBody, which is still used as a geometric summary.
inline void updateFacePenetratedBodiesFromTraceStates(Network* water, const std::vector<Network*>& objects) {
  const auto faces = water->getBoundaryFaces();
  for (const auto& f : faces) {
    f->penetratedBody = nullptr;
    const bool all_point_faces_neumann = std::ranges::all_of(f->getPoints(), [&f](const auto& p) { return isNeumannBoundaryState(p, f); });
    if (!all_point_faces_neumann) {
      for (const auto& net : objects)
        if (net->InsideQ(f->centroid)) {
          f->penetratedBody = net;
          break;
        }
    }
  }
}

/*DOC_EXTRACT 0_1_1_2_BOUNDARY_CONDITIONS

## 多重節点

NOTE: 面の向き$\bf n$がカクッと不連続に変わる節点には，$\phi$は同じでも，隣接面にそれぞれ対して異なる$\phi_n$を計算できるようにする

NOTE: $\bf n$が不連続に変化する節点まわりの要素は，自分のために用意された$\phi_n$を選択し補間に用いなければならない

これを多重節点という．

多重節点を導入すると，未知変数idは，節点idだけではなく，節点と面の組みのidとなる．

### 境界値問題の未知変数ID 多重節点との区別

* `isNeumannBieDofKey`と`isDirichletBieDofKey`は，節点と面の組みが，境界値問題の未知変数かどうかを判定する．

* `pf2ID`は，節点と面の組みを未知変数IDに変換する．多重節点でない場合は，`{p,nullptr}`が変数のキーとなり，多重節点の場合は，与えられた`{p,f}`が変数のidとなる．

*/

/* ---------------------------------------------------------------------------
   Layer 2: BieDofKey — BIE unknown の代表キーとしてどの (entity, face) 組を
                        使うべきかを判定する（「規約」）

   境界条件判定の3層構造:

     Layer 1: BoundaryState   — 接触状態としてどうなっているか
                                （BEM_node_face_state.hpp: isNeumannBoundaryState 等）
     Layer 2: BieDofKey       — BIE unknown の代表キーとしてどの組を使うべきか（ここ）
     Layer 3: ActiveBieDof    — 実際に index >= 0 の DOF が存在するか
                                （Network.hpp: findActiveBieDof 等）

   利用可能タイミング:
     setContactFaces() 完了後。Layer 1 と同じ情報源を使うが、
     isMultipleNode と nullptr 規約を加味して DOF 割り当てルールを返す。
     dof.index は参照しない。「割り当てるべきか」であり「割り当てられているか」ではない。

   DOF 割り当て規約:
     - Dirichlet DOF は常に (entity, nullptr) に1つだけ割り当てられる。
       → isDirichletBieDofKey(entity, nullptr) == true
     - Neumann DOF は isMultipleNode の場合のみ per-face (entity, f) に割り当てられる。
       → isNeumannBieDofKey(entity, f) == true  (f != nullptr, isMultipleNode)
     - 非 isMultipleNode の Neumann は (entity, nullptr) に割り当てられる。
       → isNeumannBieDofKey(entity, nullptr) == true

   旧名: isNeumannID_BEM / isDirichletID_BEM
--------------------------------------------------------------------------- */
inline bool isNeumannBieDofKey(const auto p, const netF* f) {
  if (!p)
    return false;
  if (p->isMultipleNode)
    return (f != nullptr) && isNeumannBoundaryState(p, f);
  return (f == nullptr) && hasAnyNeumannBoundaryState(p);
};
inline bool isNeumannBieDofKey(const std::tuple<netP*, netF*>& PF) { return isNeumannBieDofKey(std::get<0>(PF), std::get<1>(PF)); };
// Dirichlet DOF は常に (entity, nullptr) に割り当てられる（per-face にならない）。
// そのため isNeumannBieDofKey のように特定 face の BoundaryState を問うのではなく、
// いずれかの face が Dirichlet かどうか（hasAnyDirichletBoundaryState）で判定する。
inline bool isDirichletBieDofKey(const auto p, const auto f) {
  if (!p)
    return false;
  return (f == nullptr) && hasAnyDirichletBoundaryState(p);
};
inline bool isDirichletBieDofKey(const std::tuple<netP*, netF*>& PF) { return isDirichletBieDofKey(std::get<0>(PF), std::get<1>(PF)); };

/*DOC_EXTRACT 0_3_BEM_utilities

## 多重節点を考慮したIDの設定方法

*/

//@ pf2IDは，setNodeFaceIndicesを実行せずとも使える．pf2Indexは，setNodeFaceIndicesを実行してから使う．
inline std::tuple<networkPoint*, networkFace*> pf2ID(const networkPoint* p, const networkFace* f) {
  if (isNeumannBieDofKey(p, f) || isDirichletBieDofKey(p, f))
    return {const_cast<networkPoint*>(p), const_cast<networkFace*>(f)};
  else
    return {const_cast<networkPoint*>(p), nullptr};
}

inline std::tuple<networkPoint*, networkFace*> pf2ID(const std::tuple<networkPoint*, networkFace*>& pf) { return pf2ID(std::get<0>(pf), std::get<1>(pf)); }

inline std::vector<std::tuple<networkPoint*, networkFace*>> p2AllIDs(const networkPoint* p) {
  std::vector<std::tuple<networkPoint*, networkFace*>> ret;
  bool nullptr_found = false;
  for (const auto& f : p->getBoundaryFaces()) {
    auto PF = pf2ID(p, f);
    if (nullptr_found && std::get<1>(PF) != nullptr)
      ret.emplace_back(PF);
    else if (!nullptr_found)
      ret.emplace_back(PF);
    if (std::get<1>(PF) == nullptr)
      nullptr_found = true;
  }
  return ret;
};

/*

   係数行列を作成する場合（LU分解など）：
   pf2Index(p0, integ_f)は，積分の重みとに掛かる節点上のある量を指定のに使われる．
   これは，係数行列を作成する際に使うことになる．
   多重節点の場合でも適切にIDを返す．

   係数行列を作成する必要がない場合（GMRESなど）：
   もし，p->f2Index.at(f)が存在するなら，それは多重節点として扱われる．その値を使う．
   つまり，

*/

inline int pf2Index(const networkPoint* p, networkFace* f) {
  auto lookup = [&](networkFace* key) -> int {
    const auto* d = p->findActiveBieDof(key);
    if (d)
      return d->index;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "pf2Index: dof not found");
  };
  if (f == nullptr || !p->isMultipleNode || getVertexNodeFaceBoundaryType(p, f) == NodeFaceBoundaryType::Dirichlet)
    return lookup(nullptr);
  else
    return lookup(f);
};

/*

   pf2Index(p0, integ_f)を使えるように刷るためには，setNodeFaceIndicesを実行する必要がある．

*/

// --- 辺中点の多重節点ID/Index関数（頂点の pf2ID / pf2Index に対応） ---
// isNeumannBieDofKey / isDirichletBieDofKey は auto 化済みのため networkLine にも使用可能

inline std::tuple<networkLine*, networkFace*> lf2ID(const networkLine* l, const networkFace* f) {
  if (isNeumannBieDofKey(l, f) || isDirichletBieDofKey(l, f))
    return {const_cast<networkLine*>(l), const_cast<networkFace*>(f)};
  return {const_cast<networkLine*>(l), nullptr};
}

inline int lf2Index(const networkLine* l, networkFace* f) {
  auto lookup = [&](networkFace* key) -> int {
    const auto* d = l->findActiveBieDof(key);
    if (d)
      return d->index;
    return -1; // 非活性中点（hybrid境界辺など）
  };
  if (f == nullptr || !l->isMultipleNode || getEdgeNodeFaceBoundaryType(l, f) == NodeFaceBoundaryType::Dirichlet)
    return lookup(nullptr);
  return lookup(f);
}

// Returns 6 DOF indices for a true quadratic face: [p0, p1, p2, l0_mid, l1_mid, l2_mid]
// Topology: l0=edge(p0,p1), l1=edge(p1,p2), l2=edge(p2,p0)
// Matches TriShape<6> node ordering: N0..N2=vertices, N3=mid(p0,p1), N4=mid(p1,p2), N5=mid(p2,p0)
inline std::array<int, 6> getQuadDOFIndices(const networkFace* f) {
  const auto& [p0, p1, p2] = f->getPoints();
  const auto& [l0, l1, l2] = f->Lines;
  return {pf2Index(p0, const_cast<networkFace*>(f)),
          pf2Index(p1, const_cast<networkFace*>(f)),
          pf2Index(p2, const_cast<networkFace*>(f)),
          lf2Index(l0, const_cast<networkFace*>(f)),
          lf2Index(l1, const_cast<networkFace*>(f)),
          lf2Index(l2, const_cast<networkFace*>(f))};
}

/* -------------------------------------------------------------------------- */

inline std::size_t setNodeFaceIndices(const std::vector<Network*>& objects) {
  // Indexing policy:
  // - Assign indices in a spatially coherent order (wave-front from top corner).
  // - Improves locality for sparse/ILU preconditioners.

  // ---- helpers ----

  // Reference point: (min_x, min_y, max_z) — top-left-front corner
  double min_x = 1e30, min_y = 1e30, max_z = -1e30;
  for (const auto* water : objects)
    for (const auto* p : water->getBoundaryPoints()) {
      min_x = std::min(min_x, std::get<0>(p->X));
      min_y = std::min(min_y, std::get<1>(p->X));
      max_z = std::max(max_z, std::get<2>(p->X));
    }
  const Tddd ref_p = {min_x, min_y, max_z};

  auto dist_sq = [&](const Tddd& X) {
    double dx = std::get<0>(X) - ref_p[0], dy = std::get<1>(X) - ref_p[1], dz = std::get<2>(X) - ref_p[2];
    return dx * dx + dy * dy + dz * dz;
  };

  auto sort_faces = [&](auto faces) {
    std::stable_sort(faces.begin(), faces.end(), [&](const networkFace* a, const networkFace* b) { return dist_sq(a->centroid) < dist_sq(b->centroid); });
    return faces;
  };

  // Prune stale dof entries (faces deleted by remesh)
  auto pruneStaleDofs = [](auto* entity, const auto& alive) {
    for (auto it = entity->dofs.begin(); it != entity->dofs.end();)
      (it->first && !alive.count(it->first)) ? it = entity->dofs.erase(it) : ++it;
  };

  // ---- alive faces set ----
  std::unordered_set<networkFace*> alive_faces;
  for (const auto* water : objects)
    for (auto* f : water->getBoundaryFaces())
      alive_faces.emplace(f);

  // ---- collect points ----
  std::unordered_set<networkPoint*> unique_points;
  for (const auto* water : objects)
    for (auto* p : water->getBoundaryPoints())
      unique_points.emplace(p);
  std::vector<networkPoint*> points(unique_points.begin(), unique_points.end());
  std::stable_sort(points.begin(), points.end(), [&](const auto* a, const auto* b) { return dist_sq(a->X) < dist_sq(b->X); });

  for (auto* p : points)
    pruneStaleDofs(p, alive_faces);
  for (auto* p : points)
    setMultipleNode(p);

  // ---- assign DOF indices (common for point and line) ----
  auto assignDofIndices = [&](auto* entity, std::size_t& idx) {
    for (auto& [ff, dd] : entity->dofs)
      dd.index = -1;

    if (isNeumannBieDofKey(entity, nullptr) || isDirichletBieDofKey(entity, nullptr))
      entity->dof(nullptr).index = idx++;

    if (entity->isMultipleNode)
      for (auto* f : sort_faces(entity->getBoundaryFaces()))
        if (f && isNeumannBieDofKey(entity, f))
          entity->dof(f).index = idx++;
  };

  std::size_t i = 0;
  for (auto* p : points)
    assignDofIndices(p, i);

  // ---- line midpoints (true quadratic only) ----
  if (use_true_quadratic_element) {
    std::unordered_set<networkLine*> unique_lines;
    for (const auto* water : objects)
      for (auto* l : water->getBoundaryLines())
        unique_lines.emplace(l);
    std::vector<networkLine*> lines(unique_lines.begin(), unique_lines.end());
    auto midpoint_pos = [](const networkLine* l) { auto [pA, pB] = l->getPoints(); return 0.5 * (pA->X + pB->X); };
    std::stable_sort(lines.begin(), lines.end(), [&](const auto* a, const auto* b) {
      return dist_sq(midpoint_pos(a)) < dist_sq(midpoint_pos(b));
    });

    for (auto* l : lines)
      pruneStaleDofs(l, alive_faces);
    for (auto* l : lines)
      setMultipleNode(l);

    for (auto* l : lines) {
      // hybrid: skip lines not touching any quadratic face
      if (use_quadratic_linear_hybrid) {
        auto bf = l->getBoundaryFaces();
        if (!std::ranges::any_of(bf, [](const auto* f) { return f && f->isTrueQuadraticElement; }))
          continue;
        if (!std::ranges::all_of(bf, [](const auto* f) { return f && f->isTrueQuadraticElement; })) {
          for (auto& [ff, dd] : l->dofs)
            dd.index = -1;
          l->midpoint_index = -1;
          continue;
        }
      }

      assignDofIndices(l, i);
      l->midpoint_index = (l->dof(nullptr).index >= 0) ? l->dof(nullptr).index : -1;
    }
  }

  // ---- verify consistency ----
  auto checkConsistency = [](const auto* entity, const char* msg) {
    bool has_multiple = std::ranges::count_if(entity->dofs, [](const auto& kv) { return kv.second.index >= 0; }) > 1;
    assert(entity->isMultipleNode == has_multiple && msg);
  };
  for (auto* p : points)
    checkConsistency(p, "isMultipleNode inconsistent for point");
  for (const auto* water : objects)
    for (auto* l : water->getBoundaryLines())
      checkConsistency(l, "isMultipleNode inconsistent for line");

  return i;
};

inline std::size_t setNodeFaceIndices(const Network* objects) { return setNodeFaceIndices(std::vector<Network*>{const_cast<Network*>(objects)}); };

/* -------------------------------------------------------------------------- */
/*  initializeNodeFaceStates                                                  */
/*  全 DOF carrier の boundary (node, face) に NodeFaceState エントリを保証。  */
/*  try_emplace: 既存エントリ（detach state 等）は上書きしない。               */
/* -------------------------------------------------------------------------- */

inline void initializeNodeFaceStates(Network* water) {
  auto ensureDofs = [](auto* node) {
    node->dofs.try_emplace(nullptr);
    for (auto* f : node->getBoundaryFaces())
      node->dofs.try_emplace(f);
  };
  for (auto* p : water->getBoundaryPoints())
    ensureDofs(p);
  for (auto* l : water->getBoundaryLines())
    ensureDofs(l);
}

/* -------------------------------------------------------------------------- */

#include "BEM_BoundaryValues.hpp"

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_1_1_1_BOUNDARY_CONDITIONS

## 境界のタイプを決定する

<img src="./img/schematic_boundary_types_without_float.png" width="600px">

0. 流体と物体の衝突を判定し，流体節点が接触する物体面を保存しておく．

   * \ref{contact_angle}{`networkPoint::contact_angle`}
   * \ref{isInContact}{`networkPoint::isInContact`}
   * \ref{addContactFaces}{`networkPoint::addContactFaces`}

を使って接触判定を行っている．

 \ref{BEM:contact_range}{流体が構造物との接触を感知する半径}の設置も重要．

つぎに，その情報を使って，境界のタイプを次の順で決める．（物理量を与えるわけではない）

1. 面の境界条件：３節点全てが接触している流体面はNeumann面，それ以外はDirichlet面とする．CORNER面は設定しない．
   - Neumann面$\Gamma^{({\rm N})}$ : 3点接触流体面
   - Dirichlet面$\Gamma^{({\rm D})}$ : それ以外の面

2. 辺の境界条件 : 辺を含む２面がNeumann面ならNeumann辺，２面がDirichlet面ならDirichlet辺，それ以外はCORNERとする．
   - Neumann辺 : 隣接面2面がNeumann面の辺
   - Dirichlet辺 : 隣接面2面がDirichlet面の辺
   - CORNER辺 : それ以外の辺（Neumann面とDirichlet面の間にある辺）

3. 点の境界条件：点を含む面全てがNeumann面ならNeumann点，面全てがDirichlet面ならDirichlet点，それ以外はCORNERとする．
   - Neumann点 : 隣接面全てがNeumann面である点
   - Dirichlet点 : 隣接面全てがDirichlet面である点
   - CORNER点 : それ以外の点（Neumann面とDirichlet面の間にある点）

*/

inline void setRigidBodyVelocityAndAccel(Network* net, const double& RK_time) {
  // この関数は，オブジェクトの重心と姿勢に関するものであり，物体の変形には関係しない．
  // この結果，net->velocityとnet->accelerationが設定される．

  // net->velocity : 物体の重心と姿勢の速度ベクトル（6成分）
  // net->acceleration : 物体の重心と姿勢の加速度
  //! 注意 物体の変形は，節点毎にこれに加算される形で設定される．

  if (std::ranges::all_of(net->isFixed, [](const auto& b) { return b; })) {
    net->mass = 1E+20;
    net->inertia.fill(1E+20);
    net->COM.fill(0.);
    net->initial_center_of_mass.fill(0.);
  } else {
    for (auto i = 0; i < net->isFixed.size(); i++) {
      if (net->isFixed[i])
        net->inertia[i] = 1E+20;
    }
  }

  // MOIが設定されておらず、かつ回転が固定されていない浮遊剛体の場合、
  // 慣性モーメントが0.0のままとなり計算が破綻するのを防ぐためのチェックと警告。
  if (net->isFloatingBody) {
    for (int i = 3; i < 6; ++i) {
      if (net->isFixed[i] == false && net->inertia[i] == 0.0) {
        std::cerr << "Warning: MOI for rigid body '" << net->getName() << "' (DOF " << i << ") is zero, but it's a floating body. This may cause instability. Please specify a non-zero MOI in the input file." << std::endl;
      }
    }
  }

  std::string move_name_velocity;
  T6d default_acceleration = {0., 0., 0., 0., 0., 0.};

  if (net->inputJSON.find("velocity")) {
    move_name_velocity = normalizeVelocityName(net->inputJSON["velocity"][0]);
    std::cout << "move_name_velocity = " << move_name_velocity << std::endl;
    if (move_name_velocity == "update")
      std::cout << " velocity is already updated using acceleration" << std::endl;
    else if (move_name_velocity == "fixed")
      net->velocity.fill(0.);
    else if (move_name_velocity == "floating") {
      std::cout << "floatingの場合は，加速度の時間積分によってシミュレートされる" << std::endl;
      net->velocity = net->RK_Velocity.getX();
    } else {
      std::cout << "(RigidBodyObject) velocity is explicityly given as " << move_name_velocity << std::endl;
      if (!isSupportedRigidBodyVelocityName(move_name_velocity))
        throw error_message(__FILE__,
                            __PRETTY_FUNCTION__,
                            __LINE__,
                            "unknown rigid-body velocity profile '" + move_name_velocity +
                                "'. Supported names: " + supportedRigidBodyVelocityNames());
      double delta_t = 1E-5;
      if (move_name_velocity == "file") {
        double start_time = 0.;
        if (net->inputJSON["velocity"].size() >= 3)
          start_time = std::stod(net->inputJSON["velocity"][2]);
        net->velocity = net->intpMotionRigidBody.D(RK_time + start_time);
        default_acceleration = net->intpMotionRigidBody(RK_time + delta_t / 2.) - net->intpMotionRigidBody(RK_time - delta_t / 2.);
      } else {
        net->velocity = velocity(move_name_velocity, net->inputJSON["velocity"], RK_time);
        default_acceleration = velocity(move_name_velocity, net->inputJSON["velocity"], RK_time + delta_t / 2.) - velocity(move_name_velocity, net->inputJSON["velocity"], RK_time - delta_t / 2.);
      }
      std::cout << "net->velocity = " << net->velocity << std::endl;
      default_acceleration /= delta_t;
    }
  } else {
    std::cout << "指定がないので速度はゼロ" << std::endl;
    net->velocity.fill(0.);
  }

  std::cout << "setting acceleration" << std::endl;
  std::string move_name_accel;
  if (move_name_velocity == "fixed")
    net->acceleration.fill(0.);
  else if (move_name_velocity == "floating")
    std::cout << "floatingの場合は，加速度は計算する" << std::endl;
  else if (net->inputJSON.find("acceleration")) {
    move_name_accel = net->inputJSON["acceleration"][0];
    if (move_name_accel == "fixed")
      net->acceleration.fill(0.);
    else if (move_name_accel == "floating") {
      std::cout << "floatingの場合は，加速度の時間積分によってシミュレートされる" << std::endl;
      // この時点ではわからない
    } else
      net->acceleration = acceleration(move_name_accel, net->inputJSON["acceleration"], RK_time);
  } else {
    std::cout << "指定がないので加速度はdefault_acceleration" << std::endl;
    net->acceleration = default_acceleration;
  }

  for (const auto& p : net->getPoints()) {
    auto tmp = net->velocityRigidBody(p->X);
    p->velocity[0] = tmp[0];
    p->velocity[1] = tmp[1];
    p->velocity[2] = tmp[2];
  }
};

// b# ------------------------------------------------------ */
// b#      物体のノイマン境界の速度 u(t) at Neumann を設定        */
// b# ------------------------------------------------------ */

//\label{BEM:setBodyVelocity}
inline void setBodyVelocity(const std::vector<Network*>& objects) {
  for (auto net : objects) {
    std::cout << Green << "setBodyVelocity: " << colorReset << net->getName() << std::endl;
    //! 壁面の動きは，マイステップ更新することにした．この結果はphin()で参照される
    net->velocity.fill(0.);
    // net->acceleration.fill(0.);
    for (const auto& p : net->getPoints()) {
      p->velocity.fill(0.);
      // p->acceleration.fill(0.);
    }

    // 考え方：
    // RigidBody内部の点の動きは，物体の重心と姿勢の動きから決定される．
    // SoftBody内部の点の動きは，物体の重心と姿勢の動きだけではなく，各節点に与えられた速度relative_velocityが加算される．
    // relative_velocityの定義は，物体の初期位置に対しての相対速度として定義される．
    // ただ，relative_velocityが与えられたかどうかで，isRigidBodyとisSoftBodyは，区別されるので，現段階では，RigidBodyとSoftBodyの両方に対して同じコードで対応できる．
    // SoftBodyの場合も，net->RK_COMとnet->RK_Qを使う．

    // 1. まずは，net->RK_COMとnet->RK_Qを各ルンゲクッタの時刻に更新する．
    auto RK_time = net->RK_COM.gett(); //%各ルンゲクッタの時刻を使う
    setRigidBodyVelocityAndAccel(net, RK_time);

    // 2. 各節点の速度を設定する．（RigidBody，SoftBodyの両方に対応）
    if (net->inputJSON.find("relative_velocity")) {
      std::string move_name = net->inputJSON["relative_velocity"][0];
      std::cout << "move_name = " << move_name << std::endl;
      for (const auto& p : net->getPoints()) {
        // floatingの場合は，velocityRigidBody内部で参照されるnet->RK_COMとnet->RK_Qが各ルンゲクッタの時刻における物体の速度・姿勢を保持しているので，問題ないはず．
        Tddd V = velocityOnBody(net, p, RK_time);
        p->velocity[0] = V[0];
        p->velocity[1] = V[1];
        p->velocity[2] = V[2];
      }
    }
  }
}

// b# ------------------------------------------------------ */

/* -------------------------------------------------------------------------- */
/*                             f,l,pの境界条件を決定                             */
/* -------------------------------------------------------------------------- */

inline void setBoundaryTypes(Network* water, const std::vector<Network*>& objects) {
  std::cout << water->getName() << "の境界条件を決定 setBoundaryTypes" << std::endl;
  // 前提: initializeNodeFaceStates() + setContactFaces() 済み（refreshBoundaryStatesAndTypes 経由）

  for (const auto& f : water->getBoundaryFaces())
    f->penetratedBody = nullptr;

  for (const auto& p : water->getBoundaryPoints()) {
    p->penetratedBody = nullptr;
    for (const auto& net : objects)
      if (net->InsideQ(p->X))
        p->penetratedBody = net;
  }

  std::cout << "step2 節点-面状態と面の貫入情報を更新" << std::endl;
  const auto faces = water->getBoundaryFaces();
  prunePressureDetachedTraceStates(water);
  updateFacePenetratedBodiesFromTraceStates(water, objects);

  // Pressure-based detachment is applied per node-face state, i.e. (p,f) and
  // (l,f), not per face. Once detached, the node-face state stays Dirichlet
  // until geometric contact is lost and the stale detach flag gets pruned above.
  if (enable_pressure_detachment && time_step > 0) {
    int detached_count = 0;
    auto tryDetach = [&](auto* node) {
      for (auto* f : node->getBoundaryFaces()) {
        if (isNeumannBoundaryState(node, f) &&
            !node->penetratedBody &&
            pressureDetachmentEligible(node, f) &&
            node->dof(f).detach_negative_count >= detachment_consecutive_steps) {
          node->dof(f).detached_by_pressure = true;
          node->dof(f).detach_negative_count = 0;
          ++detached_count;
        }
      }
    };
    for (auto* p : water->getBoundaryPoints())
      tryDetach(p);
    for (auto* l : water->getBoundaryLines())
      tryDetach(l);
    updateFacePenetratedBodiesFromTraceStates(water, objects);
    if (detached_count > 0)
      std::cout << Magenta << "[pressure_detach] " << detached_count
                << " node-face states detached" << colorReset << std::endl;
  }

  std::cout << "step3 線の境界条件を決定" << std::endl;
  for (const auto& l : water->getLines()) {

    bool was_neumann = l->Neumann;
    bool was_corner = l->CORNER;

    l->Neumann = false;
    l->Dirichlet = false;
    l->CORNER = false;

    const auto line_bc = getEdgeNodeBoundaryTypeFromFaceStates(l);
    l->Neumann = (line_bc == LineBoundaryType::Neumann);
    l->CORNER = (line_bc == LineBoundaryType::Corner);
    l->Dirichlet = (line_bc == LineBoundaryType::Dirichlet);

    if (!l->Neumann && !l->CORNER && !l->Dirichlet)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "line with mixed or undefined boundary condition");

    // 境界属性が変わった場合、端点の CORNER 点の曲率キャッシュを invalidate
    // （CORNER 点は l->Neumann でフィルタして近傍を収集するため）
    if (l->Neumann != was_neumann || l->CORNER != was_corner) {
      auto [p0, p1] = l->getPoints();
      if (p0 && p0->CORNER) p0->geom_curvature.valid = false;
      if (p1 && p1->CORNER) p1->geom_curvature.valid = false;
    }
  }

  if (const char* env = std::getenv("BEM_LINE_NORMAL_DEBUG"); env && std::string(env) != "0") {
    std::size_t two_face_lines = 0;
    std::size_t noncorner_two_face_lines = 0;
    std::size_t flipped_noncorner = 0;
    double min_dot_noncorner = 1.0;
    for (const auto& l : water->getBoundaryLines()) {
      const auto bfs = l->getBoundaryFaces();
      if (bfs.size() != 2)
        continue;
      ++two_face_lines;
      const double d = Dot(bfs[0]->normal, bfs[1]->normal);
      if (!l->CORNER) {
        ++noncorner_two_face_lines;
        min_dot_noncorner = std::min(min_dot_noncorner, d);
        if (d < 0.0)
          ++flipped_noncorner;
      }
    }
    std::cout << Magenta << "[BEM:line-normal] " << Cyan
              << "two_face=" << Green << two_face_lines
              << Cyan << " noncorner_two_face=" << Green << noncorner_two_face_lines
              << Cyan << " flipped_noncorner(dot<0)=" << Red << flipped_noncorner
              << Cyan << " min_dot_noncorner=" << Yellow << min_dot_noncorner
              << colorReset << std::endl;
  }

  // CORNER line classification is now known, so re-build line-side contact
  // candidates using endpoint-derived check faces. This is used to detect
  // perpendicular supporting body surfaces around the waterline.
  {
    std::vector<networkLine*> corner_lines;
    for (auto* l : water->getBoundaryLines())
      if (l->CORNER)
        corner_lines.push_back(l);
    for (const auto& l : corner_lines)
      l->clearContactFaces();

    _Pragma("omp parallel for") for (const auto& l : corner_lines)
        l->addContactFaces(objects, false);

    // Debug: check CORNER-line contact candidates after the re-run.
    int total_corner = corner_lines.size();
    int with_contacts = 0;
    for (const auto& l : corner_lines) {
      auto cf = l->getContactFaces();
      if (!cf.empty())
        ++with_contacts;
    }
    std::cout << "[CORNER line re-run] total=" << total_corner
              << " with_contacts=" << with_contacts
              << " check_faces_sample=";
    if (!corner_lines.empty()) {
      auto sample = corner_lines[0]->getFacesForContactCheck();
      std::cout << sample.size() << " (own=" << corner_lines[0]->getBoundaryFaces().size() << ")"
                << " contact_range=" << corner_lines[0]->contact_range;
    }
    std::cout << std::endl;
  }

  auto _t0 = std::chrono::high_resolution_clock::now();
  std::cout << __FILE__ << ":" << __LINE__ << " step4 点の境界条件を集約" << std::endl;
  {
    const std::vector<networkPoint*> pts(water->getPoints().begin(), water->getPoints().end());
    const int n_pts = static_cast<int>(pts.size());
#pragma omp parallel for
    for (int i = 0; i < n_pts; i++) {
      auto* p = pts[i];
      const auto point_bc = getVertexNodeBoundaryTypeFromFaceStates(p);
      p->Neumann = (point_bc == NodeFaceBoundaryType::Neumann);
      p->Dirichlet = (point_bc == NodeFaceBoundaryType::Dirichlet);
      bool was_corner = p->CORNER;
      p->CORNER = (!p->Neumann && !p->Dirichlet);
      if (p->CORNER != was_corner)
        p->geom_curvature.valid = false;
    }
  }
  auto _t1 = std::chrono::high_resolution_clock::now();
  std::cout << __FILE__ << ":" << __LINE__ << " step4 point loop: "
            << std::chrono::duration<double, std::milli>(_t1 - _t0).count() << " ms" << std::endl;

  //@ ------------------------------------------ */

  for (const auto& f : faces) {
    if (use_quadratic_linear_hybrid) {
      // hybridモード: Dirichlet面とCORNER隣接Neumann面のみ2次
      bool is_dirichlet_face = std::ranges::all_of(f->getPoints(), [f](const auto* p) {
        return isDirichletBoundaryState(p, f);
      });
      const auto& [l0, l1, l2] = f->Lines;
      bool has_corner_edge = l0->CORNER || l1->CORNER || l2->CORNER;
      if (is_dirichlet_face || has_corner_edge) {
        f->isTrueQuadraticElement = true;
        f->isPseudoQuadraticElement = false;
        f->isLinearElement = false;
      } else {
        f->isTrueQuadraticElement = false;
        f->isPseudoQuadraticElement = false;
        f->isLinearElement = true;
      }
    } else if (use_true_quadratic_element) {
      f->isTrueQuadraticElement = true;
      f->isPseudoQuadraticElement = false;
      f->isLinearElement = false;
    } else {
      f->isTrueQuadraticElement = false;
      f->isPseudoQuadraticElement = use_pseudo_quadratic_element &&
                                    std::ranges::all_of(f->getPoints(), [f](const auto* p) {
                                      return isDirichletBoundaryState(p, f);
                                    });
      f->isLinearElement = !f->isPseudoQuadraticElement;
    }
  }
  auto _t2 = std::chrono::high_resolution_clock::now();
  std::cout << __FILE__ << ":" << __LINE__ << " face element-type loop: "
            << std::chrono::duration<double, std::milli>(_t2 - _t1).count() << " ms" << std::endl;

  //@ ------------------------------------------ */

  std::cout << "setBoundaryTypes終了" << std::endl;
};

/* -------------------------------------------------------------------------- */
/*  refreshBoundaryStatesAndTypes                                             */
/*  initializeNodeFaceStates → setContactFaces → setBoundaryTypes を一括実行  */
/* -------------------------------------------------------------------------- */

inline void refreshBoundaryStatesAndTypes(Network* water, const std::vector<Network*>& objects) {
  initializeNodeFaceStates(water);
  water->setContactFaces(objects);
  setBoundaryTypes(water, objects);
}

inline void refreshBoundaryStatesAndTypes(std::vector<Network*>& fluid_objects, const std::vector<Network*>& contact_objects) {
  for (auto& water : fluid_objects) {
    initializeNodeFaceStates(water);
    water->setContactFaces(contact_objects);
    setBoundaryTypes(water, contact_objects);
  }
}

/* -------------------------------------------------------------------------- */
/* CORNER edge quadratic interpolation along waterline                        */
/* -------------------------------------------------------------------------- */

// Return one Neumann-side face adjacent to a CORNER line.
// The line itself is classified from node-face states; this helper only picks a
// representative face on the Neumann side for geometric interpolation.
inline networkFace* getNeumannFaceOfCornerLine(const networkLine* l) {
  for (auto* f : l->getBoundaryFaces())
    if (isNeumannBoundaryState(l, f))
      return f;
  return nullptr;
}

// Find the neighboring CORNER point along the waterline, excluding one line.
inline networkPoint* getAdjacentCornerPoint(const networkPoint* p, const networkLine* exclude_line) {
  for (auto* l : p->getLines()) {
    if (l->CORNER && l != exclude_line) {
      auto [p0, p1] = l->getPoints();
      return (p0 == p) ? p1 : p0;
    }
  }
  return nullptr; // waterline endpoint
}

// Find the neighboring CORNER line along the waterline, excluding one line.
inline networkLine* getAdjacentCornerLine(const networkPoint* p, const networkLine* exclude_line) {
  for (auto* l : p->getLines()) {
    if (l->CORNER && l != exclude_line)
      return l;
  }
  return nullptr;
}

// Check whether the Neumann-side surface is smooth across two CORNER lines.
// Returns true when the representative Neumann-face normals are nearly parallel.
inline bool isNeumannSurfaceSmooth(const networkLine* l1, const networkLine* l2) {
  auto* nf1 = getNeumannFaceOfCornerLine(l1);
  auto* nf2 = getNeumannFaceOfCornerLine(l2);
  if (!nf1 || !nf2)
    return false;
  // Threshold: cos(30 deg) ≈ 0.866 — Neumann normals must be within 30 degrees
  return Dot(nf1->normal, nf2->normal) > std::cos(30.0 * M_PI / 180.0);
}

// Lagrange quadratic interpolation through 3 points with chord-length parameterization
// P0, P1, P2: points along the curve
// Returns: interpolated position at the midpoint of edge P1-P2
inline Tddd lagrangeQuadraticMidpoint(const Tddd& P0, const Tddd& P1, const Tddd& P2) {
  double d01 = Norm(P1 - P0);
  double d12 = Norm(P2 - P1);
  if (d01 < 1e-20 || d12 < 1e-20)
    return 0.5 * (P1 + P2); // degenerate: fall back to linear

  double t0 = 0.0;
  double t1 = d01;
  double t2 = d01 + d12;
  double t_mid = t1 + 0.5 * d12; // midpoint of P1-P2 segment

  double L0 = (t_mid - t1) * (t_mid - t2) / ((t0 - t1) * (t0 - t2));
  double L1 = (t_mid - t0) * (t_mid - t2) / ((t1 - t0) * (t1 - t2));
  double L2 = (t_mid - t0) * (t_mid - t1) / ((t2 - t0) * (t2 - t1));

  return L0 * P0 + L1 * P1 + L2 * P2;
}

// Compute the quadratic midpoint offset for a single CORNER line
// Uses average of forward and backward quadratic interpolations along the waterline
// Then projects onto the Neumann surface for C0 continuity
inline void computeCornerMidpointOffset(networkLine* l) {
  l->corner_midpoint_offset = {0., 0., 0.};
  if (!l->CORNER)
    return;

  auto [p_a, p_b] = l->getPoints();
  Tddd X_a = p_a->X;
  Tddd X_b = p_b->X;
  Tddd X_linear = 0.5 * (X_a + X_b);

  // Find adjacent CORNER points and check Neumann surface smoothness
  auto* l_prev = getAdjacentCornerLine(p_a, l);
  auto* l_next = getAdjacentCornerLine(p_b, l);

  networkPoint* p_prev = nullptr;
  networkPoint* p_next = nullptr;

  if (l_prev && isNeumannSurfaceSmooth(l, l_prev))
    p_prev = getAdjacentCornerPoint(p_a, l);

  if (l_next && isNeumannSurfaceSmooth(l, l_next))
    p_next = getAdjacentCornerPoint(p_b, l);

  // Compute quadratic midpoints
  Tddd X_mid = X_linear;
  int count = 0;
  Tddd sum = {0., 0., 0.};

  if (p_prev) {
    sum += lagrangeQuadraticMidpoint(p_prev->X, X_a, X_b);
    count++;
  }
  if (p_next) {
    sum += lagrangeQuadraticMidpoint(p_next->X, X_b, X_a);
    count++;
  }

  if (count > 0)
    X_mid = sum / static_cast<double>(count);
  else
    return; // no adjacent points available, keep linear (offset = 0)

  // Project onto Neumann surface for C0 continuity
  auto* nf = getNeumannFaceOfCornerLine(l);
  if (nf) {
    auto [p0, p1, p2] = nf->getPoints();
    T3Tddd neumann_triangle = {p0->X, p1->X, p2->X};
    X_mid = Nearest(X_mid, neumann_triangle);
  }

  l->corner_midpoint_offset = X_mid - X_linear;
}

// Compute offsets for all CORNER lines in a water network
inline void computeAllCornerMidpointOffsets(Network* water) {
  for (auto* l : water->getLines()) {
    if (l->CORNER)
      computeCornerMidpointOffset(l);
    else
      l->corner_midpoint_offset = {0., 0., 0.};
  }
}

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

`phiOnFace`は，各節点`p`における各面`f`に対するポテンシャル`phi`を設定するために使用される．
`phitOnFace`は，各節点`p`における各面`f`に対するポテンシャルの時間微分`dphi/dt`を設定するために使用される．
他も同様である．

*/

/* -------------------------------------------------------------------------- */
/*                         phinOnFace, phintOnFaceの設定                       */
/* -------------------------------------------------------------------------- */

// --- networkPoint / networkLine 統一アクセサ ---
inline const Tddd& getNodeX(const networkPoint* p) { return p->X; }
inline const Tddd& getNodeX(const networkLine* l) { return l->X_mid; }

inline double getNodeTime(const networkPoint* p) { return p->RK_X.gett(); }
inline double getNodeTime(const networkLine* l) { return l->RK_X.gett(); }

inline Tddd getNormalNeumann(const networkPoint* p) { return p->getNormalNeumann_BEM(); }
inline Tddd getNormalNeumann(const networkLine* l) { return getNormalNeumann_mid(l); }

inline double& phi_ref(auto* p) { return std::get<0>(p->phiphin); }
inline double& phin_ref(auto* p) { return std::get<1>(p->phiphin); }
inline double& phit_ref(auto* p) { return std::get<0>(p->phiphin_t); }
inline double& phint_ref(auto* p) { return std::get<1>(p->phiphin_t); }

inline double getDirichletPhit(networkPoint* p) { return p->aphiat(0.); }
inline double getDirichletPhit(networkLine* l) {
  auto [pA, pB] = l->getPoints();
  return 0.5 * (std::get<0>(pA->phiphin_t) + std::get<0>(pB->phiphin_t));
}

inline double sanitizePhiTInitialGuess(double candidate, double fallback = 0.) {
  constexpr double sentinel_threshold = 1e20;
  if (!std::isfinite(candidate) || std::abs(candidate) >= sentinel_threshold)
    return fallback;
  return candidate;
}

// --- 統一テンプレート: setPhiPhinOnFace の per-face 処理 ---
template <typename Node>
inline void setPhiPhinOnFace_run(Node* node,
                                 networkFace* const f,
                                 const std::unordered_set<networkFace*>& alive_faces,
                                 const std::unordered_map<networkFace*, double>& old_phi_on_face = {},
                                 const std::unordered_map<networkFace*, double>& old_phin_on_face = {}) {
  if (f && !alive_faces.count(f))
    return;
  auto& d = node->dof(f);
  if (isNeumannBieDofKey(node, f)) {
    if (node->absorbedBy != nullptr) {
      Tddd vel = node->absorbedBy->absorb_velocity(getNodeX(node), getNodeTime(node));
      if (f == nullptr)
        d.phin = phin_ref(node) = Dot(vel, getNormalNeumann(node));
      else
        d.phin = Dot(vel, f->normal);
    } else {
      if (f == nullptr)
        d.phin = phin_ref(node) = Dot(contactNormalVelocity(node), getNormalNeumann(node));
      else
        d.phin = Dot(contactNormalVelocity(node, f), f->normal);
    }
    if (node->isMultipleNode) {
      auto it = old_phi_on_face.find(f);
      d.phi = (it != old_phi_on_face.end()) ? it->second : phi_ref(node);
    } else
      d.phi = phi_ref(node);
    d.phi_t = 1E+30;
    d.phin_t = 1E+30;
  }
  if (isDirichletBieDofKey(node, f)) {
    d.phi = phi_ref(node);
    if (node->isMultipleNode) {
      auto it = old_phin_on_face.find(f);
      d.phin = (it != old_phin_on_face.end()) ? it->second : phin_ref(node);
    } else
      d.phin = phin_ref(node);
    d.phi_t = 1E+30;
    d.phin_t = 1E+30;
  }
}

inline void setPhiPhinOnFace(Network* water) {
  // この関数は，現在の節点/辺中点の「境界種別」の結果に従って，
  // BIE に渡す各 DOF の phi / phin を毎 RK step ごとに組み直す。
  //
  // 重要なのは，
  // 1. まず前 step の dof 値を退避する
  // 2. dof を一度クリアする
  // 3. Dirichlet なら phi，Neumann なら phin を再設定する
  // という順序で処理している点。
  //
  // 多重節点では「節点そのものの代表値(nullptr 側)」と
  // 「各隣接面ごとの値(face 側)」の両方を持ちうるため，
  // setPhiPhinOnFace_run(..., nullptr, ...) と
  // setPhiPhinOnFace_run(..., face, ...) の両方を呼ぶ。
  const auto surface_vec = water->getBoundaryFaces();
  // リメッシュ等で既に消えた面の DOF を誤って再利用しないため，
  // 現在生存している境界面集合を先に作る。
  const std::unordered_set<networkFace*> alive_faces(surface_vec.begin(), surface_vec.end());

  // 頂点・辺中点の共通処理: dofs の phi/phin を再構築
  // active DOF を持たない節点（linear 辺等）は代表値を nullptr DOF に保持して終了。
  std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorReset << std::endl;
  auto setPhiPhinForNode = [&](auto* node) {
    std::unordered_map<networkFace*, double> old_phi_on_face;
    std::unordered_map<networkFace*, double> old_phin_on_face;
    for (const auto& [f, d] : node->dofs) {
      old_phi_on_face[f] = d.phi;
      old_phin_on_face[f] = d.phin;
    }
    for (auto& [f, d] : node->dofs) {
      d.phi = 0.;
      d.phin = 0.;
      d.phi_t = 0.;
      d.phin_t = 0.;
    }

    bool has_dof_index = std::ranges::any_of(node->dofs, [](const auto& kv) { return kv.second.index >= 0; });
    if (!has_dof_index) {
      auto& d0 = node->dof(nullptr);
      d0.phi = std::get<0>(node->phiphin);
      d0.phin = std::get<1>(node->phiphin);
      return;
    }

    for (const auto& [f, d] : node->dofs)
      if (d.index >= 0)
        setPhiPhinOnFace_run(node, const_cast<networkFace*>(f), alive_faces, old_phi_on_face, old_phin_on_face);
  };

  for (auto* p : water->getPoints())
    setPhiPhinForNode(p);
  for (auto* l : water->getBoundaryLines())
    setPhiPhinForNode(l);

  // 面の phi 代表値（可視化用）
  std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorReset << std::endl;
  for (const auto& f : water->getBoundaryFaces()) {
    auto [p0, p1, p2] = f->getPoints();
    std::get<0>(f->phiphin) = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3.;
  }

  dumpDebugMidpointLineState(water, "post-setPhiPhinOnFace", time_step, -1);
};

inline void setPhiPhinOnFace(const std::vector<Network*>& objects) {
  for (const auto& water : objects)
    setPhiPhinOnFace(water);
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

// \label{BEM:setPhiPhin_t}
inline void setPhiPhin_t(std::vector<Network*> WATERS) {
  /* -------------------------------------------------------------------------- */
  // --- 統一テンプレート: setPhiPhin_t の per-face 処理 ---
  auto setPhiPhin_t_run = [](auto* node,
                             networkFace* const f,
                             const std::unordered_map<networkFace*, double>& old_phit_on_face = {},
                             const std::unordered_map<networkFace*, double>& old_phint_on_face = {}) {
    auto& d = node->dof(f);
    if (isNeumannBieDofKey(node, f)) {
      if (node->absorbedBy != nullptr) {
        Tddd gradPhi_t = node->absorbedBy->absorb_gradPhi_t(getNodeX(node), getNodeTime(node));
        if (f == nullptr)
          d.phin_t = phint_ref(node) = Dot(gradPhi_t, getNormalNeumann(node));
        else
          d.phin_t = Dot(gradPhi_t, f->normal);
      } else {
        if (f == nullptr)
          d.phin_t = phint_ref(node) = phint_Neumann(node);
        else
          d.phin_t = phint_Neumann(node, f);
      }
      const double fallback_phit = sanitizePhiTInitialGuess(phit_ref(node), 0.);
      if (node->isMultipleNode) {
        auto it = old_phit_on_face.find(f);
        d.phi_t = (it != old_phit_on_face.end())
                      ? sanitizePhiTInitialGuess(it->second, fallback_phit)
                      : fallback_phit;
      } else
        d.phi_t = fallback_phit;
    }
    if (isDirichletBieDofKey(node, f)) {
      d.phi_t = getDirichletPhit(node);
      if (f == nullptr)
        phit_ref(node) = d.phi_t;
      const double fallback_phint = sanitizePhiTInitialGuess(phint_ref(node), 0.);
      if (node->isMultipleNode) {
        auto it = old_phint_on_face.find(f);
        d.phin_t = (it != old_phint_on_face.end())
                       ? sanitizePhiTInitialGuess(it->second, fallback_phint)
                       : fallback_phint;
      } else
        d.phin_t = fallback_phint;
    }
  };

  // 頂点・辺中点の共通処理: dofs の phi_t/phin_t を再構築
  // active DOF を持たない節点（linear 辺等）はスキップされる。
  std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える" << colorReset << std::endl;
  auto setPhiPhinT_ForNode = [&](auto* node) {
    std::unordered_map<networkFace*, double> old_phit_on_face;
    std::unordered_map<networkFace*, double> old_phint_on_face;
    for (const auto& [f, d] : node->dofs) {
      old_phit_on_face[f] = d.phi_t;
      old_phint_on_face[f] = d.phin_t;
    }

    for (const auto& [f, d] : node->dofs)
      if (d.index >= 0)
        setPhiPhin_t_run(node, const_cast<networkFace*>(f), old_phit_on_face, old_phint_on_face);
  };

  for (const auto water : WATERS) {
    for (auto* p : water->getPoints())
      setPhiPhinT_ForNode(p);
    for (auto* l : water->getBoundaryLines())
      setPhiPhinT_ForNode(l);
  }

  std::cout << "setPhiPhin_t終了" << std::endl;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

// mode: false = phi/phin, true = phi_t/phin_t
void storePhiPhinCommon(const std::vector<Network*>& WATERS, const V_d& ans, bool time_derivative_mode) {

  // Accessors to select phi/phin vs phi_t/phin_t in dofs and phiphin
  auto get_phi_from_dof = [&](const NodeFaceState& d) -> double { return time_derivative_mode ? d.phi_t : d.phi; };
  auto get_phin_from_dof = [&](const NodeFaceState& d) -> double { return time_derivative_mode ? d.phin_t : d.phin; };
  auto set_phi_in_dof = [&](NodeFaceState& d, double v) { if (time_derivative_mode) d.phi_t = v; else d.phi = v; };
  auto set_phin_in_dof = [&](NodeFaceState& d, double v) { if (time_derivative_mode) d.phin_t = v; else d.phin = v; };
  auto get_phiphin_0 = [&](const BEM_DOF_Base* p) -> double { return time_derivative_mode ? std::get<0>(p->phiphin_t) : std::get<0>(p->phiphin); };
  auto set_phiphin_0 = [&](BEM_DOF_Base* p, double v) { if (time_derivative_mode) std::get<0>(p->phiphin_t) = v; else std::get<0>(p->phiphin) = v; };
  auto set_phiphin_1 = [&](BEM_DOF_Base* p, double v) { if (time_derivative_mode) std::get<1>(p->phiphin_t) = v; else std::get<1>(p->phiphin) = v; };

  // 頂点（networkPoint）・辺中点（networkLine）の共通処理。
  // BEM_DOF_Base の共通インターフェース（dofs, phiphin, findActiveBieDof 等）を使用。
  //
  // 旧実装では頂点と辺で別々のループがあり、以下の不整合があった:
  //   - 頂点の phin 代表値: ループ中の最後の Dirichlet DOF で上書き（ループ順依存、不定）
  //   - 辺の phin 代表値: nullptr DOF から取得（決定的）
  // これは頂点側が Dirichlet 多重節点に対応していなかったことを意味する。
  // 統一後は両方とも nullptr DOF から代表値を取得し、決定的な値を保証する。
  auto storeForNode = [&](auto* node) {
    // 1. BVP 解を per-face DOF に書き戻す
    //    - Dirichlet DOF: phin が未知量（BVP 解）、phi は代表値から設定
    //    - Neumann DOF: phi が未知量（BVP 解）
    for (auto& [f, d] : node->dofs) {
      if (d.index < 0 || d.index >= static_cast<int>(ans.size()))
        continue;
      if (isDirichletBieDofKey(node, f)) {
        set_phin_in_dof(d, ans[d.index]);
        set_phi_in_dof(d, get_phiphin_0(node));
      }
      if (isNeumannBieDofKey(node, f)) {
        set_phi_in_dof(d, ans[d.index]);
      }
    }

    // 2. phi 代表値: 純 Neumann 節点は面積加重平均で集約
    //    CORNER / mixed node では phiphin[0] は Dirichlet 側の phi を保持すべき。
    //    Neumann 側の BVP 解で上書きすると absorber 等の効果が失われる。
    if (node->Neumann && !hasAnyDirichletBoundaryState(node)) {
      double total_area = 0., weighted_phi = 0.;
      for (auto* f : node->getBoundaryFaces()) {
        if (!f) continue;
        total_area += f->area;
        const auto* d = node->findActiveBieDofOrDefault(f);
        double phi_val = d ? get_phi_from_dof(*d) : get_phiphin_0(node);
        weighted_phi += phi_val * f->area;
      }
      if (total_area > 0.)
        set_phiphin_0(node, weighted_phi / total_area);
    }

    // 3. phin 代表値: nullptr DOF から取得
    //    - Dirichlet 節点: (node, nullptr) が Dirichlet DOF で、phin は BIE で解かれる。
    //      代表値 = BIE 解。
    //    - CORNER 節点: (node, nullptr) は Dirichlet DOF（isDirichletBieDofKey は常に f==nullptr）。
    //      per-face Neumann DOF の phin（既知値）は代表値に反映されない。
    //    - 純 Neumann multiple node: (node, nullptr) は active DOF を持たない
    //      （isNeumannBieDofKey は f!=nullptr を要求、isDirichletBieDofKey は Dirichlet 面を要求）。
    //      この場合 findActiveBieDof(nullptr) が nullptr を返し、代表値は更新されない。
    //      phin 代表値は setPhiPhinOnFace で設定された値がそのまま残る。
    if (auto* d0 = node->findActiveBieDof(nullptr))
      set_phiphin_1(node, get_phin_from_dof(*d0));
  };

  for (const auto water : WATERS) {
    for (auto* p : water->getPoints())
      storeForNode(p);
    for (auto* l : water->getBoundaryLines())
      storeForNode(l);
  }

  for (const auto water : WATERS)
    dumpDebugMidpointLineState(water, time_derivative_mode ? "post-storePhiPhin_t" : "post-storePhiPhin", time_step, -1);

  if (const char* env = std::getenv("BEM_STORE_DEBUG"); env && std::string(env) != "0") {
    double max_abs_err = 0.0;
    int max_abs_err_row = -1;
    std::size_t point_checked = 0, mid_checked = 0;

    auto checkNode = [&](const auto* node, std::size_t& count) {
      for (const auto& [f, d] : node->dofs) {
        int i = d.index;
        if (i < 0 || i >= static_cast<int>(ans.size()))
          continue;
        double stored = 0.0;
        if (isDirichletBieDofKey(node, f))
          stored = get_phin_from_dof(d);
        else if (isNeumannBieDofKey(node, f))
          stored = get_phi_from_dof(d);
        else
          continue;
        ++count;
        const double err = std::abs(stored - ans[i]);
        if (err > max_abs_err) {
          max_abs_err = err;
          max_abs_err_row = i;
        }
      }
    };

    for (const auto water : WATERS) {
      for (const auto& p : ToVector(water->getPoints()))
        checkNode(p, point_checked);
      for (const auto* l : water->getBoundaryLines())
        checkNode(l, mid_checked);
    }

    std::cout << Magenta << "[BEM:store] " << Cyan
              << "point_rows=" << Green << point_checked
              << Cyan << " mid_rows=" << Green << mid_checked
              << Cyan << " max|stored-ans|=" << Yellow << max_abs_err
              << Cyan << " (row=" << Yellow << max_abs_err_row << Cyan << ")"
              << colorReset << std::endl;
  }
}

inline void storePhiPhin(const std::vector<Network*>& WATERS, const V_d& ans) { storePhiPhinCommon(WATERS, ans, false); }

inline void storePhiPhin_t(const std::vector<Network*>& WATERS, const V_d& ans) { storePhiPhinCommon(WATERS, ans, true); }
