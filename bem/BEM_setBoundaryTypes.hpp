#pragma once

#include "Network.hpp"
#include <algorithm>
#include <unordered_set>

//! 境界条件と関連が深いのでここで定義しておく

/*DOC_EXTRACT 0_1_1_2_BOUNDARY_CONDITIONS

## 多重節点

NOTE: 面の向き$\bf n$がカクッと不連続に変わる節点には，$\phi$は同じでも，隣接面にそれぞれ対して異なる$\phi_n$を計算できるようにする

NOTE: $\bf n$が不連続に変化する節点まわりの要素は，自分のために用意された$\phi_n$を選択し補間に用いなければならない

これを多重節点という．

多重節点を導入すると，未知変数idは，節点idだけではなく，節点と面の組みのidとなる．

### 境界値問題の未知変数ID 多重節点との区別

* `isNeumannID_BEM`と`isDirichletID_BEM`は，節点と面の組みが，境界値問題の未知変数かどうかを判定する．

* `pf2ID`は，節点と面の組みを未知変数IDに変換する．多重節点でない場合は，`{p,nullptr}`が変数のキーとなり，多重節点の場合は，与えられた`{p,f}`が変数のidとなる．

*/

inline bool isNeumannID_BEM(const auto p, const netF* f) {
  if (p->Neumann || p->CORNER) {
    if (p->isMultipleNode)
      return (f != nullptr) && f->Neumann;
    else
      return (f == nullptr);
  } else
    return false;
};
inline bool isNeumannID_BEM(const std::tuple<netP*, netF*>& PF) { return isNeumannID_BEM(std::get<0>(PF), std::get<1>(PF)); };
inline bool isDirichletID_BEM(const auto p, const auto f) { return (p->Dirichlet || p->CORNER) && (f == nullptr); };
inline bool isDirichletID_BEM(const std::tuple<netP*, netF*>& PF) { return isDirichletID_BEM(std::get<0>(PF), std::get<1>(PF)); };

/*DOC_EXTRACT 0_3_BEM_utilities

## 多重節点を考慮したIDの設定方法

*/

//@ pf2IDは，setNodeFaceIndicesを実行せずとも使える．pf2Indexは，setNodeFaceIndicesを実行してから使う．
inline std::tuple<networkPoint*, networkFace*> pf2ID(const networkPoint* p, const networkFace* f) {
  if (isNeumannID_BEM(p, f) || isDirichletID_BEM(p, f))
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
  try {
    if (f == nullptr || !p->isMultipleNode || f->Dirichlet)
      return p->f2Index.at(nullptr);
    else
      return p->f2Index.at(f);
  } catch (const std::out_of_range& e) {
    for (const auto& [f, i] : p->f2Index)
      std::cout << f << " " << i << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error");
  }
};

/*

   pf2Index(p0, integ_f)を使えるように刷るためには，setNodeFaceIndicesを実行する必要がある．

*/

// --- 辺中点の多重節点ID/Index関数（頂点の pf2ID / pf2Index に対応） ---
// isNeumannID_BEM / isDirichletID_BEM は auto 化済みのため networkLine にも使用可能

inline std::tuple<networkLine*, networkFace*> lf2ID(const networkLine* l, const networkFace* f) {
  if (isNeumannID_BEM(l, f) || isDirichletID_BEM(l, f))
    return {const_cast<networkLine*>(l), const_cast<networkFace*>(f)};
  return {const_cast<networkLine*>(l), nullptr};
}

inline int lf2Index(const networkLine* l, networkFace* f) {
  if (f == nullptr || !l->isMultipleNode || (f && f->Dirichlet))
    return l->f2Index.at(nullptr);
  return l->f2Index.at(f);
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

inline std::size_t setNodeFaceIndices(const std::vector<Network*>& objects) {
  // Indexing policy:
  // - Assign indices in a spatially coherent order.
  // - Use z-descending so free-surface (high z) IDs get earlier indices.
  // - This improves locality for sparse/ILU preconditioners and makes ordering deterministic.

  std::unordered_set<networkPoint*> unique_points;
  for (const auto* water : objects)
    for (auto* p : water->getBoundaryPoints())
      unique_points.emplace(p);

  std::vector<networkPoint*> points(unique_points.begin(), unique_points.end());

  // Compute bounding box to determine the reference point (Top corner)
  // We use (min_x, min_y, max_z) as the reference point to start indexing from the "top-left-front".
  double min_x = 1e30, min_y = 1e30, max_z = -1e30;
  for (const auto* p : points) {
    min_x = std::min(min_x, std::get<0>(p->X));
    min_y = std::min(min_y, std::get<1>(p->X));
    max_z = std::max(max_z, std::get<2>(p->X));
  }
  const std::array<double, 3> ref_p = {min_x, min_y, max_z};

  // Sort by distance from the reference point.
  // This creates a "wave-front" ordering starting from the top corner,
  // which generally improves locality compared to pure coordinate slicing.
  auto dist_sq_from_ref = [&](const std::array<double, 3>& X) {
    double dx = std::get<0>(X) - ref_p[0];
    double dy = std::get<1>(X) - ref_p[1];
    double dz = std::get<2>(X) - ref_p[2]; // dz is negative or zero, but squared
    // Weight Z direction more to prioritize depth layering (optional, currently equal weights)
    return dx * dx + dy * dy + dz * dz;
  };

  auto cmp_point = [&](const networkPoint* a, const networkPoint* b) { return dist_sq_from_ref(a->X) < dist_sq_from_ref(b->X); };

  std::stable_sort(points.begin(), points.end(), cmp_point);

  std::size_t i = 0;
  for (auto* q : points) {
    q->f2Index.clear();

    if (isNeumannID_BEM(q, nullptr) || isDirichletID_BEM(q, nullptr))
      q->f2Index[nullptr] = i++;

    auto faces = q->getBoundaryFaces();
    // Sort faces of the same node based on their centroid distance from ref point as well
    // to maintain consistency with node ordering.
    auto cmp_face = [&](const networkFace* a, const networkFace* b) { return dist_sq_from_ref(a->centroid) < dist_sq_from_ref(b->centroid); };
    std::stable_sort(faces.begin(), faces.end(), cmp_face);

    for (auto* f : faces)
      if (isNeumannID_BEM(q, f) || isDirichletID_BEM(q, f))
        q->f2Index[f] = i++;
  }

  // For true quadratic elements, assign DOF indices to edge midpoints
  if (use_true_quadratic_element) {
    std::unordered_set<networkLine*> unique_lines;
    for (const auto* water : objects)
      for (auto* l : water->getBoundaryLines())
        unique_lines.emplace(l);

    std::vector<networkLine*> lines(unique_lines.begin(), unique_lines.end());

    // Sort lines by midpoint distance from reference point for spatial coherence
    auto midpoint_dist_sq = [&](const networkLine* l) {
      auto [pA, pB] = l->getPoints();
      double mx = 0.5 * (std::get<0>(pA->X) + std::get<0>(pB->X));
      double my = 0.5 * (std::get<1>(pA->X) + std::get<1>(pB->X));
      double mz = 0.5 * (std::get<2>(pA->X) + std::get<2>(pB->X));
      double dx = mx - ref_p[0], dy = my - ref_p[1], dz = mz - ref_p[2];
      return dx * dx + dy * dy + dz * dz;
    };
    std::stable_sort(lines.begin(), lines.end(), [&](const networkLine* a, const networkLine* b) {
      return midpoint_dist_sq(a) < midpoint_dist_sq(b);
    });

    for (auto* l : lines) {
      l->f2Index.clear();
      l->isMultipleNode = l->CORNER;

      // Primary DOF (Dirichlet side for CORNER, sole DOF for non-CORNER)
      if (isDirichletID_BEM(l, nullptr) || isNeumannID_BEM(l, nullptr)) {
        l->f2Index[nullptr] = i;
        l->midpoint_index = i; // backward compat
        i++;
      }

      // CORNER: additional Neumann face DOF
      if (l->CORNER) {
        auto faces = l->getBoundaryFaces();
        std::stable_sort(faces.begin(), faces.end(), [&](const networkFace* a, const networkFace* b) {
          return dist_sq_from_ref(a->centroid) < dist_sq_from_ref(b->centroid);
        });
        for (auto* f : faces) {
          if (f && f->Neumann && isNeumannID_BEM(l, f)) {
            l->f2Index[f] = i++;
            break; // CORNER edge has exactly 1 Neumann face
          }
        }
      }
    }
  }

  return i;
};

// inline std::size_t setNodeFaceIndices(const std::vector<Network *> &objects) {
//   std::size_t i = 0;
//   for (const auto water : objects)
//     for (const auto &q : water->getBoundaryPoints()) {
//       q->f2Index.clear();
//       if (isNeumannID_BEM(q, nullptr) || isDirichletID_BEM(q, nullptr))
//         q->f2Index[nullptr] = i++;
//       for (const auto &f : q->getBoundaryFaces())
//         if (isNeumannID_BEM(q, f) || isDirichletID_BEM(q, f))
//           q->f2Index[f] = i++;
//     }
//   return i;
// };

inline std::size_t setNodeFaceIndices(const Network* objects) { return setNodeFaceIndices(std::vector<Network*>{const_cast<Network*>(objects)}); };

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
    move_name_velocity = net->inputJSON["velocity"][0];
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

inline void setIsMultipleNode(const auto& p) {
  if (p->CORNER)
    p->isMultipleNode = true;
  else {
    auto n = p->getNormal_BEM();
    p->isMultipleNode = std::ranges::any_of(p->getBoundaryFaces(), [&](const auto& f) { return !isFlat(n, f->normal, 20 * M_PI / 180.); });
  }
};

/* -------------------------------------------------------------------------- */
/*                             f,l,pの境界条件を決定                             */
/* -------------------------------------------------------------------------- */

inline void setBoundaryTypes(Network* water, const std::vector<Network*>& objects) {
  std::cout << water->getName() << "の境界条件を決定 setBoundaryTypes" << std::endl;

  water->setContactFaces(objects);

  for (const auto& f : water->getBoundaryFaces())
    f->penetratedBody = nullptr;

  for (const auto& p : water->getBoundaryPoints()) {
    p->penetratedBody = nullptr;
    for (const auto& net : objects)
      if (net->InsideQ(p->X))
        p->penetratedBody = net;
  }

  std::cout << "step2 面の境界条件を判定" << std::endl;
  const auto faces = water->getBoundaryFaces();
  for (const auto& f : faces) {
    f->Neumann = std::ranges::all_of(f->getPoints(), [&f](const auto& p) { return (p->getNearestContactFace(f) != nullptr); });
    if (!f->Neumann) {
      for (const auto& net : objects)
        if (net->InsideQ(f->centroid)) {
          f->Neumann = true;
          f->penetratedBody = net;
        }
    }
    f->Dirichlet = !f->Neumann;
  }

  std::cout << "step3 線の境界条件を決定" << std::endl;
  for (const auto& l : water->getLines()) {

    l->Neumann = false;
    l->Dirichlet = false;
    l->CORNER = false;

    l->Neumann = std::ranges::all_of(l->getBoundaryFaces(), [](const auto& f) { return f->Neumann; });
    if (l->Neumann)
      continue;

    bool has_dirichlet = std::ranges::any_of(l->getBoundaryFaces(), [](const auto& f) { return f->Dirichlet; });
    bool has_neumann = std::ranges::any_of(l->getBoundaryFaces(), [](const auto& f) { return f->Neumann; });
    l->CORNER = (has_dirichlet && has_neumann);
    if (l->CORNER)
      continue;

    bool all_dirichlet = std::ranges::all_of(l->getBoundaryFaces(), [](const auto& f) { return f->Dirichlet; });
    l->Dirichlet = all_dirichlet;
    if (l->Dirichlet)
      continue;

    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "line with mixed or undefined boundary condition");
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

  // CORNER lines need contact faces from endpoint vertices to detect perpendicular body surfaces.
  // Re-run addContactFaces for CORNER lines now that CORNER flag is set.
  {
    std::vector<networkLine*> corner_lines;
    for (auto* l : water->getBoundaryLines())
      if (l->CORNER)
        corner_lines.push_back(l);
    for (const auto& l : corner_lines)
      l->clearContactFaces();
#pragma omp parallel for
    for (const auto& l : corner_lines)
      l->addContactFaces(objects, false);

    // Debug: check CORNER line contact faces
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

  std::cout << "step4 点の境界条件を決定" << std::endl;
  for (const auto& p : water->getPoints()) {
    p->Neumann = std::ranges::all_of(p->getBoundaryFaces(), [](const auto& f) { return f->Neumann; });
    p->Dirichlet = std::ranges::all_of(p->getBoundaryFaces(), [](const auto& f) { return f->Dirichlet; });
    p->CORNER = (!p->Neumann && !p->Dirichlet);
    // step5 多重節点の判定
    setIsMultipleNode(p);
  }

  //@ ------------------------------------------ */

  for (const auto& f : faces) {
    if (use_true_quadratic_element) {
      f->isTrueQuadraticElement = true;
      f->isPseudoQuadraticElement = false;
      f->isLinearElement = false;
    } else {
      f->isTrueQuadraticElement = false;
      f->isPseudoQuadraticElement = use_pseudo_quadratic_element && f->Dirichlet;
      f->isLinearElement = !f->isPseudoQuadraticElement;
    }
  }

  //@ ------------------------------------------ */

  std::cout << "setBoundaryTypes終了" << std::endl;
};

/* -------------------------------------------------------------------------- */
/* CORNER edge quadratic interpolation along waterline                        */
/* -------------------------------------------------------------------------- */

// Get the Neumann face adjacent to a CORNER line
inline networkFace* getNeumannFaceOfCornerLine(const networkLine* l) {
  for (auto* f : l->getBoundaryFaces())
    if (f->Neumann)
      return f;
  return nullptr;
}

// Find the adjacent CORNER point of p along the waterline, excluding the given line
inline networkPoint* getAdjacentCornerPoint(const networkPoint* p, const networkLine* exclude_line) {
  for (auto* l : p->getLines()) {
    if (l->CORNER && l != exclude_line) {
      auto [p0, p1] = l->getPoints();
      return (p0 == p) ? p1 : p0;
    }
  }
  return nullptr; // waterline endpoint
}

// Find the CORNER line connecting p to its adjacent CORNER point (excluding exclude_line)
inline networkLine* getAdjacentCornerLine(const networkPoint* p, const networkLine* exclude_line) {
  for (auto* l : p->getLines()) {
    if (l->CORNER && l != exclude_line)
      return l;
  }
  return nullptr;
}

// Check if Neumann surface is smooth at a CORNER point between two CORNER lines
// Returns true if the dihedral angle is close to 180 degrees (normals nearly parallel)
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
inline double& phint_ref(auto* p) { return std::get<1>(p->phiphin_t); }

inline double getDirichletPhit(networkPoint* p) { return p->aphiat(0.); }
inline double getDirichletPhit(networkLine* l) {
  auto [pA, pB] = l->getPoints();
  return 0.5 * (std::get<0>(pA->phiphin_t) + std::get<0>(pB->phiphin_t));
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
  if (isNeumannID_BEM(node, f)) {
    if (node->absorbedBy != nullptr) {
      Tddd vel = node->absorbedBy->absorb_velocity(getNodeX(node), getNodeTime(node));
      if (f == nullptr)
        node->phinOnFace[f] = phin_ref(node) = Dot(vel, getNormalNeumann(node));
      else
        node->phinOnFace[f] = Dot(vel, f->normal);
    } else {
      if (f == nullptr)
        node->phinOnFace[f] = phin_ref(node) = Dot(contactNormalVelocity(node), getNormalNeumann(node));
      else
        node->phinOnFace[f] = Dot(contactNormalVelocity(node, f), f->normal);
    }
    if (node->isMultipleNode) {
      auto it = old_phi_on_face.find(f);
      node->phiOnFace[f] = (it != old_phi_on_face.end()) ? it->second : phi_ref(node);
    } else
      node->phiOnFace[f] = phi_ref(node);
    node->phitOnFace[f] = 1E+30;
    node->phintOnFace[f] = 1E+30;
  }
  if (isDirichletID_BEM(node, f)) {
    node->phiOnFace[f] = phi_ref(node);
    if (node->isMultipleNode) {
      auto it = old_phin_on_face.find(f);
      node->phinOnFace[f] = (it != old_phin_on_face.end()) ? it->second : phin_ref(node);
    } else
      node->phinOnFace[f] = phin_ref(node);
    node->phitOnFace[f] = 1E+30;
    node->phintOnFace[f] = 1E+30;
  }
}

inline void setPhiPhinOnFace(Network* water) {
  const auto surface_vec = water->getBoundaryFaces();
  const std::unordered_set<networkFace*> alive_faces(surface_vec.begin(), surface_vec.end());

  // b! 点
  std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える" << colorReset << std::endl;

#pragma omp parallel for
  for (const auto& p : ToVector(water->getPoints())) {
    auto old_phi_on_face = std::move(p->phiOnFace);
    auto old_phin_on_face = std::move(p->phinOnFace);
    p->phiOnFace.clear();
    p->phinOnFace.clear();
    p->phitOnFace.clear();
    p->phintOnFace.clear();

    setPhiPhinOnFace_run(p, (networkFace*)nullptr, alive_faces, old_phi_on_face, old_phin_on_face);
    for (auto* face : p->getBoundaryFaces())
      setPhiPhinOnFace_run(p, face, alive_faces, old_phi_on_face, old_phin_on_face);
  }

  // b! 面
  std::cout << Green << "RKのtime step毎に，Dirichlet面にはΦを与える．Neumann面にはΦnを与える．" << colorReset << std::endl;
  for (const auto& f : water->getBoundaryFaces()) {
    auto [p0, p1, p2] = f->getPoints();
    std::get<0>(f->phiphin) = (std::get<0>(p0->phiphin) + std::get<0>(p1->phiphin) + std::get<0>(p2->phiphin)) / 3.;
  }

  // b! 辺中点（真2次要素の場合）
  if (use_true_quadratic_element) {
    for (auto* l : water->getBoundaryLines()) {
      auto old_phi_on_face = std::move(l->phiOnFace);
      auto old_phin_on_face = std::move(l->phinOnFace);
      l->phiOnFace.clear();
      l->phinOnFace.clear();
      l->phitOnFace.clear();
      l->phintOnFace.clear();

      if (l->f2Index.empty()) {
        l->phiOnFace[nullptr] = l->phiphin[0];
        l->phinOnFace[nullptr] = l->phiphin[1];
        continue;
      }

      // Bootstrap: 最初のRKステップ前は端点平均を使う
      if (l->RK_phi.steps == 0) {
        auto [pA, pB] = l->getPoints();
        phi_ref(l) = 0.5 * (std::get<0>(pA->phiphin) + std::get<0>(pB->phiphin));
      }

      for (const auto& [f, idx] : l->f2Index)
        setPhiPhinOnFace_run(l, f, alive_faces, old_phi_on_face, old_phin_on_face);
    }
  }
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
  //    for (const auto water : WATERS)
  // #pragma omp parallel
  //       for (const auto &p : water->getPoints())
  // #pragma omp single nowait
  //          for (const auto &[F, i] : p->f2Index) {
  //             // auto [p, F] = PBF;
  //             //!!ノイマンの場合はこれでDphiDtは計算できませんよ
  //             if (isDirichletID_BEM(p, F))
  //                p->phitOnFace.at(F) = std::get<0>(p->phiphin_t) = p->aphiat(0.);
  //             else if (isNeumannID_BEM(p, F)) {
  //                for (auto &[f, phin_t] : p->phintOnFace) {
  //                   // phin_t = std::get<1>(p->phiphin_t) = (f != nullptr) ? phint_Neumann(f) : phint_Neumann(p);
  //                   phin_t = (f != nullptr) ? phint_Neumann(p, f) : phint_Neumann(p);  // \label{BEM:setphint}
  //                }
  //                std::get<1>(p->phiphin_t) = phint_Neumann(p);
  //             }
  //          }

  /* -------------------------------------------------------------------------- */
  // --- 統一テンプレート: setPhiPhin_t の per-face 処理 ---
  auto setPhiPhin_t_run = [](auto* node,
                             networkFace* const f,
                             const std::unordered_map<networkFace*, double>& old_phit_on_face = {},
                             const std::unordered_map<networkFace*, double>& old_phint_on_face = {}) {
    if (isNeumannID_BEM(node, f)) {
      if (node->absorbedBy != nullptr) {
        Tddd gradPhi_t = node->absorbedBy->absorb_gradPhi_t(getNodeX(node), getNodeTime(node));
        if (f == nullptr)
          node->phintOnFace[f] = phint_ref(node) = Dot(gradPhi_t, getNormalNeumann(node));
        else
          node->phintOnFace[f] = Dot(gradPhi_t, f->normal);
      } else {
        if (f == nullptr)
          node->phintOnFace[f] = phint_ref(node) = phint_Neumann(node);
        else
          node->phintOnFace[f] = phint_Neumann(node, f);
      }
      if (node->isMultipleNode) {
        auto it = old_phit_on_face.find(f);
        node->phitOnFace[f] = (it != old_phit_on_face.end()) ? it->second : 0.;
      } else
        node->phitOnFace[f] = 0.;
    }
    if (isDirichletID_BEM(node, f)) {
      node->phitOnFace[f] = getDirichletPhit(node);
      if (node->isMultipleNode) {
        auto it = old_phint_on_face.find(f);
        node->phintOnFace[f] = (it != old_phint_on_face.end()) ? it->second : 0.;
      } else
        node->phintOnFace[f] = 0.;
    }
  };

  // b! 点
  std::cout << Green << "RKのtime step毎に，Dirichlet点にはΦtを与える．Neumann点にはΦntを与える" << colorReset << std::endl;

  for (const auto water : WATERS)
#pragma omp parallel for
    for (const auto& p : ToVector(water->getPoints())) {
      auto old_phit_on_face = std::move(p->phitOnFace);
      auto old_phint_on_face = std::move(p->phintOnFace);
      p->phitOnFace.clear();
      p->phintOnFace.clear();

      setPhiPhin_t_run(p, (networkFace*)nullptr, old_phit_on_face, old_phint_on_face);
      for (const auto& f : p->getBoundaryFaces())
        setPhiPhin_t_run(p, f, old_phit_on_face, old_phint_on_face);
    }

  // b! 辺中点（真2次要素の場合）
  if (use_true_quadratic_element) {
    for (const auto water : WATERS) {
      for (auto* l : water->getBoundaryLines()) {
        auto old_phit_on_face = std::move(l->phitOnFace);
        auto old_phint_on_face = std::move(l->phintOnFace);
        l->phitOnFace.clear();
        l->phintOnFace.clear();

        for (const auto& [f, idx] : l->f2Index)
          setPhiPhin_t_run(l, f, old_phit_on_face, old_phint_on_face);
      }
    }
  }

  std::cout << "setPhiPhin_t終了" << std::endl;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

template <typename T1, typename T2, typename T3>
void storePhiPhinCommon(const std::vector<Network*>& WATERS, const V_d& ans, T1 phiphinProperty, T2 phiOnFaceProperty, T3 phinOnFaceProperty) {

  for (const auto water : WATERS) {
#pragma omp parallel for
    for (const auto& p : ToVector(water->getPoints()))
      for (const auto& [f, i] : p->f2Index) {
        if (isDirichletID_BEM(p, f)) {
          (p->*phinOnFaceProperty).at(f) = std::get<1>(p->*phiphinProperty) = ans[i];
          (p->*phiOnFaceProperty).at(f) = std::get<0>(p->*phiphinProperty);
        }
        if (isNeumannID_BEM(p, f)) {
          (p->*phiOnFaceProperty).at(f) = std::get<0>(p->*phiphinProperty) = ans[i];
        }
      }
  }

  //^ 隣接フェイスの面積で重み付けした phi 値の寄与を計算し，結果を phiphinProperty に格納
  for (const auto water : WATERS) {
#pragma omp parallel for
    for (const auto& p : ToVector(water->getPoints()))
      if (p->Neumann) {
        double total = 0;
        std::get<0>(p->*phiphinProperty) = 0;
        for (const auto& f : p->getBoundaryFaces()) {
          total += f->area;
          if ((p->*phiOnFaceProperty).count(f))
            std::get<0>(p->*phiphinProperty) += (p->*phiOnFaceProperty).at(f) * f->area;
          else
            std::get<0>(p->*phiphinProperty) += (p->*phiOnFaceProperty).at(nullptr) * f->area;
        }
        std::get<0>(p->*phiphinProperty) /= total;
      }
  }

  // True quadratic: write back solved midpoint values to per-face maps
  // Note: For phi_t solve, phi_mid/phin_mid are temporarily swapped with phi_t_mid/phin_t_mid,
  // so writing to phiOnFace/phinOnFace here correctly updates the phi_t values
  // (because the swap also covers the per-face maps — see BEM_solveBVP.hpp).
  if (use_true_quadratic_element) {
    for (const auto water : WATERS)
      for (auto* l : water->getBoundaryLines())
        for (const auto& [f, idx] : l->f2Index) {
          if (idx < 0 || idx >= static_cast<int>(ans.size()))
            continue;
          if (isDirichletID_BEM(l, f))
            l->phinOnFace[f] = ans[idx]; // Dirichlet DOF: phin is unknown
          if (isNeumannID_BEM(l, f))
            l->phiOnFace[f] = ans[idx]; // Neumann DOF: phi is unknown
        }
    // Update representative scalars — area-weighted average of per-face values,
    // analogous to the vertex treatment above (L982-996).
    // Without this, CORNER midpoints would keep a stale nullptr-key value
    // that diverges from the BIE-solved face-specific values.
    for (const auto water : WATERS)
      for (auto* l : water->getBoundaryLines()) {
        // phi representative: area-weighted average over boundary faces
        {
          double total_area = 0., weighted_phi = 0.;
          for (auto* f : l->getBoundaryFaces()) {
            if (f == nullptr)
              continue;
            double phi_f = (l->phiOnFace.count(f))         ? l->phiOnFace.at(f)
                           : (l->phiOnFace.count(nullptr)) ? l->phiOnFace.at(nullptr)
                                                           : l->phiphin[0];
            weighted_phi += phi_f * f->area;
            total_area += f->area;
          }
          if (total_area > 0.)
            l->phiphin[0] = weighted_phi / total_area;
          else if (l->phiOnFace.count(nullptr))
            l->phiphin[0] = l->phiOnFace.at(nullptr);
        }
        // phin representative: keep from nullptr key only.
        // At CORNER edges, phin is naturally different on each face (different normals).
        // phiphin[1] feeds back into the next RK step as Neumann BC source (line 817),
        // so area-weighted averaging would corrupt the boundary condition.
        if (l->phinOnFace.count(nullptr))
          l->phiphin[1] = l->phinOnFace.at(nullptr);
      }
  }

  if (const char* env = std::getenv("BEM_STORE_DEBUG"); env && std::string(env) != "0") {
    double max_abs_err_point = 0.0;
    int max_abs_err_point_row = -1;
    double max_abs_err_mid = 0.0;
    int max_abs_err_mid_row = -1;
    std::size_t point_checked = 0;
    std::size_t mid_checked = 0;

    for (const auto water : WATERS) {
      for (const auto& p : ToVector(water->getPoints())) {
        for (const auto& [f, i] : p->f2Index) {
          if (i < 0 || i >= static_cast<int>(ans.size()))
            continue;
          double stored = 0.0;
          bool has_row = false;
          if (isDirichletID_BEM(p, f)) {
            stored = (p->*phinOnFaceProperty).at(f);
            has_row = true;
          } else if (isNeumannID_BEM(p, f)) {
            stored = (p->*phiOnFaceProperty).at(f);
            has_row = true;
          }
          if (!has_row)
            continue;
          ++point_checked;
          const double err = std::abs(stored - ans[i]);
          if (err > max_abs_err_point) {
            max_abs_err_point = err;
            max_abs_err_point_row = i;
          }
        }
      }
    }

    if (use_true_quadratic_element) {
      for (const auto water : WATERS) {
        for (const auto* l : water->getBoundaryLines()) {
          for (const auto& [f, idx] : l->f2Index) {
            if (idx < 0 || idx >= static_cast<int>(ans.size()))
              continue;
            ++mid_checked;
            double stored = 0.0;
            if (isDirichletID_BEM(l, f))
              stored = l->phinOnFace.at(f);
            else if (isNeumannID_BEM(l, f))
              stored = l->phiOnFace.at(f);
            else
              continue;
            const double err = std::abs(stored - ans[idx]);
            if (err > max_abs_err_mid) {
              max_abs_err_mid = err;
              max_abs_err_mid_row = idx;
            }
          }
        }
      }
    }

    std::cout << Magenta << "[BEM:store] " << Cyan
              << "point_rows=" << Green << point_checked
              << Cyan << " max|stored-ans|=" << Yellow << max_abs_err_point
              << Cyan << " (row=" << Yellow << max_abs_err_point_row << Cyan << ")  "
              << "mid_rows=" << Green << mid_checked
              << Cyan << " max|stored-ans|=" << Yellow << max_abs_err_mid
              << Cyan << " (row=" << Yellow << max_abs_err_mid_row << Cyan << ")"
              << colorReset << std::endl;
  }
}

inline void storePhiPhin(const std::vector<Network*>& WATERS, const V_d& ans) { storePhiPhinCommon(WATERS, ans, &networkPoint::phiphin, &networkPoint::phiOnFace, &networkPoint::phinOnFace); }

inline void storePhiPhin_t(const std::vector<Network*>& WATERS, const V_d& ans) { storePhiPhinCommon(WATERS, ans, &networkPoint::phiphin_t, &networkPoint::phitOnFace, &networkPoint::phintOnFace); }
