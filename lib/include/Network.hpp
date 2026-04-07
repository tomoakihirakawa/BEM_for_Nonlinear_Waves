#pragma once

#define BEM // BEMのメンバー変数関数を有効化する
#include "NetworkCommon.hpp"
#include "RigidBodyDynamics.hpp"
#include "basic.hpp"
#include "basic_surface_geometry.hpp"
#include "interpolations.hpp"
#include "rootFinding.hpp"
#ifdef USE_TETGEN
#include "tetgen.h"
#endif
#include "vtkWriter.hpp"
#include <cstdlib>
#include <cstdint>
#include <concepts>
#include <functional>
#include <sstream>
#include <typeinfo>
#include <unordered_set>
#include "integrationOfODE.hpp"

class networkLine;
class networkTetra;
class networkFace;
class networkPoint;
class Network;

using netP = networkPoint;
using netF = networkFace;
using netPp = networkPoint*;
using netFp = networkFace*;
using netTp = networkTetra*;
using netL = networkLine;
using netLp = networkLine*;
using V_netLp = std::vector<networkLine*>;
using V_netPp = std::vector<networkPoint*>;
using VV_netPp = std::vector<std::vector<networkPoint*>>;
using V_netFp = std::vector<networkFace*>;
using VV_netFp = std::vector<std::vector<networkFace*>>;
using V_Netp = std::vector<Network*>;
//
using T_LL = std::array<networkLine*, 2>;
using T_LLL = std::array<networkLine*, 3>;
using T_3L = std::array<networkLine*, 3>;
using T_LLLL = std::array<networkLine*, 4>;
using T_4L = std::array<networkLine*, 4>;
using T_5L = std::array<networkLine*, 5>;
using T_6L = std::array<networkLine*, 6>;
//
using T_PP = std::array<networkPoint*, 2>;
using T_2P = std::array<networkPoint*, 2>;
using T_PPP = std::array<networkPoint*, 3>;
using T_3P = std::array<networkPoint*, 3>;
using T_4P = std::array<networkPoint*, 4>;
using T_5P = std::array<networkPoint*, 5>;
using T_6P = std::array<networkPoint*, 6>;
using T4TPPP = std::array<T_PPP, 4>;
using T4T3P = std::array<T_PPP, 4>;
//
using T_FF = std::array<networkFace*, 2>;
using T_FFF = std::array<networkFace*, 3>;
using T_3F = std::array<networkFace*, 3>;
using T_4F = std::array<networkFace*, 4>;
using T_TT = std::array<networkTetra*, 2>;

enum class TetraState {
  Air = 0,
  Water = 1,
  Inside = 2,
  Outside = 3,
};

enum class FaceBadQualityType : std::uint8_t {
  TinyFace = 0,
  SubsurfaceAltitude = 1
};

struct FaceBadQualityEvent {
  int time_step = -1;
  int rk_step = -1;
  FaceBadQualityType type = FaceBadQualityType::TinyFace;
  double value = 0.0;
  double threshold = 0.0;
};

/* ------------------------------------------------------ */

// 前方宣言を追加
bool isInContact(const networkPoint* p, const T3Tddd& f_target);
std::pair<bool, Tddd> isInContact_(const networkPoint* p, const T3Tddd& f_target);
bool isInContact(const networkPoint* p, const Tddd& pX, const T3Tddd& f_target);
bool isInContact(const Tddd& X /*base*/, const Tddd& n /*base*/, const T3Tddd& f_target, const double contact_range);
bool isInContact(const networkPoint* p, const networkFace* f_normal, const networkFace* f_target);
bool isInContact(const networkPoint* p, const networkFace* f_target);
bool isInContact(const networkPoint* p, const networkFace* f_normal, const std::unordered_set<networkFace*>& faces_target);

/* ------------------------------------------------------ */

template <class T>
V_netFp extractFaces(const std::vector<T*>& Ls) {
  std::unordered_set<networkFace*> ret;
  for (const auto& l : Ls)
    for (const auto& f : l->getPoints())
      ret.emplace(f);
  return V_netFp(ret.begin(), ret.end());
};

template <class T>
V_netPp extractPoints(const std::vector<T*>& Ls) {
  std::unordered_set<networkPoint*> ret;
  for (const auto& l : Ls)
    for (const auto& p : l->getPoints())
      ret.emplace(p);
  return V_netPp(ret.begin(), ret.end());
};

template <class T>
V_netLp extractLines(const std::vector<T*>& points) {
  std::unordered_set<networkLine*> ret;
  for (const auto& p : points)
    for (const auto& l : p->getLines())
      ret.emplace(l);
  return V_netLp(ret.begin(), ret.end());
};

/* ------------------------------------------------------ */
/*              NodeFaceState                                */
/*  Per-(node, face) state: BIE values, contact, detach   */
/* ------------------------------------------------------ */

struct NodeFaceState {
  // BIE values
  double phi = 0., phin = 0.;
  double phi_FMM = 0., phin_FMM = 0.;
  double phi_t = 0., phin_t = 0.;

  // BVP system
  int index = -1;

  // Pressure / detachment
  double pressure = 0.;
  int detach_negative_count = 0;
  bool detached_by_pressure = false;

  // Contact (face pointers only; geometry recomputed on use)
  std::vector<networkFace*> contact_opponent_faces;

  networkFace* nearestContactFace() const noexcept {
    return contact_opponent_faces.empty() ? nullptr : contact_opponent_faces.front();
  }
  bool isNeumann() const noexcept {
    return !detached_by_pressure && !contact_opponent_faces.empty();
  }
};

/* ------------------------------------------------------ */
/*                  contactAcceptanceAngle                */
/* ------------------------------------------------------ */

// 近い面は広い角度で受け入れ、遠い面は厳しく絞る
// distance ≈ 0 → near_angle_deg、distance ≈ contact_range → far_angle_deg
inline double contactAcceptanceAngle(const double distance, const double contact_range,
                                     const double near_angle_deg = 60., const double far_angle_deg = 30.) {
  if (!(contact_range > 0) || !std::isfinite(contact_range) || !std::isfinite(distance))
    return M_PI * near_angle_deg / 180.;
  const auto r = std::clamp(std::abs(distance) / contact_range, 0.0, 1.0);
  auto deg = std::clamp(near_angle_deg - (near_angle_deg - far_angle_deg) * r, far_angle_deg, near_angle_deg);
  return M_PI * deg / 180.;
}

/* ------------------------------------------------------ */
/*                  ContactDetectable                     */
/* ------------------------------------------------------ */

struct ContactDetectable {
  double contact_range = 0.;
  Network* penetratedBody = nullptr;

  // --- per-(node, face) DOF (canonical storage) ---
  std::unordered_map<networkFace*, NodeFaceState> dofs;

  // Write accessor: creates entry if absent (use for assignment)
  NodeFaceState& dof(networkFace* f) { return dofs[f]; }

  // Contact/boundary state accessor: returns any DOF entry regardless of BIE index.
  // Use for reading contact_opponent_faces, nearestContactFace(), isNeumann(), detached_by_pressure.
  const NodeFaceState* findContactState(const networkFace* f) const noexcept {
    auto it = dofs.find(const_cast<networkFace*>(f));
    return (it != dofs.end()) ? &it->second : nullptr;
  }

  /* -------------------------------------------------------------------------
     Layer 3: ActiveBieDof — 実際に index >= 0 の DOF が存在するかを確認する（「結果」）

     境界条件判定の3層構造:

       Layer 1: BoundaryState — 接触状態としてどうなっているか
                                （BEM_node_face_state.hpp: isNeumannBoundaryState 等）
       Layer 2: BieDofKey     — BIE unknown の代表キーとしてどの組を使うべきか
                                （BEM_setBoundaryTypes.hpp: isNeumannBieDofKey 等）
       Layer 3: ActiveBieDof  — 実際に index >= 0 の DOF が存在するか（ここ）

     利用可能タイミング:
       setNodeFaceIndices() 完了後（DOF index が割り当てられた後）。
       それより前に呼ぶと、全ての DOF が index < 0 であるため常に nullptr を返す。

     Layer 1, 2 は「すべきか」を問う。Layer 3 は「なっているか」を問う。
  ------------------------------------------------------------------------- */

  // index >= 0 の active な BIE DOF のみを返す。
  // contact detection で作られた stale エントリ（index < 0）は除外される。
  const NodeFaceState* findActiveBieDof(const networkFace* f) const noexcept {
    auto it = dofs.find(const_cast<networkFace*>(f));
    return (it != dofs.end() && it->second.index >= 0) ? &it->second : nullptr;
  }

  // face f の active DOF を試し、なければ nullptr（primary DOF）にフォールバック。
  const NodeFaceState* findActiveBieDofOrDefault(const networkFace* f) const noexcept {
    if (auto* d = findActiveBieDof(f))
      return d;
    return findActiveBieDof(nullptr);
  }

  virtual ~ContactDetectable() = default;

  // --- pure virtual ---
  virtual const Tddd& getPosition() const = 0;
  virtual void setXSingle(const Tddd&) = 0;
  virtual std::vector<networkFace*> getBoundaryFaces() const = 0;
  virtual Network* getNetwork() const = 0;

  // Layer 3: entity が1つ以上の active BIE DOF を持つか。
  // setNodeFaceIndices() 後に使用可能。
  bool hasActiveBieDof() const noexcept {
    return std::any_of(dofs.begin(), dofs.end(), [](const auto& kv) { return kv.second.index >= 0; });
  }

  // Override to extend faces used for isFacing check (e.g., CORNER lines include endpoint vertex faces)
  virtual std::vector<networkFace*> getFacesForContactCheck() const { return getBoundaryFaces(); }

  // Unified contact face detection for both networkPoint and networkLine
  void addContactFaces(const std::vector<Network*>&, bool);

  // --- accessors (derived from dofs) ---
  std::vector<networkFace*> getContactFaces() const noexcept {
    std::vector<networkFace*> ret;
    for (const auto& [f, d] : this->dofs)
      if (!d.contact_opponent_faces.empty())
        for (auto* cf : d.contact_opponent_faces)
          if (std::ranges::find(ret, cf) == ret.end())
            ret.push_back(cf);
    return ret;
  };

  networkFace* getNearestContactFace(const networkFace* const f) const noexcept {
    const auto* d = findContactState(f);
    return d ? d->nearestContactFace() : nullptr;
  };

  void clearContactFaces() {
    for (auto& [f, d] : this->dofs)
      d.contact_opponent_faces.clear();
  };
};

/* ------------------------------------------------------ */
/*              BEM_DOF_Base                               */
/*  Shared BEM degree-of-freedom state for both           */
/*  networkPoint (vertex DOF) and networkLine (midpoint)  */
/* ------------------------------------------------------ */

struct BEM_DOF_Base : public ContactDetectable {
#if defined(BEM)
  // --- BIE unknown / solution ---
  Tdd phiphin = {0., 0.};        // {phi, phin}
  Tdd phiphin_t = {0., 0.};      // {phi_t, phin_t} (time derivative)
  Tdd phiphin_copy = {0., 0.};   // {phi, phin}
  Tdd phiphin_t_copy = {0., 0.}; // {phi_t, phin_t} (time derivative)
  double diag_coeff_BEM = 0.;    // diagonal coefficient (alpha)
  double phi_tmp = 0.;

  // --- boundary classification ---
  bool isMultipleNode = false;
  bool CORNER = false;

  // --- velocity fields ---
  Tddd u_potential_BEM = {0., 0., 0.};
  Tddd u_omega_VPM = {0., 0., 0.};
  Tddd u_total = {0., 0., 0.};
  Tddd u_reloc = {0., 0., 0.}; // velocity including smoothing correction (set by setNodeVelocity)
  Tddd X_reloc = {0., 0., 0.}; // predicted position after relocation (= RK_X target with u_reloc)

  // --- node relocation / surface ---
  Tddd vecToSurface = {0., 0., 0.};
  Tddd clungSurface = {0., 0., 0.};
  networkFace* relocation_face = nullptr;
  Tdd relocation_param = {0., 0.};

  // --- absorption ---
  Tddd u_absorbed = {0., 0., 0.};
  double phi_absorbed = 0.;
  double absorb_gamma = 0.;
  Network* absorbedBy = nullptr;
  double signed_distance = 0.;

  // --- RK4 integrators (shared by vertex and midpoint DOFs) ---
  RungeKutta<double> RK_phi;
  RungeKutta<Tddd> RK_X;

  // --- curvature (mesh density control) ---
  surface_geometry::LocalFrame local_frame;                   // 局所座標系
  surface_geometry::PrincipalCurvatureResult geom_curvature;  // 主曲率・主方向
  virtual void buildLocalFrame() = 0;

  // --- debug counters ---
  int debug_direction_info_count = 0;
  int debug_contact_faces_count = 0;
  int debug_body_vertices_count = 0;
  int debug_isInContact_pass_count = 0;

  void copy_phiphin() {
    this->phiphin_copy = this->phiphin;
    this->phiphin_t_copy = this->phiphin_t;
  }

  // --- material derivative of phi ---
  // Uses virtual getPosition() (Line→X_mid, Point→X)
  double DphiDt(const Tddd& U_modified, const double pressure, const Tddd& offset = {0., 0., 0.}, const Tddd& g_center = {0., 0., -1E+20}) const noexcept {
    const auto& pos = this->getPosition();
    double h = (offset[0] == 0. && offset[1] == 0. && offset[2] == 0. && g_center[0] == 0. && g_center[1] == 0. && g_center[2] == -1E+20)
                   ? std::get<2>(pos)
                   : Dot(pos - offset, Normalize(pos - g_center));
    return -0.5 * Dot(this->u_potential_BEM, this->u_potential_BEM) - _GRAVITY_ * h - pressure / _WATER_DENSITY_ + Dot(U_modified, this->u_potential_BEM);
  }

  double DphiDt_damped(const std::array<double, 2> gamma_phi, const Tddd& U_modified, const double pressure, const Tddd& offset = {0., 0., 0.}, const Tddd& g_center = {0., 0., -1E+20}) const noexcept {
    auto [gamma, reference_phi] = gamma_phi;
    return -gamma * (std::get<0>(this->phiphin) - reference_phi) + DphiDt(U_modified, pressure, offset, g_center);
  }
#endif
};

/* ------------------------------------------------------ */
/*     *                                                  */
/*    / \                                                 */
/*   *===*                                                */
/*    \ /                                                 */
/*     *                                                  */
/* ------------------------------------------------------ */

class networkLine : public target4FMM, public CoordinateBounds, public BEM_DOF_Base {

public:
  // std::array<double, 3> X_on_curve;

  /* ----------------------------------------------- */
  /*                 MOORING LINE                    */
  /* ----------------------------------------------- */

public:
  double diameter = 0.;               //! [m]
  double natural_length = 0.;         //! [m]
  double stiffness = 0.;              //! [N/m]
  double damping = 0.;                //! [N/(m/s^2)]
  double weight_per_unit_length = 0.; //! [kg/m]

  /* ------------------------------------------------ */
public:
#ifdef DEM
public:
  double tension;
#endif

public:
  // CORNER is inherited from BEM_DOF_Base
  bool Dirichlet;
  bool Neumann;
  Tddd corner_midpoint_offset = {0., 0., 0.}; // quadratic correction from linear midpoint for CORNER edges
  bool isIntxn();
#if defined(BEM)
  V_d interpoltedX;
  int midpoint_index = -1;   // DOF index for BIE system (-1 = not assigned)
  Tddd X_mid = {0., 0., 0.}; // absolute midpoint position

  void setXSingle(const Tddd& xyz_IN) override {
    this->X_mid = xyz_IN;
    this->Xtarget = xyz_IN;
  }
#endif

  void buildLocalFrame() override;  // 実装は Network.cpp（networkPoint の定義が必要）
public:
  //---------------------------------
  Network* network;
  Network* getNetwork() const override { return this->network; };

  // bool status;

  // bool getStatus() const noexcept { return this->status; };
  // void setStatus(const bool TorF) { this->status = TorF; };
  // bool isStatus(const bool TorF) const noexcept { return this->status ==
  // TorF; };

  std::vector<networkTetra*> Tetras;

  V_netFp Faces;

  networkPoint* Point_A;
  networkPoint* Point_B;
  // なぜこの型なのか？　いつかうunordered_setに変更する．

protected:
  bool BoundaryQ() const;
  bool ComputationalBoundaryQ() const;

  //---------------------------------
public:
  networkLine(Network* network_IN, netP* sPoint_IN, netP* ePoint_IN);
  // コピーコンストラクタ
  networkLine(const networkLine* l) : CoordinateBounds(l->getBounds()), Point_A(nullptr), Point_B(nullptr) {
    std::cout << "copying networkLine ...";
    this->network = l->getNetwork();
    this->Faces = l->getBoundaryFaces();

    this->Dirichlet = l->Dirichlet;
    this->Neumann = l->Neumann;
    this->CORNER = l->CORNER;

    auto tmpL = l->getPoints();
    set(std::get<0>(tmpL), std::get<1>(tmpL));

    setBoundsSingle();
    std::cout << " done" << std::endl;
  };
  ~networkLine();
  //---------------------------------
  // bool setBounds();
  void setBoundsSingle();
  //---------------------------------
  void set(networkPoint* const sPoint_IN, networkPoint* const ePoint_IN) {
    this->Point_A = sPoint_IN;
    this->Point_B = ePoint_IN;
  };
  void set(networkFace* const Face0_IN, networkFace* const Face1_IN) { this->Faces = {Face0_IN, Face1_IN}; };
  void set(networkFace* const Face0_IN) { this->Faces = {Face0_IN}; };
  //---------------------------------
  netPp Split(const Tddd& midX = {1E+80, 1E+80, 1E+80});
  bool canFlip(const double) const;
  bool Flip(bool);
  bool checkTopology() const;
  bool flipIfIllegal();
  bool flipIfBetter(const double min_degree_to_flat = M_PI / 180., const double min_inner_angle = M_PI / 180., const int min_n = 5);
  bool flipIfTopologicallyBetter(const double min_degree_of_line = M_PI / 180., const double min_degree_of_face = M_PI / 180, const int s_meanIN = 6);
  void divideIfIllegal();
  bool isAdjacentFacesFlat(const double) const;
  bool islegal() const;
  bool isGoodForQuadInterp() const noexcept {
    if (this->CORNER)
      return false;
    else
      return true;
  };
  bool isGoodForQuadInterp_Geo() const noexcept {
    // 線の中心位置を決めるために，線が２次補間で近似できるか，
    // 周辺の三角形の状況から判断する
    if (this->Neumann && !this->isAdjacentFacesFlat(M_PI / 3.))
      return false;
    if (this->CORNER)
      return false;
    if (!this->islegal())
      return false;
    return true;
  };
  bool isMergeable() const;
  netPp Collapse();                              // targetX を内部で決定（中点）
  netPp Collapse(const Tddd& externalTargetX);   // targetX を外部から指定
  netPp mergeIfMergeable(); // deleteしていない方のpointを返す
  //---------------------------------
  //  netP* operator()(netP* a){return this->Point2Point[a];};
  template <class T>
  T* getTheOther(const std::vector<T*>& PorF, const T* const a) const noexcept {
    if (PorF.size() < 2)
      return nullptr;
    else if (a == PorF[0])
      return PorF[1];
    else if (a == PorF[1])
      return PorF[0];
    else
      return nullptr;
  };

  netP* operator()(const netP* const a) const noexcept {
    if (this->Point_A == a)
      return this->Point_B;
    else if (this->Point_B == a)
      return this->Point_A;
    else
      return nullptr;
  };

  netF* operator()(const netF* const a) const noexcept {
    if (this->Faces.size() < 2)
      return nullptr;
    else if (a == this->Faces[0])
      return this->Faces[1];
    else if (a == this->Faces[1])
      return this->Faces[0];
    else
      return nullptr;
  };

  // VV_d getLocations() const { return obj3D::extractX(this->getPoints()); };
  T2Tddd getLocationsTuple() const;

  // V_d getNormal() const;
  Tddd getNormal() const;

  //---------------------------------
  double length() const;
  //------------------
  // V_netPp getPoints() const { return {this->Point_A, this->Point_B}; };
  T_PP getPoints(const networkPoint* const p) const noexcept {
    if (p == this->Point_A)
      return {this->Point_A, this->Point_B};
    else
      return {this->Point_B, this->Point_A};
  };
  T_PP getPoints() const noexcept { return {this->Point_A, this->Point_B}; };
  T_PP getPointsTuple() const noexcept { return {this->Point_A, this->Point_B}; };
  const V_netFp& getFaces() const noexcept { return this->Faces; };
  V_netFp getBoundaryFaces(networkFace* const f_excluded) const;
  std::vector<networkFace*> getBoundaryFaces() const override { return getBoundaryFaces(nullptr); }
  const Tddd& getPosition() const override { return this->X_mid; }
  std::vector<networkFace*> getFacesForContactCheck() const override;
  void setContactRange();

  V_netFp getFacesExcept(const netFp f) const noexcept { return TakeExcept(this->Faces, f); };
  V_netFp getFacesExcept(const V_netFp& fs) const noexcept { return TakeExcept(this->Faces, fs); };
  V_netFp getFaces(const bool TorF) const;
  //------------------
  // 派生クラスのクラス名で，選ばれる関数
  // int Find(netP *p_IN) const { return network::find(this->getPoints(), p_IN);
  // };
  bool replace(const netP* oldP, netP* newP);  // implemented in networkLine.cpp
  bool erase(const netP* p_IN) { return replace(p_IN, nullptr); };
  bool Replace(netP* oldP, netP* newP);
  //------------------
  int Find(netF* f_IN) const noexcept { return network::find(this->Faces, f_IN); };
  bool erase(netF* const f_IN) { return ::erase_element(this->Faces, f_IN); };
  bool add(netF* f_IN) { return ::add_element(this->Faces, f_IN); };
  bool replace(netF* oldF, netF* newF) {
    auto ret = network::myswitch(this->Faces, oldF, newF);
    this->Faces = DeleteDuplicates(this->Faces); // 重複を削除する
    return ret;
  };
  void Delete();
};

/* -------------------------------------------------------------------------- */

/* ------------------------------------------------------ */
/*    \ /                                                 */
/*   --@--                                                */
/*    / \                                                 */
/* ------------------------------------------------------ */

#include "integrationOfODE.hpp"

class networkPoint : public target4FMM, public CoordinateBounds, public CRS, public BEM_DOF_Base {

  /* ------------------------------ MOORING LINE ------------------------------*/
  /*
     例えば，端部には速度が与えられたとしても，内部の節点の速度は，運動方程式に従って，隣の節点との間の張力や重力，抗力そして相対速度に依存する減衰力によって決まる．
     * RK_X_sub
     * RK_velocity_sub
     の二つを使う
  */

public:
  int index = -1;
  // contact_range is now in ContactDetectable base class
  void setContactRange(const std::vector<Network*>& nets);

  std::array<double, 3> getGravitationalForce() { return this->mass * _GRAVITY3_; }

  std::array<double, 3> getTension() {
    std::array<double, 3> force, v, relative_velocity;
    force.fill(0);
    double disp, strain;
    for (const auto& l : this->getLines()) {
      v = (*l)(this)->RK_X_sub.getX() - this->RK_X_sub.getX();
      disp = Norm(v) - l->natural_length;
      strain = disp / l->natural_length;
      if (disp > 0.)
        force += l->stiffness * strain * Normalize(v);
      relative_velocity = this->RK_velocity_sub.getX() - (*l)(this)->RK_velocity_sub.getX();
      force -= l->damping * relative_velocity; // / this->RK_velocity_sub.dt_fixed;
    }
    return force;
  }

  std::array<double, 3> getDragForce(const double Cd = 0.3) {
    std::array<double, 3> drag_force, relative_velocity, mean_v, fluid_velocity, normalized_relative_velocity, A2B;
    drag_force.fill(0);
    fluid_velocity.fill(0);
    double A;
    for (const auto& l : this->getLines()) {
      // mean_v = 0.5 * ((*l)(this)->RK_velocity_sub.getX() +
      // this->RK_velocity_sub.getX());
      mean_v = this->RK_velocity_sub.getX();
      relative_velocity = fluid_velocity - mean_v;
      // A = M_PI * std::pow(l->diameter, 2);
      A2B = (*l)(this)->RK_X_sub.getX() - this->RK_X_sub.getX();
      normalized_relative_velocity = Normalize(relative_velocity);
      A = l->diameter * Norm(A2B) * 0.5;                                                           //! half area
      A *= Norm(normalized_relative_velocity - Dot(normalized_relative_velocity, Normalize(A2B))); //! projected area
      drag_force += 0.5 * _WATER_DENSITY_ * Dot(relative_velocity, relative_velocity) * Cd * A * normalized_relative_velocity;

      if (!isFinite(A2B) || !isFinite(Normalize(A2B)) || !isFinite(normalized_relative_velocity) || !isFinite(relative_velocity) || !isFinite(drag_force)) {
        std::cout << "A2B = " << A2B << std::endl;
        std::cout << "Normalize(A2B) = " << Normalize(A2B) << std::endl;
        std::cout << "normalized_relative_velocity = " << normalized_relative_velocity << std::endl;
        std::cout << "relative_velocity = " << relative_velocity << std::endl;
        std::cout << "drag_force = " << drag_force << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }
    }
    return drag_force;
  }

  std::array<double, 3> getForce() { return getTension() + getDragForce() + getGravitationalForce(); }

  /* --------------------------------------------------------------------------
   */

  std::array<std::array<double, 3>, 3> deformationGradient() {
    std::array<std::array<double, 3>, 3> F;
    F.fill({0., 0., 0.});
    for (const auto& l : this->getLines()) {
      auto [p, q] = l->getPoints();
      auto X = p->RK_X_sub.getX() - this->RK_X_sub.getX();
      auto Y = q->RK_X_sub.getX() - this->RK_X_sub.getX();
      F[0] += X;
      F[1] += Y;
      F[2] += Cross(X, Y);
    }
    return F;
  }

  /* --------------------------------------------------------------------------
   */

public:
  std::array<double, 3> X_last;
  netPp tmpPoint;
  V_netLp Lines;
  std::vector<networkFace*> Faces;
  std::vector<networkTetra*> Tetras;

  bool MemberQ(const networkFace* f_IN) const noexcept { return std::find(Faces.begin(), Faces.end(), f_IN) != Faces.end(); };

  Network* network;
  Network* getNetwork() const override { return this->network; };
  // Network *getStorage() const { return this->storage; };

  const V_netLp& getLines() const noexcept { return this->Lines; };

  int Find(netL* l_IN) { return network::find(this->Lines, l_IN); };

  bool erase(networkLine* l) { geom_curvature.valid = false; return ::erase_element(this->Lines, l); };
  bool add(networkLine* l) { geom_curvature.valid = false; return ::add_element(this->Lines, l); };
  bool erase(networkFace* f) { geom_curvature.valid = false; return ::erase_element(this->Faces, f); };
  bool add(networkFace* f) { geom_curvature.valid = false; return ::add_element(this->Faces, f); };

public:
  // 今のところBEMのためのもの
  RungeKutta<Tddd> RK_X_sub;
  RungeKutta<T3Tddd> RK_defGrad;
  T3Tddd defGrad = {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}};
  T3Tddd RightCauchyGreen = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
  T3Tddd GreenLagrangeStrain = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
  T3Tddd SecondPiolaKirchhoffStress = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
  T3Tddd FirstPiolaKirchhoffStress = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
  T3Tddd CauchyStress = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
  T3Tddd velocityGrad = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
  Tddd stressDiv = {0., 0., 0.};

  //
  // 今のところSoft bodyのためのもの
  RungeKutta<Tddd> RK_velocity, RK_velocity_sub;
  RungeKutta<T6d> RK_generalized_velocity;
  RungeKutta<Tddd> RK_force;
  //
  // 今のところSPHのためのもの
  RungeKutta<Tddd> RK_U;
  RungeKutta<double> RK_rho;
  RungeKutta<double> RK_P;
  //
  // なぜか，SPHでルンゲクッタがうまくいかないので（多分dtがDUDtに入ってくるから整合性が取れない）リープフロッグを使ってみる
  LeapFrog<Tddd> LPFG_X;
  LeapFrog<double> LPFG_rho;
  /* --------------------------------------------------------------------------
   */
public:
  V_netLp getLinesCORNER() const noexcept {
    V_netLp ret;
    ret.reserve(this->Lines.size());
    for (const auto& l : this->Lines)
      if (l->CORNER)
        ret.emplace_back(l);
    return ret;
  };
  V_netLp getLinesNeumann() const noexcept {
    V_netLp ret;
    ret.reserve(this->Lines.size());
    for (const auto& l : this->Lines)
      if (l->Neumann)
        ret.emplace_back(l);
    return ret;
  };
  V_netLp getLinesDirichlet() const noexcept {
    V_netLp ret;
    ret.reserve(this->Lines.size());
    for (const auto& l : this->Lines)
      if (l->Dirichlet)
        ret.emplace_back(l);
    return ret;
  };

  bool replace(const netL* oldL, netL* newL) {
    bool ret = false;
    for (auto& l : this->Lines)
      if (l == oldL) {
        l = newL;
        ret = true;
      }
    this->Lines = DeleteDuplicates(this->Lines); // 重複を削除する
    return ret;
  };

  bool replace(netF* oldF, netF* newF) {
    bool ret = false;
    for (auto& f : this->Faces)
      if (f == oldF) {
        f = newF;
        ret = true;
      }
    this->Faces = DeleteDuplicates(this->Faces); // 重複を削除する
    return ret;
  };

  void sortLinesByLength() {
    std::sort(this->Lines.begin(), this->Lines.end(), [](const netLp l0, const netLp l1) { return (l0->length() < l1->length()); });
  };
  /* ------------------------------------------------------ */
  /*                       For Physics                      */
  /* ------------------------------------------------------ */
  Tddd initialX; // 必ず設定される
  //
  T6d force = {0., 0., 0., 0., 0., 0.};
  T6d inertia = {0., 0., 0., 0., 0., 0.};
  T6d velocity = {0., 0., 0., 0., 0., 0.};
  T6d acceleration = {0., 0., 0., 0., 0., 0.};
  // 応力
  Tddd stress = {0., 0., 0.};
  double mass = 0.;
  Tddd forceTranslational() const noexcept { return {force[0], force[1], force[2]}; };
  Tddd forceRotational() const noexcept { return {force[3], force[4], force[5]}; };
  Tddd velocityTranslational() const noexcept { return {velocity[0], velocity[1], velocity[2]}; };
  Tddd velocityRotational() const noexcept { return {velocity[3], velocity[4], velocity[5]}; };
  Tddd accelTranslational() const noexcept { return {acceleration[0], acceleration[1], acceleration[2]}; };
  Tddd accelRotational() const noexcept { return {acceleration[3], acceleration[4], acceleration[5]}; };
  /* ------------------------------------------------------ */
  T6d& F = this->force;
  T6d& I = this->inertia;
  T6d& A = this->acceleration;
  T6d& V = this->velocity;
  /* ------------------------------------------------------ */
  double density = 0, density_ = 0;
  double& rho = this->density;
  double& rho_ = this->density_;
  double volume = 0., volume_ = 0.;
  double radius = 0.;
  double pressure = 0.;
  /* ------------------------------------------------------ */
  Tddd Fxyz() const noexcept { return {force[0], force[1], force[2]}; };
  Tddd Txyz() const noexcept { return {force[3], force[4], force[5]}; };
  Tddd Vxyz() const noexcept { return {velocity[0], velocity[1], velocity[2]}; };
  Tddd Wxyz() const noexcept { return {velocity[3], velocity[4], velocity[5]}; };
  /* ------------------------------------------------------ */

  Tddd signed_distance_vector;

  // geom_curvature is now in BEM_DOF_Base (shared by networkPoint and networkLine)

  //! SPH
  std::vector<networkPoint*> points_in_SML;

#ifdef DEM
public:
  V_netPp contactP;
  V_netPp neighborP; //! for sph
  void setParticle(double volume_IN, double densityIN) {
    this->density = densityIN;
    this->volume = volume_IN;
    this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
    this->mass = this->volume * this->density;
  };

  // 2023/05/16
  //\label{SPH:auxiliaryPoints}
  std::array<networkPoint*, 0> auxiliaryPoints;
  networkPoint* surfacePoint = nullptr;
  networkPoint* auxPoint = nullptr;
  double W;
  /////////////////////////
  // 物性
  double mu_SPH = _WATER_MU_10deg_; /*粘性係数：水は約0.001016 Pa.s*/
  //
  double dt_CFL;
  double total_weight, total_weight_, total_weight__, total_N;
  Tddd normal_SPH, v_to_surface_SPH;
  networkFace* mirroring_face;
  double d_empty_center;
  double pn_SPH;
  bool pn_is_set = false;
  bool isSurface = false, isSurface_next = false;
  bool isSurface_tmp = false, isSurface_next_tmp = false;
  bool isSurface_last_tmp = true;
  int isSurface_count = 0;
  bool isNearSurface = false;
  bool isNeumannSurface = false;
  bool isInsideOfBody = false;
  bool isCaptured = false, isCaptured_ = false;
  bool isChecked = false;
  bool isFluid = false, isFirstWallLayer = false;
  bool isAuxiliary = false;
  bool hasAuxiliary() { return this->isFluid && this->isCaptured && (this->isSurface || this->isSurface_next); };
  bool isFreeFalling = false;
  // double radius_SPH;
  double particle_spacing;
  double min_particle_spacing = 0;
  double min_particle_spacing_next = 0;

  double SML() const { return C_SML * this->particle_spacing; };
  double SML_next() const { return C_SML_next * this->particle_spacing; };
  double SML_grad() const { return 1.2 * C_SML * this->particle_spacing; };
  double SML_grad_next() const { return 1.2 * C_SML_next * this->particle_spacing; };

  networkPoint* nearest_wall_p_next = nullptr;
  networkPoint* nearest_wall_p = nullptr;
  double C_SML, C_SML_next;
  double number_density_SPH;
  double density_interpolated_SPH;
  double pressure_SPH, pressure_SPH_, pressure_SPH__;
  int pressure_equation_index = 99;
  double& p_SPH = this->pressure_SPH;
  double& p_SPH_ = this->pressure_SPH_;
  double& p_SPH__ = this->pressure_SPH__;
  double p_SPH_smoothed = 0.;
  double p_EISPH = 0, p_SPH_last = 0;
  double dp_SPH, dp_SPH_;
  double p_SPH_SPP;
  double DPDt_SPH;
  double div_U_error;
  int checked_points_in_radius_of_fluid_SPH, checked_points_in_radius_SPH, checked_points_SPH;
  int checked_points_in_radius_of_fluid_SPH_next, checked_points_in_radius_SPH_next, checked_points_SPH_next;
  double pressure_Tait(const double rho, double C0 = 1466.) const {
    // double C0 = 1466.; //[m/s]
    // C0 /= 5;
    // double C0 = Norm(this->U_SPH); //[m/s]
    // double r = 7.15;
    double r = 7.;
    double rho_w = 1000.;
    double B = rho_w * C0 * C0 / r;
    return B * (std::pow(rho / rho_w, r) - 1.);
  };
  double pressure_Tait() { return pressure_Tait(this->density); };
  double getStaticPressure() {
    double rho_w = 1000.;
    double g = 9.81;
    return -rho_w * g * std::get<2>(getXtuple());
  };
  /////////////////////////
  T3Tddd tensorproduct_grad_Uij, tensorproduct_grad_Uij_next;
  Tddd lap_U = {0., 0., 0.};
  T3Tddd grad_U = {{ {0., 0., 0.},
                     {0., 0., 0.},
                     {0., 0., 0.} }};
  T3Tddd Mat1, Mat2, Mat3, Mat_B;
  std::unordered_map<networkPoint*, double> map_p_grad;
  std::vector<std::tuple<networkPoint*, double>> vector_p_grad;
  // std::unordered_map<networkPoint *, double> extra_storage;
  // std::vector<std::tuple<networkPoint *, double>> extra_storage_vector;
  // void setLapU(const V_d &lap_U_IN) { this->lap_U = lap_U_IN; };
  void setLapU(const Tddd& lap_U_IN) { this->lap_U = lap_U_IN; };
  void setDensityVolume(const double& den, const double& v) {
    // 質量は保存
    this->density = den;
    this->volume = v;
    this->mass = v * den;
    this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
  };
  void setDensity(const double& den) {
    // 質量は保存
    this->density = den;
    this->volume = this->mass / this->density;
    this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
  };
  void setDensity_ConstantVolume(const double& den) {
    // 質量は保存
    this->density = den;
    this->mass = this->volume * this->density;
    this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
  };
  void setVolume(const double& v) {
    // 質量は保存
    this->volume = v;
    this->density = this->mass / this->volume;
    this->radius = std::pow(this->volume / (4. * M_PI / 3.), 1 / 3.);
  };
  double div_U = 0, div_tmpU = 0, div_tmpU_ = 0, PoissonRHS = 0, div_U_next = 0, DrhoDt_SPH_next = 0;
  Tddd grad_div_U = {0., 0., 0.}, grad_div_U_ = {0., 0., 0.};
  Tddd gradP_SPH = {0., 0., 0.}, gradP_SPH_ = {0., 0., 0.};
  Tddd gradP_W_SPH = {0., 0., 0.}, Integral_gradP_W_SPH = {0., 0., 0.};
  std::unordered_map<networkPoint*, Tddd> grad_coeff;
  std::unordered_map<networkPoint*, Tddd> grad_coeff_next;
  //////////////////////////
  netFp face_org;
  double a_viscosity;
  Tddd viscosity_term; // nu*laplacian(U)
  Tddd U_SPH = {0., 0., 0.}, U_SPH_ = {0., 0., 0.}, marker_X = {0., 0., 0.}, marker_U = {0., 0., 0.};
  T3Tddd nabla_otimes_U;
  Tddd U_XSPH = {0., 0., 0.}, U_XSPH_next = {0., 0., 0.}, U_next = {0., 0., 0.};
  InterpolationLagrange<std::array<double, 3>>* interp_U_lag = nullptr;
  InterpolationBspline<std::array<double, 3>>* interp_U_Bspline = nullptr;
  std::vector<double> vec_time_SPH;
  std::vector<Tddd> vec_U_SPH;
  double tmp_density;
  Tddd pre_U_SPH;
  Tddd mu_lap_rho_g_SPH;
  Tddd interp_normal, interp_normal_original, interp_normal_original_modified, interp_normal_original_choped, interp_normal_original_next_choped;
  Tddd interp_normal_water, interp_normal_water_next;
  Tddd intp_normal_Eigen;
  Tddd interp_normal_rigid, interp_normal_rigid_next;
  Tddd X_next;
  double volume_next, mass_next;
  Tddd COM_SPH, vec2COM, vec2COM_next;
  double intp_density = _WATER_DENSITY_, intp_density_next = _WATER_DENSITY_, ddr_intp_density = _WATER_DENSITY_;
  double totalMass_SPH;
  Tddd interp_normal_next, interp_normal_original_next;
  std::vector<Tddd> vector_to_polygon;
  std::vector<Tddd> vector_to_polygon_next;
  Tddd cg_neighboring_particles_SPH;
  Tddd b_vector;
  //
  std::array<std::array<double, 3>, 3> grad_corr_M = _I3_, grad_corr_M_next = _I3_;
  std::array<std::array<double, 3>, 3> laplacian_corr_M = _I3_, laplacian_corr_M_next = _I3_;
  std::array<std::array<double, 3>, 3> inv_grad_corr_M = _I3_, inv_grad_corr_M_next = _I3_;
  std::array<std::array<double, 3>, 3> Eigenvectors_of_M = _I3_, Eigenvectors_of_M1 = _I3_, Eigenvectors_of_M_next = _I3_;
  std::array<double, 3> Eigenvalues_of_M, Eigenvalues_of_M1, Eigenvalues_of_M_next, Eigenvalues_of_M1_next;
  //
  double var_Eigenvalues_of_M = 0.;
  double min_Eigenvalues_of_M = 0.;
  double min_Eigenvalues_of_M_next = 0.;
  Tddd min_Eigenvector_of_M = {0., 0., 0.};
  Tddd min_Eigenvector_of_M_next = {0., 0., 0.};
  double max_Eigenvalues_of_M = 0.;
  double var_Eigenvalues_of_M1 = 0.;
  double min_Eigenvalues_of_M1 = 0.;
  double max_Eigenvalues_of_M1 = 0.;
  double var_Eigenvalues_of_M_next = 0.;
  // ダミー粒子としての情報

  std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> v_reeDW, v_reeDW_next;
  std::array<std::array<std::array<double, 3>, 3>, 3> v_eeDW, v_eeDW_next;
  std::array<std::array<std::array<double, 3>, 3>, 3> v_rrDW, v_rrDW_next;

  /* ------------------- 多段の時間発展スキームのため ------------------- */
  Tddd DUDt_SPH = {0., 0., 0.};
  Tddd DUDt_update_SPH = {0., 0., 0.};
  Tddd DUDt_modify_SPH = {0., 0., 0.};
  Tddd DUDt_modify_SPH_2 = {0., 0., 0.};
  Tddd ViscousAndGravityForce_, tmp_ViscousAndGravityForce_;
  Tddd ViscousAndGravityForce, tmp_ViscousAndGravityForce;
  Tddd repulsive_force_SPH;
  double DrhoDt_SPH = 0.;
#endif

  std::array<double, 3> U_FVM = {0., 0., 0.};

  /* --------------------------------------------------------------------------
   */
  std::tuple<networkFace* /*補間に使った三角形の頂点*/, T_PPP, /*補間に使った三角形の頂点*/
             Tdd /*パラメタt0,t1*/, double /*深さ方向距離*/, double /*粒子間隔*/>
      particlize_info;
  /* ------------------------ 境界条件 ------------------------ */
public:
  double b_diff_RHS_FMM = 0; //! FMMのため
  double b_RHS_Direct = 0;   //! FMMのため
  double b_RHS_FMM = 0;      //! FMMのため
  // CORNER is inherited from BEM_DOF_Base
  bool Dirichlet = false;
  bool Neumann = false;

  /* --------------------------------------------------------------------------
   */

  std::map<Network*, int> net_depth;
  int minDepthFromCORNER, minDepthFromCORNER_;             // remeshのために導入
  int minDepthFromMultipleNode, minDepthFromMultipleNode_; // remeshのために導入
  //
  bool isCorner() const noexcept { return this->CORNER; };
  bool isDirichlet() const noexcept { return this->Dirichlet; };
  bool isNeumann() const noexcept { return this->Neumann; };
  void setC() {
    this->Dirichlet = false;
    this->Neumann = false;
    this->CORNER = true;
  };
  void unsetC() { this->CORNER = false; };
  //
  void setDirichlet() {
    this->Dirichlet = true;
    this->Neumann = false;
    this->CORNER = false;
  };
  void unsetDirichlet() { this->Dirichlet = false; };
  //
  void setNeumann() {
    this->Dirichlet = false;
    this->Neumann = true;
    this->CORNER = false;
  };
  void unsetNeumann() { this->Neumann = false; };
  void setD() { this->setDirichlet(); };
  void unsetD() { this->unsetDirichlet(); };
  void setN() { this->setNeumann(); };
  void unsetN() { this->unsetNeumann(); };

#ifdef BEM
  /* ----------------------------------- BEM ---------------------------------- */
public:
  std::array<double, 2> wGPhin_wGnPhi_FMM = {0., 0.}; // spherical harmonics のチェックに利用2024/07/03
  std::array<double, 2> wGPhin_wGnPhi_near = {0., 0.}, wGPhin_wGnPhi_far = {0., 0.};
  using T_PBF = std::tuple<netP*, bool, netF*>;
  std::unordered_map<T_PBF, std::unordered_map<T_PBF, Tdd>> IGIGn;
  double almost_solid_angle = 0.;
  double phin_Dirichlet = 0;

  std::unordered_map<networkFace*, double> phiOnFace_copy, phinOnFace_copy;
  std::unordered_map<networkFace*, double> phitOnFace_copy, phintOnFace_copy;

  Tddd vecToSurface_BUFFER = {0., 0., 0.}, vecToSurface_BUFFER_BUFFER = {0., 0., 0.};
  double pressure_BEM = 0.;

  double height(const Tddd& offset = {0., 0., 0.}, const Tddd& g_center = {0., 0., -1E+20}) const noexcept {
    if (offset[0] == 0. && offset[1] == 0. && offset[2] == 0. && g_center[0] == 0. && g_center[1] == 0. && g_center[2] == -1E+20)
      return std::get<2>(this->X);
    else
      return Dot(this->X - offset, Normalize(this->X - g_center));
  };

  double aphiat(const double pressure, const Tddd& offset = {0., 0., 0.}, const Tddd& g_center = {0., 0., -1E+20}) const noexcept {
    if (offset[0] == 0. && offset[1] == 0. && offset[2] == 0. && g_center[0] == 0. && g_center[1] == 0. && g_center[2] == -1E+20)
      return -0.5 * Dot(this->u_potential_BEM, this->u_potential_BEM) - _GRAVITY_ * std::get<2>(this->X) - pressure / _WATER_DENSITY_;
    else
      return -0.5 * Dot(this->u_potential_BEM, this->u_potential_BEM) - _GRAVITY_ * height(offset, g_center) - pressure / _WATER_DENSITY_;
  };

  using BEM_DOF_Base::DphiDt;
  double DphiDt(const double pressure, const Tddd& offset = {0., 0., 0.}, const Tddd& g_center = {0., 0., -1E+20}) const noexcept { return DphiDt(this->u_potential_BEM, pressure, offset, g_center); };
#endif

  // メッシュの衝突を対し噛める際に使用2022/03/10
  // bool isThereAnyFacingFace(const networkFace *const f, const double rad =
  // 1E-10) const;

private:
  // b% ------------------------------------------------------ */
  // b%       面に対する鏡像位置の粒子．一つの面に対して一点決まる　        */
  // b% ------------------------------------------------------ */
  std::unordered_map<networkFace*, networkPoint*> map_Face_MirrorPoint;

public:
  void makeMirroredPoints(const Buckets<networkFace*>& B_face, const double mirroring_distance);
  void clearMirroredPoints() {
    for (const auto& [f, p] : this->map_Face_MirrorPoint)
      delete p;
    this->map_Face_MirrorPoint.clear();
  };
  //% ------------------------------------------------------ */
  //%                      接触の判別用　                      */
  //% ------------------------------------------------------ */
  // f_nearestContactFaces, ContactFaces, contact_range, penetratedBody,
  // getContactFaces, getNearestContactFace(s), clearContactFaces
  // are now in ContactDetectable base class.

public:
  const Tddd& getPosition() const override { return this->X; }
  std::vector<networkFace*> getBoundaryFaces() const override;

  std::unordered_set<networkPoint*> ContactPoints;
  std::unordered_map<Network*, std::unordered_set<networkPoint*>> map_Net_ContactPoints;

  const std::unordered_set<networkPoint*>& getContactPoints() const noexcept { return this->ContactPoints; };
  const std::unordered_set<networkPoint*>& getContactPoints(Network* const net) const noexcept {
    auto it = this->map_Net_ContactPoints.find(net);
    if (it != this->map_Net_ContactPoints.end())
      return it->second;
    else
      return this->map_Net_ContactPoints.at(nullptr);
  };
  std::unordered_set<networkPoint*> getContactPoints(const std::vector<Network*>& net) const noexcept {
    std::unordered_set<networkPoint*> ret;
    for (const auto& n : net) {
      auto tmp = this->getContactPoints(n);
      ret.insert(tmp.begin(), tmp.end());
    }
    return ret;
  };
  void clearContactPoints() {
    this->ContactPoints.clear();
    this->map_Net_ContactPoints.clear();
    this->map_Net_ContactPoints[nullptr] = {};
  };
  networkPoint(Network* network_IN, const Tddd& xyz_IN, networkLine* xline_IN = nullptr, networkFace* xface_IN = nullptr);
  ~networkPoint();
  /* --------------------------------------------------------------------------
   */
  NewtonRaphson<double> NR_double;
  /* --------------------------------------------------------------------------
   */
  V_d getFaceAreas() const;
  std::unordered_set<networkLine*> getLinesAround() const;
  std::unordered_set<networkLine*> getLinesOppsoite() const;

  /*for a cross point*/
  //--------------
  // void set(const V_d& xyz_IN){object3D::setBounds(xyz_IN);};
  void setX(const V_d& xyz_IN);
  void setX(const Tddd& xyz_IN);
  void setXSingle(const Tddd& xyz_IN) override {
    // this->pre_X = this->X;
    this->CoordinateBounds::setBounds(xyz_IN);
    this->Xtarget = xyz_IN;
    this->geom_curvature.valid = false;
    for (auto* q : this->getNeighbors())
      if (q) q->geom_curvature.valid = false;
  };

  void buildLocalFrame() override {
    this->local_frame = {};
    Tddd normal = this->getNormalAreaAveraged();
    double n_len = Norm(normal);
    if (!(n_len > 1e-10))
      return;
    Tddd n_axis = normal / n_len;
    double best_proj_len = -1.0;
    Tddd best_u_dir = {0, 0, 0};
    for (auto* q : this->getNeighbors()) {
      if (!q) continue;
      Tddd d = q->X - this->X;
      Tddd proj = d - Dot(d, n_axis) * n_axis;
      double pl = Norm(proj);
      if (pl > best_proj_len) {
        best_proj_len = pl;
        best_u_dir = proj;
      }
    }
    if (!(best_proj_len > 1e-10))
      return;
    Tddd u_axis = best_u_dir / best_proj_len;
    Tddd v_axis = Cross(n_axis, u_axis);
    this->local_frame = {u_axis, v_axis, n_axis, true};
  }
  // bool setXcarefully(const V_d &xyz_IN);
  // void resetXinfo();
  // void setBounds();
  void setBoundsSingle() { CoordinateBounds::setBounds(this->X); };
  // V_netPp getXNeighbors() const;
  V_netPp getNeighbors() const noexcept {
    V_netPp ret(this->Lines.size());
    int i = 0;
    for (const auto& l : this->Lines)
      ret[i++] = (*l)(this);
    return ret;
  };

  int Find(const netP* lookforobj) const noexcept {
    for (int i = 0; i < this->Lines.size(); i++)
      if ((*(this->Lines[i]))(this) == lookforobj)
        return i;
    return -1;
  };

  netL* getLineBetween(const netP* const p) const noexcept {
    for (const auto& l : this->Lines)
      if (p == (*l)(this))
        return l;
    return nullptr;
  };

  bool isOpen() const noexcept {
    for (const auto& l : this->Lines)
      if (l->getBoundaryFaces().size() < 2)
        return true;
    return false;
  };
  bool isClosed() const noexcept {
    for (const auto& l : this->Lines)
      if (l->getBoundaryFaces().size() < 2)
        return false;
    return true;
  };

  V_netPp getNeighborsSort() const;
  V_netPp getNeighborsSort2() const;
  V_netPp getNeighborsSort(bool TorF);

  V_netFp getFaces(networkLine* line) const;

  V_d getAngles() const;
  V_d getAngles(networkLine* const base_line) const;

  Tddd normal_to_be_preserved = {0., 0., 0.};

#define use_binary_search
  /**
  2021/09/02unordered_setを使うよう修正した
  将来的にはunordered setを返す関数に修正すべき
   **/

  /* --------------------------------- SURFACE --------------------------------
   */

  bool BoundaryQ() const noexcept;
  bool ComputationalBoundaryQ() const noexcept;

  // getBoundaryFaces() is declared above as override of ContactDetectable
  std::vector<networkLine*> getBoundaryLines() const;
  std::vector<networkPoint*> getNeighborPointsOnSurfaces() const;

  /* --------------------------------------------------------------------------
   */

  const std::vector<networkFace*>& getFaces() const noexcept { return this->Faces; };
  std::vector<networkTetra*> getTetras() const noexcept;
  std::vector<networkFace*> setFaces(const V_netLp& this_lines);
  std::vector<networkFace*> setFaces();
  std::unordered_set<networkFace*> getFacesFromLines() const noexcept {
    std::unordered_set<networkFace*> ret;
    for (const auto& l : this->Lines)
      for (const auto& f : l->getBoundaryFaces())
        ret.emplace(f);
    return ret;
  };
  std::unordered_set<networkFace*> getFaces0() const noexcept {
    std::unordered_set<networkFace*> ret;
    for (const auto& l : this->Lines)
      for (const auto& f : l->getBoundaryFaces())
        ret.emplace(f);
    return ret;
  };
  std::unordered_set<networkFace*> getFaces1() const noexcept {
    std::unordered_set<networkFace*> ret;
    for (const auto& p : this->getNeighbors())
      ret.merge(p->getFaces0());
    return ret;
  };
  V_netFp getFacesNeumann() const;
  V_netFp getFacesDirichlet() const;

  V_netFp getFacesSort(networkLine* const l) const;
  V_netFp getFacesSort() const;

  // V_d getNormal() const override;
  Tddd getNormalTuple() const;
  Tddd getNormalDirichlet() const;
  Tddd getNormalNeumann() const;
  /* ------------------------------------------------------ */
  Tddd getNormal_BEM() const;
  Tddd getNormalDirichlet_BEM() const;
  Tddd getNormalNeumann_BEM() const;
  /* ------------------------------------------------------ */
  // Tddd getNormal_BEM_Buffer() const;
  // Tddd getNormalDirichlet_BEM_Buffer() const;
  // Tddd getNormalNeumann_BEM_Buffer() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalSplineKernelAveraged() const;
  // Tddd getNormalDirichletSplineKernelAveraged() const;
  // Tddd getNormalNeumannSplineKernelAveraged() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalSubAreaAveraged() const;
  // Tddd getNormalDirichletSubAreaAveraged() const;
  // Tddd getNormalNeumannSubAreaAveraged() const;
  /* ------------------------------------------------------ */
  Tddd getNormalAreaAveraged() const;
  Tddd getNormalDirichletAreaAveraged() const;
  Tddd getNormalNeumannAreaAveraged() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalAreaAveraged_Buffer() const;
  // Tddd getNormalDirichletAreaAveraged_Buffer() const;
  // Tddd getNormalNeumannAreaAveraged_Buffer() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalInscribedCircleAreaAveraged() const;
  // Tddd getNormalDirichletInscribedCircleAreaAveraged() const;
  // Tddd getNormalNeumannInscribedCircleAreaAveraged() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalAngleAveraged() const;
  // Tddd getNormalDirichletAngleAveraged() const;
  // Tddd getNormalNeumannAngleAveraged() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalOptimum() const;
  // Tddd getNormalDirichletOptimum() const;
  // Tddd getNormalNeumannOptimum() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalQuadInterpAngleAveraged() const;
  // Tddd getNormalNeumannQuadInterpAngleAveraged() const;
  // Tddd getNormalDirichletQuadInterpAngleAveraged() const;
  /* ------------------------------------------------------ */
  // Tddd getNormalArithmeticAveraged() const;
  // Tddd getNormalDirichletArithmeticAveraged() const;
  // Tddd getNormalNeumannArithmeticAveraged() const;

  void Delete();
  //--------------

  /*SolidAngle_detail

  SolidAngle_detail*/
  /*SolidAngle_detail_code*/
  double getSolidAngle() const;
  double getSolidAngleBuffer() const;
  double getMinimalSolidAngle() const;
  // double getSolidAngle(bool TorF);
  // double getSolidAngle(const V_netFp &faces);
  /*SolidAngle_detail_code*/
};

inline double localEdgeLength(const networkPoint* p) {
  if (!p)
    return 0.0;
  double sum = 0.0;
  std::size_t count = 0;
  for (const auto* l : p->getLines()) {
    if (!l)
      continue;
    const double len = l->length();
    if (!(len > 0.0) || !std::isfinite(len))
      continue;
    sum += len;
    ++count;
  }
  return (count > 0) ? (sum / static_cast<double>(count)) : 0.0;
}

inline double localEdgeLength(const networkLine* l) {
  if (!l)
    return 0.0;
  const auto [p0, p1] = l->getPoints();
  std::unordered_set<networkLine*> adjacent_lines;
  if (p0)
    for (const auto& L : p0->getBoundaryLines())
      if (L != l)
        adjacent_lines.insert(L);
  if (p1)
    for (const auto& L : p1->getBoundaryLines())
      if (L != l)
        adjacent_lines.insert(L);
  if (!adjacent_lines.empty()) {
    double sum = 0.0;
    std::size_t count = 0;
    for (const auto& L : adjacent_lines) {
      const double len = L->length();
      if (!(len > 0.0) || !std::isfinite(len))
        continue;
      sum += len;
      ++count;
    }
    if (count > 0)
      return sum / static_cast<double>(count);
  }
  const double len = l->length();
  if (len > 0.0 && std::isfinite(len))
    return len;
  return std::max(localEdgeLength(p0), localEdgeLength(p1));
}

inline bool isEdge(const networkPoint* p) {
  return p->getLines().empty() || std::ranges::any_of(p->getLines(), [&](const auto& l) { return (!l || l->getFaces().size() <= 1); });
};

template <>
inline double windingNumber(const Tddd& X, const std::vector<networkPoint*>& V_Points) { return 1.; };

template <std::size_t N>
std::array<Tddd, N> ToX(const std::array<networkPoint*, N>& ps) {
  std::array<Tddd, N> ret;
  for (int i = 0; i < N; i++)
    ret[i] = ps[i]->X;
  return ret;
};

inline T2Tddd ToX(const networkLine* const l) { return l->getLocationsTuple(); };
inline const Tddd& ToX(const networkPoint* const p) { return p->X; };

// Representative position for ALE/contact template unification.
// networkPoint: vertex position (X).  networkLine: quadratic midpoint position (X_mid).
inline const Tddd& getPosition(const networkPoint* const p) { return p->X; }
inline const Tddd& getPosition(const networkLine* const l) { return l->X_mid; }

template <typename T>
std::vector<Tddd> ToX(const std::unordered_set<T>& ps) {
  std::vector<Tddd> ret(ps.size());
  int i = 0;
  for (const auto& p : ps)
    ret[i++] = ToX(p);
  return ret;
};

template <typename T, std::size_t N>
std::vector<std::array<Tddd, N>> ToX(const std::unordered_set<std::array<T, N>>& ps) {
  std::vector<std::array<Tddd, N>> ret(ps.size());
  int i = 0;
  for (const auto& p : ps)
    ret[i++] = ToX(p);
  return ret;
};

template <typename T, std::size_t N>
std::vector<std::array<Tddd, N>> ToX(const std::vector<std::array<T, N>>& ps) {
  std::vector<std::array<Tddd, N>> ret(ps.size());
  int i = 0;
  for (const auto& p : ps)
    ret[i++] = ToX(p);
  return ret;
};

template <typename T>
std::vector<Tddd> ToX(const std::vector<T>& ps) {
  std::vector<Tddd> ret(ps.size());
  int i = 0;
  for (const auto& p : ps)
    ret[i++] = ToX(p);
  return ret;
};

template <std::size_t N>
std::array<double, N> ToPhi(const std::array<networkPoint*, N>& ps) {
  std::array<double, N> ret;
  for (int i = 0; i < N; i++)
    ret[i] = std::get<0>(ps[i]->phiphin);
  return ret;
};

template <std::size_t N>
std::array<double, N> ToPhin(const std::array<networkPoint*, N>& ps) {
  std::array<double, N> ret;
  for (int i = 0; i < N; i++)
    ret[i] = std::get<0>(ps[i]->phiphin);
  return ret;
};

//@ ------------------------ 抽出用関数 ----------------------- */
//@ --------------------------------------------------------- */
// std::unordered_set<networkFace *> extFaces_(const V_netPp &ps) {
//    std::unordered_set<networkFace *> ret;
//    for (const auto &p : ps)
//       for (const auto &f : p->getBoundaryFaces())
//          ret.emplace(f);
//    return ret;
// };

// 注意：unordered_setからデータを取る場合．順番は保障されない
inline std::vector<Tddd> extX(const V_netPp& ps) {
  std::vector<Tddd> ret;
  ret.reserve(ps.size());
  for (const auto& p : ps)
    ret.emplace_back(p->X);
  return ret;
};
// std::vector<Tddd> extX(const std::unordered_set<networkPoint *> &ps)
// {
// 	std::vector<Tddd> ret;
// 	for (const auto &p : ps)
// 		ret.emplace_back(p->X);
// 	return ret;
// };
// std::vector<Tddd> extXBuffer(const std::unordered_set<networkPoint *> &ps) {
//    std::vector<Tddd> ret;
//    for (const auto &p : ps)
//       ret.emplace_back(p->getXBuffer());
//    return ret;
// };
// std::vector<Tddd> extXBuffer(const std::vector<networkPoint *> &ps) {
//    std::vector<Tddd> ret;
//    for (const auto &p : ps)
//       ret.emplace_back(p->getXBuffer());
//    return ret;
// };
// std::vector<Tddd> extNormals(const V_netPp &ps) {
//    std::vector<Tddd> ret;
//    ret.reserve(ps.size());
//    for (const auto &p : ps)
//       ret.emplace_back(p->getNormalTuple());
//    return ret;
// };

// #ifdef BEM
// double extPhi(const networkPoint *const p) {
//    return std::get<0>(p->phiphin);
// };

// double extPhin(const networkPoint *const p) {
//    return std::get<1>(p->phiphin);
// };

// V_d extPhi(const V_netPp &ps) {
//    V_d ret(ps.size());
//    int i = 0;
//    for (const auto &p : ps)
//       ret[i++] = std::get<0>(p->phiphin);
//    return ret;
// };

// V_d extPhin(const V_netPp &ps) {
//    V_d ret(ps.size());
//    int i = 0;
//    for (const auto &p : ps)
//       ret[i++] = std::get<1>(p->phiphin);
//    return ret;
// };
// std::vector<Tdd> extPhiphin(const V_netPp &ps) {
//    std::vector<Tdd> ret(ps.size());
//    int i = 0;
//    for (const auto &p : ps)
//       ret[i++] = p->phiphin;
//    return ret;
// };

// #endif

//@ ------------------------------------------------------ */

/*networkPoint_code*/
// double distance(const networkPoint *p0, const networkPoint *p1) { return
// Norm(p0->X - p1->X); };

/* -------------------------------------------------------------------------- */

inline netP* Point(const netL* const l0, const netL* const l1) {
  auto [a, b] = l0->getPoints();
  auto [p, q] = l1->getPoints();
  if (a && (a == p || a == q))
    return a;
  else if (b && (b == p || b == q))
    return b;
  else
    return nullptr;
};
inline T_3P Points(const netL* const l0, const netL* const l1, const netL* const l2) { return {Point(l2, l0), Point(l0, l1), Point(l1, l2)}; };

inline netL* Line(const netP* const p0, const netP* const p1) {
  for (const auto& l : p0->Lines)
    if (p1 == (*l)(p0))
      return l;
  return nullptr;
};
inline T_3L Lines(const netP* const p0, const netP* const p1, const netP* const p2) { return {Line(p0, p1), Line(p1, p2), Line(p2, p0)}; };
/* -------------------------------------------------------------------------- */
inline netL* link(netP* const p0, netP* const p1, Network* const net) {
  // 2022/09/05
  // もし，点が既に繋がっていた場合，それをつなげる線を返す
  // もし，繋がっていない場合，新たな線を生成し，そのポインターを返す．
  // もし，矛盾した繋ぎ方をしていた場合や，つなげることができない点が与えられた場合は例外処理をする
  try {
    if (!p0 || !p1 || p1 == p0) {
      std::stringstream ss;
      ss << p0 << "<--cannot link-->" << p1;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    }

    auto l0 = Line(p0, p1);
    auto l1 = Line(p1, p0);

    if (l0 && l1) {
      if (l0 == l1)
        return l0;
      else
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Points are linked but by differrent lines");
    } else if (!l0 && l1) {
      p0->add(l1);
      return l1; // only the other side was linked
    } else if (l0 && !l1) {
      p1->add(l0);
      return l0; // only the other side was linked
    } else {
      // std::cout <<
      // "lineはコンストラクタで引数のオブジェクトを自身のオブジェクトリスト(Points,Faces)に保存し，さらにオブジェクトのLineリストに自身thisを保存する"
      // << std::endl;
      return new networkLine(net, p0, p1);
    }
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  }
};
inline std::array<netL*, 3> link(netP* const p0, netP* const p1, netP* const p2, Network* const net) { return {link(p0, p1, net), link(p1, p2, net), link(p2, p0, net)}; };
inline netL* link(netP* const p0, netP* const p1) { return link(p0, p1, p0->getNetwork()); };

/* -------------------------------------------------------------------------- */
class pathInfo {
public:
  V_netFp face;
  T3Tddd xyz; /* {start_xyz, cwp_xyz, crosspoint_xyz} */
  double incidentAngle;
  double r;
  pathInfo(V_netFp face_IN, const T3Tddd& xyz_IN, const double incidentAngle_IN, const double r_IN) : face(face_IN), xyz(xyz_IN), incidentAngle(incidentAngle_IN), r(r_IN) {};
  ~pathInfo() {
    //    std::cout << Red << "destructed" << colorReset << std::endl;
  }
  void info() {
    std::cout << Blue << "         face :" << face << std::endl;
    std::cout << Blue << "          xyz :" << xyz << std::endl;
    std::cout << Blue << "incidentAngle :" << incidentAngle << std::endl;
    std::cout << Blue << "            r :" << r << colorReset << std::endl;
  };
};

struct DodecaPoints;

/* =====================================================*/
/*     @       */
/*    / \      */
/*   @-->@     */
/* =====================================================*/

class networkFace : public Triangle {

  /*

  面の積分点が，核（カーネル）に基づくN体相互作用系のソース点(sources)と対応するので，sourcesは面のクラスが所有すべきである．

  */

public:
  std::vector<std::shared_ptr<source4FMM<target4FMM>>> sources;

  void clearSources() { sources.clear(); }

public:
  int index = -1;
  int component_id = -1; // multi-part rigid body: index into Network::component_names. -1 = single part

  /* --------------------------------------------------------------------------
   */

  std::array<double, 3> gradP1(const std::function<double(networkPoint*)>& func) {
    std::array<double, 3> V012 = {func(this->Points[0]), func(this->Points[1]), func(this->Points[2])};
    std::array<std::array<double, 3>, 3> X012 = {this->Points[0]->X, this->Points[1]->X, this->Points[2]->X};
    return ::gradP1(X012, V012);
  };

  /* --------------------------------------------------------------------------
   */

public:
  T_LLL Lines = {nullptr, nullptr, nullptr}; // 2月28日(月)導入
  T_TT Tetras = {nullptr, nullptr};

  Network* network;
  Network* getNetwork() const noexcept { return this->network; };

  void setLines(const T_LLL& ls) { this->Lines = ls; };

public:
  bool MemberQ(const networkLine* const l) const noexcept { return MemberQ_(this->Lines, l) && MemberQ_(this->Lines, l) && MemberQ_(this->Lines, l); };
  bool MemberQ(const networkPoint* const p) const noexcept { return MemberQ_(this->Points, p); };
  bool AllMemberQ(const T_3P& points) const noexcept { return MemberQ_(this->Points, std::get<0>(points)) && MemberQ_(this->Points, std::get<1>(points)) && MemberQ_(this->Points, std::get<2>(points)); };
  bool AnyMemberQ(const T_3P& points) const noexcept { return MemberQ_(this->Points, std::get<0>(points)) || MemberQ_(this->Points, std::get<1>(points)) || MemberQ_(this->Points, std::get<2>(points)); };
  // Tdd grid_pull_factor;
  int grid_pull_depth;

  Tddd normalVelocityRigidBody(const Tddd& X) const;
  //% ==================================================================
  //%                        接続情報の事前計算
  //% ==================================================================

  std::array<std::shared_ptr<DodecaPoints>, 3> dodecaPoints;
  void setDodecaPoints();

  // # ==================================================================
  // #                               BEM
  // # ==================================================================

  bool is_active_face = false;
  int minDepthToCORNER;
  Network* penetratedBody = nullptr;

  // Pressure-based detachment
  double pressure_neumann = 0.;      // 面平均Neumann側圧力（前ステップ）
  int detach_negative_count = 0;     // 負圧連続ステップ数
  bool detached_by_pressure = false; // 圧力剥離されたフラグ

  // # ------------------------------------------------------------------
  // # BEMの積分を効率的に行うためにあらかじめ積分情報を作成しておく
  // # ------------------------------------------------------------------

  using linear_triangle_integration_info = std::tuple<Tdd /*0: 2D parameter {[0,1], [0,1]} (integration variables)*/, double /*1: gaussian weight (integration weight)*/, Tddd /*2: 3D parameter {xi0=[0,1], xi1=[0,1-xi0]} to move on a triangle and assscociated with the 2D parameter*/, Tddd /*3: 3d position vector using {xi0,xi1,xi2}*/, Tddd /*4: cross product dX/dxi0 x dX/dxi1 at the point*/, double /*5: the norm of the cross product*/>;
  using pseudo_quadratic_triangle_integration_info = std::tuple<Tdd /*0: 2D parameter {[0,1], [0,1]} (integration variables)*/, double /*1: gaussian weight (integration weight)*/, std::array<T6d, 4> /*3: shape function for the quadratic element*/, Tddd /*4: 3d position vector using {xi0,xi1,xi2}*/, Tddd /*5: cross product dX/dxi0 x dX/dxi1 at the point*/, double /*6: the norm of the cross product*/>;

  std::vector<std::unordered_map<networkPoint*, std::vector<linear_triangle_integration_info>>> map_Point_LinearIntegrationInfo_vector;
  std::vector<std::unordered_map<networkPoint*, std::vector<pseudo_quadratic_triangle_integration_info>>> map_Point_PseudoQuadraticIntegrationInfo_vector;

  using BEM_IGIGn_info_type = std::tuple<networkPoint*, networkFace*, double, double>;
  std::unordered_map<networkPoint*, std::vector<BEM_IGIGn_info_type>> map_Point_BEM_IGIGn_info_init;

  void setIntegrationInfo();

  bool isPseudoQuadraticElement = false;
  bool isTrueQuadraticElement = false;
  bool isLinearElement = false;
  std::vector<FaceBadQualityEvent> bad_quality_history;

  // # ------------------------------------------------------------------
  // # True quadratic 要素用の事前計算データ
  // #
  // # true_quadratic 要素では、3頂点 + 3辺中点の6自由度を持つ。
  // # BIE の直接積分 (setDirectIntegration) で繰り返し必要になる
  // # 幾何量（曲面上の位置・クロス積）と形状関数値を、
  // # 各求積ルールの求積点で事前に計算して保存する。
  // #
  // # 幾何は2次補間（6節点: 3頂点 + 3辺中点座標）で評価する。
  // # 求積点の数は vector で管理し、ルール変更に柔軟に対応する。
  // #
  // # ハイブリッド重み分配ルール（DodecaPoints::N6_new の else 分岐と同構造）:
  // #   - Dirichlet面: 全辺中点が独立DOF（full quadratic）
  // #   - Neumann面:
  // #     - CORNER辺の中点: 独立DOF（重み保持）
  // #     - 非CORNER辺の中点: 両端頂点に 0.5 ずつ再分配（実質linear）
  // #   CORNER辺を0〜2本持つNeumann面が存在しうる。
  // # ------------------------------------------------------------------

  /// 求積点での事前計算データ（幾何・形状関数・重み）
  struct TrueQuadQP {
    Tddd X;                   ///< 曲面上の位置（2次幾何補間）
    Tddd cross;               ///< クロス積 dX/dt0 × dX/dt1（法線方向・局所Jacobian）
    std::array<double, 6> N6; ///< 形状関数値（ハイブリッド再分配適用済み）
    double w;                 ///< 求積重み
  };

  /// 6DOF のグローバルインデックス [0-2]=頂点, [3-5]=辺中点
  /// 非活性な辺中点（非CORNER Neumann辺）は負値: 再分配済みのため直接参照しない
  std::array<int32_t, 6> true_quad_dof_indices = {-1, -1, -1, -1, -1, -1};

  /// 遠方非隣接ソース用の求積点（例: Dunavant p2m設定点数）
  std::vector<TrueQuadQP> true_quad_qp_far;
  /// 近傍非隣接ソース用の求積点（例: GW5x5 25点）
  std::vector<TrueQuadQP> true_quad_qp_near;
  /// 辺中点Duffy用の求積点（辺番号 0/1/2 ごとに 2×25=50点）
  /// 辺分割はターゲット非依存（辺番号のみ依存）なので事前計算可能
  /// 注: 頂点Duffyは頂点並べ替えがターゲット依存のため事前計算不可
  std::array<std::vector<TrueQuadQP>, 3> true_quad_qp_duffy;

  /// 辺中点DOFの活性状態を返す
  /// midpoint_index >= 0 の辺のみ中点DOFとして有効。
  /// PLPLPL辺順序: l0=p0-p1, l1=p1-p2, l2=p2-p0
  /// 戻り値: {l0活性, l1活性, l2活性}
  std::array<bool, 3> trueQuadMidpointActive() const noexcept {
    auto [p0, l0, p1, l1, p2, l2] = this->PLPLPL;
    (void)p0;
    (void)p1;
    (void)p2;
    return {l0->midpoint_index >= 0, l1->midpoint_index >= 0, l2->midpoint_index >= 0};
  }

  /// ハイブリッド true_quadratic 形状関数を計算
  /// 非活性な辺中点の重みを両端頂点に 0.5 ずつ再分配
  /// TriShape<6> 順序: N0=p0, N1=p1, N2=p2, N3=mid(l0), N4=mid(l1), N5=mid(l2)
  std::array<double, 6> trueQuadN6(double t0, double t1) const noexcept {
    const auto mid_active = trueQuadMidpointActive();
    auto [N0, N1, N2, N3, N4, N5] = TriShape<6>(t0, t1);
    std::array<double, 6> N6 = {N0, N1, N2, N3, N4, N5};
    // l0 (p0-p1間): 非活性ならN3を頂点p0,p1に再分配
    if (!mid_active[0]) {
      N6[0] += 0.5 * N6[3];
      N6[1] += 0.5 * N6[3];
      N6[3] = 0.;
    }
    // l1 (p1-p2間): 非活性ならN4を頂点p1,p2に再分配
    if (!mid_active[1]) {
      N6[1] += 0.5 * N6[4];
      N6[2] += 0.5 * N6[4];
      N6[4] = 0.;
    }
    // l2 (p2-p0間): 非活性ならN5を頂点p2,p0に再分配
    if (!mid_active[2]) {
      N6[2] += 0.5 * N6[5];
      N6[0] += 0.5 * N6[5];
      N6[5] = 0.;
    }
    return N6;
  }

  Tdd phiphin = {0., 0.};
  Tddd phinTuple = {0., 0., 0.};

  void clearBadQualityHistory() { bad_quality_history.clear(); }

  void appendBadQualityEvent(const FaceBadQualityEvent& event, const std::size_t max_history = 16) {
    if (!bad_quality_history.empty()) {
      auto& last = bad_quality_history.back();
      if (last.time_step == event.time_step &&
          last.rk_step == event.rk_step &&
          last.type == event.type) {
        last.value = std::min(last.value, event.value);
        last.threshold = event.threshold;
        return;
      }
    }
    bad_quality_history.emplace_back(event);
    while (bad_quality_history.size() > max_history)
      bad_quality_history.erase(bad_quality_history.begin());
  }

  // # ========================================================================

public:
  T_PPP Points;
  std::tuple<networkPoint*, networkLine*, networkPoint*, networkLine*, networkPoint*, networkLine*> PLPLPL;
  //% ------------------------------------------------------ */
  std::unordered_set<networkPoint*> ParametricPoints;
  //! ------------------------------------------------------ */
  std::unordered_set<networkPoint*> ContactPoints;
  std::unordered_map<Network*, std::unordered_set<networkPoint*>> map_Net_ContactPoints = {{nullptr, {}}};

public:
  //% ------------------------------------------------------ */
  //%           パーティクライズで作られるParametricPoints        */
  //% ------------------------------------------------------ */
  const std::unordered_set<networkPoint*>& getParametricPoints() const noexcept { return this->ParametricPoints; };
  void clearParametricPoints() {
    for (const auto& p : this->ParametricPoints)
      delete p;
    this->ParametricPoints.clear();
  };
  void addParametricPoints(networkPoint* const p) { this->ParametricPoints.emplace(p); };
  void addParametricPoints(const V_netPp& P) { this->ParametricPoints.insert(P.begin(), P.end()); };
  void addParametricPoints(const std::unordered_set<networkPoint*>& P) { this->ParametricPoints.insert(P.begin(), P.end()); };
  //! ------------------------------------------------------ */
  //!                          接触の判別                      */
  //! ------------------------------------------------------ */
  const std::unordered_set<networkPoint*>& getContactPoints() const noexcept { return this->ContactPoints; };
  const std::unordered_set<networkPoint*>& getContactPoints(Network* const net) const noexcept {
    auto it = this->map_Net_ContactPoints.find(net);
    if (it != this->map_Net_ContactPoints.end())
      return it->second;
    else
      return this->map_Net_ContactPoints.at(nullptr);
  };
  const std::unordered_set<networkPoint*> getContactPoints(const std::vector<Network*>& nets) const noexcept {
    std::unordered_set<networkPoint*> ret;
    for (const auto& n : nets) {
      auto points = getContactPoints(n);
      ret.insert(points.begin(), points.end());
    }
    return ret;
  };
  void clearContactPoints() {
    this->ContactPoints.clear();
    this->map_Net_ContactPoints.clear();
    this->map_Net_ContactPoints[nullptr] = {};
  };
  //! ------------------------------------------------------ */
  // void reverseNormal() {
  //    // std::reverse(this->Lines.begin(), this->Lines.end());
  //    this->Lines = Reverse(this->Lines);
  //    this->setGeometricProperties(ToX(this->setPoints()));  //
  //    setBoundsは，setPointsFromCurrentLines()を実行する．
  // };
  /* ------------------------------------------------------ */
  Tddd mirror(const Tddd& v) const noexcept {
    // 与えられたvの鏡像を返す
    return v - 2. * Dot(v, this->normal) * this->normal;
  };
  Tddd mirrorPosition(const Tddd& v) const noexcept {
    // 与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
    auto u = std::get<0>(this->Points)->X;
    return v + 2. * Dot(u - v, this->normal) * this->normal;
  };
  Tddd mirrorPosition(const Tddd& v, const double scale) const noexcept {
    // 与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
    return v + scale * Dot(this->X - v, this->normal) * this->normal;
  };
  Tddd mirrorPosition(const networkPoint* p, const double scale = 2.) const noexcept {
    // 与えられた位置ベクトルに対して，鏡像位置を返す．ただし，最短距離にできる鏡像位置．
    return p->X + scale * Dot(this->X - p->X, this->normal) * this->normal;
  };

  // networkFace(Network *network_IN, const T_LLL &Lines_IN, T_3P);
  networkFace(Network* network_IN, networkPoint* const p0, networkLine* const l0, networkPoint* const p1, networkLine* const l1, networkPoint* const p2, networkLine* const l2);
  networkFace(Network* network_IN, networkPoint* const p0, networkPoint* const p1, networkPoint* const p2);
  /*コピーコンストラクタ*/
  networkFace(const netFp f);
  ~networkFace();

  std::unordered_set<networkPoint*> particlize(const double dx, const V_d& depth_list);

  V_netFp getNeighbors() const noexcept {
    V_netFp ret;
    networkFace* f;
    std::ranges::for_each(this->Lines, [&](const auto& l) {
      if (f = (*l)(this)) {
        ret.emplace_back(f);
      }; });
    return ret;
  };

  netL* getLineBetween(const networkFace* const f) const noexcept {
    auto [l0, l1, l2] = this->Lines;
    if (l0 != nullptr && f == (*l0)(this))
      return l0;
    else if (l1 != nullptr && f == (*l1)(this))
      return l1;
    else if (l2 != nullptr && f == (*l2)(this))
      return l2;
    else
      return nullptr;
  };

public:
  T_PPP getPointsFromLines(const T_LLL& lines) const {
    //! p0 <-> line[0] <-> p1 <-> line[1] <-> p2 <-> line[2] を 作成
    const auto [l0, l1, l2] = lines;
    const auto l0_points = l0->getPoints();
    const auto l1_points = l1->getPoints();
    const auto l2_points = l2->getPoints();
    auto A = Intersection(l0_points, l2_points);
    auto B = Intersection(l1_points, l0_points);
    auto C = Intersection(l2_points, l1_points);
    if ((A.size() != 1 || B.size() != 1 || C.size() != 1)) {
      std::stringstream ss;
      ss << "\nthis:" << this << "\nlines:" << lines;
      ss << "\nl0_points:" << l0_points << "\nl1_points:" << l1_points << "\nl2_points:" << l2_points << "\n";
      ss << "\ngetPoints():" << this->getPoints() << "\n";
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Intersection size is not 1. " + ss.str());
    } else
      return {*A.begin(), *B.begin(), *C.begin()};
  };

  T_PPP getPointsFromLines() const noexcept { return getPointsFromLines(this->Lines); };

  void setPoints(networkPoint* const p0, networkLine* const l0, networkPoint* const p1, networkLine* const l1, networkPoint* const p2, networkLine* const l2) {
    this->PLPLPL = {p0, l0, p1, l1, p2, l2};
    this->Points = {p0, p1, p2};
    this->Lines = {l0, l1, l2};
  };

  void reverse() {
    auto [p0, l0, p1, l1, p2, l2] = this->PLPLPL;
    this->setPoints(p1, l0, p0, l2, p2, l1);
  };

  T_PPP setPoints(netL* l0, netL* l1, netL* l2) {
    try {
      const auto [p0, p1, p2] = this->Points = getPointsFromLines({l0, l1, l2});
      this->PLPLPL = {p0, l0, p1, l1, p2, l2};
      return this->Points;
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    }
  };

  T_PPP setPoints(const T_LLL& lines) {
    const auto [l0, l1, l2] = lines;
    return setPoints(l0, l1, l2);
  }

  T_PPP syncPLPLPL() {
    try {
      const auto [p0, p1, p2] = this->Points = getPointsFromLines(this->Lines);
      this->PLPLPL = {p0, std::get<0>(this->Lines), p1, std::get<1>(this->Lines), p2, std::get<2>(this->Lines)};
      return this->Points;
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    }
  };

  T_PPP setPoints() { return setPoints(this->Lines); };

  bool setGeometricProperties(const T3Tddd& toX_this_points) {
    /**
     *! 依存関係の明示方法
     * 下のように使う．これによって，setPointsFromLinesした後に，これを実行する必要があることを印象付けることができる．
     * f->setGeometricProperties(f->setPointsFromLines())
     **/
    //@ Pointsに従って，値を設定する．
    //@
    // 線の持つ点のインターセクションをチェックすることで，面の持つ点を間接的に取得する．
    try {
      /**
      @ networkFacesの持つ
      @ this->Points = {p0,p1,p2}
      @ this->Lines = {l0,l1,l2}
      @ の関係:
      @ 		    p2
      @         /\
      @ 		   /a2\
      @       /    \
      @  l2  /      \ l1
      @ 	   /a0    a1\
      @    -------------
      @  p0      l0      p1
      */
      /*
      b# setBoundsの伝播．おそらくトポロジカルな変更はなく，値の再計算のみでいい
      @ 点の座標変更：PointのsetBounds()を実行
      @   V
      @  PointのsetBounds()   !!!!!!
      Facesをどうにかして決めれないか？-->全てのsetBounds()が終わった後に，取り込めばいい．
      @   |
      $   |      線のつなぎ変え：LineのsetBounds()を実行
      $   |       V
      $   +-----> LineのsetBounds()
      !           |
      !           +--->FaceのsetBounds() -(LinesとPointsは確定済み)->
      normal,area,angles
      *
      @  PointのsetBoundsSingle() ---->
      CoordinateBounds::setBounds(this->X);　伝播なし $  LineのsetBoundsSingle()
      ----> CoordinateBounds::setBounds(this->X);　伝播なし
      */
      CoordinateBounds::setBounds(toX_this_points);
      Triangle::setProperties(toX_this_points);
      if (!isFinite(toX_this_points) || !isFinite(this->area) || !isFinite(this->normal) || !isFinite(this->angles)) {
        std::stringstream ss;
        ss << "線が更新されておらずエラーになる可能性がある" << std::endl;
        ss << "全ての線が更新されている必要があるため" << std::endl;
        ss << "this->Points = " << this->Points << std::endl;
        ss << "this->Lines = " << this->Lines << std::endl;
        ss << "toX_this_points = " << toX_this_points << std::endl;
        ss << "this->area = " << this->area << std::endl;
        ss << "this->normal = " << this->normal << std::endl;
        ss << "this->angles = " << this->angles << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      } else
        return true;
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
  };

  bool setGeometricProperties() { return setGeometricProperties(ToX(this->Points)); };

  //------------------------
  void Delete();
  //------------------------
  // bool Replace(netL *oldL, netL *newL, netF *newF = nullptr) { return
  // network::replace(this, this->Lines, oldL, newL, newF); };
  bool Replace(netL* oldL, netL* newL) {
    // switchでないと，順番に意味のあるFaceではおかしくなるので注意
    if (this->replace(oldL, newL) && oldL->erase(this) && newL->add(this))
      return true;
    else
      return false;
  };

  double getSubArea(const networkPoint* const p) const noexcept {
    auto [p0, p1, p2] = this->getPoints(p);
    auto a = (p0->X + p1->X) / 2.;
    auto b = (p0->X + p2->X) / 2.;
    auto c = (p0->X + p1->X + p2->X) / 3.;
    return (TriangleArea(p0->X, a, b) + TriangleArea(a, c, b));
  };
  double getInscribedCircleArea() const noexcept {
    double a = std::get<0>(this->Lines)->length();
    double b = std::get<1>(this->Lines)->length();
    double c = std::get<2>(this->Lines)->length();
    // Hellon
    return std::sqrt((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) / 16.);
  };
  //////////////////////
  // このfaceの保存状況に従って，lに対する前後のpointへのインデックスを取得できる
  Tiii point_indicies(const netL* l) const {
    try {
      // for (auto i = 0; i < 3; i++)
      // 	if ((l->getPoints()[0] == this->Points[i] /*back*/ &&
      // std::get<1>(l->getPoints()) == this->Points[(i + 1) % 3] /*front*/) ||
      // 		(std::get<1>(l->getPoints()) == this->Points[i] /*back*/
      // && l->getPoints()[0] == this->Points[(i + 1) % 3] /*front*/))
      // return {i, (i + 1) % 3, (i + 2) % 3};

      auto p0 = std::get<0>(l->getPoints());
      auto p1 = std::get<0>(l->getPoints());
      if ((p0 == std::get<0>(this->Points) /*back*/ && p1 == std::get<1>(this->Points) /*front*/) || (p1 == std::get<0>(this->Points) /*back*/ && p0 == std::get<1>(this->Points) /*front*/))
        return {0, 1, 2};
      else if ((p0 == std::get<1>(this->Points) /*back*/ && p1 == std::get<2>(this->Points) /*front*/) || (p1 == std::get<1>(this->Points) /*back*/ && p0 == std::get<2>(this->Points) /*front*/))
        return {1, 2, 0};
      else if ((p0 == std::get<2>(this->Points) /*back*/ && p1 == std::get<0>(this->Points) /*front*/) || (p1 == std::get<2>(this->Points) /*back*/ && p0 == std::get<0>(this->Points) /*front*/))
        return {2, 0, 1};

      std::stringstream ss;
      ss << "このlineを基準としてインデックスをつくれない：この線はこの面のいっぺんではない"
         << " setBoundsを忘れていませんか？ input line l = " << l;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    }
  };
  // 与えられたpをindex[0]として，this->Pointsのインデックスを返す
  Tiii point_indicies(const networkPoint* const p) const {
    try {
      if (p == std::get<0>(this->Points))
        return {0, 1, 2};
      else if (p == std::get<1>(this->Points))
        return {1, 2, 0};
      else if (p == std::get<2>(this->Points))
        return {2, 0, 1};
      // for (auto i = 0; i < 3; i++)
      // 	if (p == this->Points[i])
      // 		return {i, (i + 1) % 3, (i + 2) % 3};
      std::stringstream ss;
      ss << "このpointを基準としてインデックスをつくれない：この線はこの面のいっぺんではない"
         << " setBoundsを忘れていませんか？ input point p = " << p;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    }
  };

  const Tddd& getAngles() const { return this->angles; };
  Tddd getAngles(const netL* l) const {
    try {
      if (l == std::get<0>(this->Lines))
        return this->angles;
      else if (l == std::get<1>(this->Lines))
        return RotateLeft(this->angles, 1);
      else
        return RotateLeft(this->angles, 2);
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
  };
  Tddd getAngles(const networkPoint* const p) const {
    try {
      auto [i, j, k] = point_indicies(p);
      auto [a0, a1, a2] = this->angles;
      if (i == 0)
        return this->angles;
      else if (i == 1)
        return {a1, a2, a0};
      else
        return {a2, a0, a1};
      // return {this->angles[i], /*ここにこの線lが位置する*/ this->angles[j],
      // this->angles[k]};
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
  };
  double getAngle(const networkPoint* p) const {
    if (p == std::get<0>(this->Points))
      return std::get<0>(this->angles);
    else if (p == std::get<1>(this->Points))
      return std::get<1>(this->angles);
    else if (p == std::get<2>(this->Points))
      return std::get<2>(this->angles);
    else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point is not found");
  };

  /* ------------------------------------------------------ */

  /*
  PointsはsetBoundsと同時にかならず，更新される．
  そのため，問題なくりようできる．
  ただし，位置関係を決めておかないと問題が生じる．
  */
  // const V_netPp &getPoints() const { return this->Points; };

  /* --------------------------------- SURFACE --------------------------------
   */

  bool BoundaryQ() const noexcept { return this->Tetras[0] == nullptr || this->Tetras[1] == nullptr; };

  bool ComputationalBoundaryQ() const noexcept;
  networkTetra* getInsideTetra() const noexcept;
  networkTetra* getOutsideTetra() const noexcept;

  /* ------------------------------------------------------ */
  //% タプル
  /* memo
  @ this->Points = {p0,p1,p2}
  @ this->Lines = {l0,l1,l2}
  @ [p0,p1,p2] = getPoints()
  @ [l0,l1,l2] = getLines()
  @ の関係:
  @ 		   p2
  @         /\
  @        /a2\
  @   l2  /    \ l1
  @      /a0  a1\
  @      ----------
  @   p0     l0    p1
  */
  /* --------------------------------------------------------------------------
   */
  const std::array<networkPoint*, 3>& getPoints() const noexcept { return this->Points; };

  std::array<networkPoint*, 3> getPoints(const networkPoint* const p) const noexcept {
    if (p == std::get<0>(this->Points))
      return this->Points;
    else if (p == std::get<1>(this->Points))
      return {std::get<1>(this->Points), std::get<2>(this->Points), std::get<0>(this->Points)};
    else
      return {std::get<2>(this->Points), std::get<0>(this->Points), std::get<1>(this->Points)};
  };
  std::array<networkPoint*, 3> getPoints(const networkLine* const l) const noexcept {
    if (l == std::get<1>(this->PLPLPL))
      return this->Points;
    else if (l == std::get<3>(this->PLPLPL))
      return {std::get<1>(this->Points), std::get<2>(this->Points), std::get<0>(this->Points)};
    else
      return {std::get<2>(this->Points), std::get<0>(this->Points), std::get<1>(this->Points)};
  };
  /* --------------------------------------------------------------------------
   */
  const T_LLL& getLines() const noexcept { return this->Lines; };
  T_LLL getLines(const networkLine* const l) const noexcept {
    if (l == std::get<0>(this->Lines))
      return this->Lines;
    else if (l == std::get<1>(this->Lines))
      return RotateLeft(this->Lines,
                        1); //{std::get<1>(this->Lines), std::get<2>(this->Lines),
                            // std::get<0>(this->Lines)};
    else
      return RotateLeft(this->Lines, 2);
  };
  /* --------------------------------------------------------------------------
   */
  const T_TT& getTetras() const noexcept { return this->Tetras; };
  /* --------------------------------------------------------------------------
   */
  bool replace(netL* const oldL, netL* const newL) {
    bool found = false;
    if (std::get<0>(this->Lines) == oldL) {
      std::get<0>(this->Lines) = newL;
      found = true;
    } else if (std::get<1>(this->Lines) == oldL) {
      std::get<1>(this->Lines) = newL;
      found = true;
    } else if (std::get<2>(this->Lines) == oldL) {
      std::get<2>(this->Lines) = newL;
      found = true;
    }
    if (found) {
      if (std::get<1>(this->PLPLPL) == oldL)
        std::get<1>(this->PLPLPL) = newL;
      else if (std::get<3>(this->PLPLPL) == oldL)
        std::get<3>(this->PLPLPL) = newL;
      else if (std::get<5>(this->PLPLPL) == oldL)
        std::get<5>(this->PLPLPL) = newL;
    }
    return found;
  };
  bool replace(netP* const oldP, netP* const newP);  // implemented in networkFace.cpp
  std::vector<networkFace*> getNeighbors(const networkFace* const obj) const noexcept { return {(*std::get<0>(this->Lines))(obj), (*std::get<1>(this->Lines))(obj), (*std::get<2>(this->Lines))(obj)}; };
  /* ------------------------------------------------------ */
  std::tuple<T_PPP, T_LLL> getPointsLinesTuple(const networkLine* const l) const { return {this->getPoints(l), getLines(l)}; };
  /* ------------------------------------------------------ */
  netP* getPointFront(const netL* l) const noexcept { return std::get<1>(getPoints(l)); };
  netP* getPointBack(const netL* l) const noexcept { return std::get<0>(getPoints(l)); };
  netP* getPointOpposite(const netL* l) const noexcept { return std::get<2>(getPoints(l)); };
  /* ------------------------------------------------------ */
  netL* getLineFront(const networkPoint* p) const noexcept { return std::get<0>(getLinesTupleFrom(p)); };
  netL* getLineBack(const networkPoint* p) const noexcept { return std::get<2>(getLinesTupleFrom(p)); };
  netL* getLineOpposite(const networkPoint* p) const noexcept { return std::get<1>(getLinesTupleFrom(p)); };
  /* ------------------------------------------------------ */

  netL* getLine(const netL* l, int j = 0) const {
    auto it = std::find(this->Lines.begin(), this->Lines.end(), l);
    if (it != this->Lines.end()) {
      int j_3 = j % 3; // should be [2,1,0,-1,-2]
      j = (j_3 == -1 ? 2 : (j_3 == -2 ? 1 : j_3));
      return this->Lines[(std::distance(Lines.begin(), it) + j) % 3];
    } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Line not found");
  };

  netL* getLineFront(const netL* l) const noexcept { return getLine(l, 1); };
  netL* getLineBack(const netL* l) const noexcept { return getLine(l, -1); };
  netF* getFaceFront(const netL* l) { return (*getLine(l, 1))(this); };
  netF* getFaceBack(const netL* l) { return (*getLine(l, -1))(this); };
  T_LLL getLinesTupleFrom(const networkPoint* const p) const {
    try {
      if (p == std::get<0>(this->Points))
        return this->Lines;
      else if (p == std::get<1>(this->Points))
        return RotateLeft(this->Lines, 1);
      else if (p == std::get<2>(this->Points))
        return RotateLeft(this->Lines, 2);
      else
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
  };
  //------------------------
  // netF *getFaceFront(const networkPoint *p) const { return
  // (*(getLinesFrom(p)[0]))(this); }; netF *getFaceBack(const networkPoint *p)
  // const { return (*(getLinesFrom(p)[0]))(this); }; netF
  // *getFaceOpposite(const networkPoint *p) const { return
  // (*(getLinesFrom(p)[0]))(this); };
  netF* getFaceFront(const networkPoint* p) const noexcept { return (*(std::get<0>(getLinesTupleFrom(p))))(this); };
  netF* getFaceBack(const networkPoint* p) const noexcept { return (*(std::get<2>(getLinesTupleFrom(p))))(this); };
  netF* getFaceOpposite(const networkPoint* p) const noexcept { return (*(std::get<1>(getLinesTupleFrom(p))))(this); };

  /* --------------------------------------------------------------------------
   */

  template <size_t N>
  std::array<networkPoint*, N> getPoints() const {
    static_assert(N == 3 || N == 6, "Unsupported getPoints. Only 3 or 6 are supported.");
    if constexpr (N == 3) {
      return ToArray(this->Points);
    } else if constexpr (N == 6) {
      auto [l0, l1, l2] = this->getLines();
      auto [p0_f0, p1_f0, p2_f0] = (*l1)(this)->getPoints(l1); // f0
      auto [p0_f1, p1_f1, p2_f1] = (*l2)(this)->getPoints(l2); // f1
      auto [p0_f2, p1_f2, p2_f2] = (*l0)(this)->getPoints(l0); // f2
      return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
    }
  }
  //
  template <size_t N, bool networkLine::* BoolMember = nullptr>
  std::array<networkPoint*, N> getPoints(const networkPoint* P = nullptr) const {
    static_assert(N == 3 || N == 6, "Unsupported getPoints. Only 3 or 6 are supported.");
    if constexpr (N == 3) {
      if (P == std::get<0>(this->Points) || P == nullptr)
        return ToArray(this->Points);
      else if (P == std::get<1>(this->Points))
        return ToArray(RotateLeft(this->Points, 1));
      else
        return ToArray(RotateLeft(this->Points, 2));
    } else if constexpr (N == 6) {
      /*              0
       *              0  p2_f0
       *            /  \
       *   p0_f0  3/-l1-\5  p1_f0
       *         /  \  /  \
       *  1    1/---4(p)---\2   2　p2_f2
       */
      auto [l0, l1, l2] = this->getLinesTupleFrom(P);
      auto [p0_f0, p1_f0, p2_f0] = (*l1)(this)->getPoints(l1); // f0
      auto [p0_f1, p1_f1, p2_f1] = (*l2)(this)->getPoints(l2); // f1
      auto [p0_f2, p1_f2, p2_f2] = (*l0)(this)->getPoints(l0); // f2
      if (BoolMember) {
        if (l0->*BoolMember)
          p2_f2 = nullptr;
        if (l1->*BoolMember)
          p2_f0 = nullptr;
        if (l2->*BoolMember)
          p2_f1 = nullptr;
      }
      return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
    }
  }

  template <size_t N, bool networkLine::* BoolMember = nullptr>
  std::array<networkPoint*, N> getPoints(const networkLine* L = nullptr) const {
    static_assert(N == 3 || N == 6, "Unsupported getPoints. Only 3 or 6 are supported.");
    if constexpr (N == 3) {
      return ToArray(this->getPoints(L));
    } else if constexpr (N == 6) {
      auto [l1, l2, l0] = this->getLines(L);
      auto [p0_f0, p1_f0, p2_f0] = (*l1)(this)->getPoints(l1); // f0
      auto [p0_f1, p1_f1, p2_f1] = (*l2)(this)->getPoints(l2); // f1
      auto [p0_f2, p1_f2, p2_f2] = (*l0)(this)->getPoints(l0); // f2
      if (BoolMember) {
        if (l0->*BoolMember)
          p2_f2 = nullptr;
        if (l1->*BoolMember)
          p2_f0 = nullptr;
        if (l2->*BoolMember)
          p2_f1 = nullptr;
      }
      return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
    }
  }

  /* --------------------------------------------------------------------------
   */

  std::tuple<networkPoint*, networkLine*, networkPoint*, networkLine*, networkPoint*, networkLine*> getPointsAndLines(const networkPoint* origin) const {
    auto [P0, L0, P1, L1, P2, L2] = this->PLPLPL;
    if (P0 == origin || nullptr == origin) {
      return this->PLPLPL;
    } else if (P1 == origin) {
      return {P1, L1, P2, L2, P0, L0};
    } else if (P2 == origin) {
      return {P2, L2, P0, L0, P1, L1};
    } else
      return this->PLPLPL;
  };

  std::tuple<networkPoint*, networkLine*, networkPoint*, networkLine*, networkPoint*, networkLine*> getPointsAndLinesFromNearest(const networkPoint* origin) const {
    auto [P0, L0, P1, L1, P2, L2] = this->PLPLPL;
    if (P0 == origin || nullptr == origin) {
      return this->PLPLPL;
    } else if (P1 == origin || (Norm(origin->X - P0->X) > Norm(origin->X - P1->X) && Norm(origin->X - P2->X) > Norm(origin->X - P1->X))) {
      return {P1, L1, P2, L2, P0, L0};
    } else if (P2 == origin || (Norm(origin->X - P0->X) > Norm(origin->X - P2->X) && Norm(origin->X - P1->X) > Norm(origin->X - P2->X))) {
      return {P2, L2, P0, L0, P1, L1};
    } else
      return this->PLPLPL;
  };

  std::tuple<networkPoint*, networkLine*, networkPoint*, networkLine*, networkPoint*, networkLine*> getPointsAndLines(const networkLine* origin) const {
    auto [P0, L0, P1, L1, P2, L2] = this->PLPLPL;
    if (L0 == origin || nullptr == origin) {
      return this->PLPLPL;
    } else if (L1 == origin) {
      return {P1, L1, P2, L2, P0, L0};
    } else if (L2 == origin) {
      return {P2, L2, P0, L0, P1, L1};
    } else
      return this->PLPLPL;
  };

  /* --------------------------------------------------------------------------
   */

  std::tuple<netPp, netPp, netPp, netPp, netPp, netPp> get6PointsTuple(const networkLine* l /*基準*/) const {
    try {
      auto [l0, l1, l2] = this->getLines(l); // 修正した2022/03/21
      // if (l1 == l)
      // {
      // 	l0 = std::get<1>(this->Lines);
      // 	l1 = std::get<2>(this->Lines);
      // 	l2 = std::get<0>(this->Lines);
      // }
      // else if (l2 == l)
      // {
      // 	l0 = std::get<2>(this->Lines);
      // 	l1 = std::get<0>(this->Lines);
      // 	l2 = std::get<1>(this->Lines);
      // }

      /*
       *                  0  (p2_f0)
       *                 / \
       *                /f0 \
       *               /     \
       *    (p0_f0)  3/--l-l0-\5 (p0_f2)
       *             / \     / \
       *            /   l1  l2  \
       *           / f1  \ /  f2 \
       * (p2_f1) 1/-------4-------\2  (p2_f2)
       *               (p0_f1)
       */
      auto [p0_f0, p1_f0, p2_f0] = (*l0)(this)->getPoints(l0); // f0
      auto [p0_f1, p1_f1, p2_f1] = (*l1)(this)->getPoints(l1); // f1
      auto [p0_f2, p1_f2, p2_f2] = (*l2)(this)->getPoints(l2); // f2
      return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
  };

  std::tuple<netPp, netPp, netPp, netPp, netPp, netPp> get6PointsTuple(const networkPoint* p) const {
    try {
      if (p == nullptr || !MemberQ_(this->Points, p)) {
        p = std::get<0>(this->Points);
      }
      auto [l0, l1, l2] = this->getLinesTupleFrom(p);
      /*              0
       *              0
       *            /  \
       *          3/-l1-\5
       *         /  \  /  \
       *  1    1/---4(p)---\2   2
       */
      return (*l1)(this)->get6PointsTuple(l1);
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
  };

  //@ 2024/02/02本格的に使い始める
  std::array<networkPoint*, 6> get6Points(const networkPoint* p) const {
    if (p == nullptr || !MemberQ_(this->Points, p))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "p is nullptr or not found in this->Points");
    /*
     *                  0  (p2_f0)
     *                 / \
     *                /f0 \
     *               /     \
     *    (p0_f0)  3/--l-l0-\5 (p0_f2)
     *             / \     / \
     *            /   l1  l2  \
     *           / f1  \ /  f2 \
     * (p2_f1) 1/-------4-------\2  (p2_f2)
     *               (p0_f1)
     */
    auto [l0, l1, l2] = this->getLinesTupleFrom(p);
    auto [p0_f0, p1_f0, p2_f0] = (*l0)(this)->getPoints(l0); // f0
    auto [p0_f1, p1_f1, p2_f1] = (*l1)(this)->getPoints(l1); // f1
    auto [p0_f2, p1_f2, p2_f2] = (*l2)(this)->getPoints(l2); // f2

    if (l0->CORNER)
      p2_f0 = nullptr;
    if (l1->CORNER)
      p2_f1 = nullptr;
    if (l2->CORNER)
      p2_f2 = nullptr;
    return {p2_f0, p2_f1, p2_f2, p0_f0, p0_f1, p0_f2};
  };

  //-------------------------
  Tddd getMeanX() const { return Mean(getLocationsTuple()); };

  T3Tddd getLocationsTuple() const { return {std::get<0>(this->Points)->X, std::get<1>(this->Points)->X, std::get<2>(this->Points)->X}; };
  T3Tddd getXVertices() const { return {std::get<0>(this->Points)->X, std::get<1>(this->Points)->X, std::get<2>(this->Points)->X}; };

  netFp Split(netLp DivL, netLp newDivL, netLp newMidL, int type);

  /* ------------------------------------------------------ */
  bool isFacing(const networkFace* const F, const double rad = 1E-10) const { return isFlat(this->normal, -F->normal, rad) || isFlat(this->normal, F->normal, rad); };
};

/* -------------------------------------------------------------------------- */
/*                                 クラスに属さない関数                                 */
/* -------------------------------------------------------------------------- */

bool isEdge(const networkFace* const f);

/* -------------------------------------------------------------------------- */
T3Tddd ToX(const networkFace* const f);
std::vector<T3Tddd> ToX(const std::vector<networkFace*>& fs);
std::vector<T3Tddd> ToX(const std::unordered_set<networkFace*>& fs);
/* -------------------------------------------------------------------------- */
Tddd TriangleNormal(const networkFace* const f);
//@ ------------------------ 抽出用関数など ----------------------- */
std::vector<Tddd> extX(const std::unordered_set<networkFace*>& fs);
std::vector<T3Tddd> extVertices(const std::unordered_set<networkFace*>& fs);
std::vector<T2Tddd> extX(const std::unordered_set<networkLine*>& ls);
std::vector<Tddd> extNormals(const V_netFp& fs);

//@ ------------------------------------------------------ */

std::vector<netL*> link(const V_netPp& obj, Network* net);
netL* unlink(netP* obj, netP* obj_);

//@ ------------------------------------------------------ */
//@                         extract                        */
//@ ------------------------------------------------------ */

std::unordered_set<networkPoint*> extPointsCORNER_(const std::vector<networkPoint*>& ps);
std::unordered_set<networkPoint*> extPointsCORNER_(const std::unordered_set<networkPoint*>& ps);
std::unordered_set<networkPoint*> extPointsCornerOrNeumann_(const std::vector<networkPoint*>& ps);
std::unordered_set<networkPoint*> extPointsCornerOrNeumann_(const std::unordered_set<networkPoint*>& ps);

// b! ------------------------------------------------------ */
// b! ---------------------- extLines ---------------------- */
// b! ------------------------------------------------------ */
/* ------------------ for unordered_set ----------------- */
std::unordered_set<networkLine*> extLinesCORNER_(const std::unordered_set<networkFace*>& fs);
std::unordered_set<networkLine*> extLines_(const std::unordered_set<networkFace*>& fs);
std::unordered_set<networkLine*> extLinesCORNER_(const std::unordered_set<networkPoint*>& ps);
std::unordered_set<networkLine*> extLinesCORNER_(const networkPoint* p);
std::unordered_set<networkLine*> extLines_(const std::unordered_set<networkPoint*>& ps);
/* --------------------- for vector --------------------- */
std::unordered_set<networkLine*> extLinesCORNER_(const std::vector<networkFace*>& fs);
std::unordered_set<networkLine*> extLines_(const std::vector<networkFace*>& fs);
std::unordered_set<networkLine*> extLinesCORNER_(const std::vector<networkPoint*>& ps);
std::unordered_set<networkLine*> extLines_(const std::vector<networkPoint*>& ps);
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */
// b! ------------------------------------------------------ */
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
V_d extLength(const std::unordered_set<networkLine*>& ls);
V_d extLength(const V_netLp& ls);
Tddd extLength(const std::array<networkLine*, 3>& ls);
V_d extAreas(const V_netFp& fs);
V_d extractAreas(const V_netFp& fs);

template <class T>
std::vector<Tddd> extractXtuple(const std::unordered_set<T*>& object) {
  std::vector<Tddd> ret(object.size());
  int i = 0;
  for (auto it = object.begin(); it != object.end(); ++it)
    ret[i++] = (*it)->X;
  return ret;
};
template <class T>
std::vector<Tddd> extractXtuple(const std::vector<T*>& object) {
  std::vector<Tddd> ret(object.size());
  for (auto i = 0; i < object.size(); i++)
    ret[i] = object[i]->X;
  return ret;
};

template <class T>
VV_d extractX(const std::unordered_set<T*>& object) {
  VV_d ret(object.size(), V_d(3));
  int i = 0;
  for (auto it = object.begin(); it != object.end(); ++it)
    ret[i++] = ToVector((*it)->X);
  return ret;
};
template <class T>
VV_d extractX(const std::vector<T*>& object) {
  VV_d ret(object.size(), V_d(3));
  for (auto i = 0; i < object.size(); i++)
    ret[i] = ToVector(object[i]->X);
  return ret;
};
// 2021/09/06追加
std::vector<Tddd> extXtuple(const V_netPp& points);
std::vector<Tddd> extXtuple(const V_netFp& points);
T3Tddd extractXtuple(networkFace const* f);
std::vector<std::vector<Tddd>> extractXtuple(const std::vector<networkLine*>& lines);
/* -------------------------------------------------------------------------- */

#include "searcher.hpp"

/* -------------------------------------------------------------------------- */
/*      @     */
/*     /|\    */
/*    / | \   */
/*   /  |  \  */
/*  @---@---@ */
/* -------------------------------------------------------------------------- */

inline double det3(const Tddd& u, const Tddd& v, const Tddd& w) noexcept { return CrossDot(u, v, w); }

inline double orient3d(const Tddd& a, const Tddd& b, const Tddd& c, const Tddd& d) noexcept { return det3(a - d, b - d, c - d); }

inline double insphere(const Tddd& a, const Tddd& b, const Tddd& c, const Tddd& d, const Tddd& e) noexcept {
  const auto ea = a - e;
  const auto eb = b - e;
  const auto ec = c - e;
  const auto ed = d - e;

  const double ea2 = Dot(ea, ea);
  const double eb2 = Dot(eb, eb);
  const double ec2 = Dot(ec, ec);
  const double ed2 = Dot(ed, ed);

  const double det = det3(ea, eb, ec) * ed2 - det3(ea, eb, ed) * ec2 + det3(ea, ec, ed) * eb2 - det3(eb, ec, ed) * ea2;

  const double orient = orient3d(a, b, c, d);
  return orient > 0. ? det : -det;
}

class networkTetra : public Tetrahedron, public ElasticBodyDynamics {
public:
  Network* network;
  // Network connectivities
  T_4P Points;
  T_6L Lines;
  T_4F Faces;
  TetraState tetra_state = TetraState::Water;

  double mass = 0;
  networkTetra(Network* const network_IN, const T_4P& points, const T_6L& lines, const T_4F& faces);
  ~networkTetra();

  bool IsolatedQ() const noexcept {
    return std::ranges::all_of(Faces, [](auto* f) { return f->BoundaryQ(); });
  }

  networkFace* faceOpposite(const networkPoint* point) const noexcept {
    for (auto* face : Faces) {
      if (!face)
        continue;
      const auto pts = face->getPoints();
      if (std::ranges::any_of(pts, [](auto* fp) { return fp == nullptr; }))
        continue;
      if (!std::ranges::any_of(pts, [&](auto* fp) { return fp == point; }))
        return face;
    }
    return nullptr;
  }

  /* --------------------------------------------------------------------------
   */

  std::array<double, 3> grad(const std::function<double(networkPoint*)>& func) {
    auto [p0, p1, p2, p3] = Points;
    std::array<double, 4> values = {func(p0), func(p1), func(p2), func(p3)}, ans;
    std::array<std::array<double, 4>, 4> positions = {{{1., p0->X[0], p0->X[1], p0->X[2]}, {1., p1->X[0], p1->X[1], p1->X[2]}, {1., p2->X[0], p2->X[1], p2->X[2]}, {1., p3->X[0], p3->X[1], p3->X[2]}}};
    lapack_svd_solve(positions, ans, values);
    return {ans[1], ans[2], ans[3]};
  };

  std::array<std::array<double, 3>, 3> grad(const std::function<std::array<double, 3>(networkPoint*)>& func) {
    auto [p0, p1, p2, p3] = Points;
    std::array<std::array<double, 4>, 4> positions = {{{1., p0->X[0], p0->X[1], p0->X[2]}, {1., p1->X[0], p1->X[1], p1->X[2]}, {1., p2->X[0], p2->X[1], p2->X[2]}, {1., p3->X[0], p3->X[1], p3->X[2]}}};
    auto ans = Dot(Inverse(positions), std::array<std::array<double, 3>, 4>{{func(p0), func(p1), func(p2), func(p3)}});
    return {ans[1], ans[2], ans[3]};
  };

  // std::array<std::array<double, 3>, 3> grad(const
  // std::function<std::array<double, 3>(networkPoint *)> &func) {
  //    return {grad([&func](networkPoint *p) -> double { return func(p)[0]; }),
  //            grad([&func](networkPoint *p) -> double { return func(p)[1]; }),
  //            grad([&func](networkPoint *p) -> double { return func(p)[2];
  //            })};
  // };

  double div(const std::function<std::array<double, 3>(networkPoint*)>& func) {
    //@ func returns {f0,f1,f2}
    //@ div . {f0,f1,f2} = f0x + f1y + f2z
    //@ which is tr(grad(func))
    return Trace(this->grad(func));
  };

  std::array<double, 3> div(const std::function<std::array<std::array<double, 3>, 3>(networkPoint*)>& func) {
    //@ func returns {f0,f1,f2}
    //@ div . {f0,f1,f2} = f0x + f1y + f2z
    //@ which is tr(grad(func))
    return {this->div([&func](networkPoint* p) -> std::array<double, 3> { return func(p)[0]; }),
            this->div([&func](networkPoint* p) -> std::array<double, 3> { return func(p)[1]; }),
            this->div([&func](networkPoint* p) -> std::array<double, 3> { return func(p)[2]; })};
  };

  std::array<networkPoint*, 4> getPoints(const networkPoint* p) const {
    if (p == std::get<1>(this->Points))
      return {Points[1], Points[2], Points[3], Points[0]};
    else if (p == std::get<2>(this->Points))
      return {Points[2], Points[3], Points[0], Points[1]};
    else if (p == std::get<3>(this->Points))
      return {Points[3], Points[0], Points[1], Points[2]};
    else
      return this->Points;
  };

  bool DelaunayQ() {
    auto [p0, p1, p2, p3] = this->Points;
    return insphere(p0->X, p1->X, p2->X, p3->X, Mean(p0->X, p1->X, p2->X, p3->X)) < 0;
  };
};

/* -------------------------------------------------------------------------- */

netF* genFace(Network* const net, netL* const l0, netL* const l1, netL* const l2);
std::tuple<bool, networkTetra*> genTetra(Network* const net, netP* const p0, netP* const p1, netP* p2, netP* p3);
std::tuple<bool, networkTetra*> genTetra(Network* const net, const T_4P& abcd);
Tddd Nearest(const Tddd& X, const networkFace* f);
double Distance(const Tddd& X, const networkFace* f);
std::tuple<Tddd, networkFace*> Nearest_(const Tddd& X, const std::unordered_set<networkFace*>& faces);
std::tuple<Tddd, networkFace*> Nearest_(const Tddd& X, const std::vector<networkFace*>& faces);
std::tuple<Tddd, networkFace*> Nearest_(const Tddd& X, const double& r, const std::unordered_set<networkFace*>& faces);
Tddd Nearest(const Tddd& X, const std::vector<networkFace*>& faces);
Tddd Nearest(const Tddd& X, const std::unordered_set<networkFace*>& faces);
Tddd Nearest(const networkPoint* p, const std::unordered_set<networkFace*>& faces);
double DirectedHausdorffDistance(const V_netFp& query_faces, const V_netFp& target_faces, int triangle_subdivision = 6);
double DirectedHausdorffDistance(const V_netFp& query_faces, const V_Netp& target_networks, int triangle_subdivision = 6);

/* -------------------------------------------------------------------------- */

#ifdef tetgenH

inline tetgenio generate_tetgenio_input(const auto& surfaces) {
  std::cout << Red << "generate_tetgenio_input" << colorReset << std::endl;

  /* -------- distribute index ----------*/

  // surfacesが参照する全頂点を一意に集める（tetgenへ渡す点集合）
  std::unordered_set<networkPoint*> tmp;
  for (const auto& f : surfaces)
    for (const auto& p : f->getPoints())
      tmp.insert(p);

  // 重複排除済み集合から順序付きベクタへ移し替え
  const std::vector<networkPoint*> points(std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));

  // tetgen用のpointインデックスを割り当てる
  int count = 0;
  for (const auto& p : points)
    p->index = count++;
  // こちらは面のインデックス割り当て（後続で参照しやすくするため）
  count = 0;
  for (const auto& f : surfaces)
    f->index = count++;

  /* -------- tetgenio in の作成 -------- */

  tetgenio in;

  // pointlistに3次元座標を詰める（tetgen形式：xyzxyz...）
  in.pointlist = (REAL*)malloc(points.size() * 3 * sizeof(REAL));
  if (!in.pointlist) {
    printf("Memory allocation error for in.pointlist\n");
    throw std::runtime_error("Memory allocation error for in.pointlist");
  }
  // 事前に割り当てたインデックスを利用
  for (const auto& p : points) {
    in.pointlist[p->index * 3] = p->X[0];
    in.pointlist[p->index * 3 + 1] = p->X[1];
    in.pointlist[p->index * 3 + 2] = p->X[2];
  }

  // 面集合をtetgenのfacet配列へコピー（各三角形→単一ポリゴン）
  in.facetlist = (tetgenio::facet*)malloc(surfaces.size() * sizeof(tetgenio::facet));
  if (!in.facetlist) {
    printf("Memory allocation error for in.facetlist\n");
    free(in.pointlist);
    throw std::runtime_error("Memory allocation error for in.facetlist");
  }

  // 面データを in.facetlist に保存
  for (const auto surface : surfaces) {
    auto f = &in.facetlist[surface->index];    // tetgenのfacet配列の該当面
    in.init(f);                                // facetの初期化
    f->numberofpolygons = 1;                   // 各facetは単一ポリゴン（三角形）で構成される
    f->polygonlist = new tetgenio::polygon[1]; // ポリゴン配列の確保
    auto p = &f->polygonlist[0];               // ポリゴンデータへのポインタ
    in.init(p);                                // ポリゴンの初期化
    p->numberofvertices = 3;                   // 三角形なので頂点数は3
    const auto& [p0, p1, p2] = surface->getPoints();
    p->vertexlist = new int[3]; // 頂点インデックス配列の確保
    p->vertexlist[0] = p0->index;
    p->vertexlist[1] = p1->index;
    p->vertexlist[2] = p2->index;
  }

  // tetgenへ渡すメッシュサイズ情報
  in.numberofpoints = points.size();
  in.numberoffacets = surfaces.size();

  std::cout << Green << "generate_tetgenio_input: done" << colorReset << std::endl;

  return in;
}

// 点群からtetgenio入力を生成する関数
inline tetgenio generate_tetgenio_input_from_points(const std::vector<networkPoint*>& points) {
  std::cout << Red << "generate_tetgenio_input_from_points" << colorReset << std::endl;

  tetgenio in;
  in.firstnumber = 0;

  // pointlistに3次元座標を詰める
  in.numberofpoints = points.size();
  in.pointlist = (REAL*)malloc(in.numberofpoints * 3 * sizeof(REAL));
  if (!in.pointlist) {
    throw std::runtime_error("Memory allocation failed for pointlist");
  }

  // インデックスを割り当てつつ座標をコピー
  int count = 0;
  for (const auto& p : points) {
    p->index = count;
    in.pointlist[count * 3 + 0] = p->X[0];
    in.pointlist[count * 3 + 1] = p->X[1];
    in.pointlist[count * 3 + 2] = p->X[2];
    count++;
  }

  // 点群のみの場合、facetlist（面情報）は作成しない
  in.numberoffacets = 0;
  in.facetlist = nullptr;

  std::cout << Green << "generate_tetgenio_input_from_points: done" << colorReset << std::endl;
  return in;
}

// [追加] 点群と面群の両方からtetgenio入力を生成する関数 (PLC用)
inline tetgenio generate_tetgenio_input_mixed(const std::vector<networkPoint*>& points, const std::vector<networkFace*>& faces) {
  std::cout << Red << "generate_tetgenio_input_mixed" << colorReset << std::endl;

  tetgenio in;
  in.firstnumber = 0;

  // 1. Point List (全点)
  in.numberofpoints = points.size();
  in.pointlist = (REAL*)malloc(in.numberofpoints * 3 * sizeof(REAL));
  if (!in.pointlist) {
    throw std::runtime_error("Memory allocation failed for pointlist");
  }

  // インデックスを割り当てつつ座標をコピー
  int count = 0;
  for (const auto& p : points) {
    p->index = count;
    in.pointlist[count * 3 + 0] = p->X[0];
    in.pointlist[count * 3 + 1] = p->X[1];
    in.pointlist[count * 3 + 2] = p->X[2];
    count++;
  }

  // 2. Facet List (面)
  in.numberoffacets = faces.size();
  in.facetlist = (tetgenio::facet*)malloc(in.numberoffacets * sizeof(tetgenio::facet));
  if (!in.facetlist && in.numberoffacets > 0) {
    throw std::runtime_error("Memory allocation failed for facetlist");
  }

  for (size_t i = 0; i < faces.size(); ++i) {
    auto f = faces[i];
    tetgenio::facet* F = &in.facetlist[i];
    tetgenio::init(F);
    F->numberofpolygons = 1;
    F->polygonlist = (tetgenio::polygon*)malloc(F->numberofpolygons * sizeof(tetgenio::polygon));
    F->numberofholes = 0;
    F->holelist = nullptr;

    tetgenio::polygon* P = &F->polygonlist[0];
    tetgenio::init(P);
    P->numberofvertices = 3;
    P->vertexlist = (int*)malloc(P->numberofvertices * sizeof(int));

    auto [p0, p1, p2] = f->getPoints();
    P->vertexlist[0] = p0->index;
    P->vertexlist[1] = p1->index;
    P->vertexlist[2] = p2->index;
  }

  std::cout << Green << "generate_tetgenio_input_mixed: done" << colorReset << std::endl;
  return in;
}

#endif

/* -------------------------------------------------------------------------- */
/*                                    接続関係                                    */
/* -------------------------------------------------------------------------- */

networkFace* isFace(const networkPoint* const a, const networkPoint* const b, const networkPoint* const c);
networkLine* ConnectedQ(const networkFace* const a, const networkFace* const b);
networkLine* ConnectedQ(const networkLine* const l, const networkPoint* const a);
networkLine* ConnectedQ(const networkPoint* const a, const networkLine* const l);
networkLine* ConnectedQ(const networkPoint* const a, const networkPoint* const b);
T_3L ConnectedQ(const networkPoint* const a, const networkPoint* const b, const networkPoint* const c);
T_3L ConnectedQ(const T_3P& abc);

/* ========================================================================== */
/*   Networkは持っているPointsから情報を取り出すメソッドを提供する */
/*   @-@-@                                                                    */
/*   |\|/|                                                                    */
/*   @-@-@                                                                    */
/*   |/|\|                                                                    */
/*   @-@-@                                                                    */
/*  ========================================================================= */

class MooringLine;
class Network : public CoordinateBounds, public RigidBodyDynamics {
public:
  virtual ~Network(); // = default;  // Virtual destructor

  vtkPolygonWriter<networkPoint*> vtpWriter;

public:
  Network* surfaceNet;

  // b#
  // --------------------------------------------------------------------------
  // */ b#                                 octreeに関する */ b#
  // --------------------------------------------------------------------------
  // */
public:
  octree<networkFace*>* octreeOfFaces;
  octree<networkPoint*>* octreeOfPoints;
  // Geometry/topology signature cache for setGeometricProperties() skipping.
  // Kept here so any caller can rely on consistent detection.
  bool setGeometricProperties();
  void genOctreeOfFaces(const Tii& depthlimit, const int objnum) {
    std::cout << this->getName() << "genOctreeOfFaces(" << depthlimit << ", " << objnum << ")";
    this->setGeometricPropertiesForce();
    if (this->octreeOfFaces)
      delete this->octreeOfFaces;
    this->octreeOfFaces = new octree(this->scaledBounds(expand_bounds), depthlimit, objnum, this->Faces);
    Print(this->getName() + "->octreeOfFaces->deleteOuside()");
    this->octreeOfFaces->deleteOuside();
    Print(this->getName() + "->octreeOfFaces->setNeighbors()");
    this->octreeOfFaces->setNeighbors();
  };

  void genOctreeOfPoints(const Tii& depthlimit, const int objnum) {
    Print("octreeOfPointsを生成");
    this->setGeometricPropertiesForce();
    if (this->octreeOfPoints)
      delete this->octreeOfPoints;
    this->octreeOfPoints = new octree(this->scaledBounds(expand_bounds), depthlimit, objnum, this->Points);
    this->octreeOfPoints->deleteOuside();
  };

  /* --------------------------------------------------------------------------*/

  std::tuple<bool, octree<networkFace*>*, networkFace*> isInside_MethodOctree(const Tddd& X, const double small = 1E-13) const {
    //! セルに含まれるかどうかではなく，ポリゴンの内部にあるかどうかで判定するための関数
    if (this->octreeOfFaces) {
      auto cells = this->octreeOfFaces->getIntersectInside(X);
      if (cells.empty())
        return {false, nullptr, nullptr}; //! どのセルにもXは入らない
      else {
        auto cell = cells[0];
        if (cell->faces_.empty()) {
          //% cellがfaceを持たない場合，
          //% まずはセルの頂点における回転数を調べることで，
          //% ポリゴンの内部外部のどちらに位置するか調べる
          if (std::ranges::all_of(cell->WNs, [](const auto& w_num) { return w_num > 0.6; }) /*should be inside*/)
            return {true, cell, nullptr};
          else if (std::ranges::all_of(cell->WNs, [](const auto& w_num) { return w_num < 0.1; }) /*should be outside*/)
            return {false, cell, nullptr};
          else {
            //@ 構造物の内部外部に位置するか回転数で判断するのが困難な場合
            std::unordered_set<networkFace*> faces;
            for (const auto& N0 : cell->neighbors)
              for (const auto& nei : N0->neighbors)
                faces.insert(nei->faces_.begin(), nei->faces_.end());
            if (faces.empty())
              return {false, cell, nullptr};
            else {
              auto [nearest, f] = Nearest_(X, faces);
              return {(Dot(nearest - X, f->normal) >= 0. || Norm(nearest - X) < small), cell, f};
            }
          }
        } else {
          //% faceを持つ場合，最も近い面とdotが同じかどうか
          auto [nearest, f] = Nearest_(X, cell->faces_);
          return {(Dot(nearest - X, f->normal) >= 0. || Norm(nearest - X) < small), cell, f};
        }
      }
    } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "octreeOfFaces has not been set");
  };

  Tddd interpolateVector(const Tddd& X) const {
    auto [isinside, cell, _] = isInside_MethodOctree(X);
    if (cell)
      return cell->Interpolate(X, cell->vectors);
    else
      return {0., 0., 0.};
  };

  // --------------------------------------------------------------------------

public:
  RungeKutta<Tddd> RK_COM;
  RungeKutta<T4d> RK_Q;
  RungeKutta<T6d> RK_Velocity;

public:
  bool IGNORE;
  // Tdd grid_pull_factor;
  int grid_pull_depth;
  /* ------------------------------------------------------ */
  /*                          体積の計算                      */
  /* ------------------------------------------------------ */
  double getVolume() const;
  double getVolume(std::function<std::array<double, 3>(const networkPoint*)> func_X) const;
  double GaussIntegral2(const Tddd& X) const;
  double GaussIntegral(const Tddd& X) const;
  double windingNumber(const Tddd& X) const;
  bool InsideQ(const Tddd& X) const;
  std::vector<double> windingNumber(const std::vector<Tddd>& Xs) const;
  bool all_of_isInside(const std::vector<Tddd>& Xs) const;
  //! ------------------------------------------------------ */
  //!                          接触の判別                      */
  //! ------------------------------------------------------ */

  std::unordered_set<networkFace*> getContactFacesOfPoints() const;
  std::unordered_set<networkPoint*> getContactPointsOfPoints() const;
  std::unordered_set<networkPoint*> getContactPointsOfPoints(const std::vector<Network*>& nets) const;

  // b$ =================================================================== */
  // b$       ３次元空間分割に関する．　バケツ．Faces,Points,ParametricPoints
  // b$ =================================================================== */

public:
  Buckets<networkFace*> BucketFaces, BucketSurfaces;
  Buckets<networkPoint*> BucketPoints;
  Buckets<networkTetra*> BucketTetras;

  const double expand_bounds = 1.5;

  /* --------------------------------------------------------------------------*/

  void makeBuckets(const double spacing = 1E+20);

  void makeBucketTetras(double spacing = 1E+20);

  /*面は点と違って，複数のバケツ（セル）と接することがある*/
  void makeBucketFaces(double spacing = 1E+20);

  double last_makeBucketPoints_spacing = 0.;

  void makeBucketPoints(double spacing = 1E+20);

  /* --------------------------------------------------------------------------*/

  std::vector<networkFace*> ContactFaces;

  void setContactFaces(const std::vector<Network*>& nets) {

    std::cout << this->getName() << " setContactFaces()" << std::endl;
    auto points = ToVector(this->Points);
    /*
     *  p->contact_range
     *  p->ContactFaces
     *  this->ContactFaces
     を設定する．
     */

    this->setGeometricPropertiesForce();

    for (const auto& p : points) {
      p->setContactRange(nets);
      p->clearContactFaces();
    }

    _Pragma("omp parallel for") for (const auto& p : points)
        p->addContactFaces(nets, false);
    // Lines の接触判定（端点の contact_range 設定後に実行）
    auto lines = ToVector(this->getLines());
    for (const auto& l : lines) {
      l->setContactRange();
      l->clearContactFaces();
    }

    _Pragma("omp parallel for") for (const auto& l : lines)
        l->addContactFaces(nets, false);

    std::vector<networkFace*> valid_faces;
    for (const auto& net : nets) {
      const auto boundary_faces = net->getBoundaryFaces();
      valid_faces.insert(valid_faces.end(), boundary_faces.begin(), boundary_faces.end());
    }
    auto is_valid_face = [&](networkFace* f) {
      return f != nullptr && std::find(valid_faces.begin(), valid_faces.end(), f) != valid_faces.end();
    };
    std::vector<networkFace*> faces;
    auto append_unique_face = [&](networkFace* f) {
      if (!is_valid_face(f))
        return;
      if (std::find(faces.begin(), faces.end(), f) == faces.end())
        faces.emplace_back(f);
    };
    for (const auto& p : points)
      for (auto* cf : p->getContactFaces())
        append_unique_face(cf);
    for (const auto& l : lines)
      for (auto* cf : l->getContactFaces())
        append_unique_face(cf);
    this->ContactFaces = std::move(faces);
  };

  // b$

  // ----------------------------------------------------------------------------------------- */

  /*
  Nearestは，外部のNearest関数と似たように，Networkの境界面上の最も近い点を返す．
  ただし，Network内の面に限定して探索を行う．
  */
  std::tuple<networkFace*, std::array<double, 3>> Nearest(const std::array<double, 3>& X, const std::function<Tddd(const networkPoint*)>& position = [](const networkPoint* p) { return p->X; }) const {
    double min_distance = 1E+20, r;
    std::tuple<networkFace*, std::array<double, 3>> ret = {nullptr, {1E+20, 1E+20, 1E+20}};
    T3Tddd X0X1X2;
    Tddd X_nearest;
    for (const auto& f : this->getBoundaryFaces()) {
      auto [p0, p1, p2] = f->getPoints();
      X0X1X2 = {position(p0), position(p1), position(p2)};
      if (min_distance > (r = (Norm((X_nearest = ::Nearest(X, X0X1X2)) - X)))) {
        min_distance = r;
        ret = {f, X_nearest};
      }
    }
    return ret;
  }

  // query_faces 上に一様サンプル点を置き，この Network の境界面への最近点距離の最大値を返す．
  double DirectedHausdorffDistance(const V_netFp& query_faces, const int triangle_subdivision = 6) const {
    return ::DirectedHausdorffDistance(query_faces, this->getBoundaryFaces(), triangle_subdivision);
  }

  // @ ============================================================== */
  // @                 RigidBodyDynamicsを使ったメソッド                 */
  // @ ============================================================== */

  void RigidBodyUpdatePoints() {
    for (const auto& p : this->getPoints())
      p->setXSingle(this->rigidTransformation(p->initialX));
    // p->setXSingle(this->quaternion.Rv(p->initialX - this->ICOM) -
    // (p->initialX - this->ICOM) + trans + p->initialX);
    this->setGeometricPropertiesForce();
  };

  void calcPhysicalProperties() {
    this->force.fill(0.);
    this->inertia.fill(0.);
    this->mass = 0.;
    for (const auto& p : this->Points) {
      // for (auto i = 0; i < 3; ++i)
      // 	this->force[i] += p->force[i];
      std::get<0>(this->force) += std::get<0>(p->force);
      std::get<1>(this->force) += std::get<1>(p->force);
      std::get<2>(this->force) += std::get<2>(p->force);
      this->mass += p->mass; // total mass
    }
    std::get<0>(this->inertia) = std::get<1>(this->inertia) = std::get<2>(this->inertia) = this->mass; // total mass

    // this->center_of_mass = {0, 0, 0};
    // //! -------------------------------------- */
    // for (const auto &p : this->Points)
    // 	this->center_of_mass += p->X;
    // this->center_of_mass /= this->mass;
    // //! -------------------------------------- */
    Tddd tmp = {0, 0, 0};
    for (const auto& p : this->Points) {
      tmp = p->X - this->center_of_mass;
      std::get<3>(this->force) += std::get<0>(p->force) * std::get<0>(tmp);
      std::get<4>(this->force) += std::get<1>(p->force) * std::get<1>(tmp);
      std::get<5>(this->force) += std::get<2>(p->force) * std::get<2>(tmp);
      //
      std::get<3>(this->inertia) += p->mass * std::pow(std::get<0>(tmp), 2.);
      std::get<4>(this->inertia) += p->mass * std::pow(std::get<1>(tmp), 2.);
      std::get<5>(this->inertia) += p->mass * std::pow(std::get<2>(tmp), 2.);
    }
  };

  //* ------------------- 物体全体にかかる力の計算方法２つ -------------------
  //*/

  // まずは，点にかかる力を計算してから使うこと．
  void sumForceOfPoints() {
    //%this->center_of_massが設定してある必要がある
    //%ネットワークの各点にかかっている力の総和を計算する
    this->force.fill(0.);
    for (const auto& p : this->Points) {
      auto F = p->force;
      auto r = (p->X - this->center_of_mass);
      auto [Fx, Fy, Fz, Fx_, Fy_, Fz_] = F;
      auto [Nx, Ny, Nz] = Cross(r, Tddd{Fx, Fy, Fz}); // モーメント
      std::get<0>(this->force) += Fx;
      std::get<1>(this->force) += Fy;
      std::get<2>(this->force) += Fz;
      std::get<3>(this->force) += Nx;
      std::get<4>(this->force) += Ny;
      std::get<5>(this->force) += Nz;
    }
  };
  void integrateForceOnFace() {
    /*(*checked by mathematica*)
    f0 = {f0x, f0y, f0z}
    f1 = {f1x, f1y, f1z}
    f2 = {f2x, f2y, f2z}
    x0 = {x0x, x0y, x0z}
    x1 = {x1x, x1y, x1z}
    x2 = {x2x, x2y, x2z}
    f[x0_, x1_] := Dot[{f0, f1, f2}, {x0, x1, 1 - x0 - x1}]
    X[z0_, z1_] := Dot[{x0, x1, x2}, {z0, z1, 1 - z0 - z1}]
    CForm@Flatten[
    Join[Simplify@
            Integrate[f[z0, z1], {\[Xi]0, 0, 1}, {z1, 0, 1 - z0}],
    Simplify@
            Integrate[
            Cross[X[z0, z1], f[z0, z1]], {z0, 0, 1}, {z1, 0, 1 - z0}]]]
    */
    this->force.fill(0.);
    for (const auto& f : this->Faces) {
      auto [p0, p1, p2] = f->getPoints();
      auto [f0x, f0y, f0z, f0x_, f0y_, f0z_] = p0->force;
      auto [f1x, f1y, f1z, f1x_, f1y_, f1z_] = p1->force;
      auto [f2x, f2y, f2z, f2x_, f2y_, f2z_] = p2->force;
      auto [x0x, x0y, x0z] = p0->X - this->center_of_mass;
      auto [x1x, x1y, x1z] = p1->X - this->center_of_mass;
      auto [x2x, x2y, x2z] = p2->X - this->center_of_mass;
      this->force += T6d{(f0x + f0y + f0z) / 6.,
                         (f1x + f1y + f1z) / 6.,
                         (f2x + f2y + f2z) / 6.,
                         (f2z * x1x + f2z * x1y + 2 * f2z * x1z + f2x * (2 * x1x + x1y + x1z) + f2y * (x1x + 2 * x1y + x1z) - 2 * f1x * x2x - f1y * x2x - f1z * x2x - f1x * x2y - 2 * f1y * x2y - f1z * x2y - f1x * x2z - f1y * x2z - 2 * f1z * x2z) / 24.,
                         (-(f2z * x0x) - f2z * x0y - 2 * f2z * x0z - f2x * (2 * x0x + x0y + x0z) - f2y * (x0x + 2 * x0y + x0z) + 2 * f0x * x2x + f0y * x2x + f0z * x2x + f0x * x2y + 2 * f0y * x2y + f0z * x2y + f0x * x2z + f0y * x2z + 2 * f0z * x2z) / 24.,
                         (f1z * x0x + f1z * x0y + 2 * f1z * x0z + f1x * (2 * x0x + x0y + x0z) + f1y * (x0x + 2 * x0y + x0z) - 2 * f0x * x1x - f0y * x1x - f0z * x1x - f0x * x1y - 2 * f0y * x1y - f0z * x1y - f0x * x1z - f0y * x1z - 2 * f0z * x1z) / 24.};
    }
  };
  //    void calcAccelFromForce() {
  //       //! 慣性を設定しておく必要がある．
  //       this->acceleration = this->force / this->inertia;
  //    };

  //@ ------------------------------------------------------ */

  // BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM
  /*Network_BEM_detail
  境界要素法用のNetworkには，以下が設定されている．

  * bool BASE;
  * bool CORNER;
  * bool Dirichlet;
  * bool Neumann;

  `networkPoint`も`this->getNetwork()->setB()`などとして，自身の境界条件判断できる．
  Network_BEM_detail*/
  /*Network_BEM_code*/
public:
  bool isFluid;
  bool isRigidBody;
  bool isSoftBody;
  bool isFloatingBody;
  bool isAbsorber;

  // Multi-part rigid body: component_id (on each face) indexes into this vector
  std::vector<std::string> component_names;

  std::function<Tddd(const Tddd&, double)> absorb_velocity = [](const Tddd&, double) -> Tddd { return {0., 0., 0.}; };
  std::function<Tddd(const Tddd&, double)> absorb_gradPhi_t = [](const Tddd&, double) -> Tddd { return {0., 0., 0.}; };
  std::function<double(const Tddd& X, const double t)> absorb_phi = [](const Tddd& X, const double t) -> double { return 0.; };
  std::function<double(const Tddd& X, const double t)> absorb_eta = [](const Tddd& X, const double t) -> double { return 0.; };
  std::function<double(double signed_distance)> absorb_gamma = [](double) -> double { return 0.; };
  WaterWaveTheory water_wave_theory;
  RandomWaterWaveTheory random_water_wave_theory;

  // Initial condition lambdas (set via "initial_condition" in Fluid JSON; nullptr = not set)
  std::function<double(const Tddd& X, const double t)> ic_phi = nullptr;
  std::function<double(const Tddd& X, const double t)> ic_eta = nullptr;

#ifdef BEM
public:
  double _current_time_;
  //
  std::array<bool, 6> isFixed = {false, false, false, false, false, false};
  JSON inputJSON;
  std::tuple<std::string, double> velocity_name_start;
  const double move_amplitude = 0.4;
  T6d velocityPredefined();
  Tddd translationPredefined();
  V_d nabla_phi;
  V_d meanX;
#endif
  /*Network_BEM_code*/
  // BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM

private:
  using V_netLp = std::vector<networkLine*>;
  std::string name = "no_name";
  std::string filename = "no_name";
  // 接続している点をLinesをとおして素早く参照できる
protected:
  std::unordered_set<networkTetra*> Tetras;
  std::unordered_set<networkPoint*> Points;
  std::unordered_set<networkFace*> Faces;
  std::vector<networkPoint*> Points_vector;
  std::vector<networkFace*> Faces_vector;
  V_netPp PointGarbage;
  V_netFp FaceGarbage;

  // 接続を持ち得ないnetworkPoint
  std::unordered_set<networkPoint*> IsolatedPoints;

  // ============================
  // Coordinate Scaling for FMM numerical stability
  // ============================
  double scale_factor_ = 1.0;
  bool is_scaled_ = false;
  // 元座標は各点の X_temporary に保存（反復順序に依存しない安全な方法）

public:
  // Compute characteristic length (3/4 of bounding box diagonal)
  // FMM uses rho = |r - r_center|, with rho_max = diagonal/2 (center to corner).
  // Using 3/4 of diagonal provides good numerical stability while keeping rho normalized.
  double computeCharacteristicLength() const {
    auto [xrange, yrange, zrange] = this->bounds;
    double Lx = std::abs(std::get<1>(xrange) - std::get<0>(xrange));
    double Ly = std::abs(std::get<1>(yrange) - std::get<0>(yrange));
    double Lz = std::abs(std::get<1>(zrange) - std::get<0>(zrange));
    return std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
  }

  static bool scalingDebugEnabled() {
    if (const char* env = std::getenv("BEM_SCALE_DEBUG"))
      return std::string(env) != "0";
    return false;
  }

  void debugPrintScalingState(const std::string& tag) const {
    if (!scalingDebugEnabled())
      return;
    std::ostringstream oss;
    oss << Magenta << "[scale] " << Cyan << tag << colorReset
        << " net=" << (this->name.empty() ? "<unnamed>" : this->name)
        << " is_scaled=" << is_scaled_
        << " scale_factor=" << scale_factor_
        << " charL=" << computeCharacteristicLength()
        << " points=" << Points.size()
        << " lines=" << Lines.size()
        << " faces=" << Faces.size()
        << " bounds=" << static_cast<const CoordinateBounds&>(*this);
    if (!Points.empty()) {
      const auto* p = *Points.begin();
      oss << " sampleX=" << p->X;
    }
    std::cout << oss.str() << std::endl;
  }

  // Apply coordinate scaling: store original coordinates in X_temporary, scale by 1/L
  void applyScaling(double L) {
    if (is_scaled_) {
      throw std::runtime_error("Network::applyScaling: already scaled");
    }
    if (L <= 0) {
      throw std::runtime_error("Network::applyScaling: L must be positive");
    }
    scale_factor_ = L;
    // Store original coordinates in X_temporary and scale
    for (auto& p : this->Points) {
      p->X_temporary = p->X;
      p->setXSingle(p->X / L);
    }
    // Scale edge midpoint positions (true quadratic geometry)
    for (auto& l : this->Lines) {
      l->setXSingle(l->X_mid / L);
    }
    is_scaled_ = true;
    this->setGeometricPropertiesForce();
  }

  // Remove coordinate scaling: restore original coordinates from X_temporary
  void removeScaling() {
    if (!is_scaled_) {
      return; // Nothing to do
    }
    // Restore original coordinates from X_temporary
    for (auto& p : this->Points) {
      p->setXSingle(p->X_temporary);
    }
    // Restore edge midpoint positions
    for (auto& l : this->Lines) {
      l->setXSingle(l->X_mid * scale_factor_);
    }
    scale_factor_ = 1.0;
    is_scaled_ = false;
    this->setGeometricPropertiesForce();
  }

  double getScaleFactor() const { return scale_factor_; }
  bool isScaled() const { return is_scaled_; }

public:
  std::unordered_set<networkLine*> Lines;
  bool MemberQ(networkPoint* const& p_IN) const;
  bool MemberQ(networkFace* const& f_IN) const;
  void add(netPp const p_IN);
  void add(netLp const l_IN);
  void add(netFp const f_IN);
  void add(netTp const t_IN);
  void add(const V_netPp& ps_IN);
  void add(const V_netFp& fs_IN);
  bool erase(networkFace* f);
  void erase(const V_netFp& fs);
  bool erase(networkPoint* p);
  bool erase(networkLine* l);

  // パッチ貼り付け: patchA 領域を patchB で置き換える
  // 外周点の対応は座標比較（coord_eps）で取る
  bool replacePatch(Network& patchB, double coord_eps = 1e-10);

  //% ------------------------- obj ------------------------ */
  //%
  // NetwrokObjは，Networkとしてコンストラクトするために，下のコンストラクタように修正した．*/
  // \label{Network::constructor}
  // Netwrokコンストラク
  Network(const std::string& filename = "file_name_is_not_given", const std::string& name_IN = "no_name");
  //% ------------------------------------------------------ */
  const std::string& getName() const { return this->name; };
  const std::string& getFileName() const { return this->filename; };
  void setName(const std::string nameIN) { this->name = nameIN; };
  // V_netPp linkXPoints(Network &water, Network &obj);
  /* ------------------------------------------------------ */
  void displayStates() {
    std::cout << "/* -------------------------------------------------------------------------- */" << colorReset << std::endl;
    const int size = 15;
    V_i num(size, 0), pointsLines(size, 0);
    // 点の持つ線の数をカウント
    for (const auto& p : this->Points)
      for (auto i = 0; i < size; i++)
        if (p->Lines.size() == (size_t)i)
          pointsLines[i]++;
    V_i connection(3);
    for (const auto& l : this->getLines()) {
      if ((l->getFaces()).size() == 0)
        connection[0]++;
      else if ((l->getFaces()).size() == 1)
        connection[1]++;
      else if ((l->getFaces()).size() == 2)
        connection[2]++;
    }
    std::cout << Magenta << std::setw(25) << "this->getName(): " << this->getName() << ", " << this << std::endl;
    std::cout << Magenta << std::setw(25) << "   Lines.size(): " << getLines().size() << std::endl;
    std::cout << Magenta << std::setw(25) << "  Points.size(): " << Points.size() << std::endl;
    std::cout << Magenta << std::setw(25) << "   Faces.size(): " << Faces.size() << std::endl;
    std::cout << Magenta << std::setw(25) << "   this->bounds: " << this->bounds << std::endl;
    std::cout << Green << std::setw(25) << "                  " << GridVector(Subdivide(0., (double)size - 1., size - 1), 6) << std::endl;
    std::cout << Magenta << std::setw(25) << "Lines of points : " << GridVector(pointsLines, 6) << std::endl;
    std::cout << Magenta << std::setw(25) << " Connection of faces : " << GridVector(connection, 6);
    if (connection[0] == 0 && connection[1] == 0 && connection[2] != 0)
      std::cout << Blue << " 全ての線が2つの面と接続しているため，格子は閉じた面を形成している" << colorReset;
    std::cout << colorReset << std::endl;
    std::cout << "/* -------------------------------------------------------------------------- */" << colorReset << std::endl;
  };
  double rotation_offset;
  V_d translation_offset;
  double rotation;
  V_d translation;
  void rotate(const double theta, const Tddd& axis) {
    for (auto& p : this->getPoints())
      p->X = Dot(RotationMatrix(theta, axis), p->X);
    this->setGeometricPropertiesForce();
  };
  /* ------------------------------------------------------ */
  /*             ネットワークの姿勢に関する変数とメソッド          */
  /* ------------------------------------------------------ */
  void rotate(const Quaternion& Q, const Tddd& barycenter = {0., 0., 0.});
  void scale(const double ratio, const Tddd& center = {0, 0, 0});
  void scale(const Tddd& ratio, const Tddd& center = {0, 0, 0});
  void translate(const Tddd& translation);
  void translateFromInitialX(const Tddd& translation);
  void resetInitialX();
  void applyTransformations(const JSON& json); // Implemented in a .cpp file not provided
  const std::unordered_set<networkPoint*>& getPoints() const;
  template <typename F>
  const std::vector<networkPoint*> getPoints(F&& cond) const {
    std::vector<networkPoint*> ret;
    ret.reserve(this->Points.size());
    for (const auto& p : this->Points)
      if (cond(p))
        ret.emplace_back(p);
    return ret;
  }
  std::vector<networkPoint*> getBoundaryPoints() const;
  /* ------------------------------ 境界条件に関する設定関数 ------------------------------ */

  std::unordered_set<networkPoint*> getPointsCORNER() const;
  std::unordered_set<networkPoint*> getPointsNeumann() const;
  std::unordered_set<networkPoint*> getPointsDirichlet() const;
  void setMinDepthFromCORNER();

  /* --------------------------------------------------------------------------*/

  const std::unordered_set<networkFace*>& getFaces() const;
  const std::vector<networkFace*>& getFacesVector() const;
  const std::vector<networkPoint*>& getPointsVector() const;
  // Update Faces_vector based on the contents of Faces
  void setFacesVector();
  // Update Points_vector based on the contents of Points
  void setPointsVector();
  const std::unordered_set<networkTetra*>& getTetras() const;
  const std::vector<networkTetra*> getIsolatedTetras() const;
  std::unordered_set<networkPoint*> getParametricPoints() const;
  Tddd getMeanX() const;
  void getLocations(VV_d& ret) const;
  std::unordered_set<networkLine*> getLines() const;
  std::vector<networkFace*> getActiveFaces() const {
    std::vector<networkFace*> ret;
    ret.reserve(this->Faces.size());
    for (const auto& f : this->Faces)
      if (f && f->is_active_face)
        ret.emplace_back(f);
    if (!ret.empty())
      return ret;
    return this->getBoundaryFaces();
  }
  std::vector<networkFace*> getBoundaryFaces() const;
  std::vector<networkLine*> getBoundaryLines() const;
  std::vector<networkFace*> getComputationalBoundaryFaces() const;
  std::vector<networkLine*> getComputationalBoundaryLines() const;
  std::vector<networkPoint*> getComputationalBoundaryPoints() const;
  std::vector<networkFace*> getInteriorFaces() const;
  std::vector<networkLine*> getInteriorLines() const;
  std::vector<networkPoint*> getInteriorPoints() const;
  template <typename F>
  std::vector<networkLine*> getLines(F&& cond) const {
    std::vector<networkLine*> ret;
    ret.reserve(this->Lines.size());
    for (const auto& l : this->Lines)
      if (cond(l))
        ret.emplace_back(l);
    return ret;
  };

  std::unordered_set<networkLine*> getLinesUO() const;
  V_netLp getLinesIntxn() const;
  //======================================
  /*setXPoints_detail
  `setXPoints`は，thisとtargetの`networkFace`と`networkLine`の干渉を検知して，干渉点として`networkPoint`を新たに作成する．

  -
  干渉点`networkPoint`は，干渉した`networkFace`と`networkLine`の`XPoints`に`push`される．
  -
  干渉点として作成される`networkPoint`には，`xline`と`xface`が与えられる．これらは，デフォルトで`xline=nullptr`，`xface=nullptr`となっている．

  setXPoints_detail*/

  // void setXPoints(Network &target, Network *interactionNet);

  /**
  p->setFaces()もf->setFaces()も,自身のp->Lines,f->Linesを元に，p->Faces,f->Facesを決定し保存する．
  p->Linesとf->Linesが正しく設定してあるかチェックする．
  特に，flipやdivideの後には，p->Lines,f->Linesが正しく設定されていない可能性があるので，要注意．
  */

  void setGeometricPropertiesForce();
  void computePrincipalCurvatures(bool force_all = true);
  std::size_t topologySignature() const;
  std::size_t geometrySignature(double inv_eps = 1e12) const;

  // 局所パッチのコピー: 辺または点の近傍を別の Network にコピーする
  networkLine* copyLocalPatch(Network& patch, networkLine* l, int ring_depth = 1) const;
  networkPoint* copyLocalPatch(Network& patch, networkPoint* p, int ring_depth = 1) const;

private:
  void assignPointFaceIndices();
  void setGeometricPropertiesImpl();
  std::size_t computeTopologySignature() const;
  std::size_t computeGeometrySignature(double inv_eps = 1e12) const;
  mutable std::size_t topology_signature_cached = 0;
  mutable std::size_t geometry_signature_cached = 0;
  mutable bool geometry_signature_valid = false;

  /**
  pointのFacesとLinesの関係は整合性があるか？faceのFacesとLinesの関係は整合性があるか？をチェック．
  setGeometricPropertiesForce()を実行していれば，この整合性は保たれるはずではある．
  */

  std::array<bool, 4> validateConectivity();

  /* --------------------------------------------------------------------------*/

public:
  V_netPp setPoints(const std::vector<Tddd>& v_IN);
  V_netPp setPoints(const VV_d& v_IN);
  /* --------------------------------------------------------------------------
   */
public:
  void setFaces(const std::vector<T3Tddd>& v_IN, const double resolution = 1E-10);
  void setFaces(const VV_i& f_v, const V_netPp& points);
  void setLines(const VV_i& l_v, const V_netPp& points);
  void mergeObjFile(const std::string& filename, int component_id);
  /* -------------------------------------------------------------------------- */

public:
  void DeleteParticles();
  void DeleteIsolatedTetras();
  void DeleteInteriorTetras();

  /* ------------------------------------------------------ */
public:
  std::vector<MooringLine*> mooringLines;

#ifdef tetgenH

  void genTetra(const tetgenio& out) {
    std::cout << Red << "genTetra" << colorReset << std::endl;
    int index;
    Tddd X;
    std::array<networkPoint*, 4> points;
    // tetgenが出力した全四面体をループ
    for (int i = 0; i < out.numberoftetrahedra; i++) {
      points = {nullptr, nullptr, nullptr, nullptr};
      // 四面体を構成する4つの頂点をループ
      for (int j = 0; j < 4; j++) {
        // 頂点のインデックスと座標を取得
        index = out.tetrahedronlist[i * 4 + j];
        X = {out.pointlist[index * 3], out.pointlist[index * 3 + 1], out.pointlist[index * 3 + 2]};
        // 座標を元に、対応するnetworkPointオブジェクトを空間ハッシュ(Bucket)から高速に検索
        for (auto& p : this->BucketPoints.getData(X))
          if (Norm(p->X - X) < 1E-12) {
            points[j] = p; // 見つかったポインタを格納
            break;
          }
      }
      if (points[0] && points[1] && points[2] && points[3])
        ::genTetra(this, points);
    }
    std::cout << Green << "Total Tetra : " << this->Tetras.size() << colorReset << std::endl;
  }

  /* --------------------------------------------------------------------------*/

  bool tetrahedralize_(const std::string command, tetgenio& in) {

    printf("Loaded %d vertices and %d faces from Network.\n", in.numberofpoints, in.numberoffacets);

    // 2. tetgenの動作を制御するオプションを設定
    tetgenbehavior b;
    b.parse_commandline(const_cast<char*>(command.c_str()));

    // 3. tetgenライブラリを実行し、結果を出力データ(out)に格納
    tetgenio out;
    ::tetrahedralize(&b, &in, &out);
    printf("Generated %d tetrahedra.\n", out.numberoftetrahedra);

    /* --------------------------------------------------------------------------*/

    // 4. tetgenが追加した新しい頂点（シュタイナー点）をNetworkに追加．既存の点と重複しないように、空間ハッシュ(Bucket)を使って高速にチェックする
    this->makeBucketPoints(this->getScale() / 50.);
    std::vector<networkPoint*> addedPoints;
    int count = 0;
    for (int i = 0; i < out.numberofpoints; i++) {
      std::array<double, 3> X = {out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2]};
      // 同じ座標の点がまだNetworkに存在しないかチェック
      if (std::ranges::none_of(this->BucketPoints.getData(X), [&](auto& p) { return Norm(p->X - X) < 1E-15; })) {
        auto p = new networkPoint(this, X);
        this->BucketPoints.add(X, p);
        addedPoints.push_back(p);
        count++;
      }
    }
    std::cout << "added points : " << count << std::endl;
    this->setGeometricPropertiesForce();
    genTetra(out);
    this->setGeometricPropertiesForce();

    // 5. 不要な面を削除（両側に四面体が存在しない面）
    auto faces = ToVector(this->getFaces());
    for (auto& f : faces)
      if (f->Tetras[0] == nullptr && f->Tetras[1] == nullptr)
        delete f;

    this->setGeometricPropertiesForce();
    this->makeBuckets();
    this->checkConnectivity();

    for (auto& p : addedPoints) {
      auto& phiphin = p->phiphin;
      auto& phiphin_t = p->phiphin_t;
      phiphin = Tdd{0, 0};
      phiphin_t = Tdd{0, 0};
      for (auto& q : p->getNeighborPointsOnSurfaces()) {
        phiphin += q->phiphin;
        phiphin_t += q->phiphin_t;
      }
      phiphin /= p->getNeighborPointsOnSurfaces().size();
      phiphin_t /= p->getNeighborPointsOnSurfaces().size();
    }

    return true;
  }

  bool tetrahedralize(const std::string command, const auto& surfaces, const bool delete_tetra_first = true) {
    //! delete all existing tetrahedra
    if (delete_tetra_first) {
      auto existingTetras = ToVector(this->getTetras());
      for (auto& t : existingTetras)
        delete t;
    }
    tetgenio in = generate_tetgenio_input(surfaces);
    return tetrahedralize_(command, in);
  }

  /* --------------------------------------------------------------------------*/

  bool tetrahedralize(const std::string command = "p/0.01Dq1.YSM") { return this->tetrahedralize(command, this->getBoundaryFaces()); }

  /* --------------------------------------------------------------------------*/
  /*                     点群から四面体を生成するメソッド                      */
  /* --------------------------------------------------------------------------*/

  bool tetrahedralizeFromPoints(const std::string command = "") {

    // 既存の四面体を削除
    auto existingTetras = ToVector(this->getTetras());
    for (auto& t : existingTetras)
      delete t;

    // 1. 現在の点群からtetgen用の入力データ(in)を生成
    // unordered_setからvectorへ変換して渡す
    auto pointsVec = ToVector(this->getPoints());
    tetgenio in = generate_tetgenio_input_from_points(pointsVec);
    return tetrahedralize_(command, in);
  }

  /* --------------------------------------------------------------------------*/
  /*           閉じていない面と点群から四面体を生成するメソッド (PLC)          */
  /* --------------------------------------------------------------------------*/
  bool tetrahedralizePLC(const std::string command = "") {

    // 既存の四面体を削除
    auto existingTetras = ToVector(this->getTetras());
    for (auto& t : existingTetras)
      delete t;

    // 1. 全点と全面を取得
    // Network内の全ての点と面を使用します
    auto pointsVec = ToVector(this->getPoints());
    auto facesVec = ToVector(this->getFaces());

    // 2. tetgenio入力生成 (Mixed)
    tetgenio in = generate_tetgenio_input_mixed(pointsVec, facesVec);
    return tetrahedralize_(command, in);
  }

#endif

  void checkConnectivity() const;

public:
  bool erase_element(networkPoint* const p);
  bool erase_element(networkLine* const l);
  bool erase_element(networkFace* const f);
  bool erase_element(networkTetra* const t);

  size_t improveTetrahedraDelaunay(const double eps = 1E-12);
  struct Flip23Candidate {
    networkFace* f = nullptr;
    networkTetra *t0 = nullptr, *t1 = nullptr;
    networkPoint *a = nullptr, *b = nullptr, *c = nullptr, *p = nullptr, *q = nullptr;
    double smaller_is_better{std::numeric_limits<double>::infinity()};
    double s_before{0.};
    double s_after{0.};
  };

  struct Flip32Candidate {
    networkLine* edge = nullptr;
    std::array<networkTetra*, 3> tetras{};
    networkPoint *p = nullptr, *q = nullptr, *a = nullptr, *b = nullptr, *c = nullptr;
    double smaller_is_better{std::numeric_limits<double>::infinity()};
    double s_before{0.};
    double s_after{0.};
  };
  Flip32Candidate improveTetrahedraDelaunay32(double eps = 1E-12, std::function<double(double, double)> quality_metric = nullptr);
  void improve(Flip32Candidate& best32);

  /* --------------------------------------------------------------------------
   */

  Flip23Candidate improveTetrahedraDelaunay23(double eps = 1E-12, std::function<double(const double, const double)> quality_metric = nullptr);

  void improve(Flip23Candidate& best23);
};

//% ==========================================================================
//% ==========================================================================

#include "my_vtk.hpp"
// #include "networkPoint.cpp"

// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM
#ifdef BEM
inline V_Netp takeNetworks(const V_netPp& ps) {
  V_Netp ret({});
  for (auto& p : ps)
    network::add(ret, p->getNetwork());
  return ret;
};

#endif
// BEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEMBEM

//@ ------------------------------------------------------ */
//@                          外部関数                       */
//@ ------------------------------------------------------ */

// !
// 注意することは，deleteは接続関係を全て解除することを伴い，再構築が必要になること．
// !
// flip等の操作では，できる限りdeleteせずにreplaceなどを多用し，最後に要らなくなったものをdeleteするようにする

inline void networkFace::Delete() {
  try {
    // 1. this->Linesからthisを削除する
    for (auto& l : this->Lines)
      if (l)
        l->erase(this);

    // 2. this->Pointsからthisを削除する
    for (auto& p : this->Points)
      if (p)
        p->erase(this);

    // 3. this->Tetrasからthisを削除する. 四面体を全て削除する
    auto [t0, t1] = this->Tetras;
    if (t0)
      delete t0; // 注意 ~networkTetraは内部でthis->Tetrasを改変する
    if (t1)
      delete t1; // 注意 ~networkTetraは内部でthis->Tetrasを改変する

    // 4. Networkからthisを削除する
    if (!this->network->erase(this))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "is not a member of " + this->network->getName());
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  }
};

/* -------------------------------------------------------------------------- */

/*!
pointのDelete()と同時に，pointが接続している辺のDelete()を必ず伴うべきだ．
*/

inline void networkPoint::Delete() {
  try {
    if (!this->network->erase(this)) {
      std::stringstream ss;
      ss << "この点は，networkに含まれていない" << std::endl;
      ss << "network name : " << this->network->getName() << std::endl;
      if (this->network->MemberQ(this))
        ss << this->network->getName() << "に含まれている" << std::endl;
      else
        ss << this->network->getName() << "にも含まれていない" << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    } else {
      // Delete adjacent faces first to clear line/point references safely.
      auto faces = this->getFaces();
      for (auto* f : faces)
        if (f)
          delete f;

      // Then delete all incident lines (faces should already be detached).
      auto lines = this->getLines();
      for (auto* l : lines)
        if (l)
          delete l;
    }
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  }
};

/* -------------------------------------------------------------------------- */

inline void networkLine::Delete() {
  if (this->Point_A) {
    this->Point_A->erase(this);
    this->Point_A = nullptr;
  }
  if (this->Point_B) {
    this->Point_B->erase(this);
    this->Point_B = nullptr;
  }

  if (!this->network->erase(this))
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "is not a member of " + this->network->getName());

  bool found = false;
  do {
    found = false;
    auto copy_Faces = this->Faces;
    for (auto f : copy_Faces)
      if (f) {
        delete f;
        found = true;
      }
  } while (found);
};

/* -------------------------------------------------------------------------- */

///////////////////////////////////////////////////////////////
inline networkPoint::~networkPoint() { this->Delete(); };
inline networkFace::~networkFace() { this->Delete(); };
inline networkLine::~networkLine() { this->Delete(); };
///////////////////////////////////////////////////////////////
inline Network::~Network() {
  if (octreeOfFaces)
    delete this->octreeOfFaces;
  //
  try {
    // std::cout << Red << this->getName() << colorReset << std::endl;
    // std::cout << "destroying Tetras" << std::endl;
    {
      auto tmp = this->Tetras;
      for (const auto& t : tmp)
        delete t;
    }
    // std::cout << "destroying Faces" << std::endl;
    {
      auto tmp = this->Faces;
      for (const auto& f : tmp)
        delete f;
    }
    // std::cout << "destroying Lines" << std::endl;
    {
      auto tmp = this->Lines;
      for (const auto& l : tmp)
        delete l;
    }
    // std::cout << "destroying Points" << std::endl;
    {
      auto tmp = this->Points;
      for (const auto& p : tmp)
        delete p;
    }
    if (this->Points.size() == 0 && this->Faces.size() == 0 && this->Lines.size() == 0 && this->Tetras.size() == 0) {
      // std::cout << this->getName() << " is deleted successfully" <<
      // std::endl;
    } else {
      std::stringstream message;
      message << "Points.size() = " << this->Points.size() << std::endl;
      message << "Faces.size() = " << this->Faces.size() << std::endl;
      message << "Lines.size() = " << this->Lines.size() << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message.str());
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << colorReset << std::endl;
  };
};

/*XNetwork_code*/

/*networkFace::divide_detail
divideの操作を整理して理解するために，networkFaceにdivideの一部操作を任せている．
ここでnetworkFaceは，与えられた線と自身がすでに所持していた線を，適切に接続し直す．
点と線の接続はnetworkFaceに行わせず，netowrkLineに行わせる．
networkFace::divide_detail*/
//
/*networkFace::divide_code*/
inline netFp networkFace::Split(netLp DivL /*this*/, netLp newDivL, netLp newMidL, int type) {
  if (!MemberQ_(this->Lines, DivL))
    throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "divL is not found in this->Lines"));
  auto oldF = this;
  auto newF = new networkFace(this);
  auto fL = oldF->getLineFront(DivL);
  auto bL = oldF->getLineBack(DivL);
  //
  //          /     \                     /     \
  //         /       \                   /       \
  //        /         \                 /         \
  //       /    / \    \               /    / \    \
  // bL<--/--- /   \----\-->fL   bL<--/----/   \----\-->fL
  //   --/--->/oldF \<---\--         /    /newF \    \
  //    /     --|A---     \         /     --|----     \
  //   /        V|         \       /        V          \
  //    -------DivL-------         --------DivL---------
  //
  if (oldF->replace(fL,
                    newMidL)) { // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
    std::stringstream ss;
    // ss << "oldF->getLines() = " << oldF->getLines();
    ss << ", oldF = " << oldF << ", fL = " << fL << ", newMidL = " << newMidL;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
  }

  // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
  if (newF->replace(bL, newMidL)) {
    std::stringstream ss;
    ss << "newF = " << newF << ", bL = " << bL << ", newMidL = " << newMidL;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
  }
  // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
  if (fL->replace(oldF, newF))
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

  //          /     \                     /     \
  //         /       \                   /       \
  //        /         \                 /         \
  //       /    / \    \               /    / \    \
  // bL<--/----/   \<---\-- newMidL<--/--  /   \----\-->fL
  //   --/--->/oldF \----\->         /    /newF \<---\---
  //    /     --|A---     \         /     --|----     \
  //   /        V|         \       /        V          \
  //    -------DivL-------         --------DivL---------
  //
  newMidL->set(oldF, newF);

  //          /     \                     /     \
  //         /       \                   /       \
  //        /         \                 /         \
  //       /    / \    \               /    / \    \
  // bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
  //  ---/--->/oldF \<---\--       --/--->/newF \<---\---
  //    /     --|A---     \         /     --|----     \
  //   /        V|         \       /        V          \
  //    -------DivL-------         --------DivL---------
  if (type == 1) {
    // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
    if (newF->replace(DivL, newDivL))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    newDivL->set(newF);
    //          /     \                     /     \
    //         /       \                   /       \
    //        /         \                 /         \
    //       /    / \    \               /    / \    \
    // bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
    //  ---/--->/oldF \<---\--       --/--->/newF \<---\---
    //    /     --|A---     \         /     --|A---     \
    //   /        V|         \       /        V|         \
    //    -------DivL/*this*/---     -------newDivL------
  } else if (type == 2) {
    // oldFのfLをnewMidLに繋ぎ直し，fLはnewFと繋げる
    if (oldF->Replace(DivL, newDivL))
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    DivL->add(newF);
    //          /     \                     /     \
    //         /       \                   /       \
    //        /         \                 /         \
    //       /    / \    \               /    / \    \
    // bL<--/----/   \----\-->newMidL<--/----/   \----\-->fL
    //  ---/--->/oldF \<---\--       --/--->/newF \<---\---
    //    /     --|A---     \         /     --|A---     \
    //   /        V|         \       /        V|         \
    //    ------newDivL------         ----DivL/*this*/----
  }
  return newF;
};

inline std::unordered_map<netP*, Tdd>& operator+=(std::unordered_map<netP*, Tdd>& ii_dd, const std::unordered_map<netP*, Tdd>& jj_dd) {
  std::unordered_map<netP*, Tdd>::iterator it;
  for (auto& [jj, dd] : jj_dd)
    if ((it = ii_dd.find(jj)) != ii_dd.end())
      it->second += dd;
    else
      ii_dd[jj] = dd;
  return ii_dd;
};

inline std::map<netP*, Tdd>& operator+=(std::map<netP*, Tdd>& ii_dd, const std::map<netP*, Tdd>& jj_dd) {
  std::map<netP*, Tdd>::iterator it;
  for (auto& [jj, dd] : jj_dd)
    if ((it = ii_dd.find(jj)) != ii_dd.end())
      it->second += dd;
    else
      ii_dd[jj] = dd;
  return ii_dd;
};

inline std::map<netP*, Tdd>& operator+=(std::map<netP*, Tdd>& ii_dd, const std::tuple<netP*, Tdd>& jj_dd) {
  std::map<netP*, Tdd>::iterator it;
  if ((it = ii_dd.find(std::get<0>(jj_dd))) != ii_dd.end())
    it->second += std::get<1>(jj_dd);
  else
    ii_dd[std::get<0>(jj_dd)] = std::get<1>(jj_dd);
  return ii_dd;
};

/* -------------------------------------------------------------------------- */

inline T4T3P ToT4T3P(const T_4P& p0123) {
  auto [p0, p1, p2, p3] = p0123;
  return {{{p0, p1, p2}, {p0, p1, p3}, {p0, p2, p3}, {p1, p2, p3}}};
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

// struct PseudoQuadPatch {
//    /*this nodes forms a quadratic element like below
//    0
//    | \
//    |   \
//    |     \
//    3(p2)-l1-5(p1)
//    | \    |   \
//    |  l2  l0   \
//    |   \  |     \
//    1 ----4(p0 --- 2
//    */
//    const std::array<networkPoint *, 3> p_of_face;
//    const std::array<networkLine *, 3> l_of_face;
//    const std::array<bool, 3> ignore;
//    const std::array<networkPoint *, 6> points;
//    PseudoQuadPatch(const networkPoint *origin, const networkFace *f)
//        : p_of_face(f->getPoints(origin)),
//          l_of_face(f->getLinesTupleFrom(std::get<0>(p_of_face))),
//          ignore({!std::get<1>(l_of_face)->CORNER,
//                  !std::get<2>(l_of_face)->CORNER,
//                  !std::get<0>(l_of_face)->CORNER}),
//          points({f->getPointOpposite(std::get<1>(l_of_face)) /*0*/,
//                  f->getPointOpposite(std::get<2>(l_of_face)) /*1*/,
//                  f->getPointOpposite(std::get<0>(l_of_face)) /*2*/,
//                  std::get<2>(p_of_face) /*3*/,
//                  std::get<0>(p_of_face) /*4*/,
//                  std::get<1>(p_of_face) /*5*/}){};
// };

struct PseudoQuadPatch {
  networkFace* f;
  const std::tuple<networkPoint*, networkLine*, networkPoint*, networkLine*, networkPoint*, networkLine*> PLPLPL;
  std::array<bool, 3> available; // availability of the lines (face)
  std::array<networkFace*, 3> faces;
  std::array<networkPoint*, 6> points;
  std::function<bool(const networkLine*)> conditionChecker;
  std::array<std::tuple<networkPoint*, networkFace*>, 6> points_faces;

  /*       0
         /   \
        / f1  \
       /       \
   3(p2)--(l1)--5(p1)
     /  \       / \
    /  (l2) fc l0  \
   / f2   \   /  f0 \
  1------ 4(p0)------2

  points_facesは，以下のように，点と点が所属する面を保持する（多重節点を使う際に役立つ）．

  {{point 0, face 1},
  {point 1, face 2},
  {point 2, face 0},
  {point 3, face center},
  {point 4, face center},
  {point 5, face center}}

  与えられた点または線が，面のPLPLPLのメンバーなら，
  その点または線を先頭にするように，PLPLPLを設定し，保持する．

  */

  // PseudoQuadPatch(const networkPoint *origin_point, const networkFace *f) : f(f),
  // PLPLPL(f->getPointsAndLinesFromNearest(origin_point)) { initialize(); };
  // PseudoQuadPatch(const networkLine *origin_line, const networkFace *f) : f(f),
  // PLPLPL(f->getPointsAndLines(origin_line)) { initialize(); };

  PseudoQuadPatch(networkFace* f, const networkPoint* origin_point = nullptr, std::function<bool(const networkLine*)> checker = [](const networkLine* line) { return !line->CORNER && line->getBoundaryFaces().size() == 2; }) : f(f), PLPLPL(f->getPointsAndLinesFromNearest(origin_point)), conditionChecker(checker) { initialize(); };

  PseudoQuadPatch(networkFace* f, const networkLine* origin_line = nullptr, std::function<bool(const networkLine*)> checker = [](const networkLine* line) { return !line->CORNER && line->getBoundaryFaces().size() == 2; }) : f(f), PLPLPL(f->getPointsAndLines(origin_line)), conditionChecker(checker) { initialize(); };

  void initialize() {
    auto [p0, l0, p1, l1, p2, l2] = this->PLPLPL;

    this->available = {conditionChecker(l1), conditionChecker(l2), conditionChecker(l0)};

    auto fs0 = l0->getBoundaryFaces();
    auto fs1 = l1->getBoundaryFaces();
    auto fs2 = l2->getBoundaryFaces();

    // if (fs0.size() == 2 && fs1.size() == 2 && fs2.size() == 2)

    if (fs0.empty() || fs1.empty() || fs2.empty()) {
      //! この辺に対しては，定義できない．
      std::cout << "l0 :" << l0 << " l0->Faces:" << l0->getFaces() << std::endl;
      std::cout << "l1 :" << l1 << " l1->Faces:" << l1->getFaces() << std::endl;
      std::cout << "l2 :" << l2 << " l2->Faces:" << l2->getFaces() << std::endl;

      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "fs0.empty() || fs1.empty() || fs2.empty()");
    }

    networkFace *f0, *f1, *f2;
    {
      f0 = fs0.size() == 1 ? fs0[0] : (fs0[0] == f ? fs0[1] : fs0[0]);
      f1 = fs1.size() == 1 ? fs1[0] : (fs1[0] == f ? fs1[1] : fs1[0]);
      f2 = fs2.size() == 1 ? fs2[0] : (fs2[0] == f ? fs2[1] : fs2[0]);
    }
    // else throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
    // "fs0.size() != 2 || fs1.size() != 2 || fs2.size() != 2");

    this->points = {f1->getPointOpposite(l1) /*0*/, f2->getPointOpposite(l2) /*1*/, f0->getPointOpposite(l0) /*2*/, p2 /*3*/, p0 /*4*/, p1 /*5*/};

    this->faces = {f1 /*0*/, f2 /*1*/, f0};

    this->points_faces[0] = {this->points[0], f1};
    this->points_faces[1] = {this->points[1], f2};
    this->points_faces[2] = {this->points[2], f0};
    this->points_faces[3] = {this->points[3], f};
    this->points_faces[4] = {this->points[4], f};
    this->points_faces[5] = {this->points[5], f};
  };
};

/*
線形補間するとして，
pの重みがわかる必要がある．

線形三角補間の場合，割り振りは，{t0,t1,1-t0-t1}である．
2次三角補間の場合，割り振りは，{t0,t1,1-t0-t1,t0*t1,t0*(1-t0-t1),t1*(1-t0-t1)}である．


最後は，
{t0,t1,1-t0-t1}

{t0,t1,1-t0-t1,t0*t1,t0*(1-t0-t1),t1*(1-t0-t1)}

(t0)*
(t1)*
(1-t0-t1) ->

p2 = (Dot({-1, -1, 0, 2, 4, 4},{p0,p1,p2,p3,p4,p5}) + Dot({-1, -1, 0, 2, 4,
4},{p0,p1,p2,p3,p4,p5}))/16


今ある関数でうまくできる？

N = {-1, -1, 0, 2, 4, 4}/8

2次補間の重みを，さらに節点上の補間の重みを通して，節点上の補間の重みを求める．

t0 * p0の変数
t1 * p1の変数
(1-t0-t1) * p2の変数
[t0*t1] * (Dot(N,{p0,p1,p2,p3,p4,p5}) + Dot(N,{p0,p1,p2,p3,p4,p5}))/2
[t0*(1-t0-t1)] * (Dot(N,{p0,p1,p2,p3,p4,p5}) + Dot(N,{p0,p1,p2,p3,p4,p5}))/2
[t1*(1-t0-t1)] * (Dot(N,{p0,p1,p2,p3,p4,p5}) + Dot(N,{p0,p1,p2,p3,p4,p5}))/2

これをみると，
{t0, t1, 1-t0-t1, N * t0*t1, N * t0*(1-t0-t1), N * t1*(1-t0-t1)}
となっていることがわかる．N = {-1, -1, 0, 2, 4, 4}/8．

言い換えると，
t0は，p0の重み
t1は，p1の重み
1-t0-t1は，p2の重み
N *
t0*t1は，p3の重みであるが，p3は存在しないので，これはさらに分配され，f_l0とfの節点に分配される．
N *
t0*(1-t0-t1)は，p4の重みであるが，p4は存在しないので，これはさらに分配され，f_l1とfの節点に分配される．
N *
t1*(1-t0-t1)は，p5の重みであるが，p5は存在しないので，これはさらに分配され，f_l2とfの節点に分配される．

*/

inline bool useOppositeFace(const networkLine* line, const double alpha = M_PI / 3.) {
  auto faecs = line->getBoundaryFaces();
  if (faecs.size() != 2)
    return false;
  auto [p0, p1] = line->getPoints();
  auto f0 = faecs[0];
  auto f1 = faecs[1];

  double signed_angle = VectorAngle(f0->normal, f1->normal, p1->X - p0->X);
  return std::abs(signed_angle) < alpha;
  // auto [a0, a1, f0_oppositeX] = f0->getPoints(p0);
  // if (f0_oppositeX == p1)
  //    std::swap(a1, f0_oppositeX);
  // auto [b0, b1, f1_oppositeX] = f1->getPoints(p0);
  // if (f1_oppositeX == p1)
  //    std::swap(b1, f1_oppositeX);
  // double base_length = Norm(p0->X - p1->X);

  // std::array<double, 3> Nmid = {5. / 4., 3 / 4., 1 - 2 * 3 / 4.};
  // std::array<double, 3> Nopposite = {-1., 1., 1 - (-1. + 1.)};
  // Norm(f0_oppositeX->X - (-f1_oppositeX->X + p0->X + p1->X));

  // auto f0_ideal_X = (-f1_oppositeX->X + p0->X + p1->X);
  // auto f1_ideal_X = (-f0_oppositeX->X + p0->X + p1->X);

  // Norm(f0_oppositeX->X - RotationMatrix(signed_angle, f0_ideal_X - p0->X) +
  // p0->X); Norm(f1_oppositeX->X - RotationMatrix(signed_angle, f1_ideal_X -
  // p0->X) + p0->X);

  // return (Norm(f0_oppositeX->X - (-f1_oppositeX->X + p0->X + p1->X)) /
  // base_length < alpha) &&
  //        (Norm(f1_oppositeX->X - (-f0_oppositeX->X + p0->X + p1->X)) /
  //        base_length < alpha);
};

template <typename T>
struct is_std_array_of_6 : std::false_type {
};
template <typename U>
struct is_std_array_of_6<std::array<U, 6>> : std::true_type {
};
template <typename T>
inline constexpr bool is_std_array_of_6_v = is_std_array_of_6<T>::value;

struct DodecaPoints {

  networkFace* f; // given center face
  const std::function<bool(const networkLine*)> conditionChecker;
  const PseudoQuadPatch quadpoint;
  const std::tuple<networkPoint*, networkLine*, networkPoint*, networkLine*, networkPoint*, networkLine*> PLPLPL;
  const networkLine *l0, *l1, *l2;
  const std::array<bool, 3> available; // availability of the lines (face)
  networkFace *f0, *f1, *f2;

  const PseudoQuadPatch quadpoint_l0;
  const PseudoQuadPatch quadpoint_l1;
  const PseudoQuadPatch quadpoint_l2;

  const T3Tdd shape_map_center = {{{0., 0.5} /*quad 4 -> linear 0*/, {0.5, 0.} /*quad 5 -> linear 1*/, {0.5, 0.5} /*quad 3 -> linear 2*/}};
  //! これとN3の内積を新なパラメタとして利用すると，(t0,t1)=(1,0)で{0.,
  //! 0.5}に，(t0,t1)=(0,1)で{0.5, 0.}に，t0=t1=0で{0.5, 0.5}になる．
  const T3Tdd shape_map_l0_face = {{{0.5, 0.} /*0*/, {0., 0.5} /*1*/, {0., 0.} /*2*/}};
  const T3Tdd shape_map_l1_face = {{{0., 0.} /*0*/, {0.5, 0.} /*1*/, {0., 0.5} /*2*/}};
  const T3Tdd shape_map_l2_face = {{{0., 0.5} /*0*/, {0., 0.} /*1*/, {0.5, 0.} /*2*/}};

  const PseudoQuadPatch* getPseudoQuadPatch(const int i) {
    if (i == 0)
      return &this->quadpoint_l0;
    else if (i == 1)
      return &this->quadpoint_l1;
    else if (i == 2)
      return &this->quadpoint_l2;
    else
      return &this->quadpoint;
  };

  DodecaPoints(
      networkFace* f, const networkPoint* origin_point = nullptr, std::function<bool(const networkLine*)> checker = [](const networkLine* line) { return !line->CORNER && line->getBoundaryFaces().size() == 2; })
      : f(f), conditionChecker(checker), quadpoint(f, origin_point, checker), PLPLPL(quadpoint.PLPLPL), l0(std::get<1>(PLPLPL)), l1(std::get<3>(PLPLPL)), l2(std::get<5>(PLPLPL)), available(quadpoint.available), f0(l0->getBoundaryFaces().size() == 1 ? l0->getBoundaryFaces()[0] : (l0->getBoundaryFaces()[0] == f ? l0->getBoundaryFaces()[1] : l0->getBoundaryFaces()[0])),
        f1(l1->getBoundaryFaces().size() == 1 ? l1->getBoundaryFaces()[0] : (l1->getBoundaryFaces()[0] == f ? l1->getBoundaryFaces()[1] : l1->getBoundaryFaces()[0])), f2(l2->getBoundaryFaces().size() == 1 ? l2->getBoundaryFaces()[0] : (l2->getBoundaryFaces()[0] == f ? l2->getBoundaryFaces()[1] : l2->getBoundaryFaces()[0])), quadpoint_l0(f0, l0, checker), quadpoint_l1(f1, l1, checker), quadpoint_l2(f2, l2, checker) {};
  /*
  これをみると，
  {t0, t1, 1-t0-t1, N * t0*t1, N * t0*(1-t0-t1), N * t1*(1-t0-t1)}
  となっていることがわかる．N = {-1, -1, 0, 2, 4, 4}/8．
  言い換えると，
  t0は，p0の重み
  t1は，p1の重み
  1-t0-t1は，p2の重み
  N *
  t0*t1は，p3の重みであるが，p3は存在しないので，これはさらに分配され，f_l0とfの節点に分配される．
  N *
  t0*(1-t0-t1)は，p4の重みであるが，p4は存在しないので，これはさらに分配され，f_l1とfの節点に分配される．
  N *
  t1*(1-t0-t1)は，p5の重みであるが，p5は存在しないので，これはさらに分配され，f_l2とfの節点に分配される．
  */

  /* --------------------------------------------------------------------------
   */

  /*

  `findMininum`は，節点のデータから補間された関数上での最小値を求める関数である．
  なので，滑らかな関数として補間できない場合は，この関数は使えない．

  節点の空間座標を補間した表面上で，ある別の関数の最小値を探したい場合，
  対象となる'別の関数'の微分も用意しておく必要がある．

  */

  // #define USE_BROYDEN_METHOD_FOR_FIND_MINIMUM
  // #define USE_NEWTON_METHOD_FOR_FIND_MINIMUM
#define USE_SUBDIVISION_METHOD_FOR_FIND_MINIMUM

  template <std::size_t N>
  std::tuple<Tdd, std::array<double, N>> findMinimum(const std::function<std::array<double, N>(networkPoint*)>& conversion, const auto& minimizing_function, const std::array<double, 2> t0_range = {0., 1.}) {

    Tdd first_guess = {0.25, 0.25};
    std::array<double, N> first_point = this->interpolate(first_guess, conversion), tmp;

    // if (first_point > (tmp = this->interpolate(0., 1., conversion))) {
    //    first_guess = {0., 1.};
    //    first_point = tmp;
    // }
    // if (first_point > (tmp = this->interpolate(1., 0., conversion))) {
    //    first_guess = {1., 0.};
    //    first_point = tmp;
    // }
    // if (first_point > (tmp = this->interpolate(0., 0., conversion))) {
    //    first_guess = {0., 0.};
    //    first_point = tmp;
    // }

    std::array<double, N> F, dFdt0, dFdt1, dFdt0t0, dFdt0t1, dFdt1t1;

#ifdef USE_BROYDEN_METHOD_FOR_FIND_MINIMUM
    const Tdd dt0t1 = {0., 1E-10};
    auto grad_FF = [&](const double t0, const double t1) -> Tdd {
      F = this->interpolate(t0, t1, conversion);
      dFdt0 = this->D_interpolate<1, 0>(t0, t1, conversion);
      dFdt1 = this->D_interpolate<0, 1>(t0, t1, conversion);
      return Tdd{Dot(dFdt0, F), Dot(dFdt1, F)};
    };
    BroydenMethod<Tdd> BM(first_guess, first_guess + dt0t1);
    for (auto i = 0; i < 50; i++) {
      BM.updateGoodBroyden(grad_FF(BM.X[0], BM.X[1]), grad_FF(BM.X[0] - BM.dX[0], BM.X[1] - BM.dX[1]));
      BM.X[0] = std::clamp(BM.X[0], t0_range[0], t0_range[1]);
      BM.X[1] = std::clamp(BM.X[1], 0., 1. - BM.X[0]);
      if (Norm(BM.dX) < 1e-13 && i > 3) {
        // std::cout << "iteration : " << i << "\n";
        break;
      }
    }
    return {BM.X, this->interpolate(BM.X[0], BM.X[1], conversion)};
#elif defined(USE_NEWTON_METHOD_FOR_FIND_MINIMUM)

    NewtonRaphson<Tdd> NR(first_guess);
    for (auto i = 0; i < 50; i++) {
      auto [t0, t1] = NR.X;

      F = this->interpolate(t0, t1, conversion);
      dFdt0 = this->D_interpolate<1, 0>(t0, t1, conversion);
      dFdt1 = this->D_interpolate<0, 1>(t0, t1, conversion);
      Tdd grad = {Dot(dFdt0, F), Dot(dFdt1, F)};

      dFdt0t1 = this->D_interpolate<1, 1>(t0, t1, conversion);
      dFdt0t0 = this->D_interpolate<2, 0>(t0, t1, conversion);
      dFdt1t1 = this->D_interpolate<0, 2>(t0, t1, conversion);
      T2Tdd Hessian = {Tdd{Dot(dFdt0t0, F) + Dot(dFdt0, dFdt0), Dot(dFdt0t1, F) + Dot(dFdt1, dFdt0)}, Tdd{Dot(dFdt0t1, F) + Dot(dFdt0, dFdt1), Dot(dFdt1t1, F) + Dot(dFdt1, dFdt1)}};

      NR.update(grad, Hessian);

      NR.X[0] = std::clamp(NR.X[0], t0_range[0], t0_range[1]);
      NR.X[1] = std::clamp(NR.X[1], 0., 1. - NR.X[0]);
      if (Norm(NR.dX) < 1e-13 && i > 2)
        // std::cout << "iteration : " << i << "\n";
        break;
    }
    return {NR.X, this->interpolate(NR.X[0], NR.X[1], conversion)};

    /* --------------------------------------------------------------------------
     */
#elif defined(USE_SUBDIVISION_METHOD_FOR_FIND_MINIMUM)
    double min = 1E10;
    Tdd t0t1_optimum = {0.25, 0.25};
    double t0_min = 0;
    double t0_max = 1;
    double t1_min = 0;
    double t1_max = 1;
    double div = 3, TMP;
    auto samples = convertions(conversion);
    double t0, t1;
    int count = 0;
    double eps = 1E-10;
    for (auto i = 0; i < 10; i++) {
      double dt0 = (t0_max - t0_min) / div;
      double dt1 = (t1_max - t1_min) / div;
      for (int j = 0; j <= div; j++) {
        t0 = t0_min + j * dt0;
        for (int k = 0; k <= div; k++) { // Adjust based on your constraints
          t1 = t1_min + k * dt1;
          if (t0 + t1 > 1)
            break;
          TMP = minimizing_function(this->Dot4(N6_new(t0, t1), samples));
          if (min > TMP) {
            min = TMP;
            t0t1_optimum = {t0, t1};
          }
        }
      }

      t0_min = std::clamp(t0t1_optimum[0] - 0.5 * dt0, t0_range[0], t0_range[1] - eps);
      t0_max = std::clamp(t0t1_optimum[0] + 0.5 * dt0, t0_min, t0_range[1]);
      t1_min = std::clamp(t0t1_optimum[1] - 0.5 * dt1, 0., 1. - eps);
      t1_max = std::clamp(t0t1_optimum[1] + 0.5 * dt1, t1_min, 1.);
    }

    // std::cout << "count : " << count << "\n";
    return {t0t1_optimum, this->interpolate(t0t1_optimum[0], t0t1_optimum[1], conversion)};

#endif
  };

  /* --------------------------------------------------------------------------
   */

  /*
  保存された，this->quadpoint.pointsを使った擬2次補間
  */

  std::array<T6d, 4> N6_new(const double t0, const double t1) const {
    const auto [N0, N1, N2, N3, N4, N5] = TriShape<6>(t0, t1);
    T6d Nc = {0., 0., 0., N2, N0, N1};
    T6d Nl0 = {0., 0., 0., 0., 0., 0.};
    T6d Nl1 = {0., 0., 0., 0., 0., 0.};
    T6d Nl2 = {0., 0., 0., 0., 0., 0.};
    if (this->conditionChecker(this->l0)) {
      Nc += 0.5 * N3 * TriShape<6>(0.25, 0.25, this->quadpoint.available);
      Nl0 = 0.5 * N3 * TriShape<6>(0.25, 0.25, this->quadpoint_l0.available);
    } else
      Nc += T6d{0., 0., 0., 0., 0.5 * N3, 0.5 * N3};

    if (this->conditionChecker(this->l1)) {
      Nc += 0.5 * N4 * TriShape<6>(0.5, 0.25, this->quadpoint.available);
      Nl1 = 0.5 * N4 * TriShape<6>(0.25, 0.25, this->quadpoint_l1.available);
    } else
      Nc += T6d{0., 0., 0., 0.5 * N4, 0., 0.5 * N4};

    if (this->conditionChecker(this->l2)) {
      Nc += 0.5 * N5 * TriShape<6>(0.25, 0.5, this->quadpoint.available);
      Nl2 = 0.5 * N5 * TriShape<6>(0.25, 0.25, this->quadpoint_l2.available);
    } else
      Nc += T6d{0., 0., 0., 0.5 * N5, 0.5 * N5, 0.};

    return {Nc, Nl0, Nl1, Nl2};
  };

  // CORNER edge quadratic offset at parametric coordinates
  // Returns {0,0,0} when all offsets are zero (inactive)
  Tddd corner_offset(const double t0, const double t1) const {
    const auto [N0, N1, N2, N3, N4, N5] = TriShape<6>(t0, t1);
    Tddd offset = {0., 0., 0.};
    if (!this->conditionChecker(this->l0) && this->l0->CORNER)
      offset += N3 * this->l0->corner_midpoint_offset;
    if (!this->conditionChecker(this->l1) && this->l1->CORNER)
      offset += N4 * this->l1->corner_midpoint_offset;
    if (!this->conditionChecker(this->l2) && this->l2->CORNER)
      offset += N5 * this->l2->corner_midpoint_offset;
    return offset;
  };

  template <int i, int j>
  std::array<T6d, 4> N(const double t0, const double t1) const { return N6_new<i, j>(t0, t1); };

  template <int i, int j>
  std::array<T6d, 4> D_N6_new(const double t0, const double t1) const {
    const auto [DN0, DN1, DN2, DN3, DN4, DN5] = D_TriShape_Quadratic<i, j>(t0, t1);
    T6d Nc = {0., 0., 0., DN2, DN0, DN1};
    T6d Nl0 = {0., 0., 0., 0., 0., 0.}, Nl1 = {0., 0., 0., 0., 0., 0.}, Nl2 = {0., 0., 0., 0., 0., 0.};
    if (this->conditionChecker(this->l0)) {
      Nc += 0.5 * DN3 * TriShape<6>(0.25, 0.25, this->quadpoint.available);
      Nl0 = 0.5 * DN3 * TriShape<6>(0.25, 0.25, this->quadpoint_l0.available);
    } else
      Nc += T6d{0., 0., 0., 0., 0.5 * DN3, 0.5 * DN3};

    if (this->conditionChecker(this->l1)) {
      Nc += 0.5 * DN4 * TriShape<6>(0.5, 0.25, this->quadpoint.available);
      Nl1 = 0.5 * DN4 * TriShape<6>(0.25, 0.25, this->quadpoint_l1.available);
    } else
      Nc += T6d{0., 0., 0., 0.5 * DN4, 0., 0.5 * DN4};

    if (this->conditionChecker(this->l2)) {
      Nc += 0.5 * DN5 * TriShape<6>(0.25, 0.5, this->quadpoint.available);
      Nl2 = 0.5 * DN5 * TriShape<6>(0.25, 0.25, this->quadpoint_l2.available);
    } else
      Nc += T6d{0., 0., 0., 0.5 * DN5, 0.5 * DN5, 0.};
    return {Nc, Nl0, Nl1, Nl2};
  };

  template <int i, int j>
  std::array<T6d, 4> D_N(const double t0, const double t1) const { return D_N6_new<i, j>(t0, t1); };

  std::array<double, 3> X(const double t0, const double t1) const {
    auto [Nc, N0, N1, N2] = N6_new(t0, t1);
    return Dot(Nc, ToX(this->quadpoint.points)) + Dot(N0, ToX(this->quadpoint_l0.points)) + Dot(N1, ToX(this->quadpoint_l1.points)) + Dot(N2, ToX(this->quadpoint_l2.points));
  };

  std::array<double, 3> X(const std::array<double, 2>& t0t1) const { return X(std::get<0>(t0t1), std::get<1>(t0t1)); };

  template <int i, int j>
  std::array<double, 3> D_X(const double t0, const double t1) const {
    auto [Nc, N0, N1, N2] = D_N6_new<i, j>(t0, t1);
    return Dot(Nc, ToX(this->quadpoint.points)) + Dot(N0, ToX(this->quadpoint_l0.points)) + Dot(N1, ToX(this->quadpoint_l1.points)) + Dot(N2, ToX(this->quadpoint_l2.points));
  };

  std::array<double, 3> cross(const double t0, const double t1) const { return Cross(D_X<1, 0>(t0, t1), D_X<0, 1>(t0, t1)); };

  // 最適化版: X, N6_new, crossを一度に計算（重複計算を削減）
  std::tuple<Tddd, std::array<T6d, 4>, Tddd> X_N6_cross(const double t0, const double t1) const {
    // N6_newとD_N6_newを一度に計算
    const auto [N0, N1, N2, N3, N4, N5] = TriShape<6>(t0, t1);
    const auto [DN0_t0, DN1_t0, DN2_t0, DN3_t0, DN4_t0, DN5_t0] = D_TriShape_Quadratic<1, 0>(t0, t1);
    const auto [DN0_t1, DN1_t1, DN2_t1, DN3_t1, DN4_t1, DN5_t1] = D_TriShape_Quadratic<0, 1>(t0, t1);

    // N6_new計算
    T6d Nc = {0., 0., 0., N2, N0, N1};
    T6d Nl0 = {0., 0., 0., 0., 0., 0.};
    T6d Nl1 = {0., 0., 0., 0., 0., 0.};
    T6d Nl2 = {0., 0., 0., 0., 0., 0.};

    // D_N6_new<1,0>計算
    T6d DNc_t0 = {0., 0., 0., DN2_t0, DN0_t0, DN1_t0};
    T6d DNl0_t0 = {0., 0., 0., 0., 0., 0.};
    T6d DNl1_t0 = {0., 0., 0., 0., 0., 0.};
    T6d DNl2_t0 = {0., 0., 0., 0., 0., 0.};

    // D_N6_new<0,1>計算
    T6d DNc_t1 = {0., 0., 0., DN2_t1, DN0_t1, DN1_t1};
    T6d DNl0_t1 = {0., 0., 0., 0., 0., 0.};
    T6d DNl1_t1 = {0., 0., 0., 0., 0., 0.};
    T6d DNl2_t1 = {0., 0., 0., 0., 0., 0.};

    // 条件に応じて重みを分配
    const auto shape_center = TriShape<6>(0.25, 0.25, this->quadpoint.available);
    const auto shape_l0 = TriShape<6>(0.25, 0.25, this->quadpoint_l0.available);
    const auto shape_l1 = TriShape<6>(0.25, 0.25, this->quadpoint_l1.available);
    const auto shape_l2 = TriShape<6>(0.25, 0.25, this->quadpoint_l2.available);
    const auto shape_center_l1 = TriShape<6>(0.5, 0.25, this->quadpoint.available);
    const auto shape_center_l2 = TriShape<6>(0.25, 0.5, this->quadpoint.available);

    if (this->conditionChecker(this->l0)) {
      Nc += 0.5 * N3 * shape_center;
      Nl0 = 0.5 * N3 * shape_l0;
      DNc_t0 += 0.5 * DN3_t0 * shape_center;
      DNl0_t0 = 0.5 * DN3_t0 * shape_l0;
      DNc_t1 += 0.5 * DN3_t1 * shape_center;
      DNl0_t1 = 0.5 * DN3_t1 * shape_l0;
    } else {
      Nc += T6d{0., 0., 0., 0., 0.5 * N3, 0.5 * N3};
      DNc_t0 += T6d{0., 0., 0., 0., 0.5 * DN3_t0, 0.5 * DN3_t0};
      DNc_t1 += T6d{0., 0., 0., 0., 0.5 * DN3_t1, 0.5 * DN3_t1};
    }

    if (this->conditionChecker(this->l1)) {
      Nc += 0.5 * N4 * shape_center_l1;
      Nl1 = 0.5 * N4 * shape_l1;
      DNc_t0 += 0.5 * DN4_t0 * shape_center_l1;
      DNl1_t0 = 0.5 * DN4_t0 * shape_l1;
      DNc_t1 += 0.5 * DN4_t1 * shape_center_l1;
      DNl1_t1 = 0.5 * DN4_t1 * shape_l1;
    } else {
      Nc += T6d{0., 0., 0., 0.5 * N4, 0., 0.5 * N4};
      DNc_t0 += T6d{0., 0., 0., 0.5 * DN4_t0, 0., 0.5 * DN4_t0};
      DNc_t1 += T6d{0., 0., 0., 0.5 * DN4_t1, 0., 0.5 * DN4_t1};
    }

    if (this->conditionChecker(this->l2)) {
      Nc += 0.5 * N5 * shape_center_l2;
      Nl2 = 0.5 * N5 * shape_l2;
      DNc_t0 += 0.5 * DN5_t0 * shape_center_l2;
      DNl2_t0 = 0.5 * DN5_t0 * shape_l2;
      DNc_t1 += 0.5 * DN5_t1 * shape_center_l2;
      DNl2_t1 = 0.5 * DN5_t1 * shape_l2;
    } else {
      Nc += T6d{0., 0., 0., 0.5 * N5, 0.5 * N5, 0.};
      DNc_t0 += T6d{0., 0., 0., 0.5 * DN5_t0, 0.5 * DN5_t0, 0.};
      DNc_t1 += T6d{0., 0., 0., 0.5 * DN5_t1, 0.5 * DN5_t1, 0.};
    }

    // 座標データを取得
    const auto X_c = ToX(this->quadpoint.points);
    const auto X_l0 = ToX(this->quadpoint_l0.points);
    const auto X_l1 = ToX(this->quadpoint_l1.points);
    const auto X_l2 = ToX(this->quadpoint_l2.points);

    // X計算
    Tddd X = Dot(Nc, X_c) + Dot(Nl0, X_l0) + Dot(Nl1, X_l1) + Dot(Nl2, X_l2);

    // D_X計算
    Tddd dXdt0 = Dot(DNc_t0, X_c) + Dot(DNl0_t0, X_l0) + Dot(DNl1_t0, X_l1) + Dot(DNl2_t0, X_l2);
    Tddd dXdt1 = Dot(DNc_t1, X_c) + Dot(DNl0_t1, X_l0) + Dot(DNl1_t1, X_l1) + Dot(DNl2_t1, X_l2);

    // CORNER edge quadratic correction (offset is {0,0,0} when not active)
    if (!this->conditionChecker(this->l0) && this->l0->CORNER) {
      X += N3 * this->l0->corner_midpoint_offset;
      dXdt0 += DN3_t0 * this->l0->corner_midpoint_offset;
      dXdt1 += DN3_t1 * this->l0->corner_midpoint_offset;
    }
    if (!this->conditionChecker(this->l1) && this->l1->CORNER) {
      X += N4 * this->l1->corner_midpoint_offset;
      dXdt0 += DN4_t0 * this->l1->corner_midpoint_offset;
      dXdt1 += DN4_t1 * this->l1->corner_midpoint_offset;
    }
    if (!this->conditionChecker(this->l2) && this->l2->CORNER) {
      X += N5 * this->l2->corner_midpoint_offset;
      dXdt0 += DN5_t0 * this->l2->corner_midpoint_offset;
      dXdt1 += DN5_t1 * this->l2->corner_midpoint_offset;
    }

    return {X, {Nc, Nl0, Nl1, Nl2}, Cross(dXdt0, dXdt1)};
  };

  /*
  ６点の値を与えた場合，標準の2次補間を使えるようにする．
  */

  std::array<double, 3> X(const double t0, const double t1, const std::array<std::array<double, 3>, 6>& samples) const { return Dot(TriShape<6>(t0, t1), samples); };
  std::array<double, 3> X(const std::array<double, 2>& t0t1, const std::array<std::array<double, 3>, 6>& samples) const { return Dot(TriShape<6>(std::get<0>(t0t1), std::get<1>(t0t1)), samples); };

  template <int i, int j>
  std::array<double, 3> D_X(const double t0, const double t1, const std::array<std::array<double, 3>, 6>& samples) const { return Dot(D_TriShape<6, i, j>(t0, t1), samples); };

  std::array<double, 3> cross(const double t0, const double t1, const std::array<std::array<double, 3>, 6>& samples) const { return Cross(D_X<1, 0>(t0, t1, samples), D_X<0, 1>(t0, t1, samples)); };

  /* --------------------------------------------------------------------------
   */

  template <typename T>
  T Dot4(const std::array<T6d, 4>& N, const std::array<std::array<T, 6>, 4>& samples) const {
    // 変更なし
    return Dot(std::get<0>(N), std::get<0>(samples)) + Dot(std::get<1>(N), std::get<1>(samples)) + Dot(std::get<2>(N), std::get<2>(samples)) + Dot(std::get<3>(N), std::get<3>(samples));
  };

  template <typename ConversionFunc>
  auto conversions(ConversionFunc conversion) const {
    // 変更なし
    using T = std::invoke_result_t<ConversionFunc, networkPoint*>;
    std::array<std::array<T, 6>, 4> samples;
    for (auto i = 0; i < 6; ++i) {
      std::get<0>(samples)[i] = conversion(this->quadpoint.points[i]);
      std::get<1>(samples)[i] = conversion(this->quadpoint_l0.points[i]);
      std::get<2>(samples)[i] = conversion(this->quadpoint_l1.points[i]);
      std::get<3>(samples)[i] = conversion(this->quadpoint_l2.points[i]);
    }
    return samples;
  };

  template <typename ConversionFunc>
  auto interpolate(const double t0, const double t1, ConversionFunc conversion) const {
    using ReturnT = std::invoke_result_t<ConversionFunc, networkPoint*>;
    if constexpr (is_std_array_of_6_v<ReturnT>) {
      auto samples = conversion(this->f->getPoints()[0]); // 引数のPointはダミー
      return Dot(TriShape<6>(t0, t1), samples);
    } else {
      auto result = this->Dot4(N6_new(t0, t1), conversions(conversion));
      if constexpr (std::is_same_v<ReturnT, Tddd>)
        result += corner_offset(t0, t1);
      return result;
    }
  };

  //! ----- interpolate　強制的に与えられたサンプル点を使って補間する関数
  //! (完全な2次補間用) -----

  // template <typename ConversionFunc>
  // auto interpolate(const double t0, const double t1, const
  // std::array<std::array<double, 3>, 6> &samples, ConversionFunc conversion)
  // const {
  //    using ReturnT = std::invoke_result_t<ConversionFunc, networkPoint *>;
  //    if constexpr (is_std_array_of_6_v<ReturnT>) {
  //       auto samples = conversion(this->f->getPoints()[0]);  //
  //       引数のPointはダミー return Dot(TriShape<6>(t0, t1), samples);
  //    } else {
  //       std::array<std::array<T, 3>, 6> conv;
  //       std::get<0>(conv)[i] = conversion(samples_points[0]);
  //       std::get<1>(conv)[i] = conversion(samples_points[1]);
  //       std::get<2>(conv)[i] = conversion(samples_points[2]);
  //       std::get<3>(conv)[i] = conversion(samples_points[3]);
  //       std::get<4>(conv)[i] = conversion(samples_points[4]);
  //       std::get<5>(conv)[i] = conversion(samples_points[5]);
  //       return Dot(D_TriShape<6, i, j>(t0, t1), conv);
  //    }
  // };

  //! -----------------------------------------------------------------

  template <typename ConversionFunc>
  auto interpolate(const Tdd t0t1, ConversionFunc conversion) const { return interpolate(std::get<0>(t0t1), std::get<1>(t0t1), conversion); };

  // ベクトルや配列を補間するバージョンも同様に if constexpr で分岐
  template <std::size_t N, typename ConversionFunc>
  auto interpolate(const std::array<Tdd, N>& t0t1, ConversionFunc conversion) const {
    using ReturnT = std::invoke_result_t<ConversionFunc, networkPoint*>;

    if constexpr (is_std_array_of_6_v<ReturnT>) {
      using T_elem = typename ReturnT::value_type;
      std::array<T_elem, N> ret;
      auto samples = conversion(this->f->getPoints()[0]);
      for (auto i = 0; i < t0t1.size(); ++i)
        ret[i] = Dot(TriShape<6>(std::get<0>(t0t1[i]), std::get<1>(t0t1[i])), samples);
      return ret;
    } else {
      using T_elem = ReturnT;
      std::array<T_elem, N> ret;
      auto samples = conversions(conversion);
      for (auto i = 0; i < t0t1.size(); ++i) {
        ret[i] = this->Dot4(N6_new(std::get<0>(t0t1[i]), std::get<1>(t0t1[i])), samples);
        if constexpr (std::is_same_v<T_elem, Tddd>)
          ret[i] += corner_offset(std::get<0>(t0t1[i]), std::get<1>(t0t1[i]));
      }
      return ret;
    }
  };
  // 1つ目の関数：std::vector<std::array<Tdd, N>> を引数に取るバージョン
  template <std::size_t N, typename ConversionFunc>
  auto interpolate(const std::vector<std::array<Tdd, N>>& t0t1, ConversionFunc conversion) const {
    using ReturnT = std::invoke_result_t<ConversionFunc, networkPoint*>;

    if constexpr (is_std_array_of_6_v<ReturnT>) {
      // 【新機能】6節点の値を直接与える場合
      using T_elem = typename ReturnT::value_type; // Tdddなどを抽出
      std::vector<std::array<T_elem, N>> ret(t0t1.size());
      auto samples = conversion(this->f->getPoints()[0]); // データを一度だけ取得
      for (auto i = 0; i < t0t1.size(); ++i)
        for (auto j = 0; j < N; ++j)
          ret[i][j] = Dot(TriShape<6>(std::get<0>(t0t1[i][j]), std::get<1>(t0t1[i][j])), samples);
      return ret;
    } else {
      // 【既存機能】隣接面を考慮する場合
      using T_elem = ReturnT; // Tdddなどを抽出
      std::vector<std::array<T_elem, N>> ret(t0t1.size());
      auto samples = conversions(conversion);
      for (auto i = 0; i < t0t1.size(); ++i)
        for (auto j = 0; j < N; ++j) {
          ret[i][j] = this->Dot4(N6_new(std::get<0>(t0t1[i][j]), std::get<1>(t0t1[i][j])), samples);
          if constexpr (std::is_same_v<T_elem, Tddd>)
            ret[i][j] += corner_offset(std::get<0>(t0t1[i][j]), std::get<1>(t0t1[i][j]));
        }
      return ret;
    }
  };

  // 2つ目の関数：std::vector<Tdd> を引数に取るバージョン
  template <typename ConversionFunc>
  auto interpolate(const std::vector<Tdd>& t0t1, ConversionFunc conversion) const {
    using ReturnT = std::invoke_result_t<ConversionFunc, networkPoint*>;

    if constexpr (is_std_array_of_6_v<ReturnT>) {
      // 【新機能】6節点の値を直接与える場合
      using T_elem = typename ReturnT::value_type;
      std::vector<T_elem> ret(t0t1.size());
      auto samples = conversion(this->f->getPoints()[0]);
      // ★バグ修正：ループの上限をt0t1.size()に変更
      for (auto i = 0; i < t0t1.size(); ++i)
        ret[i] = Dot(TriShape<6>(std::get<0>(t0t1[i]), std::get<1>(t0t1[i])), samples);
      return ret;
    } else {
      // 【既存機能】隣接面を考慮する場合
      using T_elem = ReturnT;
      std::vector<T_elem> ret(t0t1.size());
      auto samples = conversions(conversion);
      // ★バグ修正：ループの上限をt0t1.size()に変更
      for (auto i = 0; i < t0t1.size(); ++i)
        ret[i] = this->Dot4(N6_new(std::get<0>(t0t1[i]), std::get<1>(t0t1[i])), samples);
      return ret;
    }
  };

  template <int i, int j, typename ConversionFunc>
  auto D_interpolate(const double t0, const double t1, ConversionFunc conversion) const {
    using ReturnT = std::invoke_result_t<ConversionFunc, networkPoint*>;
    if constexpr (is_std_array_of_6_v<ReturnT>) {
      auto samples = conversion(this->f->getPoints()[0]); // 引数のPointはダミー
      return Dot(D_TriShape<6, i, j>(t0, t1), samples);
    } else {
      return Dot4(D_N6_new<i, j>(t0, t1), conversions(conversion));
    }
  };

  template <int i, int j, typename ConversionFunc>
  auto D_interpolate(const Tdd t0t1, ConversionFunc conversion) const { return D_interpolate<i, j>(std::get<0>(t0t1), std::get<1>(t0t1), conversion); };

  template <int i, int j, typename ConversionFunc>
  auto D_interpolate(const std::vector<Tdd> t0t1, ConversionFunc conversion) const {
    using ReturnT = std::invoke_result_t<ConversionFunc, networkPoint*>;

    if constexpr (is_std_array_of_6_v<ReturnT>) {
      using T_elem = typename ReturnT::value_type;
      std::vector<T_elem> ret(t0t1.size());
      auto samples = conversion(this->f->getPoints()[0]);
      for (auto k = 0; k < t0t1.size(); ++k)
        ret[k] = Dot(D_TriShape<6, i, j>(std::get<0>(t0t1[k]), std::get<1>(t0t1[k])), samples);
      return ret;
    } else {
      using T_elem = ReturnT;
      std::vector<T_elem> ret(t0t1.size());
      auto samples = conversions(conversion);
      for (auto k = 0; k < t0t1.size(); ++k)
        ret[k] = this->Dot4(D_N6_new<i, j>(std::get<0>(t0t1[k]), std::get<1>(t0t1[k])), samples);
      return ret;
    }
  };
};

/* -------------------------------------------------------------------------- */
template <typename U, typename T>
typename std::enable_if<std::is_same<T, Network*>::value || std::is_same<T, std::vector<networkLine*>>::value || std::is_same<T, std::unordered_set<networkLine*>>::value || std::is_same<T, std::vector<networkFace*>>::value || std::is_same<T, std::unordered_set<networkFace*>>::value || std::is_same<T, std::vector<networkTetra*>>::value || std::is_same<T, std::unordered_set<networkTetra*>>::value>::type vtkPolygonWrite(std::ofstream& ofs, const T& uo_f,
                                                                                                                                                                                                                                                                                                                                                                                                                                     const std::unordered_map<networkPoint*, U>& data_double = {}) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uo_f.size());
  for (const auto& f : uo_f) {
    auto abc = f->getPoints();
    vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
    vtp.addPolygon(abc);
  }
  if (!data_double.empty())
    vtp.addPointData("data", data_double);
  vtp.write(ofs);
}

template <typename face_list>
void vtkPolygonWrite(std::ofstream& ofs, const face_list& uo_f, const std::vector<std::tuple<std::string, std::unordered_map<networkPoint*, std::variant<double, std::array<double, 3>>>>>& name_uo_data) {
  vtkUnstructuredGridWriter<networkPoint*> vtu;
  vtu.reserve(uo_f.size());
  for (const auto& f : uo_f) {
    auto abc = f->getPoints();
    vtu.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
    vtu.addGrid(abc);
  }
  for (auto i = 0; i < name_uo_data.size(); ++i) {
    std::unordered_map<networkPoint*, double> double_map;
    std::unordered_map<networkPoint*, std::array<double, 3>> array_map;

    for (const auto& [key, value] : std::get<1>(name_uo_data[i])) {
      std::visit(
          [&](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, double>) {
              double_map[key] = arg;
            } else if constexpr (std::is_same_v<T, std::array<double, 3>>) {
              array_map[key] = arg;
            }
          },
          value);
    }

    if (!double_map.empty()) {
      vtu.addPointData(std::get<0>(name_uo_data[i]), double_map);
    }
    if (!array_map.empty()) {
      vtu.addPointData(std::get<0>(name_uo_data[i]), array_map);
    }
  }
  vtu.write(ofs);
};

inline void vtkPolygonWrite(std::ofstream& ofs, const T_4F& tuple_f) {
  vtkPolygonWriter<networkPoint*> vtp;
  std::ranges::for_each(tuple_f, [&](const auto& f) {
    auto abc = f->getPoints();
    vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
    vtp.addPolygon(abc); });
  vtp.write(ofs);
};

inline void vtkPolygonWrite(std::ofstream& ofs, const std::unordered_set<networkLine*>& uoL) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoL.size());
  for (const auto& l : uoL) {
    auto [a, b] = l->getPoints();
    vtp.add(a, b);
    vtp.addLine(a, b);
  }
  vtp.write(ofs);
};

template <typename U>
void vtkPolygonWrite(std::ofstream& ofs, const std::unordered_set<networkLine*>& uoL, const std::unordered_map<networkPoint*, U>& data_double) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoL.size());
  for (const auto& l : uoL) {
    auto [a, b] = l->getPoints();
    vtp.add(a, b);
    vtp.addLine(a, b);
  }
  vtp.addPointData("data", data_double);
  vtp.write(ofs);
};

// for const std::vector<std::tuple<std::string, std::unordered_map<T, U>>
template <typename U>
void vtkPolygonWrite(std::ofstream& ofs, const std::unordered_set<networkLine*>& uoL, const std::vector<std::tuple<std::string, std::unordered_map<networkPoint*, U>>>& name_uo_data = {}) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoL.size());
  for (const auto& l : uoL) {
    auto [a, b] = l->getPoints();
    vtp.add(a, b);
    vtp.addLine(a, b);
  }
  for (auto i = 0; i < name_uo_data.size(); ++i)
    vtp.addPointData(std::get<0>(name_uo_data[i]), std::get<1>(name_uo_data[i]));
  vtp.write(ofs);
};

// for std::variant<double, std::array<double, 3>>

inline void vtkPolygonWrite(std::ofstream& ofs, const std::unordered_set<networkPoint*>& uoP, const std::vector<std::tuple<std::string, std::unordered_map<networkPoint*, std::variant<double, std::array<double, 3>>>>>& name_uo_data = {}) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoP.size());
  for (const auto& p : uoP) {
    vtp.add(p);
  }

  for (const auto& [name, data] : name_uo_data) {
    std::unordered_map<networkPoint*, double> double_map;
    std::unordered_map<networkPoint*, std::array<double, 3>> array_map;

    for (const auto& [key, value] : data) {
      std::visit(
          [&](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, double>) {
              double_map[key] = arg;
            } else if constexpr (std::is_same_v<T, std::array<double, 3>>) {
              array_map[key] = arg;
            }
          },
          value);
    }

    if (!double_map.empty()) {
      vtp.addPointData(name, double_map);
    }
    if (!array_map.empty()) {
      vtp.addPointData(name, array_map);
    }
  }

  vtp.write(ofs);
}

inline void vtkUnstructuredGridWrite(std::ofstream& ofs, const std::unordered_set<networkFace*>& uoF) {
  vtkUnstructuredGridWriter<networkPoint*> vtu;
  vtu.reserve(uoF.size());
  for (const auto& f : uoF) {
    auto abc = f->getPoints();
    vtu.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
    vtu.addGrid(abc);
  };
  vtu.write(ofs);
};

inline void vtkPolygonWrite(std::ofstream& ofs, const std::unordered_set<networkFace*>& uoF) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoF.size());
  for (const auto& f : uoF) {
    auto abc = f->getPoints();
    vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
    vtp.addPolygon(abc);
  };
  vtp.write(ofs);
};

inline void vtkPolygonWrite(std::ofstream& ofs, const std::vector<networkFace*>& uoF, const std::vector<std::tuple<std::string, std::unordered_map<networkPoint*, double>>>& name_uo_data) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoF.size());

  for (const auto& f : uoF) {
    auto abc = f->getPoints();
    vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
    vtp.addPolygon(abc);
  }

  for (const auto& [name, data] : name_uo_data) {
    vtp.addPointData(name, data);
  }

  vtp.write(ofs);
}

inline void vtkPolygonWrite(std::ofstream& ofs, const std::unordered_set<networkTetra*>& uoTet) {
  // vtuがいいならそうなるように修正しなければならない．
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoTet.size());
  for (const auto& tet : uoTet)
    std::ranges::for_each(tet->Faces, [&](const auto& f) {
      auto abc = f->getPoints();
      vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
      vtp.addPolygon(abc); });
  // for (const auto &tet : uoTet) {
  //    vtp.add(std::get<0>(tet->Points),
  //            std::get<1>(tet->Points),
  //            std::get<2>(tet->Points),
  //            std::get<3>(tet->Points));
  //    vtp.addPolygon(tet->Points);
  // };
  vtp.write(ofs);
};

inline void vtkPolygonWrite(std::ofstream& ofs, const std::vector<networkTetra*>& uoTet) {
  vtkPolygonWriter<networkPoint*> vtp;
  vtp.reserve(uoTet.size());
  for (const auto& tet : uoTet)
    std::ranges::for_each(tet->Faces, [&](const auto& f) {
      auto abc = f->getPoints();
      vtp.add(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc));
      vtp.addPolygon(abc); });
  vtp.write(ofs);
};
template <typename T>
typename std::enable_if<std::is_same<T, std::vector<networkTetra*>>::value || std::is_same<T, std::unordered_set<networkTetra*>>::value>::type vtkUnstructuredGridWrite(std::ofstream& ofs, const T& uoTet) {
  vtkUnstructuredGridWriter<networkPoint*> vtu;
  vtu.reserve(uoTet.size());
  for (const auto& tet : uoTet) {
    auto abcd = tet->Points;
    vtu.add(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd));
    vtu.addGrid(abcd);
  }
  vtu.write(ofs);
};

// Correctly define the `DataVariant`
using DataVariant3 = std::variant<std::unordered_map<networkPoint*, double>, std::unordered_map<networkPoint*, Tddd>, std::unordered_map<networkPoint*, std::variant<double, Tddd>>>;

// Function Definition
template <typename T>
typename std::enable_if<std::is_same<T, std::vector<networkTetra*>>::value || std::is_same<T, std::unordered_set<networkTetra*>>::value>::type vtkUnstructuredGridWrite(std::ofstream& ofs, const T& uoTet, const std::vector<std::tuple<std::string, DataVariant3>>& name_uo_data) {
  if (uoTet.empty())
    return;

  vtkUnstructuredGridWriter<networkPoint*> vtu;
  vtu.reserve(uoTet.size());

  for (const auto& tet : uoTet) {
    auto abcd = tet->Points;
    vtu.add(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd));
    vtu.addGrid(abcd);
  }

  for (const auto& [name, data] : name_uo_data) {
    std::unordered_map<networkPoint*, double> double_map;
    std::unordered_map<networkPoint*, Tddd> array_map;

    std::visit(
        [&](auto&& arg) {
          using VType = std::decay_t<decltype(arg)>;

          if constexpr (std::is_same_v<VType, std::unordered_map<networkPoint*, double>>) {
            double_map = arg;
          } else if constexpr (std::is_same_v<VType, std::unordered_map<networkPoint*, Tddd>>) {
            array_map = arg;
          } else if constexpr (std::is_same_v<VType, std::unordered_map<networkPoint*, std::variant<double, Tddd>>>) {
            for (const auto& [key, value] : arg) {
              std::visit(
                  [&](auto&& val) {
                    if constexpr (std::is_same_v<std::decay_t<decltype(val)>, double>) {
                      double_map[key] = val;
                    } else if constexpr (std::is_same_v<std::decay_t<decltype(val)>, Tddd>) {
                      array_map[key] = val;
                    }
                  },
                  value);
            }
          }
        },
        data);

    if (!double_map.empty()) {
      vtu.addPointData(name, double_map);
    }
    if (!array_map.empty()) {
      vtu.addPointData(name, array_map);
    }
  }

  vtu.write(ofs);
}

/*

重要：
------------- networkPoint ------------
std::vector<networkLine *> networkPoint::Lines
std::vector<networkFace *> networkPoint::Faces

------------- networkLine ------------
networkLine::Point_A
networkLine::Point_B
std::vector<networkFace *> networkLine::Face

------------- networkFace ------------
std::array<networkPoint *,3> networkFace::Points
std::array<networkLine *,3> networkFace::Lines

-------------

l: f  -> F の入れ替えでは
f->erase(l)
F->add(l)
とすることができ．

l or f :p->Pの入れ替えでは，
p->erase(l or f)
P->add(l or f)
とすることができる．

*/

/* -------------------------------------------------------------------------- */

void dual_replace(netP* p, netL* line_before, netL* line_after);
void dual_replace(netP* p, netF* face_before, netF* face_after);

/* -------------------------------------------------------------------------- */

void dual_replace(netL* l, netP* point_before, netP* point_after);
void dual_replace(netL* l, netF* face_before, netF* face_after, netL* l_face_before, netL* l_face_after);

/* -------------------------------------------------------------------------- */

void dual_replace(netF* f, netL* line_before, netL* line_after);
void dual_replace(netF* f, netL* line_before, netL* line_after, netF* f_line_before);
void dual_replace(netF* f, netL* line_before, netL* line_after, netF* f_line_before, netF* f_line_after);
void dual_replace(netF* f, netP* point_before, netP* point_after);
void dual_replace(netF* f, netP* point_before, netP* point_after, netF* f_point_before);

/* -------------------------------------------------------------------------- */

inline bool networkPoint::BoundaryQ() const noexcept {
  return std::any_of(this->Faces.begin(), this->Faces.end(), [](const auto& f) { return f->BoundaryQ(); });
}

inline bool networkPoint::ComputationalBoundaryQ() const noexcept {
  return std::any_of(this->Faces.begin(), this->Faces.end(), [](const auto& f) { return f->ComputationalBoundaryQ(); });
}

inline Tddd tetraCentroid(const networkTetra* t) {
  auto [p0, p1, p2, p3] = t->Points;
  return (p0->X + p1->X + p2->X + p3->X) * 0.25;
}

inline bool isInsideTetra(const networkFace* face, const networkTetra* tetra) {
  if (!face || !tetra) return false;
  return Dot(tetraCentroid(tetra) - face->X, face->normal) < 0.0;
}

inline bool isOutsideTetra(const networkFace* face, const networkTetra* tetra) {
  if (!face || !tetra) return false;
  return Dot(tetraCentroid(tetra) - face->X, face->normal) > 0.0;
}

inline bool networkFace::ComputationalBoundaryQ() const noexcept {
  auto* t0 = this->Tetras[0];
  auto* t1 = this->Tetras[1];
  if (!t0 && !t1)
    return false;  // 両面 nullptr → 孤立面
  if (!t0 || !t1) {
    // 片面のみ tetra あり
    auto* t = t0 ? t0 : t1;
    if (t->tetra_state == TetraState::Inside)
      return true;   // Inside tetra に接する蓋面 → 計算境界
    if (t->tetra_state == TetraState::Outside)
      return false;  // Outside tetra の蓋面 → 計算境界ではない
    return true;     // Water/Air 等（boundary tetra 未生成時）→ BoundaryQ() と同じ
  }
  // 両面に tetra あり → Inside/Outside が異なれば計算境界
  bool t0_inside = Dot(tetraCentroid(t0) - this->X, this->normal) < 0.0;
  bool t1_inside = Dot(tetraCentroid(t1) - this->X, this->normal) < 0.0;
  return t0_inside != t1_inside;
}

inline networkTetra* networkFace::getInsideTetra() const noexcept {
  for (auto* t : this->Tetras) {
    if (!t) continue;
    if (Dot(tetraCentroid(t) - this->X, this->normal) < 0.0)
      return t;
  }
  return nullptr;
}

inline networkTetra* networkFace::getOutsideTetra() const noexcept {
  for (auto* t : this->Tetras) {
    if (!t) continue;
    if (Dot(tetraCentroid(t) - this->X, this->normal) > 0.0)
      return t;
  }
  return nullptr;
}

#include "MooringLine.hpp"
#include "NetworkUtility.hpp"
// #include "networkFace.cpp"
// #include "networkLine.cpp"
// #include "networkTetra.cpp"
