#pragma once

#include "basic.hpp"
#include "interpolations.hpp"

// @ ============================================================== */
// @                        剛体の力学に関する                         */
// @ ============================================================== */
//@ メッシュが関わらない剛体の運動を表すクラス
struct RigidBodyDynamics {

  RigidBodyDynamics() : quaternion(Quaternion()) {};
  virtual ~RigidBodyDynamics() = default;

  T6d forced_velocity = {0., 0., 0., 0., 0., 0.};
  T6d forced_acceleration = {0., 0., 0., 0., 0., 0.};

  InterpolationBspline<std::array<double, 6>> intpMotionRigidBody;
  //* ------------------------------------------------------ */
  //*                     運動を表す量                         */
  //* ------------------------------------------------------ */
  T6d force = {0., 0., 0., 0., 0., 0.};
  T6d velocity = {0., 0., 0., 0., 0., 0.}; // = {velocity,angular velocity}
  T6d acceleration = {0., 0., 0., 0., 0., 0.};
  InterpolationLagrange<T6d> interp_accel;
  Tddd velocityTranslational() const { return {std::get<0>(velocity), std::get<1>(velocity), std::get<2>(velocity)}; };
  Tddd velocityRotational() const { return {std::get<3>(velocity), std::get<4>(velocity), std::get<5>(velocity)}; };
  Tddd accelTranslational() const { return {std::get<0>(acceleration), std::get<1>(acceleration), std::get<2>(acceleration)}; };
  Tddd accelRotational() const { return {std::get<3>(acceleration), std::get<4>(acceleration), std::get<5>(acceleration)}; };

  //! 固定された空間座標におけるベクトルであることを頭に入れておくこと．
  //! 回転，移動をする物体の座標系ではないので，固定座標にとって，回転前後でinertiaは書き換える必要がある．
  //! inertiaの慣性モーメントはそのまま固定座標における回転行列をかけて，更新すればいい
  T6d &F = this->force;
  T6d &A = this->acceleration;
  T6d &V = this->velocity;
  //* ------------------------------------------------------ */
  //*                   位置や姿勢を表す量                      */
  //* ------------------------------------------------------ */
  Tddd center_of_mass = {0., 0., 0.}; // 現在の座標を表す
  Quaternion quaternion;              // 現在の姿勢を表す，初期のクォータニオンは固定で{1,0,0,0}
  Tddd &COM = this->center_of_mass;
  Quaternion &Q = this->quaternion;
  //* ------------------------------------------------------ */
  //*                   　　 不変の量                          */
  //* ------------------------------------------------------ */
  Tddd initial_center_of_mass = {0., 0., 0.}; // 初期のの座標を表す
  Tddd &ICOM = this->initial_center_of_mass;
  double mass = 1E+20;
  T6d inertia = {1E+20, 1E+20, 1E+20, 1E+20, 1E+20, 1E+20};
  T6d &I = this->inertia;

  std::tuple<double, double, double, double, double, double, T3Tddd, T3Tddd> getInertiaGC(double alpha = 1.0) // Global coordinate
  {
    auto [mx, my, mz, Ix, Iy, Iz] = this->inertia * alpha;
    auto R = this->quaternion.Rv();
    // Add a check for very large or zero inertia values to prevent NaN/inf
    if (Ix <= 1e-12 || Iy <= 1e-12 || Iz <= 1e-12) {
      std::cout << this->inertia << std::endl;
      std::cerr << "Warning: small moment of inertia components (Ix, Iy, Iz): " << Ix << ", " << Iy << ", " << Iz << ". IG_inv may contain inf or NaN values." << std::endl;
    } else if (Ix > 1e18 || Iy > 1e18 || Iz > 1e18) {
      std::cerr << "Warning: Extremely large moment of inertia components (Ix, Iy, Iz): " << Ix << ", " << Iy << ", " << Iz << ". IG_inv may contain inf or NaN values." << std::endl;
    }
    T3Tddd IG = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
    T3Tddd IG_inv = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
    for (auto i = 0; i < 3; ++i) {
      for (auto j = 0; j < 3; ++j) {
        IG[i][j] += R[0][i] * R[0][j] * Ix;
        IG[i][j] += R[1][i] * R[1][j] * Iy;
        IG[i][j] += R[2][i] * R[2][j] * Iz;
        IG_inv[i][j] += R[0][i] * R[0][j] / Ix;
        IG_inv[i][j] += R[1][i] * R[1][j] / Iy;
        IG_inv[i][j] += R[2][i] * R[2][j] / Iz;
      }
    }
    return {mx, my, mz, Ix, Iy, Iz, IG, IG_inv};
  };

  T6d getInertiaBC() { return this->inertia; }; // Body coordinate
  Tddd getMass3D() { return {std::get<0>(this->inertia), std::get<1>(this->inertia), std::get<2>(this->inertia)}; };

  /* ------------------------------------------------------ */

  Tddd velocityRigidBody(const Tddd &X) const { return velocityTranslational() + Cross(velocityRotational(), X - this->COM); }; //! \label{velocityRigidBody}
  Tddd accelRigidBody(const Tddd &X) const {
    Tddd r = X - this->COM;
    Tddd w = velocityRotational();
    return accelTranslational() + Cross(accelRotational(), r) - CrossCross(w, r, w) /*Cross(w, Cross(w, r))*/;
  }; //! \label{accelRigidBody}
  void calcAccelFromForce() {
    //! 慣性を設定しておく必要がある．
    this->acceleration = this->force / this->inertia;
  };
  std::array<double, 3> rigidTransformation(const std::array<double, 3> &initial_position) const {
    auto rotation = Dot(this->quaternion.Rv(), initial_position - this->initial_center_of_mass);
    auto translation = this->center_of_mass;
    return rotation + translation;
  };
};

// @ ============================================================== */
// @                        弾性体の力学に関する                         */
// @ ============================================================== */
//@ メッシュが関わらない弾性体の運動を表すクラス
struct ElasticBodyDynamics {

  ElasticBodyDynamics() = default;
  virtual ~ElasticBodyDynamics() = default;

  //* ------------------------------------------------------ */
  //*                     運動を表す量                         */
  //* ------------------------------------------------------ */
  T6d force = {0., 0., 0., 0., 0., 0.};
  T6d velocity = {0., 0., 0., 0., 0., 0.}; // = {velocity, angular velocity}
  T6d acceleration = {0., 0., 0., 0., 0., 0.};

  //* ------------------------------------------------------ */
  //*                   位置や変形を表す量                      */
  //* ------------------------------------------------------ */
  Tddd center_of_mass = {0., 0., 0.}; // 現在の重心座標
  Quaternion quaternion;              // 回転（必要なら）

  //* ------------------------------------------------------ */
  //*                   　　 材料特性                          */
  //* ------------------------------------------------------ */
  double mass = 1.0;
  double youngs_modulus = 1.0e+9; // ヤング率
  double poisson_ratio = 0.3;     // ポアソン比
};
