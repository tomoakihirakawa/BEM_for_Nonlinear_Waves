#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Network.hpp"
#include "integrationOfODE.hpp"

/*DOC_EXTRACT 1_LumpedCable

## 集中質量ケーブル — `Network` の派生クラス

`LumpedCable` は `Network` を派生し、一本のケーブル（橋梁 stay, 海洋係留、
弦、ロープ）を「集中質量点 + 線形ばね要素」の鎖として表現する。旧名は
`MooringLine` で、今でも `MooringLine` は `LumpedCable` へのエイリアス
として利用できる（後方互換のため）。

`networkLine` には `natural_length`, `stiffness`, `damping`,
`weight_per_unit_length` の 4 パラメタ、`LumpedCable` 自身は `total_length`
を持つ。`natural_length` は `total_length / (n_points - 1)` で決まる。

```cpp
const int n_points = 100;
const double total_length = 522.; //! [m]
std::array<double, 3> X_anchor   = {r * std::cos(0.), r * std::sin(0.), h};
std::array<double, 3> X_fairlead = {0., 0., 0.};
auto cable = new LumpedCable(X_anchor, X_fairlead, total_length, n_points);
```

物性はまとめて設定する：

```cpp
cable->setDensityStiffnessDampingDiameter(mass_per_length, stiffness, damp, diam);
```

第 1 引数は **単位長あたり質量 [kg/m]**（`weight_per_unit_length` メンバに代入
されるが、実態は質量密度）。Phase 0.5 (2026-04-12) で引数名を `density` →
`mass_per_length` にリネームして意味を一致させた。

典型的な使い方：フェアリード位置が既知で、そこの張力を知りたい。前時刻の
形状から RK4 で時間発展させ、現時刻の形状・張力・力を得る。動く節点
（フェアリード等）には境界条件を与え、自由な節点は張力・重力・抗力から
運動方程式を解いて進める。

### 時間積分メソッドの 3 つの層（Phase 0.5 で整理）

1. **`step(t, dt, BC)`** — 純粋な RK4 1 ステップ進行。内部で dt ランプを行わず、
   `DragForceCoefficient` はインスタンス値をそのまま使い、stdout には何も
   出さない。BEM 時間ループ、または `LumpedCableSystem::advanceRKStage()`
   から呼ばれる。
2. **`simulate(t, dt, BC, silent=false)`** — 既存の dt ランプ付きサブステップ
   ループ。ユーザが与えた 1 回の dt に対し、内部で 1e-8..1e-3 の幅でマイクロ
   ステップを積み上げながら進む。`setEquilibriumState()` の内部で使用。
   BEM 時間ループから直接呼ぶ用途には dt ランプが非効率なので、代わりに
   `step()` を使うこと。
3. **`setEquilibriumState(BC, tol, max_iters)`** — 疑似緩和による静的平衡
   探索。内部で `DragForceCoefficient` を一時的に 1000 に上げ、`simulate()`
   を収束するまで反復する。Phase 0.5 で収束判定 break が正しく動くように
   修正した（旧版はコメントアウトで常に 100 回固定）。

*/

/* -------------------------------------------------------------------------- */

class LumpedCable : public Network {
  public:
   // default destructor
   ~LumpedCable() = default;

   double total_length;
   networkPoint* firstPoint = nullptr;
   networkPoint* lastPoint = nullptr;

   LumpedCable(const std::array<double, 3> X0,
               const std::array<double, 3> X1,
               const double total_length,
               const int n)
       : Network(), total_length(total_length) {
      std::vector<networkPoint*> points;
      for (int i = 0; i < n; ++i) {
         auto p = new networkPoint(this, X0 + (X1 - X0) * i / (n - 1.));
         points.push_back(p);
         if (i == 0) firstPoint = p;
         if (i == n - 1) lastPoint = p;
      }
      for (auto i = 0; i < points.size() - 1; ++i)
         new networkLine(this, points[i], points[i + 1]);

      setNaturalLength();
   }
   double natural_length() const { return total_length / (this->getPoints().size() - 1); };

   //! total_length [m]
   //! weight_per_unit_length [kg/m]  (実際は kg/m の mass_per_length)
   //! natural_length [m]

   void setNaturalLength() {
      for (auto& l : this->getLines()) l->natural_length = natural_length();
   };

   // Phase 0.5: arg renamed `density` → `mass_per_length`.
   // Semantics unchanged: value is stored in `networkLine::weight_per_unit_length`
   // and used as kg/m in mass accumulation below. Callers pass positional args
   // so this rename is source-compatible.
   void setDensityStiffnessDampingDiameter(const double mass_per_length,
                                           const double stiffness,
                                           const double damp,
                                           const double diam) {
      for (auto& l : this->getLines()) {
         l->weight_per_unit_length = mass_per_length;  // kg/m (misnamed historically; name kept for now)
         l->stiffness = stiffness;
         l->damping = damp;
         l->diameter = diam;
      }
      for (auto& p : this->getPoints()) {
         p->mass = 0;
         for (auto& l : p->getLines())
            p->mass += l->natural_length / 2. * l->weight_per_unit_length;
      }
   };

   /* -------------------------------------------------------------------------- */
   /*                                     CFL                                    */
   /* -------------------------------------------------------------------------- */
   //! function to get dt automatically
   double get_dt(const double max_dt) {

      const double C_VELOCITY = 0.01, C_ACCEL = 1;
      double dt_cfl = max_dt, dt_tmp, norm;

      for (const auto& p : this->getPoints())
         for (const auto& l : p->getLines()) {
            //! CFL based on velocity
            norm = (l->stiffness / p->mass) * Norm(Tddd{p->velocity[0], p->velocity[1], p->velocity[2]});
            if (norm != 0.) {
               dt_tmp = C_VELOCITY * l->length() / norm;
               if (dt_cfl > dt_tmp)
                  dt_cfl = dt_tmp;
            }
            // //! CFL based on acceleration
            norm = (l->stiffness / p->mass) * std::max(0.1 * _GRAVITY_, 0.1 * Norm(Tddd{p->acceleration[0], p->acceleration[1], p->acceleration[2]}));
            if (norm != 0) {
               dt_tmp = C_ACCEL * std::sqrt(l->length() / norm);
               if (dt_cfl > dt_tmp)
                  dt_cfl = dt_tmp;
            }
         }

      return dt_cfl;
   }

   // double DragForceCoefficient = 0.3;
   double DragForceCoefficient = 2.5;  // Palm2016はだいたいこのくらい

   // Ambient fluid around the cable (water for mooring, air for bridge stays).
   // Default = water density preserves legacy mooring behavior.
   double FluidDensity = _WATER_DENSITY_;  // [kg/m^3]

   // Wind / ambient-fluid velocity field U(X, t). nullptr = 0 (still fluid).
   // Used by the 3 LumpedCable getDragForce call sites below.
   std::function<Tddd(const Tddd&, double)> wind_field;

   // Rayleigh mass-proportional damping [1/s]. Subtracts alpha_M * mass * v
   // from the net force at each RK4 stage. Default 0 => no effect in dynamic
   // simulations. solveEquilibrium()/setEquilibriumState() install an RAII
   // guard that raises this to ~2 * fundamental frequency so pseudo-relaxation
   // damps the low-frequency pendulum mode exponentially (quadratic drag
   // alone gives only 1/t convergence near v=0). No physical effect outside
   // the guarded region.
   double alpha_M = 0.0;

   /* -------------------------------------------------------------------------- */
   /*         RAII guard that temporarily overrides DragForceCoefficient         */
   /* -------------------------------------------------------------------------- */
   // Used by setEquilibriumState to safely crank Cd up to ~1000 for pseudo
   // relaxation, with automatic restoration on scope exit (including early
   // return / exception). Replaces the manual save/restore that was easy to
   // skip on an error path.
   struct DragForceCoefficientGuard {
      LumpedCable* cable;
      double saved;
      DragForceCoefficientGuard(LumpedCable* c, double temporary_value)
          : cable(c), saved(c->DragForceCoefficient) {
         cable->DragForceCoefficient = temporary_value;
      }
      ~DragForceCoefficientGuard() {
         cable->DragForceCoefficient = saved;
      }
      DragForceCoefficientGuard(const DragForceCoefficientGuard&) = delete;
      DragForceCoefficientGuard& operator=(const DragForceCoefficientGuard&) = delete;
   };

   // Same pattern for FluidDensity — pseudo-relaxation should use water density
   // regardless of target fluid so equilibrium convergence rate is fluid-agnostic.
   struct FluidDensityGuard {
      LumpedCable* cable;
      double saved;
      FluidDensityGuard(LumpedCable* c, double temporary_value)
          : cable(c), saved(c->FluidDensity) {
         cable->FluidDensity = temporary_value;
      }
      ~FluidDensityGuard() { cable->FluidDensity = saved; }
      FluidDensityGuard(const FluidDensityGuard&) = delete;
      FluidDensityGuard& operator=(const FluidDensityGuard&) = delete;
   };

   // Temporarily raise alpha_M (Rayleigh mass-proportional damping) during
   // pseudo-relaxation. Critical-damping the fundamental pendulum mode cuts
   // equilibrium step counts by 10-50x for typical bridge/mooring cables.
   struct AlphaMGuard {
      LumpedCable* cable;
      double saved;
      AlphaMGuard(LumpedCable* c, double temporary_value)
          : cable(c), saved(c->alpha_M) {
         cable->alpha_M = temporary_value;
      }
      ~AlphaMGuard() { cable->alpha_M = saved; }
      AlphaMGuard(const AlphaMGuard&) = delete;
      AlphaMGuard& operator=(const AlphaMGuard&) = delete;
   };

   // Suppress wind during equilibrium: gravity-only relaxation converges
   // to a deterministic catenary; wind-on equilibrium is a user opt-in.
   struct WindFieldGuard {
      LumpedCable* cable;
      std::function<Tddd(const Tddd&, double)> saved;
      WindFieldGuard(LumpedCable* c, std::function<Tddd(const Tddd&, double)> temporary_value)
          : cable(c), saved(std::move(c->wind_field)) {
         cable->wind_field = std::move(temporary_value);
      }
      ~WindFieldGuard() { cable->wind_field = std::move(saved); }
      WindFieldGuard(const WindFieldGuard&) = delete;
      WindFieldGuard& operator=(const WindFieldGuard&) = delete;
   };

   /* -------------------------------------------------------------------------- */
   /*          Pure single RK4 step — no warmup, no stdout, no Cd override       */
   /* -------------------------------------------------------------------------- */
   /*DOC_EXTRACT 2a_LumpedCable_step

   `step` 関数は、指定された `dt` を 1 回の RK4 4 段で一気に進める。
   `simulate` のような dt ランプも、stdout 出力も、`DragForceCoefficient`
   の一時書き換えもしない。BEM 時間ループのように、呼び出し側の RK ステージ
   と同期して `dt` を与えられる用途に最適。

   `RK_velocity_sub`, `RK_X_sub` はいずれも通常の「RK4 substate」として使い、
   関数終了時には **実座標は変更されない**（`applyMooringSimulationResult`
   を呼ぶまで `p->X`, `p->velocity` は元のまま）。次時刻の位置・速度は
   `p->RK_X_sub.get_x()`, `p->RK_velocity_sub.get_x()` で取得可能。

   */
   void step(const double current_time,
             const double dt,
             const std::function<void(networkPoint*)>& setBoundaryCondition) {
      auto points = ToVector(this->getPoints());
      // Initialize RK substates for this single step
      for (auto& p : points) {
         p->RK_velocity_sub.initialize(dt, current_time, p->velocityTranslational(), 4);
         p->RK_X_sub.initialize(dt, current_time, p->X, 4);
      }
      // Run the 4 RK4 stages
      std::array<double, 3> a;
      while (true) {
         for (auto& p : points) {
            Tddd U_wind = this->wind_field ? this->wind_field(p->X, current_time) : Tddd{0., 0., 0.};
            a = (p->getTension() + p->getDragForce(this->DragForceCoefficient, this->FluidDensity, U_wind) + p->getGravitationalForce()) / p->mass;
            if (this->alpha_M != 0.)
               a -= this->alpha_M * p->RK_velocity_sub.get_x();
            std::get<0>(p->acceleration) = std::get<0>(a);
            std::get<1>(p->acceleration) = std::get<1>(a);
            std::get<2>(p->acceleration) = std::get<2>(a);
            setBoundaryCondition(p);
         }
         for (auto& p : points) {
            p->RK_X_sub.push(p->RK_velocity_sub.get_x());
            p->RK_velocity_sub.push(p->accelTranslational());
         }
         if ((*points.begin())->RK_X_sub.finished)
            break;
      }
   }

   /* -------------------------------------------------------------------------- */
   /*DOC_EXTRACT 2_LumpedCable

   `simulate` 関数は，`netwrokPoint`が持つ`getForce`関数を用いて力を計算する．
   内部では dt を `get_dt` と step カウントに基づいて 1e-8 〜 1e-3 dt の
   マイクロステップにランプアップしながら積み上げて進む。`setEquilibriumState`
   専用の warm-up スキーム。

   Phase 0.5: `silent=true` を渡せば進捗 stdout を抑制できる。引数順は
   既存コード互換（末尾に default 付きで追加）。

   */
   /* -------------------------------------------------------------------------- */

   void simulate(const double current_time,
                 const double dt,
                 const std::function<void(networkPoint*)>& setBoundaryCondition,
                 const bool silent = false) {

      double dt_acum = 0;
      bool first = true;

      //! 実際の値は変更されない
      int i = 0, step = 0;
      auto points = ToVector(this->getPoints());
      while (dt_acum < dt) {
         double dt_cfl = get_dt(dt);
         if (step < 10)
            dt_cfl = dt * 1E-8;
         else if (step < 100)
            dt_cfl = dt * 1E-7;
         else if (step < 1000)
            dt_cfl = dt * 1E-6;
         else
            dt_cfl = std::clamp(dt_cfl, dt * 1E-4, dt * 1E-3);

         // std::cout << Red << "dt_cfl = " << dt_cfl << colorReset << std::endl;
         if (dt_acum + dt_cfl >= dt)
            dt_cfl = dt - dt_acum;
         /* -------------------------------------------------------------------------- */
         //! initialize RK
         for (auto& p : points) {
            p->RK_velocity_sub.initialize(dt_cfl, current_time, first ? p->velocityTranslational() : p->RK_velocity_sub.get_x(), 4);
            p->RK_X_sub.initialize(dt_cfl, current_time, first ? p->X : p->RK_X_sub.get_x(), 4);
         }
         first = false;
         /* -------------------------------------------------------------------------- */
         // Print("simulate ", dt_cfl, ",", dt_acum, ",", dt);
         std::array<double, 3> a;

         while (1) {
            // std::cout << "RK : " << (*this->getPoints().begin())->RK_X_sub.current_step << std::endl;
            for (auto& p : points) {
               Tddd U_wind = this->wind_field ? this->wind_field(p->X, current_time) : Tddd{0., 0., 0.};
               a = (p->getTension() + p->getDragForce(this->DragForceCoefficient, this->FluidDensity, U_wind) + p->getGravitationalForce()) / p->mass;
               if (this->alpha_M != 0.)
                  a -= this->alpha_M * p->RK_velocity_sub.get_x();
               std::get<0>(p->acceleration) = std::get<0>(a);  // accelは変更しても構わない
               std::get<1>(p->acceleration) = std::get<1>(a);  // accelは変更しても構わない
               std::get<2>(p->acceleration) = std::get<2>(a);  // accelは変更しても構わない
               setBoundaryCondition(p);
            }

            for (auto& p : points) {
               p->RK_X_sub.push(p->RK_velocity_sub.get_x());
               p->RK_velocity_sub.push(p->accelTranslational());
            }

            if ((*points.begin())->RK_X_sub.finished)
               break;
         }

         dt_acum += dt_cfl;
         if (dt_acum >= dt)
            return;

         double norm_velcoity = 0, norm_acceleration = 0;
         for (auto& p : points) {
            norm_velcoity += Norm(p->velocity);
            norm_acceleration += Norm(p->acceleration);
         }
         // show percentage (suppressed when silent)
         auto percentage = (dt_acum / dt) * 100.;
         if (!silent && percentage > 10 * i) {
            std::cout << "percentage = " << percentage << std::endl;
            std::cout << "dt = " << dt << ", dt_cfl = " << dt_cfl << ", time = " << dt_acum << ", norm_velcoity = " << norm_velcoity << ", norm_acceleration = " << norm_acceleration << std::endl;
            i++;
         }
         if (step++ > 10000000) {
            std::stringstream ss;
            ss << "step > " << step;
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
         }
      }
   }

   /* -------------------------------------------------------------------------- */
   /*DOC_EXTRACT 3_LumpedCable_setEquilibriumState

   静的平衡解を疑似緩和（dynamic relaxation）で求める。`DragForceCoefficient`
   を一時的に 1000 に引き上げ（`DragForceCoefficientGuard` により自動復元）、
   `simulate()` を反復して速度ノルムが `tol` を下回るまで走る。

   Phase 0.5 (2026-04-12) で修正した点:
   1. 収束判定 break を**復活**（旧版はコメントアウトで常に 100 回固定ループ）
   2. `tol` と `max_iters` を引数化（デフォルトは旧互換の 1e-3, 100）
   3. 収束可否を bool 返値で返す
   4. `DragForceCoefficient` 差し替えを RAII guard に（例外時も自動復元）
   5. 内部 `simulate()` 呼び出しを `silent=true` で行う

   */
   /* -------------------------------------------------------------------------- */
   bool setEquilibriumState(const std::function<void(networkPoint*)>& setBoundaryCondition,
                            const double tol = 1e-3,
                            const int max_iters = 100) {
      // Only apply pseudo-relaxation overrides when no wind is set — with
      // wind, the steady drag shifts the physical fixed point and the
      // overrides would push the cable to the wrong equilibrium.
      std::unique_ptr<DragForceCoefficientGuard> cd_guard;
      std::unique_ptr<FluidDensityGuard> rho_guard;
      std::unique_ptr<AlphaMGuard> alpha_guard;
      if (!this->wind_field) {
         cd_guard = std::make_unique<DragForceCoefficientGuard>(this, 1000.);
         rho_guard = std::make_unique<FluidDensityGuard>(this, _WATER_DENSITY_);
         // Rayleigh mass-proportional damping ~ 2x fundamental pendulum
         // frequency. Rough string estimate (T unknown here): use cable
         // length alone as omega ~ pi * sqrt(g/L). Adequate for legacy path.
         const double L = this->natural_length();
         const double omega = (M_PI) * std::sqrt(std::max(_GRAVITY_ / L, 1e-12));
         alpha_guard = std::make_unique<AlphaMGuard>(this, 2.0 * omega);
      }

      auto points = ToVector(this->getPoints());
      double n = static_cast<double>(points.size());
      double norm_total_velocity = 0;
      bool converged = false;

      for (int i = 0; i < max_iters; ++i) {
         //% 内部で計算回数を制限しているのでこのようになる
         this->simulate(0, 1., setBoundaryCondition, /*silent=*/true);
         norm_total_velocity = 0;
         for (const auto& p : points)
            norm_total_velocity += Norm(p->RK_velocity_sub.get_x());
         norm_total_velocity /= n;
         std::cout << "setEquilibriumState iter " << i
                   << "  norm_total_velocity = " << norm_total_velocity << std::endl;
         applyMooringSimulationResult();

         if (norm_total_velocity < tol) {
            converged = true;
            break;
         }
      }
      // DragForceCoefficientGuard destructor restores original Cd here.
      return converged;
   }

   void applyMooringSimulationResult() {
      for (auto& p : this->getPoints()) {
         p->setX(p->RK_X_sub.get_x());
         auto v = p->RK_velocity_sub.get_x();
         std::get<0>(p->velocity) = std::get<0>(v);
         std::get<1>(p->velocity) = std::get<1>(v);
         std::get<2>(p->velocity) = std::get<2>(v);
      }
   }
};

/* ========================================================================== */
/*                        Cable system (Phase 1, 2026-04-12)                  */
/* ========================================================================== */

/*DOC_EXTRACT 4_LumpedCableSystem

## `LumpedCableSystem` — 複数ケーブルの集約とライフサイクル管理

橋梁の 12 本 stay、海洋係留の 3-脚カテナリ、BEM 浮体連成などで「1 個の構造
物に取り付いた複数本のケーブル」を 1 つの単位として扱う。

`CableAttachment` は各端点が**世界固定**（WorldFixed — 橋梁塔・海底アンカー）
か**浮体取付**（BodyFrame — フェアリード）のどちらかを保持する。浮体取付の
場合は、その `Network*` 浮体の剛体姿勢から現在の world 位置を動的に計算する。
これにより、`setEquilibriumState` ループと `step(t, dt, BC)` ループの両方で、
端点が**浮体と共に動く**挙動が得られる。

使用例:

```cpp
// 橋梁ケーブル (両端 WorldFixed、BEM なしで独立計算):
LumpedCableSystem sys;
for (const auto& cable_spec : parsed_cables) {
    sys.addCable(cable_spec.name,
                 CableAttachment::worldFixed(cable_spec.point_a),
                 CableAttachment::worldFixed(cable_spec.point_b),
                 cable_spec.cable_length,
                 cable_spec.n_segments,
                 {cable_spec.mass_per_length, cable_spec.EA,
                  cable_spec.damping, cable_spec.diameter});
}
sys.solveEquilibrium();

// BEM 浮体係留 (アンカー world 固定、フェアリード body 追従):
floater->cable_system = std::make_unique<LumpedCableSystem>();
CableAttachment anchor = CableAttachment::worldFixed({500, 0, -58});
CableAttachment fairlead = CableAttachment::onBody(floater, {0, 0, 0});
floater->cable_system->addCable("lineA", anchor, fairlead, 522., 100, props);
// ... 時間ループ内 ...
floater->cable_system->advanceRKStage(t, dt);
if (floater->RK_Q.finished)
    floater->cable_system->commitRKStep();
auto [F, T] = floater->cable_system->forceOnBody(floater);
```

*/

/* -------------------------------------------------------------------------- */
/*                                POD: Properties                             */
/* -------------------------------------------------------------------------- */
struct CableProperties {
   double mass_per_length = 0.0;  // [kg/m] 単位長あたり質量
   double EA = 0.0;               // [N]    軸剛性
   double damping = 0.0;          // [N·s/m] networkLine::damping 相当
   double diameter = 0.0;         // [m]    抗力射影面積用
   double EI = 0.0;               // [N·m²] 曲げ剛性 (将来用、現状未使用)
};

/* -------------------------------------------------------------------------- */
/*                                POD: Attachment                             */
/* -------------------------------------------------------------------------- */
/*
 * Represents one end of a cable. Either world-fixed (e.g. bridge tower, seabed
 * anchor) or attached to a body (e.g. floater fairlead). The body-attached
 * case uses the cable-endpoint node's own initial position as the body-frame
 * anchor for rigid transformation — this mirrors the existing
 * `nextPositionOnBody()` logic from BEM_calculateVelocities.hpp.
 */
struct CableAttachment {
   enum Kind { WorldFixed, BodyFrame };
   Kind kind = WorldFixed;

   // For WorldFixed: the fixed world position.
   // For BodyFrame: the INITIAL world position (used by addCable to seed
   // the straight-line initializer and to derive body_offset at t=0).
   std::array<double, 3> world_position{0., 0., 0.};

   // For BodyFrame: the attached body. The cable endpoint node itself is
   // stored in `node` by LumpedCableSystem::addCable — it is the
   // networkPoint whose `initialX` and `RK_X` drive the rigid transform.
   Network* body = nullptr;

   // Filled in by LumpedCableSystem::addCable once the cable is built.
   networkPoint* node = nullptr;

   // -------- Factory helpers --------
   static CableAttachment worldFixed(const std::array<double, 3>& pos) {
      CableAttachment a;
      a.kind = WorldFixed;
      a.world_position = pos;
      return a;
   }
   static CableAttachment onBody(Network* b,
                                  const std::array<double, 3>& initial_world_pos) {
      CableAttachment a;
      a.kind = BodyFrame;
      a.body = b;
      a.world_position = initial_world_pos;
      return a;
   }

   // -------- Current world position (reflecting body's current RK pose) --------
   // Mirrors `nextPositionOnBody` from bem/core/BEM_calculateVelocities.hpp:287
   // with all three branches (SoftBody, relative_velocity, RigidBody) + the
   // isFixed axis-lock handling. For WorldFixed, just returns `world_position`.
   std::array<double, 3> currentWorldPosition() const {
      if (kind == WorldFixed) return world_position;
      if (!body || !node) return world_position;  // defensive fallback

      // Branch A: SoftBody or "relative_velocity" JSON flag
      if (body->isSoftBody || body->inputJSON.find("relative_velocity")) {
         auto v = node->velocityTranslational();
         for (int i = 0; i < 3; ++i)
            if (body->isFixed[i]) v[i] = 0;
         if (body->isFixed.size() == 1 && body->isFixed[0])
            v.fill(0.0);
         return node->RK_X.getX(v);
      }

      // Branch B: RigidBody — rigid-transform the node's initialX via the
      // body's current predicted COM and quaternion.
      if (body->isRigidBody) {
         auto velocity = body->velocityTranslational();
         auto rotation = body->velocityRotational();
         for (int i = 0; i < 3; ++i)
            if (body->isFixed[i]) velocity[i] = 0;
         for (int i = 3; i < 6; ++i)
            if (body->isFixed[i]) rotation[i - 3] = 0;
         auto COM_next = body->RK_COM.getX(velocity);
         auto Q_next = Quaternion(body->RK_Q.getX(
             AngularVelocityTodQdt(rotation, body->RK_Q.getX())));
         return rigidTransformation(body->ICOM, COM_next, Q_next.Rv(), node->initialX);
      }

      // Branch C: Fallback — node's current (actual) position
      return node->X;
   }
};

/* -------------------------------------------------------------------------- */
/*                              LumpedCableSystem                             */
/* -------------------------------------------------------------------------- */
class LumpedCableSystem {
  public:
   // SNAPSHOT callback signature: (iter, max_vel, positions_by_cable_name)
   using SnapshotCallback = std::function<void(
       int /*iter*/,
       double /*max_vel*/,
       const std::map<std::string, std::vector<std::array<double, 3>>>& /*positions*/)>;

   ~LumpedCableSystem() {
      for (auto* c : _cables_view) delete c;
   }

   // ------------------------------------------------------------------
   // Build
   // ------------------------------------------------------------------
   LumpedCable* addCable(const std::string& name,
                         CableAttachment end_a,
                         CableAttachment end_b,
                         double natural_length,
                         int n_points,
                         const CableProperties& props) {
      // Straight-line initial positions from each attachment's world_position
      const auto X0 = end_a.world_position;
      const auto X1 = end_b.world_position;

      auto* cable = new LumpedCable(X0, X1, natural_length, n_points);
      cable->setName(name);
      cable->setDensityStiffnessDampingDiameter(
          props.mass_per_length, props.EA, props.damping, props.diameter);

      // Wire the attachment's node pointer to the actual cable endpoint.
      end_a.node = cable->firstPoint;
      end_b.node = cable->lastPoint;

      _cables_view.push_back(cable);
      _names.push_back(name);
      _end_a.push_back(end_a);
      _end_b.push_back(end_b);
      _props.push_back(props);
      return cable;
   }

   // ------------------------------------------------------------------
   // Query
   // ------------------------------------------------------------------
   const std::vector<LumpedCable*>& cables() const { return _cables_view; }
   size_t size() const { return _cables_view.size(); }

   LumpedCable* cableByName(const std::string& name) const {
      for (size_t i = 0; i < _names.size(); ++i)
         if (_names[i] == name) return _cables_view[i];
      return nullptr;
   }
   const std::string& nameOf(size_t i) const { return _names.at(i); }
   const CableAttachment& endA(size_t i) const { return _end_a.at(i); }
   const CableAttachment& endB(size_t i) const { return _end_b.at(i); }

   // ------------------------------------------------------------------
   // Static equilibrium (bridge-cable / floater initialisation use case)
   // ------------------------------------------------------------------
   /*
    * Ported from cable/cable_solver.cpp:237-319 — the self-owned RK4
    * equilibrium driver that bypasses LumpedCable::simulate()'s warmup
    * ramp. Uses a fixed CFL=1.0 dt based on the axial wave speed, inflates
    * DragForceCoefficient to 1000 (via RAII) for pseudo-relaxation, and
    * stops as soon as (step > 1000 && max_vel < tol) across all cables.
    * Returns true on convergence, false if max_steps reached.
    */
   bool solveEquilibrium(double tol = 0.01,
                         int max_steps = 500000,
                         int snapshot_interval = 10000,
                         SnapshotCallback snapshot_cb = {}) {
      if (_cables_view.empty()) return true;

      // 1. RAII guards for pseudo-relaxation. Two regimes:
      //   (a) No wind on any cable: inflate Cd=1000 and use water density for
      //       velocity-proportional damping so relaxation is fast and
      //       fluid-agnostic. Equilibrium is drag-free (v=0 => drag=0) so
      //       these overrides do not shift the final equilibrium.
      //   (b) Any cable has wind_field != nullptr: the physical equilibrium
      //       includes a steady drag term 0.5*rho*Cd*A*|U_wind|^2 that shifts
      //       the fixed point. Applying the Cd/rho override here would move
      //       the cable to the wrong fixed point. Keep user's settings
      //       unchanged — convergence may be slower but correct.
      bool any_wind = false;
      for (auto* c : _cables_view) {
         if (c->wind_field) { any_wind = true; break; }
      }
      std::vector<std::unique_ptr<LumpedCable::DragForceCoefficientGuard>> cd_guards;
      std::vector<std::unique_ptr<LumpedCable::FluidDensityGuard>> rho_guards;
      std::vector<std::unique_ptr<LumpedCable::AlphaMGuard>> alpha_guards;
      if (!any_wind) {
         cd_guards.reserve(_cables_view.size());
         rho_guards.reserve(_cables_view.size());
         alpha_guards.reserve(_cables_view.size());
         for (size_t i = 0; i < _cables_view.size(); ++i) {
            auto* c = _cables_view[i];
            cd_guards.emplace_back(std::make_unique<LumpedCable::DragForceCoefficientGuard>(c, 1000.0));
            rho_guards.emplace_back(std::make_unique<LumpedCable::FluidDensityGuard>(c, _WATER_DENSITY_));
            // Rayleigh alpha_M: critical-damp the SLOWEST mode, which
            // dominates convergence time. For a suspended cable, this is
            // typically the pendulum/sway mode omega_p ~ sqrt(g/L) (not the
            // stretched-string fundamental omega_s = pi/L * sqrt(T/mu), which
            // is usually much higher). Over-damping omega_p with large
            // alpha_M would slow convergence, not speed it up. Higher modes
            // remain under-damped here but are damped fast by the inflated
            // quadratic drag (Cd=1000).
            const double L = c->natural_length();
            const double omega_p = std::sqrt(std::max(_GRAVITY_ / L, 1e-12));
            const double alpha = 2.0 * omega_p;  // ~critical for pendulum mode
            alpha_guards.emplace_back(std::make_unique<LumpedCable::AlphaMGuard>(c, alpha));
         }
      }

      // 2. Per-cable CFL-based dt (NOT shared — each cable uses its own
      //    native dt = natural_length / wave_speed at CFL=1.0). Sharing a
      //    single min-dt across all cables causes longer cables to run with
      //    sub-optimal dt and converge to a slightly different fixed point
      //    at the same step count, producing 0.2% tension differences in
      //    multi-cable mode vs the per-cable single-cable runs. Per-cable
      //    dt restores byte-identity with the legacy cable_solver.cpp
      //    self-owned RK4 driver for each individual cable.
      std::vector<double> dt_cfl_per_cable(_cables_view.size(), 0.0);
      for (size_t i = 0; i < _cables_view.size(); ++i) {
         auto* c = _cables_view[i];
         double wave_speed = std::sqrt(_props[i].EA / _props[i].mass_per_length);
         dt_cfl_per_cable[i] = 1.0 * c->natural_length() / wave_speed;  // CFL=1.0
      }

      // 3. Per-cable ordered point list (firstPoint → lastPoint traversal).
      std::vector<std::vector<networkPoint*>> ordered_per_cable;
      ordered_per_cable.reserve(_cables_view.size());
      for (auto* c : _cables_view)
         ordered_per_cable.push_back(orderedPointsOf(c));

      // 5. Main RK4 relaxation loop.
      // Per-cable convergence lock: once a cable's max_vel drops below tol
      // (after the warm-up grace period of 1000 steps), it stops being
      // integrated. This restores byte-identical equivalence with the
      // legacy single-cable cable_solver.cpp driver — each cable freezes at
      // exactly the iteration it would have terminated at in isolation.
      // The outer loop ends once every cable is converged.
      std::vector<bool> converged_per_cable(_cables_view.size(), false);
      // Per-cable max velocity from the last RK step (used for snapshot
      // emission and the convergence check). Reset each step.
      std::vector<double> max_vel_per_cable(_cables_view.size(), 0.0);

      bool all_converged = false;
      for (int step = 0; step < max_steps; ++step) {
         double max_vel_all = 0.0;

         // --- For each not-yet-converged cable, run one RK4 step ---
         for (size_t ci = 0; ci < _cables_view.size(); ++ci) {
            if (converged_per_cable[ci]) {
               // Locked: contribute its frozen max_vel (already below tol)
               // to the global max but skip integration entirely.
               if (max_vel_per_cable[ci] > max_vel_all)
                  max_vel_all = max_vel_per_cable[ci];
               continue;
            }
            auto& ordered = ordered_per_cable[ci];
            auto* cable = _cables_view[ci];
            const double dt_cfl = dt_cfl_per_cable[ci];

            // Initialize RK4 for this step
            for (auto& p : ordered) {
               p->RK_velocity_sub.initialize(dt_cfl, 0, p->velocityTranslational(), 4);
               p->RK_X_sub.initialize(dt_cfl, 0, p->X, 4);
            }

            // 4 RK4 stages
            while (true) {
               for (auto& p : ordered) {
                  // During equilibrium wind_field is RAII-guarded to nullptr,
                  // so U_wind is always zero here. Kept in the call for symmetry.
                  Tddd U_wind = cable->wind_field ? cable->wind_field(p->X, 0.) : Tddd{0., 0., 0.};
                  auto a = (p->getTension()
                            + p->getDragForce(cable->DragForceCoefficient, cable->FluidDensity, U_wind)
                            + p->getGravitationalForce()) / p->mass;
                  if (cable->alpha_M != 0.)
                     a -= cable->alpha_M * p->RK_velocity_sub.get_x();
                  std::get<0>(p->acceleration) = std::get<0>(a);
                  std::get<1>(p->acceleration) = std::get<1>(a);
                  std::get<2>(p->acceleration) = std::get<2>(a);
                  // BC: pin firstPoint and lastPoint
                  if (p == cable->firstPoint || p == cable->lastPoint) {
                     p->acceleration.fill(0);
                     p->velocity.fill(0);
                  }
               }
               for (auto& p : ordered) {
                  p->RK_X_sub.push(p->RK_velocity_sub.get_x());
                  p->RK_velocity_sub.push(p->accelTranslational());
               }
               if (ordered[0]->RK_X_sub.finished) break;
            }

            // Commit RK4 result to actual positions/velocities for this cable
            double max_vel_this = 0.0;
            for (auto& p : ordered) {
               p->setX(p->RK_X_sub.get_x());
               auto v = p->RK_velocity_sub.get_x();
               std::get<0>(p->velocity) = std::get<0>(v);
               std::get<1>(p->velocity) = std::get<1>(v);
               std::get<2>(p->velocity) = std::get<2>(v);
               double vnorm = Norm(v);
               if (vnorm > max_vel_this) max_vel_this = vnorm;
            }
            max_vel_per_cable[ci] = max_vel_this;
            if (max_vel_this > max_vel_all) max_vel_all = max_vel_this;

            // Per-cable lock: same trigger as legacy single-line driver
            // (`step > 1000 && max_vel < tol`). Once locked, the cable
            // keeps its current state for the rest of the run.
            if (step > 1000 && max_vel_this < tol) {
               converged_per_cable[ci] = true;
            }
         }

         // --- Snapshot emission ---
         bool is_snapshot = (step % snapshot_interval == 0) || (step > 1000 && max_vel_all < tol);
         if (is_snapshot && snapshot_cb) {
            std::map<std::string, std::vector<std::array<double, 3>>> positions_by_cable;
            for (size_t ci = 0; ci < _cables_view.size(); ++ci) {
               auto& ordered = ordered_per_cable[ci];
               std::vector<std::array<double, 3>> pos;
               pos.reserve(ordered.size());
               for (auto* p : ordered) pos.push_back(p->X);
               positions_by_cable[_names[ci]] = std::move(pos);
            }
            snapshot_cb(step, max_vel_all, positions_by_cable);
         }

         // --- Global convergence check ---
         all_converged = std::all_of(converged_per_cable.begin(),
                                      converged_per_cable.end(),
                                      [](bool b) { return b; });
         if (all_converged) break;
      }
      // DragForceCoefficientGuard destructors restore the original Cd here.
      return all_converged;
   }

   // ------------------------------------------------------------------
   // BEM time-loop stepping (two-phase API)
   // ------------------------------------------------------------------
   /*
    * BEM calls advanceRKStage() at each RK stage (passing the outer stage dt),
    * and commitRKStep() exactly once when the body's RK accumulator reports
    * `.finished`. Replaces the inlined mooring block in
    * main_time_domain.cpp:904-926.
    */
   void advanceRKStage(double current_time, double dt) {
      if (_cables_view.empty() || dt <= 0.0) return;

      for (size_t ci = 0; ci < _cables_view.size(); ++ci) {
         auto* cable = _cables_view[ci];
         const auto& a = _end_a[ci];
         const auto& b = _end_b[ci];

         // Target world positions for each endpoint at the next RK sub-state.
         auto target_a = a.currentWorldPosition();
         auto target_b = b.currentWorldPosition();

         // Current (pre-step) positions
         const auto cur_a = cable->firstPoint->X;
         const auto cur_b = cable->lastPoint->X;

         // Velocity BC — for WorldFixed it's zero (target == current, for a
         // static endpoint), for BodyFrame it's (target - current) / dt.
         std::array<double, 3> v_a = {(target_a[0] - cur_a[0]) / dt,
                                       (target_a[1] - cur_a[1]) / dt,
                                       (target_a[2] - cur_a[2]) / dt};
         std::array<double, 3> v_b = {(target_b[0] - cur_b[0]) / dt,
                                       (target_b[1] - cur_b[1]) / dt,
                                       (target_b[2] - cur_b[2]) / dt};

         auto bc = [&](networkPoint* p) {
            if (p == cable->firstPoint) {
               p->acceleration.fill(0);
               p->velocity[0] = v_a[0];
               p->velocity[1] = v_a[1];
               p->velocity[2] = v_a[2];
            } else if (p == cable->lastPoint) {
               p->acceleration.fill(0);
               p->velocity[0] = v_b[0];
               p->velocity[1] = v_b[1];
               p->velocity[2] = v_b[2];
            }
         };

         cable->step(current_time, dt, bc);
      }
   }

   void commitRKStep() {
      for (size_t ci = 0; ci < _cables_view.size(); ++ci) {
         auto* cable = _cables_view[ci];
         cable->applyMooringSimulationResult();
         // Snap the pinned endpoints back to the exact attachment position
         // to undo any drift from the velocity-based BC approximation.
         cable->firstPoint->setX(_end_a[ci].currentWorldPosition());
         cable->lastPoint->setX(_end_b[ci].currentWorldPosition());
      }
   }

   // ------------------------------------------------------------------
   // Force feedback to a floating body
   // ------------------------------------------------------------------
   /*
    * Returns the total force and moment that all fairleads attached to
    * `body` exert on the body. Moment is taken about body->COM.
    *
    * Handles the symmetric case where either end_a or end_b may be the
    * body-attached end. Sign matches the existing BEM behaviour in
    * BEM_solveBVP.hpp:1324-1340 (cable-side `getForce()` summed directly).
    */
   std::pair<std::array<double, 3>, std::array<double, 3>>
   forceOnBody(const Network* body) const {
      std::array<double, 3> F{0., 0., 0.};
      std::array<double, 3> T{0., 0., 0.};

      for (size_t ci = 0; ci < _cables_view.size(); ++ci) {
         networkPoint* fairlead = nullptr;
         if (_end_a[ci].kind == CableAttachment::BodyFrame && _end_a[ci].body == body)
            fairlead = _cables_view[ci]->firstPoint;
         else if (_end_b[ci].kind == CableAttachment::BodyFrame && _end_b[ci].body == body)
            fairlead = _cables_view[ci]->lastPoint;
         if (!fairlead) continue;

         auto f = fairlead->getForce();
         F = F + f;
         T = T + Cross(fairlead->X - body->COM, f);
      }
      return {F, T};
   }

  private:
   // Walk the cable from firstPoint to lastPoint in topological order.
   // Duplicate of cable_solver.cpp::getOrderedPoints, kept private here.
   static std::vector<networkPoint*> orderedPointsOf(const LumpedCable* cable) {
      std::vector<networkPoint*> ordered;
      auto current = cable->firstPoint;
      networkPoint* prev = nullptr;
      while (current) {
         ordered.push_back(current);
         if (current == cable->lastPoint) break;
         networkPoint* next = nullptr;
         for (auto l : current->getLines()) {
            auto other = (*l)(current);
            if (other != prev) {
               next = other;
               break;
            }
         }
         prev = current;
         current = next;
      }
      return ordered;
   }

   std::vector<LumpedCable*> _cables_view;
   std::vector<std::string> _names;
   std::vector<CableAttachment> _end_a;
   std::vector<CableAttachment> _end_b;
   std::vector<CableProperties> _props;
};

/* -------------------------------------------------------------------------- */

// -----------------------------------------------------------------------------
// Backward-compat alias — kept in this header (not just in the MooringLine.hpp
// shim) so that any translation unit including Network.hpp (which now pulls in
// LumpedCable.hpp) transparently sees `MooringLine` as a synonym. This lets
// the staged migration touch BEM / cable call sites incrementally. Removal
// plan: once all call sites have been updated to `LumpedCable`, delete both
// this alias and the MooringLine.hpp shim.
// -----------------------------------------------------------------------------
using MooringLine = LumpedCable;
