#pragma once

#include "basic.hpp"
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

/*DOC_EXTRACT ODE::Runge-Kutta

### Runge-Kutta

4次のルンゲクッタの１回の計算で溜まる誤差は$O({\Delta t}^5)$となる．
しかし，加速度を4階も計算する必要がある．
このように，ルンゲクッタを使って２階微分方程式を解く場合，
２階微分方程式を２つの1階微分方程式にわけて考え，互いに独立した２つのルンゲクッタを用意し，それぞれ現時刻の微分を使って更新する．
後退オイラーのように次時刻の流速を使って位置を更新するということはできない．

\ref{ODE:RungeKutta4}{4次のRunge-Kutta}の場合，次のようになる．

$$
\begin{align*}
k_1 &= \frac{dx}{dt}(t_n, x_n)\\
k_2 &= \frac{dx}{dt}(t_n + \frac{\Delta t}{2}, x_n + \frac{\Delta t}{2} k_1)\\
k_3 &= \frac{dx}{dt}(t_n + \frac{\Delta t}{2}, x_n + \frac{\Delta t}{2} k_2)\\
k_4 &= \frac{dx}{dt}(t_n + \Delta t, x_n + \Delta t k_3)\\
x_{n+1} &= x_n + \frac{\Delta t}{6} (k_1 + 2 k_2 + 2 k_3 + k_4)
\end{align*}
$$

\ref{ODE:RungeKutta}{RungeKuttaのクラス}

*/

// \label{ODE:RungeKutta}
template <typename T> struct RungeKuttaCommon {
  double dt_fixed = 0.0; // 基本となるステップ間隔,固定
  /* ----------------------- 毎回変わる値 ----------------------- */
  double dt = 0.0;    // 次回計算する必要がある値になる増分,回数ごとに変化
  std::vector<T> _dX; // これにプッシュしていく，c1k1,c2k2,c3k3,c4k4,,,//このベクトルの長さで何回目かを判断数し，必要な回数まできたら，finished=trueとなる
  T dX;               //_dXやdtを基に計算される
  /* ------------------------- 初期値 ------------------------ */
  double t_init = 0.0; // 初期値,固定
  T Xinit;
  /* ------------------------ 積分結果 ------------------------ */
  bool finished = false; // 全てそろって積分が評価できるようになったらtrue
                         // 最終結果はdXに代入される;
  void displayStatus() {
    std::cout << "dt : " << red << this->dt << colorReset << std::endl;
    // std::cout << "_dX : " << red << this->_dX << colorReset << std::endl;
    std::cout << "dX : " << red << this->dX << colorReset << std::endl;
  };
  int steps = 0;
  int current_step = 0;
  RungeKuttaCommon() : Xinit{}, dX{} {};
  RungeKuttaCommon(const double dt_IN, const double t0, const T &X0, int stepsIN) : dt_fixed(dt_IN), dt(0), t_init(t0), Xinit(X0), _dX({}), dX(X0), steps(stepsIN), current_step(0) {
    finished = false;
    _dX.reserve(stepsIN);
  };
  ~RungeKuttaCommon() {};

  void initialize(const double dt_IN, const double t0, const T &X0, int stepsIN) {
    this->dt_fixed = dt_IN;
    this->dt = 0;
    this->t_init = t0;
    this->Xinit = X0;
    this->_dX = {};
    this->dX = X0;
    this->steps = stepsIN;
    if (stepsIN < 1 || stepsIN > 4)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "stepsINが1~4の整数でなければならない");
    this->_dX.reserve(stepsIN);
    this->current_step = 0;
    this->finished = false;
  };

  T getXinit() const { return this->Xinit; };

  T getX() const {
    if (this->steps == 0)
      return this->Xinit;
    if (this->current_step == 0)
      return this->Xinit;
    else
      return this->Xinit + this->dX;
  };
  T get_x() const { return getX(); };

  T getdX() const { return this->dX; };
  double gett() const { return this->t_init + this->dt; };
  double get_t() const { return gett(); };

  double getTimeAtCurrentStep() const {
    switch (this->steps) {
    case 1:
      return this->t_init;

    case 2:
      switch (current_step) {
      case 0:
        return this->t_init;
      case 1:
        return this->t_init + dt_fixed; // Heun: 2nd stage at t + dt
      default:
        return this->t_init + dt_fixed;
      }

    case 3:
      switch (current_step) {
      case 0:
        return this->t_init;
      case 1:
        return this->t_init + dt_fixed / 2.0;
      case 2:
        return this->t_init + dt_fixed;
      default:
        return this->t_init + dt_fixed;
      }

    case 4:
      switch (current_step) {
      case 0:
        return this->t_init;
      case 1:
        return this->t_init + dt_fixed / 2.0;
      case 2:
        return this->t_init + dt_fixed / 2.0;
      case 3:
        return this->t_init + dt_fixed;
      case 4:
        return this->t_init + dt_fixed;
      default:
        return this->t_init + dt_fixed;
      }

    default:
      std::stringstream ss;
      ss << std::to_string(this->current_step) << "/" << std::to_string(this->steps);
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "current_step/steps=" + ss.str());
    }
  }

  double getTimeAtNextStep() const {
    switch (this->steps) {
    case 1:
      return this->t_init + dt_fixed;

    case 2:
      switch (current_step) {
      case 0:
        return this->t_init + dt_fixed; // Heun: next evaluation at t + dt
      default:
        return this->t_init + dt_fixed;
      }

    case 3:
      switch (current_step) {
      case 0:
        return this->t_init + dt_fixed / 2.0;
      case 1:
        return this->t_init + dt_fixed;
      default:
        return this->t_init + dt_fixed;
      }

    case 4:
      switch (current_step) {
      case 0:
        return this->t_init + dt_fixed / 2.0;
      case 1:
        return this->t_init + dt_fixed / 2.0;
      case 2:
        return this->t_init + dt_fixed;
      case 3:
        return this->t_init + dt_fixed;
      default:
        return this->t_init + dt_fixed;
      }

    default:
      std::stringstream ss;
      ss << std::to_string(this->current_step) << "/" << std::to_string(this->steps);
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "current_step/steps=" + ss.str());
    }
  }

  double getdt() const {
    switch (this->steps) {
    case 1:
      return dt_fixed;

    case 2:
      switch (current_step) {
      case 0:
        return dt_fixed; // Heun uses full dt between stages in this implementation
      default:
        return dt_fixed;
      }

    case 3:
      switch (current_step) {
      case 0:
        return dt_fixed / 2.0;
      case 1:
        return dt_fixed;
      default:
        return dt_fixed;
      }

    case 4:
      switch (current_step) {
      case 0:
        return dt_fixed / 2.0;
      case 1:
        return dt_fixed / 2.0;
      case 2:
        return dt_fixed;
      case 3:
        return dt_fixed;
      default:
        return dt_fixed;
      }

    default:
      std::stringstream ss;
      ss << std::to_string(this->current_step) << "/" << std::to_string(this->steps);
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "current_step/steps=" + ss.str());
    }
  }

  bool repush(const T &dXdt_IN) {
    if (current_step != 0) {
      _dX.pop_back();
      current_step--;
    }
    return this->push(dXdt_IN);
  };

  bool push(const T &dXdt_IN) {
    _dX.push_back(dXdt_IN * dt_fixed);

    switch (this->steps) {
    case 1:
      switch (current_step++) {
      case 0:
        dX = dXdt_IN * (dt = dt_fixed);
        return this->finished = true;
      default:
        break;
      }
      break;

    case 2:
      switch (current_step++) {
      case 0:
        dX = dXdt_IN * (dt = dt_fixed);
        return this->finished = false;
      case 1:
        dX = (_dX[0] + _dX[1]) / 2.0;
        dt = dt_fixed;
        return this->finished = true;
      default:
        break;
      }
      break;

    case 3:
      switch (current_step++) {
      case 0:
        dX = dXdt_IN * (dt = dt_fixed / 2.0);
        return finished = false;
      case 1:
        dX = dXdt_IN * (dt = dt_fixed);
        return finished = false;
      case 2:
        dX = (_dX[0] + 4.0 * _dX[1] + _dX[2]) / 6.0;
        dt = dt_fixed;
        return finished = true;
      default:
        break;
      }
      break;

    case 4:
      switch (current_step++) {
      case 0:
      case 1:
        dX = dXdt_IN * (dt = dt_fixed / 2.0);
        return finished = false;
      case 2:
        dX = dXdt_IN * (dt = dt_fixed);
        return finished = false;
      case 3:
        dX = (_dX[0] + 2.0 * _dX[1] + 2.0 * _dX[2] + _dX[3]) / 6.0;
        dt = dt_fixed;
        return finished = true;
      default:
        break;
      }
      break;

    default:
      break;
    }

    // If none of the cases are matched, throw an error
    std::stringstream ss;
    ss << std::to_string(this->current_step) << "/" << std::to_string(this->steps);
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "current_step/steps=" + ss.str());
  }

  T get_U0_for_SPH() const {
    if (current_step == 0)
      return this->Xinit;
    else if (current_step == 1)
      return this->Xinit;
    else if (current_step == 2)
      return this->Xinit;
    else if (current_step == 3)
      return this->Xinit + (_dX[0] + 2. * _dX[1] + 2. * _dX[2]) / 6.;
    else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "currecnt_stepが4以上になっている．");
  };

  /* ------------------------------------------------------ */
  /*            微分を与えて，次の時刻に使うべきXを返す            */
  /* ------------------------------------------------------ */
  T getX(const T &dXdt_IN) const {
    if (this->steps == 0)
      return this->Xinit;
    switch (this->steps) {
    case 1:
      return this->Xinit + dXdt_IN * dt_fixed;

    case 2:
      switch (current_step) {
      case 0:
        return this->Xinit + dXdt_IN * dt_fixed;
      default:
        return this->Xinit + (_dX[0] + dXdt_IN * dt_fixed) / 2.0;
      }

    case 3:
      switch (current_step) {
      case 0:
        return this->Xinit + dXdt_IN * dt_fixed / 2.0;
      case 1:
        return this->Xinit + dXdt_IN * dt_fixed;
      default:
        return this->Xinit + (_dX[0] + 4.0 * _dX[1] + dXdt_IN * dt_fixed) / 6.0;
      }

    case 4:
      switch (current_step) {
      case 0:
      case 1:
        return this->Xinit + dXdt_IN * dt_fixed / 2.0;
      case 2:
        return this->Xinit + dXdt_IN * dt_fixed;
      default:
        return this->Xinit + (_dX[0] + 2.0 * _dX[1] + 2.0 * _dX[2] + dXdt_IN * dt_fixed) / 6.0;
      }

    default:
      std::stringstream ss;
      ss << std::to_string(this->current_step) << "/" << std::to_string(this->steps);
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "current_step/steps=" + ss.str());
    }
  }

  /* -------------------------------------------------------------------------- */

  //! QUERY: THE VECTOR TO BE ADDED TO THE ORIGINAL VECTOR TO REACH THE GIVEN POSITION AT THE NEXT TIME STEP

  // dXdt_IN == ((position - this->Xinit)*6.0 - (_dX[0] + 2.0 * _dX[1] + 2.0 * _dX[2])) dt_fixed

  T getVectorToReachAtNextTimeQ(const T &position) const {
    switch (this->steps) {
    case 1:
      return (position - this->Xinit) / dt_fixed;

    case 2:
      switch (current_step) {
      case 0:
        return (position - this->Xinit) / dt_fixed;
      default:
        return ((position - this->Xinit) * 2.0 - _dX[0]) / dt_fixed;
      }
    case 3:
      switch (current_step) {
      case 0:
        return (position - this->Xinit) * 2.0 / dt_fixed;
      case 1:
        return (position - this->Xinit) / dt_fixed;
      default:
        return ((position - this->Xinit) * 6.0 - (_dX[0] + 4.0 * _dX[1])) / dt_fixed;
      }
    case 4:
      switch (current_step) {
      case 0:
      case 1:
        return (position - this->Xinit) * 2.0 / dt_fixed;
      case 2:
        return (position - this->Xinit) / dt_fixed;
      default:
        return ((position - this->Xinit) * 6.0 - (_dX[0] + 2.0 * _dX[1] + 2.0 * _dX[2])) / dt_fixed;
      }

    default:
      std::stringstream ss;
      ss << std::to_string(this->current_step) << "/" << std::to_string(this->steps);
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "current_step/steps=" + ss.str());
    }
  }

  /* -------------------------------------------------------------------------- */

  T getdXdt(const T &dXdt_IN) const {
    if (this->steps == 1)
      return dXdt_IN; //! 何もプッシュしていない初期状態
    else if (this->steps == 2) {
      switch (current_step) {
      case 0: //! 何もプッシュしていない初期状態
        return dXdt_IN;
      default:
        return (_dX[0] / dt_fixed + dXdt_IN) / 2.0;
      }
    } else if (this->steps == 3) {
      switch (current_step) {
      case 0: //! 何もプッシュしていない初期状態
        return dXdt_IN;
      case 1: // -> {f2, f1(t,v(t))}
        return dXdt_IN;
      default:
        return (_dX[0] / dt_fixed + 4. * _dX[1] / dt_fixed + dXdt_IN) / 6.;
      }
    } else if (this->steps == 4) {
      switch (current_step) {
      case 0: //! 何もプッシュしていない初期状態
        return dXdt_IN;
      case 1: // -> {f2, f1(t,v(t))}
        return dXdt_IN;
      case 2: // -> {f3, f2, f1(t,v(t))}
        return dXdt_IN;
      default: // -> {f3, f2, f1(t,v(t))}
        return (_dX[0] / dt_fixed + 2. * _dX[1] / dt_fixed + 2. * _dX[2] / dt_fixed + dXdt_IN) / 6.;
      }
    } else {
      std::stringstream ss;
      ss << std::to_string(this->current_step) << "/" << std::to_string(this->steps);
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "必要以上に微分をプッシュしている．current_step/steps=" + ss.str());
    }
  };
};

template <typename T> struct RungeKutta : public RungeKuttaCommon<T> {
  RungeKutta(const double dt_IN, const double t0, const T &X0, int stepsIN) : RungeKuttaCommon<T>(dt_IN, t0, X0, stepsIN) {};
  RungeKutta() : RungeKuttaCommon<T>() {};
};
template <> struct RungeKutta<std::vector<double>> : public RungeKuttaCommon<std::vector<double>> {
  RungeKutta(const double dt_IN, const double t0, const std::vector<double> &X0, int stepsIN) : RungeKuttaCommon<std::vector<double>>(dt_IN, t0, X0, stepsIN) { this->dX = V_d(X0.size(), 0.); };
  RungeKutta() : RungeKuttaCommon<std::vector<double>>() {};
};
template <> struct RungeKutta<double> : public RungeKuttaCommon<double> {
  RungeKutta(const double dt_IN, const double t0, const double &X0, int stepsIN) : RungeKuttaCommon<double>(dt_IN, t0, X0, stepsIN) { this->dX = 0.; };
  RungeKutta() : RungeKuttaCommon<double>() {};
};
template <std::size_t N> struct RungeKutta<std::array<double, N>> : public RungeKuttaCommon<std::array<double, N>> {
  RungeKutta(const double dt_IN, const double t0, const std::array<double, N> &X0, int stepsIN) : RungeKuttaCommon<std::array<double, N>>(dt_IN, t0, X0, stepsIN) { this->dX.fill(0.); };
  RungeKutta() : RungeKuttaCommon<std::array<double, N>>() {};
};

/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT ODE::LeapFrog

### LeapFrog

リープフロッグの１回の計算で溜まる誤差は$O({\Delta t}^3)$となる．
時間間隔$\Delta t$が変化する場合でも使える形でプログラムしている（\ref{ODE:LeapFrog}{LeapFrogのクラス}）．
$\Delta t$が変化する場合，"半分蹴って-移動-半分蹴って"，"半分蹴って-移動-半分蹴って"の手順を繰り返す．
\ref{ODE:LeapFrog}{LeapFrogのクラス}

*/
// \label{ODE:LeapFrog}
template <typename T> class LeapFrog {
public:
  LeapFrog() : is_first(true), finished(false) {};
  LeapFrog(double dt, double t0, const T &x0, const T &v0) : dt(dt), t(t0), x(x0), v(v0), is_first(true), finished(false) {};

  void initialize(double dt, double t0, const T &x0, const T &v0) {
    this->dt = dt;
    this->t = t0;
    this->x = x0;
    this->v = v0;
    this->is_first = true;
    this->finished = false;
  }

  bool is_first, finished;
  T v_full;

  void push(const T &a) {
    double half_dt = 0.5 * dt;
    if (is_first) {
      t_old = t;
      x_old = x;
      v_old = v;
      a_old = a;
      v += half_dt * a; // half-step update of v
      x += dt * v;      // full-step update of x
      t += dt;
      is_first = false;
      finished = false;
    } else {
      t_old = t;
      x_old = x;
      v_old = v;
      a_old = a;
      v += half_dt * a; // half-step update of v
      is_first = true;
      finished = true;
    }
  }

  // The Leapfrog method (also known as the midpoint method)
  void repush(const T &a) {
    double half_dt = 0.5 * dt;
    if (!is_first) {
      v = v_old + half_dt * a; // half-step update of v
      x = x_old + dt * v;      // full-step update of x
      t = t_old + dt;
      is_first = false;
      finished = false;
    } else {
      v = v_old + half_dt * a; // half-step update of v
      is_first = true;
      finished = true;
    }
  }

  double get_dt() const { return dt; }
  double get_t() const { return t; }
  const T &get_x() const { return x; }
  const T &get_v() const { return v; }

  T get_v(const T &a) {
    double half_dt = 0.5 * dt;
    return v_old + half_dt * a; // half-step update of v
  };

  T get_x(const T &a) {
    if (is_first) {
      double half_dt = 0.5 * dt;
      auto V = v + half_dt * a; // half-step update of v
      return x + dt * V;        // full-step update of x
    } else {
      return x;
    }
  };

  T get_x_if_v(const T &v) {
    if (is_first) {
      return x + dt * v; // full-step update of x
    } else {
      return x;
    }
  };

private:
  double dt, t, t_old;
  T x, v, a_old, v_old, x_old;
};

/* -------------------------------------------------------------------------- */

template <typename T> class VelocityVerlet {
public:
  VelocityVerlet(double dt, double t0, const T &x0, const T &v0, const T &a0) : dt(dt), t(t0), x(x0), v(v0), a(a0) {};

  void push(const T &new_a) {
    v += dt * (a + new_a) * 0.5; // this should be perfomed after x is updated though this sometimes stabilizes the calculation
    x += dt * v + 0.5 * dt * dt * new_a;
    t += dt;
    a = new_a;
  }

  double get_t() const { return t; }
  const T &get_x() const { return x; }
  const T &get_v() const { return v; }
  const T &get_a() const { return a; }

private:
  double dt, t;
  T x, v, a;
};

/* -------------------------------------------------------------------------- */

template <typename T> class Beeman {
public:
  Beeman(double dt, double t0, const T &x0, const T &v0, const T &a0, const T &old_a) : dt(dt), t(t0), x(x0), v(v0), a(a0), old_a(old_a) {};
  Beeman(double dt, double t0, const T &x0, const T &v0, const T &a0) : dt(dt), t(t0), x(x0), v(v0), a(a0), old_a(a0) {};

  void push(const T &new_a) {
    x += v * dt + (2.0 / 3) * a * dt * dt - (1.0 / 6) * old_a * dt * dt;
    T old_v = v;
    v += (1.0 / 3) * new_a * dt + (5.0 / 6) * a * dt - (1.0 / 6) * old_a * dt;
    t += dt;
    old_a = a;
    a = new_a;
  }

  double get_t() const { return t; }
  const T &get_x() const { return x; }
  const T &get_v() const { return v; }
  const T &get_a() const { return a; }

private:
  double dt, t;
  T x, v, a, old_a;
};

/* -------------------------------------------------------------------------- */

template <typename T> class Gear {
public:
  Gear(double dt, double t0, const T &x0, const T &v0, const T &a0) : dt(dt), t(t0), x(x0), v(v0), a(a0) {};

  void push(const T &new_a) {
    T pred_x = x + v * dt + 0.5 * a * dt * dt;
    T pred_v = v + a * dt;
    T error = new_a - a;
    x = pred_x + (1.0 / 6) * error * dt * dt;
    v = pred_v + 0.5 * error * dt;
    t += dt;
    a = new_a;
  }

  double get_t() const { return t; }
  const T &get_x() const { return x; }
  const T &get_v() const { return v; }
  const T &get_a() const { return a; }

private:
  double dt, t;
  T x, v, a;
};
