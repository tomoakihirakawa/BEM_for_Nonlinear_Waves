#pragma once

#include "basic_arithmetic_vector_operations.hpp"
#include "basic_vectors.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

template <typename T> struct NewtonRaphson_Common {
  T X, dX;
  NewtonRaphson_Common(const T &Xinit) : X(Xinit), dX(Xinit) {};
  void initialize(const T &Xin) { X = Xin; };
};

template <typename T> struct NewtonRaphson : public NewtonRaphson_Common<T> {
  NewtonRaphson(const T &Xinit) : NewtonRaphson_Common<T>(Xinit) {};
};

template <> struct NewtonRaphson<V_d> : public NewtonRaphson_Common<V_d> {
  NewtonRaphson(const V_d &Xinit) : NewtonRaphson_Common<V_d>(Xinit) {};
  void update(const V_d &F, const VV_d &dFdx) {
    lapack_svd lu(dFdx, dX, -F);
    X += dX;
  };

  void update(const V_d &F, const VV_d &dFdx, const double a) {
    lapack_svd lu(dFdx, dX, -F);
    X += a * dX;
  };

  // pass lambda function to constrain the solution and update X and dX
  void constrains(const std::function<V_d(V_d)> &constraint) {
    /*
    X^* = X + dX
    constrained(X^*) = X^* + dX* = X + dX + dX^*
    constrained(X^*) - (X + dX) = dX^*

    constrained(X^*) = X + dX + dX^* = X + dX^{c}
    dX^{c} = dX + dX^*
    */

    //@ tedious way to update X and dX
    // auto X_ast = X;
    // auto dX_ast = constraint(X) - X_ast;
    // auto dX_c = dX + dX_ast;
    // X = constraint(X);
    // dX = dX_c;

    //@ concise way to update X and dX
    dX -= X;
    dX += (X = constraint(X));
  };
};
/* ------------------------------------------------------ */
template <> struct NewtonRaphson<double> : public NewtonRaphson_Common<double> {
  NewtonRaphson(const double Xinit = 0) : NewtonRaphson_Common<double>(Xinit) {};
  void update(const double F, const double dFdx) {
    if (std::abs(dFdx) > 1E-20)
      X += (dX = F / (-dFdx));
  };
  void update(const double F, const double dFdx, const double a) {
    if (std::abs(dFdx) > 1E-20)
      X += a * (dX = F / (-dFdx));
  };
};
template <> struct NewtonRaphson<Tdd> : public NewtonRaphson_Common<Tdd> {
  NewtonRaphson(const Tdd &Xinit) : NewtonRaphson_Common<Tdd>(Xinit) {};
  void update(const Tdd &F, const T2Tdd &dFdx) {
    lapack_svd(dX, dFdx, -F);
    X += dX;
  };
  void update(const Tdd &F, const T2Tdd &dFdx, const double a) {
    lapack_svd(dX, dFdx, -F);
    X += a * dX;
  };
};
template <> struct NewtonRaphson<Tddd> : public NewtonRaphson_Common<Tddd> {
  NewtonRaphson(const Tddd &Xinit) : NewtonRaphson_Common<Tddd>(Xinit) {};
  void update(const Tddd &F, const T3Tddd &dFdx) {
    lapack_svd(dX, dFdx, -F);
    X += dX;
  };

  void constrains(const std::function<Tddd(const Tddd &)> &constraint) {
    dX -= X;
    dX += (X = constraint(X));
  };

  void update(const Tddd &F, const T3Tddd &dFdx, const std::function<Tddd(const Tddd &)> &constraint) {
    lapack_svd(dX, dFdx, -F);
    X += dX;
    constraint(X);
  };
};
template <> struct NewtonRaphson<T4d> : public NewtonRaphson_Common<T4d> {
  NewtonRaphson(const T4d &Xinit) : NewtonRaphson_Common<T4d>(Xinit) { this->ans = V_d(4, 0.); };
  V_d ans;
  void update(const T4d &F, const T4T4d &dFdx) {
    lapack_svd(dX, dFdx, -F);
    X += dX;
    // ludcmp lu(ToVector(dFdx));
    // lu.solve(ToVector(-F), ans);
    // std::get<0>(X) += (std::get<0>(dX) = ans[0]);
    // std::get<1>(X) += (std::get<1>(dX) = ans[1]);
    // std::get<2>(X) += (std::get<2>(dX) = ans[2]);
    // std::get<3>(X) += (std::get<3>(dX) = ans[3]);
  };
};

/* ------------------------------------------------------ */
/*                           利用例                        */
/* ------------------------------------------------------ */

//! 成分毎に与えること
//! v = vx*ex + vy*ey + vz*ezの場合
//! {vx,vy,vz}, {ex,ey,ez}を与える
//! ex = {1,0,0}, ey = {0,1,0}, ez = {0,0,1}でなくとも良いが，ほぼ正規直交座標系であることが望ましいだろう
inline std::array<double, 3> optimalVector(std::vector<double> Vsample, std::vector<Tddd> Directions, const Tddd &Vinit, std::vector<double> weights, std::array<double, 3> &convergence_info) {

  if (Vsample.size() == 1)
    return Vsample[0] * Directions[0];

  double mean = 0.;
  for (const auto &v : Vsample)
    mean += std::abs(v);
  mean /= Vsample.size();
  if (mean == 0)
    return {0., 0., 0.};

  const double tolerance = 1E-13 * mean;
  convergence_info = {0., 0., 0.};

  const double threshold_angle_in_rad = 10 * M_PI / 180.;

  //! 与えられている情報が不十分な場合がある．

  /* -------------------------------------------------------------------------- */
  /*                         make directions into groups                        */
  /* -------------------------------------------------------------------------- */

  std::vector<std::tuple<Tddd, std::vector<Tddd>>> direction_groups;
  for (auto &dir : Directions) {
    bool colinear = false;
    for (auto &[representative_dir, vec] : direction_groups) {
      if (isFlat(representative_dir, dir, threshold_angle_in_rad)) {
        vec.push_back(Normalize(dir));
        representative_dir = Normalize(Mean(vec));
        colinear = true;
        break;
      }
    }
    if (!colinear) {
      direction_groups.push_back({Normalize(dir), {Normalize(dir)}});
    }
  }

  if (direction_groups.size() == 1) {
    Tddd ret = {0., 0., 0.};
    for (std::size_t i = 0; i < Vsample.size(); ++i)
      ret += Vsample[i] * Directions[i];
    return ret / Vsample.size();
  } else if (direction_groups.size() == 2) {
    Vsample.push_back(0.);
    Directions.push_back(Normalize(Cross(Normalize(Mean(std::get<1>(direction_groups[0]))), Normalize(Mean(std::get<1>(direction_groups[1]))))));
    weights.push_back(Mean(weights));
  }

  for (auto &d : Directions)
    d = Normalize(d);

  /* -------------------------------------------------------------------------- */

  auto diff = [&](const Tddd &U, const std::size_t i) -> double { return Dot(U, Directions[i]) - Vsample[i]; };

  auto optimizing_function = [&](const Tddd &U) -> double {
    double S = 0;
    for (std::size_t i = 0; i < Vsample.size(); ++i)
      S += weights[i] * std::pow(diff(U, i), 2);
    return 0.5 * S;
  };

  NewtonRaphson<Tddd> NR(Vinit);
  Tddd grad;
  T3Tddd hess;
  int iteration = 0;
  for (iteration = 0; iteration < 50; ++iteration) {
    grad.fill(0.);
    for (auto &row : hess)
      row.fill(0.);
    for (std::size_t i = 0; i < Vsample.size(); ++i) {
      grad += weights[i] * Directions[i] * diff(NR.X, i);
      hess += weights[i] * TensorProduct(Directions[i], Directions[i]);
    }
    NR.update(grad, hess);
    if ((Norm(grad) < tolerance || optimizing_function(NR.X) < tolerance))
      break;
  }

  std::get<0>(convergence_info) = (double)iteration;
  std::get<1>(convergence_info) = Norm(grad);
  std::get<2>(convergence_info) = optimizing_function(NR.X);

  return NR.X;
}

inline std::array<double, 3> optimalVector(const std::vector<double> &Vsample, const std::vector<Tddd> &Directions, const Tddd &Vinit) {
  std::vector<double> weights(Vsample.size(), 1.);
  std::array<double, 3> convergence_info;
  return optimalVector(Vsample, Directions, Vinit, weights, convergence_info);
}

inline std::array<double, 3> optimalVector(const std::vector<double> &Vsample, const std::vector<Tddd> &Directions, const Tddd &Vinit, const std::vector<double> &weights) {
  std::array<double, 3> convergence_info;
  return optimalVector(Vsample, Directions, Vinit, weights, convergence_info);
}

inline std::array<double, 3> optimalVector(const std::vector<double> &Vsample, const std::vector<Tddd> &Directions, const Tddd &Vinit, std::array<double, 3> &convergence_info) {
  std::vector<double> weights(Vsample.size(), 1.);
  return optimalVector(Vsample, Directions, Vinit, weights, convergence_info);
}

template <std::size_t N> std::array<double, N> optimumVector(const std::vector<std::array<double, N>> &sample_vectors, const std::array<double, N> &init_vector, const double tolerance = 1E-12) {
  std::array<NewtonRaphson<double>, N> NRs;
  for (std::size_t i = 0; i < N; ++i)
    NRs[i].X = init_vector[i];
  std::array<double, N> Fs, dFs;
  double w, drdx;
  for (auto j = 0; j < 500; ++j) {
    Fs.fill(0);
    dFs.fill(0);
    for (const auto &vec : sample_vectors) {
      w = 1;
      drdx = -w;
      for (std::size_t i = 0; i < N; ++i) {
        // Fs[i] += (w * (vec[i] - NRs[i].X)) * drdx;  //<- d/dx (d*d)
        // dFs[i] += drdx * drdx;
        // use std::fma
        Fs[i] = std::fma(w * (vec[i] - NRs[i].X), drdx, Fs[i]);
        dFs[i] = std::fma(drdx, drdx, dFs[i]);
      }
    }
    bool converged = true;
    for (std::size_t i = 0; i < N; ++i) {
      NRs[i].update(Fs[i], dFs[i]);
      if (std::abs(NRs[i].dX) >= tolerance)
        converged = false;
    }
    if (converged)
      break;
  }
  std::array<double, N> result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = NRs[i].X;
  return result;
}

// template <std::size_t N>
// std::array<double, N> optimumVector(const std::vector<std::array<double, N>> &sample_vectors,
//                                     const std::array<double, N> &init_vector,
//                                     const std::vector<double> &weights,
//                                     const double tolerance = 1E-12) {
//    if (weights.size() != sample_vectors.size())
//       throw std::runtime_error("The size of the weights vector must match the size of the sample_vectors vector.");

//    double max_weight = *std::max_element(weights.begin(), weights.end());
//    std::vector<double> normalized_weights(weights.size());
//    std::transform(weights.begin(), weights.end(), normalized_weights.begin(), [&](double w) { return w / max_weight; });

//    std::array<NewtonRaphson<double>, N> NRs;
//    for (std::size_t i = 0; i < N; ++i)
//       NRs[i].X = init_vector[i];

//    std::array<double, N> Fs, dFs;
//    double w, drdx;
//    for (int iteration = 0; iteration < 500; ++iteration) {
//       Fs.fill(0);
//       dFs.fill(0);
//       for (size_t i = 0; i < sample_vectors.size(); ++i) {
//          w = normalized_weights[i];
//          drdx = -w;
//          for (size_t j = 0; j < N; ++j) {
//             // Fs[j] += w * (sample_vectors[i][j] - NRs[j].X) * drdx;
//             // dFs[j] += drdx * drdx;
//             // use std::fma
//             Fs[j] = std::fma(w * (sample_vectors[i][j] - NRs[j].X), drdx, Fs[j]);
//             dFs[j] = std::fma(drdx, drdx, dFs[j]);
//          }
//       }

//       bool converged = true;
//       for (std::size_t i = 0; i < N; ++i) {
//          NRs[i].update(Fs[i], dFs[i]);
//          if (std::abs(NRs[i].dX) >= tolerance) converged = false;
//       }
//       if (converged) break;
//    }

//    std::array<double, N> result;
//    for (std::size_t i = 0; i < N; ++i) result[i] = NRs[i].X;
//    return result;
// }

// template <std::size_t N>
// std::array<double, N> optimumVector(const std::vector<std::array<double, N>> &sample_vectors,
//                                     const std::array<double, N> &init_vector,
//                                     const std::vector<double> &weights,
//                                     const std::function<std::array<double, N>(const std::array<double, N> &)> &constraint,
//                                     const double tolerance = 1E-12) {
//    if (weights.size() != sample_vectors.size())
//       throw std::runtime_error("The size of the weights vector must match the size of the sample_vectors vector.");

//    double max_weight = *std::max_element(weights.begin(), weights.end());
//    std::vector<double> normalized_weights(weights.size());
//    std::transform(weights.begin(), weights.end(), normalized_weights.begin(), [&](double w) { return w / max_weight; });

//    std::array<NewtonRaphson<double>, N> NRs;
//    for (std::size_t i = 0; i < N; ++i)
//       NRs[i].X = init_vector[i];

//    std::array<double, N> Fs, dFs;
//    double w, drdx;
//    for (int iteration = 0; iteration < 500; ++iteration) {
//       Fs.fill(0);
//       dFs.fill(0);
//       for (size_t i = 0; i < sample_vectors.size(); ++i) {
//          w = normalized_weights[i];
//          drdx = -w;
//          for (size_t j = 0; j < N; ++j) {
//             // Fs[j] += w * (sample_vectors[i][j] - NRs[j].X) * drdx;
//             // dFs[j] += drdx * drdx;
//             // use std::fma
//             Fs[j] = std::fma(w * (sample_vectors[i][j] - NRs[j].X), drdx, Fs[j]);
//             dFs[j] = std::fma(drdx, drdx, dFs[j]);
//          }
//       }

//       bool converged = true;
//       for (std::size_t i = 0; i < N; ++i) {
//          NRs[i].update(Fs[i], dFs[i]);
//       }
/*
//       std::array<NewtonRaphson<double>, N> NRs;としたので，
//       直接的なconstraintの適用は難しくなった．
//       */

//       std::array<double, N> tmpX;
//       for (std::size_t i = 0; i < N; ++i)
//          tmpX[i] = NRs[i].X;

//       std::array<double, N> constraintX = constraint(tmpX);

//       for (std::size_t i = 0; i < N; ++i) {
//          NRs[i].dX -= NRs[i].X;
//          NRs[i].dX += (NRs[i].X = constraintX[i]);
//          if (iteration > 5 || std::abs(NRs[i].dX) >= tolerance) converged = false;
//       }

//       if (converged) break;
//    }

//    std::array<double, N> result;
//    for (std::size_t i = 0; i < N; ++i) result[i] = NRs[i].X;
//    return result;
// }

inline double optimumValue(const std::vector<double> &sample_values, const double init_value, std::vector<double> weights, const double tolerance = 1E-12) {
  if (weights.size() != sample_values.size())
    throw std::runtime_error("The size of the weights vector must match the size of the sample_vectors vector.");

  double m = 0;
  for (const auto &max : weights)
    if (m < std::abs(max))
      m = std::abs(max);
  weights /= m;
  NewtonRaphson<double> NRs;
  NRs.X = init_value;

  double Fs, dFs;
  double w, drdx;
  for (int iteration = 0; iteration < 500; ++iteration) {
    Fs = 0;
    dFs = 0;
    for (size_t i = 0; i < sample_values.size(); ++i) {
      w = weights[i];
      drdx = -w;
      Fs += w * (sample_values[i] - NRs.X) * drdx;
      dFs += drdx * drdx;
    }

    bool converged = true;
    NRs.update(Fs, dFs);
    if (std::abs(NRs.dX) >= tolerance)
      converged = false;
    if (converged)
      break;
  }

  return NRs.X;
}

/* -------------------------------------------------------------------------- */

struct DispersionRelation {
  /*DOC_EXTRACT dispersionRelation

  w = \sqrt{|{\bf g}| k  \tanh{h k}}  = \sqrt{|{\bf g}| \frac{2\pi}{L}  \tanh{h \frac{2\pi}{L}}}

  */
  double w;
  double T;
  double k;
  double L;
  double h;

  DispersionRelation() : w(0), h(0), T(0), k(0), L(0) {};

  DispersionRelation(const double wIN, const double hIN) { set_w_h(wIN, hIN); };

  void set_T_h(const double TIN, const double hIN) { set_w_h(2 * M_PI / TIN, hIN); };

  void set_w_h(const double wIN, const double hIN) {
    this->w = wIN;
    this->T = 2 * M_PI / this->w;
    this->h = hIN;
    bool found = false;
    const std::vector<double> init_L = {0.5, 1., 2., 5., 10., 20., 50., 100., 200., 500., 1000., 2000.};
    for (const double l : init_L) {
      NewtonRaphson nr(2. * M_PI / l);
      for (auto i = 0; i < 30; i++) {
        nr.update(omega(nr.X, h) - w, domegadk(nr.X, h));
        if (std::abs(omega(nr.X, h) - w) < 1E-10) {
          found = true;
          this->k = std::abs(nr.X);
          break;
        }
      }
      if (found)
        break;
    }
    this->L = 2 * M_PI / this->k;
  };

  void set_L_h(const double LIN, const double hIN) {
    this->L = LIN;
    this->h = hIN;
    this->k = 2 * M_PI / this->L;
    this->w = omega(this->k, this->h);
    this->T = 2 * M_PI / this->w;
  };

  double omega(const double k, const double h) { return Sqrt(_GRAVITY_ * k * Tanh(h * k)); };

  double domegadk(const double k, const double h) { return (_GRAVITY_ * (h * k * Power(Sech(h * k), 2) + Tanh(h * k))) / (2. * Sqrt(_GRAVITY_ * k * Tanh(h * k))); };
};

inline double cosh_kzh_cosh_kh(const double k, const double h, const double z) noexcept {
  const double kh = k * h;
  const double kz = k * z;
  if (std::abs(kz) < 40)
    return std::cosh(kz) * (1 + std::tanh(kz) * std::tanh(kh));
  else
    return 0.0;
}

inline double sinh_kzh_cosh_kh(const double k, const double h, const double z) noexcept {
  const double kh = k * h;
  const double kz = k * z;

  if (std::abs(kz) < 40)
    return std::cosh(kz) * (std::tanh(kz) + std::tanh(kh));
  else
    return 0.0;
}

struct WaterWaveTheory {

  double h;
  double L;
  double T;
  double w;
  double k;
  double c;

  double phase_shift = 0.; //[rad]

  double theta = 0.; //[rad]波の向き

  WaterWaveTheory() : h(0), L(0), T(0), w(0), k(0), c(0) {};

  void set_T_h(const double TIN, const double hIN) { set_w_h(2 * M_PI / TIN, hIN); };

  void set_w_h(const double wIN, const double hIN) {
    this->w = wIN;
    this->T = 2 * M_PI / this->w;
    this->h = hIN;
    bool found = false;
    const std::vector<double> init_L = {0.5, 1., 2., 5., 10., 20., 50., 100., 200., 500., 1000., 2000.};
    for (const double l : init_L) {
      NewtonRaphson nr(2. * M_PI / l);
      for (auto i = 0; i < 30; i++) {
        nr.update(omega(nr.X, h) - w, dwdk(nr.X, h));
        if (std::abs(omega(nr.X, h) - w) < 1E-10) {
          found = true;
          this->k = std::abs(nr.X);
          break;
        }
      }
      if (found)
        break;
    }
    this->L = 2 * M_PI / this->k;
  };

  void set_L_h(const double LIN, const double hIN) {
    this->L = LIN;
    this->h = hIN;
    this->k = 2 * M_PI / this->L;
    this->w = omega(this->k, this->h);
    this->T = 2 * M_PI / this->w;
  };

  double omega(const double k, const double h) { return Sqrt(_GRAVITY_ * k * Tanh(h * k)); };

  double dwdk(const double k, const double h) { return (_GRAVITY_ * (h * k * Power(Sech(h * k), 2) + Tanh(h * k))) / (2. * Sqrt(_GRAVITY_ * k * Tanh(h * k))); };

  double clampZ(const double z, const double k) const {
    // constexpr double LOG_LIMIT = 12.0; // e^12 ≈ 1.6e5
    // const double k_eff = std::max(k, 1e-12);
    // const double z_upper = std::min(0.0, LOG_LIMIT / k_eff);
    // return std::clamp(z, std::numeric_limits<double>::lowest(), z_upper);
    return z;
  }

  double phi(const std::array<double, 3> &X, const double t, const double A, const double bottom_z) const {
    auto [x, y, z] = X;
    z = clampZ(z - h - bottom_z, k); //! distance from the bottom
    // return A * _GRAVITY_ / w * cosh_kzh_cosh_kh(k,h,z) * std::sin(k * x - w * t);
    double kx = k * std::cos(theta);
    double ky = k * std::sin(theta);
    return A * _GRAVITY_ / w * cosh_kzh_cosh_kh(k, h, z) * std::sin(kx * x + ky * y - w * t + phase_shift);
  };

  double phi(const std::array<double, 3> &X, const double t) const { return this->phi(X, t, A, bottom_z); };

  double A;
  double bottom_z = 0;

  std::array<double, 3> gradPhi(const std::array<double, 3> &X, const double t) const {
    auto [x, y, z] = X;
    double kx = k * std::cos(theta);
    double ky = k * std::sin(theta);
    z = clampZ(z - h - bottom_z, k); //! distance from the bottom
    return {A * _GRAVITY_ * kx / w * cosh_kzh_cosh_kh(k, h, z) * std::cos(kx * x + ky * y - w * t + phase_shift), A * _GRAVITY_ * ky / w * cosh_kzh_cosh_kh(k, h, z) * std::cos(kx * x + ky * y - w * t + phase_shift), A * _GRAVITY_ * k / w * sinh_kzh_cosh_kh(k, h, z) * std::sin(kx * x + ky * y - w * t + phase_shift)};
    // return {-A * w * std::cosh(k * (z + h)) / std::sinh(k * h) * std::cos(k * x - w * t),
    //         0.,
    //         A * w * std::sinh(k * (z + h)) / std::sinh(k * h) * std::sin(k * x - w * t)};
  };

  std::array<double, 3> gradPhi_t(const std::array<double, 3> &X, const double t, const double A) const {
    auto [x, y, z] = X;
    double kx = k * std::cos(theta);
    double ky = k * std::sin(theta);
    z = clampZ(z - h - bottom_z, k); //! distance from the bottom
    return {w * A * _GRAVITY_ * kx / w * cosh_kzh_cosh_kh(k, h, z) * std::sin(kx * x + ky * y - w * t + phase_shift), w * A * _GRAVITY_ * ky / w * cosh_kzh_cosh_kh(k, h, z) * std::sin(kx * x + ky * y - w * t + phase_shift), -w * A * _GRAVITY_ * k / w * sinh_kzh_cosh_kh(k, h, z) * std::cos(kx * x + ky * y - w * t + phase_shift)};
  };

  std::array<double, 3> gradPhi_t(const std::array<double, 3> &X, const double t) const { return this->gradPhi_t(X, t, A); };

  double eta(const std::array<double, 3> &X, const double t, const double A) const {
    auto [x, y, z] = X;
    double kx = k * std::cos(theta);
    double ky = k * std::sin(theta);
    return A * std::cos(kx * x + ky * y - w * t + phase_shift) + h + bottom_z;
  };

  double eta(const std::array<double, 3> &X, const double t) const { return this->eta(X, t, A); };

  double etaZeroAtRest(const std::array<double, 3> &X, const double t) const {
    auto [x, y, z] = X;
    double kx = k * std::cos(theta);
    double ky = k * std::sin(theta);
    return A * std::cos(kx * x + ky * y - w * t + phase_shift);
  };
};

enum class SpectrumType { BRETSCHNEIDER_MITSUYA, JONSWAP };
inline std::ostream &operator<<(std::ostream &os, SpectrumType type) {
  switch (type) {
  case SpectrumType::BRETSCHNEIDER_MITSUYA:
    os << "BRETSCHNEIDER_MITSUYA";
    break;
  case SpectrumType::JONSWAP:
    os << "JONSWAP";
    break;
  default:
    os << "UNKNOWN";
    break;
  }
  return os;
}

struct RandomWaterWaveTheory {

  /* ------------------------------------------------------ */
  double H13; // 有義波高
  double T13; // 有義周期
  double L13;
  double h;        // 水深
  double bottom_z; // 底面のz座標

  static constexpr std::size_t N = 1000; // 分割数
  double f_min, f_max, df;
  std::array<std::shared_ptr<WaterWaveTheory>, N> waves;
  std::mt19937 gen;
  std::uniform_real_distribution<double> random_phase;

  // JONSWAP parameters
private:
  SpectrumType spectrum_type = SpectrumType::BRETSCHNEIDER_MITSUYA;

public:
  double gamma = 3.3; // Peak enhancement factor for JONSWAP (typical value)
  double betaJ = 0.;

  /* ------------------------------------------------------ */

  // 修正ブレットシュナイダー・光易型スペクトル
  double spectrum_bretschneider_mitsuya(double f) const { return 0.205 * std::pow(H13, 2.) * std::pow(T13, -4.) * std::pow(f, -5.) * std::exp(-0.75 * std::pow(T13 * f, -4.)); }

  // JONSWAP spectrum
  double spectrum_jonswap(double f) const {
    double Tp = T13 / (1 - 0.132 * std::pow(gamma + 0.2, -0.559)); // T13 to peak period Tp
    double fp = 1.0 / Tp;
    double sigma = (f <= fp) ? 0.07 : 0.09;
    return betaJ * std::pow(H13, 2.) * std::pow(Tp, -4.) * std::pow(f, -5.) * std::exp(-1.25 * std::pow(Tp * f, -4.)) * std::pow(gamma, std::exp(-0.5 * std::pow(-(Tp * f - 1) / sigma, 2.)));
  }

  double spectrum(double f) const {
    switch (spectrum_type) {
    case SpectrumType::JONSWAP:
      return spectrum_jonswap(f);
    case SpectrumType::BRETSCHNEIDER_MITSUYA:
    default:
      return spectrum_bretschneider_mitsuya(f);
    }
  }

  RandomWaterWaveTheory() : H13(0), T13(0), h(0), bottom_z(0), gen(std::random_device()()), random_phase(0.0, 2.0 * M_PI) {};

  RandomWaterWaveTheory(double H13, double T13, double h, double bottom_z) : H13(H13), T13(T13), h(h), bottom_z(bottom_z), gen(std::random_device{}()), random_phase(0.0, 2.0 * M_PI) {
    DispersionRelation disp;
    disp.set_T_h(T13, h);
    this->L13 = disp.L;
    f_min = 0.5 / T13;
    f_max = 2.0 / T13;
    df = (f_max - f_min) / N;
    for (int i = 0; i < N; i++) {
      auto a_wave = std::make_shared<WaterWaveTheory>();
      double f = f_min + i * df;
      double T_wave = 1.0 / f;
      a_wave->A = std::sqrt(2 * spectrum(f) * df);
      a_wave->bottom_z = bottom_z;
      a_wave->set_T_h(T_wave, h);
      a_wave->phase_shift = random_phase(gen); // ランダムな位相シフトを設定
      if (!isFinite(a_wave->eta({0, 0, 0}, 0)))
        throw std::runtime_error("eta is not finite");
      waves[i] = a_wave;
    }
  }

  void setSpectrumType(SpectrumType type) {
    this->spectrum_type = type;
    this->betaJ = 0.0624 / (0.23 + 0.0336 * gamma - 0.185 / (1.9 + gamma)) * (1.094 - 0.01915 * std::log(gamma));
    DispersionRelation disp;
    disp.set_T_h(T13, h);
    this->L13 = disp.L;
    f_min = 0.5 / T13;
    f_max = 3.0 / T13;
    df = (f_max - f_min) / N;
    for (int i = 0; i < N; i++) {
      auto a_wave = std::make_shared<WaterWaveTheory>();
      double f = f_min + i * df;
      double T_wave = 1.0 / f;
      a_wave->A = std::sqrt(2 * spectrum(f) * df);
      a_wave->bottom_z = bottom_z;
      a_wave->set_T_h(T_wave, h);
      a_wave->phase_shift = random_phase(gen); // ランダムな位相シフトを設定
      if (!isFinite(a_wave->eta({0, 0, 0}, 0)))
        throw std::runtime_error("eta is not finite");
      waves[i] = a_wave;
    }
  };

  std::array<double, 3> gradPhi(const std::array<double, 3> &X, const double t) const {
    std::array<double, 3> ret{0., 0., 0.};
    for (const auto &wave : waves) {
      auto grad = wave->gradPhi(X, t);
      ret[0] += grad[0];
      ret[1] += grad[1];
      ret[2] += grad[2];
    }
    return ret;
  };

  std::array<double, 3> gradPhi_t(const std::array<double, 3> &X, const double t) const {
    std::array<double, 3> ret{0., 0., 0.};
    for (const auto &wave : waves) {
      auto grad = wave->gradPhi_t(X, t);
      ret[0] += grad[0];
      ret[1] += grad[1];
      ret[2] += grad[2];
    }
    return ret;
  };

  double phi(const std::array<double, 3> &X, const double t) const {
    double sum = 0.0;
    for (const auto &wave : waves)
      sum += wave->phi(X, t);
    return sum;
  }

  double eta(const std::array<double, 3> &X, const double t) const {
    double sum = 0.0;
    for (const auto &wave : waves)
      sum += wave->etaZeroAtRest(X, t);
    return sum + h + bottom_z;
  }

  friend std::ostream &operator<<(std::ostream &os, const RandomWaterWaveTheory &theory);
};

inline std::ostream &operator<<(std::ostream &os, const RandomWaterWaveTheory &theory) {
  std::ios::fmtflags f(os.flags());
  os << std::fixed << std::setprecision(6) << "RandomWaterWaveTheory {\n"
     << "  spectrum_type = " << theory.spectrum_type << "\n"
     << "  gamma         = " << theory.gamma << "\n"
     << "  betaJ         = " << theory.betaJ << "\n"
     << "  H13           = " << theory.H13 << "  (significant wave height)\n"
     << "  T13           = " << theory.T13 << "  (mean period)\n"
     << "  L13           = " << theory.L13 << "\n"
     << "  h             = " << theory.h << "\n"
     << "  bottom_z      = " << theory.bottom_z << "\n"
     << "  N             = " << theory.N << "  (number of components)\n"
     << "  f_min         = " << theory.f_min << "\n"
     << "  f_max         = " << theory.f_max << "\n"
     << "  df            = " << theory.df << "\n";
  os.flags(f);
  return os;
}
/* -------------------------------------------------------------------------- */

// \label{newton:LighthillRobot}
struct LighthillRobot {
  double L;
  double w;
  double k;
  double c1;
  double c2;
  int n; // node + 1 (head node is dummy)

  LighthillRobot(double L, double w, double k, double c1, double c2, int n) : L(L), w(w), k(k), c1(c1), c2(c2), n(n + 1) {};

  auto yLH(const double x, const double t) { return (c1 * x / L + c2 * std::pow(x / L, 2)) * std::sin(k * (x / L) - w * t); };

  auto X_RB(const std::array<double, 2> &a, const double q) {
    double r = L / n;
    return a + r * std::array<double, 2>{std::cos(q), std::sin(q)};
  };

  auto f(const std::array<double, 2> &a, const double q, const double t) {
    auto [x, y] = X_RB(a, q);
    return yLH(x, t) - y;
  };

  auto ddx_yLH(const double x, const double t) {
    const double a = k * (x / L) - w * t;
    return (c1 / L + 2 * c2 * x / std::pow(L, 2)) * std::sin(a) + (c1 / L * x + c2 * std::pow(x / L, 2)) * std::cos(a) * k / L;
  };

  auto ddq_f(const double q, const double t) {
    const double r = L / n;
    const double x = r * std::cos(q);
    return -r * std::sin(q) * ddx_yLH(x, t) - r * std::cos(q);
  };

  V_d getAngles(const double t) {
    std::vector<double> Q(n, 0.);
    std::array<double, 2> a{{0., 0.}};
    double error = 0, F, scale = 0.3; //\label{LighthillRobot:scale}
    NewtonRaphson nr(0.);
    for (auto i = 0; i < Q.size(); i++) {
      nr.initialize(std::atan(ddx_yLH(std::get<0>(a), t)));
      error = 0;
      for (auto k = 0; k < 100; k++) {
        F = f(a, nr.X, t);
        nr.update(F * F * 0.5, F * ddq_f(nr.X, t), scale);
        if ((error = std::abs(F)) < 1E-10)
          break;
      }
      Q[i] = nr.X;
      a = X_RB(a, Q[i]);
    }
    return Q;
  };

  std::vector<std::array<double, 2>> anglesToX(const V_d &Q) {
    std::array<double, 2> a = {0., 0.};
    std::vector<std::array<double, 2>> ret;
    ret.reserve(Q.size() + 1);
    ret.push_back(a);
    for (auto i = 0; i < Q.size(); i++)
      ret.push_back(a = X_RB(a, Q[i]));
    return ret;
  };
};
