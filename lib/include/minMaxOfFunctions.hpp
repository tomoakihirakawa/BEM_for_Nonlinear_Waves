#ifndef minMaxOfFunctions_H
#define minMaxOfFunctions_H
#pragma once

#include "basic_linear_systems.hpp"
#include <functional>
// #include "svd.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

/* -------------------------------------------------------------------------- */
/*                               Gradient Method                              */
/* -------------------------------------------------------------------------- */

struct GradientMethod {
public:
  VV_d A;
  VV_d org_A;
  int count;

  bool use_diagonal_scaling;
  V_d diagonal_scaling_vector;
  GradientMethod(const VV_d &A_IN) : org_A(A_IN), A(A_IN), use_diagonal_scaling(false) {};

  void diagonal_scaling() {
    this->use_diagonal_scaling = true;
    this->diagonal_scaling_vector = V_d(org_A.size());
    //
    int s = org_A.size();
    for (auto i = 0; i < s; ++i)
      this->diagonal_scaling_vector[i] = 1. / org_A[i][i];
    for (auto i = 0; i < s; ++i)
      this->A[i] = this->org_A[i] * this->diagonal_scaling_vector[i];
  };

  V_d solve(V_d b, const V_d &x_init = {}, double eps = 1E-10) {
    V_d x(b.size());
    if (use_diagonal_scaling)
      b *= this->diagonal_scaling_vector;

    for (auto i = 0; i < b.size(); ++i) {
      if (i >= x_init.size())
        x[i] = 0.;
      else if (!isFinite(x_init[i]))
        x[i] = 0.;
      else
        x[i] = x_init[i];
    }
    count = 0;
    double norm;
    V_d p = b - Dot(A, x); // 修正ベクトルp
    norm = Norm(p);
    while (!(norm < eps)) {
      //  std::cout << "count = " << count << ", norm = " << norm << std::endl;
      p = b - Dot(A, x += Dot(p, p) / Dot(Dot(A, p), p) * p);
      if (!isFinite(norm = Norm(p)))
        return x;
      else if (count++ > 10000)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    };
    //   std::cout << "count = " << count << ", norm = " << norm << std::endl;
    return x;
  };
  V_d solve(V_d b, double eps) { return solve(b, V_d{}, eps); };
  /* ------------------------------------------------------ */
  V_d solveCG(V_d b, const V_d &x_init = {}, double eps = 1E-10) {
    if (use_diagonal_scaling)
      b *= this->diagonal_scaling_vector;
    V_d x(b.size());
    for (auto i = 0; i < b.size(); i++) {
      if (i >= x_init.size())
        x[i] = 0.;
      else if (!isFinite(x_init[i]))
        x[i] = 0.;
      else
        x[i] = x_init[i];
    }
    count = 0;
    V_d r = b - Dot(A, x); // 残差ベクトル
    double norm = Norm(r);
    //   std::cout << "count = " << count << ", norm = " << norm << std::endl;
    if (norm < eps)
      return x;
    //
    V_d p = r;
    //
    double alpha = Dot(r, p) / Dot(Dot(A, p), p);
    x += alpha * p; // xを修正
    r = b - Dot(this->A, x /*x=x0+alpha*p0*/);
    // check
    double beta;
    norm = Norm(r);
    while (!(norm < eps)) {
      //  std::cout << "count = " << count << ", norm = " << norm << std::endl;
      beta = -Dot(r, Dot(A, p)) / Dot(Dot(A, p), p);
      p = r + beta * p; // 修正ベクトルは，最急降下方向を少し変更した物
      alpha = Dot(r, p) / Dot(Dot(A, p), p);
      x += alpha * p;                      // xを修正
      r = b - Dot(A, x /*x=x0+alpha*p0*/); // 残差ベクトルr（最急降下方向）
      if (!isFinite(norm = Norm(r)))
        return x;
      else if (count++ > 1000)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    };
    //   std::cout << "count = " << count << ", norm = " << norm << std::endl;
    return x;
  };
  V_d solveCG(V_d b, double eps) { return solveCG(b, V_d{}, eps); };

  /* ------------------------------------------------------ */
  GradientMethod() {}
  void initialize(const V_d &Xin) {
    X = Xin;
    isFirst = true;
  }
  V_d X, dX, dFdX_last, s, s_last;
  bool isFirst;

  V_d update(const V_d &dFdX /*direction*/, const double a = 1.) { return X -= (dX = a * dFdX); };

  V_d getFletcherReevesCD(const V_d &dFdX) {
    auto dX = -dFdX;
    auto dX_last = -dFdX_last;
    double beta = Dot(dX, dX) / Dot(dX_last, dX_last);
    this->dFdX_last = dFdX;
    if (isFirst) {
      s_last = V_d(dFdX.size(), 0.);
      isFirst = false;
    }
    return s = s_last = -dFdX + beta * s_last;
  };

  V_d getPolakRibiereCD(const V_d &dFdX) {
    auto dX = -dFdX;
    auto dX_last = -dFdX_last;
    double beta = Dot(dX, dX - dX_last) / Dot(dX_last, dX_last);
    this->dFdX_last = dFdX;
    if (isFirst) {
      s_last = V_d(dFdX.size(), 0.);
      isFirst = false;
    }
    return s = s_last = -dFdX + beta * s_last;
  };

  // V_d update(const V_d &F, const V_d &dF_X) {
  //    V_d X_last = X, m(X.size(), 0.);
  //    double min = 1E+10;
  //    int size = 50;
  //    double d, f, S, extend = 2;
  //    for (auto k = 0; k < size; ++k) {
  //       d = 2. * extend * step / (size - 1);
  //       S = d * k - extend * step;
  //       f = F(X_last + S * m);
  //       if (k == 0 || min >= f) {
  //          min = f;
  //          step = S;
  //       }
  //    }
  //    X = X_last + step * m;
  //    if (this->isFirst) {
  //       dF_X_last = dF_X;
  //       m = dF_X;
  //    } else {
  //       double a = -Dot(dF_X, dF_X - dF_X_last) / Dot(dF_X_last, dF_X_last);
  //       m = dF_X + a * m;
  //    }
  //    dF_X_last = dF_X;
  //    return X;
  // };

  V_d minimize(const std::function<double(V_d)> &F, const std::function<V_d(V_d)> &dF, V_d Xinit, double initial_step = 1E-12, int max_loop = 100) {
    V_d X = Xinit, X_last = Xinit, m(X.size(), 0.);
    double s = initial_step;
    for (auto i = 0; i < max_loop; ++i) {
      X_last = X;
      double min = 1E+10;
      int size = 50;
      double d, f, S, extend = 2;
      for (auto k = 0; k < size; ++k) {
        d = 2. * extend * s / (size - 1);
        S = d * k - extend * s;
        f = F(X_last + S * m);
        if (k == 0 || min >= f) {
          min = f;
          s = S;
        }
      }
      X = X_last + s * m;
      /* ------------------------------------------------------ */
      auto dF_X = dF(X);
      auto dF_Xlast = dF(X_last);
      double a = -Dot(dF_X, dF_X - dF_Xlast) / Dot(dF_Xlast, dF_Xlast);
      // double a = -Dot(dF_X, dF_X - dF_Xlast) / Dot(m, dF_X - dF_Xlast);
      /* ------------------------------------------------------ */
      m = dF(X) + a * m;
      // std::cout << "s = " << s << std::endl;
      std::cout << "i = " << i << std::endl;
      std::cout << "X = " << X << std::endl;
      std::cout << "F(X) = " << F(X) << std::endl;
      std::cout << "Norm(m) = " << Norm(m) << std::endl;
      if (Norm(m) < 1E-8)
        break;
    };
    return X;
  };
};

/* -------------------------------------------------------------------------- */
/*                               Broyden Method                               */
/* -------------------------------------------------------------------------- */

/*DOC_EXTRACT Broydens_Method

準ニュートン法は，ヤコビ行列を近似し，ニュートン法を適用する．近似ヤコビ行列を計算する指針として，セカント条件を満たすようにする．

$$
{\bf J}_k \cdot \Delta {\bf x}_k = \Delta {\bf f}_k
$$

また，

$$
{\bf J}_k = {\bf J}_{k-1} + \frac{(\Delta {\bf f}_k - {\bf J}_{k-1} \cdot \Delta {\bf x}_k) \otimes \Delta {\bf x}_k}{\Delta {\bf x}_k \cdot \Delta {\bf x}_k}
$$

$\Delta {\bf f}_k = {\bf f}_k - {\bf f}_{k-1}$．

*/

// template <typename T>
// struct BroydenMethod;

// template <typename T>
//    requires is_std_array<T>::value
// struct BroydenMethod<T> {
//    T X, dX;
//    std::array<T, std::tuple_size<T>::value> J, Inv_J;

//    // Default constructor
//    BroydenMethod() = default;

//    BroydenMethod(const T &Xin, const T &Xin_) : X(Xin), dX(Xin_ - Xin), J(TensorProduct(Xin, Xin)) {
//       IdentityMatrix(J);
//       Inv_J = J;
//    }

//    void initialize(const T &Xin, const T &dXin) {
//       X = Xin;
//       dX = dXin;
//       J = TensorProduct(Xin, Xin);
//       IdentityMatrix(J);
//       Inv_J = J;
//    }

//    void initialize(const T &Xin) {
//       X = Xin;
//       J = TensorProduct(Xin, Xin);
//       IdentityMatrix(J);
//       Inv_J = J;
//    }

//    void updateGoodBroyden(const T &F, const T &F_, const double alpha = 1.) {
//       auto dot = Dot(dX, dX);
//       auto dF = F - F_;
//       /* -------------------------------------------------------------------------- */
//       // if (dot != static_cast<double>(0.))
//       //    J += TensorProduct((dF - Dot(J, dX)), dX) / dot;
//       // // use lapack_lu
//       // lapack_lu(J, dX, -alpha * F);
//       // X += dX;
//       // X += (dX = -alpha * Dot(Inverse(J), F));  // inverseが計算できるかは未検証
//       /* -------------------------------------------------------------------------- */
//       // good Broyden's method
//       if (Dot(dF, dF) != static_cast<double>(0.))
//          Inv_J += TensorProduct((dX - Dot(Inv_J, dF)), Dot(dX, Inv_J)) / Dot(dX, Dot(Inv_J, dF));
//       X += (dX = -alpha * Dot(Inv_J, F));
//    }

//    void updateBFGS(const T &F, const T &F_, const double alpha = 1.) {
//       auto dot = Dot(dX, dX);
//       auto dF = F - F_;
//       auto s = dX;
//       auto y = dF;
//       double sy = Dot(s, y);
//       if (sy != static_cast<double>(0.))
//          J += TensorProduct(y, y) / Dot(y, s) - Dot(TensorProduct(Dot(J, s), s), Transpose(J)) / Dot(s, Dot(J, s));
//       // X += (dX = -alpha * Dot(Inv_J, F));
//       // Solve(J, dX, -alpha * F);
//       lapack_lu(J, dX, -alpha * F);
//       // lapack_svd_solve(J, dX, -alpha * F);

//       X += dX;
//    }
// };

// // Deduction guide
// template <typename T1, typename T2>
// BroydenMethod(T1, T2) -> BroydenMethod<std::decay_t<T1>>;

using Vd = std::vector<double>;
using VVd = std::vector<std::vector<double>>;

// template <>
// struct BroydenMethod<Vd> {
//    Vd X, dX;
//    VVd J, Inv_J;

//    // Default constructor
//    BroydenMethod() = default;

//    BroydenMethod(const Vd &Xin, const Vd &Xin_) : X(Xin), dX(Xin_ - Xin), J(VVd(Xin.size(), Vd(Xin.size()))) {
//       IdentityMatrix(J);
//       Inv_J = J;
//    }

//    void initialize(const Vd &Xin) { X = Xin; }

//    // void update(const Vd &F, const Vd &F_, const double alpha = 1.) {
//    //    auto dot = Dot(dX, dX);
//    //    if (dot != static_cast<double>(0.))
//    //       J += TensorProduct((F - F_ - Dot(J, dX)), dX) / dot;

//    //    auto dF = F - F_;
//    //    // good Broyden's method
//    //    if (Dot(dF, dF) != static_cast<double>(0.))
//    //       Inv_J += TensorProduct((dX - Dot(Inv_J, dF)), Dot(dX, Inv_J)) / Dot(dX, Dot(Inv_J, dF));

//    //    X += (dX = -alpha * Dot(Inv_J, F));
//    //    // X += (dX = -alpha * Dot(Inverse(J), F));
//    // }

//    void update(const Vd &F, const Vd &F_, const double alpha = 1.) {
//       // auto dot = Dot(dX, dX);
//       // auto dF = F - F_;

//       // auto s = dX;
//       // auto y = dF;
//       // double sy = Dot(s, y);
//       // if (sy != static_cast<double>(0.))
//       //    J += TensorProduct(y, y) / Dot(y, s) - Dot(TensorProduct(Dot(J, s), s), Transpose(J)) / Dot(s, Dot(J, s));

//       // if (sy != static_cast<double>(0.))
//       //    Inv_J += Dot(Dot(s, y) + Dot(Inv_J, y), TensorProduct(s, s)) / std::pow(Dot(s, y), 2) - TensorProduct(Dot(Inv_J, y), s) + Dot(TensorProduct(s, y), Inv_J) / Dot(s, y);

//       // // X += (dX = -alpha * Dot(Inv_J, F));
//       // X += (dX = -alpha * Dot(Inverse(J), F));
//       /* -------------------------------------------------------------------------- */
//       auto dot = Dot(dX, dX);
//       auto dF = F - F_;
//       // good Broyden's method
//       if (Dot(dF, dF) != static_cast<double>(0.))
//          Inv_J += TensorProduct((dX - Dot(Inv_J, dF)), Dot(dX, Inv_J)) / Dot(dX, Dot(Inv_J, dF));
//       X += (dX = -alpha * Dot(Inv_J, F));
//    }
// };

#include <array>
#include <type_traits> // for std::decay_t
#include <vector>

// (Dot, TensorProduct, IdentityMatrixなどのヘルパー関数は、
// この後 templated になっている必要があります)

// プライマリテンプレート（中身は定義しない）
template <typename T> struct vector_traits;

// std::array に対する特殊化
template <typename T_val, std::size_t N> struct vector_traits<std::array<T_val, N>> {
  using value_type = T_val;
  // ベクトル型が std::array<T, N> なら、行列型は std::array<std::array<T, N>, N>
  using matrix_type = std::array<std::array<T_val, N>, N>;

  // サイズを取得する静的関数
  static constexpr std::size_t size(const std::array<T_val, N> &) { return N; }
};

// std::vector に対する特殊化
template <typename T_val> struct vector_traits<std::vector<T_val>> {
  using value_type = T_val;
  // ベクトル型が std::vector<T> なら、行列型は std::vector<std::vector<T>>
  using matrix_type = std::vector<std::vector<T_val>>;

  // サイズを取得する静的関数
  static std::size_t size(const std::vector<T_val> &vec) { return vec.size(); }
};

template <typename VectorType> struct BroydenMethod {
  // vector_traits を使って型を定義
  using traits = vector_traits<VectorType>;
  using MatrixType = typename traits::matrix_type;
  using value_type = typename traits::value_type;

  VectorType X, dX;
  MatrixType J, Inv_J;
  MatrixType I;

  // Default constructor
  BroydenMethod() = default;

  // 統一されたコンストラクタ
  BroydenMethod(const VectorType &Xin, const VectorType &Xin_plus_dX) {
    initialize(Xin, Xin_plus_dX - Xin); // 差分をdXとして渡す
  }

  // 初期化関数も統一
  void initialize(const VectorType &Xin, const VectorType &dXin) {
    X = Xin;
    dX = dXin;
    const auto dim = traits::size(Xin);

    // if constexpr を使って、型に応じた行列の初期化を行う
    if constexpr (std::is_same_v<VectorType, std::vector<value_type>>) {
      // vector の場合はリサイズが必要
      J.resize(dim, std::vector<value_type>(dim));
      Inv_J.resize(dim, std::vector<value_type>(dim));
    }
    // std::array は何もしなくても正しいサイズで生成される

    IdentityMatrix(J); // IdentityMatrixはテンプレート化されている必要がある
    I = Inv_J = J;
  }

  // --- 更新メソッド (array版からコピーしてきて、型を一般化する) ---

  void updateGoodBroyden(const VectorType &F, const VectorType &F_, const double alpha = 1.) {
    auto dF = F - F_;
    double denominator = Dot(dX, Dot(Inv_J, dF));
    if (std::abs(denominator) > 1e-20)
      Inv_J += TensorProduct((dX - Dot(Inv_J, dF)), Dot(dX, Inv_J)) / denominator;
    X += (dX = -alpha * Dot(Inv_J, F));
  }

  void updateBFGS(const VectorType &F, const VectorType &F_, const double alpha = 1.) {
    auto dot = Dot(dX, dX);
    auto dF = F - F_;
    auto s = dX;
    auto y = dF;
    double sy = Dot(s, y);
    if (sy != static_cast<double>(0.))
      J += TensorProduct(y, y) / Dot(y, s) - Dot(TensorProduct(Dot(J, s), s), Transpose(J)) / Dot(s, Dot(J, s));
    // X += (dX = -alpha * Dot(Inv_J, F));
    // Solve(J, dX, -alpha * F);
    // lapack_lu(J, dX, -alpha * F);
    lapack_svd_solve(J, dX, -alpha * F);

    X += dX;
  }

  // 推奨される効率的な逆BFGS法の実装
  void updateInverseBFGS(const VectorType &F, const VectorType &F_, const double alpha = 1.) {
    auto dF = F - F_;
    auto s = dX;
    auto y = dF;
    double s_dot_y = Dot(s, y);

    if (std::abs(s_dot_y) < 1e-12) {
      X += (dX = -alpha * Dot(Inv_J, F));
      return;
    }

    auto s_outer_y = TensorProduct(s, y);

    auto term1 = I - s_outer_y / s_dot_y;
    Inv_J = Dot(term1, Dot(Inv_J, Transpose(term1))); // Dot, Transposeは行列にも対応している必要がある

    const double MIN_CURVATURE = 1e-12;
    if (std::abs(s_dot_y) <= MIN_CURVATURE)
      Inv_J += 0.1 * TensorProduct(s, s) / s_dot_y;
    X += (dX = -alpha * Dot(Inv_J, F));
  }
};

// Deduction guide (これはそのまま使える)
template <typename T1, typename T2> BroydenMethod(T1, T2) -> BroydenMethod<std::decay_t<T1>>;

/* -------------------------------------------------------------------------- */
/*                            Anderson Acceleration                           */
/* -------------------------------------------------------------------------- */

template <typename VectorType> struct AndersonAcceleration {
  using traits = vector_traits<VectorType>;
  using value_type = typename traits::value_type;

  int m_max; // Max history size
  int iter;
  std::vector<VectorType> X_hist; // History of X
  std::vector<VectorType> F_hist; // History of Residuals F(X) = g(X) - X

  AndersonAcceleration(int history_size = 5) : m_max(history_size), iter(0) {}

  void reset() {
    iter = 0;
    X_hist.clear();
    F_hist.clear();
  }

  // Returns the next guess for X
  VectorType compute_next(const VectorType &X, const VectorType &F) {
    X_hist.push_back(X);
    F_hist.push_back(F);

    if (X_hist.size() > m_max + 1) {
      X_hist.erase(X_hist.begin());
      F_hist.erase(F_hist.begin());
    }

    int m_k = X_hist.size() - 1;
    if (m_k == 0) {
      return X + F;
    }

    // Build Least Squares problem to find gamma
    // Minimize || f_k - sum gamma_i (f_{k-j} - f_{k-j-1}) ||
    // Let dF_i = F_curr - F_prev[i] (diff from current to history i)

    std::vector<VectorType> dF(m_k);
    std::vector<VectorType> dX(m_k);

    auto &F_curr = F_hist.back();
    auto &X_curr = X_hist.back();

    for (int i = 0; i < m_k; ++i) {
      int idx_prev = X_hist.size() - 2 - i;
      dF[i] = F_curr - F_hist[idx_prev];
      dX[i] = X_curr - X_hist[idx_prev];
    }

    // Solve LS: min || F_curr - sum gamma_i dF_i ||
    // Normal equations: (dF^T dF) gamma = dF^T F_curr

    VV_d A(m_k, V_d(m_k));
    V_d b(m_k);

    for (int i = 0; i < m_k; ++i) {
      b[i] = Dot(dF[i], F_curr);
      for (int j = i; j < m_k; ++j) {
        double val = Dot(dF[i], dF[j]);
        A[i][j] = val;
        A[j][i] = val;
      }
    }

    V_d gamma(m_k);
    try {
      ludcmp_parallel solver(A);
      solver.solve(b, gamma);
    } catch (...) {
      // Fallback if singular
      return X + F;
    }

    VectorType X_next = X_curr + F_curr;
    for (int i = 0; i < m_k; ++i) {
      X_next -= gamma[i] * (dX[i] + dF[i]);
    }

    return X_next;
  }
};

/* -------------------------------------------------------------------------- */
/*                            Aitken Acceleration                             */
/* -------------------------------------------------------------------------- */

template <typename VectorType> struct AitkenAcceleration {
  using traits = vector_traits<VectorType>;
  using value_type = typename traits::value_type;

  double omega;
  VectorType X_prev;
  VectorType F_prev;
  bool is_first;

  AitkenAcceleration(double omega_init = 1.0) : omega(omega_init), is_first(true) {}

  void reset() {
    is_first = true;
    omega = 1.0;
  }

  // Returns the next guess for X
  VectorType compute_next(const VectorType &X, const VectorType &F) {
    if (is_first) {
      X_prev = X;
      F_prev = F;
      is_first = false;
      return X + omega * F;
    }

    auto dX = X - X_prev; // Delta x_{k-1}
    auto dF = F - F_prev; // Delta f_k
    double dF_norm2 = Dot(dF, dF);

    if (dF_norm2 > 1e-20) {
      // omega_k = -omega_{k-1} * (Delta x_{k-1} . Delta f_k) / |Delta f_k|^2
      omega = -omega * Dot(dX, dF) / dF_norm2;
    }

    // std::cout << "Aitken omega: " << omega << std::endl;

    X_prev = X;
    F_prev = F;

    return X + omega * F;
  }
};

#endif