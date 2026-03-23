#ifndef interpolations_H
#define interpolations_H

#include "basic.hpp"
#include <cmath>
#include <functional>
#include <numeric> //transform_reduce
#include <type_traits>
#include <vector>
// #include "svd.hpp"
// 必要な微分，ラプラシアンなどを定義し，SPHを計算してみる．

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

/* -------------------------------- B-spline -------------------------------- */

template <typename T> struct InterpolationBspline {
  bool is_set = false;
  int K;
  V_d abscissas, sample, a, q;
  InterpolationBspline() {}
  InterpolationBspline(int K_IN, V_d abscissas, std::vector<T> sample) : K(K_IN), abscissas(abscissas) {}
};

template <> struct InterpolationBspline<double> {
  bool is_set = false;
  int K;
  V_d abscissas, sample, a, q;

  InterpolationBspline() {}

  InterpolationBspline(int K_IN, V_d abscissas, V_d sample) : K(K_IN), abscissas(abscissas), sample(sample) { initialize(); }

  void set(int K_IN, V_d abscissas, V_d sample) {
    this->K = K_IN;
    this->abscissas = abscissas;
    this->sample = sample;
    initialize();
  }

private:
  void initialize() {
    q = OpenUniformKnots(abscissas, K);
    auto A = Bspline(abscissas, q, K);
    a.resize(A.size());
    lapack_lu lu(A);
    lu.solve(sample, a);
    this->is_set = true;
  }

public:
  double operator()(const double x) const { return Dot(Bspline(x, q, K), a); }
  double D(const double x) const { return Dot(D_Bspline(x, q, K), a); }

  std::vector<double> DN(const double x) { return D_Bspline(x, q, K); };
};

template <std::size_t N> struct InterpolationBspline<std::array<double, N>> {
  bool is_set = false;
  int K;
  V_d abscissas;
  std::vector<std::array<double, N>> sample;
  std::array<V_d, N> sampleT, a;
  V_d q;

  InterpolationBspline() {}

  InterpolationBspline(int K_IN, V_d abscissas, std::vector<std::array<double, N>> sample_IN) : K(K_IN), abscissas(abscissas), sample(sample_IN), sampleT(Transpose(sample_IN)) { initialize(); }

  void set(int K_IN, V_d abscissas, std::vector<std::array<double, N>> sample_IN) {
    this->K = K_IN;
    this->abscissas = abscissas;
    this->sample = sample_IN;
    this->sampleT = Transpose(sample_IN);
    initialize();
  }

private:
  void initialize() {
    // std::cout << "InterpolationBspline::initialize()" << std::endl;
    q = OpenUniformKnots(abscissas, K);
    auto A = Bspline(abscissas, q, K);
    lapack_lu lu(A);
    // std::cout << "InterpolationBspline::lapack_lu" << stdå::endl;
    for (size_t i = 0; i < sampleT.size(); ++i) {
      a[i].resize(A.size());
      lu.solve(sampleT[i], a[i]);
    }
    this->is_set = true;
  }

public:
  std::array<double, N> operator()(const double x) const {
    std::array<double, N> ret;
    for (std::size_t i = 0; i < N; ++i)
      ret[i] = Dot(Bspline(x, q, K), a[i]);
    return ret;
  }

  std::array<double, N> D(const double x) const {
    std::array<double, N> ret;
    for (std::size_t i = 0; i < N; ++i)
      ret[i] = Dot(D_Bspline(x, q, K), a[i]);
    return ret;
  }
};

/* ---------------------- ラグランジュ補間 ---------------------- */

template <typename T> struct InterpolationLagrange {
  std::vector<double> abscissas;
  std::vector<T> values;
  std::vector<double> denominator;

  InterpolationLagrange() {};

  InterpolationLagrange(const std::vector<double> abscissas) : abscissas(abscissas) {};

  InterpolationLagrange(const std::vector<double> abscissas, const std::vector<T> values) : abscissas(abscissas), values(values) {
    if (abscissas.size() != values.size()) {
      throw std::invalid_argument("Size of abscissas and values vectors must be the same");
    }
    this->set();
  };

  void pop() {
    abscissas.erase(abscissas.begin());
    values.erase(values.begin());
    this->set();
  };

  void push(const double abscissa, const T value) {
    abscissas.push_back(abscissa);
    values.push_back(value);
    this->set();
  };

  // size
  int size() const { return abscissas.size(); };

  void set() {
    denominator.resize(abscissas.size(), 1.);
    for (auto i = 0; i < abscissas.size(); ++i)
      for (auto j = 0; j < abscissas.size(); ++j)
        if (i != j)
          denominator[i] *= (abscissas[i] - abscissas[j]);
  };

  T operator()(const double x) {
    T ret;
    double N = 1;
    ret *= 0.;
    for (auto i = 0; i < abscissas.size(); ++i) {
      for (auto j = 0; j < abscissas.size(); ++j) {
        if (i != j)
          N *= (x - this->abscissas[j]) / (this->abscissas[i] - this->abscissas[j]);
      }
      ret += N * this->values[i];
      N = 1;
    }
    return ret;
  };

  T integrate(const double x_beg, const double x_end) {

  };

  T D(const double x) {
    T ret;
    ret *= 0.;
    for (auto i = 0; i < abscissas.size(); ++i) {
      for (auto j = 0; j < abscissas.size(); ++j) {
        if (i != j) {
          T temp = this->values[i] / (this->abscissas[i] - this->abscissas[j]);
          for (auto k = 0; k < abscissas.size(); ++k) {
            if (k != i && k != j) {
              temp *= (x - this->abscissas[k]) / (this->abscissas[i] - this->abscissas[k]);
            }
          }
          ret += temp;
        }
      }
    }
    return ret;
  }

  std::vector<double> N(const double x) {
    std::vector<double> ret(abscissas.size(), 1.);
    for (auto i = 0; i < abscissas.size(); ++i) {
      for (auto j = 0; j < abscissas.size(); ++j)
        if (i != j)
          ret[i] *= (x - this->abscissas[j]) / (this->abscissas[i] - this->abscissas[j]);
      // ret[i] *= this->values[i];
    }
    return ret;
  };

  std::vector<double> DN(const double x) {
    std::vector<double> ret(abscissas.size(), 0.);
    for (auto i = 0; i < abscissas.size(); ++i) {
      for (auto j = 0; j < abscissas.size(); ++j) {
        if (i != j) {
          T temp = 1. / (this->abscissas[i] - this->abscissas[j]);
          for (auto k = 0; k < abscissas.size(); ++k) {
            if (k != i && k != j) {
              temp *= (x - this->abscissas[k]) / (this->abscissas[i] - this->abscissas[k]);
            }
          }
          ret[i] += temp;
        }
      }
    }
    return ret;
  };
};

#include "kernelFunctions.hpp"
//* ------------------------------------------------------ */
//* ------------------------------------------------------ */
//*          Radial Basis Function Interpolations          */
//* ------------------------------------------------------ */
//* ------------------------------------------------------ */

/*DOC_EXTRACT RBF

## 放射関数補間

距離$r=\left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\|$を引数とする
放射基底関数$\phi(r_i)$に重み$w_i$を掛け合わせて構築した
補間関数$f\left( \mathbf{x} \right)=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\phi \left( \left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\| \right)}$
を放射関数補間という．

### 重み$w_i$の見積もり

重み$w_i$の決定には，サンプル点$A=\left\{ {{\mathbf{a}}_{0}},{{\mathbf{a}}_{1}},...,{{\mathbf{a}}_{N-1}} \right\}$
における値$Y=\left\{ {{y}_{0}},{{y}_{1}},...,{{y}_{N-1}} \right\}$
を使い，補間関数$f$も各サンプル点$A$において値$Y$となる方程式を$w_i$について解く：

$$
\left( \begin{matrix}
   {{w}_{0}}  \\
   \vdots   \\
   {{w}_{N-1}}  \\
\end{matrix} \right)={{\left( \begin{matrix}
   \phi \left( \left\| {{\mathbf{a}}_{0}}-{{\mathbf{a}}_{0}} \right\| \right) & \cdots  & \phi \left( \left\| {{\mathbf{a}}_{0}}-{{\mathbf{a}}_{N-1}} \right\| \right)  \\
   \vdots  & \ddots  & \vdots   \\
   \phi \left( \left\| {{\mathbf{a}}_{N-1}}-{{\mathbf{a}}_{0}} \right\| \right) & \cdots  & \phi \left( \left\| {{\mathbf{a}}_{N-1}}-{{\mathbf{a}}_{N-1}} \right\| \right)  \\
\end{matrix} \right)}^{-1}}\left( \begin{matrix}
   {{y}_{0}}  \\
   \vdots   \\
   {{y}_{N-1}}  \\
\end{matrix} \right)
$$

### 放射基底関数$\phi$

#### 多重二乗（multiquadric RBF）

放射基底関数として多重二乗（multiquadric），
$\phi \left( r \right)={{\left( {{\left( \varepsilon r \right)}^{2}}+1 \right)}^{\frac{1}{2}}}$
がよく使われる．

#### 逆多重二乗（inverse multiquadric RBF）

$\phi \left( r \right)={{\left( {{\left( \varepsilon r \right)}^{2}}+1 \right)}^{-\frac{1}{2}}}$

#### ガウシアン（Gaussian RBF）

$\phi \left( r \right)={{e}^{-{{\left( \varepsilon r \right)}^{2}}}}$

### 補間関数の微分

放射関数補間の微分を少し変形すると，

$\nabla f\left( \mathbf{x} \right)=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\nabla \phi \left( \left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\| \right)}=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\nabla {{r}_{i}}\frac{\partial \phi \left( {{r}_{i}} \right)}{\partial {{r}_{i}}}}$

さらに，計算すると，

$$
\begin{align}
  & {{r}_{i}}=\left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\|={{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{1/2}} \\
 & \frac{\partial {{r}_{i}}}{\partial {{\mathbf{x}}_{k}}}=\frac{1}{2}{{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{-\frac{1}{2}}}\left( \frac{\partial }{\partial {{\mathbf{x}}_{k}}}\sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right) \\
 & =\frac{1}{2}{{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{-\frac{1}{2}}}\left( \sum\limits_{j=0}^{M=2}{2\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}\cdot {{\mathbf{e}}_{k}} \right) \\
 & ={{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{-\frac{1}{2}}}\overbrace{\left( {{\mathbf{x}}_{k}}-{{\mathbf{a}}_{ik}} \right)}^{\text{scaler}}=\frac{\overbrace{\left( {{\mathbf{x}}_{k}}-{{\mathbf{a}}_{ik}} \right)}^{\text{scaler}}}{{{r}_{i}}}
\end{align}
$$

なので，$\nabla {{r}_{i}}=\overbrace{\left( \mathbf{x}-{{\mathbf{a}}_{i}} \right)}^{\text{vecotr}}/{{r}_{i}}$であり，

$\nabla f\left( \mathbf{x} \right)=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\frac{\mathbf{x}-{{\mathbf{a}}_{i}}}{{{r}_{i}}}\frac{\partial \phi \left( {{r}_{i}} \right)}{\partial {{r}_{i}}}}$

である．分母がゼロになる可能性があるが，放射基底関数の微分でキャンセルされる．

#### 多重二乗

$\phi \left( r \right)={{\left( {{\left( \varepsilon r \right)}^{2}}+1 \right)}^{\frac{1}{2}}},\frac{\partial \phi }{\partial r}\left( r \right)=\frac{\varepsilon^2 r}{\phi \left( r \right)}$

なので，次のように分母を消すことができる．

$\nabla f\left( \mathbf{x} \right)=\varepsilon^2 \sum\limits_{i=0}^{N-1}{{{w}_{i}}\frac{\mathbf{x}-{{\mathbf{a}}_{i}}}{\phi \left( {{r}_{i}} \right)}}$

#### 逆多重二乗

$$
\begin{align}
  & \phi \left( r \right)={{\left( {{\left( \varepsilon r \right)}^{2}}+1 \right)}^{-\frac{1}{2}}} \\
 & \frac{\partial \phi }{\partial r}\left( r \right)=-{{\varepsilon }^{2}}r{{\left( {{\left( \varepsilon r \right)}^{2}}+1 \right)}^{-1}} \\
 & \nabla f=\sum\limits_{i=0}^{N-1}{-{{\varepsilon }^{2}}\left( \mathbf{x}-{{\mathbf{a}}_{i}} \right){{\phi }^{2}}\left( r \right)} \\
\end{align}
$$

#### ガウシアン

$$
\begin{align}
  & \phi \left( r \right)={{e}^{-{{\left( \varepsilon r \right)}^{2}}}} \\
 & \frac{\partial \phi }{\partial r}\left( r \right)=-2{{\varepsilon }^{2}}r{{e}^{-{{\left( \varepsilon r \right)}^{2}}}} \\
 & \nabla f=\sum\limits_{i=0}^{N-1}{-2{{\varepsilon }^{2}}{{e}^{-{{\left( \varepsilon r \right)}^{2}}}}\left( \mathbf{x}-{{\mathbf{a}}_{i}} \right)} \\
\end{align}
$$

### 最適なパラメタ${\varepsilon}$

サンプル点の平均的な間隔を${s}$とした場合，${\varepsilon = 1/s}$とパラメタをとるとよい．

*/

class InterpolationVectorRBF {
  // using V_d = std::vector<double>;
  // using VV_d = std::vector<std::vector<double>>;
  // using VVV_d = std::vector<std::vector<std::vector<double>>>;

  /* -------------- 各positionでサンプリングしたスカラーの補間 ------------- */
private:
  VV_d P; // 行 x 列 = givenサンプル数　x パラメタの成分数
  // Pは，サンプルの値がわかっている場所のパラメタ
  VV_d F; // 行 x 列 = givenサンプル数分の式の数　x givenサンプル数
  // Fは，核関数にパラメタを代入して作る行列
  VV_d invF;
  std::function<double(V_d, V_d)> phi;           // RBF basis function passed as a lambda function
  std::function<V_d(V_d, V_d)> grad_phi;         // derivative of the RBF basis function passed as a lambda function
  std::function<double(V_d, V_d)> laplacian_phi; // RBF basis function passed as a lambda function

  int parameter_dim;
  int returning_dim;
  double scale;

public:
  /* -------------- 各positionでサンプリングしたベクトルの補間 ------------- */
private:
  VV_d VV_weights; // 行 x 列 = パラメタの成分数 x givenサンプル数
  // weightは，invF.VV_valuesで求まる
  VV_d VV_values; // values of position, {{nx,ny,nz},{nx,ny,nz},{nx,ny,nz},...}

public:
  InterpolationVectorRBF() { /*これを使う場合setupをかならず行うこと*/ };
  void set(const VV_d &P_IN, const VV_d &VV_IN) {
    this->P = P_IN;
    this->VV_values = VV_IN;
    this->parameter_dim = P_IN[0].size();
    this->returning_dim = VV_IN[0].size();
    //
    if (P_IN.size() != VV_IN.size()) {
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    }
    this->scale = RBFscale(P_IN);
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    //! 関数を代入していからweight計算
    this->VV_weights = weight(P_IN, VV_IN); // VV_INと同じサイズの行列になる
  };
  // デフォルトの補間基底関数，多重二乗を使う場合
  InterpolationVectorRBF(const VV_d &P_IN, const VV_d &VV_IN) : P(P_IN), VV_values(VV_IN), parameter_dim(P_IN[0].size()), returning_dim(VV_IN[0].size()) {
    if (P_IN.size() != VV_IN.size()) {
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    }
    this->scale = RBFscale(P_IN);
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    //! 関数を代入していからweight計算
    this->VV_weights = weight(P_IN, VV_IN); // VV_INと同じサイズの行列になる
  };
  InterpolationVectorRBF(const VV_d &P_IN, const VV_d &VV_IN, const V_d &target) : P(P_IN), VV_values(VV_IN), parameter_dim(P_IN[0].size()), returning_dim(VV_IN[0].size()) {
    if (P_IN.size() != VV_IN.size()) {
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    }
    this->scale = RBFscale(P_IN, target);
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    //! 関数を代入していからweight計算
    this->VV_weights = weight(P_IN, VV_IN); // VV_INと同じサイズの行列になる
  };
  InterpolationVectorRBF(const VV_d &P_IN, const VV_d &VV_IN, const VV_d &P_for_Derivative_IN, const VV_d &VV_for_Derivative_IN, const VV_d &normals_IN) : P(P_IN), VV_values(VV_IN), parameter_dim(P_IN[0].size()), returning_dim(VV_IN[0].size()) {
    if (P_IN.size() != VV_IN.size()) {
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    }
    this->scale = RBFscale(Join(P_IN, P_for_Derivative_IN));
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    //! 関数を代入していからweight計算
    this->VV_weights = weight(P_IN, VV_IN, P_for_Derivative_IN, VV_for_Derivative_IN, normals_IN); // VV_INと同じサイズの行列になる
  };
  /* ------------------------ メソッド ------------------------ */
public:
  V_d operator()(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    V_d ret(this->returning_dim, 0.);
    for (auto i = 0; i < this->VV_weights.size(); i++) {
      /**
       *@ this->VV_weights = {wx,wy,wz}
       *@ P = {samples_x,samples_y, samples_z}
       */
      ret += this->VV_weights[i] * phi(x, this->P[i]);
    }
    return ret;
  };

  // ヤコビアンの計算のため，2021/07/12に追加
  VV_d grad(const V_d &x) const {
    /* 復習
     * f(t) = w_i * f_i(t)
     * -> ddt f(t) = w_i * ddt f_i(t)これを各成分毎に実行すればgradとなる．各成分をt_jとすれば，
     * -> ddt_j f(t_j) = w_ij * ddt_j f_i(t_j)となる．
     */
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    int j = 0;
    VV_d ret(this->parameter_dim, V_d(this->returning_dim, 0.));
    for (auto i = 0; i < this->VV_weights.size(); i++) {
      j = 0;
      for (const auto &dfdp /*各成分毎に*/ : grad_phi(x, this->P[i] /*サンプル点番号iをしてい->パラメタベクトルを返す*/))
        ret[j++] += this->VV_weights[i] * dfdp;
    }
    return ret; // parameter_dim x returning_dimということに注意
  };

  // ヤコビアンの計算のため，2021/07/12に追加
  V_d gradOfDot(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    V_d ret(this->parameter_dim, 0);
    VV_d gd = this->grad(x);
    for (auto i = 0; i < gd.size(); i++)
      ret[i] = 2. * Dot((*this)(x), gd[i]);
    return ret;
  };

  double div(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    double ret = 0.;
    for (auto i = 0; i < this->VV_weights.size(); i++) {
      ret += Dot(this->VV_weights[i], grad_phi(x, this->P[i]));
    }
    return ret;
  };

  //! 各要素に対するラプラシアン
  V_d laplacian(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    V_d ret(this->returning_dim, 0.);
    for (auto i = 0; i < VV_weights.size(); i++)
      ret += this->VV_weights[i] * laplacian_phi(x, P[i]);
    return ret;
  };

  double RBFscale(const VV_d &sample) {
    V_d r;
    for (auto i = 0; i < sample.size(); i++)
      for (auto j = i + 1; j < sample.size(); j++)
        r.push_back(Norm(sample[i] - sample[j]));

    std::sort(r.begin(), r.end(), [](const auto &lhs, const auto &rhs) { return lhs < rhs; });
    auto s = 1 + (int)((double)r.size() / 2.);
    V_d v(s);
    for (auto i = 0; i < s; i++)
      v[i] = r[i];
    return Mean(v);
  };

  double RBFscale(const VV_d &sample, const V_d &target) {
    V_d r;
    double tmp;
    for (auto i = 0; i < sample.size(); i++)
      if ((tmp = Norm(sample[i] - target)) > 1E-10 /*except itself*/)
        r.emplace_back(tmp);
    std::sort(r.begin(), r.end(), [](const auto lhs, const auto rhs) { return lhs < rhs; });
    auto s = 1 + (int)((double)r.size() / 3.);
    V_d v(s);
    for (auto i = 0; i < s; i++)
      v[i] = r[i];
    return Mean(v);
  };

  V_d normV(const V_d &x /*is {x,y,z}*/, const VV_d &P) const {
    V_d ret(P.size());
    std::transform(P.begin(), P.end(), ret.begin(), [this, &x](const auto &a) { return this->phi(x, a); });
    return ret;
  };

  VV_d normM(const VV_d &P) const {
    VV_d R(P.size(), V_d(P.size(), 0));
    std::transform(P.begin(), P.end(), R.begin(), [this, &P](const auto &a) { return normV(a, P); });
    return R;
  };

  // この重みのサイズは？->V_vOFx
  VV_d weight(const VV_d &P, const VV_d &V_vOFx /*{{n0,n1,n2},{m0,m1,m2},{a0,a1,a2},..}*/) {
    try {
      VV_d ret(0);
      this->F = normM(P);
      lapack_lu lu(normM(P) /*未知変数ベクトル側の係数行列*/);
      this->invF = lu.inverse();
      V_d w(P.size());
      for (const auto &vOFx : Transpose(V_vOFx) /*各成分*/) {
        /**
         *@ x,y,z成分毎の重みw
         *@ {wx,wy,wz}={Dot(invF,Ax),Dot(invF,Ay),Dot(invF,Az)}
         */
        lu.solve(vOFx, w);
        ret.emplace_back(w);
      }
      return Transpose(ret); // 各点
    } catch (const error_message &e) {
      throw;
    }
  };
  VV_d weight(const VV_d &P, const VV_d &V_vOFx, const VV_d &P_for_Derivative, const VV_d &V_derivative_OFx, const VV_d &normals /*{{n0,n1,n2},{m0,m1,m2},{a0,a1,a2},..}*/) {
    try {
      VV_d ret(0);
      auto s = P.size(), k = P_for_Derivative.size();
      /*
      <---s-->
      +---------------+ ---> j
      |       |       | phi
      +---------------+
      |       |       | grad_phi
      +---------------+
      */
      // phi
      this->F = VV_d(s + k, V_d(s + k));
      for (auto i = 0; i < P.size(); i++) {
        auto a = P[i];
        auto &f = this->F[i];
        for (auto j = 0; j < P.size(); j++)
          f[j] = this->phi(a, P[j]);
        for (auto j = 0; j < P_for_Derivative.size(); j++)
          f[s + j] = this->phi(a, P_for_Derivative[j]);
      }
      for (auto i = 0; i < P_for_Derivative.size(); i++) {
        auto n = normals[i];
        auto &f = this->F[s + i];
        auto a = P_for_Derivative[s + i];
        for (auto j = 0; j < P.size(); j++)
          f[j] = Dot(n, this->grad_phi(a, P[j]));
        for (auto j = 0; j < P_for_Derivative.size(); j++)
          f[s + j] = Dot(n, this->grad_phi(a, P_for_Derivative[j]));
      }
      //
      lapack_lu lu(this->F /*未知変数ベクトル側の係数行列*/);
      this->invF = lu.inverse();
      V_d w(P.size() + P_for_Derivative.size());
      for (const auto &vOFx : Transpose(Join(V_vOFx, V_derivative_OFx)) /*各成分*/) {
        /**
         *@ x,y,z成分毎の重みw
         *@ {wx,wy,wz}={Dot(invF,Ax),Dot(invF,Ay),Dot(invF,Az)}
         */
        lu.solve(vOFx, w);
        ret.emplace_back(w);
      }
      return Transpose(ret); // 各点
    } catch (const error_message &e) {
      throw;
    }
  };

  /* ------------------------- パラメタを介して使う場合 ------------------------- */
  /*
   * 復習
   * f_ij  w_j = v_i
   * -> w_i = (f^-1)_ij v_j これであっている．
   * -> f(t) = w_i f_i(t)
   * -> f(t) = (f^-1)_ij v_j f_i(t)
   * -> f(t) =  V.((invF^T).F(t))
   * -> f(t) =  V.(F(t).invF) <-これで計算できる
   * Vはgivenでサンプル点なので，F(t).invF=N(t)で形状関数と考えることができる．
   *
   * gradは注意：
   * 間違いなく，
   * ddt(F(t)ベクトル).invF=ddt(N(t))　なので，これを各成分行うことで
   * gradを計算することができる．
   * gradNを使った計算も成分毎に考えよう．gradをベクトルとしてそのまま入れ込むと，間違える．
   * -> ddt(f(t)) =  V.(ddt(F(t).invF))
   */
  V_d N(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    V_d PHI(P.size(), 0.);
    for (auto i = 0; i < P.size(); i++)
      PHI[i] = phi(x, P[i]);
    return Dot(PHI, this->invF);
  };
  V_d N(const double t0, const double t1) const { return this->N({t0, t1}); };

  VV_d gradN(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    VV_d tmp(P.size());
    for (auto i = 0; i < P.size(); i++)
      for (const auto &gp : grad_phi(x, P[i]))
        tmp[i].emplace_back(gp);
    return Dot(Transpose(tmp), this->invF);
  };

  VV_d gradN(const double t0, const double t1) const { return this->gradN({t0, t1}); };

  V_d cross(const double t0, const double t1) const {
    auto dxdt = Dot(this->gradN(t0, t1), this->VV_values);
    return Cross(dxdt[0], dxdt[1]);
  };

  double J(const double t0, const double t1) const {
    auto dxdt = Dot(this->gradN(t0, t1), this->VV_values);
    return std::sqrt(Dot(dxdt[0], dxdt[0]) * Dot(dxdt[1], dxdt[1]) - std::pow(Dot(dxdt[0], dxdt[1]), 2.));
  };
};

struct InterpolationRBF {
  using V_d = std::vector<double>;
  using VV_d = std::vector<std::vector<double>>;
  using VVV_d = std::vector<std::vector<std::vector<double>>>;

  /* -------------- 各positionでサンプリングしたスカラーの補間 ------------- */
  // private:
  V_d w;  // weight
  VV_d A; // position
  V_d V;  // values of position
  VV_d F;
  std::function<double(V_d, V_d)> phi;           // RBF basis function passed as a lambda function
  std::function<V_d(V_d, V_d)> grad_phi;         // derivative of the RBF basis function passed as a lambda function
  std::function<double(V_d, V_d)> laplacian_phi; // derivative of the RBF basis function passed as a lambda function

  int parameter_dim;
  double scale;

public:
  InterpolationRBF(const VV_d &A_IN, const V_d &V_IN, const std::function<double(V_d, V_d)> &phi_IN, const std::function<V_d(V_d, V_d)> &grad_phi_IN) : phi(phi_IN), grad_phi(grad_phi_IN), A(A_IN), V(V_IN), parameter_dim(A_IN[0].size()) {
    this->scale = RBFscale(A_IN);
    this->w = weight(A_IN, V_IN);
  };
//!
#define old_derivative_RBF
#if defined(old_derivative_RBF)
  InterpolationRBF(const VV_d &A_IN, const V_d &V_IN, const V_d &V2_IN, const VV_d &dotter_IN, double scaleIN = -1.) : A(A_IN), V(V_IN), parameter_dim(A_IN[0].size()) {
    if (scaleIN < 0.)
      this->scale = RBFscale(A_IN);
    else
      this->scale = scaleIN;
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    // this->w = weight(A_IN, V_IN);
    //! ------------------ 微分の値を考慮して重み関数を計算 ------------------ */
    VV_d M(2 * A.size(), V_d(A.size()));
    V_d row(A.size());
    for (auto i = 0; i < A.size(); ++i) { //! each row
      for (auto j = 0; j < A.size(); ++j)
        row[j] = (this->phi(A[i], A[j])); //{,,,}
      M[i] = row;
    }
    for (auto i = 0; i < A.size(); ++i) { //! each row
      for (auto j = 0; j < A.size(); ++j)
        row[j] = Dot(dotter_IN[i], this->grad_phi(A[i], A[j])); //{,,,,,,,}
      M[i + A.size()] = row;
    }
    V_d vOFx(0);
    for (const auto &v : V_IN)
      vOFx.emplace_back(v);

    this->F = M;
    this->w.resize(A.size(), 0.);
    lapack_svd lu(M);
    std::cout << M << std::endl;
    lu.solve(vOFx, this->w);
  };
  //@ ------------------------------------------------------ */
  std::vector<std::tuple<Tddd, double, double>> EQ;
  std::vector<std::tuple<Tddd, Tddd, double, double>> EQ_dot;
  InterpolationRBF(std::vector<Tddd> kernel_locations, std::vector<std::tuple<Tddd, double, double>> eq) : parameter_dim(3), EQ(eq) {
    //! ------------------ 微分の値を考慮して重み関数を計算 ------------------ */
    VV_d M(0);
    V_d RES(0);
    for (const auto &eqn : eq) {
      auto [X, res, s] = eqn;
      RES.emplace_back(res);
      V_d row;
      for (const auto &kernelX : kernel_locations)
        row.emplace_back(kernel_MQ(ToVector(X), /*横*/ ToVector(kernelX), 1. / s));
      M.emplace_back(row);
    }
    // Aは各関数の位置で，phi(.,A[i])，gradphi(.,A[i])とかになる
    this->A.clear();
    std::transform(kernel_locations.begin(), kernel_locations.end(), std::back_inserter(this->A), [](const auto &kernelX) { return ToVector(kernelX); });
    this->F = M;
    this->w = RES;
    lapack_lu lu(M);
    lu.solve(RES, this->w);
  };
  //
  double calc(const V_d &x) const {
    double ret(0.);
    for (auto i = 0; i < EQ.size(); ++i) {
      auto [X, res, s] = EQ[i];
      ret += w[i] * kernel_MQ(x, this->A[i], 1. / s);
    }
    return ret;
  };
  //
  V_d calc_grad(const V_d &x) const {
    V_d ret(x.size(), 0.);
    for (auto i = 0; i < EQ.size(); i++) {
      auto [X, res, s] = EQ[i];
      ret += this->w[i] * grad_kernel_MQ(x, this->A[i], 1. / s);
    };
    return ret;
  };
  /* ------------------------------------------------------ */
  /* ------------------------------------------------------ */
  InterpolationRBF(const VV_d &A_IN, const V_d &V_IN, const VV_d &A2_IN, const V_d &V2_IN, const VV_d &dotter_IN, double scaleIN = -1.) : A(A_IN), V(V_IN), parameter_dim(A_IN[0].size()) {
    if (scaleIN < 0.)
      this->scale = RBFscale(A_IN);
    else
      this->scale = scaleIN;
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    // this->w = weight(A_IN, V_IN);
    //! ------------------ 微分の値を考慮して重み関数を計算 ------------------ */
    auto A2 = A2_IN;
    VV_d M(A2_IN.size() + A_IN.size(), V_d(A_IN.size() + A2.size()));
    V_d row(A.size() + A2.size());
    for (auto i = 0; i < A.size(); ++i) { //! each row
      for (auto j = 0; j < A.size(); ++j)
        row[j] = (this->phi(A[i], A[j])); //{,,,}
      for (auto j = 0; j < A2.size(); ++j)
        row[A.size() + j] = (this->phi(A[i], A2[j])); //{,,,,,,,}
      M[i] = row;
    }
    for (auto i = 0; i < A2.size(); ++i) { //! each row
      for (auto j = 0; j < A.size(); ++j)
        row[j] = Dot(dotter_IN[i], this->grad_phi(A2[i], /*横*/ A[j])); //{,,,}
      for (auto j = 0; j < A2.size(); ++j)
        row[A.size() + j] = Dot(dotter_IN[i], this->grad_phi(A2[i], /*横*/ A2[j])); //{,,,,,,,}
      M[A.size() + i] = row;
    }
    // VV_d M(0);
    // auto A2 = A2_IN;
    // for (const auto &x : A)
    // { //!each row
    // 	V_d row;
    // 	for (const auto &a : A)
    // 		row.emplace_back(this->phi(x, a)); //{,,,}
    // 	for (const auto &a : A2)
    // 		row.emplace_back(this->phi(x, a)); //{,,,,,,,}
    // 	M.emplace_back(row);
    // }
    // for (auto i = 0; i < A2.size(); i++)
    // { //!each row
    // 	V_d row;
    // 	V_d d = dotter_IN[i];
    // 	V_d x = A2[i];
    // 	for (const auto &a : A)
    // 		row.emplace_back(Dot(d, this->grad_phi(x, a))); //{,,,}
    // 	for (const auto &a : A2)
    // 		row.emplace_back(Dot(d, this->grad_phi(x, a))); //{,,,,,,,}
    // 	M.emplace_back(row);
    // }

    V_d vOFx(0);
    for (const auto &v : V_IN)
      vOFx.emplace_back(v);
    for (const auto &v : V2_IN)
      vOFx.emplace_back(v);

    // std::cout << "dotter_IN = " << dotter_IN << std::endl;
    // std::cout << "M" << M << std::endl;
    //!
    this->F = M;
    this->w.resize(A.size() + A2.size(), 0.);
    lapack_lu lu(M);
    // SVD lu(M);
    lu.solve(vOFx, this->w);

    this->A = Join(A, A2);
  };
#else
  InterpolationRBF(const VV_d &A_IN, const V_d &V_IN, const VV_d &P_for_Derivative, const V_d &V_derivative_OFx, const VV_d &normals /*{{n0,n1,n2},{m0,m1,m2},{a0,a1,a2},..}*/) : A(A_IN), V(V_IN), parameter_dim(A_IN[0].size()) {
    this->scale = RBFscale(Join(A_IN, P_for_Derivative));
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    this->w = weight(A_IN, V_IN, P_for_Derivative, V_derivative_OFx, normals);
  };
#endif
  // ディフォルトの補間基底関数，多重二乗を使う場合
  InterpolationRBF(const VV_d &A_IN, const V_d &V_IN, const double scaleIN = -1.) : A(A_IN), V(V_IN), parameter_dim(A_IN[0].size()) {
    if (scaleIN < 0.)
      this->scale = RBFscale(A_IN);
    else
      this->scale = scaleIN;
    //////////////////////////////////////
    /////////// multiquadric /////////////
    //////////////////////////////////////
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };
    ////////////////////////////////////
    /////////// thin plate /////////////
    ////////////////////////////////////
    // this->phi = [this](const V_d &x, const V_d &a) {
    //   auto r = Norm(x - a);
    //   auto e = 1. / this->scale;
    //   if (1E-8 > e * r)
    //     return 0.;
    //   else
    //     return e * r * e * r * std::log(e * r);
    // };
    // this->grad_phi = [this](const V_d &x, const V_d &a) {
    //   auto r = Norm(x - a);
    //   auto e = 1. / this->scale;
    //   if (1E-13 > e * r)
    //     return 0. * (x - a);
    //   else
    //     return e * e * (1. + 2. * std::log(e * r)) * (x - a);
    // };
    ///////////////////////////////////
    this->w = weight(A_IN, V_IN);
  };
  //
  InterpolationRBF(const VV_d &A_IN, const V_d &V_IN, const V_d &target /*補間点が決まっている場合*/) : A(A_IN), V(V_IN), parameter_dim(A_IN[0].size()) {
    this->scale = RBFscale(A_IN, target);

    //////////////////////////////////////
    /////////// multiquadric /////////////
    //////////////////////////////////////
    this->phi = [this](const V_d &x, const V_d &a) { return kernel_MQ(x, a, 1. / this->scale); };
    this->grad_phi = [this](const V_d &x, const V_d &a) { return grad_kernel_MQ(x, a, 1. / this->scale); };
    this->laplacian_phi = [this](const V_d &x, const V_d &a) { return laplacian_kernel_MQ(x, a, 1. / this->scale); };

    ////////////////////////////////////
    /////////// thin plate /////////////
    ////////////////////////////////////
    // this->phi = [this](const V_d &x, const V_d &a) {
    //   auto r = Norm(x - a);
    //   auto e = 1. / this->scale;
    //   if (1E-8 > e * r)
    //     return 0.;
    //   else
    //     return e * r * e * r * std::log(e * r);
    // };
    // this->grad_phi = [this](const V_d &x, const V_d &a) {
    //   auto r = Norm(x - a);
    //   auto e = 1. / this->scale;
    //   if (1E-13 > e * r)
    //     return 0. * (x - a);
    //   else
    //     return e * e * (1. + 2. * std::log(e * r)) * (x - a);
    // };
    ///////////////////////////////////
    this->w = weight(A_IN, V_IN);
  };
  //

  /* ------------------------ メソッド ------------------------ */
public:
  double operator()(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    // return std::transform_reduce(
    // 	this->w.cbegin(), this->w.cend(), this->A.cbegin(), 0.,
    // 	[](const double a, const double d)
    // 	{ return a + d; },
    // 	[x, this](const auto &ww, const auto &aa)
    // 	{ return ww * phi(x, aa); });

    double ret(0.);
    for (auto i = 0; i < w.size(); i++)
      ret += w[i] * phi(x, A[i] /*kernel location*/);
    return ret;
  };

  // V_d nabla(const V_d &x) const {
  // 	if (this->parameter_dim != x.size()) {
  // 		std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
  // 		throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
  // 	}
  // 	// return std::transform_reduce(this->w.cbegin(),this->w.cend(),this->A.cbegin(),V_d(x.size(),0.),
  // 	//                              [](const auto& acc, const auto& res){return acc + res;},
  // 	//                              [&x, this](const auto&  w_, const auto& a){return w_ * grad_phi(x,a);});
  // 	//----------
  // 	V_d ret(x.size(), 0.);
  // 	for (auto i = 0; i < w.size(); i++)
  // 		ret += w[i] * grad_phi(x, A[i]);
  // 	return ret;
  // };

  V_d grad(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    V_d ret(x.size(), 0.);
    // return std::transform_reduce(
    // 	this->w.cbegin(), this->w.cend(), this->A.cbegin(), ret,
    // 	[](const auto &a, const auto &d)
    // 	{ return a + d; },
    // 	[x, this](const auto &ww, const auto &aa)
    // 	{ return ww * grad_phi(x, aa); });
    for (auto i = 0; i < w.size(); i++)
      ret += this->w[i] * grad_phi(x, A[i]);
    return ret;
  };

  double laplacian(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }
    // return std::transform_reduce(
    // 	this->w.cbegin(), this->w.cend(), this->A.cbegin(), 0.,
    // 	[](const double a, const double d)
    // 	{ return a + d; },
    // 	[x, this](const auto &ww, const auto &aa)
    // 	{ return ww * laplacian_phi(x, aa); });
    double ret = 0.;
    for (auto i = 0; i < w.size(); i++)
      ret += this->w[i] * laplacian_phi(x, A[i]);
    return ret;
  };

  //////////

  double RBFscale(const VV_d &sample) {
    V_d r;
    for (auto i = 0; i < sample.size(); i++)
      for (auto j = i + 1; j < sample.size(); j++)
        r.push_back(Norm(sample[i] - sample[j]));

    std::sort(r.begin(), r.end(), [](const auto &lhs, const auto &rhs) { return lhs < rhs; });
    // return Mean(V_d(r.begin(), std::next(r.begin(), 1 + (int)((double)r.size() / 2.))));
    auto s = 2 + (int)((double)r.size() / 3.);
    V_d v(s);
    for (auto i = 1; i < s; i++)
      v[i] = r[i];

    return Mean(v);
  };

  double RBFscale(const VV_d &sample, const V_d &target) {
    V_d r;
    double tmp;
    for (auto i = 0; i < sample.size(); i++)
      if ((tmp = Norm(sample[i] - target)) > 1E-10 /*except itself*/)
        r.emplace_back(tmp);

    std::sort(r.begin(), r.end(), [](const auto lhs, const auto rhs) { return lhs < rhs; });

    // return Mean(V_d(r.begin(), std::next(r.begin(), 1 + (int)((double)r.size() / 2.))));
    auto s = 1 + (int)((double)r.size() / 3.);

    V_d v(s);
    for (auto i = 0; i < s; i++)
      v[i] = r[i];

    return Mean(v);
  };

  V_d normV(const V_d &x /*is {x,y,z}*/, const VV_d &A) const {
    V_d ret(A.size());
    std::transform(A.begin(), A.end(), ret.begin(), [this, &x](const auto &a) { return this->phi(x, a); });
    return ret;
  };

  VV_d normM(const VV_d &A) const {
    VV_d R(A.size(), V_d(A.size(), 0));
    std::transform(A.begin(), A.end(), R.begin(), [this, &A](const auto &a) { return normV(a, A); });
    return R;
  };

  V_d weight(const VV_d &A, const V_d &vOFx) const {
    V_d w(A.size());
    lapack_lu lu(normM(A));
    lu.solve(vOFx, w);
    return w;
  };

  V_d weight(const VV_d &P, const V_d &V_vOFx, const VV_d &P_for_Derivative, const V_d &V_derivative_OFx, const VV_d &normals /*{{n0,n1,n2},{m0,m1,m2},{a0,a1,a2},..}*/) {
    try {
      VV_d ret(0);
      auto s = P.size(), k = P_for_Derivative.size();
      /*
      <---s-->
      +---------------+ ---> j
      |       |       | phi
      +---------------+
      |       |       | grad_phi
      +---------------+
      */
      // phi
      auto F = VV_d(s + k, V_d(s + k));
      for (auto i = 0; i < P.size(); i++) {
        auto a = P[i];
        auto &f = F[i];
        for (auto j = 0; j < P.size(); j++)
          f[j] = this->phi(a, P[j]);
        for (auto j = 0; j < P_for_Derivative.size(); j++)
          f[s + j] = this->phi(a, P_for_Derivative[j]);
      }
      for (auto i = 0; i < P_for_Derivative.size(); i++) {
        auto n = normals[i];
        auto &f = F[s + i];
        auto a = P_for_Derivative[s + i];
        for (auto j = 0; j < P.size(); j++)
          f[j] = Dot(n, this->grad_phi(a, P[j]));
        for (auto j = 0; j < P_for_Derivative.size(); j++)
          f[s + j] = Dot(n, this->grad_phi(a, P_for_Derivative[j]));
      }
      //
      lapack_lu lu(F /*未知変数ベクトル側の係数行列*/);
      V_d w(s + k);
      lu.solve(Join(V_vOFx, V_derivative_OFx), w);
      return w;
    } catch (const error_message &e) {
      throw;
    }
  };
};

/*RBF_interp_code*/
//@ ------------------------------------------------------ */
// template <typename T, typename U>
// struct InterpolationRBF_Common
// {
// 	VV_d M;
// 	std::vector<std::tuple<T, double, double>> kernel_X_s_w;
// 	std::vector<std::tuple<T, U>> EQs;
// 	/*
// 		EQs:
// 		カーネルの座標．補間点からカーネルの座標までの距離rを計算するためにつかう
// 		カーネルが通る値
// 	*/
// 	//! ------------------ 微分の値を考慮して重み関数を計算 ------------------ */
// 	InterpolationRBF_Common(const std::vector<std::tuple<T, double>> &kernel_X_s_IN,
// 							const std::vector<std::tuple<T, U>> &EQs_IN)
// 		: kernel_X_s_w(kernel_X_s_IN.size()),
// 		  M(kernel_X_s_IN.size(), V_d(kernel_X_s_IN.size())),
// 		  EQs(EQs_IN)
// 	{
// 		int i = 0;
// 		for (const auto &X_s : kernel_X_s_IN)
// 		{
// 			std::get<0>(kernel_X_s_w[i]) = std::get<0>(X_s);
// 			std::get<1>(kernel_X_s_w[i++]) = std::get<1>(X_s);
// 		}
// 		//
// 		std::vector<U> b(EQs_IN.size());
// 		V_d row(kernel_X_s_IN.size());
// 		for (auto i = 0; i < EQs.size(); ++i)
// 		{
// 			auto [sample_X, res] = EQs[i];
// 			b[i] = res;
// 			for (auto j = 0; j < kernel_X_s_IN.size(); ++j)
// 			{
// 				auto [kernel_X, s] = kernel_X_s_IN[j];
// 				row[j] = kernel_MQ(sample_X, kernel_X, 1. / s);
// 			}
// 			M[i] = row;
// 		}
// 		// Aは各関数の位置で，phi(.,A[i])，gradphi(.,A[i])とかになる
// 		lapack_lu lu(M);
// 		V_d w(kernel_X_s_w.size());
// 		lu.solve(b, w);
// 		for (auto i = 0; i < kernel_X_s_w.size(); ++i)
// 			std::get<2>(kernel_X_s_w[i]) = w[i];
// 	};
// 	//! ------------------------------------------------------ */
// 	U operator()(const T &x) const
// 	{
// 		U ret(0.);
// 		for (const auto &X_s_w : kernel_X_s_w)
// 			ret += std::get<2>(X_s_w) * kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));
// 		return ret;
// 	};
// 	//! ------------------------------------------------------ */
// };
//@ ------------------------------------------------------ */
// template <typename T, typename U>
// struct InterpolationRBF_Common
// {
// 	VV_d M;
// 	std::vector<std::tuple<T, double, double>> kernel_X_s_w;
// 	std::vector<std::tuple<T, U>> EQs;
// 	/*
// 		EQs:
// 		カーネルの座標．補間点からカーネルの座標までの距離rを計算するためにつかう
// 		カーネルが通る値
// 	*/
// 	//! ------------------ 微分の値を考慮して重み関数を計算 ------------------ */
// 	InterpolationRBF_Common(const std::vector<std::tuple<T, double>> &kernel_X_s_IN,
// 							const std::vector<std::tuple<T, U>> &EQs_IN)
// 		: kernel_X_s_w(kernel_X_s_IN.size()), M(kernel_X_s_IN.size(), V_d(kernel_X_s_IN.size())), EQs(EQs_IN){};
// };
//@ ------------------------------------------------------ */
template <typename T, typename U> struct InterpolationRBF_Common {
  VV_d M;
  std::vector<U> b;
  std::vector<std::tuple<T, double, U>> kernel_X_s_w;
  std::vector<std::tuple<T, U>> EQs;
  InterpolationRBF_Common(const std::vector<std::tuple<T, double>> &kernel_X_s_IN, const std::vector<std::tuple<T, U>> &EQs_IN) : kernel_X_s_w(kernel_X_s_IN.size()), M(kernel_X_s_IN.size(), V_d(kernel_X_s_IN.size())), b(EQs_IN.size()), EQs(EQs_IN) { set_M_b(kernel_X_s_IN, EQs_IN); };
  InterpolationRBF_Common() {};
  void set_M_b(const std::vector<std::tuple<T, double>> &kernel_X_s_IN, const std::vector<std::tuple<T, U>> &EQs_IN) {
    try {
      // std::cout << __PRETTY_FUNCTION__ << "initialize" << std::endl;
      this->M.resize(kernel_X_s_IN.size(), V_d(kernel_X_s_IN.size()));
      this->b.resize(EQs_IN.size());
      this->kernel_X_s_w.resize(kernel_X_s_IN.size());
      this->EQs = EQs_IN;
      //* ----------------------- Mとbの計算 ----------------------- */
      int i = 0;
      for (const auto &X_s : kernel_X_s_IN) {
        std::get<0>(kernel_X_s_w[i]) = std::get<0>(X_s);
        std::get<1>(kernel_X_s_w[i++]) = std::get<1>(X_s);
      }
      //
      V_d row(kernel_X_s_IN.size());
      for (auto i = 0; i < EQs.size(); ++i) {
        auto [sample_X, res] = EQs[i];
        b[i] = res;
        for (auto j = 0; j < kernel_X_s_IN.size(); ++j) {
          auto [kernel_X, s] = kernel_X_s_IN[j];
          row[j] = kernel_MQ(sample_X, kernel_X, 1. / s);
        }
        M[i] = row;
      }
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
  };
  ;
};
//@ ------------------------------------------------------ */
template <typename T, typename U> struct InterpolationRBF_ : public InterpolationRBF_Common<T, U> {
  InterpolationRBF_(const std::vector<std::tuple<T, U>> &kernel_X_s_IN, const std::vector<std::tuple<T, U>> &EQs_IN) : InterpolationRBF_Common<T, U>(kernel_X_s_IN, EQs_IN) {};
};
//@ ------------------------------------------------------ */
template <> struct InterpolationRBF_<Tdd, Tddd> : public InterpolationRBF_Common<Tdd, Tddd> {
  InterpolationRBF_() : InterpolationRBF_Common<Tdd, Tddd>() {};
  InterpolationRBF_(const std::vector<std::tuple<Tdd, double>> &kernel_X_s_IN, const std::vector<std::tuple<Tdd, Tddd>> &EQs_IN) : InterpolationRBF_Common<Tdd, Tddd>(kernel_X_s_IN, EQs_IN) { calculate_weight(); };
  //! ------------------------------------------------------ */
  void initialize(const std::vector<std::tuple<Tdd, double>> &kernel_X_s_IN, const std::vector<std::tuple<Tdd, Tddd>> &EQs_IN) {
    InterpolationRBF_Common<Tdd, Tddd>::set_M_b(kernel_X_s_IN, EQs_IN);
    calculate_weight();
  };
  //! ------------------------------------------------------ */
  void calculate_weight() {
    lapack_lu lu(this->M);
    V_d w(kernel_X_s_w.size());
    V_d b_(kernel_X_s_w.size());
    //! wは，各kernel_X_s_wの最後に保存している
    /* -------------------- */
    for (auto i = 0; i < b.size(); ++i)
      b_[i] = std::get<0>(b[i]);
    lu.solve(b_, w);
    for (auto i = 0; i < kernel_X_s_w.size(); ++i)
      std::get<0>(std::get<2>(kernel_X_s_w[i])) = w[i];
    /* -------------------- */
    for (auto i = 0; i < b.size(); ++i)
      b_[i] = std::get<1>(b[i]);
    lu.solve(b_, w);
    for (auto i = 0; i < kernel_X_s_w.size(); ++i)
      std::get<1>(std::get<2>(kernel_X_s_w[i])) = w[i];
    /* -------------------- */
    for (auto i = 0; i < b.size(); ++i)
      b_[i] = std::get<2>(b[i]);
    lu.solve(b_, w);
    for (auto i = 0; i < kernel_X_s_w.size(); ++i)
      std::get<2>(std::get<2>(kernel_X_s_w[i])) = w[i];
  };
  //! ------------------------------------------------------ */
  Tddd /*!*/ operator()(const Tdd /*!*/ &x) const {
    Tddd /*!*/ ret = {0, 0, 0};
    for (const auto &X_s_w : kernel_X_s_w)
      ret += std::get<2>(X_s_w) * kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));
    return ret;
  };
  //! ------------------------------------------------------ */
  T3Tdd grad(const Tdd &x) const {
    T3Tdd ret{{{0, 0}, {0, 0}, {0, 0}}};
    for (const auto &X_s_w : kernel_X_s_w) {
      std::get<0>(ret) += std::get<0>(std::get<2>(X_s_w)) * grad_kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));
      std::get<1>(ret) += std::get<1>(std::get<2>(X_s_w)) * grad_kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));
      std::get<2>(ret) += std::get<2>(std::get<2>(X_s_w)) * grad_kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));
    }
    /*
     * this returns {{dx/dt0,dx/dt1},{dy/dt0,dy/dt1},{dz/dt0,dz/dt1}}
     */
    return ret;
  };
  //! ------------------------------------------------------ */
  Tdd grad_total_t0t1(const Tdd &x) const {
    return Total(this->grad(x)); // t0成分とt1成分をそれぞれ足し合わせる
  };
  //! ------------------------------------------------------ */
  Tddd grad_total_xyz(const Tdd &x) const {
    auto [dxdt0_dxdt1, dydt0_dydt1, dzdt0_dzdt1] = this->grad(x);
    return {Total(dxdt0_dxdt1), Total(dydt0_dydt1), Total(dzdt0_dzdt1)};
  };
};
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
template <> struct InterpolationRBF_<Tddd, double> : public InterpolationRBF_Common<Tddd, double> {
  InterpolationRBF_(const std::vector<std::tuple<Tddd, double>> &kernel_X_s_IN, const std::vector<std::tuple<Tddd, double>> &EQs_IN) : InterpolationRBF_Common<Tddd, double>(kernel_X_s_IN, EQs_IN) {
    lapack_lu lu(this->M);
    V_d w(kernel_X_s_w.size());
    V_d b_(kernel_X_s_w.size());
    lu.solve(b, w);
    for (auto i = 0; i < kernel_X_s_w.size(); ++i)
      std::get<2>(kernel_X_s_w[i]) = w[i];
  };
  //! ------------------------------------------------------ */
  double /*!*/ operator()(const Tddd /*!*/ &x) const {
    double /*!*/ ret(0.);
    for (const auto &X_s_w : kernel_X_s_w)
      ret += std::get<2>(X_s_w) * kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));
    return ret;
  };
  //! ------------------------------------------------------ */
  Tddd grad(const Tddd &x) const {
    Tddd ret = {0, 0, 0};
    for (const auto &X_s_w : kernel_X_s_w)
      ret += std::get<2>(X_s_w) * grad_kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));
    return ret;
  };
  //! ------------------------------------------------------ */
  double laplacian(const Tddd &x) const {
    double ret = 0.;
    for (const auto &X_s_w : kernel_X_s_w)
      ret += std::get<2>(X_s_w) * laplacian_kernel_MQ(x, std::get<0>(X_s_w), 1. / std::get<1>(X_s_w));

    return ret;
  };
};
//@ ------------------------------------------------------ */
class InterpolationIDW {
  using V_d = std::vector<double>;
  using VV_d = std::vector<std::vector<double>>;
  using VVV_d = std::vector<std::vector<std::vector<double>>>;

private:
  V_d w;                                 // weight
  VV_d A;                                // position
  V_d V;                                 // values of position
  VV_d VV;                               // values of position
  std::function<double(V_d, V_d)> phi;   // RBF basis function passed as a lambda function
  std::function<V_d(V_d, V_d)> grad_phi; // derivative of the RBF basis function passed as a lambda function
  int parameter_dim;
  double scale;
  double p; // std::pow

public:
  // ディフォルトの補間基底関数，多重二乗を使う場合
  InterpolationIDW(const VV_d &A_IN, const V_d &V_IN) : A(A_IN), V(V_IN), VV({}), parameter_dim(A_IN[0].size()) {};
  //
  InterpolationIDW(const VV_d &A_IN, const V_d &V_IN, const double p_IN /*std::pow*/) : A(A_IN), V(V_IN), VV({}), parameter_dim(A_IN[0].size()), p(p_IN) {};
  //
  InterpolationIDW(const VV_d &A_IN, const VV_d &VV_IN, const double p_IN = 1. /*std::pow*/) : A(A_IN), V({}), VV(VV_IN), parameter_dim(A_IN[0].size()), p(p_IN) {};
  //
  double operator()(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }

    double ret(0), tmp, sum(0.);
    for (auto i = 0; i < this->A.size(); i++) {
      if (Norm(this->A[i] - x) < 1E-15) {
        return this->V[i];
      };
      tmp = 1. / std::pow(Norm(this->A[i] - x), this->p);
      sum += tmp;
      ret += tmp * this->V[i];
    }
    return ret / sum;
  };

  V_d calcVV(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }

    V_d invA({});
    V_d ret(this->VV[0].size(), 0.);
    double tmp, sum(0.);
    for (auto i = 0; i < this->A.size(); i++) {
      tmp = 1. / std::pow(Norm(this->A[i] - x), this->p);
      sum += tmp;
      ret += tmp * this->VV[i];
    }
    return ret / sum;
  };

  V_d nabla(const V_d &x) const {
    if (this->parameter_dim != x.size()) {
      std::string message = "The size must be " + std::to_string(parameter_dim) + ". Given argument size is " + std::to_string(x.size());
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message);
    }

    double W = 0.;
    for (auto i = 0; i < this->A.size(); i++)
      W += std::pow(Norm(this->A[i] - x), -this->p);

    auto f = (*this)(x);

    V_d ret(x.size(), 0.);
    for (auto i = 0; i < A.size(); i++) {
      ret += std::pow(Norm(x - A[i]), -p / 2. - 1.) * (x - A[i]) * (V[i] - f);
    }

    return -ret / W * p;
  };
};

// /* -------------------------------------------------------- */
// /*                    三角面のための線形補間　３点               */
// /* -------------------------------------------------------- */
// class interpolationTriangleLinear
// {
// private:
// 	VV_d s;

// public:
// 	//   0
// 	// 1   2
// 	const V_d dndt0 = {1., 0., -1.};
// 	const V_d dndt1 = {0., 1., -1.};
// 	interpolationTriangleLinear(const VV_d &s_IN) : s(s_IN){};
// 	void set(const VV_d &s_IN) { this->s = s_IN; };
// 	V_d N(const double &t0, const double &t1) const { return {t0, t1, 1 - t0 - t1}; };
// 	V_d dNdt0(const double &t0, const double &t1) const { return dndt0; };
// 	V_d dNdt1(const double &t0, const double &t1) const { return dndt1; };
// 	V_d operator()(const double &t0, const double &t1) const
// 	{
// 		return s[0] * t0 + s[1] * t1 + s[2] * (1 - t0 - t1);
// 	};
// 	V_d cross(const double &t0, const double &t1) const
// 	{
// 		return Cross(s[0] - s[2], s[1] - s[2]);
// 	};
// 	VV_d div(const double &t0, const double &t1) const
// 	{
// 		return {s[0] - s[2], s[1] - s[2]};
// 	};
// 	double J(const double &t0, const double &t1)
// 	{
// 		VV_d a = {s[0] - s[2], s[1] - s[2]};
// 		return std::sqrt(Dot(a[0], a[0]) * Dot(a[1], a[1]) -std::pow(Dot(a[0], a[1]), 2.));
// 	};
// };
//@ -------------------------------------------------------- */
//@                    三角面のための線形補間　３点               */
//@ -------------------------------------------------------- */
using Tddd = Tddd;
struct interpolationTriangleLinearTuple {
  T3Tddd s; // sample points
  //   0  (t0,t1)=(1,0)
  // 1   2 (t0,t1)=(0,0)
  /*
   * 全体(0,1,2)をパラメタ[0,1],[0,1]で線形補間する
   */
  interpolationTriangleLinearTuple() {};
  // ~interpolationTriangleLinear() { std::cout << "interpolationTriangleLinear deleted" << std::endl; };
  const Tddd dndt0 = {1., 0., -1.};
  const Tddd dndt1 = {0., 1., -1.};
  Tddd N(const double &t0, const double &t1) const { return {t0, t1, 1 - t0 - t1}; };
  Tddd dNdt0(const double &t0, const double &t1) const { return dndt0; };
  Tddd dNdt1(const double &t0, const double &t1) const { return dndt1; };
  interpolationTriangleLinearTuple(const T3Tddd &s_IN) : s(s_IN) {};
  /* ------------------------------------------------------ */
  void set(const T3Tddd &s_IN) { this->s = s_IN; };

  Tddd operator()(const double t0, const double t1) const { return Dot(this->N(t0, t1), this->s); };
  Tddd dXdt0(const double t0, const double t1) const { return Dot(this->dNdt0(t0, t1), this->s); };
  Tddd dXdt1(const double t0, const double t1) const { return Dot(this->dNdt1(t0, t1), this->s); };
  T2Tddd div(const double t0, const double t1) { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  Tddd cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    auto dxdt0 = this->dXdt0(t0, t1);
    auto dxdt1 = this->dXdt1(t0, t1);
    return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
  };
};
////////////////////////////////////////////////////////////////////////
struct interpolationTriangleLinear3D {
  T3Tddd s; // sample points
  interpolationTriangleLinear3D(const T3Tddd &s_IN) : s(s_IN) {};
  Tddd N(const double &t0, const double &t1) const { return {t0, t1, 1 - t0 - t1}; };
  Tddd dNdt0(const double &t0, const double &t1) const { return {1., 0., -1.}; };
  Tddd dNdt1(const double &t0, const double &t1) const { return {0., 1., -1.}; };
  /* ------------------------------------------------------ */
  Tddd operator()(const double t0, const double t1) const { return Dot(this->N(t0, t1), this->s); };
  Tddd dXdt0(const double t0, const double t1) const { return Dot(this->dNdt0(t0, t1), this->s); };
  Tddd dXdt1(const double t0, const double t1) const { return Dot(this->dNdt1(t0, t1), this->s); };
  T2Tddd div(const double t0, const double t1) { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  Tddd cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    Tddd dxdt0 = this->dXdt0(t0, t1);
    Tddd dxdt1 = this->dXdt1(t0, t1);
    return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
  };
};
////////////////////////////////////////////////////////////////////////
struct interpolationTriangleLinearByFixedRange3D {
  T3Tddd s; // sample points
  interpolationTriangleLinearByFixedRange3D(const T3Tddd &s_IN) : s(s_IN) {};
  Tddd N(const double &t0, const double &t1) const { return {t0, t1 * (1. - t0), (-1. + t0) * (-1. + t1)}; };
  Tddd dNdt0(const double &t0, const double &t1) const { return {1., -t1, -1. + t1}; };
  Tddd dNdt1(const double &t0, const double &t1) const { return {0., 1. - t0, -1. + t0}; };
  /* ------------------------------------------------------ */
  Tddd operator()(const double t0, const double t1) const {
    // return Dot(this->N(t0, t1), this->s);
    return t0 * std::get<0>(this->s) + t1 * (1. - t0) * std::get<1>(this->s) + (-1. + t0) * (-1. + t1) * std::get<2>(this->s);
  };
  Tddd dXdt0(const double t0, const double t1) const {
    // return Dot(this->dNdt0(t0, t1), this->s);
    return std::get<0>(this->s) + -t1 * std::get<1>(this->s) + (-1. + t1) * std::get<2>(this->s);
  };
  Tddd dXdt1(const double t0, const double t1) const {
    // return Dot(this->dNdt1(t0, t1), this->s);
    return (1. - t0) * std::get<1>(this->s) + (-1. + t0) * std::get<2>(this->s);
  };
  T2Tddd div(const double t0, const double t1) const { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  Tddd cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    Tddd dxdt0 = this->dXdt0(t0, t1);
    Tddd dxdt1 = this->dXdt1(t0, t1);
    return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
  };
};
//! ------------------------------------------------------ */
//!                    改善した線形三角補間                   */
//! ------------------------------------------------------ */

// 共通部分：ここではNなど，を取り出したい場合，まずそれ単体でテンプレートを作って，それを別のテンプレートストラクトで継承する
template <typename T> struct interpolationTriangleLinearCommon {
  std::array<T, 3> s; // sample points
  interpolationTriangleLinearCommon(const std::array<T, 3> &s_IN) : s(s_IN) {};
  Tddd N(const double &t0, const double &t1) const { return {t0, t1, 1 - t0 - t1}; };
  Tddd dNdt0(const double &t0, const double &t1) const { return {1., 0., -1.}; };
  Tddd dNdt1(const double &t0, const double &t1) const { return {0., 1., -1.}; };
  T operator()(const double t0, const double t1) const { return Dot(this->N(t0, t1), this->s); };
  T dXdt0(const double t0, const double t1) const { return Dot(this->dNdt0(t0, t1), this->s); };
  T dXdt1(const double t0, const double t1) const { return Dot(this->dNdt1(t0, t1), this->s); };
};

template <typename T> struct interpolationTriangleLinear_ : public interpolationTriangleLinearCommon<T> {
  interpolationTriangleLinear_(const std::array<T, 3> &s_IN) : interpolationTriangleLinearCommon<T>(s_IN) {};
};

template <> struct interpolationTriangleLinear_<double> : public interpolationTriangleLinearCommon<double> {
  interpolationTriangleLinear_(const Tddd &s_IN) : interpolationTriangleLinearCommon<double>(s_IN) {};
  Tdd grad(const double t0, const double t1) const { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  double J(const double t0, const double t1) const { return this->dXdt0(t0, t1) * this->dXdt1(t0, t1); };
};

template <> struct interpolationTriangleLinear_<Tddd> : public interpolationTriangleLinearCommon<Tddd> {
  interpolationTriangleLinear_(const T3Tddd &s_IN) : interpolationTriangleLinearCommon<Tddd>(s_IN) {};
  /* ------------------------------------------------------ */
  T2Tddd div(const double t0, const double t1) const { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  Tddd cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    auto dxdt0 = this->dXdt0(t0, t1);
    auto dxdt1 = this->dXdt1(t0, t1);
    return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
  };
};
//! ------------------------------------------------------ */
//! ------------------------------------------------------ */
template <typename T> struct interpolationTriangleLinearCommon0101 {
  std::array<T, 3> s; // sample points
  interpolationTriangleLinearCommon0101(const std::array<T, 3> &s_IN) : s(s_IN) {};
  Tddd N(const double &t0, const double &t1) const { return {t0, t1 * (1. - t0), (-1. + t0) * (-1. + t1)}; };
  Tddd dNdt0(const double &t0, const double &t1) const { return {1., -t1, -1. + t1}; };
  Tddd dNdt1(const double &t0, const double &t1) const { return {0., 1. - t0, -1. + t0}; };
  T operator()(const double t0, const double t1) const { return Dot(this->N(t0, t1), this->s); };
  T dXdt0(const double t0, const double t1) const { return Dot(this->dNdt0(t0, t1), this->s); };
  T dXdt1(const double t0, const double t1) const { return Dot(this->dNdt1(t0, t1), this->s); };
};

template <typename T> struct interpolationTriangleLinear0101 : public interpolationTriangleLinearCommon0101<T> {
  interpolationTriangleLinear0101(const std::array<T, 3> &s_IN) : interpolationTriangleLinearCommon0101<T>(s_IN) {};
};

template <> struct interpolationTriangleLinear0101<double> : public interpolationTriangleLinearCommon0101<double> {
  interpolationTriangleLinear0101(const Tddd &s_IN) : interpolationTriangleLinearCommon0101<double>(s_IN) {};
  Tdd grad(const double t0, const double t1) const { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  double J(const double t0, const double t1) const { return this->dXdt0(t0, t1) * this->dXdt1(t0, t1); };
};

template <> struct interpolationTriangleLinear0101<Tddd> : public interpolationTriangleLinearCommon0101<Tddd> {
  interpolationTriangleLinear0101(const T3Tddd &s_IN) : interpolationTriangleLinearCommon0101<Tddd>(s_IN) {};
  /* ------------------------------------------------------ */
  T2Tddd div(const double t0, const double t1) const { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  Tddd cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    auto dxdt0 = this->dXdt0(t0, t1);
    auto dxdt1 = this->dXdt1(t0, t1);
    // return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
    return std::sqrt(std::fma(Dot(dxdt0, dxdt0), Dot(dxdt1, dxdt1), -std::pow(Dot(dxdt0, dxdt1), 2.)));
  };
};

//! ------------------------------------------------------ */
//! ------------------------------------------------------ */

////////////////////////////////////////////////////////////////////////
struct interpolationTriangleQuadByFixedRange3D {
  /*
  /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/interpolationLagQuadIDW20210719.nb
  のinterpQuadSimplyRangeChanged0to1がこの関数と同じ
  */
  T6Tddd s; // sample points
  bool approx_p0, approx_p1, approx_p2, no_approx;
  interpolationTriangleQuadByFixedRange3D() : approx_p0(false), approx_p1(false), approx_p2(false), no_approx(true) {};
  interpolationTriangleQuadByFixedRange3D(const T6Tddd &s_IN) : s(s_IN), approx_p0(false), approx_p1(false), approx_p2(false), no_approx(true) {};
  void noApprox() { no_approx = true; };
  void approxP0() {
    no_approx = false;
    approx_p0 = true;
  };
  void approxP1() {
    no_approx = false;
    approx_p1 = true;
  };
  void approxP2() {
    no_approx = false;
    approx_p2 = true;
  };
  T6d N(const double &t0, const double &t1) const {
    if (no_approx)
      return {t0 * (-1 + 2 * t0), (-1 + t0) * t1 * (1 + 2 * (-1 + t0) * t1), (-1 + t0) * (1 + 2 * t0 * (-1 + t1) - 2 * t1) * (-1 + t1), -4 * (-1 + t0) * t0 * t1, -4 * std::pow(-1 + t0, 2) * (-1 + t1) * t1, 4 * (-1 + t0) * t0 * (-1 + t1)};
    else if (approx_p0 && approx_p1 && approx_p2)
      return {0, 0, 0, -1 + 2 * t0 + 2 * t1 - 10 * t0 * t1 + 8 * std::pow(t0, 2) * t1, 1 - 2 * t0 + 4 * t0 * t1 - 4 * std::pow(t0, 2) * t1, 1 - 2 * t1 + 6 * t0 * t1 - 4 * std::pow(t0, 2) * t1};
    else if (approx_p0 && approx_p1)
      return {0,
              0,
              (1 - t0 - (1 - t0) * t1) * (-1 + 2 * (1 - t0 - (1 - t0) * t1)),
              -t0 + 2 * std::pow(t0, 2) - t1 + 5 * t0 * t1 - 4 * std::pow(t0, 2) * t1 + 2 * std::pow(t1, 2) - 4 * t0 * std::pow(t1, 2) + 2 * std::pow(t0, 2) * std::pow(t1, 2),
              t0 - 2 * std::pow(t0, 2) + 3 * t1 - 7 * t0 * t1 + 4 * std::pow(t0, 2) * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2),
              3 * t0 - 2 * std::pow(t0, 2) + t1 - 5 * t0 * t1 + 4 * std::pow(t0, 2) * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2)};
    else if (approx_p0 && approx_p2)
      return {0, (1 - t0) * t1 * (-1 + 2 * (1 - t0) * t1), 0, -1 + 2 * t0 + 3 * t1 - 3 * t0 * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2), 1 - 2 * t0 + t1 - t0 * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2), 1 - 3 * t1 + 3 * t0 * t1 + 2 * std::pow(t1, 2) - 4 * t0 * std::pow(t1, 2) + 2 * std::pow(t0, 2) * std::pow(t1, 2)};
    else if (approx_p1 && approx_p2)
      return {t0 * (-1 + 2 * t0), 0, 0, -1 + 3 * t0 - 2 * std::pow(t0, 2) + 2 * t1 - 2 * t0 * t1, 1 - 3 * t0 + 2 * std::pow(t0, 2), 1 + t0 - 2 * std::pow(t0, 2) - 2 * t1 + 2 * t0 * t1};
    else if (approx_p0)
      return {0, (1 - t0) * t1 * (-1 + 2 * (1 - t0) * t1), (1 - t0 - (1 - t0) * t1) * (-1 + 2 * (1 - t0 - (1 - t0) * t1)), -t0 + 2 * std::pow(t0, 2) + 4 * t0 * t1 - 4 * std::pow(t0, 2) * t1, t0 - 2 * std::pow(t0, 2) + 4 * t1 - 8 * t0 * t1 + 4 * std::pow(t0, 2) * t1 - 4 * std::pow(t1, 2) + 8 * t0 * std::pow(t1, 2) - 4 * std::pow(t0, 2) * std::pow(t1, 2), 3 * t0 - 2 * std::pow(t0, 2) - 4 * t0 * t1 + 4 * std::pow(t0, 2) * t1};
    else if (approx_p1)
      return {t0 * (-1 + 2 * t0),
              0,
              (1 - t0 - (1 - t0) * t1) * (-1 + 2 * (1 - t0 - (1 - t0) * t1)),
              -t1 + 5 * t0 * t1 - 4 * std::pow(t0, 2) * t1 + 2 * std::pow(t1, 2) - 4 * t0 * std::pow(t1, 2) + 2 * std::pow(t0, 2) * std::pow(t1, 2),
              3 * t1 - 7 * t0 * t1 + 4 * std::pow(t0, 2) * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2),
              4 * t0 - 4 * std::pow(t0, 2) + t1 - 5 * t0 * t1 + 4 * std::pow(t0, 2) * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2)};
    else
      return {t0 * (-1 + 2 * t0),
              (1 - t0) * t1 * (-1 + 2 * (1 - t0) * t1),
              0,
              -1 + 3 * t0 - 2 * std::pow(t0, 2) + 3 * t1 - 3 * t0 * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2),
              1 - 3 * t0 + 2 * std::pow(t0, 2) + t1 - t0 * t1 - 2 * std::pow(t1, 2) + 4 * t0 * std::pow(t1, 2) - 2 * std::pow(t0, 2) * std::pow(t1, 2),
              1 + t0 - 2 * std::pow(t0, 2) - 3 * t1 + 3 * t0 * t1 + 2 * std::pow(t1, 2) - 4 * t0 * std::pow(t1, 2) + 2 * std::pow(t0, 2) * std::pow(t1, 2)};
  };

  T6d dNdt0(const double &t0, const double &t1) const { return {-1 + 4 * t0, 2 * (-1 + t0) * std::pow(t1, 2) + t1 * (1 + 2 * (-1 + t0) * t1), (1 + 2 * t0 * (-1 + t1) - 2 * t1) * (-1 + t1) + 2 * (-1 + t0) * std::pow(-1 + t1, 2), -4 * (-1 + t0) * t1 - 4 * t0 * t1, -8 * (-1 + t0) * (-1 + t1) * t1, 4 * (-1 + t0) * (-1 + t1) + 4 * t0 * (-1 + t1)}; };
  T6d dNdt1(const double &t0, const double &t1) const { return {0, 2 * std::pow(-1 + t0, 2) * t1 + (-1 + t0) * (1 + 2 * (-1 + t0) * t1), (-1 + t0) * (1 + 2 * t0 * (-1 + t1) - 2 * t1) + (-1 + t0) * (-2 + 2 * t0) * (-1 + t1), -4 * (-1 + t0) * t0, -4 * std::pow(-1 + t0, 2) * (-1 + t1) - 4 * std::pow(-1 + t0, 2) * t1, 4 * (-1 + t0) * t0}; };
  /* ------------------------------------------------------ */
  Tddd operator()(const double t0, const double t1) const { return Dot(this->N(t0, t1), this->s); };
  Tddd dXdt0(const double t0, const double t1) const { return Dot(this->dNdt0(t0, t1), this->s); };
  Tddd dXdt1(const double t0, const double t1) const { return Dot(this->dNdt1(t0, t1), this->s); };
  T2Tddd div(const double t0, const double t1) const { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  Tddd cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    Tddd dxdt0 = this->dXdt0(t0, t1);
    Tddd dxdt1 = this->dXdt1(t0, t1);
    return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
  };
};
struct interpolationTriangleQuadByFixedRange3D_Center {
  T6Tddd s; // sample points
  interpolationTriangleQuadByFixedRange3D_Center(const T6Tddd &s_IN) : s(s_IN) {};
  T6d N(const double &t0, const double &t1) const { return {((-1 + t0) * t0) / 2., (t0 * (1 + t0 * (-1 + t1)) * (-1 + t1)) / 2., (t0 * t1 * (-1 + t0 * t1)) / 2., t0 * (1 + t0 * (-1 + t1)), -((1 + t0 * (-1 + t1)) * (-1 + t0 * t1)), t0 - std::pow(t0, 2) * t1}; };
  T6d dNdt0(const double &t0, const double &t1) const { return {-0.5 + t0, ((1 + 2 * t0 * (-1 + t1)) * (-1 + t1)) / 2., (t1 * (-1 + 2 * t0 * t1)) / 2., 1 + 2 * t0 * (-1 + t1), -1 - 2 * t0 * (-1 + t1) * t1, 1 - 2 * t0 * t1}; };
  T6d dNdt1(const double &t0, const double &t1) const { return {0, (t0 * (1 + 2 * t0 * (-1 + t1))) / 2., (t0 * (-1 + 2 * t0 * t1)) / 2., std::pow(t0, 2), std::pow(t0, 2) * (1 - 2 * t1), -std::pow(t0, 2)}; };
  /* ------------------------------------------------------ */
  Tddd operator()(const double t0, const double t1) const { return Dot(this->N(t0, t1), this->s); };
  Tddd dXdt0(const double t0, const double t1) const { return Dot(this->dNdt0(t0, t1), this->s); };
  Tddd dXdt1(const double t0, const double t1) const { return Dot(this->dNdt1(t0, t1), this->s); };
  T2Tddd div(const double t0, const double t1) const { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  Tddd cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    Tddd dxdt0 = this->dXdt0(t0, t1);
    Tddd dxdt1 = this->dXdt1(t0, t1);
    return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
  };
};
/////////////////////////////////////////////////////////////////////////
class interpolationTriangle {
public:
  VV_d s; // sample points

  interpolationTriangle() {};
  // ~interpolationTriangle() { std::cout << "interpolationTriangle deleted" << std::endl; };
  interpolationTriangle(const VV_d &s_IN) : s(s_IN) {};
  virtual V_d N(const double &t0, const double &t1) const = 0;
  virtual V_d dNdt0(const double &t0, const double &t1) const = 0;
  virtual V_d dNdt1(const double &t0, const double &t1) const = 0;
  /* ------------------------------------------------------ */
  void set(const VV_d &s_IN) { this->s = s_IN; };

  V_d operator()(const double t0, const double t1) const { return Dot(this->N(t0, t1), this->s); };
  V_d dXdt0(const double t0, const double t1) const { return Dot(this->dNdt0(t0, t1), this->s); };
  V_d dXdt1(const double t0, const double t1) const { return Dot(this->dNdt1(t0, t1), this->s); };
  VV_d div(const double t0, const double t1) { return {this->dXdt0(t0, t1), this->dXdt1(t0, t1)}; };
  V_d cross(const double t0, const double t1) const { return Cross(this->dXdt0(t0, t1), this->dXdt1(t0, t1)); };
  double J(const double t0, const double t1) const {
    V_d dxdt0 = this->dXdt0(t0, t1);
    V_d dxdt1 = this->dXdt1(t0, t1);
    return std::sqrt(Dot(dxdt0, dxdt0) * Dot(dxdt1, dxdt1) - std::pow(Dot(dxdt0, dxdt1), 2.));
  };
};
//@ -------------------------------------------------------- */
//@                    三角面のための線形補間　３点               */
//@ -------------------------------------------------------- */
class interpolationTriangleLinear : public interpolationTriangle {
public:
  //   0
  // 1   2
  /*
   * 全体(0,1,2)をパラメタ[0,1],[0,1]で線形補間する
   */
  interpolationTriangleLinear() {};
  // ~interpolationTriangleLinear() { std::cout << "interpolationTriangleLinear deleted" << std::endl; };
  const V_d dndt0 = {1., 0., -1.};
  const V_d dndt1 = {0., 1., -1.};
  interpolationTriangleLinear(const VV_d &s_IN) : interpolationTriangle(s_IN) {};
  V_d N(const double &t0, const double &t1) const override { return {t0, t1, 1 - t0 - t1}; };
  V_d dNdt0(const double &t0, const double &t1) const override { return dndt0; };
  V_d dNdt1(const double &t0, const double &t1) const override { return dndt1; };
};
/* -------------------------------------------------------- */
/*                    三角面のための線形補間　３点               */
/* -------------------------------------------------------- */
class interpolationTriangleLinearByFixedRange : public interpolationTriangleLinear {
public:
  //   0
  // 1   2
  /*
   * 全体(0,1,2)をパラメタ[0,1],[0,1]で線形補間する
   */
  interpolationTriangleLinearByFixedRange() {};
  interpolationTriangleLinearByFixedRange(const VV_d &s_IN) : interpolationTriangleLinear(s_IN) {};
  V_d N(const double &t0, const double &t1) const override { return {t0, t1 * (1. - t0), (-1. + t0) * (-1. + t1)}; };
  V_d dNdt0(const double &t0, const double &t1) const override { return {1., -t1, -1. + t1}; };
  V_d dNdt1(const double &t0, const double &t1) const override { return {0., 1. - t0, -1. + t0}; };
};
//@ ----------------------------------------------------------------- */
//@                            2次補間 6点                             */
//@ ----------------------------------------------------------------- */
// 範囲の変更なしの，よく使われる2次の補間
class interpolationTriangleQuad : public interpolationTriangle {
  //    0
  //  3   5
  // 1  4  2
  /*
   * 全体(0,1,2)をパラメタ[0,1],[0,1]で補間する
   */
public:
  interpolationTriangleQuad() {};
  interpolationTriangleQuad(const VV_d &s_IN) : interpolationTriangle(s_IN) {};
  V_d N(const double &t0, const double &t1) const override {
    double t2 = 1 - t0 - t1;
    return {t0 * (2 * t0 - 1), t1 * (2 * t1 - 1), t2 * (2 * t2 - 1), 4 * t0 * t1, 4 * t1 * t2, 4 * t0 * t2};
  };
  V_d dNdt0(const double &t0, const double &t1) const override { return {-1. + 4. * t0, 0., 1. - 4. * (1. - t0 - t1), 4. * t1, -4. * t1, -4. * t0 + 4. * (1. - t0 - t1)}; };
  V_d dNdt1(const double &t0, const double &t1) const override { return {0., -1. + 4. * t1, 1. - 4. * (1. - t0 - t1), 4. * t0, 4. * (1. - t0 - t1) - 4. * t1, -4. * t0}; };
};
/* ----------------------------------------------------------------- */
/*                 三角面のためのラグランジュ2次補間 6点                   */
/* ----------------------------------------------------------------- */
class interpolationTriangleQuadByFixedRangeCenter : public interpolationTriangle {
  //    0
  //  3   5
  // 1  4  2
  /*
   * 中央の三角形(3,4,5)を[0,1],[0,1]で補間するようにしたもの
   */
public:
  interpolationTriangleQuadByFixedRangeCenter() {};
  interpolationTriangleQuadByFixedRangeCenter(const VV_d &s_IN) : interpolationTriangle(s_IN) {};
  V_d N(const double &t0, const double &t1) const override {
    // CForm[FullSimplify[shape2FixedRangeCenter[t0, t1]]]
    return {((-1 + t0) * t0) / 2., (t0 * (1 + t0 * (-1 + t1)) * (-1 + t1)) / 2., (t0 * t1 * (-1 + t0 * t1)) / 2., t0 * (1 + t0 * (-1 + t1)), -((1 + t0 * (-1 + t1)) * (-1 + t0 * t1)), t0 - std::pow(t0, 2) * t1};
    // return {
    // 	((-1 + t0) * t0) / 2., ((-2 + t0) * t0 * (1 + t0 * (-1 + 2 * t1)) * (3 - 4 * t1 + t0 * (-1 + 2 * t1))) / 8., ((-2 + t0) * t0 * (-1 + t0 * (-1 + 2 * t1)) * (1 - 4 * t1 + t0 * (-1 + 2 * t1))) / 8.,
    // 	-((-2 + t0) * t0 * (1 + t0 * (-1 + 2 * t1))) / 2., -(pow(-2 + t0, 2) * (-1 + t0 * (-1 + 2 * t1)) * (1 + t0 * (-1 + 2 * t1))) / 4., ((-2 + t0) * t0 * (-1 + t0 * (-1 + 2 * t1))) / 2.};
  };
  V_d dNdt0(const double &t0, const double &t1) const override {
    // CForm[FullSimplify[D[shape2FixedRangeCenter[x, t1], x] /. x -> t0]]
    return {-0.5 + t0, ((1 + 2 * t0 * (-1 + t1)) * (-1 + t1)) / 2., (t1 * (-1 + 2 * t0 * t1)) / 2., 1 + 2 * t0 * (-1 + t1), -1 - 2 * t0 * (-1 + t1) * t1, 1 - 2 * t0 * t1};
  };
  V_d dNdt1(const double &t0, const double &t1) const override {
    // CForm[FullSimplify[D[shape2FixedRangeCenter[t0, x], x] /. x -> t1]]
    return {0, (t0 * (1 + 2 * t0 * (-1 + t1))) / 2., (t0 * (-1 + 2 * t0 * t1)) / 2., std::pow(t0, 2), std::pow(t0, 2) * (1 - 2 * t1), -pow(t0, 2)};
  };
};
/* ----------------------------------------------------------------- */
/*                 三角面のためのラグランジュ2次補間 12点                  */
/* ----------------------------------------------------------------- */
class interpolationCenterTriangleQuadIDW12 : public interpolationTriangle {
  //     7  0  6
  //   8  3   5  11
  //     1  4  2
  //      9  10
  /*
   * 中央の三角形(3,4,5)を[0,1],[0,1]で補間するようにしたもの
   */
public:
  interpolationCenterTriangleQuadIDW12() {};
  interpolationCenterTriangleQuadIDW12(const VV_d &s_IN) : interpolationTriangle(s_IN) {};
  V_d N(const double &t0, const double &t1) const override {
    V_d ret = {t0 * (-3 + 2 * (-1 + t1) * t1 + t0 * (9 + (-1 + t1) * t1) + (-9 + t1 * (7 + t1 * (5 + 12 * (-2 + t1) * t1))) * std::pow(t0, 2) + (3 - (-1 + t1) * t1 * (-3 + 17 * (-1 + t1) * t1)) * std::pow(t0, 3) - (-1 + t1) * t1 * (-1 + (-1 + t1) * t1) * (1 + 6 * (-1 + t1) * t1) * std::pow(t0, 4)),
               t0 * (1 + t0 * (-1 + t1)) * (-3 + t1 + t0 * (6 + t1 * (-9 + t1 * (-10 + t1 * (-9 + 4 * t1)))) + (-3 + t1 * (7 + 2 * t1 * (-4 + t1 * (17 + (-7 + t1) * t1)))) * std::pow(t0, 2) + 2 * (7 - 3 * t1) * std::pow(t1, 2) + t1 * std::pow(t0, 3) * (1 + t1 - 10 * std::pow(t1, 2) + std::pow(t1, 3) + std::pow(t1, 4))),
               t0 * (-1 + t0 * t1) * (-6 + t0 * (18 + t1 * (-40 + t1 * (13 + (7 - 4 * t1) * t1))) + (-18 + t1 * (47 + 2 * t1 * (-15 + t1 * (-1 + t1 * (2 + t1))))) * std::pow(t0, 2) + t1 * (11 + 4 * t1 - 6 * std::pow(t1, 2)) + (-1 + t1) * (6 + (-6 + t1) * t1) * std::pow(t0, 3) * (-1 + t1 + std::pow(t1, 2))),
               t0 * ((7 + 2 * t0 * (8 + t0)) * t1 * std::pow(-1 + t0, 2) - 6 * std::pow(-1 + t0, 3) - (-1 + t0) * (5 + 2 * t0 * (18 + t0 * (-5 + 3 * t0))) * std::pow(t1, 2) + (-6 + t0 * (-41 + t0 * (72 + (6 - 19 * t0) * t0))) * std::pow(t1, 3) + t0 * (2 + t0 * (-15 + t0 * (-59 + 48 * t0))) * std::pow(t1, 4) + 3 * (2 + (11 - 9 * t0) * t0) * std::pow(t0, 2) * std::pow(t1, 5) + 2 * (-1 + t0) * std::pow(t0, 3) * std::pow(t1, 6)),
               -4 * t0 * (9 + 4 * (-1 + t1) * t1) + (36 + (-1 + t1) * t1 * (81 + 16 * (-1 + t1) * t1)) * std::pow(t0, 2) + 2 * (-6 + (-1 + t1) * t1 * (-32 + 9 * (-1 + t1) * t1)) * std::pow(t0, 3) - (-1 + t1) * t1 * (-9 + (-1 + t1) * t1 * (59 + 4 * (-1 + t1) * t1)) * std::pow(t0, 4) + t1 * (-2 + t1 * (27 + t1 * (-48 + t1 * (19 - 2 * (-3 + t1) * t1)))) * std::pow(t0, 5) + 12 * (1 + t1 - std::pow(t1, 2)),
               t0 * (12 + t1 + t0 * (-24 + t1 * (51 + t1 * (-80 + t1 * (33 + 2 * t1)))) + (12 + t1 * (-71 + t1 * (140 - 3 * t1 * (24 + t1 * (-5 + 2 * t1))))) * std::pow(t0, 2) - t1 * (-21 + t1 * (20 + t1 * (-3 + 2 * t1) * (-20 + t1 * (12 + t1)))) * std::pow(t0, 3) + (-1 + t1) * t1 * (2 + t1 * (17 + t1 * (-40 + t1 * (17 + 2 * t1)))) * std::pow(t0, 4) + (-13 + 6 * t1) * std::pow(t1, 2)),
               (1 + t0 * (-1 + t1)) * t1 * (-7 + t0 * (9 + t0 * (-2 + 3 * (-1 + t1) * t1))) * std::pow(t0, 2) * std::pow(-1 + t1, 2),
               (-1 + t1) * (-1 + t0 * t1) * (-7 + t0 * (9 + t0 * (-2 + 3 * (-1 + t1) * t1))) * std::pow(t0, 2) * std::pow(t1, 2),
               -((-1 + t0) * t0 * (1 + t0 * (-1 + t1)) * (-1 + t0 * t1) * (5 - 2 * t1 + t0 * (2 + t1 * (-7 + 2 * t1))) * std::pow(t1, 2)),
               -(t0 * (1 + t0 * (-1 + t1)) * t1 * (5 - 2 * t1 + t0 * (2 + t1 * (-7 + 2 * t1))) * std::pow(-1 + t0, 2)),
               -(t0 * (-1 + t1) * (-1 + t0 * t1) * (3 + 2 * t1 + t0 * (-3 + t1 * (3 + 2 * t1))) * std::pow(-1 + t0, 2)),
               -((-1 + t0) * t0 * (1 + t0 * (-1 + t1)) * (-1 + t0 * t1) * (3 + 2 * t1 + t0 * (-3 + t1 * (3 + 2 * t1))) * std::pow(-1 + t1, 2))};
    return ret / (12 * (1 + t1 + t0 * (-2 + t1 - std::pow(t1, 2)) - std::pow(t1, 2) + std::pow(t0 + t0 * (-1 + t1) * t1, 2)));
  };
  V_d dNdt0(const double &t0, const double &t1) const override {
    V_d ret = {-3 - 2 * t0 * (-1 + (-1 + t1) * t1) * (9 + (-1 + t1) * t1) + t1 * (-5 + t1 * (3 - 2 * (-2 + t1) * t1)) + (-42 - (-1 + t1) * t1 * (1 + (-1 + t1) * t1 * (-55 + 38 * (-1 + t1) * t1))) * std::pow(t0, 2) + 2 * (24 + (-1 + t1) * t1 * (-1 + 2 * (-1 + t1) * t1) * (-23 + 11 * (-1 + t1) * t1)) * std::pow(t0, 3) + (-27 + (-1 + t1) * t1 * (-47 + (-1 + t1) * t1 * (102 + (-1 + t1) * t1 * (13 + 42 * (-1 + t1) * t1)))) * std::pow(t0, 4) -
                   2 * (-3 + (-1 + t1) * t1 * (-5 + (-1 + t1) * t1 * (30 + (-1 + t1) * t1 * (17 + 5 * (-1 + t1) * t1)))) * std::pow(t0, 5) - 3 * (-1 + t1) * t1 * (-1 + (-1 + t1) * t1) * (1 + 6 * (-1 + t1) * t1) * std::pow(t0, 6) * std::pow(1 + (-1 + t1) * t1, 2),
               (-1 + (-1 + t1) * t1) * (3 + t1 * (-1 + 2 * t1 * (-7 + 3 * t1))) + (-42 + t1 * (67 + t1 * (93 + t1 * (5 - t1 * (-1 + 4 * t1) * (28 + 3 * (-6 + t1) * t1))))) * std::pow(t0, 2) + 2 * (24 + t1 * (-65 + t1 * (53 + t1 * (-149 + t1 * (56 + t1 * (100 + t1 * (-91 - 4 * (-7 + t1) * t1))))))) * std::pow(t0, 3) + 2 * (3 + t1 * (-11 + t1 * (41 + t1 * (-159 + t1 * (288 + t1 * (-354 + t1 * (317 + t1 * (-197 + t1 * (83 + t1 * (-19 + 2 * t1)))))))))) * std::pow(t0, 5) +
                   2 * t0 * (9 + t1 * (-4 + t1 * (-45 + t1 + (32 + t1 * (-13 + 2 * t1)) * std::pow(t1, 2)))) + std::pow(t0, 4) * (-27 + t1 * (98 + t1 * (-215 + t1 * (578 + t1 * (-665 + t1 * (478 + t1 * (-310 + t1 * (161 - 50 * t1 + 6 * std::pow(t1, 2))))))))) + 3 * t1 * std::pow(t0, 6) * (-1 + std::pow(t1, 2) * (11 - 11 * t1 + std::pow(t1, 3))) * std::pow(1 + (-1 + t1) * t1, 2),
               6 - 5 * t1 + 2 * t0 * (-1 + (-1 + t1) * t1) * (18 + t1 * (-34 + t1 * (2 + t1 * (3 + 2 * t1)))) + (84 + t1 * (-96 + t1 * (-151 + t1 * (103 + t1 * (83 + t1 * (-68 + 3 * t1 * (-3 + 4 * t1))))))) * std::pow(t0, 2) - 2 * (48 + t1 * (-64 + t1 * (-53 + t1 * (11 + t1 * (109 + t1 * (-82 + t1 * (7 + 4 * (-1 + t1) * t1))))))) * std::pow(t0, 3) +
                   2 * (-6 + t1 * (-12 + t1 * (100 + t1 * (-190 + t1 * (214 + t1 * (-169 + t1 * (86 + t1 * (-23 + t1 * (2 + t1 * (-1 + 2 * t1)))))))))) * std::pow(t0, 5) + std::pow(t0, 4) * (54 + t1 * (-53 + t1 * (-144 + t1 * (163 + t1 * (-34 + t1 * (45 + t1 * (-79 + t1 * (23 + 4 * t1 - 6 * std::pow(t1, 2))))))))) + (-21 + t1 * (13 + 2 * (5 - 3 * t1) * t1)) * std::pow(t1, 2) +
                   3 * (-1 + t1) * t1 * (6 + (-6 + t1) * t1) * std::pow(t0, 6) * (-1 + t1 + std::pow(t1, 2)) * std::pow(1 + (-1 + t1) * t1, 2),
               6 + 13 * t1 - 2 * t0 * (-1 + (-1 + t1) * t1) * (-18 + t1 * (2 + t1 * (31 + t1 * (-41 + 2 * t1)))) + (84 + t1 * (-32 + t1 * (-312 + t1 * (265 + t1 * (214 + t1 * (-179 + 4 * (11 - 3 * t1) * t1)))))) * std::pow(t0, 2) + 2 * (-48 + t1 * (76 + t1 * (119 + t1 * (-147 + t1 * (10 + t1 * (-163 + t1 * (201 + 4 * (-19 + t1) * t1))))))) * std::pow(t0, 3) +
                   (54 + t1 * (-139 + t1 * (-8 + t1 * (-100 + t1 * (266 + t1 * (128 + t1 * (-324 + t1 * (160 + t1 * (-31 + 6 * t1))))))))) * std::pow(t0, 4) - 2 * (-1 + t1) * (-6 + t1 * (10 + t1 * (12 + t1 * (94 + t1 * (-177 + t1 * (214 + t1 * (-185 + t1 * (100 + t1 * (-35 + 2 * t1))))))))) * std::pow(t0, 5) + (6 + t1 * (-8 + t1 * (-11 + 6 * t1))) * std::pow(t1, 2) + 3 * (-1 + t1) * t1 * (-2 + t1 * (4 + (-1 + t1) * t1 * (-23 + 2 * t1))) * std::pow(t0, 6) * std::pow(1 + (-1 + t1) * t1, 2),
               4 * (-1 + (-1 + t1) * t1) * (3 + (-1 + t1) * t1) - 2 * t0 * (-1 + (-1 + t1) * t1) * (24 + (-1 + t1) * t1 * (57 + 4 * (-1 + t1) * t1)) + (-72 + t1 * (266 + t1 * (-65 + 3 * t1 * (-116 + t1 * (13 - 18 * (-3 + t1) * t1))))) * std::pow(t0, 2) + 4 * (12 + (-1 + t1) * t1 * (79 + 2 * (-1 + t1) * t1 * (-27 + (-1 + t1) * t1 * (23 + 2 * (-1 + t1) * t1)))) * std::pow(t0, 3) +
                   2 * (-6 + (-1 + t1) * t1 * (-66 + (-1 + t1) * t1 * (160 + (-1 + t1) * t1 * (19 + 20 * (-1 + t1) * t1)))) * std::pow(t0, 4) - 2 * (-1 + t1) * t1 * (-1 + (-1 + t1) * t1 * (145 + (-1 + t1) * t1 * (155 + (-1 + t1) * t1 * (63 + 4 * (-1 + t1) * t1)))) * std::pow(t0, 5) - 3 * (-1 + t1) * t1 * (-2 + (-1 + t1) * t1 * (-25 + 2 * (-1 + t1) * t1)) * std::pow(t0, 6) * std::pow(1 + (-1 + t1) * t1, 2),
               -2 * t0 * (-1 + (-1 + t1) * t1) * (-24 + t1 * (51 + t1 * (-80 + t1 * (33 + 2 * t1)))) + (72 + t1 * (-280 + t1 * (385 + t1 * (209 + t1 * (-441 + t1 * (167 + 4 * t1 * (-10 + 3 * t1))))))) * std::pow(t0, 2) + (-48 + 2 * t1 * (196 + t1 * (-361 + t1 * (153 + t1 * (-170 + t1 * (329 + t1 * (-219 + 4 * t1 * (11 + t1)))))))) * std::pow(t0, 3) + (12 + t1 * (-231 + t1 * (416 + t1 * (-132 + t1 * (232 + t1 * (-564 + t1 * (432 + t1 * (-128 + (23 - 6 * t1) * t1)))))))) * std::pow(t0, 4) -
                   2 * t1 * (-29 + t1 * (6 + t1 * (211 + t1 * (-466 + t1 * (580 + t1 * (-504 + t1 * (297 + t1 * (-108 + t1 * (17 + 2 * t1))))))))) * std::pow(t0, 5) - (-1 + (-1 + t1) * t1) * (12 + t1 - 13 * std::pow(t1, 2) + 6 * std::pow(t1, 3)) + 3 * (-1 + t1) * t1 * (2 + t1 * (17 + t1 * (-40 + t1 * (17 + 2 * t1)))) * std::pow(t0, 6) * std::pow(1 + (-1 + t1) * t1, 2),
               t0 * t1 * std::pow(-1 + t1, 2) * (14 * (-1 + (-1 + t1) * t1) + t0 * (62 + t1 * (20 + t1 * (-62 + 21 * t1))) + (92 + t1 * (-93 + 2 * t1 * (30 + t1 * (-41 - 11 * (-3 + t1) * t1)))) * std::pow(t0, 3) + 2 * (-19 + t1 * (28 + t1 * (-20 + t1 * (8 + t1 * (4 + 3 * (-2 + t1) * t1))))) * std::pow(t0, 4) + 2 * std::pow(t0, 2) * (-54 + t1 * (20 + t1 * (17 + t1 - 6 * std::pow(t1, 2)))) + 3 * std::pow(t0, 5) * (2 + t1 + 3 * (-2 + t1) * std::pow(t1, 2)) * std::pow(1 + (-1 + t1) * t1, 2)),
               t0 * (-1 + t1) * std::pow(t1, 2) * (t0 * (-41 + t1 * (-41 + t1 * (-1 + 21 * t1))) + 2 * (22 + t1 * (33 + t1 * (16 + t1 * (-23 + 6 * t1)))) * std::pow(t0, 2) - (21 + t1 * (65 + 2 * t1 * (-5 + t1 * (19 + 11 * (-2 + t1) * t1)))) * std::pow(t0, 3) + (4 + 2 * t1 * (16 + t1 * (-13 + t1 * (24 + t1 * (-19 - 3 * (-4 + t1) * t1))))) * std::pow(t0, 4) + 14 * (1 + t1 - std::pow(t1, 2)) + 3 * t1 * (-2 + 3 * (-1 + t1) * t1) * std::pow(t0, 5) * std::pow(1 + (-1 + t1) * t1, 2)),
               -(pow(t1, 2) * (5 + (-3 + t1) * t1 * (-1 + 2 * t1) - 2 * t0 * (-1 + (-1 + t1) * t1) * (-8 + t1 * (-3 + 2 * t1)) + (14 + t1 * (64 - t1 * (3 + t1 * (57 + 4 * (-7 + t1) * t1)))) * std::pow(t0, 2) + (4 + 2 * t1 * (-49 + t1 * (14 + (-2 + t1) * t1 * (-1 + 4 * (-3 + t1) * t1)))) * std::pow(t0, 3) + (-11 + t1 * (71 + t1 * (-43 + t1 * (71 + t1 * (-117 + t1 * (73 + t1 * (-19 + 2 * t1))))))) * std::pow(t0, 4) -
                               2 * (-2 + t1 * (6 + t1 * (14 + t1 * (-31 + t1 * (38 + t1 * (-37 + t1 * (24 + t1 * (-11 + 2 * t1)))))))) * std::pow(t0, 5) + 3 * (-1 + t1) * t1 * (2 + t1 * (-7 + 2 * t1)) * std::pow(t0, 6) * std::pow(1 + (-1 + t1) * t1, 2))),
               -((-1 + t0) * t1 * (-5 + t1 * (-3 + (7 - 2 * t1) * t1) + 3 * t0 * (7 + t1 * (5 + t1 * (-9 + 2 * t1))) + (-27 + 4 * (-3 + t1) * (-1 + t1) * t1 * (-2 + (-2 + t1) * t1)) * std::pow(t0, 2) + (5 - 2 * t1 * (-27 + t1 * (10 + (-2 + t1) * t1 * (-13 + 4 * t1)))) * std::pow(t0, 3) + (12 + t1 * (-81 + t1 * (131 + t1 * (-127 + t1 * (89 + t1 * (-41 + (13 - 2 * t1) * t1)))))) * std::pow(t0, 4) + 3 * (-1 + t1) * (2 + t1 * (-7 + 2 * t1)) * std::pow(t0, 5) * std::pow(1 + (-1 + t1) * t1, 2))),
               -((-1 + t0) * (-1 + t1) * (3 + 3 * t0 * (5 + 2 * t1) * (-1 + (-1 + t1) * t1) - (21 + 2 * t1 * (7 + t1 * (-2 + (-2 + t1) * t1 * (9 + 4 * t1)))) * std::pow(t0, 3) + std::pow(t0, 4) * (6 + t1 * (15 + t1 * (-27 + t1 * (9 + t1 * (-9 + t1 * (-5 + t1 - 2 * std::pow(t1, 2))))))) + std::pow(t0, 2) * (27 + 4 * (-1 + t1) * t1 * (2 + t1) * (-3 + std::pow(t1, 2))) - t1 * (-5 + t1 + 2 * std::pow(t1, 2)) + 3 * t1 * (-3 + t1 * (3 + 2 * t1)) * std::pow(t0, 5) * std::pow(1 + (-1 + t1) * t1, 2))),
               -(pow(-1 + t1, 2) * ((-5 + t0 * (5 + 9 * t0 * (1 + t0))) * t1 * std::pow(-1 + t0, 3) - 3 * (-1 + 2 * t0) * std::pow(-1 + t0, 4) - std::pow(-1 + t0, 2) * (1 + 9 * t0 * (-2 + t0 + 4 * std::pow(t0, 3))) * std::pow(t1, 2) + (-1 + t0) * (2 + t0 * (-4 + t0 * (11 + t0 * (-41 + 6 * t0 * (-3 + 11 * t0))))) * std::pow(t1, 3) + t0 * (-4 + t0 * (8 + t0 * (-34 + t0 * (33 + 8 * (8 - 9 * t0) * t0)))) * std::pow(t1, 4) +
                                    (4 + t0 * (-8 + t0 * (-1 + 6 * t0 * (-4 + 7 * t0)))) * std::pow(t0, 2) * std::pow(t1, 5) - (-8 + t0 * (5 + 6 * t0 * (1 + t0))) * std::pow(t0, 3) * std::pow(t1, 6) + (-2 + (10 - 9 * t0) * t0) * std::pow(t0, 4) * std::pow(t1, 7) + 2 * (-2 + 3 * t0) * std::pow(t0, 5) * std::pow(t1, 8)))};
    return ret / (12 * std::pow(1 + t1 + t0 * (-2 + t1 - std::pow(t1, 2)) - std::pow(t1, 2) + std::pow(t0 + t0 * (-1 + t1) * t1, 2), 2));
  };
  V_d dNdt1(const double &t0, const double &t1) const override {
    V_d ret = {-(t0 * (-1 + 2 * t1) * (1 + t0 * (-3 + t0 + 2 * t0 * (-1 + t1) * t1 * (-15 + 7 * (-1 + t1) * t1) + (6 - 4 * (-1 + t1) * t1 * (-25 + (-1 + t1) * t1)) * std::pow(t0, 2) + (-9 - (-1 + t1) * t1 * (120 + (-1 + t1) * t1 * (25 + 12 * (-1 + t1) * t1))) * std::pow(t0, 3) + (1 + (-1 + t1) * t1) * (-1 + 3 * (-1 + t1) * t1 * (-3 + 2 * (-1 + t1) * t1 * (3 + (-1 + t1) * t1))) * std::pow(t0, 5) + std::pow(t0, 4) * (5 - 6 * t1 * (10 - 11 * t1 + (5 + 2 * (-3 + t1) * t1) * std::pow(t1, 3)))))),
               t0 * (2 * (-11 + t0 * (12 + 7 * t0)) * t1 * std::pow(-1 + t0, 3) - (-4 + t0 * (5 + t0)) * std::pow(-1 + t0, 4) + (-3 + 2 * t0 * (21 + t0 * (16 + 3 * t0 * (-23 + 6 * t0)))) * std::pow(-1 + t0, 2) * std::pow(t1, 2) - 2 * (-1 + t0) * (-6 + t0 * (-5 + t0 * (13 + t0 * (130 + t0 * (-187 + 46 * t0))))) * std::pow(t1, 3) + (6 + t0 * (-11 + 2 * t0 * (-47 + 3 * t0 * (-2 + t0 * (86 + 17 * (-5 + t0) * t0))))) * std::pow(t1, 4) +
                     2 * t0 * (2 + t0 * (27 - t0 * (17 + t0 * (175 - 202 * t0 + 30 * std::pow(t0, 2))))) * std::pow(t1, 5) + (-12 + t0 * (30 + t0 * (135 + (-196 + t0) * t0))) * std::pow(t0, 2) * std::pow(t1, 6) + 4 * (-2 + t0 * (-9 + t0 * (20 + 3 * t0))) * std::pow(t0, 3) * std::pow(t1, 7) - 3 * (-2 + t0 * (9 + 2 * t0)) * std::pow(t0, 4) * std::pow(t1, 8) + 2 * (2 + t0) * std::pow(t0, 5) * std::pow(t1, 9)),
               t0 * (pow(-1 + t0, 4) * (-17 + 6 * std::pow(t0, 2)) - 4 * t1 * std::pow(-1 + t0, 3) * (1 + (-17 + 9 * t0) * std::pow(t0, 2)) + (3 + t0 * (-22 + t0 * (104 + t0 * (-182 + 57 * t0)))) * std::pow(-1 + t0, 2) * std::pow(t1, 2) - 2 * (-1 + t0) * (6 + t0 * (5 + t0 * (-1 + 2 * t0) * (15 + t0 * (-38 + 3 * t0)))) * std::pow(t1, 3) + (-6 + t0 * (-9 + t0 * (4 + t0 * (12 + t0 * (49 + (16 - 69 * t0) * t0))))) * std::pow(t1, 4) +
                     2 * t0 * (2 + t0 * (-9 + t0 * (-11 + t0 * (20 + t0 * (-50 + 57 * t0))))) * std::pow(t1, 5) + (12 + t0 * (26 + t0 * (-51 + (56 - 85 * t0) * t0))) * std::pow(t0, 2) * std::pow(t1, 6) + 4 * (-2 + t0 * (3 + t0 * (2 + 9 * t0))) * std::pow(t0, 3) * std::pow(t1, 7) - 3 * (2 + t0 * (3 + 4 * t0)) * std::pow(t0, 4) * std::pow(t1, 8) + 2 * (2 + t0) * std::pow(t0, 5) * std::pow(t1, 9)),
               t0 * (-2 * (11 + 2 * t0 * (21 + t0 * (-14 + 3 * t0))) * t1 * std::pow(-1 + t0, 3) + (1 + 2 * t0 * (2 + t0)) * std::pow(-1 + t0, 4) - std::pow(-1 + t0, 2) * (6 + t0 * (54 + t0 * (-311 + 78 * t0 + 51 * std::pow(t0, 2)))) * std::pow(t1, 2) + 2 * (-1 + t0) * (6 + t0 * (49 + t0 * (42 + t0 * (-104 + t0 * (-167 + 138 * t0))))) * std::pow(t1, 3) + (6 - t0 * (-53 + t0 * (69 + 2 * t0 * (137 + t0 * (86 + 3 * t0 * (-161 + 83 * t0)))))) * std::pow(t1, 4) +
                     2 * t0 * (-2 + t0 * (20 + t0 * (121 + t0 * (18 + t0 * (-445 + 264 * t0))))) * std::pow(t1, 5) - 2 * (6 + t0 * (45 + t0 * (3 + 10 * t0 * (-26 + 17 * t0)))) * std::pow(t0, 2) * std::pow(t1, 6) + 4 * (2 + t0 * (-6 + t0 * (-41 + 33 * t0))) * std::pow(t0, 3) * std::pow(t1, 7) + 3 * (2 + (15 - 13 * t0) * t0) * std::pow(t0, 4) * std::pow(t1, 8) + 4 * (-1 + t0) * std::pow(t0, 5) * std::pow(t1, 9)),
               -(t0 * (-1 + 2 * t1) *
                 (16 + t0 * (-77 + 4 * (-2 + t1) * (-1 + t1) * t1 * (1 + t1) + 2 * t0 * (73 + (-1 + t1) * t1 * (-22 + 9 * (-1 + t1) * t1)) - 2 * (67 + (-1 + t1) * t1 * (-115 + 2 * (-1 + t1) * t1 * (-5 + 2 * (-1 + t1) * t1))) * std::pow(t0, 2) + 2 * (28 - (-1 + t1) * t1 * (173 + 2 * (-1 + t1) * t1 * (38 + 3 * (-1 + t1) * t1))) * std::pow(t0, 3) + (-5 + 2 * (-1 + t1) * t1 * (109 + 2 * (-1 + t1) * t1 * (38 + (-1 + t1) * t1 * (3 + (-1 + t1) * t1)))) * std::pow(t0, 4) +
                             2 * (1 + (-1 + t1) * t1) * (-1 + (-1 + t1) * t1 * (-24 + (-1 + t1) * t1 * (3 + (-1 + t1) * t1))) * std::pow(t0, 5)))),
               t0 * (-2 * (-1 + t0 * (-69 + 5 * t0 * (7 + 3 * t0))) * t1 * std::pow(-1 + t0, 3) - (11 + t0 * (-17 + 2 * t0)) * std::pow(-1 + t0, 4) + (1 + 9 * t0) * (6 + t0 * (-20 + t0 * (-19 + 23 * t0))) * std::pow(-1 + t0, 2) * std::pow(t1, 2) - 2 * (-1 + t0) * (6 + t0 * (49 + t0 * (-2 + t0 * (46 + t0 * (-367 + 232 * t0))))) * std::pow(t1, 3) + (-6 + t0 * (-33 + t0 * (49 + 2 * t0 * (67 + t0 * (251 + 47 * t0 * (-13 + 6 * t0)))))) * std::pow(t1, 4) -
                     2 * t0 * (2 + t0 * (16 + t0 * (65 + t0 * (84 + t0 * (-401 + 210 * t0))))) * std::pow(t1, 5) + 2 * (6 + t0 * (17 + t0 * (3 + 2 * t0 * (-74 + 43 * t0)))) * std::pow(t0, 2) * std::pow(t1, 6) + 4 * (2 + t0 * (6 + (13 - 9 * t0) * t0)) * std::pow(t0, 3) * std::pow(t1, 7) + 3 * (-2 + (-3 + t0) * t0) * std::pow(t0, 4) * std::pow(t1, 8) + 4 * (-1 + t0) * std::pow(t0, 5) * std::pow(t1, 9)),
               (-1 + t1) * std::pow(t0, 2) *
                   (7 + 7 * t1 * (-3 + (-1 + t1) * t1) + 2 * t0 * (-15 + t1 * (52 + t1 * (-6 + t1 * (-15 + 7 * t1)))) - 2 * (-25 + t1 * (95 + t1 * (-42 + t1 * (-4 + t1 * (2 + t1))))) * std::pow(t0, 2) + (15 + t1 * (-53 + t1 * (-5 + t1 * (99 + t1 * (-119 + t1 * (77 + 3 * (-7 + t1) * t1)))))) * std::pow(t0, 4) + 2 * (-1 + t1) * (1 + (-1 + t1) * t1) * std::pow(t0, 5) * (1 + (-13 + 3 * t1 * (4 + (-2 + t1) * t1)) * std::pow(t1, 2)) -
                    2 * std::pow(t0, 3) * (20 + t1 * (-78 + t1 * (41 + 5 * t1 - 7 * std::pow(t1, 3) + 6 * std::pow(t1, 4))))),
               t1 * std::pow(t0, 2) *
                   (-7 * (2 + t1 * (-2 + (-2 + t1) * t1)) + 2 * t0 * (23 + t1 * (-23 + t1 * (-9 + t1 * (-13 + 7 * t1)))) + 2 * (-27 + t1 * (12 + t1 * (32 + t1 * (14 + (-7 + t1) * t1)))) * std::pow(t0, 2) + 2 * (13 + t1 * (20 + t1 * (-76 + t1 * (55 + t1 * (-55 + (29 - 6 * t1) * t1))))) * std::pow(t0, 3) + 2 * t1 * (1 + (-1 + t1) * t1) * std::pow(t0, 5) * (3 + t1 + (-17 + 3 * t1 * (6 + (-3 + t1) * t1)) * std::pow(t1, 2)) -
                    std::pow(t0, 4) * (4 + t1 * (38 + t1 * (-96 + t1 * (78 + t1 * (-56 + 14 * t1 + 3 * std::pow(t1, 3))))))),
               -((-1 + t0) * t0 * t1 *
                 (-10 + t1 + (-18 + t1 * (-45 + 2 * t1 * (1 + t1) * (15 + 2 * (-4 + t1) * t1))) * std::pow(t0, 2) + (-2 + t1 * (77 - 2 * t1 * (46 + t1 * (-38 + t1 * (43 + t1 * (-21 + 4 * t1)))))) * std::pow(t0, 3) + t1 * (1 + (-1 + t1) * t1) * (-6 + (-2 + t1) * t1 * (-19 + t1 * (21 + t1 * (-9 + 4 * t1)))) * std::pow(t0, 5) - 2 * (-2 + t1) * std::pow(t1, 2) + t0 * (26 + t1 + (6 + t1 * (-13 + 4 * t1)) * std::pow(t1, 2)) +
                  std::pow(t0, 4) * (4 + t1 * (-28 + t1 * (8 + t1 * (23 + t1 * (-12 + 13 * t1 - 2 * std::pow(t1, 3)))))))),
               -(t0 * std::pow(-1 + t0, 2) * (5 + t1 * (-4 + 3 * t1) + t0 * (-13 + 8 * t1) + (9 + 2 * t1 * (7 + t1 * (-17 + (5 - 2 * t1) * t1))) * std::pow(t0, 2) + (1 - 2 * t1 * (18 + t1 * (-35 + t1 * (23 + 2 * (-6 + t1) * t1)))) * std::pow(t0, 3) + (1 + (-1 + t1) * t1) * (-2 + t1 * (16 + t1 * (-21 + t1 * (-1 + 5 * t1)))) * std::pow(t0, 4))),
               t0 * std::pow(-1 + t0, 2) * (4 + t1 * (-2 + 3 * t1) - t0 * (5 + 8 * t1) + (-5 + 2 * t1 * (20 + t1 * (-14 + (3 - 2 * t1) * t1))) * std::pow(t0, 2) + (9 + 2 * t1 * (-21 + t1 * (18 + t1 * (-5 + 2 * t1 * (1 + t1))))) * std::pow(t0, 3) + (1 + (-1 + t1) * t1) * (-3 + t1 * (9 + t1 * (6 + t1 * (-19 + 5 * t1)))) * std::pow(t0, 4)),
               -((-1 + t0) * t0 * (-1 + t1) *
                 (-7 + t1 * (-3 + 2 * (-1 + t1) * t1) + t0 * (24 + t1 * (10 + t1 * (-9 + t1 * (-3 + 4 * t1)))) - (27 + t1 * (29 + 2 * t1 * (-20 + t1 * (3 + 2 * (-2 + t1) * t1)))) * std::pow(t0, 2) + (7 + t1 * (61 + 2 * t1 * (-40 + t1 * (4 + t1 * (2 + (3 - 4 * t1) * t1))))) * std::pow(t0, 3) + (6 + t1 * (-60 + t1 * (93 + t1 * (-35 + (-1 + t1) * t1 * (17 + 2 * (-6 + t1) * t1))))) * std::pow(t0, 4) +
                  (-1 + t1) * (1 + (-1 + t1) * t1) * (3 + t1 * (-15 + t1 * (6 + t1 * (11 + t1 * (-3 + 4 * t1))))) * std::pow(t0, 5)))};
    return ret / (12 * std::pow(1 + t1 + t0 * (-2 + t1 - std::pow(t1, 2)) - std::pow(t1, 2) + std::pow(t0 + t0 * (-1 + t1) * t1, 2), 2));
  };
};

#endif