#ifndef InterpolationRBF_H
#define InterpolationRBF_H
/*RBF_interp_detail

### 放射関数補間

距離\[r=\left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\|\]を引数とする
放射基底関数$\phi(r_i)$に重み$w_i$を掛け合わせて構築した
補間関数\[f\left( \mathbf{x} \right)=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\phi \left( \left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\| \right)}\]
を放射関数補間という．

#### 重み$w_i$の見積もり

重み$w_i$の決定には，サンプル点\[A=\left\{ {{\mathbf{a}}_{0}},{{\mathbf{a}}_{1}},...,{{\mathbf{a}}_{N-1}} \right\}\]
における値\[Y=\left\{ {{y}_{0}},{{y}_{1}},...,{{y}_{N-1}} \right\}\]
を使い，補間関数$f$も各サンプル点$A$において値$Y$となる方程式を$w_i$について解く：

$\[\left( \begin{matrix}
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
\end{matrix} \right)\]$

#### 放射基底関数$\phi$

##### 多重二乗（multiquadric）

放射基底関数として多重二乗（multiquadric），
\[\phi \left( r \right)={{\left( {{\left( \varepsilon r \right)}^{2}}+1 \right)}^{\frac{1}{2}}}\]
がよく使われる．

#### 補間関数の微分

放射関数補間の微分を少し変形すると，

$\[\nabla f\left( \mathbf{x} \right)=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\nabla \phi \left( \left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\| \right)}=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\nabla {{r}_{i}}\frac{\partial \phi \left( {{r}_{i}} \right)}{\partial {{r}_{i}}}}\]$

さらに，計算すると，

\[\begin{align}
  & {{r}_{i}}=\left\| \mathbf{x}-{{\mathbf{a}}_{i}} \right\|={{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{1/2}} \\
 & \frac{\partial {{r}_{i}}}{\partial {{\mathbf{x}}_{k}}}=\frac{1}{2}{{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{-\frac{1}{2}}}\left( \frac{\partial }{\partial {{\mathbf{x}}_{k}}}\sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right) \\
 & =\frac{1}{2}{{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{-\frac{1}{2}}}\left( \sum\limits_{j=0}^{M=2}{2\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}\cdot {{\mathbf{e}}_{k}} \right) \\
 & ={{\left( \sum\limits_{j=0}^{M=2}{{{\left( \mathbf{x}-{{\mathbf{a}}_{ij}} \right)}^{2}}} \right)}^{-\frac{1}{2}}}\overbrace{\left( {{\mathbf{x}}_{k}}-{{\mathbf{a}}_{ik}} \right)}^{\text{scaler}}=\frac{\overbrace{\left( {{\mathbf{x}}_{k}}-{{\mathbf{a}}_{ik}} \right)}^{\text{scaler}}}{{{r}_{i}}}
\end{align}\]

なので，\[\nabla {{r}_{i}}=\overbrace{\left( \mathbf{x}-{{\mathbf{a}}_{i}} \right)}^{\text{vecotr}}/{{r}_{i}}\]であり，

$\[\nabla f\left( \mathbf{x} \right)=\sum\limits_{i=0}^{N-1}{{{w}_{i}}\frac{\mathbf{x}-{{\mathbf{a}}_{i}}}{{{r}_{i}}}\frac{\partial \phi \left( {{r}_{i}} \right)}{\partial {{r}_{i}}}}\]$

である．そのため，分母がゼロになる可能性がある．
ただ，放射基底関数$\phi$が多重二乗 (Multiquadric)であれば，

$\[\phi \left( r \right)={{\left( {{\left( \varepsilon r \right)}^{2}}+1 \right)}^{\frac{1}{2}}},\frac{\partial \phi }{\partial r}\left( r \right)=\frac{\varepsilon^2 r}{\phi \left( r \right)}\]$

なので，次のように分母を消すことができる．

$\[\nabla f\left( \mathbf{x} \right)=\varepsilon^2 \sum\limits_{i=0}^{N-1}{{{w}_{i}}\frac{\mathbf{x}-{{\mathbf{a}}_{i}}}{\phi \left( {{r}_{i}} \right)}}\]$

RBF_interp_detail*/
/*RBF_interp_code*/

#include <cmath>
#include <functional>
#include <vector>
#include "basic.hpp"

class InterpolationRBF {
   using V_d = std::vector<double>;
   using VV_d = std::vector<std::vector<double>>;
   using VVV_d = std::vector<std::vector<std::vector<double>>>;

  private:
   V_d w;                                // weight
   VV_d A;                               // position
   V_d V;                                // values of position
   std::function<double(V_d, V_d)> phi;  // RBF basis function passed as a lambda function
   std::function<V_d(V_d, V_d)> dphid;   // derivative of the RBF basis function passed as a lambda function
   int argument_size;
   double scale;

  public:
   InterpolationRBF(const VV_d &A_IN, const V_d &V_IN,
                    const std::function<double(V_d, V_d)> &phi_IN,
                    const std::function<V_d(V_d, V_d)> &dphid_IN)
       : phi(phi_IN), dphid(dphid_IN), A(A_IN), V(V_IN), argument_size(A_IN[0].size()) {
      this->w = weight(A_IN, V_IN);
      this->scale = RBFscale(A_IN);
   };
   InterpolationRBF(const VV_d &A_IN, const V_d &V_IN)
       : A(A_IN), V(V_IN), argument_size(A_IN[0].size()) {
      this->scale = RBFscale(A_IN);
      this->phi = [this](const V_d &x, const V_d &a) {auto r=Norm(x-a); auto e=1./this->scale; return std::sqrt((e*r)*(e*r) + 1.); };
      this->dphid = [this](const V_d &x, const V_d &a) {auto r=Norm(x-a); auto e=1./this->scale; return e*e*(x-a)/std::sqrt((e*r)*(e*r) + 1.); };
      this->w = weight(A_IN, V_IN);
   };
   double operator()(const V_d &x) const {
      if (this->argument_size != x.size()) {
         std::string message = "The size must be " + std::to_string(argument_size) + ". Given argument size is " + std::to_string(x.size());
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message));
      }

      return std::transform_reduce(
          this->w.cbegin(), this->w.cend(), this->A.cbegin(), 0.,
          [](const auto &acc, const auto &res) { return acc + res; },
          [&x, this](const auto &w_, const auto &a) { return w_ * phi(x, a); });
   };

   V_d nabla(const V_d &x) const {
      if (this->argument_size != x.size()) {
         std::string message = "The size must be " + std::to_string(argument_size) + ". Given argument size is " + std::to_string(x.size());
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, message));
      }

      return std::transform_reduce(
          this->w.cbegin(), this->w.cend(), this->A.cbegin(), V_d(x.size(), 0.),
          [](const auto &acc, const auto &res) { return acc + res; },
          [&x, this](const auto &w_, const auto &a) { return w_ * dphid(x, a); });
   };

   //////////

   double RBFscale(const VV_d &sample) {
      V_d r;
      for (auto i = 0; i < sample.size(); i++)
         for (auto j = i + 1; j < sample.size(); j++)
            r.push_back(Norm(sample[i] - sample[j]));

      std::sort(r.begin(), r.end(), [](const auto &lhs, const auto &rhs) { return lhs < rhs; });

      auto s = 1 + (int)((double)r.size() / 3.);
      V_d v(s);
      for (auto i = 0; i < s; i++)
         v[i] = r[i];

      return Mean(v);
   };

   //////////

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

   V_d weight(VV_d A, V_d vOFx) const {
      V_d w(A.size());
      ludcmp lu(normM(A));
      lu.solve(vOFx, w);
      return w;
   };
};

/*RBF_interp_code*/
#endif
