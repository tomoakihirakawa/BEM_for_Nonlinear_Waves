#ifndef kernelFunctions_H
#define kernelFunctions_H

#include "basic_vectors.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

//* ------------------------------------------------------ */
//*                  Kernel functions 核関数                */
//* ------------------------------------------------------ */

//! ------------------ Multiquadric 多重二乗 ----------------- */
inline double kernel_MQ(const V_d &x, const V_d &a, const double e) { return std::sqrt(std::pow(e * Norm(x - a), 2.) + 1.); };
inline V_d grad_kernel_MQ(const V_d &x, const V_d &a, const double e) { return (x - a) * (e * e / std::sqrt(std::pow(e * Norm(x - a), 2.) + 1.)); };
inline double laplacian_kernel_MQ(const V_d &x, const V_d &a /*補間点と考える*/, const double e) {
  double exyz = std::pow(e * Norm(x - a), 2.);
  return e * e * (3. + 2. * exyz) / std::pow(1. + exyz, 1.5);
};
/* ------------------------------------------------------ */
inline double kernel_MQ(const Tdd &x, const Tdd &a, const double e) { return std::sqrt(std::pow(e * Norm(x - a), 2.) + 1.); };
inline Tdd grad_kernel_MQ(const Tdd &x, const Tdd &a, const double e) { return (x - a) * (e * e / std::sqrt(std::pow(e * Norm(x - a), 2.) + 1.)); };
inline double laplacian_kernel_MQ(const Tdd &x, const Tdd &a, const double e) {
  double exyz = std::pow(e * Norm(x - a), 2.);
  return e * e * (3. + 2. * exyz) / std::pow(1. + exyz, 1.5);
};
inline double kernel_MQ(const Tddd &x, const Tddd &a, const double e) { return std::sqrt(std::pow(e * Norm(x - a), 2.) + 1.); };
inline Tddd grad_kernel_MQ(const Tddd &x, const Tddd &a, const double e) { return (x - a) * (e * e / std::sqrt(std::pow(e * Norm(x - a), 2.) + 1.)); };
inline double laplacian_kernel_MQ(const Tddd &x, const Tddd &a, const double e) {
  double exyz = std::pow(e * Norm(x - a), 2.);
  return e * e * (3. + 2. * exyz) / std::pow(1. + exyz, 1.5);
};
/* ------------------------------------------------------ */
inline double kernel_TPS(const Tdd &x, const Tdd &a, const double e) {
  double r = Norm(x - a);
  if (r < 1E-15)
    return 0.;
  return std::pow(r, 2.) * std::log(r * e);
};
inline Tdd grad_kernel_TPS(const Tdd &x, const Tdd &a, const double e) {
  double r = Norm(x - a);
  if (r < 1E-15)
    return {0., 0.};
  return -(a - x) * (1. + 2. * std::log(e * r));
};
inline double laplacian_kernel_TPS(const Tdd &x, const Tdd &a, const double e) {
  double r = Norm(x - a);
  if (r < 1E-15)
    return 0.;
  return 5. + 6. * std::log(r * e);
};
inline double kernel_TPS(const Tddd &x, const Tddd &a, const double e) {
  double r = Norm(x - a);
  if (r < 1E-15)
    return 0.;
  return std::pow(r, 2.) * std::log(r * e);
};
inline Tddd grad_kernel_TPS(const Tddd &x, const Tddd &a, const double e) {
  double r = Norm(x - a);
  if (r < 1E-15)
    return _ZEROS3_;
  return -(a - x) * (1. + 2. * std::log(e * r));
};
inline double laplacian_kernel_TPS(const Tddd &x, const Tddd &a, const double e) {
  double r = Norm(x - a);
  if (r < 1E-15)
    return 0.;
  return std::pow(r, 2.) * std::log(r * e);
};

//! --------------------------------- ５次スプライン -------------------------------- */
// \label{SPH:w_Bspline5}
/*DOC_EXTRACT 2_0_0_kernelFunctions

## 5次スプライン関数

３次元のシミュレーションで用いられる，5次スプライン関数：

$$
W_{5}(q,h) = \frac{2187}{40\pi h^3}
\begin{cases}
(1-q)^5 - 6(\frac{2}{3}-q)^5 + 15(\frac{1}{3}-q)^5 & (0 \leq q < \frac{1}{3}) \\
(1-q)^5 - 6(\frac{2}{3}-q)^5 & (\frac{1}{3} \leq q < \frac{2}{3}) \\
(1-q)^5 & (\frac{2}{3} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
$$

その勾配は，

$$
\nabla W_{5}(q,h) = \frac{2187}{40\pi h^4}
\begin{cases}
-5(1-q)^4 + 30(\frac{2}{3}-q)^4 - 75(\frac{1}{3}-q)^4 & (0 \leq q < \frac{1}{3}) \\
-5(1-q)^4 + 30(\frac{2}{3}-q)^4 & (\frac{1}{3} \leq q < \frac{2}{3}) \\
-5(1-q)^4 & (\frac{2}{3} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
$$

ヘッシアンは，

$$
\nabla^2 W_{5}(q,h) = \frac{2187}{40\pi h^5}
\begin{cases}
20(1-q)^3 - 120(\frac{2}{3}-q)^3 + 300(\frac{1}{3}-q)^3 & (0 \leq q < \frac{1}{3}) \\
20(1-q)^3 - 120(\frac{2}{3}-q)^3 & (\frac{1}{3} \leq q < \frac{2}{3}) \\
20(1-q)^3 & (\frac{2}{3} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
$$

*/
// double w_Bspline5(const double r, const double h) {
//    constexpr double a = 2187. / (40. * M_PI);
//    constexpr double one_third = 1.0 / 3.0;
//    constexpr double two_thirds = 2.0 / 3.0;
//    const double q = r / h;
//    const double c = a / (h * h * h);
//    if (q > 1.)
//       return 0;
//    else if (q < one_third)
//       return (std::pow(1. - q, 5) - 6. * std::pow(two_thirds - q, 5) + 15. * std::pow(one_third - q, 5)) * c;
//    else if (q < two_thirds)
//       return (std::pow(1. - q, 5) - 6. * std::pow(two_thirds - q, 5)) * c;
//    else
//       return (std::pow(1. - q, 5)) * c;
// };

// Tddd grad_w_Bspline5(const Tddd &xi, const Tddd &xj, const double h) {
//    constexpr double a = 2187. / (40. * M_PI);
//    constexpr double one_third = 1.0 / 3.0;
//    constexpr double two_thirds = 2.0 / 3.0;
//    constexpr std::array<double, 3> zeros = _ZEROS3_;

//    const double r = Norm(xi - xj);
//    const double q = r / h;
//    const Tddd grad_q = (xi - xj) / (r * h);
//    const double c = a / (h * h * h);

//    if (q > 1. || r < 1E-14)
//       return zeros;
//    else if (q < one_third)
//       return -grad_q * h * (5. * std::pow(1 - q, 4) - 30. * std::pow(two_thirds - q, 4) + 75. * std::pow(one_third - q, 4)) * c;
//    else if (q < two_thirds)
//       return -grad_q * h * (5. * std::pow(1 - q, 4) - 30. * std::pow(two_thirds - q, 4)) * c;
//    else
//       return -grad_q * h * (5. * std::pow(1 - q, 4)) * c;
// };
/* -------------------------------------------------------------------------- */

inline double w_Spiky(double r, const double h) {
  if (r < h) {
    return 15. / (M_PI * std::pow(h, 6)) * std::pow(h - r, 3);
  } else
    return 0.0;
}

inline Tddd grad_w_Spiky(const Tddd &xi, const Tddd &xj, const double h) {
  Tddd rij = xi - xj;
  double r = Norm(rij);
  if (r / h < 1E-20 || r / h > 1)
    return _ZEROS3_;
  else {
    double q = h - r;
    return (45.0 / (M_PI * std::pow(h, 6)) * std::pow(q, 2)) * (-rij / r);
  }
}

/* -------------------------------------------------------------------------- */

inline double w_Bspline4_(double r, const double h) {
  constexpr double factor = 1.0 / (20.0 * M_PI);
  const double q = r / h;
  if (q > 2.5)
    return 0.;

  const double a = factor / (h * h * h);
  const double tmp = 2.5 - q;
  const double term1 = tmp * tmp * tmp * tmp;
  if (q < 0.5)
    return a * (term1 - 5. * std::pow(1.5 - q, 4) + 10. * std::pow(0.5 - q, 4));
  else if (q < 1.5)
    return a * (term1 - 5. * std::pow(1.5 - q, 4));
  else
    return a * term1;
};

inline double w_Bspline4(double r, const double h) { return w_Bspline4_(r, 0.4 * h); };

inline Tddd grad_w_Bspline4_(const Tddd &xi, const Tddd &xj, const double h) {
  constexpr double factor = 1.0 / (20.0 * M_PI);
  const Tddd xi_xj = xi - xj;
  const double r = Norm(xi_xj);
  if (r < 1E-20)
    return _ZEROS3_;

  const double q = r / h;
  if (q > 2.5)
    return _ZEROS3_;

  const Tddd a = (-4. * factor / (r * h * h * h * h)) * xi_xj;
  const double tmp = 2.5 - q;
  const double term1 = tmp * tmp * tmp;
  if (q < 0.5)
    return a * (term1 - 5. * std::pow(1.5 - q, 3) + 10. * std::pow(0.5 - q, 3));
  else if (q < 1.5)
    return a * (term1 - 5. * std ::pow(1.5 - q, 3));
  else
    return a * (term1);
};

inline Tddd grad_w_Bspline4(const Tddd &xi, const Tddd &xj, const double h) { return grad_w_Bspline4_(xi, xj, 0.4 * h); };

/* -------------------------------------------------------------------------- */

inline double w_Bspline5(double q, const double &h) {
  constexpr double a = 2187. / (40. * M_PI);
  constexpr double one_third = 1.0 / 3.0;
  constexpr double two_thirds = 2.0 / 3.0;

  if ((q /= h) > 1.)
    return 0;
  else if (q < one_third)
    return (std::pow(1 - q, 5) - 6. * std::pow(two_thirds - q, 5) + 15. * std::pow(one_third - q, 5)) * a / (h * h * h);
  else if (q < two_thirds)
    return (std::pow(1 - q, 5) - 6. * std::pow(two_thirds - q, 5)) * a / (h * h * h);
  else
    return (std::pow(1 - q, 5)) * a / (h * h * h);
};

inline Tddd grad_w_Bspline5(const Tddd &xi, const Tddd &xj, const double h) {
  constexpr double a = 2187. / (40. * M_PI);
  constexpr double one_third = 1.0 / 3.0;
  constexpr double two_thirds = 2.0 / 3.0;

  const double r = Norm(xi - xj);
  const double q = r / h;

  if (q > 1. || r < 1E-20)
    return _ZEROS3_;
  else {
    const Tddd grad_q = (xi - xj) / (r * h);
    const double c = a / (h * h * h);
    if (q < one_third) {
      if (r * h * h * h * h == 0.0)
        return _ZEROS3_;
      else
        return grad_q * -(5 * std::pow(1. - q, 4) - 30. * std::pow(two_thirds - q, 4) + 75. * std::pow(one_third - q, 4)) * c;
    } else if (q < two_thirds)
      return grad_q * -(5 * std::pow(1. - q, 4) - 30. * std::pow(two_thirds - q, 4)) * c;
    else
      return grad_q * -(5 * std::pow(1. - q, 4)) * c;
  }
};

inline double Dot_grad_w_Bspline5(const Tddd &xi, const Tddd &xj, const double h) {
  const Tddd Xij = xi - xj;
  const double r = Norm(Xij);
  const double q = r / h;
  if (q > 1. || r < 1E-20)
    return 0.;
  else
    return Dot(Xij / (r * r + 1E-4 * h * h), grad_w_Bspline5(xi, xj, h));
};

// Tddd Dot_grad_w_Bspline5(const Tddd &xi, const Tddd &xj, const double h, const std::array<Tddd, 3> &M) {
//    const Tddd Xij = xi - xj;
//    const double r = Norm(Xij);
//    const double q = r / h;
//    if (q > 1. || r < 1E-20)
//       return _ZEROS3_;
//    else
//       return Dot((Xij / (r * r)), grad_w_Bspline5(xi, xj, h));
//    // return Dot((Xij / (r * r)) * grad_w_Bspline5(xi, xj, h), M);
// };

//! --------------------------------- 3次スプライン -------------------------------- */

// \label{SPH:w_Bspline3}
/*DOC_EXTRACT 2_0_0_kernelFunctions

## 3次スプライン関数

３次元のシミュレーションで用いられる，3次スプライン関数：

$$
W_{3}(q,h) = \frac{8}{\pi h^3}
\begin{cases}
(1-6q^2+6q^3) & (0 \leq q < \frac{1}{2}) \\
2(1-q)^3 & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
$$

その勾配は，

$$
\nabla W_{3}(q,h) =
-{\nabla q}\frac{48}{\pi h^4}
\begin{cases}
-2q+3q^2 & (0 \leq q < \frac{1}{2}) \\
(1-q)^2 & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
,\quad
\nabla q = \frac{\nabla r}{h} = \frac{{\bf x}_i - {\bf x}_j}{r h}
$$

ヘッシアンは，

$$
\nabla\otimes\nabla W_{3}(q,h) =
-\frac{48}{\pi h^4}
\begin{cases}
-2\nabla \otimes{(q\nabla q)}+3\nabla \otimes{(q^2\nabla q)} & (0 \leq q < \frac{1}{2}) \\
\nabla \otimes ((1-q)^2\nabla q) & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
,\quad
\nabla q = \frac{\nabla r}{h} = \frac{{\bf x}_i - {\bf x}_j}{r h}\\
= -\frac{48}{\pi h^4}
\begin{cases}
-2\nabla \otimes{(\frac{r}{h}\frac{{\bf x}_i - {\bf x}_j}{r})}+3\nabla \otimes{(\frac{r^2}{h^2}\frac{{\bf x}_i - {\bf x}_j}{r})} & (0 \leq q < \frac{1}{2}) \\
\nabla \otimes{(\frac{r^2}{h^2}\frac{{\bf x}_i - {\bf x}_j}{r})} & (\frac{1}{2} \leq q < 1) \\
0 & (q \geq 1)
\end{cases}
$$


*/

#include <numbers>

inline double w_Linear(const double r, const double h) {
  const double q = r / h;
  if (q > 1.)
    return 0.;
  else
    return 1. - q / h;
};

inline double w_Bspline3(const double r, const double h) {
  const double q = r / h;
  if (q > 1.)
    return 0.;
  else if (q < 0.5)
    return (8. + 48. * (q - 1.) * std::pow(q, 2)) / (M_PI * h * h * h);
  else
    return 16. * std::pow((1. - q) / h, 3) / M_PI;
};

inline std::array<double, 3> grad_w_Bspline3(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
  const double r = Norm(xi - xj);
  const double q = r / h;
  if (q > 1. || q < 1E-20)
    return _ZEROS3_;
  else if (q < 0.5)
    // return (xi - xj) * (-96. + 144. * q) * q / dinom;
    return q * (-96. + 144. * q) * ((xi - xj) / r) / (std::pow(h, 4) * M_PI);
  else
    // return -48. * (xi - xj) * std::pow(1. - q, 2) / dinom;
    return (-48. * std::pow(-1. + q, 2) * ((xi - xj) / r)) / (std::pow(h, 4) * M_PI);
  //
  // const auto r = Norm(xi - xj);
  // const double q = r / h;
  // const std::array<double, 3> grad_q = (xi - xj) / (r * h);
  // const double dinom = h * h * h;
  // const double c = 48. / M_PI / dinom;
  // if (q > 1. || dinom < 1E-14)
  //    return _ZEROS3_;
  // else if (q < 0.5)
  //    return -c * (-2 * q + 3 * q * q) * grad_q;
  // else
  //    return -c * std::pow(1. - q, 2) * grad_q;
};

// double ddr_w_Bspline3(const double &r, const double &h) {
//    const double q = r / h;
//    const double dqdr = 1. / h;
//    if (q > 1. || r < 1E-20)
//       return 0.;
//    else if (q < 0.5)
//       return 8. * (-12. * q + 18. * q * q) * dqdr / (M_PI * h * h * h);
//    else
//       return 8. * 2. * 3 * std::pow(1. - q, 2) * (-dqdr) / (M_PI * h * h * h);
// };

// double ddr2_w_Bspline3(const double &r, const double &h) {
//    const double q = r / h;
//    const double dqdr = 1. / h;
//    if (q > 1. || r < 1E-20)
//       return 0.;
//    else if (q < 0.5)
//       return 8. * ((-12. + 36. * q) * dqdr) / (M_PI * h * h * h);
//    else
//       return 8. * 2. * 3 * 2 * std::pow(1. - q, 1) * (-dqdr) * (-dqdr) / (M_PI * h * h * h);
// };

// std::array<std::array<double, 3>, 3> H_grad_w_Bspline3(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
//    const auto r = Norm(xi - xj);
//    const auto grad_r = (xi - xj) / Norm(xi - xj);
//    const double q = r / h;
//    const std::array<double, 3> grad_q = (xi - xj) / (r * h);
//    std::array<std::array<double, 3>, 3> I = {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}};
//    const auto H_q = TensorProduct((xi - xj) / h, grad_r / (r * r)) + (I / h) / r;
//    const double dinom = h * h * h;
//    const double c = 8. / M_PI / dinom;
//    if (q > 1. || dinom < 1E-14)
//       return {{_ZEROS3_, _ZEROS3_, _ZEROS3_}};
//    else if (q < 0.5)
//       return c * (-6. * 2 * q + 6. * 3 * q * q) * grad_q;
//    else
//       return c * 2. * 3 * std::pow(1. - q, 2) * (-grad_q);
// };
//

inline double Dot_grad_w_Bspline3(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
  const std::array<double, 3> Xij = xi - xj;
  const double r = Norm(Xij);
  const double q = r / h;
  if (q > 1. || q < 1E-20)
    return 0.;
  else
    return Dot(Xij / (r * r + 1E-4 * h * h), grad_w_Bspline3(xi, xj, h));
};

// double Dot_grad_w_Bspline3_Dot_Modified(const Tddd &xi, const Tddd &xj, const double h, const std::array<Tddd, 3> &M) {
//    const std::array<double, 3> Xij = xi - xj;
//    const double r = Norm(Xij);
//    const double q = r / h;
//    if (q > 1. || r < 1E-20)
//       return 0.;
//    else
//       // return Dot(Xij / (r * r), grad_w_Bspline3(xi, xj, h));
//       // return Dot(Dot(Xij / (r * r), M), Dot(grad_w_Bspline3(xi, xj, h), M));
//       // return Dot(Xij / (r * r), Dot(grad_w_Bspline3(xi, xj, h), M));
//       return Dot(Xij / (r * r), grad_w_Bspline3(xi, xj, h));
// };

// Tddd grad_w_Bspline3_Dot(const Tddd &xi, const Tddd &xj, const double h, const std::array<Tddd, 3> &M) {
//    const Tddd Xij = xi - xj;
//    const double r = Norm(Xij);
//    const double q = r / h;
//    if (q > 1. || r < 1E-20)
//       return _ZEROS3_;
//    else
//       // return Dot(Xij / (r * r), grad_w_Bspline3(xi, xj, h));
//       // return Dot(Dot(Xij / (r * r), M), Dot(grad_w_Bspline3(xi, xj, h), M));
//       // return Dot(Xij / (r * r), Dot(grad_w_Bspline3(xi, xj, h), M));
//       return Dot((Xij / (r * r)) * grad_w_Bspline3(xi, xj, h), M);
// };

/* -------------------------------------------------------------------------- */

// Wendland kernel
#include <cmath>

inline double w_Wendland(const double r, const double h) {
  if (r <= h) {
    constexpr double _alpha_ = 21.0 / (16.0 * M_PI);
    const double alpha = _alpha_ / (h * h * h);
    double q = r / h;
    return 8 * alpha * std::pow(1.0 - q, 4) * (4.0 * q + 1.0);
  }
  return 0.;
}
// Gradient of Wendland kernel
// Tddd is std::array<double, 3>
inline std::array<double, 3> grad_w_Wendland_(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
  const double r = Norm(xi - xj);
  if (r <= h && r > 1E-20) {
    constexpr double _alpha_ = 21.0 / (16.0 * M_PI);
    const double alpha = _alpha_ / (h * h * h);
    const double q = r / h;
    const auto dqdx = (xi - xj) / (r * h);
    return 8. * alpha * 4. * std::pow(1.0 - q, 3) * (-dqdx) * (4.0 * q + 1.0) + 8 * alpha * std::pow(1.0 - q, 4) * (4.0 * dqdx);
  }
  return _ZEROS3_;
}

inline const auto &grad_w_Wendland = grad_w_Wendland_;

// Tddd grad_w_Wendland(const Tddd &xi, const Tddd &xj, const double h) {
//    return grad_w_Wendland_(xi, xj, h);
// }

// double Dot_grad_w_Wendland_Dot(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
//    const double r = Norm(xi - xj);
//    return Dot((xi - xj) / (r * r), grad_w_Wendland(xi, xj, h));
// };

// Tddd grad_w_Wendland_Dot(const Tddd &xi, const Tddd &xj, const double h, const std::array<Tddd, 3> &M) {
//    const double r = Norm(xi - xj);
//    return Dot(((xi - xj) / (r * r)) * grad_w_Wendland(xi, xj, h), M);
// };
/* -------------------------------------------------------------------------- */

// const auto &w_Bspline = w_Bspline5;
// const auto &grad_w_Bspline = grad_w_Bspline5;
// const auto &Dot_grad_w_Bspline_Dot = Dot_grad_w_Bspline5_Dot;
// const auto &Dot_grad_w_Bspline_Dot_Modified = Dot_grad_w_Bspline5_Dot_Modified;

// const auto &w_Bspline = w_Bspline3;
// const auto &grad_w_Bspline = grad_w_Bspline3;
// const auto &Dot_grad_w_Bspline_Dot = Dot_grad_w_Bspline3_Dot;
// const auto &ddr_w_Bspline = ddr_w_Bspline3;
// const auto &Dot_grad_w_Bspline_Dot_Modified = Dot_grad_w_Bspline3_Dot_Modified;

// make switching function depending on the smmothign length
// w_Bspline
// h < 2.6 -> w_Bspline3
// h >2.6 -> w_Bsplin5

const double between_3or5 = 2.5;
const std::array<double, 2> between_3_4_5 = {2.7, 2.7};

// #define USE_WENDLAND_KERNEL

// double w_Bspline(const double &r, const double &h) {
//    // #ifdef USE_WENDLAND_KERNEL
//    return w_Wendland(r, h);
//    // #else
//    // if (h < std::get<0>(between_3_4_5))
//    // return w_Bspline3(r, h);
//    // else
//    // if (h < std::get<1>(between_3_4_5))
//    // return w_Bspline4(r, h);
//    // else
//    // return w_Bspline5(r, h);

//    // if (h < std::get<0>(between_3_4_5))
//    //    return w_Bspline4(r, h);
//    // else
//    //    return w_Wendland(r, h);

//    // #endif
// };

// const auto &w_Bspline = w_Bspline4;
inline const auto &w_Bspline = w_Wendland;

// std::array<double, 3> grad_w_Bspline(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
//    return grad_w_Bspline4(xi, xj, h);
// };

inline std::array<double, 3> grad_w_Bspline(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
  // if (Norm(xj - xi) == 0.)
  //    return _ZEROS3_;
  // else
  // return grad_w_Bspline4(xi, xj, h);
  return grad_w_Wendland_(xi, xj, h);
  // return grad_w_Bspline4(xi, xj + 1E-4 * h * Normalize(xj - xi), h);
}

inline std::array<double, 3> grad_w_Bspline_for_dot_grad(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
  return grad_w_Wendland_(xi, xj, h);
  // return grad_w_Bspline4(xi, xj, h);
  // return grad_w_Bspline4(xi, xj + 1E-4 * h * Normalize(xj - xi), h);
}

// std::array<double, 3> grad_w_Bspline(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h, const std::array<Tddd, 3> &M) {
// #ifdef USE_WENDLAND_KERNEL
//    return Dot(grad_w_Wendland(xi, xj, h), M);
// #else
//    return Dot(M, grad_w_Bspline(xi, xj, h));
// #endif
// };

inline double Dot_grad_w_Bspline(const std::array<double, 3> &xi, const std::array<double, 3> &xj, const double h) {
  const std::array<double, 3> Xij = xi - xj;
  const double r = Norm(Xij);
  const double q = r / h;
  if (q > 1. || q < 1E-20)
    return 0.;
  else
    return Dot(Xij / (h * h * (q * q + 1E-4)), grad_w_Bspline_for_dot_grad(xi, xj, h));
};

inline double Dot_grad_w_Bspline(const std::array<double, 3> &xi, const std::array<double, 3> &xj, double h, const std::array<Tddd, 3> &M) {
  const std::array<double, 3> Xij = xi - xj;
  const double r = Norm(Xij);
  const double q = r / h;
  if (q > 1. || q < 1E-20)
    return 0.;
  else {
    // #ifdef USE_WENDLAND_KERNEL
    // return Total(Dot(M, Xij / (r * r) * grad_w_Wendland(xi, xj, h)));
    return Total(Total(M * TensorProduct(Xij / (h * h * (q * q + 1E-4)), grad_w_Bspline_for_dot_grad(xi, xj, h))));
    // auto I = M;
    // IdentityMatrix(I);
    // return Total(Total(I * TensorProduct(Xij / (r * r), grad_w_Bspline(xi, xj, h))));
    // return Dot(Xij / (r * r), grad_w_Bspline(xi, xj, h, M));
    // #endif
    // if (h < between_3or5)
    //    // return Total(Dot(Xij / (r * r) * grad_w_Bspline3(xi, xj, h), M));
    //    return ;
    // else
    //    return Dot(Xij / (r * r), Dot(grad_w_Bspline5(xi, xj, h), M));
  }

  // return Total(Dot(M /*need to be changed*/, grad_w_Bspline(xi, xj, h, M)));
};

#endif
