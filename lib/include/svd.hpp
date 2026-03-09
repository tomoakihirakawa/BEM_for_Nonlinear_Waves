#ifndef svd_H
#define svd_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

using Int = int;
using Doub = double;
using VecDoub = std::vector<double>;
using VecDoub_I = VecDoub;
using VecDoub_O = VecDoub;
using MatDoub = std::vector<VecDoub>;
using MatDoub_I = MatDoub;
using MatDoub_O = MatDoub;

template <class T>
inline T SQR(const T a) { return a * a; }

template <class T>
inline const T &MAX(const T &a, const T &b) {
   return b > a ? (b) : (a);
}

inline float MAX(const double &a, const float &b) {
   return b > a ? (b) : float(a);
}

inline float MAX(const float &a, const double &b) {
   return b > a ? float(b) : (a);
}

template <class T>
inline const T &MIN(const T &a, const T &b) {
   return b < a ? (b) : (a);
}

inline float MIN(const double &a, const float &b) {
   return b < a ? (b) : float(a);
}

inline float MIN(const float &a, const double &b) {
   return b < a ? float(b) : (a);
}

template <class T>
inline T SIGN(const T &a, const T &b) {
   return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const float &a, const double &b) {
   return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const double &a, const float &b) {
   return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

struct SVD {
   Int m, n;
   MatDoub u, v;
   VecDoub w;
   Doub eps, tsh;
   SVD(MatDoub_I &a) : m(a.size()), n(a[0].size()), u(a), v(n, std::vector<double>(n, 0.)), w(n) {
      try {
         eps = std::numeric_limits<Doub>::epsilon();
         decompose();
         reorder();
         tsh = 0.5 * std::sqrt(m + n + 1.) * w[0] * eps;
      } catch (std::exception &e) {
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }
   }

   VecDoub solve(const VecDoub_I &b, Doub thresh);
   void solve(const VecDoub_I &b, VecDoub_O &x, Doub thresh);
   void solve(const MatDoub_I &b, MatDoub_O &x, Doub thresh);

   Int rank(Doub thresh);
   Int nullity(Doub thresh);
   MatDoub range(Doub thresh);
   MatDoub nullspace(Doub thresh);

   Doub inv_condition() {
      return (w[0] <= 0. || w[n - 1] <= 0.) ? 0. : w[n - 1] / w[0];
   }

   void decompose();
   void reorder();
   Doub pythag(const Doub a, const Doub b);
};
inline void SVD::solve(const VecDoub_I &b, VecDoub_O &x, Doub thresh = -1.) {
   Int i, j, jj;
   Doub s;
   if (b.size() != m || x.size() != n)
      throw("SVD::solve bad sizes");
   VecDoub tmp(n);
   tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.) * w[0] * eps);
   // #ifdef _OPENMP
   // #pragma omp parallel for
   // #endif
   for (j = 0; j < n; j++) {
      s = 0.0;
      if (w[j] > tsh) {
         for (i = 0; i < m; i++)
            s += u[i][j] * b[i];
         s /= w[j];
      }
      tmp[j] = s;
   }
   // #ifdef _OPENMP
   // #pragma omp parallel for
   // #endif
   for (j = 0; j < n; j++) {
      s = 0.0;
      for (jj = 0; jj < n; jj++)
         s += v[j][jj] * tmp[jj];
      x[j] = s;
   }
}

inline VecDoub SVD::solve(const VecDoub_I &b, Doub thresh = -1.) {
   Int i, j, jj;
   Doub s;
   VecDoub_O x(n);
   if (b.size() != m || x.size() != n)
      throw("SVD::solve bad sizes");
   VecDoub tmp(n);
   tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.) * w[0] * eps);
   // #ifdef _OPENMP
   // #pragma omp parallel for
   // #endif
   for (j = 0; j < n; j++) {
      s = 0.0;
      if (w[j] > tsh) {
         for (i = 0; i < m; i++)
            s += u[i][j] * b[i];
         s /= w[j];
      }
      tmp[j] = s;
   }
   // #ifdef _OPENMP
   // #pragma omp parallel for
   // #endif
   for (j = 0; j < n; j++) {
      s = 0.0;
      for (jj = 0; jj < n; jj++)
         s += v[j][jj] * tmp[jj];
      x[j] = s;
   }
   return x;
}

void SVD::solve(const MatDoub_I &b, MatDoub_O &x, Doub thresh = -1.) {
   int i, j, m = b[0].size();
   if (b.size() != n || x.size() != n || b[0].size() != x[0].size())
      throw("SVD::solve bad sizes");
   VecDoub xx(n);
   for (j = 0; j < m; j++) {
      for (i = 0; i < n; i++)
         xx[i] = b[i][j];
      solve(xx, xx, thresh);
      for (i = 0; i < n; i++)
         x[i][j] = xx[i];
   }
}
Int SVD::rank(Doub thresh = -1.) {
   Int j, nr = 0;
   tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.) * w[0] * eps);
   for (j = 0; j < n; j++)
      if (w[j] > tsh)
         nr++;
   return nr;
}

Int SVD::nullity(Doub thresh = -1.) {
   Int j, nn = 0;
   tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.) * w[0] * eps);
   for (j = 0; j < n; j++)
      if (w[j] <= tsh)
         nn++;
   return nn;
}

MatDoub SVD::range(Doub thresh = -1.) {
   Int i, j, nr = 0;
   MatDoub rnge(m, std::vector<double>(rank(thresh), 0.));
   for (j = 0; j < n; j++) {
      if (w[j] > tsh) {
         for (i = 0; i < m; i++)
            rnge[i][nr] = u[i][j];
         nr++;
      }
   }
   return rnge;
}

MatDoub SVD::nullspace(Doub thresh = -1.) {
   Int j, jj, nn = 0;
   MatDoub nullsp(n, std::vector<double>(nullity(thresh), 0.));
   for (j = 0; j < n; j++) {
      if (w[j] <= tsh) {
         for (jj = 0; jj < n; jj++)
            nullsp[jj][nn] = v[jj][j];
         nn++;
      }
   }
   return nullsp;
}
void SVD::decompose() {
   try {
      bool flag;
      Int i, its, j, jj, k, l, nm;
      Doub anorm, c, f, g, h, s, scale, x, y, z;
      VecDoub rv1(n);
      g = scale = anorm = 0.0;
      for (i = 0; i < n; i++) {
         auto &U = u[i];
         rv1[i] = scale * g;
         g = s = scale = 0.0;
         if (i < m) {
            for (k = i; k < m; k++)
               scale += std::abs(u[k][i]);
            if (scale != 0.0) {
               for (k = i; k < m; k++) {
                  u[k][i] /= scale;
                  s += u[k][i] * u[k][i];
               }
               f = U[i];
               g = -SIGN(std::sqrt(s), f);
               h = f * g - s;
               U[i] = f - g;
               for (j = i + 1; j < n; j++) {
                  for (s = 0.0, k = i; k < m; k++)
                     s += u[k][i] * u[k][j];
                  f = s / h;
                  for (k = i; k < m; k++)
                     u[k][j] += f * u[k][i];
               }
               for (k = i; k < m; k++)
                  u[k][i] *= scale;
            }
         }
         w[i] = scale * g;
         g = s = scale = 0.0;
         if (i + 1 <= m && i + 1 != n) {
            for (k = i + 1; k < n; k++)
               scale += std::abs(U[k]);
            if (scale != 0.0) {
               for (k = i + 1; k < n; k++) {
                  U[k] /= scale;
                  s += U[k] * U[k];
               }
               f = U[i + 1];
               g = -SIGN(std::sqrt(s), f);
               h = f * g - s;
               U[i + 1] = f - g;
               for (k = i + 1; k < n; k++)
                  rv1[k] = U[k] / h;
               for (j = i + 1; j < m; j++) {
                  for (s = 0.0, k = i + 1; k < n; k++)
                     s += u[j][k] * U[k];
                  for (k = i + 1; k < n; k++)
                     u[j][k] += s * rv1[k];
               }
               for (k = i + 1; k < n; k++)
                  U[k] *= scale;
            }
         }
         anorm = MAX(anorm, (std::abs(w[i]) + std::abs(rv1[i])));
      }
      for (i = n - 1; i >= 0; i--) {
         if (i < n - 1) {
            if (g != 0.0) {
               for (j = l; j < n; j++)
                  v[j][i] = (u[i][j] / u[i][l]) / g;
               for (j = l; j < n; j++) {
                  for (s = 0.0, k = l; k < n; k++)
                     s += u[i][k] * v[k][j];
                  for (k = l; k < n; k++)
                     v[k][j] += s * v[k][i];
               }
            }
            for (j = l; j < n; j++)
               v[i][j] = v[j][i] = 0.0;
         }
         v[i][i] = 1.0;
         g = rv1[i];
         l = i;
      }
      for (i = MIN(m, n) - 1; i >= 0; i--) {
         l = i + 1;
         g = w[i];
         for (j = l; j < n; j++)
            u[i][j] = 0.0;
         if (g != 0.0) {
            g = 1.0 / g;
            for (j = l; j < n; j++) {
               for (s = 0.0, k = l; k < m; k++)
                  s += u[k][i] * u[k][j];
               f = (s / u[i][i]) * g;
               for (k = i; k < m; k++)
                  u[k][j] += f * u[k][i];
            }
            for (j = i; j < m; j++)
               u[j][i] *= g;
         } else
            for (j = i; j < m; j++)
               u[j][i] = 0.0;
         ++u[i][i];
      }
      for (k = n - 1; k >= 0; k--) {
         for (its = 0; its < 30; its++) {
            flag = true;
            for (l = k; l >= 0; l--) {
               nm = l - 1;
               if (l == 0 || std::abs(rv1[l]) <= eps * anorm) {
                  flag = false;
                  break;
               }
               if (std::abs(w[nm]) <= eps * anorm)
                  break;
            }
            if (flag) {
               c = 0.0;
               s = 1.0;
               for (i = l; i < k + 1; i++) {
                  f = s * rv1[i];
                  rv1[i] = c * rv1[i];
                  if (std::abs(f) <= eps * anorm)
                     break;
                  g = w[i];
                  h = pythag(f, g);
                  w[i] = h;
                  h = 1.0 / h;
                  c = g * h;
                  s = -f * h;
                  for (j = 0; j < m; j++) {
                     y = u[j][nm];
                     z = u[j][i];
                     u[j][nm] = y * c + z * s;
                     u[j][i] = z * c - y * s;
                  }
               }
            }
            z = w[k];
            if (l == k) {
               if (z < 0.0) {
                  w[k] = -z;
                  for (j = 0; j < n; j++)
                     v[j][k] = -v[j][k];
               }
               break;
            }
            if (its == 29)
               throw("no convergence in 30 svdcmp iterations");
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;
            for (j = l; j <= nm; j++) {
               i = j + 1;
               g = rv1[i];
               y = w[i];
               h = s * g;
               g = c * g;
               z = pythag(f, h);
               rv1[j] = z;
               c = f / z;
               s = h / z;
               f = x * c + g * s;
               g = g * c - x * s;
               h = y * s;
               y *= c;
               for (jj = 0; jj < n; jj++) {
                  x = v[jj][j];
                  z = v[jj][i];
                  v[jj][j] = x * c + z * s;
                  v[jj][i] = z * c - x * s;
               }
               z = pythag(f, h);
               w[j] = z;
               if (z) {
                  z = 1.0 / z;
                  c = f * z;
                  s = h * z;
               }
               f = c * g + s * y;
               x = c * y - s * g;
               for (jj = 0; jj < m; jj++) {
                  y = u[jj][j];
                  z = u[jj][i];
                  u[jj][j] = y * c + z * s;
                  u[jj][i] = z * c - y * s;
               }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
         }
      }
   } catch (std::exception &e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   }
}

void SVD::reorder() {
   Int i, j, k, s, inc = 1;
   Doub sw;
   VecDoub su(m), sv(n);
   do {
      inc *= 3;
      inc++;
   } while (inc <= n);
   do {
      inc /= 3;
      for (i = inc; i < n; i++) {
         sw = w[i];
         for (k = 0; k < m; k++)
            su[k] = u[k][i];
         for (k = 0; k < n; k++)
            sv[k] = v[k][i];
         j = i;
         while (w[j - inc] < sw) {
            w[j] = w[j - inc];
            for (k = 0; k < m; k++)
               u[k][j] = u[k][j - inc];
            for (k = 0; k < n; k++)
               v[k][j] = v[k][j - inc];
            j -= inc;
            if (j < inc)
               break;
         }
         w[j] = sw;
         for (k = 0; k < m; k++)
            u[k][j] = su[k];
         for (k = 0; k < n; k++)
            v[k][j] = sv[k];
      }
   } while (inc > 1);
   for (k = 0; k < n; k++) {
      s = 0;
      for (i = 0; i < m; i++)
         if (u[i][k] < 0.)
            s++;
      for (j = 0; j < n; j++)
         if (v[j][k] < 0.)
            s++;
      if (s > (m + n) / 2) {
         for (i = 0; i < m; i++)
            u[i][k] = -u[i][k];
         for (j = 0; j < n; j++)
            v[j][k] = -v[j][k];
      }
   }
}

Doub SVD::pythag(const Doub a, const Doub b) {
   Doub absa = std::abs(a), absb = std::abs(b);
   return (absa > absb ? absa * std::sqrt(1.0 + SQR(absb / absa)) : (absb == 0.0 ? 0.0 : absb * std::sqrt(1.0 + SQR(absa / absb))));
}

#endif