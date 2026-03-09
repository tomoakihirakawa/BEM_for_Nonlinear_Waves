// scw.hpp
#ifndef SCW_H
#define SCW_H

#include <cmath>
#include <string>
#include <vector>
#include "basic.hpp"

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_i = std::vector<int>;
using VV_i = std::vector<std::vector<int>>;

using V_s = std::vector<std::string>;
using VV_s = std::vector<std::vector<std::string>>;

class SCW {
   double g_accel = 9.81;

  public:
   std::string filename_mode, filename_data;
   VV_d mode;
   V_d data;
   int max, l;
   double h, D, eps, theta, B, G, Qmax, LL, w, kappa, q, p;  // scale
   VV_d a;

   // ファイル名を入力
   SCW(const std::string &data_filename, const std::string &mode_filename) : mode({}), data({}) {
      auto data_IN = Load(data_filename, {","});
      auto mode_IN = Load(mode_filename, {","});
      for (const auto &d : data_IN)
         for (const auto &b : d)
            this->data.push_back(atof(b.c_str()));
      for (const auto &m : mode_IN)
         this->mode.push_back(V_d{atof(m[0].c_str()), atof(m[1].c_str()), atof(m[2].c_str())});
      set();
   };

   // 直接データを入力
   SCW(const V_d &data_IN, const VV_d &mode_IN) : data(data_IN), mode(mode_IN) {
      set();
   };

   //
   void set() {
      this->LL = 1.;  // 2.*M_PI;
      this->max = (int)data[0];
      this->l = (int)data[1];
      this->D = data[4];
      this->eps = data[5];
      this->theta = data[6];
      this->q = std::cos(data[6] / 180. * M_PI);
      this->p = std::sin(data[6] / 180. * M_PI);
      this->B = data[19];
      this->G = data[20];
      this->Qmax = data[21];
      this->kappa = 2. * M_PI / LL;
      this->h = D / kappa;  // 有次元の水深
      this->w = std::sqrt(g_accel * kappa / G);
      this->a.resize((int)data[0] + 1, V_d((int)data[0] + 1, 0.));
      for (const auto &kja : mode)
         a[(int)kja[0]][(int)kja[1]] = kja[2];
   };
   //====================================================
   double eta(const double t, const double y) {  // only for standing wave case
      // LL = 2*M_PI;// arbitraly
      // double kappa = 2.*M_PI/LL;
      // nondimensionalize
      double T = -w * t;
      double Y = kappa * q * y;
      return ETA(T, Y) / kappa;
   };
   V_d eta(const V_d &t, const V_d &y) {
      V_d ret(t.size(), 0.);
      for (auto i = 0; i < t.size(); i++)
         ret[i] = eta(t[i], y[i]);
      return ret;
   };
   //--------------------------------------------------
   double phi_without_mean_velocity(const double t, const double y, const double z) {
      // double kappa = 2*M_PI/LL;
      // nondimensionalize
      double T = -w * t;
      double Y = kappa * q * y;
      double Z = kappa * z;
      double dim_phi = w / pow(kappa, 2);  // phiの次元
      double dim_t = -1. / w;              // phiの次元
      return PHI_without_mean_velocity(T, Y, Z) * dim_phi;
   };
   //--------
   double phi(const double t, const double y, const double z) {
      // double kappa = 2*M_PI/LL;
      // nondimensionalize
      double T = -w * t;
      double Y = kappa * q * y;
      double Z = kappa * z;
      double dim_phi = w / pow(kappa, 2);  // phiの次元
      double dim_t = -1. / w;              // phiの次元
      return PHI(T, Y, Z) * dim_phi;
   };

   V_d phi(const V_d &t, const V_d &y, const V_d &z) {
      V_d ret(t.size());
      for (auto i = 0; i < t.size(); i++)
         ret[i] = phi(t[i], y[i], z[i]);
      return ret;
   };
   //--------------------------------------------------
   double phi_tn(const double t, const double y, const double z, const int n) {
      double T = -w * t;
      double Y = kappa * q * y;
      double Z = kappa * z;
      double dim_phi = w / pow(kappa, 2);  // phiの次元
      double dim_t = -1. / w;              // phiの次元
      return PHI_TN(T, Y, Z, n) * dim_phi / pow(dim_t, n);
   };
   V_d phi_tn(const V_d &t, const V_d &y, const V_d &z, const int n) {
      V_d ret(t.size());
      for (auto i = 0; i < t.size(); i++)
         ret[i] = phi_tn(t[i], y[i], z[i], n);
      return ret;
   };
   VV_d phi_tn_mat(const V_d &t, const V_d &y, const V_d &z, const int N) {
      VV_d ret(N + 1, V_d(t.size(), 0.));
      for (auto n = 0; n < N + 1; n++)
         for (auto l = 0; l < t.size(); l++)
            ret[n][l] = phi_tn(t[l], y[l], z[l], n);
      return ret;
   };
   //--------------------------------------------------
   double ETA(const double T, const double Y) {
      double ret(0.), akj, alp, hyp, hypz, hypzz, tmp1, tmp2, tmp3, df, f, delf(1.);
      double Gt, Gy, Gz, Gzt, Gyz, Gzz;
      int i = 0;
      do {
         Gt = 0.;
         Gy = 0.;
         Gz = 0.;
         Gzt = 0.;
         Gyz = 0.;
         Gzz = 0.;
         for (int k = 0; k < max; k++)
            for (int j = 2 - k % 2; j < max; j++) {
               akj = a[j][k];
               alp = std::sqrt(j * j * p * p + k * k * q * q);
               hyp = cosh(alp * ret) + sinh(alp * ret) * tanh(alp * D);
               hypz = alp * (sinh(alp * ret) + cosh(alp * ret) * tanh(alp * D));
               hypzz = alp * alp * (cosh(alp * ret) + sinh(alp * ret) * tanh(alp * D));
               tmp1 = akj * std::cos(k * Y) * std::cos(j * T);
               tmp2 = akj * std::sin(k * Y) * std::sin(j * T);
               tmp3 = akj * std::cos(k * Y) * std::sin(j * T);
               Gt += j * hyp * tmp1;
               Gy += -k * hyp * tmp2;
               Gz += hypz * akj * tmp3;
               Gzt += j * hypz * tmp1;
               Gyz += -k * hypz * tmp2;
               Gzz += hypzz * tmp3;
            }
         f = -Gt + B + (pow(p * Gt, 2.) + pow(q * Gy, 2.) + Gz * Gz) / 2. + G * ret;
         df = -Gzt + p * p * Gt * Gzt + q * q * Gy * Gyz + Gz * Gzz + G;
         delf = -f / df;
         ret += delf;
      } while (std::abs(delf) > 1E-14);
      return ret;
   };

   V_d ETA(const V_d &T, const V_d &Y) {
      V_d ret(T.size(), 0.);
      for (auto i = 0; i < T.size(); i++)
         ret[i] = ETA(T[i], Y[i]);
      return ret;
   };
   //--------------------------------------------------
   double PHI_without_mean_velocity(const double T, const double Y, const double Z) {
      double ret(0.), akj, alp, hyp;
      for (auto k = 0; k < max; k++)
         for (auto j = 2 - k % 2; j < max; j++) {
            akj = a[j][k];
            alp = std::sqrt(j * j * p * p + k * k * q * q);
            hyp = std::cosh(alp * Z) + std::sinh(alp * Z) * std::tanh(alp * D);
            ret += akj * hyp * std::sin(j * T) * std::cos(k * Y);
         }
      return ret;
   };
   //--------------------------------------------------
   double PHI(const double T, const double Y, const double Z) {
      double ret(0.), akj, alp, hyp;
      for (auto k = 0; k < max; k++)
         for (auto j = 2 - k % 2; j < max; j++) {
            akj = a[j][k];
            alp = std::sqrt(j * j * p * p + k * k * q * q);
            hyp = std::cosh(alp * Z) + std::sinh(alp * Z) * std::tanh(alp * D);
            ret += akj * hyp * std::sin(j * T) * std::cos(k * Y);
         }
      return ret - 1 / pow(p, 2) * T /*mean_velocity*/;
   };
   //-----------
   V_d PHI(const V_d &T, const V_d &Y, const V_d &Z) {
      V_d ret(T.size(), 0.);
      for (auto i = 0; i < T.size(); i++)
         ret[i] = PHI(T[i], Y[i], Z[i]);
      return ret;
   };
   //--------------------------------------------------
   double PHI_TN(const double T, const double Y, const double Z, const int N) {
      double ret(0.), akj, alp, hyp;
      ;

      for (auto k = 0; k < max + 1; k++)
         for (auto j = 2 - k % 2; j < max + 1; j++) {
            akj = a[j][k];
            alp = std::sqrt(j * j * p * p + k * k * q * q);
            hyp = std::cosh(alp * Z) + std::sinh(alp * Z) * std::tanh(alp * D);
            switch (N) {
               case 0:
                  ret = ret + akj * hyp * std::sin(j * T) * std::cos(k * Y);
               case 1:
                  ret = ret + j * akj * hyp * std::cos(j * T) * std::cos(k * Y);
               case 2:
                  ret = ret + -j * j * akj * hyp * std::sin(j * T) * std::cos(k * Y);
               case 3:
                  ret = ret + -j * j * j * akj * hyp * std::cos(j * T) * std::cos(k * Y);
               case 4:
                  ret = ret + j * j * j * j * akj * hyp * std::sin(j * T) * std::cos(k * Y);
               case 5:
                  ret = ret + j * j * j * j * j * akj * hyp * std::cos(j * T) * std::cos(k * Y);
               case 6:
                  ret = ret + -j * j * j * j * j * j * akj * hyp * std::sin(j * T) * std::cos(k * Y);
               case 7:
                  ret = ret + -j * j * j * j * j * j * j * akj * hyp * std::cos(j * T) * std::cos(k * Y);
            }
         }
      return ret;
   };
   double PHI_T(const double T, const double Y, const double Z) {
      return PHI_TN(T, Y, Z, 1);
   };
   V_d PHI_TN(const V_d &T, const V_d &Y, const V_d &Z, const int N) {
      V_d ret(T.size(), 0.);
      for (auto i = 0; i < T.size(); i++)
         ret[i] = PHI_TN(T[i], Y[i], Z[i], N);
      return ret;
   };
   std::vector<V_d> PHI_TN_mat(const V_d &T, const V_d &Y, const V_d &Z, const int N) {
      std::vector<V_d> ret(N, V_d(T.size(), 0.));
      for (auto n = 0; n < N; n++)
         for (auto l = 0; l < T.size(); l++)
            ret[n][l] = PHI_TN(T[l], Y[l], Z[l], n);
      return ret;
   };
};
#endif
