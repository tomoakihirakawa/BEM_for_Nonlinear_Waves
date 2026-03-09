#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>

inline double _GRAVITY_ = 9.81;
#include "TanakaSolitaryWave.hpp"

int main() {
   double h = 1.0;
   double bz = 0.0;

   // --- 1. Wave profiles for various H/h ---
   {
      double Hh_vals[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
      std::vector<TanakaSolitaryWave> waves(7);
      for (int i = 0; i < 7; ++i)
         waves[i].solve(Hh_vals[i], h, bz, 0.0, 80);

      std::ofstream ofs("tanaka_profiles.csv");
      ofs << "x";
      for (double Hh : Hh_vals) ofs << ",eta_Hh" << Hh;
      ofs << "\n";
      for (double x = -15; x <= 15; x += 0.05) {
         ofs << x;
         for (int i = 0; i < 7; ++i) {
            std::array<double, 3> X = {x, 0, 0};
            ofs << "," << (waves[i].eta(X, 0) - h - bz);
         }
         ofs << "\n";
      }
      std::cout << "Wrote tanaka_profiles.csv" << std::endl;
   }

   // --- 2. Velocity field for H/h = 0.4 ---
   {
      TanakaSolitaryWave w;
      w.solve(0.4, h, bz);

      std::ofstream ofs("tanaka_velocity.csv");
      ofs << "x,z,u,w,phi\n";
      for (double x = -8; x <= 8; x += 0.2) {
         std::array<double, 3> Xsurf = {x, 0, 0};
         double z_surf = w.eta(Xsurf, 0);
         for (double zfrac = 0.0; zfrac <= 1.001; zfrac += 0.1) {
            double z = bz + zfrac * (z_surf - bz);
            std::array<double, 3> X = {x, 0, z};
            auto v = w.gradPhi(X, 0);
            ofs << x << "," << z << "," << v[0] << "," << v[2] << "," << w.phi(X, 0) << "\n";
         }
      }
      std::cout << "Wrote tanaka_velocity.csv" << std::endl;
   }

   // --- 3. Tanaka vs KdV sech² ---
   {
      TanakaSolitaryWave w1, w2;
      w1.solve(0.3, h, bz, 0.0, 80);
      w2.solve(0.6, h, bz, 0.0, 80);

      std::ofstream ofs("tanaka_vs_kdv.csv");
      ofs << "x,tanaka_03,kdv_03,tanaka_06,kdv_06\n";
      for (double x = -15; x <= 15; x += 0.05) {
         std::array<double, 3> X = {x, 0, 0};
         double t1 = w1.eta(X, 0) - h;
         double kappa1 = std::sqrt(3.0 * 0.3 / (4.0));
         double kdv1 = 0.3 * h / (std::cosh(kappa1 * x) * std::cosh(kappa1 * x));

         double t2 = w2.eta(X, 0) - h;
         double kappa2 = std::sqrt(3.0 * 0.6 / (4.0));
         double kdv2 = 0.6 * h / (std::cosh(kappa2 * x) * std::cosh(kappa2 * x));

         ofs << x << "," << t1 << "," << kdv1 << "," << t2 << "," << kdv2 << "\n";
      }
      std::cout << "Wrote tanaka_vs_kdv.csv" << std::endl;
   }

   // --- 4. F² vs H/h ---
   {
      std::ofstream ofs("tanaka_F2.csv");
      ofs << "Hh,F2,F2_kdv,c,qc\n";
      for (double Hh = 0.05; Hh <= 0.78; Hh += 0.025) {
         TanakaSolitaryWave w;
         w.solve(Hh, h, bz, 0.0, 80);
         ofs << Hh << "," << w.F2 << "," << (1.0 + Hh) << "," << w.phase_speed << "," << w.qc << "\n";
      }
      std::cout << "Wrote tanaka_F2.csv" << std::endl;
   }

   return 0;
}
