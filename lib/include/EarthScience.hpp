#include "fundamental.hpp"

///////////////////////////////////////////////////////////
#ifndef INCL_EarthScience
   #define INCL_EarthScience

// Okada(1986), eq.(25),(26)
std::vector<double> fOkada(const std::vector<double>& U,
                           const double y,
                           const std::vector<double>& xi_eta,
                           const double del,
                           const double d,
                           const double mu,
                           const double lam) {
   double cosdel = std::cos(del);
   double sindel = std::sin(del);
   double q = y * sindel - d * cosdel;
   double tily = xi_eta[1] * cosdel + q * sindel;
   double tild = xi_eta[1] * sindel - q * cosdel;
   double X = Norm({xi_eta[0], q});
   double R = Norm({xi_eta[0], q, xi_eta[1]});
   double I5 = (std::abs(xi_eta[0]) < 10E-10) ? 0 : mu / (lam + mu) * 2 / cosdel * std::atan((xi_eta[1] * (X + q * cosdel) + X * (R + X) * sindel) / (xi_eta[0] * (R + X) * cosdel));
   double I4 = mu / (lam + mu) / cosdel * (std::log(R + tild) - sindel * std::log(R + xi_eta[1]));
   double I3 = mu / (lam + mu) * (tily / (cosdel * R + tild) - std::log(R + xi_eta[1])) + std::tan(del) * I4;
   double I2 = mu / (lam + mu) * (-std::log(R + xi_eta[1])) - I3;
   double I1 = mu / (lam + mu) * (-xi_eta[0] / (cosdel * (R + tild))) - std::tan(del) * I5;
   double tmp = (std::abs(q) < 10E-10) ? 0. : std::atan(xi_eta[0] * xi_eta[1] / (q * R));
   return {
       (-U[0] * (xi_eta[0] * q / (R * (R + xi_eta[1])) + tmp + I1 * sindel) - U[1] * (q / R - I3 * sindel * cosdel)) / (2 * M_PI),
       (-U[0] * (tily * q / (R * (R + xi_eta[1])) + q * cosdel / (R + xi_eta[1]) + I2 * sindel) - U[1] * (tily * q / (R * (R + xi_eta[0])) + cosdel * tmp - I1 * sindel * cosdel)) / (2 * M_PI),
       (-U[0] * (tild * q / (R * (R + xi_eta[1])) + q * sindel / (R + xi_eta[1]) + I4 * sindel) - U[1] * (tild * q / (R * (R + xi_eta[0])) + sindel * tmp - I5 * sindel * cosdel)) / (2 * M_PI)};
};

// Okada(1986), eq.(24)
std::vector<double> IntegrateStrikeDipOkada(const std::vector<double>& U,
                                            const std::vector<double>& xyIN,
                                            const std::vector<double>& xyIN_orign,
                                            const std::vector<double>& LW,
                                            const std::vector<double>& SoukouKeisyaDeg,
                                            const double d, const double mu, const double lam) {
   std::vector<double> SKRad = SoukouKeisyaDeg / 180. * M_PI;
   std::vector<std::vector<double>> A = {{std::cos(SKRad[0]), -std::sin(SKRad[0])}, {std::sin(SKRad[0]), std::cos(SKRad[0])}};
   std::vector<double> xy = Dot(A, xyIN - xyIN_orign);
   double x = xy[1];
   double y = -xy[0];
   double p = y * std::cos(SKRad[1]) + d * std::sin(SKRad[1]);
   std::vector<double> ret = fOkada(U, y, {x, p + LW[1]}, SKRad[1], d, mu, lam) - fOkada(U, y, {x, p}, SKRad[1], d, mu, lam) - fOkada(U, y, {x - LW[0], p + LW[1]}, SKRad[1], d, mu, lam) + fOkada(U, y, {x - LW[0], p}, SKRad[1], d, mu, lam);
   return {-ret[1], ret[0], ret[2]};
};

// #ifndef INCL_MOTION
// #include "motion.h"cd
// #endif

// class AidaFaultModel:public Motion{
//  public:
//   double w;
//   const std::vector<double> u1{0., 0.0068, 0.};
//   const std::vector<double> u2{0., 0.0046, 0.};
//   const std::vector<double> u3{0., 0.0068, 0.};
//  AidaFaultModel(const double w_IN): w(w_IN){str="AidaFaultModel";};
//   // double ret_phi_nt(const double x, const double y, const double t, const int n){
//   //   double f = w/cosh(w*(t-3.))/2.;
//   //   std::vector<double> tmp = IntegrateStrikeDipOkada(u1*f,{x,y},{-16.1039, 37.7521},{35., 35.},{-15., 20.},13.,1.,1.) + IntegrateStrikeDipOkada(u2*f,{x,y},{-23.7321, 8.88264},{35., 35.},{15., 20.},13.,1.,1.) + IntegrateStrikeDipOkada(u3*f,{x,y},{-33.0555, -25.5369},{35., 35.},{15., 20.},13.,1.,1.);
//   //   return tmp[2];
//   // };

//   double ret_phi_nt(const double r, const double yy, const double t, const int n){
//     double x = r * std::cos(160.*M_PI/180.) + 63.8479;
//     double y = r * std::sin(160.*M_PI/180.) - 31.7988;
//     double f = w/cosh(w*(t-60.))/2.;
//     if(r<10)
//       f=r/10.*r/10.*f;
//     std::vector<double> tmp = IntegrateStrikeDipOkada(u1*f,{x,y},{-16.1039, 37.7521},{35., 35.},{-15., 20.},13.,1.,1.) + IntegrateStrikeDipOkada(u2*f,{x,y},{-23.7321, 8.88264},{35., 35.},{15., 20.},13.,1.,1.) + IntegrateStrikeDipOkada(u3*f,{x,y},{-33.0555, -25.5369},{35., 35.},{15., 20.},13.,1.,1.);
//     return -tmp[2];
//   };
// };

/* ex.) PLOT 3D DIFOMATION OF SURFACE PLANE z = 0
  GNUPLOT PLOT;
  int N = 100;
  std::vector<std::vector<std::vector<double>>> data;
  for(auto i=0; i<N; i++)
    {
      std::vector<std::vector<double>> subdata;
      double x = 2.*(2.*(((double)i - N/2.)));
      for(auto j=0; j<N; j++)
        {
          double y = 2.*(2.*(((double)j - N/2.)));
          std::vector<double> xyz{x,y,0.};
          std::vector<double> u1{0., 0.0068,0.};
          std::vector<double> u2{0., 0.0046,0.};
          std::vector<double> u3{0., 0.0068,0.};
          subdata.push_back(xyz + IntegrateStrikeDipOkada(u1,{x,y},{-16.1039, 37.7521},{35., 35.},{-15., 20.},13.,1.,1.)
                            + IntegrateStrikeDipOkada(u2,{x,y},{-23.7321, 8.88264},{35., 35.},{15., 20.},13.,1.,1.)
                            + IntegrateStrikeDipOkada(u3,{x,y},{-33.0555, -25.5369},{35., 35.},{15., 20.},13.,1.,1.));

        }
      data.push_back(subdata);
    }
  int counter = 1;
  PLOT.Set({{"contour","base"}});
  PLOT.Splot(data,{{"w","l"}});
  PLOT.Clear();
*/
#endif
