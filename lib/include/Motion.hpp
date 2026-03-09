// motion.h
#ifndef INCL_MOTION
#define INCL_MOTION
#include "fundamental.hpp"
/////////////////////////////////////////////////////
struct Motion {
   virtual double ret_phi_nt(const std::vector<double>& normal_vectorIN, const double t) {
      return 0.;
   };
};
/////////////////////////////////////////////////////
struct GoringBoussinesq : public Motion {
   double H, h, t0;
   std::string name;
   std::vector<double> motion_dir_xyz;
   GoringBoussinesq(const std::vector<double>& motion_dir_xyzIN, const double WaveHeight, const double Depth, const double LeadingTime) : H(WaveHeight),
                                                                                                                                          h(Depth),
                                                                                                                                          t0(LeadingTime),
                                                                                                                                          motion_dir_xyz(motion_dir_xyzIN / Norm(motion_dir_xyzIN) /*normalized*/),
                                                                                                                                          name("GoringBoussinesq"){};

   double ret_phi_nt(const std::vector<double>& normal_vectorIN, const double t) {
      double k = std::sqrt(3. * H / (4. * h * h * h));
      double C = std::sqrt(9.81 * (h + H));
      double eta = H * pow(cosh(k * (-C * (t - t0))), -2.);
      return -C * eta / (h + eta) * Dot(motion_dir_xyz, normal_vectorIN / Norm(normal_vectorIN) /*normalized*/);
   };
};
#endif
