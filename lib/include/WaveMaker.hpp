#ifndef WaveMaker_H
#define WaveMaker_H
#pragma once

#include "fundamental.hpp"

struct WaveMaker {
   double start = -1.5;
   double move_amplitude = .005;
   double s = M_PI / 2.;
   double a = move_amplitude;
   double k = M_PI / 1.;
   auto g = 9.81;
   T6d move_dir = {1., 0., 0., 0, 0, 0};
   WaveMaker(){};
};

namespace forced_motion {
double start = -1.5;
double move_amplitude = .005;
double s = M_PI / 2.;
double a = move_amplitude;
double k = M_PI / 1.;
auto g = 9.81;
T6d move_dir = {1., 0., 0., 0, 0, 0};
#if defined(experiment_sawai)
auto h = 0.1;
auto L = 0.25;
#elif defined(experiment_Li2002)
double h = 0.3048;
double H = 0.3 * h;
double x = 0;
double c = std::sqrt(g * (H + h));
#elif defined(WenWenLi2002_MachReflection)
double h = 0.06;
double H = 0.4 * h;
double x = 0;
double c = std::sqrt(g * (H + h));
#elif defined(experiment_)
double h = 0.3048;
double H = 0.4 * h;
double x = 0;
double c = std::sqrt(g * (H + h));
#endif
Tddd translation(const double t) {
#ifdef XueAndLin2011_large_amplitude
   auto w = 3.5317;
   Tddd move_dir = {1., 0., 0.};
   move_amplitude = 0.1;
   if (t > start)
      return -move_amplitude * std::cos(w * (t - start)) * move_dir;
   else
      return {0., 0., 0.};
#elif defined(experiment_sawai)
   auto w = std::sqrt(M_PI * g / L * std::tanh(M_PI * h / L));
   Tddd move_dir = {1., 0., 0.};
   return move_amplitude * std::sin(w * t) * move_dir;
#elif defined(experiment_Li2002)
   // 使わない
   Tddd move_dir = {1., 0., 0.};
   return move_dir;
#elif defined(WenWenLi2002_MachReflection)
   Tddd move_dir = {1., 0., 0.};
   return move_dir;
#else
   /* ------------------------------------------------------ */
   Tddd move_dir = {std::cos(k * t), std::sin(k * t), 0.};
   return move_amplitude * std::exp(-t) * (std::sin(k * t - s) - std::sin(-s)) * move_dir;
      /* ------------------------------------------------------ */
      // Tddd move_dir = Normalize(Tddd{1., 1., 0.});
      // return move_amplitude * std::exp(-t) * (sin(k * t - s) - std::sin(-s)) * move_dir;
#endif
};

T6d velocity(const double t) {
#ifdef XueAndLin2011_large_amplitude
   auto w = 3.5317;
   move_amplitude = 0.1;
   if (t > start)
      return move_amplitude * w * std::sin(w * (t - start)) * move_dir;
   else
      return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_sawai)
   auto w = std::sqrt(M_PI * g / L * tanh(M_PI * h / L));
   if (t > start)
      return -move_amplitude * w * std::sin(w * (t - start)) * move_dir;
   else
      return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_Li2002)
   double kappa = std::sqrt(3. * H / (4. * h * h * h));
   double tmp = cosh(kappa * (x - c * (t - start)));
   double eta = H * std::pow(1. / tmp, 2.);
   double u = c * eta / (h + eta);
   return u * move_dir;
#elif defined(rotation_test)
   double T = .25;
   if (t > start)
      return {0., 0., 0., 20. * M_PI / 180. * std::sin(M_PI / T * (t - start)), 0., 0.};
   else
      return {0., 0., 0., 0., 0., 0.};
#elif defined(WenWenLi2002_MachReflection)
   double eta = H * std::pow(1. / cosh(std::sqrt(3. * H / (4. * std::pow(h, 3))) * (x - c * (t - start))), 2.);
   double u = c * eta / (h + eta);
   if (t > start)
      return u * move_dir;
   else
      return {0., 0., 0., 0., 0., 0.};
#else
   /* ------------------------------------------------------ */
   T6d move_dir = {std::cos(k * t), std::sin(k * t), 0., 0., 0., 0.};
   T6d ddt_move_dir = {-k * std::sin(k * t), k * std::cos(k * t), 0., 0., 0., 0.};
   // /* |U|*n_p . n_surface = phin <-- given
   auto tmp = (-move_amplitude * std::exp(-t) * (std::sin(k * t - s) - std::sin(-s)) + move_amplitude * std::exp(-t) * (std::cos(k * t - s) * k)) * move_dir;
   tmp += move_amplitude * std::exp(-t) * (std::sin(k * t - s) - std::sin(-s)) * ddt_move_dir;
   return tmp;
#endif
};

T6d acceleration(const double t) {
#ifdef XueAndLin2011_large_amplitude
   auto w = 3.5317;
   move_amplitude = 0.1;
   if (t > start)
      return move_amplitude * w * w * std::cos(w * (t - start)) * move_dir;
   else
      return {0., 0., 0., 0., 0., 0.};
#elif defined(experiment_sawai)
   auto w = std::sqrt(M_PI * g / L * tanh(M_PI * h / L));
   return -move_amplitude * w * w * std::cos(w * t) * move_dir;
#elif defined(experiment_Li2002)
   auto u = (2 * std::sqrt(3) * std::pow(c, 2) * h * H * std::sqrt(H / std::pow(h, 3)) * std::sinh(std::sqrt(3) * std::sqrt(H / std::pow(h, 3)) * (-(c * t) + x)));
   u /= pow(h + 2 * H + h * std::cosh(std::sqrt(3) * std::sqrt(H / std::pow(h, 3)) * (-(c * t) + x)), 2);
   return u * move_dir;
#elif defined(WenWenLi2002_MachReflection)
   auto u = (2 * std::sqrt(3) * std::pow(c, 2) * h * H * std::sqrt(H / std::ow(h, 3)) * std::sinh(std::sqrt(3) * std::sqrt(H / std::pow(h, 3)) * (-(c * (t - start)) + x)));
   u /= pow(h + 2 * H + h * std::cosh(std::sqrt(3) * std::sqrt(H / std::pow(h, 3)) * (-(c * (t - start)) + x)), 2);
   return u * move_dir;
#else
   return {0., 0., 0., 0., 0., 0.};
#endif
};
}  // namespace forced_motion

#endif