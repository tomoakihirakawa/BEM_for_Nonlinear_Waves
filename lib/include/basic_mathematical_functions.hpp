#ifndef basic_mathematical_functions_H
#define basic_mathematical_functions_H
#pragma once

#include "basic_alias.hpp"
#include <random>

// double RandomReal(const Tdd &minmax) {
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<double> dis(std::get<0>(minmax), std::get<1>(minmax));
//    return dis(gen);
// };

inline double RandomReal(double min, double max) {
  static std::mt19937 gen(std::random_device{}());
  static std::uniform_real_distribution<double> dis; // default [0,1)
  dis.param(std::uniform_real_distribution<double>::param_type(min, max));
  return dis(gen);
}

inline double RandomReal(const Tdd &minmax) { return RandomReal(std::get<0>(minmax), std::get<1>(minmax)); };

inline double RandomReal() { return RandomReal({0., 1.}); };

inline Tdd RandomRealTdd(double min, double max) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(min, max);
  return {{dis(gen), dis(gen)}};
};
inline Tddd RandomRealTddd(double min, double max) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(min, max);
  return {{dis(gen), dis(gen), dis(gen)}};
};
inline T3Tdd RandomRealT3Tdd(double min, double max) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(min, max);
  return {{{dis(gen), dis(gen)}, {dis(gen), dis(gen)}, {dis(gen), dis(gen)}}};
};

#endif
