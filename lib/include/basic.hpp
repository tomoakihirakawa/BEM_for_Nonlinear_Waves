#pragma once

#include "basic_constants.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>
#include <variant>
#include <vector>

/* ---------------------------------------------------------------------- */

// double real_sph(const int l, const int m, const double theta, const double psi) {
//    return std::sph_legendre(std::abs(l), std::abs(l), theta) * std::cos(psi);
// };
// double real_sph_scale_ommited(const int l, const int m, const double theta, const double psi) {
//    return real_sph(l, m, theta, psi) / std::sqrt((2. * std::abs(l) + 1.) / (4. * M_PI));
// };
/* ------------------------------------------------------ */
/*              表示方法に関するライブラリ                     */
/* ------------------------------------------------------ */
// 2021/08/02追加
using Var_SID = std::variant<std::string, int, std::size_t, double, std::vector<double>, std::vector<int>, Tddd, T6d>;
using V_Var_SID = std::vector<Var_SID>;
using VV_Var_SID = std::vector<std::vector<Var_SID>>;

inline std::string GridVector(const std::vector<int> &V, const int i = 20) {
  std::stringstream ss;
  for (const auto &v : V) {
    ss << std::setw(i) << std::right;
    ss << v;
  }
  return ss.str();
};
inline std::string GridVector(const std::vector<double> &V, const int i = 20) {
  std::stringstream ss;
  for (const auto &v : V) {
    ss << std::setw(i) << std::right;
    ss << v;
  }
  return ss.str();
};

inline std::string Grid(const V_Var_SID &V, const int i = 20) {
  std::stringstream ss;
  for (const auto &var : V) {
    ss << std::setw(i) << std::right;
    if (std::holds_alternative<int>(var))
      ss << std::get<int>(var);
    else if (std::holds_alternative<size_t>(var))
      ss << std::get<size_t>(var);
    else if (std::holds_alternative<double>(var))
      ss << std::get<double>(var);
    else if (std::holds_alternative<std::string>(var))
      ss << std::get<std::string>(var);
    else if (std::holds_alternative<Tddd>(var)) {
      auto v = std::get<Tddd>(var);
      std::stringstream s;
      s << "{" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << ","
        << "}";
      ss << s.str();
    } else if (std::holds_alternative<T6d>(var)) {
      auto v = std::get<T6d>(var);
      std::stringstream s;
      s << "{" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << "," << std::get<3>(v) << "," << std::get<4>(v) << "," << std::get<5>(v) << ","
        << "}";
      ss << s.str();
    } else if (std::holds_alternative<std::vector<double>>(var)) {
      std::vector<double> v = std::get<std::vector<double>>(var);
      std::stringstream s;
      s << "{";
      if (!v.empty()) {
        for (size_t i = 0; i < v.size() - 1; i++)
          s << v[i] << ",";
        s << *v.rbegin();
      }
      s << "}";
      ss << s.str();
    } else if (std::holds_alternative<std::vector<int>>(var)) {
      std::vector<int> v = std::get<std::vector<int>>(var);
      std::stringstream s;
      s << "{";
      if (!v.empty()) {
        for (size_t i = 0; i < v.size() - 1; i++)
          s << v[i] << ",";
        s << *v.rbegin();
      }
      s << "}";
      ss << s.str();
    }
    // else if (std::holds_alternative<std::vector<Tddd>>(var))
    // {
    // 	auto v = std::get<std::vector<Tddd>>(var);
    // 	std::stringstream s;
    // 	s << "{";
    // 	if (!v.empty())
    // 	{
    // 		for (size_t i = 0; i < v.size() - 1; i++)
    // 			s << (Tddd)v[i] << ",";
    // 		s << *v.rbegin();
    // 	}
    // 	s << "}";
    // 	ss << s.str();
    // }
  }
  return ss.str();
};

inline std::string GridGrid(const VV_Var_SID &VV, const int i = 10) {
  std::stringstream ss;
  for (const auto &V : VV) {
    ss << Grid(V, i);
    ss << std::endl;
  }
  return ss.str();
};

/* ------------------------------------------------------ */
template <typename T, typename U> std::map<T, U> KeyTake(std::map<T, U> &m, const std::vector<T> &V) {
  std::map<T, U> ret;
  for (const auto &v : V)
    if (m.find(v) != m.end())
      ret[v] = m[v];
  return ret;
};
template <typename T, typename U> std::vector<T> TakeFirst(const std::map<T, U> &m) {
  std::vector<T> ret(m.size());
  int i = 0;
  for (auto it = m.begin(); it != m.end(); ++it)
    ret[i++] = it->first;
  return ret;
};
template <typename T, typename U> std::vector<U> TakeSecond(const std::map<T, U> &m) {
  std::vector<U> ret(m.size());
  int i = 0;
  for (auto it = m.begin(); it != m.end(); ++it)
    ret[i++] = it->second;
  return ret;
};
template <typename T, typename U> std::vector<std::vector<U>> TakeSecond(const std::map<T, std::vector<U>> &m) {
  std::vector<std::vector<U>> ret(m.size());
  int i = 0;
  for (auto it = m.begin(); it != m.end(); ++it)
    ret[i++] = it->second;
  return ret;
};

///////////////////////////////////////////////////////
// 20200526
inline V_i stoi(const V_s &vec) {
  V_i ret(vec.size());
  std::transform(vec.begin(), vec.end(), ret.begin(), [](const std::string &str) { return std::stoi(str); });
  return ret;
};
inline VV_i stoi(const std::vector<V_s> &vec) {
  VV_i ret(0);
  for (const auto &v : vec)
    ret.emplace_back(stoi(v));
  return ret;
};
inline VVV_i stoi(const std::vector<std::vector<V_s>> &vec) {
  VVV_i ret(0);
  for (const auto &v : vec)
    ret.emplace_back(stoi(v));
  return ret;
};
inline std::vector<bool> stob(const V_s &vec) {
  std::vector<bool> ret(vec.size());
  std::transform(vec.begin(), vec.end(), ret.begin(), [](const std::string &str) {
    // Handle both "true"/"false" and "0"/"1" formats
    if (str == "true" || str == "1")
      return true;
    if (str == "false" || str == "0")
      return false;
    // Fallback: try boolalpha parsing
    bool b = false;
    std::istringstream(str) >> std::boolalpha >> b;
    return b;
  });
  return ret;
};
inline V_d stod(const V_s &vec) {
  V_d ret;
  ret.reserve(vec.size());
  for (const auto &str : vec) {
    ret.push_back(std::stod(str));
  }
  return ret;
}

inline VV_d stod(const std::vector<V_s> &vec) {
  VV_d ret(0);
  for (const auto &v : vec)
    ret.emplace_back(stod(v));
  return ret;
};
inline VVV_d stod(const std::vector<std::vector<V_s>> &vec) {
  VVV_d ret(0);
  for (const auto &v : vec)
    ret.emplace_back(stod(v));
  return ret;
};
///////////////////////////////////////////////////////
inline int InverseQuotientMod(const V_i &row_col, const int b) { return row_col[0] * b + row_col[1]; };
inline V_i QuotientMod(const int a, const int b) { return {(int)a / b, a % b}; };
//==========================================================

template <class T> std::vector<T> Take(const std::vector<T> &vec, const std::array<int, 3> &range = {0, -1, 1}) {
  const int start = range[0];
  int end = (range[1] < 0) ? vec.size() + range[1] : range[1];
  const int skip = range[2];

  if (start < 0 || start >= vec.size())
    throw std::out_of_range("Start index out of range");
  if (end > vec.size() || end <= start)
    throw std::out_of_range("End index out of range or not logical");
  if (skip <= 0)
    throw std::invalid_argument("Skip value must be positive");

  std::vector<T> ret;
  ret.reserve(vec.size()); // Optional optimization
  for (int i = start; i < end; i += skip) {
    ret.emplace_back(vec[i]);
  }
  return ret;
}

//==========================================================
inline VV_d Import(const std::string &fname) {
  std::ifstream fin(fname);
  VV_d mat;
  std::string str;
  double tmp;
  while (getline(fin, str)) {
    std::stringstream ss;
    ss << str;
    V_d row;
    while (ss >> tmp)
      row.emplace_back(tmp);
    mat.emplace_back(row);
  }
  fin.close();
  return mat;
};
// example
// std::vector<std::vector<double>> mat=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");

#include "basic_IO.hpp"
#include "basic_exception.hpp"
//==========================================================
template <class T, class U> bool MemberQ(const std::map<T, U> &list, const T &form) { return list.find(form) != list.end(); };
template <class T, class U> bool MemberQ(const std::map<T *, U> &list, const T *form) { return list.find(form) != list.end(); };
template <class T, class U> bool AllMemberQ(const std::map<T *, U> &list, const std::vector<T *> &form) {
  if (form.empty())
    throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "empty"));
  for (const auto &f : form)
    if (!MemberQ(list, f))
      return false;
  return true;
};
template <class T, class U> bool AnyMemberQ(const std::map<T *, U> &list, const std::vector<T *> &form) {
  if (form.empty())
    throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "empty"));
  for (const auto &f : form)
    if (MemberQ(list, f))
      return true;
  return false;
};
/* ------------------------------------------------------ */
template <class T> int Count(const std::vector<T *> &list, const T *const v) { return std::count(list.cbegin(), list.cend(), v); };
/* -------------------------------------------------------------------------- */
template <typename Container, typename T> bool MemberQ(const Container &list, const T &form) { return (std::find(list.cbegin(), list.cend(), form) != list.cend()); };
template <typename T> bool MemberQ(const std::unordered_set<T> &list, const T &form) { return list.contains(form); };
template <typename T, typename... Ts> bool MemberQ(const std::tuple<T, Ts...> &list, const T &form) {
  return std::apply([&form](const auto &...elems) { return ((elems == form) || ...); }, list);
};
template <typename T, typename... Ts> bool MemberQ(const std::tuple<T *, Ts...> &list, const T *const form) {
  return std::apply([&form](const auto &...elems) { return ((elems == form) || ...); }, list);
};
/* -------------------------------------------------------------------------- */
template <typename Container, typename T> bool MemberQ_(const Container &list, const T &form) { return (std::find(list.cbegin(), list.cend(), form) != list.cend()); };
template <typename T> bool MemberQ_(const std::unordered_set<T> &list, const T &form) { return list.contains(form); };
template <typename T, typename... Ts> bool MemberQ_(const std::tuple<T, Ts...> &list, const T &form) {
  return std::apply([&form](const auto &...elems) { return ((elems == form) || ...); }, list);
};
template <typename T, typename... Ts> bool MemberQ_(const std::tuple<T *, Ts...> &list, const T *const form) {
  return std::apply([&form](const auto &...elems) { return ((elems == form) || ...); }, list);
};
/* -------------------------------------------------------------------------- */
template <class T> bool AnyTrue(const std::vector<T> &vec) {
  for (const auto &o : vec)
    if (o)
      return true;
  return false;
}

template <class T> bool AllTrue(const std::vector<T> &vec) {
  for (const auto &o : vec)
    if (!o)
      return false;
  return true;
}

inline bool AllTrue(const T8b &v) { return std::get<0>(v) && std::get<1>(v) && std::get<2>(v) && std::get<3>(v) && std::get<4>(v) && std::get<5>(v) && std::get<6>(v) && std::get<7>(v); }
/* -------------------------------------------------------------------------- */
template <typename T> std::vector<T> TakeExcept(const std::vector<T> &list, const T &form) {
  std::vector<T> ret;
  ret.reserve(list.size());
  for (const auto &l : list)
    if (l != form)
      ret.emplace_back(l);
  return ret;
};
template <typename T> std::vector<T> TakeExcept(const std::vector<T> &list, const std::vector<T> &form) {
  std::vector<T> ret;
  ret.reserve(list.size());
  for (const auto &l : list)
    if (!MemberQ(form, l))
      ret.emplace_back(l);
  return ret;
};
//
template <typename T> std::unordered_set<T> TakeExcept(const std::unordered_set<T> &list, const T &form) {
  std::unordered_set<T> ret;
  ret.reserve(list.size());
  for (const auto &l : list)
    if (l != form)
      ret.emplace_back(l);
  return ret;
};
template <typename T> std::unordered_set<T> TakeExcept(const std::unordered_set<T> &list, const std::vector<T> &form) {
  std::unordered_set<T> ret;
  ret.reserve(list.size());
  for (const auto &l : list)
    if (!MemberQ(form, l))
      ret.emplace(l);
  return ret;
};
template <typename T> std::vector<T> TakeExceptIndex(const std::vector<T> &list, const int index) {
  std::vector<T> ret;
  ret.reserve(list.size());
  // for (auto i = 0; (i < list.size() && i < index); i++)
  //   ret.emplace_back(list[i]);
  // for (auto i = index + 1; i < list.size(); i++)
  //   ret.emplace_back(list[i]);
  for (auto i = 0; i < list.size(); i++)
    if (i != index)
      ret.emplace_back(list[i]);
  return ret;
};
template <typename T> std::vector<T> TakeExceptIndecies(const std::vector<T> &list, const V_i &form) {
  std::vector<T> ret;
  ret.reserve(list.size());
  for (auto i = 0; i < list.size(); i++) {
    if (!MemberQ(form, i))
      ret.emplace_back(list[i]);
  }
  return ret;
};
template <class T> bool AnyMemberQ(const std::vector<T> &list, const std::vector<T> &form) { // Mathematica like
  for (const auto &f : form)
    if (MemberQ(list, f))
      return true;
  return false;
};
template <class T> bool AllMemberQ(const std::vector<T> &list, const std::vector<T> &form) { // Mathematica like
  for (const auto &f : form)
    if (!MemberQ(list, f))
      return false;
  return true;
};
// double Sum(const V_d::iterator first, const V_d::iterator second) { return std::reduce(first, second, 0.); };
// double Sum(const V_d &v) { return std::reduce(v.cbegin(), v.cend(), 0.); };
// double Total(const V_d &v) { return std::reduce(v.cbegin(), v.cend(), 0.); };

inline double Sum(const V_d::iterator first, const V_d::iterator second) { return std::accumulate(first, second, 0.); };
inline double Sum(const V_d &v) { return std::accumulate(v.cbegin(), v.cend(), 0.); };
// double Total(const V_d &v) { return std::accumulate(v.cbegin(), v.cend(), 0.); };

#include "basic_arithmetic_vector_operations.hpp"

// double Mean(const V_d &v) {
//    if (v.empty())
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
//    else if (v.size() == 1)
//       return v[0];
//    else
//       return std::accumulate(v.cbegin(), v.cend(), 0.) / v.size();
// };

inline V_d Sum(const VV_d &v) {
  if (v.empty())
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
  else if (v.size() == 1)
    return v[0];

  V_d ret(v[0].size(), 0.);
  for (const auto &u : v)
    ret += u;
  return ret;
};

inline V_d Total(const VV_d &v) {
  if (v.empty())
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "An empty vector is passed");
  else if (v.size() == 1)
    return v[0];

  V_d ret(v[0].size(), 0.);
  for (const auto &u : v)
    ret += u;
  return ret;
};

inline double RMS(const V_d &v) {
  double ret(0);
  for (const auto &u : v)
    ret += u * u;
  return std::sqrt(ret);
};

// template <typename T>
// std::vector<T> RandomSample(std::vector<T> ret /*copy may change*/) {
//    std::shuffle(std::begin(ret), std::end(ret), std::default_random_engine());
//    return ret;
//    // if (!n /*if n==0*/) {
//    //    return ret;
//    // } else {
//    //    return Take(ret, {0, n});
//    // }
// };
// template <typename T>
// std::vector<T *> RandomSample(std::vector<T *> ret /*copy may change*/) {
//    std::shuffle(std::begin(ret), std::end(ret), std::default_random_engine());
//    return ret;
//    // if (!n /*if n==0*/) {

//    // } else {
//    //    return Take(ret, {0, n});
//    // }
// };

// template <typename T>
// std::vector<T> RandomSample(std::vector<T> ret) {
//    // Use the current time as a seed for the random number generator
//    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
//    std::shuffle(std::begin(ret), std::end(ret), std::default_random_engine(seed));
//    return ret;
// }

// template <typename T>
// std::vector<T *> RandomSample(std::vector<T *> ret) {
//    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
//    std::shuffle(std::begin(ret), std::end(ret), std::default_random_engine(seed));
//    return ret;
// }

template <typename V> V RandomSample(V ret) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(std::begin(ret), std::end(ret), gen);
  return ret;
}

#include "basic_vectors.hpp"

//==========================================================
// struct Timer
// {
// 	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
// 	Timer()
// 	{
// 		start = std::chrono::high_resolution_clock::now();
// 	};
// 	std::chrono::duration<double> elapsed()
// 	{
// 		end = std::chrono::high_resolution_clock::now();
// 		return end - start;
// 	};
// };
#include "lib_measurement.hpp"
//==========================================================
inline double sgn(const double x) { return (x > 0.) ? 1. : (x < 0. ? -1. : 0.); };
//==========================================================
// template <class T>
// double Mean(const std::vector<T> &v)
// {
//   if (v.size() == 1)
//     return v[0];
//   return std::accumulate(v.cbegin(), v.cend(), 0.) / v.size();
// };
// V_d Mean(const VV_d &v)
// {
//   if (v.size() == 1)
//   {
//     return v[0];
//   }
//   VV_d W = Transpose(v);
//   V_d ret(W.size());
//   for (size_t i = 0; i < ret.size(); i++)
//     ret[i] = Mean(W[i]);
//   ;
//   return ret;
// };
// double Variance(const V_d &v)
// {
//   int N = v.size();
//   double ret(0);
//   for (const auto &w : v)
//     ret = ret + w * w;
//   return ret / N - std::pow(Mean(v), 2.);
// };
//==========================================================
inline double SgLog(const double u, const double hsig, const double maxmin) { return hsig + sgn(u) * std::exp(maxmin * (std::abs(u) - 1.)); };
inline double InvSgLog(const double h, const double hsig, const double maxmin) { return sgn(h - hsig) * std::log(std::abs(h - hsig)) / maxmin + 1.; };
inline V_d SgLog(const V_d &h, const double hsig, const double maxmin) {
  V_d ret;
  std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [hsig, maxmin](double tmp) { return SgLog(tmp, hsig, maxmin); });
  return ret;
};
inline double DSgLog(const double u, const double hsig, const double maxmin) { return maxmin * sgn(u) * std::exp(maxmin * (std::abs(u) - 1.)); }; // 単調増加
inline V_d DSgLog(const V_d &h, const double hsig, const double maxmin) {
  V_d ret;
  std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [hsig, maxmin](double tmp) { return DSgLog(tmp, hsig, maxmin); });
  return ret;
};
//============================================================
inline double Sg(const double h, const double h_sig, const double beta) { return sgn(h) * std::pow(std::abs(h), beta) + h_sig; };
inline double InvSg(const double h, const double h_sig, const double beta) { return sgn(h - h_sig) * std::pow(std::abs(h - h_sig), 1. / beta); };
inline V_d Sg(const V_d &h, const double h_sig, const double beta) {
  V_d ret;
  std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [h_sig, beta](double tmp) { return Sg(tmp, h_sig, beta); });
  return ret;
};
inline double DSg(const double h, const double h_sig, const double beta) { return beta * std::pow(std::abs(h), beta - 1.); }; // 単調増加
inline V_d DSg(const V_d &h, const double h_sig, const double beta) {
  V_d ret;
  std::transform(h.cbegin(), h.cend(), std::back_inserter(ret), [h_sig, beta](double tmp) { return DSg(tmp, h_sig, beta); });
  return ret;
};
//==========================================================
inline double linspace(const V_d &v, int size, int i) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  return v[0] + (v[1] - v[0]) * i / (double)(size - 1.);
};
/* double init_position(int size, int i, V_d v, double c, double beta){ */
/*   cout << "max = " << v[1] << ", min = " << v[0] << ", d = " << v[1]-v[0] << endl; */
/*   return Sg( v[0] + (v[1] - v[0]) * i/(double)(size - 1.), c, beta);; */
/* }; */
inline double linspace(const V_d &v, int size, int i, double c, double beta) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  V_d w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
  return Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
};
inline void linspace(const V_d &v, int size, V_d &O_X) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  O_X.resize(size);
  for (int i = 0; i < size; i++)
    O_X[i] = v[0] + (v[1] - v[0]) * i / ((double)size - 1.);
};
inline void linspace(const V_d &v, int size, V_d &O_X, double c, double beta) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  V_d w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
  O_X.resize(size);
  for (int i = 0; i < size; i++)
    O_X[i] = Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
};
inline V_d linspace(const V_d &v, int size) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  V_d O_X(size);
  for (int i = 0; i < size; i++)
    O_X[i] = v[0] + (v[1] - v[0]) * i / (double)(size - 1.);
  return O_X;
};
inline V_d linspace(const V_d &v, int size, double c, double beta) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  V_d w{InvSg(v[0], c, beta), InvSg(v[1], c, beta)};
  V_d O_X(size);
  for (int i = 0; i < size; i++)
    O_X[i] = Sg(w[0] + (w[1] - w[0]) * i / (double)(size - 1.), c, beta);
  return O_X;
};
inline void linspace(const std::vector<V_d> &mat, const V_i size, V_d &O_X) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  V_d v;
  for (size_t n = 0; n < mat.size(); n++) {
    linspace(mat[n], size[n], v);
    O_X.insert(std::end(O_X), std::begin(v), std::end(v));
  }
};
////////// Join has been modified on 20200411
template <typename T> std::tuple<T, T, T, T, T, T> Join(const std::tuple<T, T, T> &a, const std::tuple<T, T, T> &b) { return {std::get<0>(a), std::get<1>(a), std::get<2>(a), std::get<0>(b), std::get<1>(b), std::get<2>(b)}; };
template <typename T> std::unordered_set<T> Join(std::unordered_set<T> ret, const std::unordered_set<T> &b, const std::unordered_set<T> &c) {
  ret.insert(b.begin(), b.end());
  ret.insert(c.begin(), c.end());
  return ret;
};

template <typename T> std::unordered_set<T> Join(std::unordered_set<T> ret, const std::unordered_set<T> &b) {
  ret.insert(b.begin(), b.end());
  return ret;
};

template <typename T> std::unordered_set<T> Join(std::unordered_set<T> ret, const std::vector<T> &b) {
  ret.insert(b.begin(), b.end());
  return ret;
};

// for two vectors
template <typename T> std::vector<T> Join(std::vector<T> ret, const std::vector<T> &b) {
  ret.reserve(ret.size() + b.size());
  ret.insert(ret.end(), b.begin(), b.end());
  return ret;
};

//! for arbitrary number of vectors
template <typename T, typename... Ts> std::vector<T> Join(std::vector<T> ret, const std::vector<T> &b, const Ts &...args) {
  ret.reserve(ret.size() + b.size());
  ret.insert(ret.end(), b.begin(), b.end());
  return Join(ret, args...);
};

template <typename T> std::vector<std::vector<T>> Join(std::vector<std::vector<T>> ret, const std::vector<std::vector<T>> &b) {
  ret.reserve(ret.size() + b.size());
  ret.insert(ret.end(), b.begin(), b.end());
  return ret;
};
template <typename T> std::unordered_set<T> Append(std::unordered_set<T> a, const T &b) {
  a.emplace(b);
  return a;
};
template <typename T> std::vector<T> Append(std::vector<T> a, const T &b) {
  a.insert(a.end(), b);
  return a;
};
template <typename T, std::size_t N> std::array<T, N + 1> Append(std::array<T, N> a, const T &b) {
  std::array<T, N + 1> ret;
  std::size_t i = 0;
  for (const auto &o : a)
    ret[i++] = o;
  ret[i] = b;
  return ret;
};
//==========================================================
inline std::vector<std::tuple<double, double>> GaussianQuadratureWeightsTuple(const int n, const double x1, const double x2) {
  double EPS = 1E-14;
  std::tuple<double, double> zeros = {0., 0.};
  std::vector<std::tuple<double, double>> XW(n, zeros);
  double z1, z, xm, xl, pp, p3, p2, p1;
  int m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);
  for (auto i = 0; i < m; i++) {
    z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (auto j = 0; j < n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    } while (std::abs(z - z1) > EPS);
    std::get<0>(XW[i]) = xm - xl * z;
    std::get<0>(XW[n - 1 - i]) = xm + xl * z;
    std::get<1>(XW[i]) = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    std::get<1>(XW[n - 1 - i]) = std::get<1>(XW[i]);
  }
  return XW;
};
/* ------------------------------------------------------ */
inline VV_d GaussianQuadratureWeights(const int n, const double x1, const double x2) {
  double EPS = 1E-14;
  VV_d XW(n, V_d(2, 0));
  double z1, z, xm, xl, pp, p3, p2, p1;
  int m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);
  for (auto i = 0; i < m; i++) {
    z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (auto j = 0; j < n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    } while (std::abs(z - z1) > EPS);
    XW[i][0] = xm - xl * z;
    XW[n - 1 - i][0] = xm + xl * z;
    XW[i][1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    XW[n - 1 - i][1] = XW[i][1];
  }
  return XW;
};

inline VV_d SingularGaussianQuadratureWeights(const int n, const double x1, const double x2, const double xi_sing, const double b /*magnitude of singularity*/) {
  VV_d gw_T = Transpose(GaussianQuadratureWeights(n, InvSg(x1, xi_sing, b), InvSg(x2, xi_sing, b)));

  return Transpose(VV_d{Sg(gw_T[0], xi_sing, b), gw_T[1] * DSg(gw_T[0], xi_sing, b)});
};
//==========================================================
inline void gauleg(const double x1, const double x2, V_d &x, V_d &w) {
#if defined FULL_DEBUG
  std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;
#endif
  const double EPS = 1.0e-14;
  double z1, z, xm, xl, pp, p3, p2, p1;
  int n = x.size();
  int m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);
  for (int i = 0; i < m; i++) {
    z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 0; j < n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    } while (std::abs(z - z1) > EPS);
    x[i] = xm - xl * z;
    x[n - 1 - i] = xm + xl * z;
    w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];
  }
};

// iはNコのデータがある場合，0<=i<=N-1となる
// i=N-1の場合，q[N-1]=1.なのでh=1.はq[N-1]=1.<= h && h < q[N]
// \label{interpolation:Bspline}
inline double Bspline(const double h, const V_d &q, const int i, const int K) {
  switch (K) {
  case 1:
    return (q[i] <= h && h < q[i + 1]) ? 1. : 0.;
  case 2:
    return ((q[i] <= h && h < q[i + 1]) ? (h - q[i]) / (q[i + 1] - q[i]) : 0.) + ((q[i + 1] <= h && h < q[i + 2]) ? (q[i + 2] - h) / (q[i + 2] - q[i + 1]) : 0.);
  case 3:
    return ((q[i] <= h && h < q[i + 1]) ? (h - q[i]) / (q[i + 2] - q[i]) * (h - q[i]) / (q[i + 1] - q[i]) : 0.) + ((q[i + 1] <= h && h < q[i + 2]) ? ((h - q[i]) / (q[i + 2] - q[i]) * (q[i + 2] - h) + (q[i + 3] - h) / (q[i + 3] - q[i + 1]) * (h - q[i + 1])) / (q[i + 2] - q[i + 1]) : 0.) + ((q[i + 2] <= h && h < q[i + 3]) ? (q[i + 3] - h) / (q[i + 3] - q[i + 1]) * (q[i + 3] - h) / (q[i + 3] - q[i + 2]) : 0.);
  default:
    return (h - q[i]) / (q[i + K - 1] - q[i]) * Bspline(h, q, i, K - 1) + (q[i + K] - h) / (q[i + K] - q[i + 1]) * Bspline(h, q, i + 1, K - 1);
  };
};

inline V_d Bspline(const V_d &h, const V_d &q, const int i, const int K) {
  // {B_0(x),B_1(x),B_2(x),B_3(x),.......}
  V_d ret(h.size());
  for (auto k = 0; k < h.size(); k++)
    ret[k] = Bspline(h[k], q, i, K);
  return ret;
};

inline V_d Bspline(const double h, const V_d &q, const int K) {
  // {B_0(h),B_1(h),B_2(h),B_3(h),.......}
  V_d ret(q.size() - K);
  for (auto i = 0; i < q.size() - K; i++)
    ret[i] = Bspline(h, q, i, K);
  return ret;
};

inline VV_d Bspline(const V_d &H, const V_d &q, const int K) {
  // {{B_0(x0),B_1(x0),B_2(x0),B_3(x0),.......},
  //  {B_0(x1),B_1(x1),B_2(x1),B_3(x1),.......},
  //  {B_0(x2),B_1(x2),B_2(x2),B_3(x2),.......}}
  VV_d ret(H.size());
  std::transform(H.begin(), H.end(), ret.begin(), [q, K](const auto h) { return Bspline(h, q, K); });
  return ret;
};

inline double D_Bspline(const double h, const V_d &q, const int i, const int K) {
  // return (K - 1.) * (1. / (q[i + K - 1] - q[i]) * Bspline(h, q, i, K - 1) - 1. / (q[i + K] - q[i + 1]) * Bspline(h, q, i + 1, K - 1));
  return (K - 1.) * (Bspline(h, q, i, K - 1) / (q[i + K - 1] - q[i]) - Bspline(h, q, i + 1, K - 1) / (q[i + K] - q[i + 1]));
};

inline double D_Bspline(const double h, const V_d &q, const int i, const int K, const int n) {
  if (K == 1 && n > 0) {
    return 0.;
  } else {
    switch (n) {
    case 0:
      return Bspline(h, q, i, K);
    default:
      return (K - 1) * (D_Bspline(h, q, i, K - 1, n - 1) / (q[i + K - 1] - q[i]) - D_Bspline(h, q, i + 1, K - 1, n - 1) / (q[i + K] - q[i + 1]));
    };
  };
};

inline V_d D_Bspline(const double h, const V_d &q, const int K) {
  // {B_0(h),B_1(h),B_2(h),B_3(h),.......}
  V_d ret(q.size() - K);
  for (auto i = 0; i < q.size() - K; i++)
    ret[i] = D_Bspline(h, q, i, K);
  return ret;
};

#include "basic_linear_systems.hpp"

// V_d PeriodicKnot(const V_d& h, const int K){
//   int s(h.size());
//   V_d q(s + K);
//   double EPS = 1.E-15;
//   for(auto i=-K; i <= -1; i++)
//     q[i] = h[i+s-1]-(s-1);
//   for(auto i=s; i <= s+K-2; i++)
//     q[i] = h[i-s+1]+(s-1);
//   for(auto i=0; i < s; i++)
//     q[i + K] = h[i];
//   return q;
// };
// センサー値の微分を計算する方法として使えないか9月25日(土)
inline V_d UniformKnots(V_d h, const int K) {
  double scale = *h.rbegin() - *h.begin();
  *h.begin() -= scale * 1E-14;
  *h.rbegin() += scale * 1E-14;
  return Subdivide(*h.begin(), *h.rbegin(), h.size() + K - 1);
};
inline V_d UniformKnots(const int size, const int K) { return UniformKnots(Subdivide(double(-1), double(1), size - 1), K); };
inline V_d OpenUniformKnots(V_d h, const int K) {
  int s(h.size());
  V_d q(s + K);
  const double EPS = 1E-14;
  *h.begin() -= EPS;
  *h.rbegin() += EPS;
  for (auto i = 0; i < K; ++i) {
    // q[i] = h[0] - EPS * (i - K + 1.);
    q[i] = std::fma(-EPS, static_cast<double>(i - K + 1.), h[0]);
    // q[i + s] = h[s - 1] + EPS * i;
    q[i + s] = std::fma(EPS, static_cast<double>(i), h[s - 1]);
  }
  for (auto i = 0; i < s - K; i++)
    q[i + K] = (h[i] + h[(i + K)]) * 0.5;

  return q;
};

inline V_d OpenUniformKnots(const int size, const int K) { return OpenUniformKnots(Subdivide(double(-1), double(1), size - 1), K); };

inline V_d PeriodicUniformKnots(const V_d &h, const int K) {
  // 2021/09/25
  int s(h.size());
  V_d q(s + K);
  for (auto i = 0; i < s; i++)
    q[i] = h[i];
  for (auto i = -2 * K; i < 0; i++)
    q[i] = h[i + s - 1] - (s - 1);
  for (auto i = s; i < s + K; i++)
    q[i] = h[i - s + 1] + (s - 1);
  return q;
};

inline V_d Bspline_knot(const int size, const int K) {
  V_d h = Subdivide(double(-1), double(1), size - 1), q(size + K);
  double EPS = 1.E-15;
  for (auto i = 0; i < K; i++)
    q[i] = h[0] - EPS * (i - K + 1.);
  for (auto i = 0; i < size - K; i++)
    q[i + K] = (h[i] + h[(i + K)]) / 2.;
  for (auto i = size; i < size + K; i++)
    q[i] = h[size - 1] + EPS * (i - size);
  // for (auto i = 0; i < K; i++)
  // {
  //   q[i] = h[0] - EPS * (i - K + 1.);
  //   q[i + size] = h[size - 1] + EPS * i;
  // }
  // for (auto i = 0; i < size - K; i++)
  //   q[i + K] = (h[i] + h[i + K]) / 2.;
  return q;
};

inline V_d Bspline_knot_periodic(const int size, const int K) {
  V_d h = Subdivide(double(-1), double(1), size - 1), q(size + 4 * K - 2);
  for (auto i = 0; i < size; i++)
    q[i] = h[i];
  for (auto i = -2 * K; i < 0; i++)
    q[i] = h[i + size - 1] - (size - 1);
  for (auto i = size; i < size + 2 * K - 1; i++)
    q[i] = h[i - size + 1] + (size - 1);
  return q;
};

inline V_d Bspline_knot(const V_d &h, const int K) {
  int s(h.size());
  V_d q(s + K);
  double EPS = 1E-15;
  for (auto i = 0; i < K; i++)
    q[i] = *h.begin() - EPS * (i - K + 1.);
  for (auto i = 0; i < s - K; i++)
    q[i + K] = (h[i] + h[(i + K)]) / 2.;
  for (auto i = s; i < s + K; i++)
    q[i] = *h.rbegin() + EPS * (i - s);
  return q;
};
inline VV_d Bspline_knot(const VV_d &h, const int K) {
  VV_d q(h.size());
  for (size_t i = 0; i < h.size(); i++)
    q[i] = Bspline_knot(h[i], K);
  return q;
};
inline V_d Bspline_vector(const double h, const int s, const int K) {
  //  --------------  s --------------->
  // {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)}
  V_d ret(s, 0.), q = OpenUniformKnots(s, K);
  for (auto j = 0; j < s; j++)
    ret[j] = Bspline(h, q, j, K);
  return ret;
}; // parametric -1 to 1
inline V_d Bspline_vector(const double h, const int s, const int K, const V_d &q) {
  V_d ret(s, 0.);
  for (auto j = 0; j < s; j++)
    ret[j] = Bspline(h, q, j, K);
  return ret;
};
inline V_d D_Bspline_vector(const double h, const int s, const int K, int n) {
  //  --------------  s --------------->
  // {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)}
  V_d ret(s, 0.), q = OpenUniformKnots(s, K);
  for (auto j = 0; j < s; j++)
    ret[j] = D_Bspline(h, q, j, K, n);
  return ret;
};
inline V_d D_Bspline_vector(const double h, const int s, const int K, int n, const V_d &q) {
  V_d ret(s, 0.);
  for (auto j = 0; j < s; j++)
    ret[j] = D_Bspline(h, q, j, K, n);
  return ret;
};
//           --------------  s  --------------->
// |        {
// |        {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)},
// |        {B0(h1),B1(h1),B2(h1),B3(h1),B4(h1)},
// |        {B0(h2),B1(h2),B2(h2),B3(h2),B4(h2)},
// h.size() {B0(h3),B1(h3),B2(h3),B3(h3),B4(h3)},
// |        {B0(h4),B1(h4),B2(h4),B3(h4),B4(h4)},
// |        {B0(h5),B1(h5),B2(h5),B3(h5),B4(h5)},
// |        {B0(h6),B1(h6),B2(h6),B3(h6),B4(h6)},
// V        }
inline VV_d Bspline_matrix(const V_d &h, const int s, const int K) {
  VV_d ret(h.size(), V_d(s, 0.));
  for (size_t i = 0; i < h.size(); i++)
    ret[i] = Bspline_vector(h[i], s, K);
  return ret;
};
inline VV_d Bspline_matrix(const V_d &h, const int s, const int K, const V_d &q) {
  VV_d ret(h.size(), V_d(s, 0.));
  for (size_t i = 0; i < h.size(); i++)
    ret[i] = Bspline_vector(h[i], s, K, q);
  return ret;
};
inline VV_d Bspline_matrix(const VV_d &h, const int s, const int K, const VV_d &q) {
  VV_d ret(h.size(), V_d(s, 0.));
  for (size_t i = 0; i < h.size(); i++)
    for (auto j = 0; j < s; j++)
      ret[i][j] = Bspline(h[i][j], q[i], s, K);
  return ret;
};
inline VV_d Bspline_matrix(const VV_d &h, const int s, const int K, const V_d &q) {
  VV_d ret(h.size(), V_d(s, 0.));
  for (size_t i = 0; i < h.size(); i++)
    for (auto j = 0; j < s; j++)
      ret[i][j] = Bspline(h[i][j], q, i, K);
  return ret;
};
inline VV_d D_Bspline_matrix(const V_d &h, const int s, const int K, const int n) {
  VV_d ret(h.size(), V_d(s, 0.));
  for (size_t i = 0; i < h.size(); i++)
    ret[i] = D_Bspline_vector(h[i], s, K, n);
  return ret;
};
inline VV_d D_Bspline_matrix(const V_d &h, const int s, const int K, const int n, const V_d &q) {
  VV_d ret(h.size(), V_d(s, 0.));
  for (size_t i = 0; i < h.size(); i++)
    ret[i] = D_Bspline_vector(h[i], s, K, n, q);
  return ret;
};
//    ---------------  s --------------->
// | {
// | {B0(h0),B1(h0),B2(h0),B3(h0),B4(h0)},
// | {B0(h1),B1(h1),B2(h1),B3(h1),B4(h1)},
// s {B0(h2),B1(h2),B2(h2),B3(h2),B4(h2)},
// | {B0(h3),B1(h3),B2(h3),B3(h3),B4(h3)},
// | {B0(h4),B1(h4),B2(h4),B3(h4),B4(h4)}
// V }
inline VV_d Bspline_matrix(const int s, const int K) { return Bspline_matrix(Subdivide(-1., 1., s - 1) /*h*/, s, K); };
inline VV_d D_Bspline_matrix(const int s, const int K, const int n) { return D_Bspline_matrix(Subdivide(-1., 1., s - 1) /*h*/, s, K, n); };
///////////////////////////////////////////////////////////
struct ParametricInterpolation {
public:
  VV_d samp3;
  V_d samp2;
  V_d q, xi, eta;
  V_d q1, q2;
  VV_d InvMatB, InvMatBT;
  int s, K;
  int s1, s2;
  ///////////// 1 parm  ///////////
  ParametricInterpolation(const V_d &samp2_IN, const int K_IN) : samp2(samp2_IN), s(samp2_IN.size()), xi(Subdivide(double(-1), double(1), samp2_IN.size() - 1)), q(OpenUniformKnots(Subdivide(double(-1), double(1), samp2_IN.size() - 1), K_IN)), K(K_IN), InvMatBT(Inverse(Transpose(Bspline_matrix(samp2_IN.size(), K_IN)))), InvMatB(Inverse(Bspline_matrix(samp2_IN.size(), K_IN))) {};
  V_d N(const double a) { return Dot(InvMatBT, Bspline_vector(a, s, K)); };
  double operator()(const double a) { return Dot(N(a), samp2); };
  //===========================
  V_d N(const double a, const int n) { return Dot(InvMatBT, D_Bspline_vector(a, s, K, n)); };
  double operator()(const double a, const int n) { return Dot(Dot(InvMatB, D_Bspline_vector(a, s, K, n)), samp2); };
  //========= specific ========
  double Basis(const double h, const int i) { return Bspline(h, q, i, K); };
  V_d Basis(const V_d &h, const int i) { return Bspline(h, q, i, K); };
  double Shape(const double h, const int i) { return Dot(InvMatBT[i], Bspline_vector(h, s, K)); };
  double DShape(const double h, const int i, const int n) { return Dot(InvMatBT[i], D_Bspline_vector(h, s, K, n)); };
  /////////// 2 parm  ////////////
  /*
  s1 is the size of a column vector : points in y direction
  s2 is the size of a row vector : points in x direction
*/
  ParametricInterpolation(const VV_d &samp3_IN, const int K_IN) : samp3(samp3_IN), s1(samp3_IN.size()) /*y*/, s2(samp3_IN[0].size()) /*x*/, q1(OpenUniformKnots(samp3_IN.size(), K_IN)), q2(OpenUniformKnots(samp3_IN[0].size(), K_IN)), K(K_IN), InvMatB(Inverse(Bspline_matrix(samp3_IN.size() /*y*/, K_IN))), InvMatBT(Inverse(Transpose(Bspline_matrix(samp3_IN[0].size() /*x*/, K_IN)))) {};
  ParametricInterpolation() {};
  void reset(const VV_d &samp3_IN, const int K_IN) {
    q1.clear();
    q2.clear();
    samp3.clear();
    InvMatB.clear();
    InvMatBT.clear();
    // insert
    samp3 = samp3_IN;
    s1 = samp3_IN.size() /*y*/;
    s2 = samp3_IN[0].size() /*x*/;
    q1 = OpenUniformKnots(samp3_IN.size(), K_IN);
    q2 = OpenUniformKnots(samp3_IN[0].size(), K_IN);
    K = K_IN;
    InvMatB = Inverse(Bspline_matrix(samp3_IN.size(), K_IN));
    InvMatBT = Inverse(Transpose(Bspline_matrix(samp3_IN[0].size(), K_IN)));
  };
  ////////////////////////////////////////////////////
  V_d Nx(const double b /*y*/) { return Dot(InvMatBT, Bspline_vector(b, s2, K, q2)); };
  V_d Ny(const double a /*x*/) { return Dot(Bspline_vector(a, s1, K, q1), InvMatB); };
  VV_d N(const V_d &ba) { return TensorProduct(Ny(ba[1]), Nx(ba[0])); };
  double operator()(const double b /*x*/, const double a /*y*/) { return Dot(Ny(a), Dot(samp3, Nx(b))); };
  double operator()(const V_d &ba) { return Dot(Ny(ba[1]), Dot(samp3, Nx(ba[0]))); };
  ////////////////// DERIVATIVES /////////////////
  V_d DNx(const double b /*y*/, const int n) { return Dot(InvMatBT, D_Bspline_vector(b, s2, K, n, q2)); };
  V_d DNy(const double a /*x*/, const int n) { return Dot(D_Bspline_vector(a, s1, K, n, q1), InvMatB); };
  VV_d DN(const V_d &ba, const V_i &n) { return TensorProduct(DNy(ba[1], n[1]), DNx(ba[0], n[0])); };
  double operator()(const double b /*x*/, const double a /*y*/, const V_i &n) { return Dot(DNy(a, n[1]), Dot(samp3, DNx(b, n[0]))); };
  double operator()(const V_d &ba, const V_i &n) { return Dot(DNy(ba[1], n[1]), Dot(samp3, DNx(ba[0], n[0]))); };
  ///////////////// specific shape function ////////////
  double Nx(const double b /*y*/, const int j) { return Dot(InvMatBT[j], Bspline_vector(b, s2, K, q2)); };
  double Ny(const double a /*x*/, const int m) {
    V_d By = Bspline_vector(a, s1, K, q1);
    double ret(0.);
    for (auto l = 0; l < s1; l++)
      ret += InvMatB[l][m] * By[l];
    return ret;
  };
  double N(const V_d &ba, const V_i &jm) { /*2D vector access*/ return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]); };
  double N(const V_d &ba, const int mj) { /*1D access*/ return Ny(ba[1], (int)(mj / s2) /*row*/) * Nx(ba[0], mj % s2 /*col*/); };
  // double Basis(const V_d& ba, const V_i &jm){
  //   return samp3[jm[1]][jm[0]] * Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]);
  // };
  double Basis(const V_d &ba, const V_i &jm) { return Bspline(ba[0], q2, jm[0], K) * Bspline(ba[1], q1, jm[1], K); };
  double Shape(const V_d &ba, const V_i &jm) { return Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]); };
};
///////////////////////////////////////////////////////////
class node_complex {
  //    *
  //    |
  // *--*--*
  //    |
  //    *
public:
  VV_d v_complex;
  std::vector<node_complex *> nodes;
  std::vector<bool> hits;
  node_complex(const VV_d &v_complex_IN) : v_complex(v_complex_IN) {};

  void connect(node_complex &node_IN) {
    node_complex *p = std::addressof(node_IN);
    for (auto &n : nodes)
      if (p == n)
        return;
    nodes.emplace_back(p);
  };
};

///////////////////////////////////////////////////////////
class glLINES {
public:
  VV_d v_complex;
  std::vector<std::vector<float>> s_complex;
  VV_i f_v_complex; // これのstd::vector<int>が１面 -> 3頂点を指定
  std::vector<bool> hits;

  VV_d samp3X, samp3Y, samp3Z;
  // sample点がない場合
  glLINES(int row = 20, int col = 20, double scale = 5.) {
    /* sampleを作成 */
    VV_d samp3X(row, V_d(col, 0));
    VV_d samp3Y(row, V_d(col, 0));
    VV_d samp3Z(row, V_d(col, 0));
    for (auto i = 0; i < row; i++) {
      for (auto j = 0; j < col; j++) {
        double x = (-1. + 2. / (col - 1) * j);
        double y = (-1. + 2. / (row - 1) * i);
        samp3X[i][j] = x * scale;
        samp3Y[i][j] = y * scale;
        samp3Z[i][j] = std::sin(2. * M_PI * x) * std::sin(2. * M_PI * y);
      }
    }
    /* gl.LINESでメッシュがかけるような形式で格納 */
    set(samp3X, samp3Y, samp3Z);
  };
  // sample点がある場合
  glLINES(const VV_d &samp3X_IN, const VV_d &samp3Y_IN, const VV_d &samp3Z_IN) : samp3X(samp3X_IN), samp3Y(samp3Y_IN), samp3Z(samp3Z_IN) {
    /* gl.LINESでメッシュがかけるような形式で格納 */
    set(samp3X_IN, samp3Y_IN, samp3Z_IN);
  };
  // verticesの格納
  void set(const VV_d &samp3X_IN, const VV_d &samp3Y_IN, const VV_d &samp3Z_IN) {
    auto col = samp3X_IN[0].size();
    auto row = samp3X_IN.size();
    v_complex.clear();
    s_complex.clear();
    f_v_complex.clear();
    /* gl.LINESでメッシュがかけるような形式で格納 */
    int l = 0;
    //    8 points
    //    0   1
    //   7*---*2
    //    |   |
    //    |   |
    //   6*---*3
    //    5   4
    for (size_t i = 0; i < row - 1; i++)
      for (size_t j = 0; j < col - 1; j++) {
        v_complex.push_back({samp3X_IN[i][j], samp3Z_IN[i][j], samp3Y_IN[i][j]});             // 0
        v_complex.push_back({samp3X_IN[i][j + 1], samp3Z_IN[i][j + 1], samp3Y_IN[i][j + 1]}); // 1
        f_v_complex.push_back({l++, l++});
        v_complex.push_back({samp3X_IN[i][j + 1], samp3Z_IN[i][j + 1], samp3Y_IN[i][j + 1]});             // 2
        v_complex.push_back({samp3X_IN[i + 1][j + 1], samp3Z_IN[i + 1][j + 1], samp3Y_IN[i + 1][j + 1]}); // 3
        f_v_complex.push_back({l++, l++});
        v_complex.push_back({samp3X_IN[i + 1][j + 1], samp3Z_IN[i + 1][j + 1], samp3Y_IN[i + 1][j + 1]}); // 4
        v_complex.push_back({samp3X_IN[i + 1][j], samp3Z_IN[i + 1][j], samp3Y_IN[i + 1][j]});             // 5
        f_v_complex.push_back({l++, l++});
        v_complex.push_back({samp3X_IN[i + 1][j], samp3Z_IN[i + 1][j], samp3Y_IN[i + 1][j]}); // 6
        v_complex.push_back({samp3X_IN[i][j], samp3Z_IN[i][j], samp3Y_IN[i][j]});             // 7
        f_v_complex.push_back({l++, l++});
      }

    s_complex.resize(v_complex.size(), std::vector<float>(4, 0)); // RGBS
    for (auto &tmp : s_complex)
      tmp = {.9, .0, .0, 1.};

    this->hits.resize(f_v_complex.size());
  };
};
///////////////////////////////////////////////////////////
struct ParametricInterpolation3D {
public:
  VV_d v_complex;
  std::vector<std::vector<float>> s_complex;
  VV_i f_v_complex; // これのstd::vector<int>が１面 -> 3頂点を指定

  VV_d samp3X, samp3Y, samp3Z;
  V_d q1, q2;
  VV_d inv_coefmat, inv_coefmatT;
  int s1, s2, K;
  ////////// 2 parm  //////////
  /*
  s1 is the size of a column vector : points in y direction
  s2 is the size of a row vector : points in x direction
*/
  ParametricInterpolation3D(const VV_d &samp3X_IN, const VV_d &samp3Y_IN, const VV_d &samp3Z_IN, const int K_IN) : K(K_IN), s1(samp3X_IN.size()) /*y*/, s2(samp3X_IN[0].size()) /*x*/, q1(OpenUniformKnots(samp3X_IN.size(), K_IN)), q2(OpenUniformKnots(samp3X_IN[0].size(), K_IN)), inv_coefmat(Inverse(Bspline_matrix(samp3X_IN.size() /*y*/, K_IN))), inv_coefmatT(Inverse(Transpose(Bspline_matrix(samp3X_IN[0].size() /*x*/, K_IN)))), samp3X(samp3X_IN), samp3Y(samp3Y_IN), samp3Z(samp3Z_IN) {
    if (samp3X_IN.size() != samp3Y_IN.size() || samp3X_IN.size() != samp3Z_IN.size()) {
      Print("sizes of input samples are differen", Red);
      Print(samp3X_IN.size(), Red);
      Print(samp3Y_IN.size(), Red);
      Print(samp3Z_IN.size(), Red);
      abort();
    }
  };
  ParametricInterpolation3D() {};
  void reset(const VV_d &samp3X_IN, const VV_d &samp3Y_IN, const VV_d &samp3Z_IN, const int K_IN) {
    if (samp3X_IN.size() != samp3Y_IN.size() || samp3X_IN.size() != samp3Z_IN.size()) {
      Print("sizes of input samples are differen", Red);
      Print(samp3X_IN.size(), Red);
      Print(samp3Y_IN.size(), Red);
      Print(samp3Z_IN.size(), Red);
      abort();
    }
    q1.clear();
    q2.clear();
    inv_coefmat.clear();
    inv_coefmatT.clear();
    samp3X.clear();
    samp3Y.clear();
    samp3Z.clear();
    // insert
    K = K_IN;
    s1 = samp3X_IN.size() /*y*/;
    s2 = samp3X_IN[0].size() /*x*/;
    q1 = OpenUniformKnots(samp3X_IN.size(), K_IN);
    q2 = OpenUniformKnots(samp3X_IN[0].size(), K_IN);
    inv_coefmat = Inverse(Bspline_matrix(samp3X_IN.size(), K_IN));
    inv_coefmatT = Inverse(Transpose(Bspline_matrix(samp3X_IN[0].size(), K_IN)));
    samp3X = samp3X_IN;
    samp3Y = samp3Y_IN;
    samp3Z = samp3Z_IN;
  };

  ////////////////////////////////////////////////////
  V_d Nx(const double b /*y*/) { return Dot(inv_coefmatT, Bspline_vector(b, s2, K, q2)); };
  V_d Ny(const double a /*x*/) { return Dot(Bspline_vector(a, s1, K, q1), inv_coefmat); };
  VV_d N(const V_d &ba) { return TensorProduct(Ny(ba[1]), Nx(ba[0])); };
  V_d operator()(const double b /*x*/, const double a /*y*/) {
    V_d shapey = Ny(a), shapex = Nx(b), ret(3, 0.);
    for (auto m = 0; m < s1; m++) {
      ret[0] += shapey[m] * std::inner_product(samp3X[m].begin(), samp3X[m].end(), shapex.begin(), 0.);
      ret[1] += shapey[m] * std::inner_product(samp3Y[m].begin(), samp3Y[m].end(), shapex.begin(), 0.);
      ret[2] += shapey[m] * std::inner_product(samp3Z[m].begin(), samp3Z[m].end(), shapex.begin(), 0.);
    }
    return ret;
  };
  V_d operator()(const V_d &ba) {
    V_d shapey = Ny(ba[1]), shapex = Nx(ba[0]), ret(3, 0.);
    for (auto m = 0; m < s1; m++) {
      ret[0] += shapey[m] * std::inner_product(samp3X[m].begin(), samp3X[m].end(), shapex.begin(), 0.);
      ret[1] += shapey[m] * std::inner_product(samp3Y[m].begin(), samp3Y[m].end(), shapex.begin(), 0.);
      ret[2] += shapey[m] * std::inner_product(samp3Z[m].begin(), samp3Z[m].end(), shapex.begin(), 0.);
    }
    return ret;
  };
  ////////////////// DERIVATIVES /////////////////
  V_d DNx(const double b /*y*/, const int n) { return Dot(inv_coefmatT, D_Bspline_vector(b, s2, K, n, q2)); };
  V_d DNy(const double a /*x*/, const int n) {
    V_d By = D_Bspline_vector(a, s1, K, n, q1);
    V_d ret(s1, 0.);
    for (auto m = 0; m < s1; m++)
      for (auto l = 0; l < s1; l++)
        ret[m] += inv_coefmat[l][m] * By[l];
    return ret;
  };
  VV_d DN(const V_d &ba, const V_i &n) {
    V_d shapey = DNy(ba[1], n[1]), shapex = DNx(ba[0], n[0]);
    VV_d ret(s1, V_d(s2, 0.));
    for (auto m = 0; m < s1; m++)
      for (auto j = 0; j < s2; j++)
        ret[m][j] = shapey[m] * shapex[j];
    return ret;
  };
  V_d operator()(const double b /*x*/, const double a /*y*/, const V_i &n) {
    V_d shapey = DNy(a, n[1]), shapex = DNx(b, n[0]), ret(3, 0.);
    for (auto m = 0; m < s1; m++)
      for (auto j = 0; j < s2; j++) {
        ret[0] += samp3X[m][j] * shapey[m] * shapex[j];
        ret[1] += samp3Y[m][j] * shapey[m] * shapex[j];
        ret[2] += samp3Z[m][j] * shapey[m] * shapex[j];
      }
    return ret;
  };
  V_d operator()(const V_d &ba, const std::vector<int> &n) {
    V_d shapey = DNy(ba[1], n[1]), shapex = DNx(ba[0], n[0]), ret(3, 0.);
    for (auto m = 0; m < s1; m++)
      for (auto j = 0; j < s2; j++) {
        ret[0] += samp3X[m][j] * shapey[m] * shapex[j];
        ret[1] += samp3Y[m][j] * shapey[m] * shapex[j];
        ret[2] += samp3Z[m][j] * shapey[m] * shapex[j];
      }
    return ret;
  };
  ///////////////// specific shape function ////////////
  double Nx(const double b /*y*/, const int j) { return std::inner_product(inv_coefmatT[j].begin(), inv_coefmatT[j].end(), Bspline_vector(b, s2, K, q2).begin(), 0.); };
  double Ny(const double a /*x*/, const int m) {
    V_d By = Bspline_vector(a, s1, K, q1);
    double ret(0.);
    for (auto l = 0; l < s1; l++)
      ret += inv_coefmat[l][m] * By[l];
    return ret;
  };
  double N(const V_d &ba, const std::vector<int> &jm) { /*2D vector access*/ return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]); };
  double N(const V_d &ba, const int mj) { /*1D access*/ return Ny(ba[1], (int)(mj / s2) /*row*/) * Nx(ba[0], mj % s2 /*col*/); };
  V_d basis(const V_d &ba, const std::vector<int> &jm) { return V_d{samp3X[jm[1]][jm[0]], samp3Y[jm[1]][jm[0]], samp3Z[jm[1]][jm[0]]} * Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]); };
  V_d cross(const V_d &ba) {
    V_d Dshapey = DNy(ba[1], 1), Dshapex = DNx(ba[0], 1), shapey = DNy(ba[1], 0), shapex = DNx(ba[0], 0), ret_dxi(3, 0.), ret_deta(3, 0.);
    for (auto m = 0; m < s1; m++)
      for (auto j = 0; j < s2; j++) {
        ret_dxi[0] += samp3X[m][j] * shapey[m] * Dshapex[j];
        ret_dxi[1] += samp3Y[m][j] * shapey[m] * Dshapex[j];
        ret_dxi[2] += samp3Z[m][j] * shapey[m] * Dshapex[j];
        ret_deta[0] += samp3X[m][j] * Dshapey[m] * shapex[j];
        ret_deta[1] += samp3Y[m][j] * Dshapey[m] * shapex[j];
        ret_deta[2] += samp3Z[m][j] * Dshapey[m] * shapex[j];
      }
    return Cross(ret_dxi, ret_deta);
  };
};
////////////////////////////////////////////////////////
struct RBF_fn {
  virtual double rbf(double r) = 0;
};
struct RBF_interp {
  int dim, n;
  VV_d pts;
  V_d vals;
  V_d w;
  RBF_fn &fn;
  bool norm;

  RBF_interp(const VV_d &ptss /* x for f(x) */, const V_d &valss /*f(x)*/, RBF_fn &func, bool nrbf = false) : dim(ptss[0].size()), n(ptss.size()), pts(ptss), vals(valss), w(n), fn(func), norm(nrbf) {
    double sum;
    VV_d rbf(n, V_d(n, 0.));
    V_d rhs(n);
    for (auto i = 0; i < n; i++) {
      sum = 0.;
      for (auto j = 0; j < n; j++)
        sum += (rbf[i][j] = fn.rbf(rad(&pts[i][0], &pts[j][0])));
      rhs[i] = norm ? sum * vals[i] : vals[i];
    }
    ludcmp lu(rbf);
    lu.solve(rhs, w);
  }
  double operator()(const V_d &pt) {
    double fval, sum = 0., sumw = 0.;
    if (pt.size() != dim)
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    for (auto i = 0; i < n; i++) {
      fval = fn.rbf(rad(&pt[0], &pts[i][0]));
      sumw += w[i] * fval;
      sum += fval;
    }
    return norm ? sumw / sum : sumw;
  }

  double rad(const double *p1, const double *p2) {
    double sum = 0.;
    for (int i = 0; i < dim; i++)
      sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    return std::sqrt(sum);
  }
};

struct RBF_multiquadric : RBF_fn {
  double r02;
  RBF_multiquadric(double scale = 1.) : r02(scale * scale) {}
  double rbf(double r) { return std::sqrt(r * r + r02); }
};
struct RBF_thinplate : RBF_fn {
  double r0;
  RBF_thinplate(double scale = 1.) : r0(scale) {}
  double rbf(double r) { return r <= 0. ? 0. : r * r * std::log(r / r0); }
};
struct RBF_gauss : RBF_fn {
  double r0;
  RBF_gauss(double scale = 1.) : r0(scale) {}
  double rbf(double r) { return exp(-0.5 * r * r / (r0 * r0)); }
};
struct RBF_inversemultiquadric : RBF_fn {
  double r02;
  RBF_inversemultiquadric(double scale = 1.) : r02(scale * scale) {}
  double rbf(double r) { return 1. / std::sqrt(r * r + r02); }
};
struct Shep_interp {
  int dim, n;
  const VV_d &pts;
  const V_d &vals;
  double pneg;

  Shep_interp(VV_d &ptss, V_d &valss, double p = 2.) : dim(ptss[0].size()), n(ptss.size()), pts(ptss), vals(valss), pneg(-p) {}

  double interp(V_d &pt) {
    double r, w, sum = 0., sumw = 0.;
    if (pt.size() != dim)
      throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ""));
    for (auto i = 0; i < n; i++) {
      if ((r = rad(&pt[0], &pts[i][0])) == 0.)
        return vals[i];
      sum += (w = std::pow(r, pneg));
      sumw += w * vals[i];
    }
    return sumw / sum;
  }

  double rad(const double *p1, const double *p2) {
    double sum = 0.;
    for (auto i = 0; i < dim; i++)
      sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    return std::sqrt(sum);
  }
};

// #ifdef INCL_lapacke
// struct LU_LAPACK{
// public:
//   char TRANS = 'N';
//   int INFO;
//   int LDA;
//   int LDB;
//   int N;
//   int NRHS = 1;
//   std::vector<int> IPIV;
//   V_d Flatten_A;

//   LU_LAPACK(const std::vector<V_d>& A/* must be square matrix */):
//     Flatten_A(Flatten(A)),
//     INFO(A.size()),
//     LDA(A.size()),
//     LDB(A.size()),
//     N(A.size()),
//     IPIV(A.size())
//   {
//     std::cout << "compute the LU factorization..." << std::endl;
//     //void LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda, lapack_int* ipiv, lapack_int *info );
//     LAPACK_dgetrf(&N, &N, &Flatten_A[0], &LDA, &IPIV[0], &INFO);
//     if(INFO)
//       std::cout << "an error occured : "<< INFO << std::endl;
//   };
//   /////////////////////////////////////////////////////////////
//   V_d solve(const V_d& B){
//     V_d ret(B);
//     // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
//     std::cout << "solving the system..."<< std::endl;
//     // void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs, const double* a, lapack_int* lda, const lapack_int* ipiv,double* b, lapack_int* ldb, lapack_int *info );
//     dgetrs_(&TRANS, &N, &NRHS, &Flatten_A[0], &LDA, &IPIV[0], &ret[0], &LDB,&INFO);
//     if(INFO)
//       {
// 	std::cout << "an error occured : "<< INFO << std::endl;
// 	abort();
//       }
//     return ret;
//   };
//   /////////////////////////////////////////////////////////////
//   void solve(const V_d& B, V_d& ret){
//     // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
//     std::cout << "solving the system..."<< std::endl;
//     // void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs, const double* a, lapack_int* lda, const lapack_int* ipiv,double* b, lapack_int* ldb, lapack_int *info );
//     dgetrs_(&TRANS, &N, &NRHS, &Flatten_A[0], &LDA,&IPIV[0], &ret[0], &LDB,&INFO);
//     if(INFO)
//       {
// 	std::cout << "an error occured : "<< INFO << std::endl;
// 	abort();
//       }
//   };
// };
// #endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//=========================================================
// void Load(const std::string& filename, const V_s& sep, std::vector<std::vector<V_d>>& mat){
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     std::cout << Red << "The file can not be opened" << reset << std::endl;
//   std::string str, S, cell;
//   double v;
//   std::size_t found, s(sep.size()), min_len, left, right;
//   bool loop;
//   int row_counter(0), col_counter(0);
//   while(!strm.eof( ) && mat.size() > row_counter)
//     {
//       std::getline(strm,str);
//       if(!str.empty())
// 	{
// 	  if((left = str.find("{"))!=std::string::npos);
// 	  if((right = str.find("}"))!=std::string::npos);
// 	  cell = str.substr(left,right);
// 	  str = str.substr(right+1);

// 	  V_d row_vector;
// 	  loop = true;
// 	  col_counter = 0;

// 	  while(loop && mat[row_counter].size() > col_counter)
// 	    {
// 	      min_len = 999999;
// 	      for(auto i=0; i<s; i++)
// 		if((found = cell.find(sep[i]))!=std::string::npos)
// 		  if(min_len > found)
// 		    {
// 		      min_len = found;
// 		      S = sep[i];
// 		    }
// 	      col_counter++;
// 	    }
// 	  row[row_counter].emplace_back(row_vector);
// 	  row_counter++;
// 	}
//     }
//   strm.close();
// };
//=========================================================
template <class Ttype> void load(std::ifstream &in, std::vector<Ttype> &vec, const int row_size) {
  vec.clear();
  if (!in) {
    std::cout << "入力ファイルが開けません" << std::endl;
  }
  std::string str;
  int i;
  Ttype tmp;
  i = 0;
  getline(in, str);
  std::cout << std::string(40, '=') << std::endl;
  std::cout << "この行は，" << str << std::endl;
  std::cout << std::string(40, '-') << std::endl;
  std::cout << "読み込み結果：" << std::endl;
  for (size_t j = 0; j < row_size; j++) {
    if (str.find(',') != std::string::npos) {
      i = str.find(',');
      //	std::cout << "\',\'は" << i << "番目にあります" << endl;
      //	std::cout << str << "　を" << i << "まで読み込みます -> " << str.substr(0,i) << " -> これを数値に変換すると ->";
      tmp = (Ttype)stod(str.substr(0, i));
      //	std::cout << std::setprecision(10) << tmp << endl;
      str = str.substr(i + 1);
      //	std::cout << i+1 << "番目以降を抜き出し文字列を置き換えます -> " << str << endl;
    } else {
      //	std::cout << "\',\'は" << "見つかりませんでした -> ";
      tmp = (Ttype)stod(str);
      //	std::cout << str << "　を" << "最後まで読み込みます -> " << std::setprecision(10) << tmp << endl;
    }
    vec.emplace_back(tmp);
    //      std::cout << tmp << "を代入 -> ";
    std::cout << "vec[" << j << "] = " << std::setprecision(5) << vec[j] << std::endl;
  }
  std::cout << std::string(40, '-') << std::endl;
};
//===========================================================
template <class Ttype> void load(std::ifstream &in, std::vector<std::vector<Ttype>> &vec, const int row_size) {
  vec.clear();
  vec.resize(row_size, std::vector<Ttype>(0));
  if (!in) {
    std::cout << "入力ファイルが開けません" << std::endl;
  }
  std::string str;
  int i;
  Ttype tmp;
  for (int line = 0;; line++) {
    i = 0;
    getline(in, str);
    if (in.eof())
      break;
    std::cout << std::string(40, '=') << std::endl;
    std::cout << "この行は，" << str << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    std::cout << "読み込み結果：" << std::endl;
    for (int j = 0; j < row_size; j++) {
      if (str.find(',') != std::string::npos) {
        i = str.find(',');
        //	std::cout << "\',\'は" << i << "番目にあります" << endl;
        //	std::cout << str << "　を" << i << "まで読み込みます -> " << str.substr(0,i) << " -> これを数値に変換すると ->";
        tmp = (Ttype)stod(str.substr(0, i));
        //	std::cout << std::setprecision(10) << tmp << endl;
        str = str.substr(i + 1);
        //	std::cout << i+1 << "番目以降を抜き出し文字列を置き換えます -> " << str << endl;
      } else {
        //	std::cout << "\',\'は" << "見つかりませんでした -> ";
        tmp = (Ttype)stod(str);
        //	std::cout << str << "　を" << "最後まで読み込みます -> " << std::setprecision(10) << tmp << endl;
      }
      vec[j].emplace_back(tmp);
      //      std::cout << tmp << "を代入 -> ";
      std::cout << "vec[" << j << "][" << line << "] = " << std::setprecision(5) << vec[j][line] << std::endl;
    }
    std::cout << std::string(40, '-') << std::endl;
  };
};
//=====================================================================
//=================== Mathematica output loader =======================
//=====================================================================

inline V_s StringSplit(const std::string &strIN, const V_s &SEP) {
  std::string str(strIN), foundFirstSep;
  V_s ret, tmp;
  size_t foundFirst;
  bool bfound;

  while (!str.empty()) {
    foundFirst = std::string::npos;
    bfound = false;

    // Find the closest separator
    for (const auto &sep : SEP) {
      size_t found = str.find(sep);
      if (!sep.empty() && found != std::string::npos && (found < foundFirst)) {
        foundFirst = found;
        foundFirstSep = sep;
        bfound = true;
      }
    }

    if (bfound) {
      tmp.emplace_back(str.substr(0, foundFirst));
      str = str.substr(foundFirst + foundFirstSep.length());
    } else {
      tmp.emplace_back(str);
      break;
    }
  }

  // Remove first and last cells if they are empty strings or contain only whitespace characters
  auto it = tmp.begin();
  while (it != tmp.end()) {
    if (it->empty() || it->find_first_not_of(" \r\n\t") == std::string::npos) {
      it = tmp.erase(it);
    } else {
      ++it;
    }
  }

  return tmp;
}

// V_s StringSplit(const std::string &strIN, const V_s &SEP) {
//    std::string str(strIN), foundFirstSep; /* expexting "1,2,3,4,5" or "1, 2,3, 4,5"*/
//    V_s ret(0), tmp(0);
//    while (!str.empty()) {
//       std::size_t foundFirst(999999);
//       bool bfound(false);
//       for (const auto &sep : SEP) {
//          std::size_t found = str.find(sep);
//          if (!sep.empty()                                       /* ignore sep, "" */
//              && found != std::string::npos                      /* if sep is found */
//              && (found < foundFirst)) /* choose closest one */  // do not use "<=" instead
//          {
//             foundFirst = found;
//             foundFirstSep = sep;
//             bfound = true;
//          }
//       }
//       if (bfound) {
//          tmp.emplace_back(str.substr(0, foundFirst));
//          str = str.substr(foundFirst + foundFirstSep.length());
//       } else {
//          tmp.emplace_back(str);

//          // remove first and last cells if they are empty strings
//          while (tmp.size() > 0)
//             if ((*tmp.begin()).empty())
//                tmp.erase(tmp.begin());
//             else if ((*tmp.rbegin()).empty() ||
//                      (*tmp.rbegin()).find_first_not_of(" ") == std::string::npos /* there is no string except spaces */ ||
//                      (*tmp.rbegin()).find_first_not_of("\r") == std::string::npos /* there is no string except \r */ ||
//                      (*tmp.rbegin()).find_first_not_of("\n") == std::string::npos /* there is no string except \n */ ||
//                      (*tmp.rbegin()).find_first_not_of("\t") == std::string::npos /* there is no string except \t */)
//                tmp.pop_back();
//             else
//                break;

//          return tmp;
//       }
//    }

//    // remove first and last cells if they are empty strings

//    while (tmp.size() > 0)
//       if ((*tmp.begin()).empty())
//          tmp.erase(tmp.begin());
//       else if ((*tmp.rbegin()).empty() ||
//                (*tmp.rbegin()).find_first_not_of(" ") == std::string::npos /* there is no string except spaces */ ||
//                (*tmp.rbegin()).find_first_not_of("\r") == std::string::npos /* there is no string except \r */ ||
//                (*tmp.rbegin()).find_first_not_of("\n") == std::string::npos /* there is no string except \n */ ||
//                (*tmp.rbegin()).find_first_not_of("\t") == std::string::npos /* there is no string except \t */)
//          tmp.pop_back();
//       else
//          break;

//    return tmp;
// };

inline std::string StringTrim(const std::string &strIN, const V_s &SEP) {
  V_s str = StringSplit(strIN, SEP);
  std::string ret;
  for (const auto &s : str)
    ret += s;
  return ret;
};
inline V_s StringTrim(const V_s &strIN, const V_s &SEP) {
  V_s strOUT(strIN.size());
  for (size_t i = 0; i < strOUT.size(); i++)
    strOUT[i] = StringTrim(strIN[i], SEP);
  return strOUT;
};
inline std::string StringJoin(const V_s &strIN, const std::string sep = "") {
  std::string ret(*strIN.begin());
  for (size_t i = 1; i < strIN.size(); i++)
    ret += sep + strIN[i];
  return ret;
};
inline std::string Directory(const std::string &f) {
  auto dir = StringSplit(f, {"/"});
  dir.pop_back();
  return "/" + StringJoin(dir, "/") + "/";
};
////////////////////////////////////////////////////////////////////////
// 2021/03/31導入
inline std::map<std::string, std::vector<std::string>> parseJSON(const std::string &str_IN) {
  enum class ParserState { Init, InKey, InValue, InArray };

  std::map<std::string, std::vector<std::string>> map_S_S;
  std::vector<std::string> SEP = {" ", "\t", "\n", "\r"};
  std::string trimmed = StringTrim(str_IN, SEP);
  std::string Lsb = "[", Rsb = "]", Col = ":", Cam = ",", Lcb = "{", Rcb = "}";
  ParserState currentState = ParserState::Init;
  std::string key = "", value = "";
  for (auto it = trimmed.begin(); it != trimmed.end(); it++) {
    switch (currentState) {
    case ParserState::Init:
      if (Lcb.find(*it) != std::string::npos)
        currentState = ParserState::InKey;
      break;
    case ParserState::InKey:
      if (Col.find(*it) != std::string::npos)
        currentState = ParserState::InValue;
      else if (*it != '"' && *it != ' ')
        key.push_back(*it);
      break;
    case ParserState::InValue:
      if (Cam.find(*it) != std::string::npos || Rcb.find(*it) != std::string::npos) {
        map_S_S[key] = StringSplit(StringTrim(value, {"\""}), {","});
        key.clear();
        value.clear();
        currentState = ParserState::InKey;
      } else if (Lsb.find(*it) != std::string::npos)
        currentState = ParserState::InArray;
      else if (*it != ' ')
        value.push_back(*it);
      break;
    case ParserState::InArray:
      if (Rsb.find(*it) != std::string::npos)
        currentState = ParserState::InValue;
      else
        value.push_back(*it);
      break;
    }
  }
  return map_S_S;
}

// std::map<std::string, std::vector<std::string>> parseJSON(const std::string &str_IN) {
//    std::map<std::string, std::vector<std::string>> map_S_S;
//    //
//    std::vector<std::string> SEP = {" ", "\t", "\n", "\r"};
//    std::string trimmed = StringTrim(str_IN, SEP);
//    std::string L = "[", R = "]", Col = ":", Cam = ",", Lb = "{", Rb = "}";
//    ParserState currentState = ParserState::Init;
//    std::string array = "", value = "";
//    std::vector<std::string> vstr(0);
//    for (auto it = trimmed.begin(); it != trimmed.end(); it++) {
//       switch (currentState) {
//          case ParserState::Init:
//             if (L.find(*it) != std::string::npos) {
//                currentState = ParserState::InArray;
//             } else if (Lb.find(*it) != std::string::npos) {
//                currentState = ParserState::InKey;
//             }
//             break;
//          case ParserState::InKey:
//             if (Col.find(*it) != std::string::npos) {
//                currentState = ParserState::InValue;
//             } else {
//                value.push_back(*it);
//             }
//             break;
//          case ParserState::InValue:
//             if (Cam.find(*it) != std::string::npos) {
//                vstr.push_back(value);
//                value.clear();
//                currentState = ParserState::InKey;
//             } else if (Rb.find(*it) != std::string::npos) {
//                vstr.push_back(value);
//                value.clear();
//             } else if (L.find(*it) != std::string::npos) {
//                currentState = ParserState::InArray;
//             } else {
//                value.push_back(*it);
//             }
//             break;
//          case ParserState::InArray:
//             if (R.find(*it) != std::string::npos) {
//                currentState = ParserState::InValue;
//                vstr.push_back(array);
//                array.clear();
//             } else {
//                array.push_back(*it);
//             }
//             break;
//       }
//    }
//    // Save the obtained std::vector<std::string> to the map.
//    for (size_t i = 0; i < vstr.size() - 1; i += 2)
//       map_S_S[StringTrim(vstr[i], {"\""})] = StringSplit(StringTrim(vstr[i + 1], {"\""}), {","});

//    return map_S_S;
// }

/* -------------------------------------------------------------------------- */

// using V_S = std::vector<std::string>;

// std::map<std::string, V_S> parseJSON(const std::string &str_IN) {
//    std::map<std::string, V_S> map_S_S;
//    //
//    V_S SEP = {" ", "\t", "\n", "\r"};
//    std::string trimed = StringTrim(str_IN, SEP);
//    std::string L = "[", R = "]", Col = ":", Cam = ",", Lb = "{", Rb = "}";
//    bool is_array = false;
//    std::string array = "", value = "";
//    V_S vstr(0);
//    for (auto it = trimed.begin(); it != trimed.end(); it++) {
//       if (is_array || L.find(*it) != std::string::npos) {
//          //[から始まるarrayの場合は，ここでキャッチする
//          if ((R.find(*it) != std::string::npos)) {
//             is_array = false;
//             vstr.emplace_back(array);
//             array.clear();
//          } else {
//             if (is_array)
//                array.push_back(*it);
//             is_array = true;
//          }
//       } else {
//          // arrayでない場合で，{},:などで始まる場合はここでキャッチする
//          if ((Lb.find(*it) != std::string::npos) ||
//              (Cam.find(*it) != std::string::npos) ||
//              (Col.find(*it) != std::string::npos) ||
//              (Rb.find(*it) != std::string::npos)) {
//             if (!value.empty()) {
//                // std::cout << red << value << std::endl;
//                vstr.emplace_back(value);
//                value.clear();
//             }
//          } else {
//             value.push_back(*it);
//          }
//       }
//    }
//    // 得られたstd::vector<std::string>をmapに保存する．
//    for (auto i = 0; i < vstr.size() - 1; i += 2)
//       map_S_S[StringTrim(vstr[i], {"\""})] = StringSplit(StringTrim(vstr[i + 1], {"\""}), {","});

//    return map_S_S;
// };

/* ------------------------------------------------------ */

/*DOC_EXTRACT basic::JSON

## C++でのJSON操作に関する実装と使用方法

### 前提条件

このドキュメントでは、`basic.hpp` と呼ばれる基本的なヘッダーファイルが含まれていると仮定します。また、サンプルのJSONデータが `./sample.json` に格納されていると仮定します。

### JSONクラス

この実装には、`JSON`という名前のC++クラスがあります。このクラスは、内部的に`std::map<std::string, std::vector<std::string>>`を持っており、JSONオブジェクトのデータを保存します。

#### コンストラクタ

このクラスにはいくつかのコンストラクタがあります:

1. ファイル名を引数として取る
   ```cpp
   JSON json("./sample.json");
   ```
2. `std::ifstream` オブジェクトを引数として取る
   ```cpp
   JSON json(std::ifstream("./sample.json"));
   ```

#### 操作

- キーを用いた値の取得
   ```cpp
   std::cout << json["translate"] << std::endl;
   ```
- キーを用いた値の設定
   ```cpp
   json["price"] = {"10."};
   ```

### 使用例

以下は、このクラスの簡単な使用例です。

```cpp
{
   std::cout << magenta << "1. ファイル名でJSONをコンストラクト" << colorReset << std::endl;
   JSON json("./sample.json");
   std::cout << json["translate"] << std::endl;
}
```

```cpp
{
   std::cout << red << "2. ifstreamでJSONをコンストラクト" << colorReset << std::endl;
   JSON json(std::ifstream("./sample.json"));
   json["price"] = {"10."};
}
```

### ファイル出力

JSONオブジェクトをファイルに出力するには、`operator<<`を使用します。

```cpp
std::ofstream os("./output.json");
os << json;
os.close();
```

### その他の機能

- `find` メソッドを使い、キーが存在するかどうかを調べることができます。
- `contains` メソッドも同様の機能を提供します。

### 注意点

1. キーが存在しない場合は、例外がスローされます。
2. 数値や真偽値は、`stod`や`stob`のような関数を使って変換する必要があります。

このドキュメントはC++でのJSON操作の基本的な概要を簡単に説明したものです。実際の使用ケースや要件に応じて、コードは適宜調整してください。

### todo

TODO: 以下のようなJSONデータを読み込めるようにする.

```
"mooring": ["simple_mooring",[10, 0.1, 0],[10, 0.1, 0]]
```

TODO: オブジェクトリテラルに対応する．

```
"name": { "asda": 1 }
```

*/

struct JSON {
  std::map<std::string, std::vector<std::string>> map_S_S;

  explicit JSON(const std::filesystem::path &path_IN) : map_S_S() {
    std::ifstream istrm(path_IN);
    if (!istrm.is_open())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not open " + std::string(path_IN));
    else {
      std::stringstream ss;
      ss << istrm.rdbuf();
      map_S_S = parseJSON(ss.str());
    }
    istrm.close();
  };

  explicit JSON(const std::string &str_IN) : map_S_S() {
    std::ifstream istrm(str_IN);
    if (!istrm.is_open())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not open " + str_IN);
    else {
      std::stringstream ss;
      ss << istrm.rdbuf();
      map_S_S = parseJSON(ss.str());
    }
    istrm.close();
  };

  explicit JSON(const std::ifstream &istrm) : map_S_S() {
    if (!istrm.is_open())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can not open");
    else {
      std::stringstream ss;
      ss << istrm.rdbuf();
      map_S_S = parseJSON(ss.str());
    }
  };

  JSON() = default;

  JSON &operator=(const JSON &other) {
    this->map_S_S = other.map_S_S;
    return *this;
  };

  std::map<std::string, std::vector<std::string>> operator()() const { return this->map_S_S; };

  std::vector<std::string> &operator[](const std::string &key) /*変更を許す*/
  {
    return this->map_S_S[key];
  };

  std::vector<std::string> at(const std::string &key) const /*変更を許さない*/
  {
    if (!this->map_S_S.contains(key)) {
      std::stringstream ss;
      ss << key << " key not found" << std::endl;
      ss << "keys: ";
      for (auto it = this->map_S_S.begin(); it != this->map_S_S.end(); it++)
        ss << it->first << " " << it->second << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    } else
      return this->map_S_S.at(key);
  };

  bool find(const std::string &key, const std::function<void(const std::vector<std::string> &)> &func) const {
    auto it = this->map_S_S.find(key);
    if (it != this->map_S_S.end()) {
      func(it->second);
      return true;
    } else
      return false;
  }
  std::vector<bool> find(const std::vector<std::string> &key) const {
    std::vector<bool> ret(key.size());
    for (auto i = 0; i < key.size(); i++)
      ret[i] = this->find(key[i]);
    return ret;
  };
  bool find(const std::string &key) const { return this->map_S_S.contains(key); };

  bool contains(const std::string &key) const { return this->map_S_S.contains(key); };

  void for_each(const std::function<void(const std::string &, const std::vector<std::string> &)> &func) const {
    for (const auto &[key, value] : this->map_S_S)
      func(key, value);
  };
};

/* ------------------------------------------------------ */
inline std::string ToString(const JSON &json) {
  auto map_S_S = json.map_S_S;
  std::stringstream stream;
  stream << "{\n";
  auto count = 0;
  for (const auto &[key, v] : map_S_S) {
    count++;
    if (v.size() == 1) {
      stream << "    "
             << "\"" << key << "\":"
             << " "
             << "\"" << v[0] << "\"";
    } else {
      stream << "    "
             << "\"" << key << "\":"
             << " "
             << "[";
      for (auto i = 0; i < v.size(); i++) {
        stream << "\"" << v[i] << "\"";
        if (i != v.size() - 1)
          stream << ",";
      }
      stream << "]";
    }
    if (map_S_S.size() != count)
      stream << ",";
    stream << "\n";
  }
  stream << "}";
  return stream.str();
}

inline std::ofstream &operator<<(std::ofstream &stream, const JSON &json) {
  stream << ToString(json);
  return stream;
};

/* ------------------------------------------------------ */

using JsonValue = std::variant<double, int, std::string, T6d, Tddd, Tdd, std::vector<double>>;
using JsonVector = std::vector<JsonValue>;

struct JSONoutput {
  std::map<std::string, JsonVector> map;

  void push(const std::string &key, const JsonValue &value) { map[key].emplace_back(value); }

  template <typename T> void set(const std::string &key, const std::vector<T> &vec_value) {
    for (const auto &value : vec_value) {
      this->push(key, value);
    }
  }

  void output(std::ofstream &os) {
    os << "{\n";
    for (const auto &[key, values] : map) {
      os << "    \"" << key << "\": [";
      for (size_t i = 0; i < values.size(); ++i) {
        std::visit(
            [&os](auto &&arg) {
              using T = std::decay_t<decltype(arg)>;
              if constexpr (std::is_same_v<T, double> || std::is_same_v<T, int>) {
                os << arg;
              } else if constexpr (std::is_same_v<T, std::string>) {
                os << "\"" << arg << "\"";
              } else if constexpr (std::is_same_v<T, T6d> || std::is_same_v<T, Tdd> || std::is_same_v<T, Tddd>) {
                os << "[";
                for (size_t j = 0; j < arg.size(); ++j) {
                  os << arg[j];
                  if (j < arg.size() - 1)
                    os << ", ";
                }
                os << "]";
              } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                os << "[";
                for (size_t j = 0; j < arg.size(); ++j) {
                  // if (j == 0) os << "[";
                  os << arg[j];
                  if (j < arg.size() - 1)
                    os << ", ";
                  // if (j == arg.size() - 1) os << "]";
                }
                os << "]";
              }
            },
            values[i]);

        if (i < values.size() - 1) {
          os << ", ";
        }
      }
      os << "]";
      os << ((&key != &std::prev(map.end())->first) ? ",\n" : "\n");
    }
    os << "}\n";
  }
};

// struct JSONoutput {
//    std::map<std::string, std::vector<T6d>> map_S_T6d;
//    std::map<std::string, std::vector<Tddd>> map_S_Tddd;
//    std::map<std::string, std::vector<double>> map_S_D;
//    std::map<std::string, std::vector<int>> map_S_I;
//    std::map<std::string, std::vector<std::string>> map_S_S;
//    JSONoutput(){};

//    const double notfinite = 1E+30;
//    const Tddd notfinite_tddd = {notfinite, notfinite, notfinite};
//    const T6d notfinite_t6d = {notfinite, notfinite, notfinite, notfinite, notfinite, notfinite};

//    void push(std::string s, const T6d d) {
//       auto D = (isFinite(d) ? d : notfinite_t6d);
//       if (this->map_S_T6d.find(s) != this->map_S_T6d.end())
//          this->map_S_T6d[s].emplace_back(D);
//       else
//          this->map_S_T6d[s] = {D};
//    };
//    void push(std::string s, const Tddd d) {
//       auto D = (isFinite(d) ? d : notfinite_tddd);
//       if (this->map_S_Tddd.find(s) != this->map_S_Tddd.end())
//          this->map_S_Tddd[s].emplace_back(D);
//       else
//          this->map_S_Tddd[s] = {D};
//    };
//    void push(std::string s, double d) {
//       auto D = (isFinite(d) ? d : notfinite);
//       if (this->map_S_D.find(s) != this->map_S_D.end())
//          this->map_S_D[s].emplace_back(D);
//       else
//          this->map_S_D[s] = {D};
//    };
//    void push(std::string s, int I) {
//       if (this->map_S_I.find(s) != this->map_S_I.end())
//          this->map_S_I[s].emplace_back(I);
//       else
//          this->map_S_I[s] = {I};
//    };
//    void push(std::string s, std::string S) {
//       if (this->map_S_S.find(s) != this->map_S_S.end())
//          this->map_S_S[s].emplace_back(S);
//       else
//          this->map_S_S[s] = {S};
//    };

//    void output(std::ofstream &os) {
//       os << "{" << std::endl;
//       {
//          // int size = this->map_S_T6d.size();
//          int i = 0;
//          for (const auto &[key, V] : this->map_S_T6d)  // for eacn name
//          {
//             os << "\"" << key << "\":[";
//             for (auto j = 0; j < V.size(); ++j)  // for eacn time step
//             {
//                auto [V0, V1, V2, V3, V4, V5] = V[j];
//                os << "[" << std::setprecision(15) << V0 << ", " << V1 << ", " << V2 << ", " << V3 << ", " << V4 << ", " << V5 << "]";
//                if (j != V.size() - 1)
//                   os << ",\n";
//                else
//                   os << "],\n";
//             }
//          }
//       }
//       {
//          // int size = this->map_S_Tddd.size();
//          // int i = 0;
//          for (const auto &[key, V] : this->map_S_Tddd) {
//             os << "\"" << key << "\":[";
//             for (auto j = 0; j < V.size(); ++j) {
//                auto [V0, V1, V2] = V[j];
//                os << "[" << std::setprecision(15) << V0 << ", " << V1 << ", " << V2 << "]";
//                if (j != V.size() - 1)
//                   os << ",\n";
//                else
//                   os << "],\n";
//             }
//          }
//       }
//       {
//          int size = this->map_S_D.size();
//          int i = 0;
//          for (const auto &[key, V] : this->map_S_D) {
//             os << "\"" << key << "\":[";
//             for (auto j = 0; j < V.size(); ++j)
//                os << std::setprecision(15) << V[j] << (j != V.size() - 1 ? "," : "]");
//             if (i == size - 1)
//                os << "\n";
//             else
//                os << ",\n";
//             i++;
//          }
//       }
//       {
//          int size = this->map_S_I.size();
//          int i = 0;
//          for (const auto &[key, V] : this->map_S_I) {
//             os << "\"" << key << "\":[";
//             for (auto j = 0; j < V.size(); ++j)
//                os << V[j] << (j != V.size() - 1 ? "," : "]");
//             if (i == size - 1)
//                os << "\n";
//             else
//                os << ",\n";
//             i++;
//          }
//       }
//       os << "\n}" << std::endl;
//    };
// };

//============================================================
//========================= Load  ============================
//============================================================
// std::vector<V_d> Load(const std::string& filename, const V_s& sep){
//   std::vector<V_d> mat;
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     std::cout << Red << filename << " can not be opened" << reset << std::endl;
//   else
//     std::cout << Blue << filename << " is opened" << reset << std::endl;

//   std::string str, S;
//   double v;
//   std::size_t found, s(sep.size()), min_len;
//   bool loop;
//   while(!strm.eof( ))
//     {
//       std::getline(strm,str);
//       if(!str.empty())
// 	{
// 	  V_d row_vector;
// 	  loop = true;
// 	  while(loop)
// 	    {
// 	      min_len = 999999;
// 	      for(auto i=0; i<s; i++)
// 		if((found = str.find(sep[i]))!=std::string::npos)
// 		  if(min_len > found)
// 		    {
// 		      min_len = found;
// 		      S = sep[i];
// 		    }

// 	      if(min_len!=999999/*if not 0*/)
// 		{
// 		  row_vector.push_back((double)stod(str.substr(0,min_len)));
// 		  str = str.substr(min_len + S.length());
// 		}
// 	      else
// 		{
// 		  row_vector.push_back((double)stod(str));
// 		  loop = false;
// 		}
// 	    }
// 	  mat.emplace_back(row_vector);
// 	}
//     }
//   strm.close();
//   return mat;
// };
//////////////////////////////
inline void Load(const std::string &filename, std::vector<V_s> &ret_mat, const V_s &SEP) {
  ret_mat.clear();
  std::ifstream strm(filename, std::ios::in);
  if (!strm) {
    std::stringstream ss;
    ss << "The file can not be opened: " << filename;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
  } else
    std::cout << Blue << filename << " is opened" << colorReset << std::endl;

  std::string read;
  // int row_counter(0);
  while (!strm.eof()) // loop until the end of file
  {
    read.clear();
    std::getline(strm, read);
    if (!StringTrim(read, {" "}).empty()) {
      ret_mat.emplace_back(StringTrim(StringSplit(read, SEP), {" ", "\t", "\n", "\r"}));
    }
  }
  strm.close();
  std::cout << Blue << filename << " is closed" << colorReset << std::endl;
};

inline std::vector<V_s> Load(const std::string &filename, const V_s &SEP) {
  std::vector<V_s> ret;
  Load(filename, ret, SEP);
  return ret;
};
////////////////////////////////
inline V_d string_to_vector_double(const std::string &strIN, const V_s &SEP) {
  V_s str = StringSplit(strIN, SEP);
  V_d ret;
  for (size_t i = 0; i < str.size(); i++)
    if (!str[i].empty())
      ret.emplace_back(stod(str[i]));

  return ret;
  // std::string str(strIN);/* expexting "1,2,3,4,5" or "1, 2,3, 4,5"*/
  // V_d ret(0);
  // std::size_t found;
  // int i(0);
  // while(!str.empty())
  //   {
  //     for(auto sep: SEP)
  //       if((found = str.find(sep))!=std::string::npos)
  // 	  {
  // 	    // std::cout << str << std::endl;
  // 	    // std::cout << found << std::endl;
  // 	    // std::cout << sep << std::endl;
  // 	    // std::cout << str << std::endl;
  // 	    ret.emplace_back(stod(
  // 			       str.substr((i==0 ? 0 : found+sep.length()), found)
  // 			       ));
  // 	    str = str.substr(found+sep.length());
  // 	    break;
  // 	  }
  // 	else
  // 	  {
  // 	    ret.emplace_back(stod(
  // 			       str.substr((i==0 ? 0 : found+sep.length()))
  // 			       ));
  // 	    return ret;
  // 	  }
  //     i++;
  //   }
  // return ret;
};
///////////////
inline V_s to_cell(const std::string &strIN) {
  std::string str(strIN); /* expexting "{1,2,3,4,5} {1,2,3,4,5} {1,2,3,4,5}"*/
  V_s ret(0);
  std::size_t right, left;
  // int i(0);
  while (!str.empty()) {
    if ((left = str.find("{")) != std::string::npos) {
      if ((right = str.find("}")) != std::string::npos) {
        std::string cell = str.substr(left, right - left + 1);
        ret.emplace_back(cell);
        str = str.substr(right + 1);
        // std::cout << str << std::endl;
      } else
        return ret;
    } else
      return ret;
  }
  return ret;
};
///////////////////////////

inline void Load(const std::string &filename, VVV_d &vvv) {
  std::ifstream strm(filename, std::ios::in);
  if (!strm) {
    std::stringstream ss;
    ss << "The file cannot be opened: " << filename;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
  }
  std::string str;
  while (std::getline(strm, str)) {
    V_s str_vec = to_cell(str); // Expecting {"{1,2,3,5,5}","{1,2,3,5,5}","{1,2,3,5,5}"}
    VV_d row_vec(0);
    for (const auto &v : str_vec) {
      std::string tmp = v.substr(v.find_first_of("{") + 1, v.find_last_of("}") - 1);
      V_d vec_doub = string_to_vector_double(tmp, {",", ", "}); // Expecting {1,2,3,5,5}
      row_vec.emplace_back(vec_doub);
    }
    vvv.emplace_back(row_vec);
  }
  strm.close();
};

// /////////////////////////////////
// void Load(const std::string& filename, std::vector<V_d>& vv)
// {
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     {
//       std::cout << Red << "The file can not be opened: " << filename << reset << std::endl;
//       std::cin.ignore();
//       abort();
//     }
//   std::string str;
//   while(!strm.eof( ))
//     {
//       std::getline(strm,str);

//       // std::cout <<  str << std::endl;
//       // str = str.substr(str.find_last_of("{")+1,  str.find_last_of("}")-1);
//       // std::cout <<  str << std::endl;

//       V_d vec_doub = string_to_vector_double(/* expecting {1,2,3,5,5} */str,{",",", "});
//       vv.emplace_back(vec_doub);
//     }
//   strm.close();
// };
//============================================================
//============================================================
template <class T> std::vector<T> PointsToSurface(T t0, T t1, const std::vector<std::vector<T>> &a) {
  T EPS = 1E-15;

  if (a.size() == 3) {
    T t2 = 1. - t1 - t0;
    T x0 = 1. - t0 - 2. * t1, x1 = 1. - t1 - 2. * t2, x2 = 1. - t2 - 2. * t0;

    if (t0 < EPS) {
      return a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
    } else if (t1 < EPS) {
      return a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
    } else if (t2 < EPS) {
      return a[2] * (1. - x2) / 2. + a[0] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
    } else {
      std::vector<T> l0 = a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
      std::vector<T> l1 = a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
      std::vector<T> l2 = a[2] * (1. - x2) / 2. + a[0] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
      std::vector<T> tmp = {1 / t0, 1 / t1, 1 / t2};
      return Dot(tmp / (1 / t0 + 1 / t1 + 1 / t2), {l0, l1, l2});
    }
  } else {
    T t2 = 1. - t0, t3 = 1. - t1;
    T x0 = 2. * t3 - 1., x1 = 2. * t0 - 1., x2 = 2. * t1 - 1., x3 = 2. * t2 - 1.;

    if (t0 < EPS) {
      return a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
    } else if (t1 < EPS) {
      return a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
    } else if (t2 < EPS) {
      return a[2] * (1. - x2) / 2. + a[3] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
    } else if (t3 < EPS) {
      return a[3] * (1. - x3) / 2. + a[0] * (1. + x3) / 2. /*+p2A3*(x3^2-1)*/;
    } else {
      std::vector<T> l0 = a[0] * (1. - x0) / 2. + a[1] * (1. + x0) / 2. /*+p2A0*(x0^2-1)*/;
      std::vector<T> l1 = a[1] * (1. - x1) / 2. + a[2] * (1. + x1) / 2. /*+p2A1*(x1^2-1)*/;
      std::vector<T> l2 = a[2] * (1. - x2) / 2. + a[3] * (1. + x2) / 2. /*+p2A2*(x2^2-1)*/;
      std::vector<T> l3 = a[3] * (1. - x3) / 2. + a[0] * (1. + x3) / 2. /*+p2A3*(x3^2-1)*/;
      std::vector<T> tmp = {1 / t0, 1 / t1, 1 / t2, 1 / t3};
      return Dot(tmp / (1 / t0 + 1 / t1 + 1 / t2 + 1 / t3), {l0, l1, l2, l3});
    }
  }
};
//============================================================
template <typename T> std::vector<T> Drop(std::vector<T> ret /*copy*/, std::vector<int> vec_n) { // 2020/03/22
  std::sort(vec_n.begin(), vec_n.end());
  for (size_t i = 0; i < vec_n.size(); i++)
    ret.erase(ret.begin() + (vec_n[i] - i));
  return ret;
};

inline int Position(std::vector<int> vecIN, const int n) { // 2020/03/22
  std::vector<int>::iterator it = std::find(vecIN.begin(), vecIN.end(), n);
  if (it == vecIN.end())
    std::cout << Red << __FILE__ << "  " << __FUNCTION__ << "  " << __LINE__ << colorReset << std::endl;

  return std::distance(vecIN.begin(), it);
};

template <class T> std::vector<T> Union(const std::vector<T> &vecs, const std::vector<T> &vecs2) { // Union without sorting!
  std::vector<T> ret = vecs;
  for (const auto &v : vecs2)
    if (!MemberQ(ret, v))
      ret.emplace_back(v);
  return ret;
};

template <class T> std::vector<T> Union(const std::vector<std::vector<T>> &vecs) { // Union without sorting!
  std::vector<T> ret(0);
  for (const auto &vec : vecs)
    for (const auto &v : vec)
      if (!MemberQ(ret, v))
        ret.emplace_back(v);

  return ret;
};

template <class Type> std::vector<std::vector<Type>> SortVectorChain(const std::vector<std::vector<Type>> &vecs) { // 2020/03/22
  std::vector<Type> v = *vecs.begin();
  std::vector<std::vector<Type>> rest = Drop(vecs, {0});
  Type knot = *v.rbegin();
  int n = 0;
  std::vector<std::vector<Type>> ret = {v};
  while (rest.size() != 0) {
    std::vector<Type> tmp = rest[n];
    // std::cout << "knot= " << knot << std::endl;
    // std::cout << "tmp= " << tmp << std::endl;
    // std::cout << "MemberQ("<< tmp << "," << knot << ") = " << MemberQ(tmp, knot) << std::endl;
    if (MemberQ(tmp, knot)) {
      // std::cout << tmp << std::endl;
      // std::cout << Green << "RotateLeft(tmp, Position(tmp, knot))= " << RotateLeft(tmp, Position(tmp, knot)) << reset << std::endl;
      //	  v.clear();
      v = RotateLeft(tmp, Position(tmp, knot));
      ret.emplace_back(v);
      knot = *v.rbegin();

      rest = Drop(rest, {n});

      n = 0;
    } else
      n++;

    // std::cout << Blue << "ret= " << ret << reset << std::endl;
  }
  return ret;
}

template <class T> std::vector<T> Reverse(std::vector<T> vec) {
  std::reverse(vec.begin(), vec.end());
  return vec;
};

inline double Subtract(const Tdd &ab) { return std::get<0>(ab) - std::get<1>(ab); };

//===========================================================
#include "basic_geometry.hpp"
#include "basic_statistics.hpp"
// //! -------------------------------------------------------------------------- */
// //!                  Operations for Elements 2023/03/11                        */
// //! -------------------------------------------------------------------------- */
// Tddd gradTangential_LinearElement(const Tddd &phi012, const T3Tddd &X012) {
//    auto [X0, X1, X2] = X012;
//    auto [phi0, phi1, phi2] = phi012;
//    auto n = TriangleNormal(X012);
//    return Cross(n, phi0 * (X2 - X1) + phi1 * (X0 - X2) + phi2 * (X1 - X0)) / (2 * TriangleArea(X012));
// };

// T3Tddd gradTangential_LinearElement(const T3Tddd &V012, const T3Tddd &X012) {
//    auto [X0, X1, X2] = X012;
//    auto [Vx012, Vy012, Vz012] = Transpose(V012);
//    return {gradTangential_LinearElement(Vx012, X012),
//            gradTangential_LinearElement(Vy012, X012),
//            gradTangential_LinearElement(Vz012, X012)};
// };

/* -------------------------------------------------------------------------- */
// b% ------------------- タプルからparticlize ------------------ */
inline std::vector<std::tuple<Tddd /*実際の座標*/, Tdd /*パラメタt0t1*/>> triangleIntoPoints(const T3Tddd &X0X1X2, const double dx) {
  /*
  |-o--o--o--o-| n = 4
  x,yはパラメタ
  */
  auto intp = [&X0X1X2](double t0, double t1) { return Dot(Tddd{t0, t1, 1 - t0 - t1}, X0X1X2); };
  auto y_list = [&intp, &dx](const double x) {
    int n = std::round(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
    // int n = std::ceil(Norm(intp(x, 0.) - intp(x, 1. - x)) /*このxでのyの長さ[0,1]*/ / dx);
    V_d ret(n);
    for (auto i = 0; i < n; ++i)
      ret[i] = ((1. - x) / n) * (i + 0.5);
    return ret;
  };
  /* ------------------------------------------------------ */
  Tddd v = intp(0., 1.) - intp(0., 0.);
  Tddd u = intp(1., 0.) - intp(0., 0.);
  double height = Norm(u - Dot(Normalize(v), u) * Normalize(v));
  int n = std::round(height / dx);
  // int n = std::ceil(height / dx);
  // std::cout << Grid({"n", n, "height", height}) << std::endl;
  double dH = 1. / n;
  V_d x_list(n);
  for (auto i = 0; i < n; i++)
    x_list[i] = dH * (i + 0.5);
  /* ------------------------------------------------------ */
  std::vector<std::tuple<Tddd /*実際の座標*/, Tdd /*パラメタt0t1*/>> ret;
  ret.reserve(x_list.size() * x_list.size());
  for (const auto &x : x_list)
    for (const auto &y : y_list(x))
      ret.push_back({intp(x, y), Tdd{x, y}});
  return ret;
};

inline std::vector<std::tuple<Tddd, Tdd>> particlize(const T3Tddd &X0X1X2, const double dx) { return triangleIntoPoints(X0X1X2, dx); };

#include "lib_spatial_partitioning.hpp"

struct Load3DFile {
  V_d eachMax;
  V_d eachMin;

  VV_d v_complex;
  std::vector<Tddd> v_complex_tuple;
  std::vector<std::vector<float>> s_complex;
  VV_i f_v_complex; // これのstd::vector<int>が１面 -> 3頂点を指定

  std::vector<std::vector<float>> s;
  VV_d v, vn;

  VV_i f_v; // start from 0
  VV_i f_t;
  VV_i f_vn;

  VV_i l_v;
  VV_i l_v_complex; // これのstd::vector<int>が１線分 -> 2頂点を指定

  void calculateMaxMin() {
    eachMax.resize(3);
    eachMin.resize(3);
    for (const auto &V : v_complex)
      for (auto i = 0; i < 3; i++) {
        if (eachMax[i] < V[i])
          eachMax[i] = V[i];
        else if (eachMin[i] > V[i])
          eachMin[i] = V[i];
      }
  };

  V_d surface(int i, double t0, double t1) {
    std::vector<int> indexes = f_v[i];
    VV_d a(indexes.size());
    for (auto n = 0; n < indexes.size(); n++)
      a[n] = v[indexes[n] - 1];
    return PointsToSurface(t0, t1, a);
  };

  Load3DFile() : eachMax(3, 0.), eachMin(3, 0.) {};

  Load3DFile(const std::string &filename) : eachMax(3, 0.), eachMin(3, 0.) {
    std::vector<V_s> read_line;
    Load(filename, read_line, {"    ", "   ", "  ", " "});
    std::cout << "Load3DFile : lines of file are splited" << std::endl;
    if (filename.contains(".obj"))
      this->load(read_line);
    else if (filename.contains(".off"))
      this->load_off(read_line);
    else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "The file format is not supported");
    std::cout << "Load3DFile : file is loaded" << std::endl;
    // std::cout << Red << f_v << colorReset << std::endl;
    // std::cout << Magenta << v << colorReset << std::endl;
    generateComplex();
    std::cout << "Load3DFile : complex is generated" << std::endl;
    calculateMaxMin();
    std::cout << "Load3DFile : max and min are calculated" << std::endl;
  };

  void generateComplex() {
    int j = 0;
    std::vector<int> tri_int{0, 0, 0};
    if (s.empty()) {
      v_complex.clear();
      v_complex_tuple.clear();
      f_v_complex.clear();
      for (const auto &Ind : f_v) {
        for (const auto &i : Ind) {
          v_complex.emplace_back(v[i]);
          v_complex_tuple.emplace_back(Tddd{v[i][0], v[i][0], v[i][0]});
        }
        tri_int = {j, j + 1, j + 2};
        f_v_complex.emplace_back(tri_int);
        j = j + 3;
      }
    } else {
      v_complex.clear();
      v_complex_tuple.clear();
      s_complex.clear();
      f_v_complex.clear();
      for (const auto &Ind : f_v) {
        //	std::vector<float> tmp = {std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),std::rand()/float(RAND_MAX),1.};
        std::vector<float> tmp = {.4, .4, .4, 1.};
        for (const auto &i : Ind) {
          v_complex.emplace_back(v[i]);
          v_complex_tuple.emplace_back(Tddd{v[i][0], v[i][1], v[i][2]});
          s_complex.emplace_back(tmp);
        }
        tri_int = {j, j + 1, j + 2};
        f_v_complex.emplace_back(tri_int);
        j = j + 3;
      }
    }

    // Handling line vertices
    for (const auto &line : l_v) {
      for (const auto &i : line) {
        v_complex.emplace_back(v[i]);
        v_complex_tuple.emplace_back(Tddd{v[i][0], v[i][1], v[i][2]});
        // You may also want to handle 's_complex' for lines, if necessary
      }
      // Assuming lines are defined by pairs of vertices
      l_v_complex.push_back({j, j + 1});
      j += 2;
    }
  };

  // stringを変換しv, f_nに格納
  void load(VV_s &read_line) {
    V_s v_t_vn;
    size_t estimated_size = read_line.size(); // Estimated size based on the number of lines
    v.reserve(estimated_size);
    f_v.reserve(estimated_size);
    f_t.reserve(estimated_size);
    f_vn.reserve(estimated_size);
    l_v.reserve(estimated_size);
    vn.reserve(estimated_size);

    std::vector<int> f_v_tmp, f_t_tmp, f_vn_tmp;
    std::vector<int> l_v_tmp;
    for (auto &line : read_line) {
      auto line0 = line[0];
      if (line0.empty() || line0 == "#")
        continue;
      else if (line0 == "v") {
        line.erase(line.begin());
        v.reserve(3);
        v.emplace_back(stod(line));
      } else if (line0 == "f") {
        line.erase(line.begin());
        f_v_tmp.clear();
        f_t_tmp.clear();
        f_vn_tmp.clear();
        for (auto &l : line) {
          v_t_vn = StringSplit(l, {"/"}); //! format may be v, v/t, v//vn, v/t/vn
          if (v_t_vn.size() == 1)
            f_v_tmp.emplace_back(stoi(v_t_vn[0]) - 1);
          else if (v_t_vn.size() == 2) {
            f_v_tmp.emplace_back(stoi(v_t_vn[0]) - 1);
            if (!v_t_vn[1].empty())
              f_t_tmp.emplace_back(stoi(v_t_vn[1]) - 1);
          } else if (v_t_vn.size() == 3) {
            f_v_tmp.emplace_back(stoi(v_t_vn[0]) - 1);
            if (!v_t_vn[1].empty())
              f_t_tmp.emplace_back(stoi((v_t_vn[1])) - 1);
            if (!v_t_vn[2].empty())
              f_vn_tmp.emplace_back(stoi(v_t_vn[2]) - 1);
          }
        }
        //! f_v_tmp.size() == 3 means triangle, 2 means line, 4 means quad
        if (f_v_tmp.size() > 0)
          f_v.emplace_back(f_v_tmp);
        if (f_t_tmp.size() > 0)
          f_t.emplace_back(f_t_tmp);
        if (f_vn_tmp.size() > 0)
          f_vn.emplace_back(f_vn_tmp);
      } else if (line0 == "l") {
        line.erase(line.begin());
        l_v_tmp.clear();
        for (auto &l : line)
          l_v_tmp.emplace_back(stoi(l) - 1);
        l_v.emplace_back(l_v_tmp);
      } else if (line0 == "vn") {
        line.erase(line.begin());
        vn.emplace_back(stod(line));
      }
      // s.resize(v.size(), std::vector<float>(4, 0));  // RGBS
      // for (auto &tmp : s)
      //    tmp = {.9, .9, .9, 1.};
    }
  };

  void load_off(VV_s &read_line) {
    bool isHeaderRead = false;
    int v_size = 0, f_size = 0;
    int v_count = 0, f_count = 0;
    for (auto &line : read_line) {
      if (line[0] == "#")
        continue;
      else if (!isHeaderRead) {
        if (line[0] != "OFF") {
          std::cerr << "Invalid OFF file header!" << std::endl;
          return;
        }
        isHeaderRead = true;
      } else if (v_count < v_size) {
        v.emplace_back(stod(line));
        v_count++;
      } else if (f_count < f_size) {
        int s = stoi(line[0]);
        std::vector<int> f_v_tmp;
        for (int i = 1; i <= s; ++i)
          f_v_tmp.emplace_back(stoi(line[i]));
        f_v.emplace_back(f_v_tmp);
        f_count++;
      } else if (v_count == 0 && f_count == 0) { // Parsing sizes
        v_size = stoi(line[0]);
        f_size = stoi(line[1]);
        // e_size is ignored as it's not needed
      }
    }
  }

  void ConstructFromString(const std::string &strIN) {
    std::stringstream strm = std::stringstream(strIN);
    std::vector<V_s> read_line;
    std::string read;
    while (std::getline(strm, read)) // loop until the end of file
    {
      if (!StringTrim(read, {" "}).empty()) {
        read_line.emplace_back(StringTrim(StringSplit(read, {"    ", "   ", "  ", " "}), {"    ", "   ", "  ", " ", "\t", "\n", "\r"}));
      }
    }
    load(read_line);
    generateComplex();
    calculateMaxMin();
    glLINES lines;
    // for (auto i = 0; i < lines.f_v_complex.size(); i++) {
    //    checkAllIntersection({lines.v_complex[lines.f_v_complex[i][0]],
    //                          lines.v_complex[lines.f_v_complex[i][1]]});
    // }
  };

  // objindex をなおそう
  std::string JSON() {
    ///////////////// vertices
    std::string str = "{\n\"vertices\": [\n";
    for (size_t i = 0; i < v.size(); i++) {
      for (size_t j = 0; j < v[i].size(); j++)
        str += std::to_string(float(v[i][j])) + ", ";
      if (i != v.size() - 1)
        str += "\n";
      if (i == v.size() - 1)
        str.pop_back();
    }
    str.pop_back();
    str += "],\n";
    /////////////////// scalars
    str += "\"scalars\": [\n";
    for (size_t i = 0; i < s.size(); i++) {
      for (size_t j = 0; j < s[i].size(); j++)
        str += std::to_string(float(s[i][j])) + ", ";
      if (i != s.size() - 1)
        str += "\n";
      if (i == s.size() - 1)
        str.pop_back();
    }
    str.pop_back();
    str += "],\n";
    ///////////////// vertices_complex
    str += "\n\"vertices_complex\": [\n";
    for (size_t i = 0; i < v_complex.size(); i++) {
      for (size_t j = 0; j < v_complex[i].size(); j++)
        str += std::to_string(float(v_complex[i][j])) + ", ";
      if (i != v_complex.size() - 1)
        str += "\n";
      if (i == v_complex.size() - 1)
        str.pop_back();
    }
    str.pop_back();
    str += "],\n";
    /////////////////// scalars_complex
    str += "\"scalars_complex\": [\n";
    for (size_t i = 0; i < s_complex.size(); i++) {
      for (size_t j = 0; j < s_complex[i].size(); j++)
        str += std::to_string(float(s_complex[i][j])) + ", ";
      if (i != s_complex.size() - 1)
        str += "\n";
      if (i == s_complex.size() - 1)
        str.pop_back();
    }
    str.pop_back();
    str += "],\n";
    ////////////////// index_complex
    str += "\"indices_complex\": [\n";
    for (size_t i = 0; i < f_v_complex.size(); i++) {
      if (f_v_complex[i].size() == 3) {
        for (size_t j = 0; j < f_v_complex[i].size(); j++)
          str += std::to_string(f_v_complex[i][j]) + ", ";
        if (i != f_v_complex.size() - 1)
          str += "\n";
        if (i == f_v_complex.size() - 1)
          str.pop_back();
      } else if (f_v_complex[i].size() == 4) {
        str += std::to_string(f_v_complex[i][0]) + ", " + std::to_string(f_v_complex[i][1]) + ", " + std::to_string(f_v_complex[i][2]) + ", ";
        str += std::to_string(f_v_complex[i][2]) + ", " + std::to_string(f_v_complex[i][3]) + ", " + std::to_string(f_v_complex[i][0]) + ", ";
        ;
        if (i != f_v_complex.size() - 1)
          str += "\n";
        if (i == f_v_complex.size() - 1)
          str.pop_back();
      }
    }
    str.pop_back();
    str += "],\n";
    ////////////////// indices
    str += "\"indices\": [\n";
    for (size_t i = 0; i < f_v.size(); i++) {
      if (f_v[i].size() == 3) {
        for (size_t j = 0; j < f_v[i].size(); j++)
          str += std::to_string(f_v[i][j]) + ", ";
        if (i != f_v.size() - 1)
          str += "\n";
        if (i == f_v.size() - 1)
          str.pop_back();
      } else if (f_v[i].size() == 4) {
        str += std::to_string(f_v[i][0]) + ", " + std::to_string(f_v[i][1]) + ", " + std::to_string(f_v[i][2]) + ", ";
        str += std::to_string(f_v[i][2]) + ", " + std::to_string(f_v[i][3]) + ", " + std::to_string(f_v[i][0]) + ", ";
        ;
        if (i != f_v.size() - 1)
          str += "\n";
        if (i == f_v.size() - 1)
          str.pop_back();
      }
    }
    str.pop_back();
    str += "]\n}";
    ////////////////////
    return str;
  };
  std::string JSON(const std::string &strIN) {
    ConstructFromString(strIN);
    return JSON();
  };
};
inline std::string obj2json(const std::string &filename) {
  Load3DFile loaded(filename);
  return loaded.JSON();
};

template <class T> std::vector<T> BarycentriCoordinate(const std::vector<std::vector<T>> &a, const T t0, const T t1) {
  switch (a.size()) {
  case 3:
    return Dot({t0, 1. - t0 - t1, t1}, a);
  case 4:
    return Dot({t0, 1. - t0, t1, 1. - t1}, a);
  }
};
//=====================================

// V_d TriangleArea(const VV_d& abcd){
//   V_d ret(abcd.size());
//   ret[0] = TriangleArea(abcd[abcd.size()-1], abcd[0], abcd[1]);
//   for(size_t i=0; i<abcd.size()-2; i++)
//     ret[i+1] = TriangleArea(abcd[i],abcd[i+1],abcd[i+2]);
//   ret[abcd.size()-1] = TriangleArea(abcd[abcd.size()-2],abcd[abcd.size()-1],abcd[0]);
//   return ret;
// };

inline std::vector<int> removePositiveMinArea(const V_d &area, VV_d &vertices, std::vector<int> &indices) {
  V_d min_area;
  std::vector<int> min_ind, k_at_min;
  for (size_t k = 0; k < indices.size(); k++) {
    if (area[k] > 0) {
      if ((min_area.empty() ? 1E+30 : min_area[0]) > area[k]) {
        min_area.insert(min_area.begin(), area[k]);
        min_ind.insert(min_ind.begin(), indices[k]);
        k_at_min.insert(k_at_min.begin(), k);
      } else {
        min_area.emplace_back(area[k]);
        min_ind.emplace_back(indices[k]);
        k_at_min.emplace_back(k);
      }
    }
  }

  int select_num = k_at_min[0]; // 20200716のところ最小のものから埋めていく他ない

  for (size_t k = 0; k < k_at_min.size(); k++) {
    if (area[k_at_min[k]] > Mean(area) * 0.01 /*too small*/) {
      select_num = k_at_min[k];
      break;
    }
  }

  if (area[select_num] < 0.01)
    Print(area[select_num], Red);

  // eraseの前にretを準備
  std::vector<int> ret = {indices[(select_num + 1) % indices.size()], indices[select_num - 1 == -1 ? indices.size() - 1 : select_num - 1]};
  // area.erase(area.begin()+select_num);
  vertices.erase(vertices.begin() + select_num);
  indices.erase(indices.begin() + select_num);
  // 結んだ点の番号

  return ret;
};

inline std::vector<int> removeMinArea(const V_d &area, VV_d &vertices, std::vector<int> &indices) {
  V_d min_area;
  std::vector<int> min_ind, k_at_min;
  for (size_t k = 0; k < indices.size(); k++) {
    if (area[k] > 0 && (min_area.empty() ? 1E+30 : min_area[0]) > std::abs(area[k])) {
      min_area.insert(min_area.begin(), area[k]);
      min_ind.insert(min_ind.begin(), indices[k]);
      k_at_min.insert(k_at_min.begin(), k);
    } else if (std::abs(area[k]) > 0) {
      min_area.emplace_back(area[k]);
      min_ind.emplace_back(indices[k]);
      k_at_min.emplace_back(k);
    }
  }

  int select_num = k_at_min[0]; // 20200716のところ最小のものから埋めていく他なない

  // eraseの前にretを準備
  std::vector<int> ret = {indices[select_num + 1 == indices.size() ? 0 : select_num + 1], indices[select_num - 1 == -1 ? indices.size() - 1 : select_num - 1]};
  // area.erase(area.begin()+select_num);
  vertices.erase(vertices.begin() + select_num);
  indices.erase(indices.begin() + select_num);
  // 結んだ点の番号

  return ret;
};

// VV_i PolygonTriangulation(const VV_d& abcd,
// 						   const V_d& n){
//   int s = abcd.size(), i=0;
//   VV_d vertices = abcd;
//   std::vector<int> indices(s);
//   for(auto& ind:indices)
//     ind = i++;
//   VV_i connectInd;
//   while(true){
//     connectInd.emplace_back(removePositiveMinArea(DirectedArea(vertices,n),
// 					       vertices,
// 					       indices));
//     if(vertices.size()==3)
//       break;
//   }
//   return connectInd;
// };

// VV_i PolygonTriangulation(const VV_d& abcd){
//   int s = abcd.size(), i=0;
//   VV_d vertices = abcd;//順番に並んでいる
//   std::vector<int> indices(s);
//   for(auto& ind:indices)
//     ind = i++;
//   VV_i connectInd;
//   while(true){
//     connectInd.emplace_back(removeMinArea(TriangleArea(vertices),
// 				       vertices,
// 				       indices));
//     if(vertices.size()==3)
//       break;
//   }
//   return connectInd;
// };

// template <class T>
// std::vector<T> DeleteDuplicates(const std::vector<T> &vec)
// {
// 	std::vector<T> ret(0);
// 	ret.reserve(vec.size());
// 	for (const auto &v : vec)
// 	{
// 		if (std::find(ret.begin(), ret.end(), v) == ret.end())
// 			ret.emplace_back(v);
// 	}
// 	return ret;
// };
/* -------------------------------------------------------------------------- */
template <typename T> std::vector<T> DeleteDuplicates(std::vector<T> V) {
  std::sort(V.begin(), V.end());
  V.erase(std::unique(V.begin(), V.end()), V.end());
  return V;
}
// template <typename T> std::vector<T *> DeleteDuplicates(const std::vector<T *> &V) {
//   std::vector<T *> ret(0);
//   for (const auto &v : V)
//     if (!std::binary_search(ret.begin(), ret.end(), v))
//       ret.insert(std::lower_bound(ret.begin(), ret.end(), v), v);
//   return ret;
// };
template <typename T> bool DuplicateFreeQ(const std::vector<T> &vec) {
  for (auto i = 0; i < vec.size(); i++)
    for (auto j = i + 1; j < vec.size(); j++)
      if (vec[i] == vec[j])
        return false;
  return true;
};

template <typename T> bool DuplicateFreeQ(const std::tuple<T, T, T> &V) {
  auto [t0, t1, t2] = V;
  return (t0 != t1 && t0 != t2 && t1 != t2);
};
template <typename T> bool DuplicateFreeQ(const std::tuple<T, T, T, T> &V) {
  auto [t0, t1, t2, t3] = V;
  return (t0 != t1 && t0 != t2 && t0 != t3 && t1 != t2 && t1 != t3 && t2 != t3);
};

template <typename... Args> bool DuplicateFreeQ(Args... args) {
  std::set<std::common_type_t<Args...>> values = {args...};
  return values.size() == sizeof...(args);
}
/* -------------------------------------------------------------------------- */
template <typename T> T Product(const std::vector<T> &vec) {
  T ret = 1;
  for (const auto &v : vec)
    ret = ret * v;
  return ret;
};
/* -------------------------------------------------------------------------- */
inline VV_d nearestPointsOfLines(const VV_d &a0a1, const VV_d &b0b1) {
  V_d a1 = a0a1[1], a0 = a0a1[0];
  V_d b1 = b0b1[1], b0 = b0b1[0];
  V_d va = a1 - a0;
  V_d vb = b1 - b0;
  V_d a0Tob0 = b0b1[0] - a0a1[0];
  V_d n = Cross(va, vb), n1 = Cross(va, n), n2 = Cross(vb, n);
  return {a0 + va * Dot(a0Tob0, n2) / Dot(va, n2), b0 + vb * Dot(-(a0Tob0), n1) / Dot(vb, n1)};
};
inline V_d midPointOfLines(const VV_d &a0a1, const VV_d &b0b1) {
  VV_d points = nearestPointsOfLines(a0a1, b0b1);
  return (points[0] + points[1]) / 2.;
};
//=============
// 20200817
inline V_d midPointOfLines(const VVV_d &vecs) {
  // vecs is {{a0,a1},{b0,b1},{c0,c1},...}
  int n = vecs.size() - 1;
  VV_d A(vecs.size(), V_d(3));
  VV_d V(vecs.size(), V_d(3));

  for (size_t i = 0; i < vecs.size(); i++) {
    A[i] = vecs[i][0];
    V[i] = vecs[i][1] - vecs[i][0];
  }

  V_d rhs(vecs.size());

  for (size_t j = 0; j < vecs.size(); j++)
    for (size_t i = 0; i < vecs.size(); i++)
      rhs[j] -= (double)((i == j) ? -n : 1.) * Dot(V[j], A[i]);

  VV_d tmp(vecs.size(), V_d(vecs.size()));

  for (size_t j = 0; j < vecs.size(); j++)
    for (size_t i = 0; i < vecs.size(); i++)
      tmp[i][j] = (double)((i == j) ? -n : 1.) * Dot(V[j], V[i]);

  auto param = Dot(rhs, Inverse(tmp));

  V_d ret{0, 0, 0};
  for (size_t i = 0; i < vecs.size(); i++)
    ret += A[i] + V[i] * param[i];

  return ret / ((double)(vecs.size()));
};

// 20200818
#include "basic_mathematical_functions.hpp"
//-----------------------
template <class T> std::vector<T> Intersection(const std::vector<T> &A, const std::vector<T> &B) {
  std::vector<T> ret;
  ret.reserve(std::min(A.size(), B.size()));
  for (const auto &b : B) {
    if (std::find_if(A.cbegin(), A.cend(), [&](const auto &a) { return a == b; }) != A.end())
      ret.emplace_back(b);
  }
  return ret;
};
template <class T> std::vector<T> Intersection(const std::unordered_set<T> &A, const std::unordered_set<T> &B) {
  std::vector<T> ret;
  ret.reserve(std::min(A.size(), B.size()));
  for (const auto &b : B) {
    if (A.find(b) != A.end())
      ret.emplace_back(b);
  }
  return ret;
};
template <typename... Args> std::vector<std::tuple<Args...>> Intersection(const std::vector<std::tuple<Args...>> &a, const std::vector<std::tuple<Args...>> &b) {
  std::vector<std::tuple<Args...>> result;
  for (const auto &x : a)
    if (std::find(b.begin(), b.end(), x) != b.end())
      result.push_back(x);
  return result;
}

template <class T> std::vector<T> Intersection(const std::vector<T> &A, const std::vector<T> &B, const std::vector<T> &C) {
  std::vector<T> ret;
  std::unordered_set<T> A_(A.begin(), A.end());
  for (const auto &a : A_)
    for (const auto &b : B)
      if (a == b)
        for (const auto &c : C)
          if (a == c)
            ret.emplace_back(c);
  return ret;
};
template <class T> std::unordered_set<T> Intersection(const std::unordered_set<T> &A, const std::unordered_set<T> &B, const std::unordered_set<T> &C) {
  std::unordered_set<T> ret;
  std::unordered_set<T> A_(A.begin(), A.end());
  for (const auto &a : A_)
    for (const auto &b : B)
      if (a == b)
        for (const auto &c : C)
          if (a == c)
            ret.emplace(c);
  return ret;
};

/* ------------------------------------------------------ */
/*                       連立1次方程式の解法                 */
/* ------------------------------------------------------ */

template <typename U, typename T> class SLE {
  // マップから行列を生成する
  //  mapと違い，行列は，行列と，それにかかるベクトルの並びが正しくなくてはならない．
  // 正しい並びで生成する方法を考える．
  //*1)indexをもつ新たなmapを生成し，それを使って，行列に割り振っていく
  //*2)pointerをソートして，そのソートした順番に従う．バイナリーサーチ
public:
  VV_d A;
  V_d b;
  std::vector<T *> keys;
  std::vector<U *> eqs_keys;
  //! ------------------------------------------------------ */
  SLE(const std::map<T *, std::map<T *, double>> &A_IN, const std::map<T *, double> &b_IN) : keys(A_IN.size(), nullptr), A(A_IN.size(), V_d(A_IN.size(), 0.)), b(b_IN.size(), 0.) {
    if (A_IN.size() != b_IN.size()) {
      std::stringstream ss;
      ss << "A.size() = " << A_IN.size() << ", b.size() = " << b.size();
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    }
    //* -------------------- 使用するキーを生成，ソート ------------------- */
    // このキーは初期に作成したあとは，固定．
    this->keys = TakeFirst(A_IN);
    std::sort(this->keys.begin(), this->keys.end());
    // この順番がマトリックスの添字となる
    // これ以降index()が使える
    //* -------------------- 行列とベクトルを生成する -------------------- */
    int i = 0;
    for (const auto &[p, m] : A_IN) // 必ず全ての点がある
    {
      auto &row = this->A[index(p)];
      for (const auto &[q, v] : m)
        if ((i = index(q)) >= 0)
          row[i] = v;
    }
    for (const auto &[p, v] : b_IN)
      this->b[index(p)] = v;
    Print("eq has made", Red);
    Print("A.size() = " + std::to_string(A.size()), Red);
    Print("b.size() = " + std::to_string(b.size()), Red);
    Print("keys.size() = " + std::to_string(keys.size()), Red);
  };
  //! ------------------------------------------------------ */
  int index(const T *key) const {
    if (std::binary_search(this->keys.cbegin(), this->keys.cend(), key)) {
      return std::distance(this->keys.cbegin(), std::lower_bound(this->keys.cbegin(), this->keys.cend(), key));
    } else
      return -1;
  };
  //! ------------------------------------------------------ */
  int eqs_index(const T *eq_key) const {
    if (std::binary_search(this->eqs_keys.cbegin(), this->eqs_keys.cend(), eq_key)) {
      return std::distance(this->eqs_keys.cbegin(), std::lower_bound(this->eqs_keys.cbegin(), this->eqs_keys.cend(), eq_key));
    } else
      return -1;
  };
  //! ------------------------------------------------------ */
  SLE(const std::vector<T *> &keys_IN /*variables*/) : keys(keys_IN) {
    //* -------------------- 使用するキーを生成，ソート ------------------- */
    std::sort(this->keys.begin(), this->keys.end());
  };
  //
  SLE(const std::vector<U *> &eqs_keys_IN, const std::vector<T *> &keys_IN) : eqs_keys(eqs_keys_IN), keys(keys_IN) {
    std::sort(this->keys.begin(), this->keys.end());
    std::sort(this->eqs_keys.begin(), this->eqs_keys.end());
  };
  //! ------------------------------------------------------ */
  // 単純に，方程式を追加したい場合，方程式に対してインデックスを割り振る必要はない．
  void add(const std::map<U *, std::map<T *, double>> &A_IN, const std::map<U *, double> &b_IN = {}) {
    //* ------------------------------------------------------ */
    auto eqn = TakeFirst(A_IN);
    std::sort(eqn.begin(), eqn.end());
    auto index_eq = [&eqn](const T *key) -> int {
      if (std::binary_search(eqn.cbegin(), eqn.cend(), key)) {
        return std::distance(eqn.cbegin(), std::lower_bound(eqn.cbegin(), eqn.cend(), key));
      } else
        return -1;
    };
    //* -------------------- 行列とベクトルを生成する -------------------- */
    int i = 0;
    VV_d newA(eqn.size(), V_d(this->keys.size(), 0.));
    V_d newb(eqn.size(), 0.);
    for (const auto &[p, m] : A_IN) // 必ず全ての点がある
    {
      auto &row = newA[index_eq(p)];
      for (const auto &[q, v] : m)
        if ((i = index(q)) >= 0)
          row[i] = v;
    }
    for (const auto &[p, v] : b_IN)
      newb[index_eq(p)] = v;
    //
    this->A = Join(this->A, newA);
    this->b = Join(this->b, newb);
  };
  //! ------------------------------------------------------ */
  // 方程式を修正する場合は，修正する方程式を識別するインデックスが必要
  // インデックスは初期に作っておく必要がある．
  void modify(const std::map<U *, std::map<T *, double>> &A_IN, const std::map<U *, double> &b_IN = {}) {
    try {
      int i = 0, j = 0;
      for (const auto &[p, m] : A_IN) // 必ず全ての点がある
      {
        j = this->eqs_index(p);
        if (j >= 0) {
          auto &row = this->A[this->eqs_index(p)];
          for (const auto &[q, v] : m)
            if ((i = this->index(q)) >= 0)
              row[i] += v; //! ここでいう修正とは，行列成分への足し算のこと
        }
      }
      for (const auto &[p, v] : b_IN) {
        j = this->eqs_index(p);
        if (j >= 0) {
          this->b[eqs_index(p)] += v;
        }
      }
    } catch (const error_message &e) {
      throw;
    }
  };
  //! ------------------------------------------------------ */
};
