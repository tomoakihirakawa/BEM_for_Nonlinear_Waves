#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <type_traits>

#if __has_include(<quadmath.h>) && !defined(_LIBCPP_VERSION)
#include <quadmath.h>
#define BEM_HAS_QUADMATH 1
#else
#define BEM_HAS_QUADMATH 0
#endif

#if __has_include(<stdfloat>)
#include <stdfloat>
#endif

#include "basic_arithmetic_array_operations.hpp"

/* -------------------------------------------------------------------------- */

// constexpr std::array<double, 2> sum_and_arounding_error(const double sum, const double input) {
//    const double t = sum + input;
//    if ((std::bit_cast<uint64_t>(sum) & 0x7FFFFFFFFFFFFFFFULL) < (std::bit_cast<uint64_t>(input) & 0x7FFFFFFFFFFFFFFFULL))
//       return {t, sum - (t - input)};
//    else
//       return {t, input - (t - sum)};
// }

// constexpr std::array<double, 2> sum_and_arounding_error(const double a, const double b) {
//    const double sum = a + b;
//    // if (std::abs(a) < std::abs(b))
//    if ((std::bit_cast<uint64_t>(a) & 0x7FFFFFFFFFFFFFFFULL) < (std::bit_cast<uint64_t>(b) & 0x7FFFFFFFFFFFFFFFULL))
//       return {sum, a - (sum - b)};
//    else
//       return {sum, b - (sum - a)};
// }

// 互換名が必要ならエイリアス
static inline __attribute__((always_inline)) std::array<double, 2> sum_and_arounding_error(const double a, const double b) noexcept {
  const double s = a + b;
  const double z = s - a;
  return {s, (a - (s - z)) + (b - z)};
}

inline constexpr std::array<double, 2> two_prod(const double a, const double b) {
  const double p = a * b;
  return {p, std::fma(a, b, -p)};
}

#if !BEM_HAS_QUADMATH

// Fallback implementation for platforms without libquadmath (e.g. Apple clang/libc++).
// This is sufficient for building and for code paths that don't require true __float128 support.
struct DoubleDouble {
  double a = 0.0;
  double b = 0.0;

  constexpr DoubleDouble() = default;
  constexpr DoubleDouble(double v) : a(v), b(0.0) {}
  constexpr DoubleDouble(double h, double l) : a(h), b(l) {}

  explicit constexpr operator double() const { return a + b; }
  constexpr operator long double() const { return static_cast<long double>(a) + static_cast<long double>(b); }

  constexpr DoubleDouble operator+(const DoubleDouble &rhs) const {
    auto [s1, e1] = sum_and_arounding_error(this->a, rhs.a);
    auto [sum, err] = sum_and_arounding_error(s1, this->b + rhs.b + e1);
    return DoubleDouble(sum, err);
  }
  constexpr DoubleDouble operator-(const DoubleDouble &rhs) const { return *this + DoubleDouble(-rhs.a, -rhs.b); }
  constexpr DoubleDouble operator-() const { return DoubleDouble(-a, -b); }

  constexpr DoubleDouble operator*(const DoubleDouble &rhs) const {
    auto [p1, e1] = two_prod(this->a, rhs.a);
    auto [s, e2] = sum_and_arounding_error(e1, this->a * rhs.b + this->b * rhs.a);
    auto [hi, lo] = sum_and_arounding_error(p1, s);
    return DoubleDouble(hi, lo + e2);
  }

  constexpr DoubleDouble operator/(const DoubleDouble &rhs) const {
    const double denom = rhs.a + rhs.b;
    const double numer = this->a + this->b;
    double q = numer / denom;
    DoubleDouble qq(q);
    // One refinement: qq = qq + (this - rhs*qq)/rhs
    DoubleDouble r = (*this) - rhs * qq;
    q += static_cast<double>(r) / denom;
    return DoubleDouble(q);
  }

  constexpr DoubleDouble &operator+=(const DoubleDouble &rhs) { return (*this = *this + rhs); }
  constexpr DoubleDouble &operator-=(const DoubleDouble &rhs) { return (*this = *this - rhs); }
  constexpr DoubleDouble &operator*=(const DoubleDouble &rhs) { return (*this = *this * rhs); }
  constexpr DoubleDouble &operator/=(const DoubleDouble &rhs) { return (*this = *this / rhs); }

  constexpr auto operator<=>(const DoubleDouble &rhs) const {
    if (auto cmp = this->a <=> rhs.a; cmp != 0)
      return cmp;
    return this->b <=> rhs.b;
  }
  constexpr bool operator==(const DoubleDouble &other) const { return a == other.a && b == other.b; }
};

inline DoubleDouble sqrt(const DoubleDouble &x) {
  const double v = static_cast<double>(x);
  if (v < 0.0)
    return DoubleDouble(NAN, NAN);
  return DoubleDouble(std::sqrt(v));
}

inline DoubleDouble Dot(const std::array<DoubleDouble, 3> &A, const std::array<DoubleDouble, 3> &B) { return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]; }

inline std::array<DoubleDouble, 3> Cross(const std::array<DoubleDouble, 3> &A, const std::array<DoubleDouble, 3> &B) {
  return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}

inline std::array<DoubleDouble, 3> Normalize(const std::array<DoubleDouble, 3> &A) {
  const DoubleDouble unity = 1.0;
  const DoubleDouble inv_norm = unity / sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
  return {A[0] * inv_norm, A[1] * inv_norm, A[2] * inv_norm};
}

template <typename T>
inline constexpr std::complex<T> Polar(const T &r, const T &theta) {
  return std::polar(r, theta);
}

inline std::ostream &operator<<(std::ostream &os, const DoubleDouble &q) { return os << static_cast<double>(q); }

#else

struct DoubleDouble {
  double a = 0.0;
  double b = 0.0;

  // コンストラクタをconstexprにする (最重要)
  constexpr DoubleDouble() = default;

  constexpr DoubleDouble(double h, double l) : a(h), b(l) {}

  constexpr DoubleDouble(const std::float128_t &q_val)
      // 1. q_valをdoubleにキャストして高位部分(a)を求める
      : a(static_cast<double>(q_val)), b(static_cast<double>(q_val - static_cast<std::float128_t>(a))) {}

  // 各演算子をconstexprにする
  constexpr DoubleDouble operator+(const DoubleDouble &rhs) const {
    auto [s1, e1] = sum_and_arounding_error(this->a, rhs.a);
    auto [sum, err] = sum_and_arounding_error(s1, this->b + rhs.b + e1);
    return DoubleDouble(sum, err);
  }

  constexpr DoubleDouble operator-(const DoubleDouble &rhs) const {
    DoubleDouble result;
    result.a = -rhs.a;
    result.b = -rhs.b;
    return *this + result;
  }

  constexpr DoubleDouble operator*(const DoubleDouble &rhs) const {
    auto [p1, e1] = two_prod(this->a, rhs.a);
    auto [s, e2] = sum_and_arounding_error(e1, this->a * rhs.b + this->b * rhs.a);
    auto [a, b] = sum_and_arounding_error(p1, s);
    return DoubleDouble(a, b + e2);
  }

  constexpr DoubleDouble operator/(const DoubleDouble &rhs) const {
    // コンパイル時評価で0除算が発生するとエラーになる
    if (std::is_constant_evaluated() && rhs.a == 0.0 && rhs.b == 0.0) {
      // コンパイル時エラーを明確にする
      // (throwは実行時エラー用)
    }
    DoubleDouble inv_rhs(1.q / static_cast<std::float128_t>(rhs));
    const DoubleDouble two(static_cast<std::float128_t>(2.0));
    inv_rhs = inv_rhs * (two - rhs * inv_rhs);
    return (*this * inv_rhs);
  }

  // DoubleDoubleクラスの定義内に追加
  DoubleDouble &operator/=(const DoubleDouble &rhs) {
    // 既に実装済みの / 演算子を使って実装するのが簡単
    *this = *this / rhs;
    return *this;
  }

  constexpr DoubleDouble operator/(const std::float128_t &rhs) const {
    // コンパイル時評価で0除算が発生するとエラーになる
    if (std::is_constant_evaluated() && rhs == 0.0q) {
      // コンパイル時エラーを明確にする
      // (throwは実行時エラー用)
    }
    return *this / DoubleDouble(rhs);
  }

  // 複合代入演算子
  constexpr DoubleDouble &operator+=(const DoubleDouble &rhs) { return (*this = *this + rhs); }
  constexpr DoubleDouble &operator-=(const DoubleDouble &rhs) {
    *this = *this - rhs;
    return *this;
  }
  constexpr DoubleDouble &operator*=(const DoubleDouble &rhs) { return (*this = *this * rhs); }

  constexpr DoubleDouble &operator*=(const std::float128_t &rhs) {
    DoubleDouble tmp(rhs);
    *this = *this * tmp;
    return *this;
  }

  // doubleへのキャスト
  explicit constexpr operator double() const { return this->a + this->b; }

  // constexpr operator std::float128_t() const {
  constexpr operator std::float128_t() const {
    // 各要素を先に__float128に変換してから足すことで、高精度を維持する
    return static_cast<std::float128_t>(this->a) + static_cast<std::float128_t>(this->b);
  }

  constexpr operator long double() const { return static_cast<long double>(static_cast<std::float128_t>(this->a) + static_cast<std::float128_t>(this->b)); }

  // 【ここに追加】宇宙船演算子 (三方比較)
  constexpr auto operator<=>(const DoubleDouble &rhs) const {
    if (auto cmp = this->a <=> rhs.a; cmp != 0)
      return cmp;
    return this->b <=> rhs.b;
  }

  // ↓ doubleやintとの比較も可能にする
  constexpr auto operator<=>(double rhs) const { return *this <=> DoubleDouble(rhs); }

  constexpr bool operator==(const DoubleDouble &other) const { return a == other.a && b == other.b; }

  // -を定義
  constexpr DoubleDouble operator-() const { return DoubleDouble(-a, -b); }
  // >を定義float128_tとの比較も可能にする
  constexpr bool operator>(const std::float128_t &rhs) const { return static_cast<std::float128_t>(*this) > rhs; }

  constexpr bool operator<(const std::float128_t &rhs) const { return static_cast<std::float128_t>(*this) < rhs; }
};

inline DoubleDouble sqrt(const DoubleDouble &a) {
  // 負の値の平方根はエラー（非数を返すなど）
  if (a < 0.0) {
    return DoubleDouble(NAN, NAN);
  }
  // 簡単な実装：一度__float128に変換して計算し、DoubleDoubleに戻す
  // これでも精度が落ちる可能性はあるが、コンパイルエラーは解決できる
  return DoubleDouble(sqrtq(static_cast<__float128>(a)));
}

inline DoubleDouble Dot(const std::array<DoubleDouble, 3> &A, const std::array<DoubleDouble, 3> &B) { return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]; }

inline std::array<DoubleDouble, 3> Cross(const std::array<DoubleDouble, 3> &A, const std::array<DoubleDouble, 3> &B) {
  DoubleDouble A1 = std::get<0>(A);
  DoubleDouble A2 = std::get<1>(A);
  DoubleDouble A3 = std::get<2>(A);
  DoubleDouble B1 = std::get<0>(B);
  DoubleDouble B2 = std::get<1>(B);
  DoubleDouble B3 = std::get<2>(B);
  return {A2 * B3 - A3 * B2, A3 * B1 - A1 * B3, A1 * B2 - A2 * B1};
}

inline std::array<DoubleDouble, 3> Normalize(const std::array<DoubleDouble, 3> &A) {
  DoubleDouble A1 = std::get<0>(A);
  DoubleDouble A2 = std::get<1>(A);
  DoubleDouble A3 = std::get<2>(A);
  DoubleDouble unity = 1.;
  DoubleDouble _norm = unity / sqrt(A1 * A1 + A2 * A2 + A3 * A3);
  return {A1 * _norm, A2 * _norm, A3 * _norm};
}

inline std::tuple<DoubleDouble, std::array<DoubleDouble, 3>> Nearest_(const std::array<DoubleDouble, 3> &X, const std::array<std::array<DoubleDouble, 3>, 2> &ab) {
  /*
  a * t + b * (1-t)
  ---------------------------
  ( a*t+b*(1-t) - X ).(a-b) = 0
  ( (a-b)*t + (b - X) ).(a-b) = 0
  t = (X-b).(a-b)/(a-b).(a-b)
  */
  // const auto a_b = std::get<0>(ab) - std::get<1>(ab);
  // const auto t = std::clamp(Dot(X - std::get<1>(ab), a_b) / Dot(a_b, a_b), 0.0, 1.0);
  // return {t, std::get<0>(ab) * t + std::get<1>(ab) * (1. - t)};

  const std::array<DoubleDouble, 3> a = std::get<0>(ab);
  const std::array<DoubleDouble, 3> b = std::get<1>(ab);
  const std::array<DoubleDouble, 3> a_b = a - b;
  const DoubleDouble dot_a_b = Dot(a_b, a_b);
  if (dot_a_b == static_cast<DoubleDouble>(1e-30))
    return {0.5, 0.5 * (a + b)};
  const DoubleDouble t = std::clamp(Dot(X - b, a_b) / dot_a_b, static_cast<DoubleDouble>(0.0), static_cast<DoubleDouble>(1.0));
  return {t, a * t + b * (static_cast<DoubleDouble>(1.0) - t)};
};

inline std::tuple<DoubleDouble, DoubleDouble, std::array<DoubleDouble, 3>, std::array<DoubleDouble, 3>> Nearest_(const std::array<DoubleDouble, 3> &X_IN, const std::array<std::array<DoubleDouble, 3>, 3> &abc) {
  const auto [a_, b_, c_] = abc;
  const std::array<DoubleDouble, 3> a = a_;
  const std::array<DoubleDouble, 3> b = b_;
  const std::array<DoubleDouble, 3> c = c_;
  const std::array<DoubleDouble, 3> X = X_IN;
  const std::array<DoubleDouble, 3> X_a = X - a;
  const std::array<DoubleDouble, 3> b_a = b - a;
  const std::array<DoubleDouble, 3> c_a = c - a;
  const std::array<DoubleDouble, 3> x = Normalize(b_a);
  const std::array<DoubleDouble, 3> z = Normalize(Cross(b - a, c - a));
  const std::array<DoubleDouble, 3> y = Cross(z, x);
  const DoubleDouble Xx = Dot(X_a, x);
  const DoubleDouble Xy = Dot(X_a, y);
  const DoubleDouble Bx = Dot(b_a, x);
  const DoubleDouble By = Dot(b_a, y);
  const DoubleDouble Cx = Dot(c_a, x);
  const DoubleDouble Cy = Dot(c_a, y);
  const DoubleDouble den = By * Cx - Bx * Cy;
  const DoubleDouble t0 = (By * Cx - Bx * Cy + (-By + Cy) * Xx + (Bx - Cx) * Xy) / den;
  const DoubleDouble t1 = -(Cy * Xx - Cx * Xy) / den;

  if (0 >= t0 && 0 >= t1)
    return {0., 0., c, z};
  else if (0 >= t0 && 1 <= t1)
    return {0., 1., b, z};
  else if (1 <= t0 && 0 >= t1)
    return {1., 0., a, z};
  else if (0 <= t0 && t0 <= static_cast<DoubleDouble>(1.) && 0 <= t1 && t1 <= static_cast<DoubleDouble>(1.) - t0)
    return {t0, t1, a * t0 + b * t1 + c * (static_cast<DoubleDouble>(1.) - t0 - t1), z};
  else if (0 >= t0) {
    auto [t1, closestX] = Nearest_(X, std::array<std::array<DoubleDouble, 3>, 2>{b, c});
    return {0., t1, closestX, z};
  } else if (0 >= t1) {
    auto [t0, closestX] = Nearest_(X, std::array<std::array<DoubleDouble, 3>, 2>{a, c});
    return {t0, 0., closestX, z};
  } else {
    auto [t, closestX] = Nearest_(X, std::array<std::array<DoubleDouble, 3>, 2>{a, b});
    return {t, static_cast<DoubleDouble>(1.) - t, closestX, z};
  }

  /* ----------------------------------- 修正後 ---------------------------------- */

  /*
        t1
        A
        |      |        t1>= 1
     -- b -----+------  t1 = 1
        | \    |
        |   \  |
        |     \|
     -- c ---- a -----  t1 = 0
        |      |        t1<=0
  t0<=0 |      | 1<=t0
      t0=0    t0=1
  */
};

inline std::tuple<DoubleDouble, DoubleDouble, std::array<DoubleDouble, 3>, std::array<DoubleDouble, 3>> DoubleDoubleNearest_(const std::array<double, 3> &X, const std::array<std::array<double, 3>, 3> &ab) {
  std::array<DoubleDouble, 3> X_dd = {DoubleDouble(X[0]), DoubleDouble(X[1]), DoubleDouble(X[2])};
  std::array<std::array<DoubleDouble, 3>, 3> ab_dd = {{
      {DoubleDouble(ab[0][0]), DoubleDouble(ab[0][1]), DoubleDouble(ab[0][2])},
      {DoubleDouble(ab[1][0]), DoubleDouble(ab[1][1]), DoubleDouble(ab[1][2])},
      {DoubleDouble(ab[2][0]), DoubleDouble(ab[2][1]), DoubleDouble(ab[2][2])},
  }};
  return Nearest_(X_dd, ab_dd);
}

inline constexpr std::float128_t operator*(const std::float128_t &lhs, const DoubleDouble &rhs) { return lhs * static_cast<std::float128_t>(rhs); }

inline constexpr std::float128_t operator*(const DoubleDouble &lhs, const std::float128_t &rhs) { return static_cast<std::float128_t>(lhs) * rhs; }

inline constexpr std::complex<std::float128_t> operator*(const std::complex<std::float128_t> &lhs, const std::complex<DoubleDouble> &rhs) {
  auto rhs_real = static_cast<std::float128_t>(rhs.real());
  auto rhs_imag = static_cast<std::float128_t>(rhs.imag()); // complex * complex の計算
  return std::complex<std::float128_t>(lhs.real() * rhs_real - lhs.imag() * rhs_imag, lhs.real() * rhs_imag + lhs.imag() * rhs_real);
}

inline constexpr std::complex<std::float128_t> operator*(const std::complex<DoubleDouble> &lhs, const std::complex<std::float128_t> &rhs) { return rhs * lhs; }

inline constexpr std::complex<std::float128_t> operator*(const std::complex<std::float128_t> &lhs, const DoubleDouble &rhs_scalar) {
  auto scalar_val = static_cast<std::float128_t>(rhs_scalar);
  return std::complex<std::float128_t>(lhs.real() * scalar_val, lhs.imag() * scalar_val);
}

inline constexpr std::complex<std::float128_t> operator*(const DoubleDouble &lhs_scalar, const std::complex<std::float128_t> &rhs) { return rhs * lhs_scalar; }

inline constexpr std::complex<std::float128_t> &operator+=(std::complex<std::float128_t> &lhs, const std::complex<DoubleDouble> &rhs) {
  lhs.real(lhs.real() + static_cast<std::float128_t>(rhs.real()));
  lhs.imag(lhs.imag() + static_cast<std::float128_t>(rhs.imag()));
  return lhs;
}

inline constexpr std::complex<DoubleDouble> &operator+=(std::complex<DoubleDouble> &lhs, const std::complex<std::float128_t> &rhs) {
  lhs.real(lhs.real() + static_cast<DoubleDouble>(rhs.real()));
  lhs.imag(lhs.imag() + static_cast<DoubleDouble>(rhs.imag()));
  return lhs;
}

inline constexpr DoubleDouble operator*(double lhs, const DoubleDouble &rhs) {
  DoubleDouble tmp;
  tmp.a = lhs;
  return tmp * rhs;
}
inline constexpr DoubleDouble operator*(const DoubleDouble &lhs, double rhs) {
  DoubleDouble tmp;
  tmp.a = rhs;
  return lhs * tmp;
}

inline std::ostream &operator<<(std::ostream &os, const DoubleDouble &q) {
  os << (q.a + q.b); // 簡単な出力：合計として表示
  return os;
}

inline constexpr std::float128_t &operator*=(std::float128_t &lhs, const DoubleDouble &rhs) { return (lhs *= static_cast<std::float128_t>(rhs)); }

inline constexpr std::complex<std::float128_t> &operator*=(std::complex<std::float128_t> &lhs, const DoubleDouble &rhs_scalar) {
  auto scalar_val = static_cast<std::float128_t>(rhs_scalar);
  lhs.real(lhs.real() * scalar_val);
  lhs.imag(lhs.imag() * scalar_val);
  return lhs;
}
template <size_t N, size_t M> constexpr void quadruple_Hadamard_product_add(const std::array<std::array<std::complex<DoubleDouble>, M>, N> &a, const std::array<std::array<std::complex<DoubleDouble>, M>, N> &b, std::array<std::array<std::complex<DoubleDouble>, M>, N> &result) {
  static_assert(N > 0 && M > 0, "N and M must be greater than 0.");

  for (std::size_t i = 0; i < N; ++i)
    for (std::size_t j = 0; j < M; ++j)
      result[i][j] += a[i][j] * b[i][j];
};

/* -------------------------------------------------------------------------- */

// long double版のヘルパー関数
inline constexpr std::array<long double, 2> sum_and_arounding_error_ld(const long double x, const long double y) {
  long double z = x + y;
  long double v = z - x;
  long double e = (x - (z - v)) + (y - v);
  return {z, e};
}

inline constexpr std::array<long double, 2> two_prod_ld(const long double a, const long double b) {
  long double p = a * b;
  // std::fma は long double のオーバーロードを持つ
  long double err = std::fma(a, b, -p);
  return {p, err};
}

// long doubleを基本とする倍々精度クラス
struct LongDoubleLongDouble {
  long double a = 0.0L;
  long double b = 0.0L;

  constexpr LongDoubleLongDouble() = default;
  // constexpr LongDoubleLongDouble(long double val) : a(val), b(0.0L) {}
  constexpr LongDoubleLongDouble(long double h, long double l) : a(h), b(l) {}

  // 他の数値型からの変換コンストラクタ
  template <typename T>
    requires std::is_arithmetic_v<T>
  constexpr LongDoubleLongDouble(const T &val) : a(static_cast<long double>(val)), b(0.0L) {}

  explicit constexpr LongDoubleLongDouble(const std::float128_t &val) {
    a = static_cast<long double>(val);
    // valからaを引く計算は高精度(float128)で行い、その結果をlong doubleにキャストする
    b = static_cast<long double>(val - static_cast<std::float128_t>(a));
  }
  // long doubleへのキャスト
  explicit constexpr operator long double() const { return this->a + this->b; }

  // 【ここに追加】 doubleへの明示的な変換演算子
  explicit constexpr operator double() const { return static_cast<double>(this->a + this->b); }

  // 演算子
  constexpr LongDoubleLongDouble operator+(const LongDoubleLongDouble &rhs) const {
    auto [s1, e1] = sum_and_arounding_error_ld(a, rhs.a);
    auto [s2, e2] = sum_and_arounding_error_ld(b, rhs.b);
    long double c = e1 + s2;
    auto [s3, e3] = sum_and_arounding_error_ld(s1, c);
    long double s4 = e2 + e3;
    auto [sum, err] = sum_and_arounding_error_ld(s3, s4);
    return LongDoubleLongDouble(sum, err);
  }

  constexpr LongDoubleLongDouble operator-() const { return LongDoubleLongDouble(-a, -b); }

  constexpr LongDoubleLongDouble operator-(const LongDoubleLongDouble &rhs) const { return *this + (-rhs); }

  constexpr LongDoubleLongDouble operator*(const LongDoubleLongDouble &rhs) const {
    auto [p1, e1] = two_prod_ld(this->a, rhs.a);
    long double p2 = this->a * rhs.b + this->b * rhs.a;
    auto [s, e2] = sum_and_arounding_error(e1, p2);
    auto [a, b] = sum_and_arounding_error(p1, s);
    return LongDoubleLongDouble(a, b + e2);
  }

  constexpr LongDoubleLongDouble operator/(const LongDoubleLongDouble &rhs) const {
    long double inv_rhs = 1.0L / static_cast<long double>(rhs);
    return *this * LongDoubleLongDouble(inv_rhs);
  }

  // 複合代入演算子
  constexpr LongDoubleLongDouble &operator+=(const LongDoubleLongDouble &rhs) {
    *this = *this + rhs;
    return *this;
  }

  constexpr LongDoubleLongDouble &operator+=(const std::float128_t &rhs) {
    LongDoubleLongDouble tmp(rhs);
    *this = *this + tmp;
    return *this;
  }

  constexpr LongDoubleLongDouble &operator-=(const LongDoubleLongDouble &rhs) {
    *this = *this - rhs;
    return *this;
  }
  constexpr LongDoubleLongDouble &operator*=(const LongDoubleLongDouble &rhs) {
    *this = *this * rhs;
    return *this;
  }
  constexpr LongDoubleLongDouble &operator/=(const LongDoubleLongDouble &rhs) {
    *this = *this / rhs;
    return *this;
  }

  // 比較演算子
  constexpr auto operator<=>(const LongDoubleLongDouble &rhs) const {
    if (auto cmp = this->a <=> rhs.a; cmp != 0) {
      return cmp;
    }
    return this->b <=> rhs.b;
  }

  template <typename T>
    requires std::is_arithmetic_v<T>
  constexpr auto operator<=>(const T &rhs) const {
    return *this <=> LongDoubleLongDouble(rhs);
  }

  constexpr bool operator==(const LongDoubleLongDouble &other) const { return a == other.a && b == other.b; }

  template <typename T>
    requires std::is_arithmetic_v<T>
  constexpr bool operator==(const T &rhs) const {
    return *this == LongDoubleLongDouble(rhs);
  }

  explicit constexpr operator std::float128_t() const {
    // 各要素を先に__float128に変換してから足すことで、高精度を維持する
    return static_cast<std::float128_t>(this->a) + static_cast<std::float128_t>(this->b);
  }
};

// --- グローバルな演算子オーバーロード ---

// スカラーとの乗算
template <typename T>
  requires std::is_arithmetic_v<T>
inline constexpr LongDoubleLongDouble operator*(const T &lhs, const LongDoubleLongDouble &rhs) {
  return LongDoubleLongDouble(lhs) * rhs;
}

template <typename T>
  requires std::is_arithmetic_v<T>
inline constexpr LongDoubleLongDouble operator*(const LongDoubleLongDouble &lhs, const T &rhs) {
  return lhs * LongDoubleLongDouble(rhs);
}

inline constexpr LongDoubleLongDouble operator*(const LongDoubleLongDouble &lhs, const std::float128_t &rhs) { return lhs * LongDoubleLongDouble(rhs); }

// std::complexとの演算
inline constexpr std::complex<LongDoubleLongDouble> operator*(const std::complex<LongDoubleLongDouble> &lhs, const std::complex<LongDoubleLongDouble> &rhs) { return std::complex<LongDoubleLongDouble>(lhs.real() * rhs.real() - lhs.imag() * rhs.imag(), lhs.real() * rhs.imag() + lhs.imag() * rhs.real()); }

template <typename T>
  requires std::is_arithmetic_v<T>
inline constexpr std::complex<LongDoubleLongDouble> operator*(const T &lhs, const std::complex<LongDoubleLongDouble> &rhs) {
  return LongDoubleLongDouble(lhs) * rhs;
}

template <typename T>
  requires std::is_arithmetic_v<T>
inline constexpr std::complex<LongDoubleLongDouble> operator*(const std::complex<LongDoubleLongDouble> &lhs, const T &rhs) {
  return lhs * LongDoubleLongDouble(rhs);
}

// std::complexの複合代入
inline constexpr std::complex<LongDoubleLongDouble> &operator+=(std::complex<LongDoubleLongDouble> &lhs, const std::complex<LongDoubleLongDouble> &rhs) {
  lhs.real(lhs.real() + rhs.real());
  lhs.imag(lhs.imag() + rhs.imag());
  return lhs;
}

// ストリーム出力
inline std::ostream &operator<<(std::ostream &os, const LongDoubleLongDouble &q) {
  os << (q.a + q.b);
  return os;
}

// 行列のアダマール積と加算
template <size_t N, size_t M> constexpr void long_double_Hadamard_product_add(const std::array<std::array<std::complex<LongDoubleLongDouble>, M>, N> &a, const std::array<std::array<std::complex<LongDoubleLongDouble>, M>, N> &b, std::array<std::array<std::complex<LongDoubleLongDouble>, M>, N> &result) {
  static_assert(N > 0 && M > 0, "N and M must be greater than 0.");

  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < M; ++j) {
      result[i][j] += a[i][j] * b[i][j];
    }
  }
};

/* -------------------------------------------------------------------------- */

template <typename T> constexpr std::complex<T> Polar(const T &r, const T &theta) {
  if constexpr (std::is_same_v<T, std::float128_t>) {
    return r * std::complex<T>(cosq(theta), sinq(theta));
  } else {
    return std::polar<T>(r, theta);
  }
}

template <typename T> constexpr std::complex<T> Polar(const double r, const T &theta) {
  if constexpr (std::is_same_v<T, std::float128_t>) {
    return static_cast<std::float128_t>(r) * std::complex<T>(cosq(theta), sinq(theta));
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    std::float128_t theta_q = static_cast<std::float128_t>(theta);
    std::float128_t cos_val = cosq(theta_q);
    std::float128_t sin_val = sinq(theta_q);
    return static_cast<DoubleDouble>(r) * std::complex<DoubleDouble>(DoubleDouble(cos_val), DoubleDouble(sin_val));
  } else {
    return std::polar(static_cast<T>(r), theta);
  }
}
inline constexpr std::complex<DoubleDouble> Polar(const DoubleDouble &r, const DoubleDouble &theta) {
  auto c = Polar(static_cast<std::float128_t>(r), static_cast<std::float128_t>(theta));
  return std::complex<DoubleDouble>(DoubleDouble(c.real()), DoubleDouble(c.imag())); // std::complex<DoubleDouble>に変換して返す
}

// template <typename T>
// std::complex<T> Polar(const int r, const T& theta) {
//    if constexpr (std::is_same_v<T, std::float128_t>) {
//       return static_cast<std::float128_t>(r) * std::complex<T>(cosq(theta), sinq(theta));
//    } else if constexpr (std::is_same_v<T, DoubleDouble>) {
//       std::float128_t theta_q = static_cast<std::float128_t>(theta);
//       std::float128_t cos_val = cosq(theta_q);
//       std::float128_t sin_val = sinq(theta_q);
//       return static_cast<DoubleDouble>(r, 0.) * std::complex<DoubleDouble>(DoubleDouble(cos_val), DoubleDouble(sin_val));
//    } else {
//       return std::polar(static_cast<T>(r), theta);
//    }
// }

/* -------------------------------------------------------------------------- */

// 実数・複素数の平方根
template <typename T> T Sqrt(const T &x) {
  if constexpr (std::is_same_v<T, std::float128_t>) {
    return sqrtq(x);
  } else {
    return std::sqrt(x);
  }
}

template <typename T> std::complex<T> Sqrt(const std::complex<T> &z) {
  if constexpr (std::is_same_v<T, std::float128_t>) {
    return csqrtq(z); // 複素数用のsqrtq
  } else {
    return std::sqrt(z);
  }
}

// 自然対数
template <typename T> T Log(const T &x) {
  if constexpr (std::is_same_v<T, std::float128_t>) {
    return logq(x);
  } else {
    return std::log(x);
  }
}

// べき乗
template <typename T> T Pow(const T &base, const T &exp) {
  if constexpr (std::is_same_v<T, std::float128_t>) {
    return powq(base, exp);
  } else {
    return std::pow(base, exp);
  }
}

inline DoubleDouble abs(const DoubleDouble &val) {
  // valが負なら符号を反転させる、というロジックで実装
  // (DoubleDoubleクラスに比較演算子<が実装されていると仮定)
  if (val < 0.0) {
    return -val;
  }
  return val;
}

template <typename T> T Abs(const T &x) {
  if constexpr (std::is_same_v<T, std::float128_t>) {
    return fabsq(x);
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    return DoubleDouble(fabsq(static_cast<std::float128_t>(x)));
  } else {
    return std::abs(x);
  }
}

template <typename T, std::size_t N> std::array<T, N> &operator*=(std::array<T, N> &lhs, const std::array<T, N> &rhs) {
  for (std::size_t i = 0; i < N; ++i)
    lhs[i] *= rhs[i];
  return lhs;
}

template <typename T, std::size_t N> std::array<T, N> operator*(std::array<T, N> lhs, const std::array<T, N> &rhs) {
  for (std::size_t i = 0; i < N; ++i)
    lhs[i] *= rhs[i];
  return lhs;
}

template <typename T, std::size_t N, std::size_t M> std::array<std::array<T, M>, N> operator*(std::array<std::array<T, M>, N> lhs, const std::array<std::array<T, M>, N> &rhs) {
  for (std::size_t i = 0; i < N; ++i)
    lhs[i] *= rhs[i];
  return lhs;
}

template <typename STREAM> auto operator<<(STREAM &stream, const std::same_as<std::float128_t> auto &q) -> STREAM & {
  char buf[128];
  quadmath_snprintf(buf, sizeof(buf), "%.36Qg", q);
  stream << buf;
  return stream;
}

template <typename STREAM> auto operator<<(STREAM &stream, const std::same_as<std::complex<std::float128_t>> auto &q) -> STREAM & {
  char buf[128];
  quadmath_snprintf(buf, sizeof(buf), "%.36Qg + %.36Qgi", std::real(q), std::imag(q));
  stream << buf;
  return stream;
}

template <typename STREAM, std::size_t N> STREAM &operator<<(STREAM &stream, const std::array<std::float128_t, N> &arr) {
  stream << "{";
  for (size_t i = 0; i < N; ++i) {
    if (i > 0)
      stream << ", ";
    stream << arr[i];
  }
  stream << "}";
  return stream;
}

template <typename STREAM, std::size_t N, std::size_t M> STREAM &operator<<(STREAM &stream, const std::array<std::array<std::float128_t, M>, N> &arr) {
  stream << "{";
  for (size_t i = 0; i < N; ++i) {
    if (i > 0)
      stream << ", ";
    stream << arr[i];
  }
  stream << "}";
  return stream;
}

#endif // !BEM_HAS_QUADMATH
