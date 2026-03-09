#ifndef lib_Fourier_H
#define lib_Fourier_H

#include <cmath>
#include <complex>
#include <execution> // C++17以上
#include <iostream>
#include <lib_quadmath.hpp>
#include <numeric> // for std::transform_reduce
#include <vector>

/* -------------------------------------------------------------------------- */

inline std::vector<std::vector<std::complex<double>>> operator*(std::vector<std::vector<std::complex<double>>> a, const std::vector<std::vector<std::complex<double>>> &b) {
  int P = b[0].size();
  for (int i = 0; i < a.size(); ++i)
    for (int j = 0; j < P; ++j)
      a[i][j] = a[i][j] * b[i][j];
  return a;
}

/* -------------------------------------------------------------------------- */

template <typename T> std::vector<std::complex<double>> DFT(const std::vector<T> &sample) {
  //! sample can be real or complex
  int s = sample.size();
  std::vector<std::complex<double>> result(s);
  std::complex<double> sum = 0;
  double c;
  for (int n = 0; n < s; ++n) {
    c = -n * 2. * M_PI / s;

    // sum = 0;
    // for (int k = 0; k <= s - 1; ++k)
    //    sum += sample[k] * std::polar(1.0, c * k);  //! sampleが複素数の場合に対応するためpolarの引数としてsampleは渡せない

    std::complex<double> sum = std::transform_reduce(std::execution::par_unseq, sample.begin(), sample.end(), std::complex<double>{0}, std::plus<>(), [&](const auto &val) {
      int k = &val - &sample[0];           //! これはポインタの差分でkを求める方法
      return val * std::polar(1.0, c * k); //! sampleが複素数の場合に対応するためpolarの引数としてsampleは渡せない
    });
    result[n] = sum / static_cast<double>(s);
  }
  return result;
}

template <typename T> std::vector<std::vector<std::complex<double>>> DFT(const std::vector<std::vector<T>> &sample2D) {
  //! sample can be real or complex
  int N = sample2D.size();    //! 行
  int M = sample2D[0].size(); //! 列
  std::vector<std::vector<std::complex<double>>> result(N, std::vector<std::complex<double>>(M));
  std::complex<double> sum = 0;
  double cx, cy;
  /*
  for (int n = 0; n < N; ++n)
     for (int m = 0; m < M; ++m) {
        sum = 0;
        cx = -n * 2 * M_PI / N;
        cy = -m * 2 * M_PI / M;
        for (int k = 0; k <= N - 1; ++k)
           for (int j = 0; j <= M - 1; ++j)
              sum += std::polar(sample2D[k][j], cx * k + cy * j);
        result[n][m] = sum / static_cast<double>(N * M);
     }
  return result;
  */
  // 上のようにまとめてDFTを行うと，O(N^2M^2)の計算量になり，計算効率が悪い．
  // このように2Dの場合は，行ごとにDFTを行い，次に列ごとにDFTを行うのが良い．
  std::vector<std::vector<std::complex<double>>> cn2D(N, std::vector<std::complex<double>>(M));
  std::vector<std::vector<std::complex<double>>> cn2Dtmp(N, std::vector<std::complex<double>>(M));
  std::vector<std::complex<double>> cn_j(N);
  // ! x方向のフーリエ変換
  for (int i = 0; i < N; ++i)
    cn2D[i] = DFT(sample2D[i]); //! 長さM
  // ! y方向のフーリエ変換 cn2D[;;][j]
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i)
      cn_j[i] = cn2D[i][j]; //! 長さN
    cn_j = DFT(cn_j);       //! 長さN
    for (int i = 0; i < N; ++i)
      cn2Dtmp[i][j] = cn_j[i];
  }
  // ! cn2Dtmp[;;][j]をcn2D[;;][j]にコピー
  return cn2Dtmp;
}

template <typename T, size_t N> std::array<std::complex<double>, N> DFT(const std::array<T, N> &sample) {
  //! sample can be real or complex
  int s = sample.size();
  double S = s;
  std::array<std::complex<double>, N> result;
  std::complex<double> sum = 0;
  double c;
  for (int n = 0; n < s; ++n) {
    sum = 0;
    c = -n * 2 * M_PI / s;
    for (int k = 0; k <= s - 1; ++k)
      sum += sample[k] * std::polar(1.0, c * k);
    result[n] = sum / S;
  }
  return result;
}

template <typename T, size_t N, size_t M> std::array<std::array<std::complex<double>, M>, N> DFT(const std::array<std::array<T, M>, N> &sample2D) {
  //! sample can be real or complex
  std::array<std::array<std::complex<double>, M>, N> result;
  std::complex<double> sum = 0;
  double cx, cy;
  std::array<std::array<std::complex<double>, M>, N> cn2D;
  std::array<std::array<std::complex<double>, M>, N> cn2Dtmp;
  std::array<std::complex<double>, N> cn_j;
  // ! x方向のフーリエ変換
  for (int i = 0; i < N; ++i)
    cn2D[i] = DFT(sample2D[i]); //! 長さM
  // ! y方向のフーリエ変換 cn2D[;;][j]
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i)
      cn_j[i] = cn2D[i][j]; //! 長さN
    cn_j = DFT(cn_j);       //! 長さN
    for (int i = 0; i < N; ++i)
      cn2Dtmp[i][j] = cn_j[i];
  }
  // ! cn2Dtmp[;;][j]をcn2D[;;][j]にコピー
  return cn2Dtmp;
}

/* -------------------------------------------------------------------------- */

template <typename T> std::vector<std::complex<double>> InverseDFT(const std::vector<T> &sample) {
  //! sample can be real or complex
  int s = sample.size();
  std::vector<std::complex<double>> result(s);
  std::complex<double> sum = 0;
  double c;
  for (int n = 0; n < s; ++n) {
    sum = 0;
    c = n * 2 * M_PI / s;
    for (int k = 0; k <= s - 1; ++k)
      sum += sample[k] * std::polar(1.0, c * k);
    result[n] = sum * (double)s;
  }
  return result;
}

template <typename T> std::vector<std::vector<std::complex<double>>> InverseDFT(const std::vector<std::vector<T>> &sample2D) {
  //! sample can be real or complex
  int N = sample2D.size();    //! 行
  int M = sample2D[0].size(); //! 列
  std::vector<std::vector<std::complex<double>>> result(N, std::vector<std::complex<double>>(M));
  std::complex<double> sum = 0;
  double cx, cy;
  std::vector<std::vector<std::complex<double>>> cn2D(N, std::vector<std::complex<double>>(M));
  // ! x方向のフーリエ変換
  for (int i = 0; i < N; ++i)
    cn2D[i] = InverseDFT(sample2D[i]); //! 長さM

  // ! y方向のフーリエ変換 cn2D[;;][j]
  std::vector<std::vector<std::complex<double>>> cn2Dtmp(N, std::vector<std::complex<double>>(M));
  std::vector<std::complex<double>> cn_j(N);
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i)
      cn_j[i] = cn2D[i][j];  //! 長さN
    cn_j = InverseDFT(cn_j); //! 長さN
    for (int i = 0; i < N; ++i)
      cn2Dtmp[i][j] = cn_j[i];
  }
  // ! cn2Dtmp[;;][j]をcn2D[;;][j]にコピー
  return cn2Dtmp;
}

// template <typename T, size_t N>
// std::array<std::complex<double>, N> InverseDFT(const std::array<T, N>& sample) {
//    //! sample can be real or complex
//    constexpr double s = static_cast<double>(N);
//    constexpr double coeff = 2. * M_PI / s;
//    std::array<std::complex<double>, N> result;
//    std::complex<double> sum = 0;
//    double c;
//    for (int n = 0; n < N; ++n) {
//       sum = 0;
//       c = n * coeff;
//       for (int k = 0; k <= s - 1; ++k)
//          sum += sample[k] * std::polar(1.0, c * k);
//       result[n] = sum * s;
//    }
//    return result;
// }

// 事前計算テーブルを生成するためのヘルパー関数を定義する名前空間
// 事前計算テーブル生成
namespace detail {
template <int N, typename T> auto make_inverse_dft_polar_table() {
  std::array<std::array<std::complex<T>, N>, N> table;

  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < N; ++k) {
      T angle = 2 * std::numbers::pi_v<T> * n * k / N;
      if constexpr (std::is_same_v<T, __float128>) {
        table[n][k] = {cosq(angle), sinq(angle)};
      } else {
        table[n][k] = {std::cos(angle), std::sin(angle)};
      }
    }
  }
  return table;
}
} // namespace detail

// Nは配列長
template <typename T, size_t N> std::array<std::complex<double>, N> InverseDFT(const std::array<T, N> &sample) {
  static constexpr auto polar_table = detail::make_inverse_dft_polar_table<N, T>();
  std::array<std::complex<double>, N> result{};
  std::complex<double> sum = 0;
  for (int n = 0; n < N; ++n) {
    sum = 0;
    const auto &polar_row = polar_table[n];
    for (int k = 0; k < N; ++k)
      sum += sample[k] * polar_row[k];
    result[n] = sum;
  }
  return result;
}

template <typename T, size_t N, size_t M> std::array<std::array<std::complex<double>, M>, N> InverseDFT(const std::array<std::array<T, M>, N> &sample2D) {
  std::array<std::array<std::complex<double>, M>, N> result;
  std::complex<double> sum = 0;
  double cx, cy;
  std::array<std::array<std::complex<double>, M>, N> cn2D;
  // ! x方向のフーリエ変換
  for (int i = 0; i < N; ++i)
    cn2D[i] = InverseDFT(sample2D[i]); //! 長さM

  // ! y方向のフーリエ変換 cn2D[;;][j]
  std::array<std::complex<double>, N> cn_j;
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i)
      cn_j[i] = cn2D[i][j];  //! 長さN
    cn_j = InverseDFT(cn_j); //! 長さN
    for (int i = 0; i < N; ++i)
      cn2D[i][j] = cn_j[i];
  }
  // ! cn2Dtmp[;;][j]をcn2D[;;][j]にコピー
  return cn2D;
}

/* -------------------------------------------------------------------------- */
template <typename T> std::vector<std::complex<double>> DiscreteConvolve(const std::vector<T> &f, const std::vector<T> &g) {
  int len = f.size() + g.size() - 1;
  std::vector<std::complex<double>> F(len, 0), G(len, 0);
  for (int i = 0; i < f.size(); ++i)
    F[i] = f[i];
  for (int i = 0; i < g.size(); ++i)
    G[i] = g[i];

  std::vector<std::complex<double>> FourierGF(len, 0);
  F = DFT(F);
  G = DFT(G);
  for (int n = 0; n < len; ++n)
    FourierGF[n] += F[n] * G[n];

  return InverseDFT(FourierGF);
}

template <typename T, size_t N> std::array<std::complex<double>, N> DiscreteConvolve(const std::array<T, N> &f, const std::array<T, N> &g) {
  int len = f.size() + g.size() - 1;
  std::array<std::complex<double>, N> F, G;
  for (int i = 0; i < f.size(); ++i)
    F[i] = f[i];
  for (int i = 0; i < g.size(); ++i)
    G[i] = g[i];

  std::array<std::complex<double>, N> FourierGF;
  F = DFT(F);
  G = DFT(G);
  for (int n = 0; n < len; ++n)
    FourierGF[n] += F[n] * G[n];

  return InverseDFT(FourierGF);
}

template <typename T> struct DiscreteConvolveClass {

  int len;
  std::vector<std::complex<double>> FourierF, FourierG, FourierGF;
  std::vector<std::complex<double>> InverseFourierGF;

  DiscreteConvolveClass(const std::vector<T> &f, const std::vector<T> &g) : len(f.size() + g.size() - 1), FourierF(len, 0), FourierG(len, 0), FourierGF(len, 0) {
    for (int i = 0; i < f.size(); ++i)
      FourierF[i] = f[i];
    for (int i = 0; i < g.size(); ++i)
      FourierG[i] = g[i];

    FourierF = DFT(FourierF);
    FourierG = DFT(FourierG);

    for (int n = 0; n < len; ++n)
      FourierGF[n] += FourierF[n] * FourierG[n];

    InverseFourierGF = InverseDFT(FourierGF);
  };
};

/* -------------------------------------------------------------------------- */
// \label{Fourier2D}
template <typename T> struct Fourier2D {
  std::size_t rows, cols; //! Matrix size + shift - 1
  std::vector<std::vector<std::complex<double>>> MatrixDFT;

  Fourier2D() {};
  Fourier2D(int size_Nk, int size_Mk, int size_Ns, int size_Ms) { reset4convolution(size_Nk, size_Mk, size_Ns, size_Ms); };

  int Nk, Mk, Ns, Ms; //! indexの最大
  void reset4convolution(int size_Nk, int size_Mk, int size_Ns, int size_Ms) {
    Nk = size_Nk - 1; //! indexの最大
    Mk = size_Mk - 1; //! indexの最大
    Ns = size_Ns - 1; //! indexの最大
    Ms = size_Ms - 1; //! indexの最大
    this->rows = Nk + Ns + 1;
    this->cols = Mk + Ms + 1; //! 1x1 1x1なら size 1必要ない
    this->MatrixDFT = std::vector<std::vector<std::complex<double>>>(rows, std::vector<std::complex<double>>(cols, std::complex<double>(0.0, 0.0)));
  };
  void reset4convolution(int N) {
    reset4convolution(N, N, N, N); //! N x N の正方行列
  };

  void reset4convolution() {
    for (auto &row : this->MatrixDFT)
      std::fill(row.begin(), row.end(), std::complex<double>(0.0, 0.0));
  };

  bool corr_row = false;
  bool corr_col = false;

  std::complex<double> fma(const std::complex<double> &a, const std::complex<double> &b, const std::complex<double> &c) {
    double real_part = std::fma(a.real(), b.real(), -a.imag() * b.imag()) + c.real();
    double imag_part = std::fma(a.real(), b.imag(), a.imag() * b.real()) + c.imag();
    return {real_part, imag_part};
  }

  // --- reverse 対応 add（shift + flip）---
  template <typename MatrixLike> void add(const MatrixLike &input, const std::array<bool, 2> &reverse = {false, false}) {
    auto MatrixPaddedDFT = DFT(fill(input, reverse));
    for (std::size_t i = 0; i < rows; ++i)
      for (std::size_t j = 0; j < cols; ++j)
        MatrixDFT[i][j] += MatrixPaddedDFT[i][j];
  }

  template <typename MatrixLike> void add(const MatrixLike &input, const std::array<bool, 2> &reverse, const MatrixLike &input2) {
    auto MatrixPaddedDFT = DFT(fill(input, reverse));
    auto MatrixPaddedDFT2 = DFT(fill(input2));
    for (std::size_t i = 0; i < rows; ++i)
      for (std::size_t j = 0; j < cols; ++j)
        MatrixDFT[i][j] += MatrixPaddedDFT[i][j] * MatrixPaddedDFT2[i][j];
    // MatrixDFT[i][j] = fma(MatrixPaddedDFT[i][j], MatrixPaddedDFT2[i][j], MatrixDFT[i][j]);
  }

  template <typename MatrixLike> void add(const MatrixLike &input, const std::array<bool, 2> &reverse, const Fourier2D<T> &inputDFT2) {
    if (inputDFT2.rows != this->rows || inputDFT2.cols != this->cols) {
      throw std::runtime_error("Input DFT size does not match the current Fourier2D size.");
    }
    auto MatrixPaddedDFT = DFT(fill(input, reverse));
    for (std::size_t i = 0; i < rows; ++i)
      for (std::size_t j = 0; j < cols; ++j) {
        MatrixDFT[i][j] += MatrixPaddedDFT[i][j] * inputDFT2.MatrixDFT[i][j];
        // MatrixDFT[i][j] = fma(MatrixPaddedDFT[i][j], inputDFT2.MatrixDFT[i][j], MatrixDFT[i][j]);
      }
  }

  template <typename MatrixLike> void add_multipy(const MatrixLike &input, const std::array<bool, 2> &reverse, const MatrixLike &input2, double a) {
    auto MatrixPaddedDFT = DFT(fill(input, reverse));
    auto MatrixPaddedDFT2 = DFT(fill(input2));
    for (std::size_t i = 0; i < rows; ++i)
      for (std::size_t j = 0; j < cols; ++j) {
        this->MatrixDFT[i][j] += MatrixPaddedDFT[i][j] * MatrixPaddedDFT2[i][j] * a;
        // this->MatrixDFT[i][j] = fma(MatrixPaddedDFT[i][j], MatrixPaddedDFT2[i][j], this->MatrixDFT[i][j]);
      }
  }

  template <typename MatrixLike> std::vector<std::vector<T>> fill(const MatrixLike &input) {
    const auto si = input.size();
    const auto sj = input[0].size();
    std::vector<std::vector<T>> fillingMatrix(std::vector<std::vector<T>>(this->rows, std::vector<T>(this->cols, 0)));
    for (std::size_t i = 0; i < si; ++i)
      for (std::size_t j = 0; j < sj; ++j)
        fillingMatrix[i][j] = input[i][j];
    return fillingMatrix;
  }

  template <typename MatrixLike> std::vector<std::vector<T>> fill(const MatrixLike &input, const std::array<bool, 2> &reverse) {
    const auto si = input.size();
    const auto sj = input[0].size();
    std::vector<std::vector<T>> fillingMatrix(std::vector<std::vector<T>>(this->rows, std::vector<T>(this->cols, 0)));
    if (reverse[0] && reverse[1]) {
      for (std::size_t i = 0; i < si; ++i)
        for (std::size_t j = 0; j < sj; ++j)
          fillingMatrix[i][j] = input[si - 1 - i][sj - 1 - j];
    } else if (reverse[0]) {
      for (std::size_t i = 0; i < si; ++i)
        for (std::size_t j = 0; j < sj; ++j)
          fillingMatrix[i][j] = input[si - 1 - i][j];
    } else if (reverse[1]) {
      for (std::size_t i = 0; i < si; ++i)
        for (std::size_t j = 0; j < sj; ++j)
          fillingMatrix[i][j] = input[i][sj - 1 - j];
    } else {
      for (std::size_t i = 0; i < si; ++i)
        for (std::size_t j = 0; j < sj; ++j)
          fillingMatrix[i][j] = input[i][j];
    }
    return fillingMatrix;
  }
};

/* ============================================================================ */
struct Quadruple {
  double hi = 0.0;
  double lo = 0.0;

  Quadruple(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}

  // 2つのdoubleの和と誤差を計算
  static std::array<double, 2> two_sum(double a, double b) {
    double s = a + b;
    double be = s - a;
    // (a - (s - be)) は a - (s - (s-a)) = a - a = 0 にはならない
    // s-be は a' に近い値であり、 a - a' で誤差の高位部分が求まる
    double err = (a - (s - be)) + (b - be);
    return {s, err};
  }

  // 2つのdoubleの積と誤差を計算 (FMAを利用)
  static std::array<double, 2> two_prod(double a, double b) {
    double p = a * b;
    double err = std::fma(a, b, -p);
    return {p, err};
  }

  // 加算
  Quadruple operator+(const Quadruple &rhs) const {
    auto [s1, e1] = two_sum(this->hi, rhs.hi);
    // ここでの加算で発生する小さな誤差は無視する（一般的な実装）
    double s2 = this->lo + rhs.lo + e1;
    auto [sum, err] = two_sum(s1, s2);
    return Quadruple(sum, err);
  }

  // 減算（タイプミスを修正）
  Quadruple operator-(const Quadruple &rhs) const { return *this + Quadruple(-rhs.hi, -rhs.lo); }

  // 乗算
  Quadruple operator*(const Quadruple &rhs) const {
    auto [p1, e1] = two_prod(this->hi, rhs.hi);
    // lo*loの項は非常に小さいため無視し、中間的な誤差も無視する（一般的な実装）
    double p2 = this->hi * rhs.lo + this->lo * rhs.hi;
    auto [sum, err] = two_sum(p1, e1 + p2);
    return Quadruple{sum, err};
  }

  // division
  Quadruple operator/(const Quadruple &rhs) const {
    if (rhs.hi == 0.0 && rhs.lo == 0.0) {
      throw std::runtime_error("Division by zero in Quadruple division.");
    }
    double inv_hi = 1.0 / rhs.hi;
    double inv_lo = -rhs.lo * inv_hi * inv_hi; // loの項は非常に小さいため無視しない
    return Quadruple(this->hi * inv_hi + this->lo * inv_lo, this->lo * inv_hi);
  }

  // 複合代入演算子
  Quadruple &operator+=(const Quadruple &rhs) {
    *this = *this + rhs;
    return *this;
  }
  Quadruple &operator-=(const Quadruple &rhs) {
    *this = *this - rhs;
    return *this;
  }
  Quadruple &operator*=(const Quadruple &rhs) {
    *this = *this * rhs;
    return *this;
  }

  // doubleへのキャスト
  operator double() const { return this->hi + this->lo; }
};

// std::complex<Quadruple> operator*(const std::complex<Quadruple>& a, const std::complex<Quadruple>& b) {
//    return std::complex<Quadruple>(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
// }

// std::complex<Quadruple>& operator*=(std::complex<Quadruple>& a, const std::complex<Quadruple>& b) {
//    a = std::complex<Quadruple>(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
//    return a;
// }

// std::complex<Quadruple> operator*(const std::complex<Quadruple>& a, const Quadruple& b) {
//    return std::complex<Quadruple>(a.real() * b - a.imag() * b, a.real() * b + a.imag() * b);
// }

// std::complex<Quadruple>& operator*=(std::complex<Quadruple>& a, const Quadruple& b) {
//    a = std::complex<Quadruple>(a.real() * b - a.imag() * b, a.real() * b + a.imag() * b);
//    return a;
// }

std::ostream &operator<<(std::ostream &os, const Quadruple &q) {
  os << (q.hi + q.lo); // 簡単な出力：合計として表示
  return os;
}

template <size_t N, size_t M> void quadruple_adamard_product_add(const std::array<std::array<std::complex<Quadruple>, M>, N> &a, const std::array<std::array<std::complex<Quadruple>, M>, N> &b, std::array<std::array<std::complex<Quadruple>, M>, N> &result) {

  static_assert(N > 0 && M > 0, "N and M must be greater than 0.");

  for (std::size_t i = 0; i < N; ++i)
    for (std::size_t j = 0; j < M; ++j)
      result[i][j] += a[i][j] * b[i][j];
};

//^ ========================================================================== */
//^                               Fourier2D                                    */
//^ ========================================================================== */

// \label{Fourier2D}
template <typename INPUT, typename TREATED, size_t Nk, size_t Mk, size_t Ns, size_t Ms> struct Fourier2DConvolution {

  static constexpr std::size_t rows = Nk + Ns + 1;
  static constexpr std::size_t cols = Mk + Ms + 1;

  std::array<std::array<std::complex<TREATED>, cols>, rows> MatrixDFT;
  std::array<std::array<std::complex<TREATED>, cols>, rows> convolution;
  std::array<std::array<std::complex<TREATED>, rows>, rows> polar_1_c_k_for_rows;
  std::array<std::array<std::complex<TREATED>, cols>, cols> polar_1_c_k_for_cols;

  void clear() { this->MatrixDFT = {}; }

  void convolve() { this->convolution = InverseDFT(this->MatrixDFT); }

  void initialize_polar_1_c_k() {
    TREATED two_pi;
    const TREATED unity = 1; // キャストは不要

    if constexpr (std::is_same_v<TREATED, __float128> || std::is_same_v<TREATED, _Float128>) {
      // float128 の場合は、quadmath.h の高精度な定数を使う
      two_pi = 2 * M_PIq;
    } else {
      // double や long double の場合は、C++20 の numbers ヘッダが最も良い
      two_pi = 2 * std::numbers::pi_v<TREATED>;
    }

    {
      const TREATED S = static_cast<TREATED>(1) / rows;
      const TREATED coeff = two_pi * S;
      for (int n = 0; n < rows; ++n)
        for (int k = 0; k < rows; ++k)
          this->polar_1_c_k_for_rows[n][k] = Polar(unity, -n * coeff * k) * S;
    }
    {
      const TREATED S = static_cast<TREATED>(1) / cols;
      const TREATED coeff = two_pi * S;
      for (int n = 0; n < cols; ++n)
        for (int k = 0; k < cols; ++k)
          this->polar_1_c_k_for_cols[n][k] = Polar(unity, -n * coeff * k) * S;
    }

    /* -------------------------- */
  }

  //% コンストラクタ
  Fourier2DConvolution() {
    this->initialize_polar_1_c_k();
    this->clear();
  }

  template <typename TYPE, size_t ROWS, size_t COLS> Fourier2DConvolution(const std::array<std::array<TYPE, COLS>, ROWS> &sample2DIN, const TREATED coeff, const bool reverse0 = false, const bool reverse1 = false) {
    this->initialize_polar_1_c_k();
    this->MatrixDFT = dft_fill(sample2DIN, reverse0, reverse1) * coeff;
  }

  template <typename TYPE, size_t ROWS, size_t COLS> Fourier2DConvolution(const std::array<std::array<TYPE, COLS>, ROWS> &sample2DIN, const bool reverse0 = false, const bool reverse1 = false) {
    this->initialize_polar_1_c_k();
    this->MatrixDFT = dft_fill(sample2DIN, reverse0, reverse1);
  }

  //% フーリエ変換のアダマール積を蓄積するためのメソッド
  template <typename TYPE, size_t ROWS, size_t COLS> void add(const std::array<std::array<TYPE, COLS>, ROWS> &sample2DIN) {
    static_assert(ROWS <= rows && COLS <= cols, "Input size exceeds the defined size of Fourier2DConvolution.");
    this->MatrixDFT = dft_fill(sample2DIN);
  }

  template <typename TYPE, size_t ROWS, size_t COLS> void add(const std::array<std::array<TYPE, COLS>, ROWS> &sample2DIN, const TREATED coeff) {
    static_assert(ROWS <= rows && COLS <= cols, "Input size exceeds the defined size of Fourier2DConvolution.");
    this->MatrixDFT = dft_fill(sample2DIN) * coeff;
  }

  /* -------------------------------------------- */

  template <typename ANYTYPE> void adamard_product_add(const Fourier2DConvolution<ANYTYPE, TREATED, Nk, Mk, Ns, Ms> &inputDFT2, const Fourier2DConvolution<ANYTYPE, TREATED, Nk, Mk, Ns, Ms> &inputDFT1) { this->MatrixDFT += inputDFT2.MatrixDFT * inputDFT1.MatrixDFT; };

  template <typename ANYTYPE> void adamard_product_add(const Fourier2DConvolution<ANYTYPE, TREATED, Nk, Mk, Ns, Ms> &inputDFT2, const Fourier2DConvolution<ANYTYPE, TREATED, Nk, Mk, Ns, Ms> &inputDFT1, const TREATED coeff) { this->MatrixDFT += inputDFT2.MatrixDFT * inputDFT1.MatrixDFT * coeff; };

  std::complex<TREATED> fma_increment(const std::complex<TREATED> &a, const std::complex<TREATED> &b, std::complex<TREATED> &c) {
    c.real(std::fma(a.real(), b.real(), std::fma(-a.imag(), b.imag(), c.real())));
    c.imag(std::fma(a.real(), b.imag(), std::fma(a.imag(), b.real(), c.imag())));
    return c;
  }

  template <typename TYPE, size_t N> std::array<std::complex<TREATED>, N> InverseDFT(const std::array<TYPE, N> &sample) {
    static constexpr auto polar_table = detail::make_inverse_dft_polar_table<N, __float128>();
    std::array<std::complex<TREATED>, N> result{};
    std::complex<TREATED> sum = 0;
    for (int n = 0; n < N; ++n) {
      sum = 0;
      const auto &polar_row = polar_table[n];
      for (int k = 0; k < N; ++k)
        sum += static_cast<std::complex<TREATED>>(sample[k]) * polar_row[k];
      // fma_increment(sample[k], polar_row[k], sum);
      result[n] = sum;
    }
    return result;
  }

  template <typename TYPE> std::array<std::array<std::complex<TREATED>, cols>, rows> InverseDFT(const std::array<std::array<TYPE, cols>, rows> &sample2D) {
    std::array<std::array<std::complex<TREATED>, cols>, rows> result;
    std::complex<TREATED> sum = 0;
    TREATED cx, cy;
    std::array<std::array<std::complex<TREATED>, cols>, rows> cn2D;
    // ! x方向のフーリエ変換
    for (int i = 0; i < rows; ++i)
      cn2D[i] = InverseDFT(sample2D[i]); //! 長さM

    // ! y方向のフーリエ変換 cn2D[;;][j]
    std::array<std::complex<TREATED>, rows> cn_j;
    for (int j = 0; j < cols; ++j) {
      for (int i = 0; i < rows; ++i)
        cn_j[i] = cn2D[i][j];  //! 長さN
      cn_j = InverseDFT(cn_j); //! 長さN
      for (int i = 0; i < rows; ++i)
        cn2D[i][j] = cn_j[i];
    }
    // ! cn2Dtmp[;;][j]をcn2D[;;][j]にコピー
    return cn2D;
  }

  template <typename TYPE, size_t COLS> std::array<std::complex<TREATED>, cols> dft1d_cols(const std::array<TYPE, COLS> &sample) {
    static_assert(COLS <= cols, "Input size exceeds the defined size of Fourier2DConvolution.");
    constexpr TREATED S = static_cast<TREATED>(cols);
    std::array<std::complex<TREATED>, cols> result;
    std::complex<TREATED> sum = 0;
    for (int n = 0; n < cols; ++n) {
      sum = 0;
      for (int k = 0; k < COLS; ++k)
        fma_increment(static_cast<std::complex<TREATED>>(sample[k]), this->polar_1_c_k_for_cols[n][k], sum);
      result[n] = sum;
    }
    return result;
  }

  template <typename TYPE, size_t ROWS> std::array<std::complex<TREATED>, rows> dft1d_rows(const std::array<TYPE, ROWS> &sample) {
    static_assert(ROWS <= rows, "Input size exceeds the defined size of Fourier2DConvolution.");
    constexpr TREATED S = static_cast<TREATED>(rows);
    std::array<std::complex<TREATED>, rows> result;
    std::complex<TREATED> sum = 0;
    for (int n = 0; n < rows; ++n) {
      sum = 0;
      for (int k = 0; k < ROWS; ++k)
        fma_increment(static_cast<std::complex<TREATED>>(sample[k]), this->polar_1_c_k_for_rows[n][k], sum);
      // sum += sample[k] * this->polar_1_c_k_for_rows[n][k];
      result[n] = sum;
    }
    return result;
  }

  template <typename TYPE, size_t ROWS, size_t COLS> std::array<std::array<std::complex<TREATED>, cols>, rows> dft_fill(const std::array<std::array<TYPE, COLS>, ROWS> &sample2DIN, bool reverse0 = false, bool reverse1 = false) {
    static_assert(ROWS <= rows && COLS <= cols, "Input size exceeds the defined size of Fourier2DConvolution.");
    std::array<std::array<std::complex<TREATED>, cols>, rows> cn2D, cn2Dtmp;
    std::array<std::complex<TREATED>, rows> cn_j_dft;
    std::array<std::complex<TREATED>, ROWS> cn_j;
    std::complex<TREATED> sum = 0;

    if (reverse0) {
      for (int i = 0; i < ROWS; ++i)
        cn2D[i] = dft1d_cols(sample2DIN[ROWS - 1 - i /*due to reverse0*/]); //! 長さM
    } else {
      for (int i = 0; i < ROWS; ++i)
        cn2D[i] = dft1d_cols(sample2DIN[i /*due to reverse0*/]); //! 長さM
    }
    // ! y方向のフーリエ変換 cn2D[;;][j]
    for (int j = 0; j < cols; ++j) {
      for (int i = 0; i < ROWS; ++i)
        cn_j[i] = cn2D[i][j];      //! 長さN
      cn_j_dft = dft1d_rows(cn_j); //! 長さN
      for (int i = 0; i < rows; ++i)
        cn2Dtmp[i][j] = cn_j_dft[i];
    }
    return cn2Dtmp;
  };
};

/* ============================================================================ */

template <typename T> std::vector<std::vector<std::complex<double>>> InverseDFT(const Fourier2D<T> &convolver2Dkernel) { return InverseDFT(convolver2Dkernel.MatrixDFT); }

inline std::vector<std::vector<double>> Re(const std::vector<std::vector<std::complex<double>>> &matrix2d) {
  std::vector<std::vector<double>> Re_matrix2d(matrix2d.size(), std::vector<double>(matrix2d[0].size(), 0));
  for (std::size_t i = 0; i < matrix2d.size(); ++i)
    for (std::size_t j = 0; j < matrix2d[0].size(); ++j)
      Re_matrix2d[i][j] = matrix2d[i][j].real();
  return Re_matrix2d;
}

template <typename T> std::vector<std::vector<std::complex<double>>> Convolution(const Fourier2D<T> &convolver2Dkernel, const Fourier2D<T> &convolver2Dshift) { return InverseDFT(convolver2Dkernel.MatrixDFT * convolver2Dshift.MatrixDFT); }

template <typename T> std::vector<std::vector<double>> ReConvolution(const Fourier2D<T> &convolver2Dkernel, const Fourier2D<T> &convolver2Dshift) {
  //! kernelDFTとshiftDFTを掛け算して，逆DFTをとる
  auto convolution_result = Convolution(convolver2Dkernel, convolver2Dshift);
  std::vector<std::vector<double>> ret(convolver2Dkernel.rows, std::vector<double>(convolver2Dkernel.cols, 0));
  for (int i = 0; i < convolver2Dkernel.rows; ++i)
    for (int j = 0; j < convolver2Dkernel.cols; ++j)
      ret[i][j] = convolution_result[i][j].real();
  return ret;
}

// \label{Convolver2D}
template <typename T> struct Convolver2D {
  std::size_t rows, cols; //! kernel size + shift - 1
  std::vector<std::vector<std::complex<double>>> kernelDFT;
  std::vector<std::vector<std::complex<double>>> shiftDFT;
  std::vector<std::vector<std::complex<double>>> convolution;

  Convolver2D() {};

  int Nk, Mk, Ns, Ms; //! indexの最大
  void reset(int size_Nk, int size_Mk, int size_Ns, int size_Ms) {
    Nk = size_Nk - 1; //! indexの最大
    Mk = size_Mk - 1; //! indexの最大
    Ns = size_Ns - 1; //! indexの最大
    Ms = size_Ms - 1; //! indexの最大
    this->rows = Nk + Ns + 1;
    this->cols = Mk + Ms + 1; //! 1x1 1x1なら size 1必要ない
    // 複素数行列をゼロ初期化
    this->kernelDFT = std::vector<std::vector<std::complex<double>>>(rows, std::vector<std::complex<double>>(cols, std::complex<double>(0.0, 0.0)));
    this->shiftDFT = std::vector<std::vector<std::complex<double>>>(rows, std::vector<std::complex<double>>(cols, std::complex<double>(0.0, 0.0)));
    this->kernelpadded = std::vector<std::vector<double>>(this->rows, std::vector<T>(this->cols, 0));
    this->shiftpadded = std::vector<std::vector<double>>(this->rows, std::vector<T>(this->cols, 0));
  };

  bool corr_row = false;
  bool corr_col = false;

  std::vector<std::vector<double>> kernelpadded;
  void addKernel(const std::vector<std::vector<T>> &kernel) {
    if (kernel.size() > this->rows || kernel[0].size() > this->cols) {
      std::cerr << "Error: Kernel size does not match." << std::endl;
      return;
    } else {
      // this->kernelpadded = std::vector<std::vector<double>>(this->rows, std::vector<T>(this->cols, 0));

      for (int i = 0; i < kernel.size(); ++i)
        for (int j = 0; j < kernel[0].size(); ++j)
          kernelpadded[i][j] = kernel[i][j];

      auto kernelpaddedDFT = DFT(kernelpadded);
      for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
          kernelDFT[i][j] += kernelpaddedDFT[i][j];
    }
  };

  std::vector<std::vector<double>> shiftpadded;
  void addShift(std::vector<std::vector<T>> shift, std::array<bool, 2> reverse = {false, false}) {
    if (shift.size() > this->rows || shift[0].size() > this->cols) {
      std::cerr << "Error: Shift size does not match." << std::endl;
      return;
    } else {
      int si = shift.size();
      int sj = shift[0].size();
      // this->shiftpadded = std::vector<std::vector<double>>(this->rows, std::vector<T>(this->cols, 0));

      for (int i = 0; i < shift.size(); ++i)
        for (int j = 0; j < shift[0].size(); ++j)
          shiftpadded[i][j] = shift[reverse[0] ? (si - 1 - i) : i][reverse[1] ? (sj - 1 - j) : j];

      auto shiftpaddedDFT = DFT(shiftpadded);
      for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
          shiftDFT[i][j] += shiftpaddedDFT[i][j];
    }
  };

  void convolve() { convolution = InverseDFT(kernelDFT * shiftDFT); };

  std::vector<std::vector<double>> getReConvolution() const {
    std::vector<std::vector<double>> ReConvolution(rows, std::vector<double>(cols, 0));
    for (int i = 0; i < rows; ++i)
      for (int j = 0; j < cols; ++j)
        ReConvolution[i][j] = convolution[i][j].real();
    return ReConvolution;
  };
};

#endif