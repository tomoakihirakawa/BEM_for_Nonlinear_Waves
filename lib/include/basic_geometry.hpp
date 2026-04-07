#ifndef basic_geometry_H
#define basic_geometry_H

#include "basic_linear_systems.hpp"
#include "basic_statistics.hpp"
#include "basic_vectors.hpp"
#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"
#include <execution>
// 面と面の干渉？？？
// 球と面の干渉チェックか．
// そのために，
// タプルを作ろうとしている，

// using Tdd = std::tuple<double, double>;
// using Tddd = std::tuple<double, double, double>;
// using T2Tddd = std::tuple<Tddd, Tddd>;
// using T3Tddd = std::tuple<Tddd, Tddd, Tddd>;
// using T3Tdd = std::tuple<Tdd, Tdd, Tdd>;

/* -------------------------------------------------------------------------- */
inline double SolidAngle_VanOosteromAandStrackeeJ1983(const Tddd& p, Tddd A, Tddd B, Tddd C) {
  // A -= p;
  // B -= p;
  // C -= p;

  // const double a = Norm(A);
  // const double b = Norm(B);
  // const double c = Norm(C);

  // const auto crossBC = Cross(B, C);
  // const double numerator = std::abs(Dot(A, crossBC));

  // const double denominator = a * b * c +
  //                            Dot(A, B) * c +
  //                            Dot(A, C) * b +
  //                            Dot(B, C) * a;

  // return 2.0 * std::atan2(numerator, denominator);
  /* -------------------------------------------------------------- */
  // The solid angle of a plane triangle
  // Van Oosterom, A. and Strackee, J. (1983)
  const auto [a0, a1, a2] = (A -= p);
  const auto [b0, b1, b2] = (B -= p);
  const auto [c0, c1, c2] = (C -= p);
  const double nA = Norm(A);
  const double nB = Norm(B);

  // double nC = Norm(C);
  // return 2. * std::atan2((-(a2 * b1 * c0) + a1 * b2 * c0 + a2 * b0 * c1 - a0 * b2
  // * c1 - a1 * b0 * c2 + a0 * b1 * c2), 				  (nB * (a0 * c0 + a1 * c1 + a2 * c2) +
  // nA * (b0 * c0 + b1 * c1 + b2 * c2) + (a0 * b0 + a1 * b1 + a2 * b2) * nC +
  // nA * nB * nC));
  return 2. * std::atan2(-(a2 * b1 * c0) + a1 * b2 * c0 + a2 * b0 * c1 - a0 * b2 * c1 - a1 * b0 * c2 + a0 * b1 * c2, b0 * c0 * nA + b1 * c1 * nA + b2 * c2 * nA + a0 * c0 * nB + a1 * c1 * nB + a2 * c2 * nB + (a0 * b0 + a1 * b1 + a2 * b2 + nA * nB) * Norm(C));
};

inline double SolidAngle_VanOosteromAandStrackeeJ1983(const Tddd& p, const T3Tddd& ABC) { return SolidAngle_VanOosteromAandStrackeeJ1983(p, std::get<0>(ABC), std::get<1>(ABC), std::get<2>(ABC)); };
inline double SolidAngle_VanOosteromAandStrackeeJ1983(const Tddd& p, const Tddd& A) { return 0.; };
//
// double SolidAngle(const Tddd &o, const Tddd &A, const Tddd &B, const Tddd &C) {
//    double c = VectorAngle(A - o, B - o), a = VectorAngle(B - o, C - o), b = VectorAngle(C - o, A - o);
//    double s = (a + b + c) * 0.5;
//    if (Between(s, {M_PI - 1E-10, M_PI + 1E-10}))
//       return 4. * M_PI / 2.;
//    else
//       return 4. * std::std::atan(std::sqrt(tan(s * 0.5) * std::tan((s - a) * 0.5) * std::tan((s - b) * 0.5) * std::tan((s - c) * 0.5)));
//    // return SolidAngle_VanOosteromAandStrackeeJ1983(o, A, B, C);
// };

inline double SolidAngle_UsingVectorAngle(const Tddd& p0, const Tddd& p1, const Tddd& p2, const Tddd& p3) {
  double c = VectorAngle(p1 - p0, p2 - p0), a = VectorAngle(p2 - p0, p3 - p0), b = VectorAngle(p3 - p0, p1 - p0);
  double s = (a + b + c) * 0.5;
  if (Between(s, {M_PI - 1E-10, M_PI + 1E-10}))
    return 4. * M_PI / 2.;
  else
    return 4. * std::atan(std::sqrt(tan(s * 0.5) * std::tan((s - a) * 0.5) * std::tan((s - b) * 0.5) * std::tan((s - c) * 0.5)));
};

inline double SolidAngle_(const Tddd& p0, const Tddd& p1, const Tddd& p2, const Tddd& p3) {

  //     1,2,3 (outward rotation)
  //          3
  // 0,3,2  / | \ 0,1,3     --- 1 ---
  //       2--|--1           \      /
  // 0,2,1  \ | /             2    0
  //          0                 \/
  //
  /* -------------------------------------------------------------------------- */
  // return SolidAngle_VanOosteromAandStrackeeJ1983(p0, p1, p2, p3);
  /* -------------------------------------------------------------------------- */
  // auto A = Normalize(p1 - p0);
  // auto B = Normalize(p2 - p0);
  // auto C = Normalize(p3 - p0);
  // return 2 * atan2(Dot(A, Cross(B, C)), (1 + Dot(A, B) + Dot(B, C) + Dot(C, A)));
  auto A = (p1 - p0);
  auto B = (p2 - p0);
  auto C = (p3 - p0);
  return 2. * std::atan2(Det(T3Tddd{A, B, C}) /*Dot(A, Cross(B, C))*/, Norm(A) * Norm(B) * Norm(C) + Dot(A, B) * Norm(C) + Dot(A, C) * Norm(B) + Dot(B, C) * Norm(A));
};

inline double SolidAngle(const Tddd& o, const Tddd& A, const Tddd& B, const Tddd& C) {
  // return SolidAngle_(o, a, b, c);
  //
  // return std::abs(SolidAngle_VanOosteromAandStrackeeJ1983(o, a, b, c));
  //
  // auto v1 = A - o;
  // auto v2 = B - o;
  // auto v3 = C - o;
  // auto cv1v2 = Normalize(Cross(v1, v2));
  // auto cv1v3 = Normalize(Cross(v1, v3));
  // auto cv2v3 = Normalize(Cross(v2, v3));
  // return std::abs(acos(Dot(cv1v2, cv1v3)) + acos(-Dot(cv1v2, cv2v3)) + acos(Dot(cv1v3, cv2v3)) - M_PI);

  double c = VectorAngle(A - o, B - o), a = VectorAngle(B - o, C - o), b = VectorAngle(C - o, A - o);
  double s = (a + b + c) * 0.5;
  if (Between(s, {M_PI - 1E-10, M_PI + 1E-10}))
    return 4. * M_PI / 2.;
  else
    return 4. * std::atan(std::sqrt(std::tan(s * 0.5) * std::tan((s - a) * 0.5) * std::tan((s - b) * 0.5) * std::tan((s - c) * 0.5)));

  // return SolidAngle_VanOosteromAandStrackeeJ1983(o, A, B, C);
};

inline double SolidAngle(const Tddd& p, const T3Tddd& ABC) {
  return SolidAngle(p, std::get<0>(ABC), std::get<1>(ABC), std::get<2>(ABC));
  // return SolidAngle_VanOosteromAandStrackeeJ1983(p, std::get<0>(ABC),
  // std::get<1>(ABC), std::get<2>(ABC));
};

inline T4d SolidAngles(const Tddd& o, const Tddd& a, const Tddd& b, const Tddd& c) {
  return {(SolidAngle_VanOosteromAandStrackeeJ1983(o, a, b, c)), (SolidAngle_VanOosteromAandStrackeeJ1983(a, b, c, o)), (SolidAngle_VanOosteromAandStrackeeJ1983(b, c, o, a)), (SolidAngle_VanOosteromAandStrackeeJ1983(c, o, a, b))};
  //
  // return {std::abs(SolidAngle_VanOosteromAandStrackeeJ1983(o, a, b, c)),
  //         std::abs(SolidAngle_VanOosteromAandStrackeeJ1983(a, b, c, o)),
  //         std::abs(SolidAngle_VanOosteromAandStrackeeJ1983(b, c, o, a)),
  //         std::abs(SolidAngle_VanOosteromAandStrackeeJ1983(c, o, a, b))};
};

inline T4d SolidAngles(const T4Tddd& oabc) {
  auto [o, a, b, c] = oabc;
  return SolidAngles(o, a, b, c);
};

inline double SolidAngle(const Tddd& o, const std::vector<Tddd>& xyz) {
  double total = 0., tmp, angle;
  int sz = xyz.size();
  double local_mean_length = 0;
  std::array<double, 3> local_mean_normal1 = {0., 0., 0.}, local_mean_normal2 = {0., 0., 0.};
  Tddd cross0, cross1;
  bool isAllFlat = true;
  for (auto i = 0; i < sz; ++i) {
    local_mean_length += Norm(xyz[i] - o);
    cross1 = Cross(xyz[i] - o, xyz[(i + 1) == sz ? 0 : (i + 1)] - o);
    local_mean_normal1 += -Normalize(cross0); /*内向き*/
    if (i != 0) {
      if (!isFlat(cross0, cross1, 1E-5))
        isAllFlat = false;
    }
    local_mean_normal2 += Normalize(xyz[i] - o);
    cross0 = cross1;
  }
  local_mean_length /= sz;
  local_mean_normal1 /= sz;
  local_mean_normal2 /= sz;
  Tddd n = isAllFlat ? Normalize(local_mean_normal1) : Normalize(local_mean_normal2);
  Tddd A, B, n0, n1;
  std::array<double, 3> point_at_outward_direction = local_mean_length * n + o;
  for (auto i = 0; i < xyz.size(); ++i) {
    A = xyz[i];
    B = xyz[(i + 1) == sz ? 0 : (i + 1)];
    tmp = SolidAngle(o, A, B, point_at_outward_direction);
    // n0 = Cross(A - o, B - o);
    // n1 = Cross(B - o, point_at_outward_direction - o);
    // if (!isFlat(n0, n1, 1E-10) && !isFlat(n0, -n1, 1E-10))
    total += tmp;
  }
  return total;
};

inline T4d TetrahedronSolidAngle_UsingVectorAngle(const Tddd& X0, const Tddd& X1, Tddd X2, Tddd X3) {
  if (Dot(TriangleNormal(X1, X2, X3), (X1 + X2 + X3) / 3. - X0) < 0)
    X2.swap(X3);
  return {(SolidAngle_UsingVectorAngle(X0, X1, X2, X3)), (SolidAngle_UsingVectorAngle(X1, X0, X3, X2)), (SolidAngle_UsingVectorAngle(X2, X0, X1, X3)), (SolidAngle_UsingVectorAngle(X3, X0, X2, X1))};
};

inline T4d TetrahedronSolidAngle(const Tddd& X0, const Tddd& X1, Tddd X2, Tddd X3) {
  if (Dot(TriangleNormal(X1, X2, X3), (X1 + X2 + X3) / 3. - X0) < 0)
    X2.swap(X3);
  return {(SolidAngle_VanOosteromAandStrackeeJ1983(X0, X1, X2, X3)), (SolidAngle_VanOosteromAandStrackeeJ1983(X1, X0, X3, X2)), (SolidAngle_VanOosteromAandStrackeeJ1983(X2, X0, X1, X3)), (SolidAngle_VanOosteromAandStrackeeJ1983(X3, X0, X2, X1))};
};

inline T4d TetrahedronSolidAngle(const T4Tddd& abcd) { return TetrahedronSolidAngle(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd)); };

inline T4d TetrahedronSolidAngle_UsingVectorAngle(const T4Tddd& abcd) { return TetrahedronSolidAngle_UsingVectorAngle(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd)); };
// % -------------------------------------------------------------------------- */
// %                                  外接球                                     */
// % -------------------------------------------------------------------------- */
inline Tddd Circumcenter(const Tddd& a, const Tddd& b_, const Tddd& c_) {
  Tddd b = b_ - a;
  Tddd c = c_ - a;
  const auto x = Normalize(b);
  const auto z = Normalize(Cross(b, c));
  const auto y = Normalize(Cross(x, z)); // 1x1
  //
  Tdd bxy = {Dot(b, x), Dot(b, y)};
  Tdd cxy = {Dot(c, x), Dot(c, y)};
  // auto [Cx, Cy] = Dot(Inverse(T2Tdd{bxy, cxy}), 0.5 * Tdd{Dot(bxy, bxy), Dot(cxy, cxy)});  // circumcenter
  //! use solve instead of inverse
  Tdd CxCy = {0., 0.};
  // Solve(T2Tdd{bxy, cxy}, CxCy, Tdd{0.5 * Dot(bxy, bxy), 0.5 * Dot(cxy, cxy)});
  lapack_lu lu(T2Tdd{bxy, cxy}, CxCy, Tdd{0.5 * Dot(bxy, bxy), 0.5 * Dot(cxy, cxy)});
  // lapack_svd_solve(T2Tdd{bxy, cxy}, CxCy, Tdd{0.5 * Dot(bxy, bxy), 0.5 * Dot(cxy, cxy)});
  auto [Cx, Cy] = CxCy;
  return a + Cx * x + Cy * y;
};
inline Tddd Circumcenter(const T3Tddd& abcd) { return Circumcenter(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd)); };
inline Tddd Circumcenter(const Tddd& a, const Tddd& b, const Tddd& c, const Tddd& d) {
  // http://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter#:~:text=For%20all%20tetrahedra%2C%20there%20exists,circumsphere%20is%20called%20the%20circumcentre.
  // double a2 = Dot(a, a);
  // return 0.5 * Dot(Inverse(T3Tddd{b - a, c - a, d - a}), Tddd{Dot(b, b) - a2, Dot(c, c) - a2, Dot(d, d) - a2});
  std::array<double, 3> b_ = b - a, c_ = c - a, d_ = d - a, CxCyCz;
  // return a + 0.5 * Dot(Inverse(T3Tddd{b_, c_, d_}), Tddd{Dot(b_, b_), Dot(c_, c_), Dot(d_, d_)});
  //! use solve
  // Solve(T3Tddd{b_, c_, d_}, CxCyCz, Tddd{Dot(b_, b_), Dot(c_, c_), Dot(d_, d_)});
  lapack_lu lu(T3Tddd{b_, c_, d_}, CxCyCz, Tddd{Dot(b_, b_), Dot(c_, c_), Dot(d_, d_)});
  // lapack_svd_solve(T3Tddd{b_, c_, d_}, CxCyCz, Tddd{Dot(b_, b_), Dot(c_, c_), Dot(d_, d_)});
  return a + 0.5 * CxCyCz;
};
inline Tddd Circumcenter(const T4Tddd& abcd) { return Circumcenter(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd)); };
/* -------------------------------------------------------------------------- */
inline double Circumradius(const Tddd& a, const Tddd& b, const Tddd& c) {
  // auto X = Circumcenter(a, b, c);
  // return (Norm(a - X) + Norm(b - X) + Norm(c - X)) / 3.;
  const double ab = Norm(a - b);
  const double bc = Norm(b - c);
  const double ca = Norm(c - a);
  const double s = (ab + bc + ca) * 0.5;
  return ab * bc * ca / (4. * std::sqrt(s * (s - ab) * (s - bc) * (s - ca)));
};
inline double Circumradius(const T3Tddd& abcd) { return Circumradius(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd)); };
inline double CircumArea(const T3Tddd& abcd) { return std::pow(Circumradius(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd)), 2) * M_PI; };
inline double Circumradius(const Tddd& a, const Tddd& b, const Tddd& c, const Tddd& d) {
  auto X = Circumcenter(a, b, c, d);
  return (Norm(a - X) + Norm(b - X) + Norm(c - X) + Norm(d - X)) / 4.;
};
inline double Circumradius(const T4Tddd& abcd) { return Circumradius(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd)); };
// % -------------------------------------------------------------------------- */
// %                                  内接球                                     */
// % -------------------------------------------------------------------------- */
inline double Inradius(Tddd p0, Tddd p1, Tddd p2) {
  p2 -= p0;
  p1 -= p0;
  p0.fill(0.);
  auto l2 = Norm(p2 - (-p1 * (Dot(p1 - p2, p1) / Dot(p1, p1)) + p1));
  return (Norm(p1) * l2) / (Norm(p1) + Norm(p2) + Norm(p1 - p2));
};
inline double Inradius(const T3Tddd& p0123) { return Inradius(std::get<0>(p0123), std::get<1>(p0123), std::get<2>(p0123)); };
inline double InArea(const T3Tddd& p0123) { return std::pow(Inradius(std::get<0>(p0123), std::get<1>(p0123), std::get<2>(p0123)), 2) * M_PI; };
inline double Inradius(Tddd p0, Tddd p1, Tddd p2, Tddd p3) {
  // see /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/四面体の内接球外接球.nb
  p3 -= p0;
  p2 -= p0;
  p1 -= p0;
  p0 = {0, 0, 0};
  double A0 = TriangleArea(p1, p2, p3);
  double A1 = TriangleArea(p0, p2, p3);
  double A2 = TriangleArea(p0, p1, p3);
  double A3 = TriangleArea(p0, p1, p2);
  return 0.5 * CrossDot(p1, p2, p3) / (A0 + A1 + A2 + A3);
  //
  // // see /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/四面体の内接球外接球.nb
  // double l3 = Norm(Dot(p3 - p0, Normalize(Cross(p1 - p0, p2 - p0))));
  // double A0 = TriangleArea(p1, p2, p3);
  // double A1 = TriangleArea(p0, p2, p3);
  // double A2 = TriangleArea(p0, p1, p3);
  // double A3 = TriangleArea(p0, p1, p2);
  // return (A3 * l3) / (A0 + A1 + A2 + A3);
};
inline double Inradius(const T4Tddd& p0123) { return Inradius(std::get<0>(p0123), std::get<1>(p0123), std::get<2>(p0123), std::get<3>(p0123)); };
/* -------------------------------------------------------------------------- */
inline Tddd Incenter(const Tddd& p0, const Tddd& p1, const Tddd& p2) {
  // see /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/四面体の内接球外接球.nb
  // https://en.wikipedia.org/wiki/Tetrahedron
  auto len = [](const Tddd& a, const Tddd& b, const Tddd& c) {
    const auto a_b = a - b;
    return Norm(c - ((a_b)*Dot(c - b, a_b) / Dot(a_b, a_b) + b));
  };
  auto l0 = len(p1, p2, p0);
  auto l1 = len(p2, p0, p1);
  auto l2 = len(p0, p1, p2);
  return (p0 / l0 + p1 / l1 + p2 / l2) * Inradius(p0, p1, p2);
};
inline Tddd Incenter(const T3Tddd& p0123) { return Incenter(std::get<0>(p0123), std::get<1>(p0123), std::get<2>(p0123)); };
inline Tddd Incenter(const Tddd& p0, const Tddd& p1, const Tddd& p2, const Tddd& p3) {
  // see /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/四面体の内接球外接球.nb
  // https://en.wikipedia.org/wiki/Tetrahedron
  double l0 = Norm(Dot(p0 - p1, Normalize(Cross(p2 - p1, p3 - p1))));
  double l1 = Norm(Dot(p1 - p2, Normalize(Cross(p3 - p2, p0 - p2))));
  double l2 = Norm(Dot(p2 - p3, Normalize(Cross(p0 - p3, p1 - p3))));
  double l3 = Norm(Dot(p3 - p0, Normalize(Cross(p1 - p0, p2 - p0))));
  return (p0 / l0 + p1 / l1 + p2 / l2 + p3 / l3) * Inradius(p0, p1, p2, p3);
};
inline Tddd Incenter(const T4Tddd& p0123) { return Incenter(std::get<0>(p0123), std::get<1>(p0123), std::get<2>(p0123), std::get<3>(p0123)); };
// % -------------------------------------------------------------------------- */
inline Tddd Centroid(const Tddd& a, const Tddd& b, const Tddd& c) { return (a + b + c) / 3.; };
inline Tddd Centroid(const T3Tddd& abcd) { return Centroid(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd)); };
inline Tddd Centroid(const Tddd& a, const Tddd& b, const Tddd& c, const Tddd& d) { return (a + b + c + d) / 4.; };
inline Tddd Centroid(const T4Tddd& abcd) { return Centroid(std::get<0>(abcd), std::get<1>(abcd), std::get<2>(abcd), std::get<3>(abcd)); };
// % -------------------------------------------------------------------------- */
inline bool InsideQ(const T3Tdd& bounds, const T3Tdd& bounds_) {
  //! bounds < bounds_ ?
  //! Is bounds inside bounds_ ?
  auto [X0, X1] = std::get<0>(bounds);
  auto [Y0, Y1] = std::get<1>(bounds);
  auto [Z0, Z1] = std::get<2>(bounds);
  auto [X0_, X1_] = std::get<0>(bounds_);
  auto [Y0_, Y1_] = std::get<1>(bounds_);
  auto [Z0_, Z1_] = std::get<2>(bounds_);
  return {X0_ < X0 && X1 < X1_ && Y0_ < Y0 && Y1 < Y1_ && Z0_ < Z0 && Z1 < Z1_};
};

inline bool InsideQ(const Tddd& X, const Tddd& Xcenter, const double& r) {
  //! point v.s. sphere
  //! point < sphere ?
  //! Is point X inside the sphere with center Xcenter and radius r ?
  return Norm(X - Xcenter) < r;
};
inline bool InsideQ(const T3Tdd& bounds, const Tddd& Xcenter, const double& r) {
  //! cube v.s. sphere
  //! cube < sphere ?
  //! Is cube bounds inside the sphere with center Xcenter and radius r ?
  auto [X0, X1] = std::get<0>(bounds);
  auto [Y0, Y1] = std::get<1>(bounds);
  auto [Z0, Z1] = std::get<2>(bounds);
  auto isInside = [&](const Tddd& X, const Tddd& Xcenter, const double& r) {
    // point v.s. sphere
    return Norm(X - Xcenter) < r;
  };
  return (isInside({X0, Y0, Z0}, Xcenter, r) && isInside({X1, Y0, Z0}, Xcenter, r) && isInside({X1, Y1, Z0}, Xcenter, r) && isInside({X0, Y1, Z0}, Xcenter, r) && isInside({X0, Y0, Z1}, Xcenter, r) && isInside({X1, Y0, Z1}, Xcenter, r) && isInside({X0, Y1, Z1}, Xcenter, r) && isInside({X1, Y1, Z1}, Xcenter, r));
};
inline bool InsideQ(const Tddd& Xcenter, const double& r, const T3Tdd& bounds) {
  auto [X0, X1] = std::get<0>(bounds);
  auto [Y0, Y1] = std::get<1>(bounds);
  auto [Z0, Z1] = std::get<2>(bounds);
  auto isInside = [&](const Tddd& X, const Tddd& Xcenter, const double& r) {
    // point v.s. sphere
    return Norm(X - Xcenter) < r;
  };
  return (!isInside({X0, Y0, Z0}, Xcenter, r) && !isInside({X1, Y0, Z0}, Xcenter, r) && !isInside({X1, Y1, Z0}, Xcenter, r) && !isInside({X0, Y1, Z0}, Xcenter, r) && !isInside({X0, Y0, Z1}, Xcenter, r) && !isInside({X1, Y0, Z1}, Xcenter, r) && !isInside({X0, Y1, Z1}, Xcenter, r) && !isInside({X1, Y1, Z1}, Xcenter, r));
};

/* -------------------------------------------------------------------------- */

inline Tddd grad_CircumradiusToInradius(const Tddd& a, const Tddd& b, const Tddd& c) {
  auto grad = [](const Tddd& a, const Tddd& b) { return (a - b) / Norm(a - b); };
  const double ab = Norm(a - b);
  const double bc = Norm(b - c);
  const double ac = Norm(a - c);
  const double den1 = -ab + bc + ac;
  const double den2 = ab - bc + ac;
  const double den3 = ab + bc - ac;
  const double denominator = den1 * den2 * den3;
  const std::array<double, 3> grad1 = -grad(a, b) + grad(a, c); // derivative of (-ab + bc + ca)
  const std::array<double, 3> grad2 = grad(a, b) + grad(a, c);  // derivative of (ab - bc + ca)
  const std::array<double, 3> grad3 = grad(a, b) - grad(a, c);  // derivative of (ab + bc - ca)
  const std::array<double, 3> grad_denominator_a = grad1 * den2 * den3 + den1 * grad2 * den3 + den1 * den2 * grad3;
  //
  const double numerator = (ab * bc * ac);
  const std::array<double, 3> grad_numerator_a = grad(a, b) * bc * ac + ab * bc * grad(a, c);
  //
  return 2. * (grad_numerator_a * denominator - numerator * grad_denominator_a) / std::pow(denominator, 2);
};

inline double differenceFromEquilateralTriangle(const Tddd& a, const Tddd& b, const Tddd& c) {
  auto weight = [&](const Tddd& a, const Tddd& b, const Tddd& c) -> double {
    Tddd Xmid = (c + b) * 0.5;
    double height = Norm(c - b) * std::sqrt(3.) * 0.5;
    Tddd vertical = Normalize(Chop(a - Xmid, c - b));
    return Norm(height * vertical + Xmid - a) / height;
  };
  return (weight(a, b, c) + weight(b, c, a) + weight(c, a, b)) / 3.;
};

// 四面体の R/r を計算
inline double CircumradiusToInradius(const std::array<double, 3>& A, const std::array<double, 3>& B, const std::array<double, 3>& C, const std::array<double, 3>& D) {
  // 体積と表面積
  double V = TetVolume(A, B, C, D);
  double S = TriangleArea(B, C, D) + TriangleArea(A, C, D) + TriangleArea(A, B, D) + TriangleArea(A, B, C);
  if (V < 1e-15 || S < 1e-15)
    return std::numeric_limits<double>::infinity();
  double r = 3.0 * V / S;

  // 外心の計算
  std::array<double, 3> u = B - A;
  std::array<double, 3> v = C - A;
  std::array<double, 3> w = D - A;

  std::array<std::array<double, 3>, 3> M = {{{2 * u[0], 2 * u[1], 2 * u[2]}, {2 * v[0], 2 * v[1], 2 * v[2]}, {2 * w[0], 2 * w[1], 2 * w[2]}}};
  std::array<double, 3> rhs = {Dot(u, u), Dot(v, v), Dot(w, w)};
  std::array<double, 3> x;

  lapack_lu lu(M, x, rhs);
  // lapack_svd_solve(M, x, rhs);

  std::array<double, 3> O = {A[0] + x[0], A[1] + x[1], A[2] + x[2]};
  double R = Norm(O - A);

  return R / r;
}

// triangle distorsion measure
inline double CircumradiusToInradius(const Tddd& a, const Tddd& b, const Tddd& c) {
  const double ab = Norm(a - b);
  const double bc = Norm(b - c);
  const double ca = Norm(c - a);
  const double den1 = (-ab + bc + ca);
  const double den2 = (ab - bc + ca);
  const double den3 = (ab + bc - ca);
  return 2. * ab * bc * ca / (den1 * den2 * den3);
  //
  // const double s = (ab + bc + ca) * 0.5;
  // const double area = std::sqrt(s * (s - ab) * (s - bc) * (s - ca));
  // const double inradius = area / s;
  // const double circumradius = ab * bc * ca / (4. * area);
  // return ab * bc * ca / (4. * ((s - ab) * (s - bc) * (s - ca)));
  /* ------------------------------- OLD version ------------------------------ */
  // return Circumradius(X0, X1, X2) / Inradius(X0, X1, X2);
  /* -------------------------------------------------------------------------- */
  // X2 -= X0;
  // X1 -= X0;
  // X0.fill(0.);
  // auto l2 = Norm(X2 - ((-X1) * Dot(X2 - X1, -X1) / Dot(-X1, -X1) + X1));
  // auto X = Circumcenter(X0, X1, X2);
  // return ((Norm(X0 - X) + Norm(X1 - X) + Norm(X2 - X))) * (Norm(X1) + Norm(X2) + Norm(X1 - X2)) / (3. * Norm(X1) * l2);
};

inline double CircumradiusToInradius(const T3Tddd& X012) {
  // return Circumradius(X012) / Inradius(X012);
  return CircumradiusToInradius(std::get<0>(X012), std::get<1>(X012), std::get<2>(X012));
};

// InradiusToCircumradius = 2 * r_in / r_out
// 正三角形で 1.0、退化で 0.0。有界 [0, 1] で品質スコアとして安全。
inline double InradiusToCircumradius(const Tddd& a, const Tddd& b, const Tddd& c) {
  const double ab = Norm(a - b);
  const double bc = Norm(b - c);
  const double ca = Norm(c - a);
  const double den1 = (-ab + bc + ca);
  const double den2 = (ab - bc + ca);
  const double den3 = (ab + bc - ca);
  return (den1 * den2 * den3) / (2. * ab * bc * ca);
}

inline double InradiusToCircumradius(const T3Tddd& X012) {
  return InradiusToCircumradius(std::get<0>(X012), std::get<1>(X012), std::get<2>(X012));
}

inline double log10_CircumradiusToInradius(Tddd X0, Tddd X1, Tddd X2) { return std::log10(Circumradius(X0, X1, X2)) - std::log10(Inradius(X0, X1, X2)); };

/* -------------------------------------------------------------------------- */
/*
M. Meyer, M. Desbrun, P. Schröder, and A. H. Barr, “Discrete
Differential-Geometry Operators for Triangulated 2-Manifolds BT  - Visualization
and Mathematics III,” Vis. Math. III, pp. 35–57, 2003.
*/

inline Tddd ToSphericalCoodrinates(const Tddd& xyz) {
  double r = Norm(xyz);
  return {r, std::atan(std::get<2>(xyz) / r), std::atan2(std::get<1>(xyz), std::get<0>(xyz))};
};

/* -------------------------------------------------------------------------- */

struct Point {
  Tddd X;
  Point(const Tddd& XIN) : X(XIN) {};
};

/* -------------------------------------------------------------------------- */

template <typename T>
struct Edge {
  Edge(T a, T b);
};

/* -------------------------------------------------------------------------- */

inline T3Tdd Distance(const T3Tdd& b, const Tddd& a) {
  auto [mmX, mmY, mmZ] = b;
  T3Tdd ret = b;
  //
  Tdd to_ax = std::get<0>(a) - mmX;
  Tdd to_ay = std::get<1>(a) - mmY;
  Tdd to_az = std::get<2>(a) - mmZ;
  //
  std::get<0>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ax)), std::abs(std::get<1>(to_ax))}); //<-Tdd
  std::get<1>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ay)), std::abs(std::get<1>(to_ay))}); //<-Tdd
  std::get<2>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_az)), std::abs(std::get<1>(to_az))}); //<-Tdd
  //
  if (Between(std::get<0>(a), mmX))
    std::get<0>(std::get<0>(ret)) = 0;
  if (Between(std::get<1>(a), mmY))
    std::get<0>(std::get<1>(ret)) = 0;
  if (Between(std::get<2>(a), mmZ))
    std::get<0>(std::get<2>(ret)) = 0;
  return ret;
};

/* -------------------------------------------------------------------------- */

// namespace geometry {
// struct Line {
//    T2Tddd X;
//    double minX, maxX, minY, maxY, minZ, maxZ;
//    T3Tdd bounds;
//    Line(const T2Tddd &XIN)
//        : X(XIN),
//          minX(Min(Tdd{std::get<0>(std::get<0>(XIN)),
//                       std::get<0>(std::get<1>(XIN))})),
//          maxX(Max(Tdd{std::get<0>(std::get<0>(XIN)),
//                       std::get<0>(std::get<1>(XIN))})),
//          minY(Min(Tdd{std::get<1>(std::get<0>(XIN)),
//                       std::get<1>(std::get<1>(XIN))})),
//          maxY(Max(Tdd{std::get<1>(std::get<0>(XIN)),
//                       std::get<1>(std::get<1>(XIN))})),
//          minZ(Min(Tdd{std::get<2>(std::get<0>(XIN)),
//                       std::get<2>(std::get<1>(XIN))})),
//          maxZ(Max(Tdd{std::get<2>(std::get<0>(XIN)),
//                       std::get<2>(std::get<1>(XIN))})),
//          bounds{{{minX, maxX}, {minY, maxY}, {minZ, maxZ}}} {};
//    Line(const Tddd &X0IN, const Tddd &X1IN)
//        : X({X0IN, X1IN}),
//          minX(Min(Tdd{std::get<0>(X0IN), std::get<0>(X1IN)})),
//          maxX(Max(Tdd{std::get<0>(X0IN), std::get<0>(X1IN)})),
//          minY(Min(Tdd{std::get<1>(X0IN), std::get<1>(X1IN)})),
//          maxY(Max(Tdd{std::get<1>(X0IN), std::get<1>(X1IN)})),
//          minZ(Min(Tdd{std::get<2>(X0IN), std::get<2>(X1IN)})),
//          maxZ(Max(Tdd{std::get<2>(X0IN), std::get<2>(X1IN)})),
//          bounds{{{minX, maxX}, {minY, maxY}, {minZ, maxZ}}} {};
// };
// /* ------------------------------------------------------ */
// struct Triangle {
//    T3Tddd X;
//    Tddd center;
//    T3Tdd bounds;
//    Tddd normal;
//    Triangle(const T3Tddd &XIN)
//        : X(XIN),
//          center(Mean(XIN)),
//          bounds{{{(Min(Tddd{std::get<0>(std::get<0>(XIN)),
//                             std::get<0>(std::get<1>(XIN)),
//                             std::get<0>(std::get<2>(XIN))})),
//                   (Max(Tddd{std::get<0>(std::get<0>(XIN)),
//                             std::get<0>(std::get<1>(XIN)),
//                             std::get<0>(std::get<2>(XIN))}))},
//                  {(Min(Tddd{std::get<1>(std::get<0>(XIN)),
//                             std::get<1>(std::get<1>(XIN)),
//                             std::get<1>(std::get<2>(XIN))})),
//                   (Max(Tddd{std::get<1>(std::get<0>(XIN)),
//                             std::get<1>(std::get<1>(XIN)),
//                             std::get<1>(std::get<2>(XIN))}))},
//                  {(Min(Tddd{std::get<2>(std::get<0>(XIN)),
//                             std::get<2>(std::get<1>(XIN)),
//                             std::get<2>(std::get<2>(XIN))})),
//                   (Max(Tddd{std::get<2>(std::get<0>(XIN)),
//                             std::get<2>(std::get<1>(XIN)),
//                             std::get<2>(std::get<2>(XIN))}))}}},
//          normal(Normalize(Cross(std::get<1>(X) - std::get<0>(X),
//                                 std::get<2>(X) - std::get<0>(X)))){};
//    Triangle(const Tddd &X0IN, const Tddd &X1IN, const Tddd &X2IN)
//        : X({X0IN, X1IN, X2IN}),
//          center((X0IN + X1IN + X2IN) / 3.),
//          bounds{{{(Min(Tddd{std::get<0>(X0IN), std::get<0>(X1IN),
//                             std::get<0>(X2IN)})),
//                   (Max(Tddd{std::get<0>(X0IN), std::get<0>(X1IN),
//                             std::get<0>(X2IN)}))},
//                  {(Min(Tddd{std::get<1>(X0IN), std::get<1>(X1IN),
//                             std::get<1>(X2IN)})),
//                   (Max(Tddd{std::get<1>(X0IN), std::get<1>(X1IN),
//                             std::get<1>(X2IN)}))},
//                  {(Min(Tddd{std::get<2>(X0IN), std::get<2>(X1IN),
//                             std::get<2>(X2IN)})),
//                   (Max(Tddd{std::get<2>(X0IN), std::get<2>(X1IN),
//                             std::get<2>(X2IN)}))}}},
//          normal(Normalize(Cross(X1IN - X0IN, X2IN - X0IN))){};
// };

// /* ------------------------------------------------------ */

// struct Sphere {
//    Tddd X;
//    double radius;
//    T3Tdd bounds;
//    Sphere(const Tddd &XIN, const double radiusIN = 0.)
//        : X(XIN),
//          radius(radiusIN),
//          bounds{{{(std::get<0>(XIN) - radiusIN), (std::get<0>(XIN) + radiusIN)},
//                  {(std::get<1>(XIN) - radiusIN), (std::get<1>(XIN) + radiusIN)},
//                  {(std::get<2>(XIN) - radiusIN),
//                   (std::get<2>(XIN) + radiusIN)}}} {};
// };
// };  // namespace geometry

/* ------------------------------------------------------ */

/*DOC_EXTRACT 0_1_coordinatebounds

### CoordinateBounds クラスについて

#### 概要
`CoordinateBounds`クラスは，
3次元の座標領域（バウンディングボックス）を管理するためのクラスである．

#### メンバ変数

| 変数名 | 型 | 説明 |
|--------|----|------|
| bounds | std::array<std::array<double,2>,3> | x, y, zそれぞれの範囲を保持する |
| X      | std::array<double,3>  | 中心座標を保持する |

#### メンバ関数

| メソッド名     | 戻り値型 | 説明 |
|--------------|----------|------|
| scaledBounds  | std::array<std::array<double,2>,3>    | 指定されたスケールでバウンディングボックスを拡大・縮小する |
| setBounds     | void     | バウンディングボックスを設定する（オーバーロードあり） |
| getXtuple     | const std::array<double,3> & | 中心座標を返す |
| getBounds     | const std::array<std::array<double,2>,3> & | バウンディングボックスの範囲を返す |
| Distance      | Tdd       | 指定座標との最小・最大距離を計算する |
| getVolume     | double    | バウンディングボックスの体積を計算する |
| getScale      | double    | バウンディングボックスのスケールを計算する |
| getCenter     | std::array<double,3>      | バウンディングボックスの中心座標を計算する |
| getVertices   | std::array<std::array<double,3>,8>    | バウンディングボックスの8つの頂点を計算する |

#### オペレーター

| オペレーター名 | 戻り値型 | 説明 |
|--------------|----------|------|
| ()            | const std::array<std::array<double,2>,3> & | バウンディングボックスの範囲を返す（関数呼び出しオペレータ） |

---

#### 有用性
CoordinateBounds クラスは、3次元空間での領域制限やクエリ処理、衝突判定などに使用できる。簡易的な操作で座標の範囲や距離、体積などを計算できるため、効率的な空間分割やデータ処理が可能である。

*/

// structをわざわざ作るのは，T3Tddではなく，coordinateboundsとして意味を具体的にした状態で持ち回りたいから．それだけ．

struct CoordinateBounds {
  T3Tdd bounds;
  // X is center
  Tddd X;
  Tddd X_temporary;
  /* -------------------------------------------------------- */
  Tdd xbounds() const { return std::get<0>(bounds); };
  Tdd ybounds() const { return std::get<1>(bounds); };
  Tdd zbounds() const { return std::get<2>(bounds); };
  /* ------------------------------------------------------ */
  T3Tdd getUniformBounds(const double scale = 1.) const {
    auto [x0, x1] = this->bounds[0];
    auto [y0, y1] = this->bounds[1];
    auto [z0, z1] = this->bounds[2];
    double L = std::max(std::abs(x1 - x0), std::abs(y1 - y0));
    L = std::max(L, std::abs(z1 - z0));
    Tddd center = {(x0 + x1) / 2., (y0 + y1) / 2., (z0 + z1) / 2.};
    return {{{center[0] - L * scale, center[0] + L * scale}, {center[1] - L * scale, center[1] + L * scale}, {center[2] - L * scale, center[2] + L * scale}}};
  };

  T3Tdd scaledBounds(const double scale) const {
    auto [xrange, yrange, zrange] = this->bounds;
    auto [x0, x1] = xrange;
    auto [y0, y1] = yrange;
    auto [z0, z1] = zrange;
    auto cX = 0.5 * (x0 + x1);
    auto cY = 0.5 * (y0 + y1);
    auto cZ = 0.5 * (z0 + z1);
    // xrange = {cX + scale * (x0 - cX), cX + scale * (x1 - cX)};
    // yrange = {cY + scale * (y0 - cY), cY + scale * (y1 - cY)};
    // zrange = {cZ + scale * (z0 - cZ), cZ + scale * (z1 - cZ)};
    xrange = {std::fma(scale, x0 - cX, cX), std::fma(scale, x1 - cX, cX)};
    yrange = {std::fma(scale, y0 - cY, cY), std::fma(scale, y1 - cY, cY)};
    zrange = {std::fma(scale, z0 - cZ, cZ), std::fma(scale, z1 - cZ, cZ)};
    return {{xrange, yrange, zrange}};
  };
  /* -------------------------------------------------------------------------- */
  void setBounds(const std::vector<Tddd>& Vxyz) {
    this->bounds = MinMaxColumns(Vxyz);
    this->X = Mean(Transpose(this->bounds));
  };
  void setBounds(const CoordinateBounds& bs) {
    this->bounds = bs.bounds;
    this->X = bs.X;
  };
  void setBounds(const T3Tddd& X012) {
    this->bounds = MinMaxColumns(X012);
    this->X = Mean(X012);
  };
  void setBounds(const Tddd& x) {
    this->bounds = {{{std::get<0>(x), std::get<0>(x)}, {std::get<1>(x), std::get<1>(x)}, {std::get<2>(x), std::get<2>(x)}}};
    this->X = x;
  };
  const Tddd& getXtuple() const { return this->X; /*面のsetBoundsでは，T3Tdddの平均がXとなるようになっている．バウンディングボックスの中心ではない.*/ };
  // const Tddd &getX() const { return this->X;
  // /*面のsetBoundsでは，T3Tdddの平均がXとなるようになっている．バウンディングボックスの中心ではない.*/
  // }; V_d getX() const { return {std::get<0>(this->X), std::get<1>(this->X),
  // std::get<2>(this->X)}; };
  const T3Tdd& getBounds() const { return this->bounds; };
  /* ------------------------------------------------------ */
  CoordinateBounds() : bounds{{{1E+20, -1E+20}, {1E+20, -1E+20}, {1E+20, -1E+20}}} {};
  CoordinateBounds(const CoordinateBounds& bs) : bounds(bs.bounds), X(bs.X) {};
  CoordinateBounds(const Tdd& minmaxX, const Tdd& minmaxY, const Tdd& minmaxZ) : bounds({minmaxX, minmaxY, minmaxZ}), X(Mean(Transpose(bounds))) {};
  CoordinateBounds(const Tddd& x) : bounds{{{std::get<0>(x), std::get<0>(x)}, {std::get<1>(x), std::get<1>(x)}, {std::get<2>(x), std::get<2>(x)}}}, X(x) {};
  CoordinateBounds(const T3Tdd& minmax) : bounds(minmax), X(Mean(Transpose(bounds))) {};
  CoordinateBounds(const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ) : bounds{{{minX, maxX}, {minY, maxY}, {minZ, maxZ}}}, X(Mean(Transpose(bounds))) {};
  CoordinateBounds(const T2Tddd& x) : bounds(MinMaxColumns(x)), X(Mean(x)) {};
  CoordinateBounds(const T3Tddd& x) : bounds(MinMaxColumns(x)), X(Mean(x)) {};
  CoordinateBounds(const T4Tddd& x) : bounds(MinMaxColumns(x)), X(Mean(x)) {};
  CoordinateBounds(const Tddd& x0, const Tddd& x1, const Tddd& x2) : bounds(MinMaxColumns(T3Tddd{x0, x1, x2})), X((x0 + x1 + x2) / 3.) {};
  CoordinateBounds(const std::vector<T3Tddd>& V_X);
  CoordinateBounds(const std::vector<Tddd>& X) : bounds(MinMaxColumns(X)), X(Mean(Transpose(bounds))) {};
  // CoordinateBounds(const geometry::Line &L) : bounds(MinMaxColumns(L.X)){};
  // CoordinateBounds(const geometry::Sphere &S) : bounds{{{std::get<0>(S.X) - S.radius, std::get<0>(S.X) + S.radius}, {std::get<1>(S.X) - S.radius, std::get<1>(S.X) + S.radius}, {std::get<2>(S.X) - S.radius, std::get<2>(S.X) + S.radius}}} {};
  // CoordinateBounds(const geometry::Triangle &T) : bounds(MinMaxColumns(T.X)){};
  /* ------------------------------------------------------ */
  Tdd Distance(const Tddd& a) const {
    auto [mmX, mmY, mmZ] = this->bounds;
    T3Tdd ret = this->bounds;
    //
    Tdd to_ax = mmX - std::get<0>(a);
    Tdd to_ay = mmY - std::get<1>(a);
    Tdd to_az = mmZ - std::get<2>(a);
    //
    std::get<0>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ax)), std::abs(std::get<1>(to_ax))}); //<-Tdd
    std::get<1>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ay)), std::abs(std::get<1>(to_ay))}); //<-Tdd
    std::get<2>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_az)), std::abs(std::get<1>(to_az))}); //<-Tdd
    //
    if (Between(std::get<0>(a), mmX))
      std::get<0>(std::get<0>(ret)) = 0;
    if (Between(std::get<1>(a), mmY))
      std::get<0>(std::get<1>(ret)) = 0;
    if (Between(std::get<2>(a), mmZ))
      std::get<0>(std::get<2>(ret)) = 0;
    auto [min, max] = Transpose(ret);
    return {Norm(min), Norm(max)};
  };
  /* ------------------------------------------------------ */
  const T3Tdd& operator()() const { return this->bounds; };
  const Tdd& getXMinMax() const { return std::get<0>(this->bounds); };
  const Tdd& getYMinMax() const { return std::get<1>(this->bounds); };
  const Tdd& getZMinMax() const { return std::get<2>(this->bounds); };
  double getVolume() const { return (std::get<1>(std::get<0>(this->bounds)) - std::get<0>(std::get<0>(this->bounds))) * (std::get<1>(std::get<1>(this->bounds)) - std::get<0>(std::get<1>(this->bounds))) * (std::get<1>(std::get<2>(this->bounds)) - std::get<0>(std::get<2>(this->bounds))); };
  double getScale() const {
    auto [min, max] = Transpose(this->bounds);
    return Norm(max - min);
  };
  Tddd getCenter() const { return Mean(Transpose(this->bounds)); };
  T8Tddd getVertices() const {
    auto [X0, X1] = std::get<0>(this->bounds);
    auto [Y0, Y1] = std::get<1>(this->bounds);
    auto [Z0, Z1] = std::get<2>(this->bounds);
    // connectivity:  {{x4, x6, x7, x5}, {x0, x2, x6, x4}, {x2, x3, x7, x6}, {x3, x1, x5, x7}, {x0, x4, x5, x1}, {x0, x1, x3, x2}}
    return T8Tddd{{{X0, Y0, Z0},   // 000, 0
                   {X0, Y1, Z0},   // 010, 1
                   {X1, Y1, Z0},   // 110, 2
                   {X1, Y0, Z0},   // 100, 3
                   {X0, Y0, Z1},   // 001, 4
                   {X0, Y1, Z1},   // 011, 5
                   {X1, Y1, Z1},   // 111, 6
                   {X1, Y0, Z1}}}; // 101, 7
  };
  //% ---------------- キャスト定義 -------------- */
  //% 型が明示され，関数よりもわかりやすい．
  operator T8Tddd() const { return this->getVertices(); };

  operator T6T4Tddd() const {
    auto [X0, X1] = std::get<0>(this->bounds);
    auto [Y0, Y1] = std::get<1>(this->bounds);
    auto [Z0, Z1] = std::get<2>(this->bounds);
    return {{{{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}}}, {{{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}}}, {{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}}}, {{{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}}}, {{{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}, {X0, Y1, Z0}}}, {{{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X1, Y1, Z0}}}}};
  };

  operator T12T2Tddd() const {
    auto [X0, X1] = std::get<0>(this->bounds);
    auto [Y0, Y1] = std::get<1>(this->bounds);
    auto [Z0, Z1] = std::get<2>(this->bounds);
    return {{/*X面*/ {{{X0, Y0, Z0}, {X1, Y0, Z0} /*to X*/}},
             {{{X0, Y0, Z0}, {X0, Y0, Z1} /*to Z*/}},
             /*X面*/ {{{X0, Y0, Z1}, {X1, Y0, Z1} /*to X*/}},
             {{{X1, Y0, Z0}, {X1, Y0, Z1} /*to Z*/}},
             /*X面*/ {{{X0, Y1, Z0}, {X1, Y0, Z0} /*to X*/}},
             {{{X0, Y0, Z0}, {X0, Y1, Z1} /*to Z*/}},
             /*X面*/ {{{X0, Y1, Z1}, {X1, Y0, Z1} /*to X*/}},
             {{{X1, Y0, Z0}, {X1, Y1, Z1} /*to Z*/}},
             /*Z面*/ {{{X0, Y0, Z0}, {X0, Y1, Z0} /*to Y*/}},
             {{{X1, Y0, Z0}, {X1, Y1, Z0} /*to Y*/}},
             /*Z面*/ {{{X0, Y0, Z1}, {X0, Y1, Z1} /*to Y*/}},
             {{{X1, Y0, Z1}, {X1, Y1, Z1} /*to Y*/}}}};
  };

  operator T12T3Tddd() const {
    auto [X0, X1] = std::get<0>(this->bounds);
    auto [Y0, Y1] = std::get<1>(this->bounds);
    auto [Z0, Z1] = std::get<2>(this->bounds);
    return T12T3Tddd{{{{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}}},
                      {{{X0, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}}},
                      {{{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}}},
                      {{{X0, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}}},
                      {{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}}},
                      {{{X0, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}}},
                      {{{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}}},
                      {{{X0, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}}},
                      {{{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}}},
                      {{{X0, Y0, Z0}, {X0, Y1, Z1}, {X0, Y1, Z0}}},
                      {{{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}}},
                      {{{X1, Y0, Z0}, {X1, Y1, Z1}, {X1, Y1, Z0}}}}};
  };
  //% -------------------------------------------- */
  std::tuple<CoordinateBounds, // 0
             CoordinateBounds, // 1
             CoordinateBounds, // 2
             CoordinateBounds, // 3
             CoordinateBounds, // 4
             CoordinateBounds, // 5
             CoordinateBounds, // 6
             CoordinateBounds> // 7
  to8Bounds(const double eps = 0.) const {
    /*
    +-----+-----+ Y1
    |     |     |
    +-----+-----+ Yc
    |     |     |
    +-----+-----+ Y0
    X0    Xc    X1
    */
    auto [X0, X1] = std::get<0>(this->bounds);
    auto [Y0, Y1] = std::get<1>(this->bounds);
    auto [Z0, Z1] = std::get<2>(this->bounds);
    auto [Xc, Yc, Zc] = this->getCenter();
    return {CoordinateBounds(Tdd{{X0 - eps, Xc + eps}}, Tdd{{Y0 - eps, Yc + eps}}, Tdd{{Z0 - eps, Zc + eps}}), CoordinateBounds(Tdd{{Xc - eps, X1 + eps}}, Tdd{{Y0 - eps, Yc + eps}}, Tdd{{Z0 - eps, Zc + eps}}), CoordinateBounds(Tdd{{X0 - eps, Xc + eps}}, Tdd{{Yc - eps, Y1 + eps}}, Tdd{{Z0 - eps, Zc + eps}}), CoordinateBounds(Tdd{{Xc - eps, X1 + eps}}, Tdd{{Yc - eps, Y1 + eps}}, Tdd{{Z0 - eps, Zc + eps}}),
            CoordinateBounds(Tdd{{X0 - eps, Xc + eps}}, Tdd{{Y0 - eps, Yc + eps}}, Tdd{{Zc - eps, Z1 + eps}}), CoordinateBounds(Tdd{{Xc - eps, X1 + eps}}, Tdd{{Y0 - eps, Yc + eps}}, Tdd{{Zc - eps, Z1 + eps}}), CoordinateBounds(Tdd{{X0 - eps, Xc + eps}}, Tdd{{Yc - eps, Y1 + eps}}, Tdd{{Zc - eps, Z1 + eps}}), CoordinateBounds(Tdd{{Xc - eps, X1 + eps}}, Tdd{{Yc - eps, Y1 + eps}}, Tdd{{Zc - eps, Z1 + eps}})};
  };
  std::tuple<CoordinateBounds, // 0
             CoordinateBounds, // 1
             CoordinateBounds, // 2
             CoordinateBounds, // 3
             CoordinateBounds, // 4
             CoordinateBounds, // 5
             CoordinateBounds, // 6
             CoordinateBounds> // 7
  to8Bounds(const Tddd& center) const {
    /*
    +-----+-----+ Y1
    |     |     |
    +-----+-----+ Yc
    |     |     |
    +-----+-----+ Y0
    X0    Xc    X1
    */
    auto [X0, X1] = std::get<0>(this->bounds);
    auto [Y0, Y1] = std::get<1>(this->bounds);
    auto [Z0, Z1] = std::get<2>(this->bounds);
    auto [Xc, Yc, Zc] = center;
    return {CoordinateBounds(Tdd{{X0, Xc}}, Tdd{{Y0, Yc}}, Tdd{{Z0, Zc}}), CoordinateBounds(Tdd{{Xc, X1}}, Tdd{{Y0, Yc}}, Tdd{{Z0, Zc}}), CoordinateBounds(Tdd{{X0, Xc}}, Tdd{{Yc, Y1}}, Tdd{{Z0, Zc}}), CoordinateBounds(Tdd{{Xc, X1}}, Tdd{{Yc, Y1}}, Tdd{{Z0, Zc}}),
            CoordinateBounds(Tdd{{X0, Xc}}, Tdd{{Y0, Yc}}, Tdd{{Zc, Z1}}), CoordinateBounds(Tdd{{Xc, X1}}, Tdd{{Y0, Yc}}, Tdd{{Zc, Z1}}), CoordinateBounds(Tdd{{X0, Xc}}, Tdd{{Yc, Y1}}, Tdd{{Zc, Z1}}), CoordinateBounds(Tdd{{Xc, X1}}, Tdd{{Yc, Y1}}, Tdd{{Zc, Z1}})};
  };
  inline bool InsideQ(const Tddd& X) const {
    if ((std::get<0>(X) < std::get<0>(std::get<0>(this->bounds)) || std::get<1>(std::get<0>(this->bounds)) < std::get<0>(X)) || (std::get<1>(X) < std::get<0>(std::get<1>(this->bounds)) || std::get<1>(std::get<1>(this->bounds)) < std::get<1>(X)) || (std::get<2>(X) < std::get<0>(std::get<2>(this->bounds)) || std::get<1>(std::get<2>(this->bounds)) < std::get<2>(X)))
      return false;
    else
      return true;
  };
};
/* -------------------------------------------------------------------------- */

inline bool InsideQ(const Tddd& X, const T3Tdd& bounds) {
  // point v.s. cube
  CoordinateBounds b(bounds);
  return b.InsideQ(X);
};
struct Sphere : public CoordinateBounds {
  Tddd center;
  double radius;
  Sphere(const Tddd& XIN, const double radiusIN = 0.) : CoordinateBounds(T3Tdd{{{(std::get<0>(XIN) - radiusIN), (std::get<0>(XIN) + radiusIN)}, {(std::get<1>(XIN) - radiusIN), (std::get<1>(XIN) + radiusIN)}, {(std::get<2>(XIN) - radiusIN), (std::get<2>(XIN) + radiusIN)}}}), center(XIN), radius(radiusIN) {};
};

inline Sphere InSphere(const Tddd& p0, const Tddd& p1, const Tddd& p2, const Tddd& p3) { return Sphere(Incenter(p0, p1, p2, p3), Inradius(p0, p1, p2, p3)); };

inline Sphere CircumSphere(const Tddd& p0, const Tddd& p1, const Tddd& p2, const Tddd& p3) { return Sphere(Circumcenter(p0, p1, p2, p3), Circumradius(p0, p1, p2, p3)); };

struct Triangle : public CoordinateBounds {
  std::array<std::array<double, 3>, 3> vertices;
  Tddd normal;
  Tddd angles;
  double area;
  // 中心点
  Tddd centroid;
  Tddd& center = this->centroid;
  //
  Tddd circumcenter /*外心*/;
  double circumradius /*外半径*/;
  //
  Tddd incenter /*内心*/;
  double inradius /*内接*/;
  //
  Triangle(const std::array<std::array<double, 3>, 3>& XIN) : CoordinateBounds(XIN), vertices(XIN), normal(TriangleNormalHP(XIN)), area(TriangleArea(XIN)), angles(TriangleAngles(XIN)), centroid(Centroid(XIN)), circumcenter(Circumcenter(XIN)), circumradius(Circumradius(XIN)), incenter(Incenter(XIN)), inradius(Inradius(XIN)) {};
  Triangle(const Tddd& X0IN, const Tddd& X1IN, const Tddd& X2IN) : CoordinateBounds(X0IN, X1IN, X2IN), vertices({X0IN, X1IN, X2IN}), normal(TriangleNormalHP(X0IN, X1IN, X2IN)), area(TriangleArea(X0IN, X1IN, X2IN)), angles(TriangleAngles(vertices)), centroid(Centroid(vertices)), circumcenter(Circumcenter(vertices)), circumradius(Circumradius(vertices)), incenter(Incenter(vertices)), inradius(Inradius(vertices)) {};

  void setProperties(const std::array<std::array<double, 3>, 3>& vertices_IN) {
    this->vertices = vertices_IN;
    this->normal = TriangleNormalHP(vertices_IN);
    this->area = TriangleArea(vertices_IN);
    this->angles = TriangleAngles(vertices_IN);
    this->centroid = Centroid(vertices_IN);
    this->circumcenter = Circumcenter(vertices_IN);
    this->circumradius = Circumradius(vertices_IN);
    this->incenter = Incenter(vertices_IN);
    this->inradius = Inradius(vertices_IN);
  };
  operator const std::array<std::array<double, 3>, 3>&() const { return vertices; };
  Tddd normal_to_be_preserved;
};

/* -------------------------------------------------------------------------- */

struct Tetrahedron : public CoordinateBounds {
  //     1,3,2
  //          3
  // 0,2,3  / | \ 0,3,1     --- 1 ---
  //       2--|--1           \      /
  // 0,1,2  \ | /             2    0
  //          0                 \/
  //
  //$ GEOMETRIC PROPERTIES
  T4Tddd vertices;
  double volume;
  Tddd centroid;
  //$ circumscribed sphere (circumsphere)
  Tddd circumcenter /*外心*/;
  double circumradius /*外半径*/;
  //$ inscribed sphere (insphere)
  Tddd incenter /*内心*/;
  double inradius /*内接*/;
  //
  T4Tddd normals;
  T4d solidangles; // いつかチェック
  Tetrahedron(const Tddd& X0IN, const Tddd& X1IN, const Tddd& X2IN, const Tddd& X3IN) : CoordinateBounds(T4Tddd{X0IN, X1IN, X2IN, X3IN}), vertices(T4Tddd{X0IN, X1IN, X2IN, X3IN}), volume(TetrahedronVolume(vertices)), centroid(Centroid(vertices)), circumcenter(Circumcenter(vertices)), circumradius(Circumradius(vertices)), incenter(Incenter(vertices)), inradius(Inradius(vertices)), normals(TetrahedronNormals(vertices)), solidangles(TetrahedronSolidAngle_UsingVectorAngle(vertices)) {};

  Tetrahedron(const T4Tddd& XIN) : CoordinateBounds(XIN), vertices(XIN), volume(TetrahedronVolume(XIN)), centroid(Centroid(XIN)), circumcenter(Circumcenter(XIN)), circumradius(Circumradius(XIN)), incenter(Incenter(XIN)), inradius(Inradius(XIN)), normals(TetrahedronNormals(XIN)), solidangles(TetrahedronSolidAngle_UsingVectorAngle(XIN)) {};

  void setProperties(const T4Tddd& vertices_IN) {
    this->vertices = vertices_IN;
    CoordinateBounds::setBounds(this->vertices);
    this->volume = TetrahedronVolume(this->vertices);
    this->centroid = Centroid(this->vertices);
    this->circumcenter = Circumcenter(this->vertices);
    this->circumradius = Circumradius(this->vertices);
    this->incenter = Incenter(this->vertices);
    this->inradius = Inradius(this->vertices);
    this->normals = TetrahedronNormals(this->vertices);
    this->solidangles = TetrahedronSolidAngle_UsingVectorAngle(this->vertices);
  };

  Tetrahedron scaled(const auto& s = 0.9) { return Tetrahedron({(std::get<0>(this->vertices) - centroid) * s + centroid, (std::get<1>(this->vertices) - centroid) * s + centroid, (std::get<2>(this->vertices) - centroid) * s + centroid, (std::get<3>(this->vertices) - centroid) * s + centroid}); };

  inline bool InsideQ(const std::array<double, 3>& X) const {
    if (CoordinateBounds::InsideQ(X)) {
      const auto [x, y, z] = X;
      const auto [X0, X1, X2, X3] = this->vertices;
      std::array<double, 4> W;
      lapack_svd_solve(std::array<std::array<double, 4>, 4>{{{X0[0], X1[0], X2[0], X3[0]}, {X0[1], X1[1], X2[1], X3[1]}, {X0[2], X1[2], X2[2], X3[2]}, {1., 1., 1., 1.}}}, W, {x, y, z, 1.});
      return 0. <= W[0] && W[0] <= 1. && 0. <= W[1] && W[1] <= 1. - W[0] && 0. <= W[2] && W[2] <= 1. - W[0] - W[1];
    }
    return false;
  };

  // T4d SolidAngles(const T4Tddd &X) {
  //    auto [X0, X1, X2, X3] = X;
  //    return {(SolidAngle_(X0, X1, X2, X3)),
  //            (SolidAngle_(X1, X0, X3, X2)),
  //            (SolidAngle_(X2, X0, X1, X3)),
  //            (SolidAngle_(X3, X0, X2, X1))};
  // };

  operator T6T2Tddd() const {
    auto [p0, p1, p2, p3] = this->vertices;
    return {{{p0, p1}, {p0, p2}, {p0, p3}, {p1, p2}, {p2, p3}, {p3, p1}}};
  };

  // operator T6Tddd() const {
  //    auto [p0, p1, p2, p3] = this->vertices;
  //    return {(p0 + p1) / 2.,
  //            (p0 + p2) / 2.,
  //            (p0 + p3) / 2.,
  //            (p1 + p2) / 2.,
  //            (p2 + p3) / 2.,
  //            (p3 + p1) / 2.};
  // };

  operator T4T3Tddd() const {
    auto [p0, p1, p2, p3] = this->vertices;
    return {{{p0, p1, p2}, {p0, p1, p3}, {p0, p2, p3}, {p1, p2, p3}}};
  };

  operator T4Tddd() const {
    auto [p0, p1, p2, p3] = this->vertices;
    return {p0, p1, p2, p3};
  };
};

/* -------------------------------------------------------------------------- */

struct Ellipsoid : public CoordinateBounds {
  Quaternion Q;   // Rotation matrix
  Tddd center;    // Center of the ellipsoid
  double a, b, c; // Semi-axes lengths in the x, y, and z directions

  // Constructor
  Ellipsoid(const double aIN, const double bIN, const double cIN, const Tddd& centerIN = {0., 0., 0.}, const Quaternion& QIN = Quaternion()) { this->setProperties(aIN, bIN, cIN, centerIN, QIN); }

  // Method to check if a point is inside the ellipsoid
  inline bool InsideQ(const Tddd& point) const {
    const auto [x, y, z] = Dot(point - center, Q.Rv());
    return (std::pow(x / this->a, 2) + std::pow(y / this->b, 2) + std::pow(z / this->c, 2)) <= 1.0;
  }

  // Function to set properties of the ellipsoid (if they change)
  void setProperties(const double aIN, const double bIN, const double cIN) {
    a = aIN;
    b = bIN;
    c = cIN;
    double max = a;
    if (b > max)
      max = b;
    if (c > max)
      max = c;
    this->setBounds(std::vector<Tddd>{this->center + max, this->center - max});
  }
  void setProperties(const double aIN, const double bIN, const double cIN, const Tddd& centerIN) {
    this->center = centerIN;
    setProperties(aIN, bIN, cIN);
  }
  void setProperties(const double aIN, const double bIN, const double cIN, const Tddd& centerIN, const Quaternion& QIN) {
    this->Q = QIN;
    setProperties(aIN, bIN, cIN, centerIN);
  }
};

/* -------------------------------------------------------------------------- */
inline std::ostream& operator<<(std::ostream& stream, const CoordinateBounds& bounds) { return (stream << bounds.bounds); };
inline CoordinateBounds operator+(const CoordinateBounds& b0, const CoordinateBounds& b1) {
  auto [minX0, maxX0] = std::get<0>(b0.bounds);
  auto [minX1, maxX1] = std::get<0>(b1.bounds);
  auto [minY0, maxY0] = std::get<1>(b0.bounds);
  auto [minY1, maxY1] = std::get<1>(b1.bounds);
  auto [minZ0, maxZ0] = std::get<2>(b0.bounds);
  auto [minZ1, maxZ1] = std::get<2>(b1.bounds);

  return CoordinateBounds(std::min(minX0, minX1), std::max(maxX0, maxX1), std::min(minY0, minY1), std::max(maxY0, maxY1), std::min(minZ0, minZ1), std::max(maxZ0, maxZ1));
};

inline CoordinateBounds& operator+=(CoordinateBounds& b0, const CoordinateBounds& b1) {
  auto [minX0, maxX0] = std::get<0>(b0.bounds);
  auto [minX1, maxX1] = std::get<0>(b1.bounds);
  auto [minY0, maxY0] = std::get<1>(b0.bounds);
  auto [minY1, maxY1] = std::get<1>(b1.bounds);
  auto [minZ0, maxZ0] = std::get<2>(b0.bounds);
  auto [minZ1, maxZ1] = std::get<2>(b1.bounds);

  b0.bounds = {{{std::min(minX0, minX1), std::max(maxX0, maxX1)}, {std::min(minY0, minY1), std::max(maxY0, maxY1)}, {std::min(minZ0, minZ1), std::max(maxZ0, maxZ1)}}};
  return b0;
};
inline CoordinateBounds::CoordinateBounds(const std::vector<T3Tddd>& V_X) {
  CoordinateBounds ret;
  for (const auto& X : V_X)
    ret += CoordinateBounds(X);
  this->bounds = ret.bounds;
};
/* ------------------------------------------------------ */
struct IntersectionSphereLine {
  /*
     /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/三角形と球の干渉.nb
     X = (p0,p1).(t,1-t) = (p0-p1)*t + p1 であることから，
     t = 0 -> X = p0
     t = 1 -> X = p1
     となる
  */
  const Tddd p0, p1;
  double t;
  Tddd X;
  double distance;
  bool isIntersecting;
  IntersectionSphereLine(const Tddd& center, const double radius, const T2Tddd& p01) : p0(std::get<0>(p01)), p1(std::get<1>(p01)), t(Dot(center - p1, p0 - p1) / Dot(p0 - p1, p0 - p1)), X(p0 * t + p1 * (1 - t)), distance(Norm(X - center)), isIntersecting(distance <= radius && 0. <= t && t <= 1.) {
    if (!isIntersecting) {
      // 線分と干渉しない場合でも，球に最も近い線分上の点を返すようにする．
      // これを実行しない場合，線分ではなく，直線上の点を返すことになる．
      if (Norm(p0 - center) < Norm(p1 - center)) {
        this->t = 1;
        this->X = p0;
        this->distance = Norm(p0 - center);
      } else {
        this->t = 0;
        this->X = p1;
        this->distance = Norm(p1 - center);
      }
    }
  };
  IntersectionSphereLine(const double radius, const T2Tddd& p01) : p0(std::get<0>(p01)), p1(std::get<1>(p01)), t(Dot(-p1, p0 - p1) / Dot(p0 - p1, p0 - p1)), X(p0 * t + p1 * (1 - t)), distance(Norm(X)), isIntersecting(distance <= radius && 0. <= t && t <= 1.) {
    if (!isIntersecting) {
      // 線分と干渉しない場合でも，球に最も近い線分上の点を返すようにする．
      // これを実行しない場合，線分ではなく，直線上の点を返すことになる．
      if (Norm(p0) < Norm(p1)) {
        this->t = 1;
        this->X = p0;
        this->distance = Norm(p0);
      } else {
        this->t = 0;
        this->X = p1;
        this->distance = Norm(p1);
      }
    }
  };
};
/* ------------------------------------------------------ */
struct IntersectionTriangles {
  /*
  /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/三角形と球の干渉.nb
  */
  T3Tddd P0, P1;
  T2Tddd L;
  bool isIntersecting;
  double eps = 0.; // 1E-13;
  double a0, b0, a1, b1, min0, max0, min1, max1, deno0, deno1;
  double I10InT, I11InT, I00InT, I01InT;
  Tddd I00, I01, I10, I11;
  double lmitMin(const double a0, const double b0) const {
    Tddd ret = {0., 0., 0.};
    if (std::abs(1 - a0) > eps) {
      if ((1 - a0) > 0)
        std::get<0>(ret) = -b0 / (1 - a0);
      else
        std::get<0>(ret) = (1 - b0) / (1 - a0);
    }

    if (std::abs(-a0) > eps) {
      if (-a0 > 0)
        std::get<1>(ret) = b0 / a0;
      else
        std::get<1>(ret) = -(1 - b0) / a0;
    }
    return Max(ret);
  };

  double lmitMax(const double a0, const double b0) const {
    Tddd ret = {1., 1., 1.};
    if (std::abs(1 - a0) > eps) {
      if ((1 - a0) > 0)
        std::get<0>(ret) = (1 - b0) / (1 - a0);
      else
        std::get<0>(ret) = -b0 / (1 - a0);
    }

    if (std::abs(-a0) > eps) {
      if (-a0 > 0)
        std::get<1>(ret) = -(1 - b0) / a0;
      else
        std::get<1>(ret) = b0 / a0;
    }
    return Min(ret); // この内最も小さいものが，上限となる
  };

  Tddd p00, p01, p02;
  Tddd p10, p11, p12;
  Tddd normalP0, normalP1;
  /* -------------------------------------------------------------- */
  IntersectionTriangles(const T3Tddd& P0_IN, const T3Tddd& P1_IN) : P0(P0_IN), P1(P1_IN), isIntersecting(false), p00(std::get<0>(P0)), p01(std::get<1>(P0)), p02(std::get<2>(P0)), p10(std::get<0>(P1)), p11(std::get<1>(P1)), p12(std::get<2>(P1)), normalP0(Normalize(Cross(p01 - p00, p02 - p00))), normalP1(Normalize(Cross(p11 - p10, p12 - p10))) {
    // if (IntersectQ(CoordinateBounds(P0), CoordinateBounds(P0)))
    // 	return;
    // auto [p00, p01, p02] = P0;
    // auto [p10, p11, p12] = P1;
    // auto normalP0 = Normalize(Cross(p01 - p00, p02 - p00));
    // auto normalP1 = Normalize(Cross(p11 - p10, p12 - p10));
    if (std::abs(std::abs(Dot(normalP0, normalP1)) - 1) < eps) {
      // Print("干渉線を定義することはできない");
      // Print("線を共有するような場合はありえる");
      return;
    }
    deno0 = Dot(normalP1, p01 - p02);
    deno1 = Dot(normalP0, p11 - p12);

    if (std::abs(deno0) < eps)
      for (auto i = 1; i < 3; ++i) {
        P0 = RotateLeft(P0, i);
        p00 = std::get<0>(P0);
        p01 = std::get<1>(P0);
        p02 = std::get<2>(P0);
        deno0 = Dot(normalP1, p01 - p02);
        if (!(std::abs(deno0) < eps))
          break;
        if (i == 2)
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }

    if (std::abs(deno1) < eps)
      for (auto i = 1; i < 3; ++i) {
        P1 = RotateLeft(P1, i);
        p10 = std::get<0>(P1);
        p11 = std::get<1>(P1);
        p12 = std::get<2>(P1);
        deno1 = Dot(normalP0, p11 - p12);
        if (!(std::abs(deno1) < eps))
          break;
        if (i == 2)
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      }

    a0 = Dot(normalP1, p00 - p02) / deno0;
    b0 = Dot(normalP1, p12 - p02) / deno0;
    a1 = Dot(normalP0, p10 - p12) / deno1;
    b1 = Dot(normalP0, p02 - p12) / deno1;

    min0 = lmitMin(a0, b0);
    max0 = lmitMax(a0, b0);
    min1 = lmitMin(a1, b1);
    max1 = lmitMax(a1, b1);

    auto t = min0;
    I00 = Dot(Tddd{{t, b0 - a0 * t, 1 - t - (b0 - a0 * t)}}, P0);
    t = max0;
    I01 = Dot(Tddd{{t, b0 - a0 * t, 1 - t - (b0 - a0 * t)}}, P0);
    t = min1;
    I10 = Dot(Tddd{{t, b1 - a1 * t, 1 - t - (b1 - a1 * t)}}, P1);
    t = max1;
    I11 = Dot(Tddd{{t, b1 - a1 * t, 1 - t - (b1 - a1 * t)}}, P1);

    auto tmp = (I10 - I00) / (I01 - I00);
    I10InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
    auto I10IsInside = Between(I10InT, {0, 1});

    tmp = (I11 - I00) / (I01 - I00);
    I11InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
    auto I11IsInside = Between(I11InT, {0, 1});

    tmp = (I00 - I10) / (I11 - I10);
    I00InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
    auto I00IsInside = Between(I00InT, {0, 1});

    tmp = (I01 - I10) / (I11 - I10);
    I01InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
    auto I01IsInside = Between(I01InT, {0, 1});

    if (I10IsInside) {
      if (I11IsInside && Norm(I10 - I11) > eps) {
        L = {I10, I11};
        isIntersecting = true;
      }
      if (I01IsInside && Norm(I10 - I01) > eps) {
        L = {I10, I01};
        isIntersecting = true;
      }
      if (I00IsInside && Norm(I10 - I00) > eps) {
        L = {I10, I00};
        isIntersecting = true;
      }
    }
    if (I11IsInside) {
      if (I00IsInside && Norm(I11 - I00) > eps) {
        L = {I11, I00};
        isIntersecting = true;
      }
      if (I01IsInside && Norm(I11 - I01) > eps) {
        L = {I11, I01};
        isIntersecting = true;
      }
    }
    if (I00IsInside && I01IsInside && Norm(I00 - I01) > eps) {
      L = {I00, I01};
      isIntersecting = true;
    }
  };
};

/* ------------------------------------------------------ */
struct IntersectionLineTriangle {
  Tddd v1, v2, v3;   // vertices of the triangle
  Tddd p0, p1;       // endpoints of the line segment
  Tddd intersection; // intersection point
  bool isIntersecting;

  IntersectionLineTriangle(const Tddd& v1, const Tddd& v2, const Tddd& v3, const Tddd& p0, const Tddd& p1) : v1(v1), v2(v2), v3(v3), p0(p0), p1(p1), isIntersecting(false) {
    // Implement Möller–Trumbore intersection algorithm here
    Tddd h, s, q;
    double a, f, u, v;
    Tddd edge1 = v2 - v1;
    Tddd edge2 = v3 - v1;
    h = Cross(p1 - p0, p0 - p1);
    a = Dot(edge1, h);

    if (a > -1e-5 && a < 1e-5)
      return; // The line is parallel to the triangle.

    f = 1.0 / a;
    s = p0 - v1;
    u = f * Dot(s, h);

    if (u < 0.0 || u > 1.0)
      return;

    q = Cross(s, edge1);
    v = f * Dot(p1 - p0, q);

    if (v < 0.0 || u + v > 1.0)
      return;

    // At this stage we can compute the intersection point
    double t = f * Dot(edge2, q);

    if (t > 1e-5) { // Ray intersection
      intersection = p0 + (p1 - p0) * t;
      isIntersecting = true;
    }
  }
};

/* ------------------------------------------------------ */
struct IntersectionSphereTriangle {
  double eps = 1E-8;
  Tddd X0, X1, X2, center;
  Tddd n;
  double n0, n1, n2;
  double x00, x01, x02, x10, x11, x12, x20, x21, x22;
  double denominator;
  double scale, t0, t1;
  double distance;
  Tddd X;
  const IntersectionSphereLine intxnL0;
  const IntersectionSphereLine intxnL1;
  const IntersectionSphereLine intxnL2;
  bool isIntersecting;
  bool isIntersectingInsideTriangle;
  bool isIntersectingEdgeTriangle;
  //
  // double F(const double t0, const double t1) {
  //    auto v = Dot({t0, t1, 1}, this->X0X1X2);
  //    return Dot(v, v);
  // };
  // double dFdt(const double t0, const double t1) {
  //    auto v = Dot({t0, t1, 1}, this->X0X1X2);
  //    auto dvdt0 = Dot({1, 0, 0}, this->X0X1X2);
  //    auto dvdt1 = Dot({0, 1, 0}, this->X0X1X2);
  //    return Dot(v, v);
  // };
  IntersectionSphereTriangle(const Tddd& centerIN, const double radius, const T3Tddd& X0X1X2)
      : center(centerIN), X0(std::get<0>(X0X1X2) - centerIN), X1(std::get<1>(X0X1X2) - centerIN), X2(std::get<2>(X0X1X2) - centerIN), x00(std::get<0>(X0)), x01(std::get<1>(X0)), x02(std::get<2>(X0)), x10(std::get<0>(X1)), x11(std::get<1>(X1)), x12(std::get<2>(X1)), x20(std::get<0>(X2)), x21(std::get<1>(X2)), x22(std::get<2>(X2)), n(Normalize(Cross(X1 - X0, X2 - X0))), n0(std::get<0>(n)), n1(std::get<1>(n)), n2(std::get<2>(n)),
        denominator(n2 * (x00 * x11 - x11 * x20 + x01 * (-x10 + x20) - x00 * x21 + x10 * x21) + n1 * (x02 * x10 - x00 * x12 - x02 * x20 + x12 * x20 + x00 * x22 - x10 * x22) + n0 * (-(x02 * x11) + x01 * x12 + x02 * x21 - x12 * x21 - x01 * x22 + x11 * x22)), t0((-(n2 * x11 * x20) + n1 * x12 * x20 + n2 * x10 * x21 - n0 * x12 * x21 - n1 * x10 * x22 + n0 * x11 * x22) / denominator),
        t1((n2 * x01 * x20 - n1 * x02 * x20 - n2 * x00 * x21 + n0 * x02 * x21 + n1 * x00 * x22 - n0 * x01 * x22) / denominator), scale((-(x02 * x11 * x20) + x01 * x12 * x20 + x02 * x10 * x21 - x00 * x12 * x21 - x01 * x10 * x22 + x00 * x11 * x22) / denominator), intxnL0(1E+20, T2Tddd{X0, X1}), intxnL1(1E+20, T2Tddd{X1, X2}), intxnL2(1E+20, T2Tddd{X2, X0}), distance(std::abs(scale)), X(scale * n + centerIN),
        isIntersecting(distance <= radius && (0 <= t0 && t0 <= 1) && (0 <= t1 && t1 <= 1) && (0 <= (1 - t0 - t1) && (1 - t0 - t1) <= 1)), isIntersectingInsideTriangle(isIntersecting), isIntersectingEdgeTriangle(false) {
    /*
     *   P = (p0, p1, p2)
     *   P.T = p0*t0 + p1*t1 + p2*t2
     *
     *                           p0 <- T=(1,0,0)
     *                          *
     *     P.(t0,1-t0,0) ->    / \   <- P.(1-t2,0,t2) = p0*(1-t2) + p2*t2
     *  =>(p0,p1).(t0,1-t0) # /   \     => (p0,p2).(1-t2,t2)
     *                       *-----*    => (p2,p0).(t2,1-t2) #
     *       T=(0,1,0) -> p1           p2 <- T=(0,0,1)
     *                     P.(0,t1,1-t1)
     *                   => (p1,p2).(t1,1-t1) #
     */
    if (!isIntersectingInsideTriangle) {
      double online_dist = 1E+10;
      if (intxnL0.isIntersecting /*distanceのチェックは不要．必ず上書きする*/) {
        this->t0 = intxnL0.t;
        this->t1 = 1 - intxnL0.t;
        this->X = intxnL0.X + center;
        online_dist = this->distance = intxnL0.distance;
        this->isIntersectingInsideTriangle = false;
        this->isIntersectingEdgeTriangle = this->isIntersecting = true;
      }
      if (intxnL1.isIntersecting && intxnL1.distance <= online_dist) {
        this->t0 = 0;
        this->t1 = intxnL1.t;
        this->X = intxnL1.X + center;
        online_dist = this->distance = intxnL1.distance;
        this->isIntersectingInsideTriangle = false;
        this->isIntersectingEdgeTriangle = this->isIntersecting = true;
      }
      if (intxnL2.isIntersecting && intxnL2.distance <= online_dist) {
        this->t0 = 1 - intxnL2.t;
        this->t1 = 0;
        this->X = intxnL2.X + center;
        online_dist = this->distance = intxnL2.distance;
        this->isIntersectingInsideTriangle = false;
        this->isIntersectingEdgeTriangle = this->isIntersecting = true;
      }
    }
    if (!isIntersectingInsideTriangle || !isIntersectingEdgeTriangle) {
      if (Norm(X0) < this->distance) {
        Norm(this->X = X0 + center);
        this->distance = Norm(X0);
      }
      if (Norm(X1) < this->distance) {
        Norm(this->X = X1 + center);
        this->distance = Norm(X1);
      }
      if (Norm(X2) < this->distance) {
        Norm(this->X = X2 + center);
        this->distance = Norm(X2);
      }
    }
  };
};
/* ------------------------------------------------------ */
struct IntersectionSphereTriangleLimitedToNormalRegion {
  const double eps = 1E-13;
  const T3Tddd P;
  const Tddd center;
  //
  const Tddd p0, p1, p2;
  const Tddd n;
  const double nx, ny, nz;
  const double p02x, p02y, p02z, p12x, p12y, p12z;
  const double determ;
  const T3Tddd mat;
  const Tddd ans;
  const double t0, t1, scale;
  const Tddd X;
  const bool isIntersecting;
  IntersectionSphereTriangleLimitedToNormalRegion(const Tddd& centerIN, const double radius, const T3Tddd& P_IN)
      : P(P_IN), center(centerIN), p0(std::get<0>(P_IN)), p1(std::get<1>(P_IN)), p2(std::get<2>(P_IN)), n(Normalize(Cross(p1 - p0, p2 - p0))), nx(std::get<0>(n)), ny(std::get<1>(n)), nz(std::get<2>(n)), p02x(std::get<0>(p0 - p2)), p02y(std::get<1>(p0 - p2)), p02z(std::get<2>(p0 - p2)), p12x(std::get<0>(p1 - p2)), p12y(std::get<1>(p1 - p2)), p12z(std::get<2>(p1 - p2)), determ(-(nz * p02y * p12x) + ny * p02z * p12x + nz * p02x * p12y - nx * p02z * p12y - ny * p02x * p12z + nx * p02y * p12z),
        mat({{{nz * p12y - ny * p12z, -(nz * p02y) + ny * p02z, p02z * p12y - p02y * p12z}, {-(nz * p12x) + nx * p12z, nz * p02x - nx * p02z, -(p02z * p12x) + p02x * p12z}, {ny * p12x - nx * p12y, -(ny * p02x) + nx * p02y, p02y * p12x - p02x * p12y}}}), ans(Dot(center - p2, mat / determ)), t0(std::get<0>(ans)), t1(std::get<1>(ans)), scale(std::get<2>(ans)), X(scale * n + center),
        isIntersecting((std::abs(determ) > eps) && (std::abs(scale) <= radius && (0 <= t0 && t0 <= 1) && (0 <= t1 && t1 <= 1) && (0 <= (1 - t0 - t1) && (1 - t0 - t1) <= 1))) {};
  Tddd getNearestX() const {
    if (this->isIntersecting)
      return this->X;
    else {
      Tddd X_ = {1E+50, 1E+50, 1E+50};
      double mindistance = 1E+50;
      auto intxnL0 = IntersectionSphereLine(center, 1E+50, T2Tddd{std::get<0>(P), std::get<1>(P)});
      auto intxnL1 = IntersectionSphereLine(center, 1E+50, T2Tddd{std::get<1>(P), std::get<2>(P)});
      auto intxnL2 = IntersectionSphereLine(center, 1E+50, T2Tddd{std::get<2>(P), std::get<0>(P)});
      if (intxnL0.isIntersecting) {
        X_ = intxnL0.X;
        mindistance = Norm(X_ - center);
      }
      if (intxnL1.isIntersecting) {
        if (Norm(intxnL1.X - center) < mindistance) {
          X_ = intxnL1.X;
          mindistance = Norm(X_ - center);
        }
      }
      if (intxnL2.isIntersecting) {
        if (Norm(intxnL2.X - center) < mindistance) {
          X_ = intxnL2.X;
          mindistance = Norm(X_ - center);
        }
      }
      if (Norm(std::get<0>(P) - center) < mindistance) {
        X_ = std::get<0>(P);
        mindistance = Norm(X_ - center);
      }
      if (Norm(std::get<1>(P) - center) < mindistance) {
        X_ = std::get<1>(P);
        mindistance = Norm(X_ - center);
      }
      if (Norm(std::get<2>(P) - center) < mindistance) {
        X_ = std::get<2>(P);
      }
      return X_;
    }
  };
};

// inline std::tuple<double, Tddd> Nearest_(const Tddd &X, const T2Tddd &ab) {
//   /*
//   a * t + b * (1-t)
//   ---------------------------
//   ( a*t+b*(1-t) - X ).(a-b) = 0
//   ( (a-b)*t + (b - X) ).(a-b) = 0
//   t = (X-b).(a-b)/(a-b).(a-b)
//   */
//   // const auto a_b = std::get<0>(ab) - std::get<1>(ab);
//   // const auto t = std::clamp(Dot(X - std::get<1>(ab), a_b) / Dot(a_b, a_b), 0.0, 1.0);
//   // return {t, std::get<0>(ab) * t + std::get<1>(ab) * (1. - t)};

//   const auto &a = std::get<0>(ab);
//   const auto &b = std::get<1>(ab);
//   const auto a_b = a - b;
//   const auto dot_a_b = Dot(a_b, a_b);
//   if (dot_a_b == 1e-30)
//     return {0.5, 0.5 * (a + b)};
//   const auto t = std::clamp(Dot(X - b, a_b) / dot_a_b, 0.0, 1.0);
//   return {t, a * t + b * (1.0 - t)};
// };

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

inline double AbsMax3(const Tddd& v) { return std::max({std::abs(v[0]), std::abs(v[1]), std::abs(v[2])}); }

inline std::tuple<double, Tddd> Nearest_(const Tddd& X, const T2Tddd& ab) {
  const auto& a = std::get<0>(ab);
  const auto& b = std::get<1>(ab);

  const Tddd d = a - b; // b -> a
  const Tddd xb = X - b;

  // --- 入力が非有限なら壊れない返しにする（必要なら assert に変更） ---
  const double sa = AbsMax3(a), sb = AbsMax3(b), sx = AbsMax3(X);
  if (!std::isfinite(sa) || !std::isfinite(sb) || !std::isfinite(sx)) {
    return {0.5, (a + b) * 0.5};
  }

  // --- スケールに応じた縮退判定（相対閾値） ---
  // 代表スケール：座標の大きさ（最低 1.0）
  const double S = std::max({sa, sb, 1.0});
  // 線分長の閾値：係数は「安全側」(用途で調整)
  const double tol = 64.0 * std::numeric_limits<double>::epsilon() * S;

  const double sd = AbsMax3(d);
  if (!(sd > tol)) { // sd==0 も tol 未満もまとめて縮退扱い
    return {0.5, (a + b) * 0.5};
  }

  // --- スケーリングで overflow/underflow 回避 ---
  // d = sd * dn,  |dn|_inf = 1
  const double inv_sd = 1.0 / sd;
  const Tddd dn = d * inv_sd;

  const double den = Dot(dn, dn); // ~ O(1)（最大でも 3 程度）
  const double num = Dot(xb, dn);

  if (!(den > 0.0) || !std::isfinite(num) || !std::isfinite(den)) {
    return {0.5, (a + b) * 0.5};
  }

  double t = num / (sd * den);
  if (!std::isfinite(t)) {
    return {0.5, (a + b) * 0.5};
  }

  // --- watertight のため端点は厳密に返す ---
  if (t <= 0.0)
    return {0.0, b};
  if (t >= 1.0)
    return {1.0, a};

  return {t, b + d * t};
}

inline Tdd Nearest_(const T2Tddd& ab, const T2Tddd& AB) {
  const auto [a, b] = ab;
  const auto [A, B] = AB;
  const auto a_b = a - b;
  const auto A_B = A - B;
  // const auto [t, tau] = Dot(Inverse(T2Tdd{{{Dot(a_b, a_b), -Dot(A_B, a_b)},
  //                                          {Dot(a_b, A_B), -Dot(A_B, A_B)}}}),
  //                           Tdd{-Dot(b, a_b) + Dot(B, a_b),
  //                               -Dot(b, A_B) + Dot(B, A_B)});
  // use lapack_lu
  std::array<double, 2> ans;
  // lapack_lu(ans, T2Tdd{{{Dot(a_b, a_b), -Dot(A_B, a_b)}, {Dot(a_b, A_B), -Dot(A_B, A_B)}}}, Tdd{-Dot(b, a_b) + Dot(B, a_b), -Dot(b, A_B) + Dot(B, A_B)});
  // Solve(ans,
  //       T2Tdd{{{Dot(a_b, a_b), -Dot(A_B, a_b)}, {Dot(a_b, A_B), -Dot(A_B, A_B)}}},
  //       Tdd{-Dot(b, a_b) + Dot(B, a_b), -Dot(b, A_B) + Dot(B, A_B)});
  // lapack_lu lu(ans,
  //              T2Tdd{{{Dot(a_b, a_b), -Dot(A_B, a_b)}, {Dot(a_b, A_B), -Dot(A_B, A_B)}}},
  //              Tdd{-Dot(b, a_b) + Dot(B, a_b), -Dot(b, A_B) + Dot(B, A_B)});
  lapack_lu lu(ans, T2Tdd{{{Dot(a_b, a_b), -Dot(A_B, a_b)}, {Dot(a_b, A_B), -Dot(A_B, A_B)}}}, Tdd{-Dot(b, a_b) + Dot(B, a_b), -Dot(b, A_B) + Dot(B, A_B)});
  const auto [t, tau] = ans;

  if (Between(t, array_0_1) && Between(tau, array_0_1))
    return {t, tau};
  auto [t0, X0] = Nearest_(a, {A, B});
  auto [t1, X1] = Nearest_(b, {A, B});
  auto [t2, X2] = Nearest_(A, {a, b});
  auto [t3, X3] = Nearest_(B, {a, b});
  auto d0 = Norm(X0 - a);
  auto d1 = Norm(X1 - b);
  auto d2 = Norm(X2 - A);
  auto d3 = Norm(X3 - B);
  if (d0 <= d1 && d0 <= d2 && d0 <= d3)
    return {1., t0};
  else if (d1 <= d2 && d1 <= d3)
    return {0., t1};
  else if (d2 <= d3)
    return {t2, 1.};
  else
    return {t3, 0.};
};
inline Tddd Nearest(const Tddd& X, const T2Tddd& ab) { return std::get<1>(Nearest_(X, ab)); };
inline std::tuple<double, double, Tddd> NearestXOnPlane_(const Tddd& X, const T3Tddd& abc) {
  // アンダースコアがついているものはパラメタも返す
  const auto [a, b, c] = abc;
  // use SolveLinearSystem
  std::array<double, 3> ans;
  // lapack_lu(ans, T3Tddd{a - c, b - c, Cross(a - c, b - c)}, X - c);
  // X.A = b
  // Solve(ans, T3Tddd{a - c, b - c, Cross(a - c, b - c)}, X - c);
  // lapack_lu lu(ans, T3Tddd{a - c, b - c, Cross(a - c, b - c)}, X - c);

  // lapack_svd_solve(Transpose(T3Tddd{a - c, b - c, Cross(a - c, b - c)}), ans, X - c);
  lapack_svd_solve(ans, T3Tddd{a - c, b - c, Cross(a - c, b - c)}, X - c);

  const auto [t0, t1, alpah] = ans;
  return {t0, t1, a * t0 + b * t1 + c * (1. - t0 - t1)};
};
inline Tddd NearestXOnPlane(const Tddd& X, const T3Tddd& abc) {
  // アンダースコアがついているものはパラメタも返す
  return std::get<2>(NearestXOnPlane_(X, abc));
};

// \label{Nearest_}
// inline std::tuple<double, double, Tddd /*to nearest*/, Tddd /*normal vec*/> Nearest_(const Tddd &X, const T3Tddd &abc) {
//   /* ----------------------------------- 修正前 ---------------------------------- */
//   // double T0, T1;
//   // const auto [t0, t1, XOnPlane] = NearestXOnPlane_(X, abc);
//   // //! a*t0 + b*t1 + c*(1-t0-t1)
//   // const auto is_inside_of_triangle_cylinder = Between(t0, array_0_1) && Between(t1, array_0_1) && Between(t0 + t1, array_0_1);
//   // const auto [a, b, c] = abc;
//   // auto [u, X0] = Nearest_(X, T2Tddd{a, b});
//   // //! a*u + b*(1-u)
//   // Tddd ret;
//   // if (is_inside_of_triangle_cylinder && (NormSquared(XOnPlane - X) < NormSquared(X0 - X))) {
//   //    ret = XOnPlane;
//   //    T0 = t0;
//   //    T1 = t1;
//   //    return {T0, T1, ret};
//   // } else {
//   //    ret = X0;
//   //    T0 = u;
//   //    T1 = 1 - u;  //! a*u + b*(1-u) + c * 0
//   // }
//   // auto [v, X1] = Nearest_(X, T2Tddd{b, c});
//   // //! b*v + c*(1-v)
//   // if (NormSquared(ret - X) > NormSquared(X1 - X)) {
//   //    ret = X1;
//   //    T0 = 0;
//   //    T1 = v;  //! a*0 + b*v + c * (1-v)
//   // }
//   // auto [w, X2] = Nearest_(X, T2Tddd{c, a});
//   // //! c*w + a*(1-w)
//   // if (NormSquared(ret - X) > NormSquared(X2 - X)) {
//   //    ret = X2;
//   //    T0 = 1 - w;
//   //    T1 = 0;  //! a*(1-w) + b*0 + c * w
//   // }
//   // return {T0, T1, ret};

//   /* ----------------------------------- 修正その２ ---------------------------------- */
//   /*
//   ２次元に置き換える方法
//   */

//   const auto [a, b, c] = abc;
//   const auto X_a = X - a;
//   const auto b_a = b - a;
//   const auto c_a = c - a;
//   const auto x = Normalize(b_a);
//   const auto z = Normalize(Cross(b_a, c_a));
//   const auto y = Cross(z, x);
//   const double Xx = Dot(X_a, x);
//   const double Xy = Dot(X_a, y);
//   const double Bx = Dot(b_a, x);
//   const double By = Dot(b_a, y);
//   const double Cx = Dot(c_a, x);
//   const double Cy = Dot(c_a, y);
//   //   const double den = By * Cx - Bx * Cy;
//   //   const double t0 = (By * Cx - Bx * Cy + (-By + Cy) * Xx + (Bx - Cx) * Xy) / den;
//   //   const double t1 = -(Cy * Xx - Cx * Xy) / den;
//   const double den = std::fma(By, Cx, -Bx * Cy);
//   const double t0 = std::fma(By, Cx, std::fma(-Bx, Cy, std::fma(-By + Cy, Xx, (Bx - Cx) * Xy))) / den;
//   const double t1 = -std::fma(Cy, Xx, -Cx * Xy) / den;

//   if (0 >= t0 && 0 >= t1)
//     return {0., 0., c, z};
//   else if (0 >= t0 && 1 <= t1)
//     return {0., 1., b, z};
//   else if (1 <= t0 && 0 >= t1)
//     return {1., 0., a, z};
//   else if (0 <= t0 && t0 <= 1. && 0 <= t1 && t1 <= 1. - t0)
//     return {t0, t1, a * t0 + b * t1 + c * (1. - t0 - t1), z};
//   else if (0 >= t0) {
//     auto [t1, closestX] = Nearest_(X, T2Tddd{b, c});
//     return {0., t1, closestX, z};
//   } else if (0 >= t1) {
//     auto [t0, closestX] = Nearest_(X, T2Tddd{a, c});
//     return {t0, 0., closestX, z};
//   } else {
//     auto [t, closestX] = Nearest_(X, T2Tddd{a, b});
//     return {t, 1. - t, closestX, z};
//   }

//   /* ----------------------------------- 修正後 ---------------------------------- */

//   /*
//         t1
//         A
//         |      |        t1>= 1
//      -- b -----+------  t1 = 1
//         | \    |
//         |   \  |
//         |     \|
//      -- c ---- a -----  t1 = 0
//         |      |        t1<=0
//   t0<=0 |      | 1<=t0
//       t0=0    t0=1
//   */

//   //! パラメタをチェックして，三角柱にあるかどうかをチェックせずとも，最近点を求めることができる．
//   // auto [a, b, c] = abc;
//   // double t0, t1, alpah;
//   // const double eps = 1E-5;
//   // if (std::abs(a[2] - b[2]) < eps && std::abs(c[2] - b[2]) < eps) {
//   //    // alpha = a[2];
//   //    Tdd ans;
//   //    lapack_svd_solve(ans, T2Tdd{Tdd{a[0] - c[0], a[1] - c[1]}, Tdd{b[0] - c[0], b[1] - c[1]}}, Tdd{X[0] - c[0], X[1] - c[1]});
//   //    t0 = ans[0];
//   //    t1 = ans[1];
//   // } else {
//   //    Tddd ans;
//   //    // lapack_svd_solve(ans, T3Tddd{a - c, b - c, Cross(a - c, b - c)}, X - c);
//   //    lapack_svd_solve(ans, T3Tddd{a - c, b - c, Cross(a - c, b - c)}, X - c);
//   //    t0 = ans[0];
//   //    t1 = ans[1];
//   //    alpah = ans[2];
//   // }
//   // // lapack_lu(ans, T3Tddd{a - c, b - c, Cross(a - c, b - c)}, X - c);
//   // if (0 >= t0 && 0 >= t1)
//   //    return {0., 0., c};
//   // else if (0 >= t0 && 1 <= t1)
//   //    return {0., 1., b};
//   // else if (1 <= t0 && 0 >= t1)
//   //    return {1., 0., a};
//   // else if (0 <= t0 && t0 <= 1. && 0 <= t1 && t1 <= 1. - t0)
//   //    return {t0, t1, a * t0 + b * t1 + c * (1. - t0 - t1)};
//   // else if (0 >= t0) {
//   //    auto [t1, closestX] = Nearest_(X, T2Tddd{b, c});
//   //    return {0., t1, closestX};
//   // } else if (0 >= t1) {
//   //    auto [t0, closestX] = Nearest_(X, T2Tddd{a, c});
//   //    return {t0, 0., closestX};
//   // } else {
//   //    auto [t, closestX] = Nearest_(X, T2Tddd{a, b});
//   //    return {t, 1. - t, closestX};
//   // }

//   // if (0 >= t0) {
//   //    auto [t1, closestX] = Nearest_(X, T2Tddd{b, c});
//   //    return {0., t1, closestX};
//   // } else if (0 >= t1) {
//   //    auto [t0, closestX] = Nearest_(X, T2Tddd{a, c});
//   //    return {t0, 0., closestX};
//   // } else if (t0 <= 1. && t1 <= 1. && t0 + t1 <= 1.) {
//   //    return {t0, t1, a * t0 + b * t1 + c * (1. - t0 - t1)};
//   // } else {
//   //    auto [t, closestX] = Nearest_(X, T2Tddd{a, b});
//   //    return {t, 1. - t, closestX};
//   // }
// };

// ...existing code...
// \label{Nearest_}
inline std::tuple<double, double, Tddd /*to nearest*/, Tddd /*normal vec*/> Nearest_(const Tddd& X, const T3Tddd& abc) {
  const auto [a, b, c] = abc;

  // 1. エッジベクトルと法線の計算
  const auto edge1 = b - a;
  const auto edge2 = c - a;
  const auto normal = Cross(edge1, edge2);
  const double det = Dot(normal, normal); // det = |normal|^2

  // 縮退三角形（面積がほぼ0）の対策
  if (det < 1.0e-20) {
    // 3辺（線分）との距離を比較して最も近いものを返す
    auto [ta, Pa] = Nearest_(X, T2Tddd{a, b}); // ta is weight for a
    auto [tb, Pb] = Nearest_(X, T2Tddd{b, c}); // tb is weight for b
    auto [tc, Pc] = Nearest_(X, T2Tddd{c, a}); // tc is weight for c
    double da = NormSquared(X - Pa);
    double db = NormSquared(X - Pb);
    double dc = NormSquared(X - Pc);

    if (da <= db && da <= dc)
      return {ta, 1.0 - ta, Pa, {0., 0., 0.}};
    if (db <= dc)
      return {0.0, tb, Pb, {0., 0., 0.}};
    return {1.0 - tc, 0.0, Pc, {0., 0., 0.}};
  }

  const auto unitNormal = normal * (1.0 / std::sqrt(det));

  // 2. 点Xを三角形を含む平面に投影 (Xproj)
  const auto diff = X - a;
  const double distToPlane = Dot(diff, unitNormal);
  const auto Xproj = X - distToPlane * unitNormal;

  // 3. 平面上の点 Xproj の重心座標 (u, v, w) を計算
  // P = a + u * edge1 + v * edge2  =>  P = (1-u-v)*a + u*b + v*c
  // 連立方程式:
  //   (edge1.edge1) * u + (edge1.edge2) * v = (Xproj-a).edge1
  //   (edge1.edge2) * u + (edge2.edge2) * v = (Xproj-a).edge2
  // ここで (Xproj-a).edge = (diff - dist*N).edge = diff.edge (Nとedgeは直交)

  const double d00 = Dot(edge1, edge1);
  const double d01 = Dot(edge1, edge2);
  const double d11 = Dot(edge2, edge2);
  const double b0 = Dot(diff, edge1);
  const double b1 = Dot(diff, edge2);

  // det = d00*d11 - d01*d01 と等価
  const double u = (d11 * b0 - d01 * b1) / det;
  const double v = (d00 * b1 - d01 * b0) / det;
  const double w = 1.0 - u - v;

  // 4. 内部判定: 重心座標がすべて正なら三角形の内部
  if (u >= 0.0 && v >= 0.0 && w >= 0.0) {
    // 戻り値: {weight_a, weight_b, Point, Normal}
    // P = w*a + u*b + v*c
    return {w, u, Xproj, unitNormal};
  }

  // 5. 外部判定: 3つの辺（線分）への最短距離を比較
  // 領域判定ロジック（Voronoi Region）は複雑でバグを生みやすいため、
  // 計算コストは僅かに増えるが、明示的に3辺を調べる方が圧倒的に堅牢。

  // Segment AB: Nearest_ returns weight for first point (a)
  auto [t_ab, P_ab] = Nearest_(X, T2Tddd{a, b});
  // Segment BC: Nearest_ returns weight for first point (b)
  auto [t_bc, P_bc] = Nearest_(X, T2Tddd{b, c});
  // Segment CA: Nearest_ returns weight for first point (c)
  auto [t_ca, P_ca] = Nearest_(X, T2Tddd{c, a});

  double d_ab = NormSquared(X - P_ab);
  double d_bc = NormSquared(X - P_bc);
  double d_ca = NormSquared(X - P_ca);

  if (d_ab <= d_bc && d_ab <= d_ca) {
    return {t_ab, 1.0 - t_ab, P_ab, unitNormal};
  } else if (d_bc <= d_ca) {
    // P_bc = t_bc * b + (1-t_bc) * c
    // weights: a=0, b=t_bc
    return {0.0, t_bc, P_bc, unitNormal};
  } else {
    // P_ca = t_ca * c + (1-t_ca) * a
    // weights: a=1-t_ca, b=0
    return {1.0 - t_ca, 0.0, P_ca, unitNormal};
  }
};
// ...existing code...

inline Tddd Nearest(const Tddd& X, const T3Tddd& abc) { return std::get<2>(Nearest_(X, abc)); };

inline Tddd Nearest(const Tddd& X, const std::vector<T3Tddd>& ABC) {
  double distance = 1E+20, tmp;
  Tddd near, ret;
  for (const auto& abc : ABC) {
    if (distance > (tmp = Norm(X - (near = Nearest(X, abc))))) {
      distance = tmp;
      ret = near;
    }
  }
  return ret;
};

inline Tddd Nearest(const Tddd& X, const T3Tdd& minmax3) {
  double distance = 1E+20, tmp;
  Tddd near, ret;
  CoordinateBounds B(minmax3);
  // std::ranges::for_each((T12T3Tddd)(B), [&](const auto &abc) {
  //    if (distance > (tmp = Norm(X - (near = Nearest(X, abc))))) {
  //       distance = tmp;
  //       ret = near;
  //    }
  // });

  for (const auto& abc : (T12T3Tddd)(B)) {
    if (distance > (tmp = Norm(X - (near = Nearest(X, abc))))) {
      distance = tmp;
      ret = near;
    }
  };

  return ret;
};

inline Tddd Nearest(const Tddd& X, const Tddd& Y) { return Y; };

inline T4d approximateNearest(const T3Tddd& XYZ, const T3Tddd& ABC) {

  auto t0t1_init = std::array<double, 2>{1 / 3., 1 / 3.};

  auto gradF = [&](const std::array<double, 2>& t0t1) -> std::array<double, 2> {
    auto [t0, t1] = t0t1;
    auto f_t0 = (XYZ[0] - XYZ[2]);
    auto f_t1 = (XYZ[1] - XYZ[2]);
    Tddd XonA = (XYZ[0] - XYZ[2]) * t0 + (XYZ[1] - XYZ[2]) * t1 + XYZ[2];
    auto [s0, s1, X, _] = Nearest_(XonA, ABC);
    auto f = XonA - X;
    return {Dot(f, f_t0), Dot(f, f_t1)};
  };

  auto s0s1 = [&](double t0, double t1) -> std::array<double, 2> {
    Tddd XonA = (XYZ[0] - XYZ[2]) * t0 + (XYZ[1] - XYZ[2]) * t1 + XYZ[2];
    auto [s0, s1, X, _] = Nearest_(XonA, ABC);
    return {s0, s1};
  };

  BroydenMethod<std::array<double, 2>> Broyden(t0t1_init, t0t1_init + 1E-10);
  double min_distance = 1E+20;
  for (auto i = 0; i < 30; ++i) {
    auto t0t1_pre = Broyden.X;
    Broyden.updateBFGS(gradF(Broyden.X), gradF(Broyden.X - Broyden.dX), 0.1);
    if (i % 2 == 0) {
      Broyden.X[0] = std::clamp(Broyden.X[0], 0., 1.);
      Broyden.X[1] = std::clamp(Broyden.X[1], 0., 1. - Broyden.X[0]);
    } else {
      Broyden.X[1] = std::clamp(Broyden.X[1], 0., 1.);
      Broyden.X[0] = std::clamp(Broyden.X[0], 0., 1. - Broyden.X[1]);
    }
    Broyden.dX = Broyden.X - t0t1_pre;
  }
  auto [s0, s1] = s0s1(Broyden.X[0], Broyden.X[1]);
  return {Broyden.X[0], Broyden.X[1], s0, s1};
  /* -------------------------------------------------------------------------- */
  // auto t0t1s0s1_init = std::array<double, 4>{1 / 3., 1 / 3., 1 / 3., 1 / 3.};

  // auto f = [&](const std::array<double, 4> &t0t1s0s1) -> Tddd {
  //    auto [t0, t1, s0, s1] = t0t1s0s1;
  //    Tddd XonA = (XYZ[0] - XYZ[2]) * t0 + (XYZ[1] - XYZ[2]) * t1 + XYZ[2];
  //    Tddd XonB = (ABC[0] - ABC[2]) * s0 + (ABC[1] - ABC[2]) * s1 + ABC[2];
  //    return XonA - XonB;
  // };
  // auto gradF = [&](const std::array<double, 4> &t0t1s0s1) -> std::array<double, 4> {
  //    auto f_t0 = (XYZ[0] - XYZ[2]);
  //    auto f_t1 = (XYZ[1] - XYZ[2]);
  //    auto f_s0 = (ABC[0] - ABC[2]);
  //    auto f_s1 = (ABC[1] - ABC[2]);
  //    auto value = f(t0t1s0s1);
  //    return {Dot(value, f_t0), Dot(value, f_t1), Dot(value, f_s0), Dot(value, f_s1)};
  // };

  // if (Norm(f(t0t1s0s1_init)) == 0.)
  //    return t0t1s0s1_init;

  // BroydenMethod<std::array<double, 4>> Broyden(t0t1s0s1_init, t0t1s0s1_init + 1E-10);
  // double min_distance = 1E+20;
  // T4d ret;
  // for (auto i = 0; i < 30; ++i) {
  //    auto t0t1s0s1_pre = Broyden.X;
  //    Broyden.updateBFGS(gradF(Broyden.X), gradF(Broyden.X - Broyden.dX), 0.1);
  //    if (i % 2 == 0) {
  //       Broyden.X[0] = std::clamp(Broyden.X[0], 0., 1.);
  //       Broyden.X[1] = std::clamp(Broyden.X[1], 0., 1. - Broyden.X[0]);
  //       Broyden.X[2] = std::clamp(Broyden.X[2], 0., 1.);
  //       Broyden.X[3] = std::clamp(Broyden.X[3], 0., 1. - Broyden.X[2]);
  //    } else {
  //       Broyden.X[1] = std::clamp(Broyden.X[1], 0., 1.);
  //       Broyden.X[0] = std::clamp(Broyden.X[0], 0., 1. - Broyden.X[1]);
  //       Broyden.X[3] = std::clamp(Broyden.X[3], 0., 1.);
  //       Broyden.X[2] = std::clamp(Broyden.X[2], 0., 1. - Broyden.X[3]);
  //    }
  //    Broyden.dX = Broyden.X - t0t1s0s1_pre;
  // }

  // return Broyden.X;

  /* -------------------------------------------------------------------------- */

  // auto t0t1s0s1_init = std::array<double, 4>{1 / 3., 1 / 3., 1 / 3., 1 / 3.};

  // auto f = [&](const T3Tddd &A, const T3Tddd &B, const std::array<double, 4> &t0t1s0s1) -> Tddd {
  //    auto [t0, t1, s0, s1] = t0t1s0s1;
  //    Tddd XonA = A[0] * t0 + A[1] * t1 * (1 - t0) + A[2] * (1 - t0 - t1 * (1 - t0));
  //    Tddd XonB = B[0] * s0 + B[1] * s1 * (1 - s0) + B[2] * (1 - s0 - s1 * (1 - s0));
  //    return XonA - XonB;
  // };
  // auto gradF = [&](const T3Tddd &A, const T3Tddd &B, const std::array<double, 4> &t0t1s0s1) -> std::array<double, 4> {
  //    auto [t0, t1, s0, s1] = t0t1s0s1;
  //    auto f_t0 = A[0] - A[1] * t1 + A[2] * (-1 + t1);
  //    auto f_t1 = A[1] * (1 - t0) + A[2] * (-(1 - t0));
  //    auto f_s0 = B[0] - B[1] * s1 + B[2] * (-1 + s1);
  //    auto f_s1 = B[1] * (1 - s0) + B[2] * (-(1 - s0));
  //    auto value = f(A, B, t0t1s0s1);
  //    return {Dot(value, f_t0), Dot(value, f_t1), Dot(value, f_s0), Dot(value, f_s1)};
  // };

  // if (Norm(f(ABC, XYZ, t0t1s0s1_init)) == 0.)
  //    return t0t1s0s1_init;

  // BroydenMethod<std::array<double, 4>> Broyden(t0t1s0s1_init, t0t1s0s1_init + 1E-10);
  // double min_distance = 1E+20;
  // for (auto i = 0; i < 30; ++i) {
  //    auto t0t1s0s1_pre = Broyden.X;
  //    Broyden.updateBFGS(gradF(XYZ, ABC, Broyden.X), gradF(XYZ, ABC, Broyden.X - Broyden.dX), 0.1);
  //    Broyden.X[0] = std::clamp(Broyden.X[0], 0., 1.);
  //    Broyden.X[1] = std::clamp(Broyden.X[1], 0., 1.);
  //    Broyden.X[2] = std::clamp(Broyden.X[2], 0., 1.);
  //    Broyden.X[3] = std::clamp(Broyden.X[3], 0., 1.);
  //    Broyden.dX = Broyden.X - t0t1s0s1_pre;
  //    if (i > 3 && Norm(Broyden.X - t0t1s0s1_pre) == 0.)
  //       break;
  // }

  // return Broyden.X;

  /* -------------------------------------------------------------------------- */
};

//\label{Nearest(const T3Tddd &XYZ, const T3Tddd &ABC)}
inline T2Tddd Nearest(const T3Tddd& XYZ, const T3Tddd& ABC) {
  Tddd X0, X1, X0_near, X1_near;
  double nearest_distance = 1E+20, t0, t1, T0, T1, t0_min = 0, t0_max = 1, t1_min = 0, t1_max = 1, distance, w = 1;
  double nearest_t0, nearest_t1;
  double div_init;
  int i, j, I, J;
  nearest_t0 = 1. / 3.;
  nearest_t1 = 1. / 3.;
  const double N = 3.;
  /* -------------------------------------------------------------------------- */
  for (int step = 0; step < 30; ++step) {
    for (i = 0; i <= N; ++i) {
      t0 = i * (t0_max - t0_min) / N + t0_min;
      for (j = 0; j <= N; ++j) {
        t1 = j * (t1_max - t1_min) / N + t1_min;
        X0 = XYZ[0] * t0 + XYZ[1] * t1 * (1. - t0) + XYZ[2] * (1. - t0 - t1 * (1. - t0));
        auto [_, __, X1, ___] = Nearest_(X0, ABC);
        if (nearest_distance > (distance = Norm(X0 - X1))) {
          nearest_distance = distance;
          nearest_t0 = t0;
          nearest_t1 = t1;
          X0_near = X0;
          X1_near = X1;
        }
      }
    }
    w *= 0.7;
    t0_min = std::clamp(nearest_t0 - w, 0., 1.);
    t0_max = std::clamp(nearest_t0 + w, 0., 1.);
    t1_min = std::clamp(nearest_t1 - w, 0., 1.);
    t1_max = std::clamp(nearest_t1 + w, 0., 1.);
  }
  return {X0_near, X1_near};
  /* -------------------------------------------------------------------------- */
  // Tddd X0, X1;
  // double nearest_distance = 1E+20, t0, t1, T0, T1, t0_min, t0_max, t1_min, t1_max, T0_min, T0_max, T1_min, T1_max, distance, w = 1;
  // double nearest_t0, nearest_t1, nearest_T0, nearest_T1;
  // double div_init;
  // int i, j, I, J;
  // auto tmp = approximateNearest(XYZ, ABC);
  // nearest_t0 = tmp[0];
  // nearest_t1 = tmp[1];
  // nearest_T0 = tmp[2];
  // nearest_T1 = tmp[3];
  // double N;
  // /* -------------------------------------------------------------------------- */
  // for (int step = 0; step < 10; ++step) {
  //    if (step == 0) {
  //       N = 5.;
  //       w = 1.;
  //    } else {
  //       N = 2.;
  //       w *= 0.75;
  //    }
  //    t0_min = std::clamp(nearest_t0 - w, 0., 1.);
  //    t0_max = std::clamp(nearest_t0 + w, 0., 1.);
  //    t1_min = std::clamp(nearest_t1 - w, 0., 1.);
  //    t1_max = std::clamp(nearest_t1 + w, 0., 1.);

  //    T0_min = std::clamp(nearest_T0 - w, 0., 1.);
  //    T0_max = std::clamp(nearest_T0 + w, 0., 1.);
  //    T1_min = std::clamp(nearest_T1 - w, 0., 1.);
  //    T1_max = std::clamp(nearest_T1 + w, 0., 1.);

  //    for (i = 0; i <= N; ++i) {
  //       t0 = i * (t0_max - t0_min) / N + t0_min;
  //       for (j = 0; j <= N; ++j) {
  //          t1 = std::clamp(j * (t1_max - t1_min) / N + t1_min, 0., 1. - t0);
  //          X0 = XYZ[0] * t0 + XYZ[1] * t1 + XYZ[2] * (1. - t0 - t1);
  //          for (I = 0; I <= N; ++I) {
  //             T0 = I * (T0_max - T0_min) / N + T0_min;
  //             for (J = 0; J <= N; ++J) {
  //                T1 = std::clamp(J * (T1_max - T1_min) / N + T1_min, 0., 1. - T0);
  //                X1 = ABC[0] * T0 + ABC[1] * T1 + ABC[2] * (1. - T0 - T1);
  //                if (nearest_distance > (distance = Norm(X0 - X1))) {
  //                   nearest_distance = distance;
  //                   nearest_t0 = t0;
  //                   nearest_t1 = t1;
  //                   nearest_T0 = T0;
  //                   nearest_T1 = T1;
  //                }
  //             }
  //          }
  //       }
  //    }
  // }
  // return {XYZ[0] * nearest_t0 + XYZ[1] * nearest_t1 + XYZ[2] * (1. - nearest_t0 - nearest_t1),
  //         ABC[0] * nearest_T0 + ABC[1] * nearest_T1 + ABC[2] * (1. - nearest_T0 - nearest_T1)};
  /* -------------------------------------------------------------------------- */
  // for (int step = 0; step < 10; ++step) {
  //    if (step == 0) {
  //       N = 5;
  //       w *= 0.75;
  //    } else {
  //       N = 2;
  //       w *= 0.75;
  //    }
  //    t0_min = std::clamp(nearest_t0 - w, 0., 1.);
  //    t0_max = std::clamp(nearest_t0 + w, 0., 1.);
  //    t1_min = std::clamp(nearest_t1 - w, 0., 1.);
  //    t1_max = std::clamp(nearest_t1 + w, 0., 1.);

  //    T0_min = std::clamp(nearest_T0 - w, 0., 1.);
  //    T0_max = std::clamp(nearest_T0 + w, 0., 1.);
  //    T1_min = std::clamp(nearest_T1 - w, 0., 1.);
  //    T1_max = std::clamp(nearest_T1 + w, 0., 1.);

  //    for (i = 0; i <= N; ++i) {
  //       t0 = i * (t0_max - t0_min) / N + t0_min;
  //       for (j = 0; j <= N; ++j) {
  //          t1 = j * (t1_max - t1_min) / N + t1_min;
  //          X0 = XYZ[0] * t0 + XYZ[1] * t1 * (1 - t0) + XYZ[2] * (1. - t0 - t1 * (1 - t0));
  //          for (I = 0; I <= N; ++I) {
  //             T0 = I * (T0_max - T0_min) / N + T0_min;
  //             for (J = 0; J <= N; ++J) {
  //                T1 = J * (T1_max - T1_min) / N + T1_min;
  //                X1 = ABC[0] * T0 + ABC[1] * T1 * (1 - T0) + ABC[2] * (1. - T0 - T1 * (1 - T0));
  //                if (nearest_distance > (distance = Norm(X0 - X1))) {
  //                   nearest_distance = distance;
  //                   nearest_t0 = t0;
  //                   nearest_t1 = t1;
  //                   nearest_T0 = T0;
  //                   nearest_T1 = T1;
  //                }
  //             }
  //          }
  //       }
  //    }
  // }
  // return {XYZ[0] * nearest_t0 + XYZ[1] * nearest_t1 * (1 - nearest_t0) + XYZ[2] * (1. - nearest_t0 - nearest_t1 * (1 - nearest_t0)),
  //         ABC[0] * nearest_T0 + ABC[1] * nearest_T1 * (1 - nearest_T0) + ABC[2] * (1. - nearest_T0 - nearest_T1 * (1 - nearest_T0))};
};

// // ２つの三角形の最短ベクトルを求める
// T2Tddd Nearest(const T3Tddd &A, const T3Tddd &B) {

// };

inline double scalefactorToReach(const T2Tddd& line, const T3Tddd& triangle) {
  auto [a, b] = line;
  auto [p0, p1, p2] = triangle;
  // オーダーが匹敵する物を選ぶ
  double log_b_a = std::log10(Norm(b - a));
  double diff0 = std::abs(std::log10(Norm(p0 - a) - log_b_a));
  double diff1 = std::abs(std::log10(Norm(p1 - a) - log_b_a));
  double diff2 = std::abs(std::log10(Norm(p2 - a) - log_b_a));
  Tddd n = TriangleNormal(p0, p1, p2);
  if (diff0 < diff1 && diff0 < diff2)
    return Dot(p0 - a, n) / Dot(b - a, n);
  else if (diff1 < diff0 && diff1 < diff2)
    return Dot(p1 - a, n) / Dot(b - a, n);
  else
    return Dot(p2 - a, n) / Dot(b - a, n);
};

// b* -------------------------------------------------------------------------- */
// b*                                 IntersectQ                                 */
// b* -------------------------------------------------------------------------- */
// \label{basic_geometry:IntersectQ}
//! cube - cube
// 暗黙の変換は，予期しない変換を行なってしまう場合が多発するので，ここでは型変換が起きないよう，CoordinateBounds同士の干渉は行わないようにする．
// bool IntersectQ(const T3Tdd &b0, const T3Tdd &b1) {
//    return !((std::get<0>(std::get<0>(b0)) > std::get<0>(std::get<0>(b1)) && std::get<0>(std::get<0>(b0)) > std::get<1>(std::get<0>(b1)) /*1のxの最大最小が，0のxの最小よりも小さい*/) ||
//             (std::get<0>(std::get<1>(b0)) > std::get<0>(std::get<1>(b1)) && std::get<0>(std::get<1>(b0)) > std::get<1>(std::get<1>(b1)) /*1のyの最大最小が，0のyの最小よりも小さい*/) ||
//             (std::get<0>(std::get<2>(b0)) > std::get<0>(std::get<2>(b1)) && std::get<0>(std::get<2>(b0)) > std::get<1>(std::get<2>(b1)) /*1のzの最大最小が，0のzの最小よりも小さい*/) ||
//             (std::get<1>(std::get<0>(b0)) < std::get<0>(std::get<0>(b1)) && std::get<1>(std::get<0>(b0)) < std::get<1>(std::get<0>(b1)) /*1のxの最大最小が，0のxの最大よりも大きい*/) ||
//             (std::get<1>(std::get<1>(b0)) < std::get<0>(std::get<1>(b1)) && std::get<1>(std::get<1>(b0)) < std::get<1>(std::get<1>(b1)) /*1のyの最大最小が，0のyの最大よりも大きい*/) ||
//             (std::get<1>(std::get<2>(b0)) < std::get<0>(std::get<2>(b1)) && std::get<1>(std::get<2>(b0)) < std::get<1>(std::get<2>(b1)) /*1のzの最大最小が，0のzの最大よりも大きい*/) /*これがtrueの場合，逆にhitなし*/);
// };

inline bool IntersectQ(const T3Tdd& b0, const T3Tdd& b1) {
  auto [x0min, x0max] = b0[0];
  auto [y0min, y0max] = b0[1];
  auto [z0min, z0max] = b0[2];
  auto [x1min, x1max] = b1[0];
  auto [y1min, y1max] = b1[1];
  auto [z1min, z1max] = b1[2];
  if (x1max < x0min || y1max < y0min || z1max < z0min || x1min > x0max || y1min > y0max || z1min > z0max)
    return false;
  else
    return true;

  // return !(std::get<1>(b0[0]) < std::get<0>(b1[0]) || std::get<1>(b1[0]) < std::get<0>(b0[0]) ||  // x軸
  //          std::get<1>(b0[1]) < std::get<0>(b1[1]) || std::get<1>(b1[1]) < std::get<0>(b0[1]) ||  // y軸
  //          std::get<1>(b0[2]) < std::get<0>(b1[2]) || std::get<1>(b1[2]) < std::get<0>(b0[2]));   // z軸
}

//! cube - point
inline bool IntersectQ(const Tddd& X, const T3Tdd& minmax3) { return !((std::get<0>(X) < std::get<0>(std::get<0>(minmax3)) || std::get<1>(std::get<0>(minmax3)) < std::get<0>(X)) || (std::get<1>(X) < std::get<0>(std::get<1>(minmax3)) || std::get<1>(std::get<1>(minmax3)) < std::get<1>(X)) || (std::get<2>(X) < std::get<0>(std::get<2>(minmax3)) || std::get<1>(std::get<2>(minmax3)) < std::get<2>(X))); };
inline bool IntersectQ(const T3Tdd& minmax3, const Tddd& X) { return !((std::get<0>(X) < std::get<0>(std::get<0>(minmax3)) || std::get<1>(std::get<0>(minmax3)) < std::get<0>(X)) || (std::get<1>(X) < std::get<0>(std::get<1>(minmax3)) || std::get<1>(std::get<1>(minmax3)) < std::get<1>(X)) || (std::get<2>(X) < std::get<0>(std::get<2>(minmax3)) || std::get<1>(std::get<2>(minmax3)) < std::get<2>(X))); };
//! sphere - point
inline bool IntersectQ(const Tddd& center, const double& r, const Tddd& a) { return Norm(a - center) <= r; };
inline bool IntersectQ(const Sphere& sp, const Tddd& a) { return Norm(a - sp.center) <= sp.radius; };
inline bool IntersectQ(const Tddd& a, const Sphere& sp) { return Norm(a - sp.center) <= sp.radius; };
//! sphere - cube
inline bool IntersectQ(const Tddd& X, const double& r, const T3Tdd& minmax3) { return IntersectQ(X, minmax3) || (r >= Norm(Nearest(X, minmax3) - X)); };
inline bool IntersectQ(const Sphere& s, const T3Tdd& minmax3) { return IntersectQ(s.center, s.radius, minmax3); };
inline bool IntersectQ(const T3Tdd& minmax3, const Sphere& s) { return IntersectQ(s.center, s.radius, minmax3); };
//! sphere - sphere
inline bool IntersectQ(const Tddd& x0, const double r0, const Tddd& x1, const double r1) { return Norm(x0 - x1) <= (r0 + r1); };
inline bool IntersectQ(const Sphere& s0, const Sphere& s1) { return Norm(s0.X - s1.X) >= (s0.radius + s1.radius); };
//! sphere - line
inline bool IntersectQ(const Tddd& center, const double radius, const T2Tddd& ab) {
  const auto a = std::get<0>(ab) - center;
  const auto b = std::get<1>(ab) - center;
  if (Norm(a) < radius || Norm(b) < radius)
    return true;
  const Tddd a2b = b - a;
  const double t = Dot(-a, a2b) / Dot(a2b, a2b);
  return (0. <= t && t <= 1. && Norm(a + a2b * t) < radius);
};
inline bool IntersectQ(const Tddd& center, const double radius, const T6T2Tddd& ab) { return IntersectQ(center, radius, std::get<0>(ab)) || IntersectQ(center, radius, std::get<1>(ab)) || IntersectQ(center, radius, std::get<2>(ab)) || IntersectQ(center, radius, std::get<3>(ab)) || IntersectQ(center, radius, std::get<4>(ab)) || IntersectQ(center, radius, std::get<5>(ab)); };

inline bool IntersectQ(const Sphere& sp, const T6T2Tddd& ab) { return IntersectQ(sp.center, sp.radius, std::get<0>(ab)) || IntersectQ(sp.center, sp.radius, std::get<1>(ab)) || IntersectQ(sp.center, sp.radius, std::get<2>(ab)) || IntersectQ(sp.center, sp.radius, std::get<3>(ab)) || IntersectQ(sp.center, sp.radius, std::get<4>(ab)) || IntersectQ(sp.center, sp.radius, std::get<5>(ab)); };

inline bool IntersectQ(const Sphere& s, const T2Tddd& ab) { return IntersectQ(s.center, s.radius, ab); };
inline bool IntersectQ(const T2Tddd& ab, const Sphere& sp) { return IntersectQ(sp, ab); };
//! sphere - triangle
inline bool IntersectQ(const Tddd& X, const double r, const T3Tddd& abc) { return r >= Norm(Nearest(X, abc) - X); };
inline bool IntersectQ(const Sphere& sp, const T3Tddd& abc) { return sp.radius > Norm(Nearest(sp.center, abc) - sp.center); };
//! line - triangle
inline std::tuple<bool, Tddd, Tddd> IntersectQ_(const T2Tddd& AB, const T3Tddd& abc) {

  if (IntersectQ(CoordinateBounds(AB).bounds, CoordinateBounds(abc).bounds)) {
    std::array<double, 3> x;
    // Solve(x, T3Tddd{a - c, b - c, A - B}, A - c);
    // lapack_lu lu(x, T3Tddd{a - c, b - c, A - B}, A - c);
    const auto [a, b, c] = abc;
    const auto [A, B] = AB;
    lapack_svd_solve(x, T3Tddd{a - c, b - c, A - B}, A - c);
    const auto [t0, t1, T] = x;
    static const Tdd range = {0., 1.};
    return {Between(T, range) && Between(t0, range) && Between(t1, range) && Between(t0 + t1, range), (1. - T) * A + T * B, {t0, t1, 1. - t0 - t1}};
    /*
    Hit = A + T * (B - A)
    Hit = a * t0 + b * t1 + c * (1-t0-t1) = (a-c) * t0 + (b-c) * t1 + c
    equalize the two equations
    {t0,t1,T} . {a-c, b-c, A-B} = A-c
    */
  } else
    return {false, {0., 0., 0.}, {0., 0., 0.}};
  /*

  $X$は，$A$と$B$を結ぶ直線上にあるとする．この$X$を頂点の重心座標で表す．
  この二つの表現を等しくするような$t_0,t_1,T$を求める．
  結果から，$T$が$0$から$1$の間にあるかどうかで，$X$が線分上にあるかどうかがわかる．
  また，$t_0,t_1,1-t_0-t_1$が$0$から$1$の間にあるかどうかで，$X$が三角形の内部にあるかどうかがわかる．

  $$
  X = (t0,t1,1-t0-t1)\cdot(a-c,b-c,c)\\
  X = A + T\cdot(B-A)
  $$
  */

  /* -------------------------------------------------------------------------- */

  /*

  $X$は，$A$と$B$を結ぶ直線上にあるとする．この$X$を頂点の重心座標で表す．
  この二つの表現を等しくするような$t_0,t_1,T$を求める．
  結果から，$T$が$0$から$1$の間にあるかどうかで，$X$が線分上にあるかどうかがわかる．
  また，$t_0,t_1,1-t_0-t_1$が$0$から$1$の間にあるかどうかで，$X$が三角形の内部にあるかどうかがわかる．

  $$
  X = (t0,t1,1-t0-t1)\cdot(a-c,b-c,c)\\
  X = A + T\cdot(B-A)
  $$
  */
};

inline bool IntersectQ(const T2Tddd& AB, const T3Tddd& abc) {
  const auto [a, b, c] = abc;
  const auto [A, B] = AB;
  // const auto [t0, t1, T] = Dot(A - c, Inverse(T3Tddd{a - c, b - c, A - B}));
  std::array<double, 3> x;
  lapack_lu(x, T3Tddd{a - c, b - c, A - B}, A - c);
  // lapack_svd_solve(x, T3Tddd{a - c, b - c, A - B}, A - c);
  // Solve(x, T3Tddd{a - c, b - c, A - B}, A - c);
  const auto [t0, t1, T] = x;
  static const Tdd range = array_0_1;
  return Between(T, range) && Between(t0, range) && Between(t1, range) && Between(t0 + t1, range);
};
inline bool IntersectQ(const T3Tddd& abc, const T2Tddd& AB) { return IntersectQ(AB, abc); };
//! cube - line
inline bool IntersectQ(const T3Tdd& minmax3, const T2Tddd& AB) {
  const auto [A, B] = AB;
  if (IntersectQ(minmax3, A) || IntersectQ(minmax3, B) || IntersectQ(minmax3, 0.5 * (A + B)))
    return true;
  const auto minmax = Transpose(minmax3);
  const auto [x0, y0, z0] = std::get<0>(minmax);
  const auto [x1, y1, z1] = std::get<1>(minmax);
  const Tddd a = {x0, y0, z0};
  const Tddd b = {x0, y1, z1};
  const Tddd d = {x0, y1, z0};
  const Tddd c = {x1, y1, z0};
  //! --------------------------------- */
  // const auto inv = Inverse(Transpose(T3Tddd{a - d, b - d, c - d}));
  // const Tddd num = Dot(inv, d - B);
  // const Tddd den = Dot(inv, A - B);
  // 変更後
  const auto adbdcd = T3Tddd{a - d, b - d, c - d};
  Tddd num, den;
  // Solve(num, adbdcd, d - B);
  // Solve(den, adbdcd, A - B);
  lapack_lu lu_num(num, adbdcd, d - B);
  // lapack_svd_solve(num, adbdcd, d - B);
  lapack_lu lu_den(den, adbdcd, A - B);
  // lapack_svd_solve(den, adbdcd, A - B);
  //
  auto [int0, int1, int2] = Transpose(T2Tddd{num / den, (1 + num) / den});
  //! --------------------------------- */
  //! 線分の成分が一致する場合や，denが０の場合にに発生する，indeterminateな場合結果を避けるための処理
  const double small = 1E-13;
  // const auto [T0, T1, T2] = Dot(inv, (A + B) * 0.5 - d);
  // 変更後
  Tddd T0T1T2;
  // Solve(T0T1T2, adbdcd, (A + B) * 0.5 - d);
  // lapack_lu lu_T0T1T2(T0T1T2, adbdcd, (A + B) * 0.5 - d);
  lapack_svd svd_T0T1T2(Transpose(adbdcd), T0T1T2, (A + B) * 0.5 - d);
  const auto [T0, T1, T2] = T0T1T2;
  //
  const auto [dx, dy, dz] = A - B;
  if (Norm(dx) < small || std::abs(std::get<0>(den)) < small)
    if (Between(T0, array_0_1))
      int0 = {0, 1};
    else
      return false;
  if (Norm(dy) < small || std::abs(std::get<1>(den)) < small)
    if (Between(T1, array_0_1))
      int1 = {0, 1};
    else
      return false;
  if (Norm(dz) < small || std::abs(std::get<2>(den)) < small)
    if (Between(T2, array_0_1))
      int2 = {0, 1};
    else
      return false;
  if ((std::get<0>(int0) < 0 && std::get<1>(int0) < 0) || (std::get<0>(int0) > 1 && std::get<1>(int0) > 1) || (std::get<0>(int1) < 0 && std::get<1>(int1) < 0) || (std::get<0>(int1) > 1 && std::get<1>(int1) > 1) || (std::get<0>(int2) < 0 && std::get<1>(int2) < 0) || (std::get<0>(int2) > 1 && std::get<1>(int2) > 1))
    return false;
  //! --------------------------------- */
  if (std::get<0>(int0) > std::get<1>(int0))
    Swap(int0);
  if (std::get<0>(int1) > std::get<1>(int1))
    Swap(int1);
  if (std::get<0>(int2) > std::get<1>(int2))
    Swap(int2);
  //! --------------------------------- */
  return Min(Tddd{std::get<1>(int0), std::get<1>(int1), std::get<1>(int2)}) >= Max(Tddd{std::get<0>(int0), std::get<0>(int1), std::get<0>(int2)});
};
inline bool IntersectQ(const T2Tddd& AB, const T3Tdd& minmax3) { return IntersectQ(minmax3, AB); };

//! cube - tringle
inline bool IntersectQ(const T3Tdd& minmax3, const T3Tddd& abc) {
  // チェックの短縮化が必要
  CoordinateBounds boundTri(abc);
  if (IntersectQ(minmax3, boundTri.bounds)) {
    if (IntersectQ(minmax3, std::get<0>(abc)) || IntersectQ(minmax3, std::get<1>(abc)) || IntersectQ(minmax3, std::get<2>(abc)) || IntersectQ(minmax3, Mean(abc)))
      return true;
    auto [a, b, c] = abc;
    return IntersectQ(minmax3, T2Tddd{a, b}) || IntersectQ(minmax3, T2Tddd{b, c}) || IntersectQ(minmax3, T2Tddd{c, a}) || std::ranges::any_of((T12T2Tddd)CoordinateBounds(minmax3), [&](const auto& AB) { return IntersectQ(abc, AB); });
  } else
    return false;
};

inline Tddd XonTriangle(const T3Tddd& abc, const T2Tddd& AB) {
  auto [a, b, c] = abc;
  auto [A, B] = AB;
  // auto [t0, t1, T] = Dot(Inverse(Transpose(T3Tddd{a - c, b - c, A - B})), A - c);
  // auto [t0, t1, T] = Dot(A - c, Inverse(T3Tddd{a - c, b - c, A - B}));
  std::array<double, 3> ans;
  lapack_lu(ans, T3Tddd{a - c, b - c, A - B}, A - c);
  // lapack_svd_solve(ans, T3Tddd{a - c, b - c, A - B}, A - c);
  // Solve(ans, T3Tddd{a - c, b - c, A - B}, A - c);
  const auto [t0, t1, T] = ans;
  //
  return A + (B - A) * T;
};

//! tetrahedron - point
inline bool IntersectQ(const T4Tddd& abcd, Tddd X) {
  //@ barycentric coordinates
  //
  // | a0, b0, c0, d0 | | t0 |   | x |
  // | a1, b1, c1, d1 | | t1 |   | y |
  // | a2, b2, c2, d2 | | t2 | = | z |
  // |  1,  1,  1,  1 | | t3 |   | 1 |
  //
  //                    | a0, a1, a2,  1 |
  // | t0, t1, t2, t3 | | b0, b1, b1,  1 |  = | x, y, y, z |
  //                    | c0, c2, c2,  1 |
  //                    | d0, d1, d2,  1 |
  //
  // -> Dot({t0,t1,t2,t3}, {{a,1},{b,1},{c,1},{d,1}}) = {X,1}
  // {a,1} and {X,1} mean {a0, a1, a2, 1} and {x, y, z, 1}
  //
  const auto [a, b, c, d] = abcd;
  // const auto [t0, t1, t2, t3] = Dot(T4d{std::get<0>(X), std::get<1>(X), std::get<2>(X), 1.},
  //                                   Inverse(T4T4d{{{std::get<0>(a), std::get<1>(a), std::get<2>(a), 1.},
  //                                                  {std::get<0>(b), std::get<1>(b), std::get<2>(b), 1.},
  //                                                  {std::get<0>(c), std::get<1>(c), std::get<2>(c), 1.},
  //                                                  {std::get<0>(d), std::get<1>(d), std::get<2>(d), 1.}}}));

  // use SolveLinearSystem
  std::array<double, 4> x;
  // lapack_lu(x, T4T4d{{{std::get<0>(a), std::get<1>(a), std::get<2>(a), 1.}, {std::get<0>(b), std::get<1>(b), std::get<2>(b), 1.}, {std::get<0>(c), std::get<1>(c), std::get<2>(c), 1.}, {std::get<0>(d), std::get<1>(d), std::get<2>(d), 1.}}}, T4d{std::get<0>(X), std::get<1>(X), std::get<2>(X), 1.});
  lapack_svd_solve(x, T4T4d{{{std::get<0>(a), std::get<1>(a), std::get<2>(a), 1.}, {std::get<0>(b), std::get<1>(b), std::get<2>(b), 1.}, {std::get<0>(c), std::get<1>(c), std::get<2>(c), 1.}, {std::get<0>(d), std::get<1>(d), std::get<2>(d), 1.}}}, T4d{std::get<0>(X), std::get<1>(X), std::get<2>(X), 1.});
  const auto [t0, t1, t2, t3] = x;

  //
  //
  // auto mean = Mean(abcd);
  // auto [a, b, c, d] = abcd - mean;
  // X -= Mean(abcd);
  // const auto [t0, t1, t2, t3] = Dot(T4d{std::get<0>(X), std::get<1>(X), std::get<2>(X), 1.},
  //                                   Inverse(T4T4d{{std::get<0>(a), std::get<1>(a), std::get<2>(a), 1.},
  //                                                 {std::get<0>(b), std::get<1>(b), std::get<2>(b), 1.},
  //                                                 {std::get<0>(c), std::get<1>(c), std::get<2>(c), 1.},
  //                                                 {std::get<0>(d), std::get<1>(d), std::get<2>(d), 1.}}));

  return (Between(t0, array_0_1) && Between(t1, {0., 1. - t0}) && Between(t2, {0., 1. - t0 - t1}));
}

//! tetrahedron - sphere
inline bool IntersectQ(const Tddd& center, const double r, const T4Tddd& abcd) {
  if (IntersectQ(abcd, center))
    return true;
  else {
    const auto [a, b, c, d] = abcd;
    return (Norm(Nearest(center, std::vector<T3Tddd>{T3Tddd{a, b, c}, T3Tddd{a, b, d}, T3Tddd{a, c, d}, T3Tddd{b, c, d}}) - center) < r);
  }
}

inline bool IntersectQ(const Sphere& sp, const Tetrahedron& t) { return IntersectQ(sp.center, sp.radius, t.vertices); }
inline bool IntersectQ(const Tetrahedron& t, const Sphere& sp) { return IntersectQ(sp.center, sp.radius, t.vertices); }

//! tetrahedron - line
inline bool IntersectQ(const T4Tddd& abcd, const T2Tddd& AB) {
  auto [a, b, c, d] = abcd;
  auto [A, B] = AB;
  return IntersectQ(abcd, A) || IntersectQ(abcd, B) || IntersectQ(T3Tddd{{a, b, c}}, AB) || IntersectQ(T3Tddd{{a, b, d}}, AB) || IntersectQ(T3Tddd{{a, c, d}}, AB) || IntersectQ(T3Tddd{{b, c, d}}, AB);
};

inline bool IntersectQ(const Tetrahedron& Tet, const T2Tddd& AB) {
  if (IntersectQ(Tet.bounds, AB))
    return IntersectQ(Tet.vertices, AB);
  else
    return false;
};

inline bool IntersectQ(const Tetrahedron& Tet0, const Tetrahedron& Tet1) {
  if (IntersectQ(Tet0.bounds, Tet1.bounds)) {
    if (IntersectQ(Tet0.incenter, Tet0.inradius, Tet1.incenter, Tet1.inradius))
      return true;
    return std::ranges::any_of((T6T2Tddd)(Tet0), [&](const auto& ab) { return IntersectQ(Tet1, ab); }) ||

           std::ranges::any_of((T6T2Tddd)(Tet1), [&](const auto& AB) { return IntersectQ(Tet0, AB); });
  } else
    return false;
};

inline bool IntersectQ(const T4Tddd& abcd, const T4Tddd& ABCD) { return IntersectQ(Tetrahedron(abcd), Tetrahedron(ABCD)); };

inline Tddd t0_t1_alpha(const T3Tddd& p0p1p2, const Tddd& X) {
  //@ ３点の張る面　と　１点　の関係を調べる
  // a = (p0,p1,p2).(t0,t1,1-t0,t1)は，三角形が張る面上で，Xに最も近い点
  // この位置aから，alpha*nだけ移動した位置にXがある．n=(nx, ny, nz)は，三角形がつくる単位法線ベクトル
  auto [p0, p1, p2] = p0p1p2;
  auto [nx, ny, nz] = Normalize(Cross(p1 - p0, p2 - p0)); // 三角形がつくる単位法線ベクトル
  auto [p02x, p02y, p02z] = p0 - p2;
  auto [p12x, p12y, p12z] = p1 - p2;
  return Dot(X - p2, T3Tddd{{{nz * p12y - ny * p12z, -(nz * p02y) + ny * p02z, p02z * p12y - p02y * p12z}, {-(nz * p12x) + nx * p12z, nz * p02x - nx * p02z, -(p02z * p12x) + p02x * p12z}, {ny * p12x - nx * p12y, -(ny * p02x) + nx * p02y, p02y * p12x - p02x * p12y}}}) / (-(nz * p02y * p12x) + ny * p02z * p12x + nz * p02x * p12y - nx * p02z * p12y - ny * p02x * p12z + nx * p02y * p12z);
};

/* -------------------------------------------------------------------------- */
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using V_i = std::vector<int>;
using VV_i = std::vector<std::vector<int>>;
/* ------------------------------------------------------ */

inline double normalDirDistanceFromTriangle(const T3Tddd& ps, const Tddd& a) { return Dot(TriangleNormal(ps), std::get<0>(ps) - a); };

inline double factorOfVectorToReachTriangle(const Tddd& p0, const Tddd& p1, const Tddd& p2, const Tddd& a, const Tddd& b) {
  // オーダーが匹敵する物を選ぶ
  double log_b_a = std::log10(Norm(b - a));
  double diff0 = std::abs(std::log10(Norm(p0 - a) - log_b_a));
  double diff1 = std::abs(std::log10(Norm(p1 - a) - log_b_a));
  double diff2 = std::abs(std::log10(Norm(p2 - a) - log_b_a));
  auto n = TriangleNormal(p0, p1, p2);
  if (diff0 < diff1 && diff0 < diff2)
    return Dot(p0 - a, n) / Dot(b - a, n);
  else if (diff1 < diff0 && diff1 < diff2)
    return Dot(p1 - a, n) / Dot(b - a, n);
  else
    return Dot(p2 - a, n) / Dot(b - a, n);
};

/* ------------------------------------------------------ */
inline Tddd vectorToInfiniteLine(const Tddd& P, Tddd A, Tddd B) {
  // Tddd BA = B - A, AP = A - P;
  // return AP - BA * Dot(AP, BA) / Dot(BA, BA);
  /* ------------------------------------------------------ */
  B -= A;
  A -= P;
  return A - B * Dot(A, B) / Dot(B, B);
};
inline double distanceToInfiniteLine(const Tddd& P, const Tddd& A, const Tddd& B) { return Norm(vectorToInfiniteLine(P, A, B)); };
/* ------------------------------------------------------ */
inline Tddd pOnSurfaceTuple(const Tddd& p0, const Tddd& p1, const Tddd& p2, const Tddd& a, const Tddd& b) {
  Tddd n = TriangleNormal(p0, p1, p2), b_a = b - a;
  return a + b_a * Dot(p0 - a, n) / Dot(b_a, n); // 分母が0の場合はあり得る
}
inline Tddd pOnSurfaceTuple(const T3Tddd& p0p1p2, const T2Tddd& ab) { return pOnSurfaceTuple(std::get<0>(p0p1p2), std::get<1>(p0p1p2), std::get<2>(p0p1p2), std::get<0>(ab), std::get<1>(ab)); }
/* -------------------------------------------------------------------------- */

inline bool isConvexPolygon(const std::vector<Tddd>& ps, const Tddd& n) {
  auto s = ps.size();
  if (s < 3)
    return false;
  else if (s == 3)
    return true;

  for (auto i = 0; i < s; ++i) {
    auto v0 = ps[i + 1] - ps[i];
    auto v1 = ps[i + 2] - ps[i + 1];
    auto angle = CrossDot(v0, v1, n);
    if (angle <= 1E-13)
      return false; // 符号が変わったらfalse
  }
  return true;
};

inline bool isConcavePolygon(const std::vector<Tddd>& ps, const Tddd& n) { return !isConvexPolygon(ps, n); };

inline bool isConvexPolygon(const std::vector<Tdd>& ps) {
  std::vector<Tddd> Ps(ps.size());
  int i = 0;
  for (const auto& v : ps)
    Ps[i++] = {std::get<0>(v), std::get<1>(v), 0.};
  return isConvexPolygon(Ps, Tddd{0., 0., 1});
};

inline bool isConcavePolygon(const std::vector<Tdd>& ps) { return !isConvexPolygon(ps); };

inline bool isConcavePolygon(const std::vector<Tddd>& ps) {
  auto s = ps.size();
  if (s < 3)
    return false;

  auto v0 = *ps.begin() - *ps.rbegin();
  auto v1 = *(ps.begin() + 1) - *ps.begin();
  auto normal = Cross(v0, v1);
  normal = normal / Norm(normal);
  double angle = MyVectorAngle(v0 /*基準*/, v1, normal);
  if (!isFinite(normal))
    return false;

  for (auto i = 0; i < s; i++) {
    // 0->1 1->2
    angle = angle * MyVectorAngle(ps[i] - ps[(s + i - 1) % s] /*基準*/, ps[(s + i + 1) % s] - ps[i], normal);
    if (!isFinite(angle))
      return false;
    if (angle < 0.)
      return true; // 符号が変わったらtrue
  }
  return false; // 符号が変わらなかったのでfalse
};

/* ------------------------------------------------------ */

template <typename T>
inline double windingNumber(const Tddd& X, const std::vector<std::array<T, 3>>& V_vertices) {
  double ret = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : ret)
#endif
  for (const auto& V : V_vertices)
    ret += SolidAngle_VanOosteromAandStrackeeJ1983(X, ToX(V));
  return ret / (4. * M_PI);
};

template <typename T>
inline double windingNumber(const Tddd& X, const std::vector<T>& V_vertices) {
  double ret = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : ret)
#endif
  for (const auto& V : V_vertices)
    ret += SolidAngle_VanOosteromAandStrackeeJ1983(X, ToX(V));
  return ret / (4. * M_PI);
};

template <>
inline double windingNumber(const Tddd& X, const std::vector<T3Tddd>& V_vertices) {
  double ret = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : ret)
#endif
  for (const auto& V : V_vertices)
    ret += SolidAngle_VanOosteromAandStrackeeJ1983(X, V);
  return ret / (4. * M_PI);
};

template <typename T>
inline T8d windingNumber(const T8Tddd& Xs, const std::vector<T>& V_vertices) {
  T8d ret = {0., 0., 0., 0., 0., 0., 0., 0.};
  // for (const auto &V : V_vertices)
  {
    for_each(ret, Xs, [&](auto& r, auto& X) {
      for (const auto& V : V_vertices)
        r += SolidAngle_VanOosteromAandStrackeeJ1983(X, ToX(V));
    });

    // std::ranges::for_each(ret, Xs, [&](auto &r, const auto &X) {
    //    // r += SolidAngle_VanOosteromAandStrackeeJ1983(X, V);
    //    r += SolidAngle_VanOosteromAandStrackeeJ1983(X, ToX(V));
    // });
  }
  return ret / (4. * M_PI);
};

template <>
inline T8d windingNumber(const T8Tddd& Xs, const std::vector<T3Tddd>& V_vertices) {
  T8d ret = {0., 0., 0., 0., 0., 0., 0, 0.};
  for (const auto& V : V_vertices)
    for_each(ret, Xs, [&](auto& r, auto& X) { r += SolidAngle_VanOosteromAandStrackeeJ1983(X, V); });
  // std::ranges::for_each(ret, Xs, [&](auto &r, const auto &X) { r += SolidAngle_VanOosteromAandStrackeeJ1983(X, V); });
  return ret / (4. * M_PI);
};
template <>
inline T8d windingNumber(const T8Tddd& Xs, const std::vector<Tddd>& V_vertices) { return {0., 0., 0., 0., 0., 0., 0, 0.}; };

inline T8d windingNumber(const T8Tddd& Xs, const std::unordered_set<Tddd>& V_vertices) { return {0., 0., 0., 0., 0., 0., 0, 0.}; };

inline std::vector<double> windingNumber(const std::vector<Tddd>& Xs, const std::vector<T3Tddd>& V_vertices) {
  std::vector<double> ret(Xs.size(), 0.);
  for (auto i = 0; i < Xs.size(); ++i) {
    auto X = Xs[i];
    auto tmp = 0;
    for (const auto& vertices : V_vertices)
      tmp += SolidAngle_VanOosteromAandStrackeeJ1983(X, vertices);
    ret[i] = tmp / (4. * M_PI);
  }
  return ret;
};

inline T8d windingNumber(const T8Tddd& Xs, const std::vector<T3Tddd>& V_vertices) {
  T8d ret = {0., 0., 0., 0., 0., 0., 0, 0.};
  for (const auto& vertices : V_vertices)
    for_each(ret, Xs, [&](auto& r, auto& X) { r += SolidAngle_VanOosteromAandStrackeeJ1983(X, vertices); });
  return ret / (4. * M_PI);
};

inline double windingNumber(const Tddd& X, const std::vector<Tddd>& V_vertices) { return 0.; };

/* -------------------------------------------------------------------------- */

inline bool InsideQ(const Tddd& X, const std::vector<T3Tddd>& V_vertices) { return (CoordinateBounds(V_vertices).InsideQ(X) && windingNumber(X, V_vertices) < 0.75); };
inline bool InsideQ(const Tddd& X, const std::vector<Tddd>& V_vertices) { return !CoordinateBounds(V_vertices).InsideQ(X); };
// bool InsideQ(const CoordinateBounds &bounds, const geometry::Sphere &s) {
//    // cube v.s. sphere
//    // cube < sphere ?
//    auto [X0, X1] = std::get<0>(bounds.bounds);
//    auto [Y0, Y1] = std::get<1>(bounds.bounds);
//    auto [Z0, Z1] = std::get<2>(bounds.bounds);
//    return (Norm(Tddd{X0, Y0, Z0} - s.X) < s.radius &&
//            Norm(Tddd{X1, Y0, Z0} - s.X) < s.radius &&
//            Norm(Tddd{X0, Y1, Z0} - s.X) < s.radius &&
//            Norm(Tddd{X1, Y1, Z0} - s.X) < s.radius &&
//            Norm(Tddd{X0, Y0, Z1} - s.X) < s.radius &&
//            Norm(Tddd{X1, Y0, Z1} - s.X) < s.radius &&
//            Norm(Tddd{X0, Y1, Z1} - s.X) < s.radius &&
//            Norm(Tddd{X1, Y1, Z1} - s.X) < s.radius);
// };
/* ------------------------------------------------------ */

inline Tddd ToX(const Tddd* X) { return *X; };
inline Tddd ToX(const Tddd X) { return X; };
inline T3Tddd ToX(const T3Tddd X) { return X; };

// b% -------------------------------------------------------------------------- */
// b%                                     八分木                                  */
// b% -------------------------------------------------------------------------- */
// 基本的な形状に対してoctreeを生成できるようにする．
template <typename T>
struct octree : public CoordinateBounds {
  /*
  top level -> depth = 0
  ---------
  depth = 1
  ---------
  last level -> depth = depthlimit = 2
  ---------
  */
  octree<T>* const parent;
  const int depth;
  const octree* top; // 最上階のデータ
  const std::vector<T> faces_only_for_top;
  std::vector<T> faces_;
  T8d WNs; // winding numbers
  bool inside;
  std::vector<octree<T>*> children;
  /* ------------------------------------------------------ */
  // objnumの数だけセル中にオブジェクトがあれば，最後までツリーを作成するというもの．もしなければ，minDepthまで作成して終了する．
  // A winding number of 0 means the point is outside the polygon
  octree(const CoordinateBounds& boundsIN, const Tii& depthlimit, const int objnum, const std::vector<T>& FACES)
      : CoordinateBounds(boundsIN), parent(nullptr), depth(0), top(this), faces_only_for_top(FACES), faces_({}), // inside(windingNumber(boundsIN.getCenter(), top->faces_only_for_top) > 0.75),
        WNs(windingNumber((T8Tddd)boundsIN, top->faces_only_for_top)), inside(std::ranges::any_of(WNs, [](const auto& w_num) { return w_num > 0.6; })), children(generateChildrenParallel(depthlimit, objnum, FACES)) {};
  octree(const CoordinateBounds& boundsIN, const Tii& depthlimit, const int objnum, const std::unordered_set<T>& FACES)
      : CoordinateBounds(boundsIN), parent(nullptr), depth(0), top(this), faces_only_for_top(std::vector<T>(FACES.begin(), FACES.end())), faces_({}), // inside(windingNumber(boundsIN.getCenter(), top->faces_only_for_top) > 0.75),
        WNs(windingNumber((T8Tddd)boundsIN, top->faces_only_for_top)), inside(std::ranges::any_of(WNs, [](const auto& w_num) { return w_num > 0.6; })), children(generateChildrenParallel(depthlimit, objnum, faces_only_for_top)) {};
  /* ------------------------------------------------------ */
  octree(const CoordinateBounds& boundsIN, const Tii& depthlimit, const int objnum, const std::vector<T>& FACES, octree<T>* const parentIN)
      : CoordinateBounds(boundsIN), parent(parentIN), depth(parentIN->depth + 1), top(parentIN->top), faces_only_for_top({}), faces_({}),
        //  inside((parentIN && FACES.empty()) ? parentIN->inside : (windingNumber(boundsIN.getCenter(), top->faces_only_for_top) > 0.75)),
        WNs((parentIN && FACES.empty()) ? parentIN->WNs : windingNumber((T8Tddd)boundsIN, top->faces_only_for_top)), inside(std::ranges::any_of(WNs, [](const auto& w_num) { return w_num > 0.6; })), children(generateChildrenParallel(depthlimit, objnum, FACES)) {};
  /* ------------------------------------------------------ */
  std::vector<octree<T>*> generateChildrenParallel(const Tii& depthlimit, const int objnum, const std::vector<T>& FACES) {
    // std::vector<T> faces_;
    faces_.reserve(FACES.size());
    // これまでの取得方法
    //  for (const auto &f : FACES)
    //     if (IntersectQ(this->bounds, ToX(f))) faces_.emplace_back(f);
    //  キューブの球に入るものを取得．SPHで内部外部のチェックがうまくいかないことがあったので，faceを多く取得することにした
    //
    for (const auto& f : FACES) {
      if (Norm(Nearest(this->X, ToX(f)) - this->X) <= this->getScale() / 2.)
        faces_.emplace_back(f);
    }
    //
    if (std::get<0>(depthlimit) >= this->depth || (std::get<1>(depthlimit) >= this->depth && objnum <= faces_.size())) {
      const auto bounds8 = to8Bounds();

#ifdef _OPENMP
      std::vector<octree<T>*> ret(8);
#pragma omp parallel sections
      {
#pragma omp section
        ret[0] = new octree(std::get<0>(bounds8), depthlimit, objnum, faces_, this);
#pragma omp section
        ret[1] = new octree(std::get<1>(bounds8), depthlimit, objnum, faces_, this);
#pragma omp section
        ret[2] = new octree(std::get<2>(bounds8), depthlimit, objnum, faces_, this);
#pragma omp section
        ret[3] = new octree(std::get<3>(bounds8), depthlimit, objnum, faces_, this);
#pragma omp section
        ret[4] = new octree(std::get<4>(bounds8), depthlimit, objnum, faces_, this);
#pragma omp section
        ret[5] = new octree(std::get<5>(bounds8), depthlimit, objnum, faces_, this);
#pragma omp section
        ret[6] = new octree(std::get<6>(bounds8), depthlimit, objnum, faces_, this);
#pragma omp section
        ret[7] = new octree(std::get<7>(bounds8), depthlimit, objnum, faces_, this);
      }
      return ret;
#else
      return {new octree(std::get<0>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<1>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<2>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<3>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<4>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<5>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<6>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<7>(bounds8), depthlimit, objnum, faces_, this)};
#endif
    } else
      return {};
  };
  std::vector<octree<T>*> generateChildren(const Tii& depthlimit, const int objnum, std::vector<T>& FACES) {
    // std::vector<T> faces_;
    faces_.reserve(FACES.size());
    for (const auto& f : FACES)
      if (IntersectQ(this->bounds, ToX(f)))
        faces_.emplace_back(f);
    if (std::get<0>(depthlimit) >= this->depth || (std::get<1>(depthlimit) >= this->depth && objnum <= faces_.size())) {
      //@ min#Objects より多くなければ分割されない
      //@ maxDepth までしか回想は作れればい．最低は０階
      // 例えば，{10,1}の場合，最大で10階まで分割されている．またセルが１個含んでいればそれ以上分割されない．
      //@ また，オブジェクトがゼロなら１つ目の条件から分割されない
      const auto bounds8 = to8Bounds();
      return {new octree(std::get<0>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<1>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<2>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<3>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<4>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<5>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<6>(bounds8), depthlimit, objnum, faces_, this), new octree(std::get<7>(bounds8), depthlimit, objnum, faces_, this)};
    } else
      return {};
  };
  /* ------------------------------------------------------ */
  octree(const CoordinateBounds& boundsIN, const std::vector<T>& FACES, const Tii& depthlimit) : CoordinateBounds(boundsIN), parent(nullptr), depth(0), top(this), faces_only_for_top(FACES), faces_(voronoiConnectivity(FACES)), children(voronoiChildren(depthlimit, faces_)) {};
  octree(const CoordinateBounds& boundsIN, const std::vector<T>& FACES, const Tii& depthlimit, octree<T>* const parentIN) : CoordinateBounds(boundsIN), parent(parentIN), depth(parentIN->depth + 1), top(parentIN->top), faces_(voronoiConnectivity(FACES)), children(voronoiChildren(depthlimit, faces_)) {};
  std::vector<T> voronoiConnectivity(const std::vector<T>& FACES) const {
    Tdd minmax, r = {1E+10, 1E+10};
    std::vector<T> ret;
    for (const auto& f : FACES) {
      minmax = Distance(ToX(f));
      if (std::get<0>(r) >= std::get<0>(minmax))
        std::get<0>(r) = std::get<0>(minmax);
      if (std::get<1>(r) >= std::get<1>(minmax))
        std::get<1>(r) = std::get<1>(minmax);
    }
    ret.reserve(FACES.size());
    for (const auto& f : FACES) {
      minmax = Distance(ToX(f));
      if (Between(std::get<0>(minmax), r) || Between(std::get<1>(minmax), r) || Between(std::get<0>(r), minmax))
        ret.emplace_back(f);
    }
    return ret;
  };
  std::vector<octree<T>*> voronoiChildren(const Tii& depthlimit, const std::vector<T>& faces) {
    if (faces.empty() || std::get<0>(depthlimit) <= this->depth || std::get<1>(depthlimit) >= faces.size()) {
      return {};
    } else {
      std::vector<octree<T>*> ret(8);
      const auto bounds8 = to8Bounds();
#pragma omp parallel sections
      {
#pragma omp section
        ret[0] = new octree<T>(std::get<0>(bounds8), faces, depthlimit, this);
#pragma omp section
        ret[1] = new octree<T>(std::get<1>(bounds8), faces, depthlimit, this);
#pragma omp section
        ret[2] = new octree<T>(std::get<2>(bounds8), faces, depthlimit, this);
#pragma omp section
        ret[3] = new octree<T>(std::get<3>(bounds8), faces, depthlimit, this);
#pragma omp section
        ret[4] = new octree<T>(std::get<4>(bounds8), faces, depthlimit, this);
#pragma omp section
        ret[5] = new octree<T>(std::get<5>(bounds8), faces, depthlimit, this);
#pragma omp section
        ret[6] = new octree<T>(std::get<6>(bounds8), faces, depthlimit, this);
#pragma omp section
        ret[7] = new octree<T>(std::get<7>(bounds8), faces, depthlimit, this);
      }
      return ret;
    }
  };
  void deleteOuside() {
    if (!this->children.empty()) {
      auto tmp = this->children;
      for (auto& c : tmp)
        c->deleteOuside();
    } else if (!this->inside)
      delete this;
  };
  ~octree() {
    auto tmp = this->children;
    for (const auto& c : tmp)
      delete c;
    if (parent)
      parent->children.erase(std::remove(parent->children.begin(), parent->children.end(), this), parent->children.end());
  };
  void getDescendants(auto& accum) const {
    for (const auto& c : this->children) {
      accum.emplace(c);
      c->getDescendants(accum);
    }
  };
  bool hasChildren() const { return !this->children.empty(); };
  /* ------------------------------------------------------ */
  template <typename U>
  bool isAllVertexInsideOf(const U& object) const {
    /*
    このセルの頂点全てが，何かsの中に入っているか？
    all vertices of this cube are in the given object ?
    */
    auto [X0, X1, X2, X3, X4, X5, X6, X7] = this->getVertices();
    return IntersectQ(object, X0) && IntersectQ(object, X1) && IntersectQ(object, X2) && IntersectQ(object, X3) && IntersectQ(object, X4) && IntersectQ(object, X5) && IntersectQ(object, X6) && IntersectQ(object, X7);
  };
  bool isAllVertexInsideOf(const Tddd& s) const { return false; };
  bool isAllVertexInsideOf(const Sphere& s) const {
    auto [X0, X1, X2, X3, X4, X5, X6, X7] = this->getVertices();
    return (Norm(X0 - s.center) <= s.radius && Norm(X1 - s.center) <= s.radius && Norm(X2 - s.center) <= s.radius && Norm(X3 - s.center) <= s.radius && Norm(X4 - s.center) <= s.radius && Norm(X5 - s.center) <= s.radius && Norm(X6 - s.center) <= s.radius && Norm(X7 - s.center) <= s.radius);
  };
  // b* -------------------------------------------------------------------------- */
  // b*                                INTERSECTIONS                               */
  // b* -------------------------------------------------------------------------- */
  std::unordered_set<octree<T>*> getIntersectAsUnorderedSet(const auto& s) const { // 交わる最深階層にあるキューブ
    std::unordered_set<octree<T>*> accum;
    accum.reserve(100000);
    for (const auto& c : this->children)
      c->getIntersect(accum, s);
    return accum;
  };
  void getIntersect(std::unordered_set<octree<T>*>& accum, const auto& s) {
    if (!this->children.empty())
      for (const auto& c : this->children)
        c->getIntersect(accum, s);
    else if (this->inside && IntersectQ(this->bounds, s))
      accum.emplace(this);
  };
  void getAllDeepest(std::unordered_set<octree<T>*>& accum) {
    if (!this->children.empty())
      for (const auto& c : this->children)
        c->getAllDeepest(accum);
    else
      accum.emplace(this);
  };
  void getAllDeepestInside(std::unordered_set<octree<T>*>& accum) {
    if (!this->children.empty())
      for (const auto& c : this->children)
        c->getAllDeepestInside(accum);
    else if (this->inside)
      accum.emplace(this);
  };
  void getAllDeepestInside(std::vector<octree<T>*>& accum) {
    if (!this->children.empty())
      for (const auto& c : this->children)
        c->getAllDeepestInside(accum);
    else if (this->inside)
      accum.emplace_back(this);
  };
  std::unordered_set<octree<T>*> getAllDeepest() {
    std::unordered_set<octree<T>*> ret;
    if (!this->children.empty())
      for (const auto& c : this->children)
        c->getAllDeepest(ret);
    else
      ret.emplace(this);
    return ret;
  };
  /* ------------------------------------------------------ */
  void getDepth(const int& d, std::unordered_set<octree<T>*>& accum) {
    if (this->depth < d && !this->children.empty())
      for (const auto& c : this->children)
        c->getDepth(d, accum);
    else if (this->depth == d)
      accum.emplace(this);
  };
  std::unordered_set<octree<T>*> getDepth(const int& d) {
    std::unordered_set<octree<T>*> ret;
    if (this->depth < d && !this->children.empty())
      for (const auto& c : this->children)
        c->getDepth(d, ret);
    else if (this->depth == d)
      ret.emplace(this);
    return ret;
  };
  /* ------------------------------------------------------ */
  std::unordered_set<octree<T>*> getAllDeepestInside() {
    std::unordered_set<octree<T>*> ret;
    if (!this->children.empty())
      for (const auto& c : this->children)
        c->getAllDeepestInside(ret);
    else if (this->inside) // 全てあるかないかなので，一つでもnullptrなら，最後の階層と考えられる．
      ret.emplace(this);
    return ret;
  };
  /* ------------------------------------------------------ */
  // ここからかんがえてみよう
  void isIntersectInside(bool& found, const auto& s) {
    if (!found && IntersectQ(this->bounds, s)) {
      if (!this->children.empty()) {
        for (const auto& c : this->children)
          c->isIntersectInside(found, s);
      } else if (this->inside)
        found = true;
    }
  };
  bool isIntersectInside(const auto& s) const { // 交わる最深階層にあるキューブ
    bool found = false;
    for (const auto& c : this->children)
      c->isIntersectInside(found, s);
    return found;
  };
  // b@ -------------------------------------------------------------------------- */
  // b@                             STL like functions                             */
  // b@ -------------------------------------------------------------------------- */
  /*
    octreeの場合は，
    最上層から最下層へ探査する．下層にいくにつれて探査すべきセルは増えていくが，もし，途中のセルにオブジェクトが入っていなければ，その時点で下層への探査を中断し別のセルの探査へ移る．

    球体内に点が存在するかどうかを調べる場合．
    * 少なくとも１つのオブジェクトを持つセルが，完全に球体に入る場合，-> trueとなる．
    * また，オブジェクトを持たないセルが，完全に球体に入る場合，そのセルはそれ以上探査する必要はない．
    このように，高速に判定が終了する場合が多い．
    球体表面付近だけ，オブジェクトが位置する場合，探査に時間がかかってしまうと考えられる．
   */
  bool findIntersect(const auto& s, const std::function<bool(const T&)>& func) const {
    if (this->faces_.empty() || !IntersectQ(this->bounds, s.bounds))
      return false;
    // 最下層の場合は，オブジェクト全てとの干渉を調べる
    if (this->children.empty() /*最下層に到達した場合*/) {
      for (const auto& f : this->faces_)
        if (func(f))
          return true;
    }
    // さらに下層を探査する
    for (const auto& c : this->children) {
      if (c->findIntersect(s, func))
        return true;
    }
    return false;
  };
  bool none_of(const Tddd& x, const double radius /*検索半径*/, const std::function<bool(const T&)>& func) const {
    Sphere sphere(x, radius);
    if (IntersectQ(this->bounds, sphere))
      for (const auto& c : this->children)
        if (c->findIntersect(sphere, func))
          return false;
    return true;
  };
  bool none_of(const Tddd& x, const std::function<bool(const T&)>& func) const {
    /*検索半径無限大*/
    Sphere sphere(x, 1E+20);
    if (IntersectQ(this->bounds, sphere))
      for (const auto& c : this->children)
        if (c->findIntersect(sphere, func))
          return false;
    return true;
  };
  // bool any_of(const Tddd &x, const double radius /*検索半径*/, const std::function<bool(const T &)> &func) const {
  //    Sphere sphere(x, radius);
  //    if (IntersectQ(this->bounds, sphere))
  //       for (const auto &c : this->children)
  //          if (c->findIntersect(sphere, func))
  //             return true;
  //    return true;
  // };
  // bool any_of(const Tddd &x, const std::function<bool(const T &)> &func) const {
  //    /*検索半径無限大*/
  //    Sphere sphere(x, 1E+20);
  //    if (IntersectQ(this->bounds, sphere))
  //       for (const auto &c : this->children)
  //          if (c->findIntersect(sphere, func))
  //             return true;
  //    return false;
  // };
  /* ------------------------------------------------------- */
  void getIntersect(std::vector<octree<T>*>& accum, const auto& s) {
    if (IntersectQ(this->bounds, s)) {
      if (!this->children.empty())
        for (const auto& c : this->children)
          c->getIntersect(accum, s);
      else
        accum.emplace_back(this);
    }
  };
  std::vector<octree<T>*> getIntersect(const auto& s) const { // 交わる最深階層にあるキューブ
    std::vector<octree<T>*> accum;
    accum.reserve(100000);
    for (const auto& c : this->children)
      c->getIntersect(accum, s);
    return accum;
  };
  /* ------------------------------------------------------- */
  void getIntersect(octree<T>*& p_cell, const Tddd& X) {
    if (!p_cell && IntersectQ(this->bounds, X)) {
      if (!this->children.empty()) {
        for (const auto& c : this->children)
          c->getIntersect(p_cell, X);
      } else if (this->inside)
        p_cell = this;
    }
  };
  octree<T>* getIntersect(const Tddd& X) const { // 交わる最深階層にあるキューブ
    octree<T>* p_cell = nullptr;
    for (const auto& c : this->children)
      c->getIntersect(p_cell, X);
    return p_cell;
  };

  /* ------------------------------------------------------- */
  void getIntersectInside(std::vector<octree<T>*>& accum, const auto& s) {
    if (IntersectQ(this->bounds, s)) {
      if (!this->children.empty()) {  //! c->c0とせずにc0とだけすることによる問題が多い
                                      // b! 今のエラーは，ここのエラーではない．sphereの干渉チェックの問題
        if (isAllVertexInsideOf(s)) { // 全ての頂点がsの内部にあるので，以降，干渉チェックは行わない．
          for (const auto& c : this->children)
            c->getAllDeepestInside(accum);
        } else {
          for (const auto& c : this->children)
            c->getIntersectInside(accum, s);
        }
      } else if (this->inside)
        accum.emplace_back(this);
    }
  };
  std::vector<octree<T>*> getIntersectInside(const auto& s) const { // 交わる最深階層にあるキューブ
    std::vector<octree<T>*> accum;
    accum.reserve(100000);
    for (const auto& c : this->children)
      c->getIntersectInside(accum, s);
    return accum;
  };
  void apply(const std::function<void(octree<T>*)>& func) const { // 交わる最深階層にあるキューブ
    for (const auto& c : this->children)
      func(c);
  };
  void apply(const std::function<void(octree<T>*)>& func, const auto& s) {
    if (IntersectQ(this->bounds, s)) {
      if (!this->children.empty()) { //! c->c0とせずにc0とだけすることによる問題が多い
        for (const auto& c : this->children)
          c->apply(func, s);
      } else if (this->inside)
        func(this);
    }
  };
  // b! -------------------------------------------------------------------------- */
  // b!                              SETTING NEIGHBORS                             */
  // b! -------------------------------------------------------------------------- */
  std::unordered_set<octree<T>*> neighbors;
  std::array<T, 8> nearest_face;
  std::unordered_set<T> nearest_faces;
  std::array<bool, 8> bools;
  T8d scalers;
  T8Tddd vectors;
  std::unordered_set<T> checked_faces_passed;
  template <typename U>
  U Interpolate(const Tddd& X, const std::array<U, 8>& c) const {
    const auto [x, y, z] = X;
    const auto [c000, c001, c010, c011, c100, c101, c110, c111] = c;
    const auto [x0, x1] = std::get<0>(this->bounds);
    const auto [y0, y1] = std::get<1>(this->bounds);
    const auto [z0, z1] = std::get<2>(this->bounds);
    return (-(((c111 * (x - x0) + c011 * (-x + x1)) * (y - y0) + (c101 * (-x + x0) + c001 * (x - x1)) * (y - y1)) * (z - z0)) + ((c110 * (x - x0) + c010 * (-x + x1)) * (y - y0) + (c100 * (-x + x0) + c000 * (x - x1)) * (y - y1)) * (z - z1)) / ((x0 - x1) * (y0 - y1) * (z0 - z1));
  };
  template <typename U>
  U Integrate(const std::array<U, 8>& c) const {
    const auto [c000, c001, c010, c011, c100, c101, c110, c111] = c;
    const auto [x0, x1] = std::get<0>(this->bounds);
    const auto [y0, y1] = std::get<1>(this->bounds);
    const auto [z0, z1] = std::get<2>(this->bounds);
    return -((c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111) * (x0 - x1) * (y0 - y1) * (z0 - z1)) / 8.;
  };
  /* -------------------------------------------------------------------------- */
  void setNeighbors() {
    /*
    最深層のセルの近傍セルを
    neighbors
    にセットする．
    */
    double distance;
    for (const auto& c : this->getAllDeepest()) {
      c->neighbors.clear();
      auto [X0, X1] = std::get<0>(c->bounds);
      auto [Y0, Y1] = std::get<1>(c->bounds);
      auto [Z0, Z1] = std::get<2>(c->bounds);
      auto d0 = (X1 - X0) * 0.2;
      auto d1 = (Y1 - Y0) * 0.2;
      auto d2 = (Z1 - Z0) * 0.2;
      for (const auto& v : top->getIntersectInside(T3Tdd{{{X0 - d0, X1 + d0}, {Y0 - d1, Y1 + d1}, {Z0 - d2, Z1 + d2}}})) {
        if (v != c)
          c->neighbors.emplace(v);
      }
    }
  };
};

// template <typename = typename std::enable_if<std::is_same<T, T3Tddd>::value>::type>
template <typename T>
inline void setVectorsToTriangle(octree<T>& tree) {
  auto tmp = tree.getAllDeepestInside();
  for (const auto& cell : tmp) {
    cell->checked_faces_passed.clear();
    cell->bools = {false, false, false, false, false, false, false, false};
    for (const auto& f : cell->faces_) {
      // 各頂点にとって最も近い点を抽出
      for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools, [&](const auto& x, auto& s, auto& v, auto& f4v, auto& b) {
        auto XonTriangle = Nearest(x, ToX(f));
        if (Norm(XonTriangle - x) <= s || !b) {
          s = Norm(XonTriangle - x);
          v = XonTriangle - x;
          f4v = f;
          b = true;
        }
      });
      // for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
      //               [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
      //                  auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
      //                  if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
      //                     s = intsp.distance;
      //                     v = intsp.X - x;
      //                     f4v = f;
      //                     b = true;
      //                  }
      //               });
    }
  }
  //
  for (int i = 0; i < 5; ++i)
    for (const auto& cell : tmp) { // cellにとって，最も近い面を，neighborsから探す
      for (const auto& nei : cell->neighbors) {
        std::ranges::for_each(nei->nearest_face, nei->bools, [&](const auto& f, const auto& B) {
          if (B && (cell->checked_faces_passed.emplace(f)).second) {
            // std::cout << "各頂点にとって最も近い点を抽出" << std::endl;
            for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools, [&](const auto& x, auto& s, auto& v, auto& f4v, auto& b) {
              auto XonTriangle = Nearest(x, ToX(f));
              if (Norm(XonTriangle - x) <= s || !b) {
                s = Norm(XonTriangle - x);
                v = XonTriangle - x;
                f4v = f;
                b = true;
              }
            });
            // for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
            //               [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
            //                  auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
            //                  if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
            //                     s = intsp.distance;
            //                     v = intsp.X - x;
            //                     f4v = f;
            //                     b = true;
            //                  }
            //               });
          }
        });
        for (const auto& f : nei->faces_)
          if (cell->checked_faces_passed.emplace(f).second) {
            // std::cout << "各頂点にとって最も近い点を抽出" << std::endl;

            for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools, [&](const auto& x, auto& s, auto& v, auto& f4v, auto& b) {
              auto XonTriangle = Nearest(x, ToX(f));
              if (Norm(XonTriangle - x) <= s || !b) {
                s = Norm(XonTriangle - x);
                v = XonTriangle - x;
                f4v = f;
                b = true;
              }
            });

            // for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
            //               [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
            //                  auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
            //                  if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
            //                     s = intsp.distance;
            //                     v = intsp.X - x;
            //                     f4v = f;
            //                     b = true;
            //                  }
            //               });
          }
      }
    }
};
template <typename T>
inline void setVectorsToTriangle(octree<T>* tree) { setVectorsToTriangle(*tree); };
#endif
