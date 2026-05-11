// BEM_remesh_common.hpp
//
// Shared declarations for the Trial remesh path.

#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include "Network.hpp"
#include "NetworkUtility.hpp"
#include "BEM_inputfile_reader.hpp"

class PVDWriter;

#ifndef BEM_REMESH_WATER_SNAPSHOT_WRITER_DEFINED
#define BEM_REMESH_WATER_SNAPSHOT_WRITER_DEFINED
using RemeshWaterSnapshotWriter =
    std::function<void(Network*, int, const std::string&, PVDWriter&)>;
#endif

// =============================================================================
// Feature edge / vertex classification (shared across all remesh methods)
// =============================================================================
//
// In this code "feature edge" unifies two independent sources of sharpness:
//   1. l->BCInterface     — boundary-condition interface (e.g. waterline, where
//                      Dirichlet meets Neumann). Not a geometric crease.
//   2. l->SharpQ(a)  — geometric crease detected by dihedral-angle > a.
//
// A vertex is classified based on the number of distinct feature-line
// directions incident to it (collinear incidents count as the same line):
//
//     Smooth  = no incident feature edges (and not on a BC BCInterface / BorderQ)
//     Feature = exactly 1 feature line passes through  → motion constrained
//               to that 1D line.
//     Corner  = 2+ distinct feature-line directions, or the vertex itself is
//               flagged as p->BCInterface (BC junction) or p->BorderQ() (mesh
//               boundary) → pinned.
//
// These rules apply to BOTH BC-interface corners (our p->BCInterface / l->BCInterface)
// AND purely geometric corners (a tank edge, a column top rim, etc.) so that
// non-BC geometric features are protected just as rigorously as BC interfaces.

enum class FeatureClass { Smooth, Feature, Corner };

// A line is a "feature" if it carries either a BC interface flag or a sharp
// dihedral. Kept inline and predicate-style so every method (and the smooth
// projector below) uses the same rule.
inline bool is_feature_edge(const networkLine* l, double feature_angle) {
  if (!l) return false;
  if (l->BCInterface) return true;
  if (l->SharpQ(feature_angle)) return true;
  return false;
}

inline bool allows_free_surface_smoothing(const networkPoint* p) {
  if (!p) return false;
  if (p->BorderQ() || p->BCInterface) return false;

  // Smoothing is intentionally limited to pure Dirichlet/free-surface points.
  // Neumann/body-side and BCInterface-neighbor points are projected/repaired by
  // the geometry pipeline; moving them here can pull the waterline/body patch
  // away from the intended geometric constraints before the next repair pass.
  const bool has_dirichlet = p->Dirichlet;
  const bool has_neumann = p->Neumann;
  if (!has_dirichlet || has_neumann) return false;

  for (auto* l : p->getBoundaryLines()) {
    if (!l) continue;
    if (l->BCInterface || l->Neumann)
      return false;
  }
  return true;
}

inline FeatureClass classify_feature_vertex(const networkPoint* p,
                                             double feature_angle,
                                             Tddd* tangent_out = nullptr) {
  if (!p) return FeatureClass::Corner;
  // BC corner / mesh border are always pinned regardless of geometry.
  if (p->BCInterface) return FeatureClass::Corner;
  if (p->BorderQ()) return FeatureClass::Corner;

  // Collect unit directions of incident feature edges.
  std::vector<Tddd> feature_dirs;
  for (auto* l : p->getBoundaryLines()) {
    if (!is_feature_edge(l, feature_angle)) continue;
    auto [pa, pb] = l->getPoints();
    if (!pa || !pb) continue;
    const networkPoint* other = (pa == p) ? pb : pa;
    if (!other) continue;
    const Tddd d = other->X - p->X;
    const double dn = Norm(d);
    if (dn < 1e-14) continue;
    feature_dirs.push_back(d / dn);
  }

  // Merge collinear (both orientations = same line).
  const double collinear_cos_thr = std::cos(10.0 * M_PI / 180.0);
  std::vector<Tddd> distinct;
  for (const auto& d : feature_dirs) {
    bool merged = false;
    for (const auto& e : distinct) {
      if (std::abs(Dot(d, e)) > collinear_cos_thr) { merged = true; break; }
    }
    if (!merged) distinct.push_back(d);
  }

  if (distinct.empty()) return FeatureClass::Smooth;
  if (distinct.size() == 1) {
    if (tangent_out) {
      Tddd tan = {0., 0., 0.};
      for (const auto& d : feature_dirs)
        tan += (Dot(d, distinct[0]) >= 0.0 ? d : -d);
      const double tn = Norm(tan);
      *tangent_out = (tn > 1e-12) ? (tan / tn) : distinct[0];
    }
    return FeatureClass::Feature;
  }
  return FeatureClass::Corner;
}

// Shared collapse-target selector.
// Returns std::nullopt when the collapse must be skipped to preserve a
// feature/corner. Otherwise returns the position the merged vertex should
// end up at.
//
//   Corner ↔ Corner : skip (can't preserve both)
//   Corner ↔ Feature: skip (would break one or the other)
//   Corner ↔ Smooth : Corner position (Corner pinned)
//   Feature ↔ Feature (on different lines): skip
//   Feature ↔ Smooth: Feature position (snap onto feature line)
//   Smooth ↔ Smooth : midpoint
inline std::optional<Tddd> feature_aware_collapse_target(
    const networkPoint* pA, const networkPoint* pB, double feature_angle) {
  if (!pA || !pB) return std::nullopt;
  const auto cA = classify_feature_vertex(pA, feature_angle);
  const auto cB = classify_feature_vertex(pB, feature_angle);
  if (cA != FeatureClass::Smooth && cB != FeatureClass::Smooth)
    return std::nullopt;
  if (cA != FeatureClass::Smooth) return pA->X;
  if (cB != FeatureClass::Smooth) return pB->X;
  return 0.5 * (pA->X + pB->X);
}

// ============================================================================
// Triangle quality primitives
// ============================================================================

// Minimum interior angle of triangle ABC, in radians.
// Returns 0 for degenerate triangles (zero-length edges).
inline double triangle_min_angle(const Tddd& A, const Tddd& B, const Tddd& C) {
  const Tddd AB = B - A, AC = C - A, BC = C - B;
  const double a = Norm(BC);
  const double b = Norm(AC);
  const double c = Norm(AB);
  if (a < 1e-30 || b < 1e-30 || c < 1e-30) return 0.0;
  auto ang_from_sides = [](double opp, double s1, double s2) {
    const double v = (s1 * s1 + s2 * s2 - opp * opp) / (2.0 * s1 * s2);
    return std::acos(std::clamp(v, -1.0, 1.0));
  };
  const double A_ang = ang_from_sides(a, b, c);
  const double B_ang = ang_from_sides(b, a, c);
  const double C_ang = ang_from_sides(c, a, b);
  return std::min({A_ang, B_ang, C_ang});
}

inline double equilateral_vertex_deviation(const Tddd& X0,
                                           const Tddd& X1,
                                           const Tddd& X2) {
  const Tddd base = X2 - X1;
  const double base_len = Norm(base);
  if (!(base_len > 1e-30) || !std::isfinite(base_len))
    return std::numeric_limits<double>::infinity();
  const Tddd Xmid = 0.5 * (X1 + X2);
  const Tddd v_raw = Chop(X0 - Xmid, base);
  const double v_norm = Norm(v_raw);
  if (!(v_norm > 1e-30) || !std::isfinite(v_norm))
    return std::numeric_limits<double>::infinity();
  const double height = 0.5 * std::sqrt(3.0) * base_len;
  if (!(height > 1e-30) || !std::isfinite(height))
    return std::numeric_limits<double>::infinity();
  const Tddd X_ideal = Xmid + height * (v_raw / v_norm);
  return Norm(X_ideal - X0) / height;
}

inline double equilateral_coordinate_deviation(const Tddd& A,
                                               const Tddd& B,
                                               const Tddd& C) {
  const double dA = equilateral_vertex_deviation(A, B, C);
  const double dB = equilateral_vertex_deviation(B, C, A);
  const double dC = equilateral_vertex_deviation(C, A, B);
  if (!std::isfinite(dA) || !std::isfinite(dB) || !std::isfinite(dC))
    return std::numeric_limits<double>::infinity();
  return (dA * dA + dB * dB + dC * dC) / 3.0;
}

// Non-Delaunay test for edge flip.
//   Edge (A,B) with opposite vertices C, D on the two adjacent triangles.
//   Returns true if the flip improves triangle quality geometrically:
//     |C - D| < |A - B|   (opposite points are closer than edge endpoints)
// This is the classical Delaunay test adapted to 3D via chord length.
inline bool is_non_delaunay(const networkLine* l) {
  if (!l) return false;
  auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1]) return false;
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB) return false;
  networkPoint* pC = nullptr;
  networkPoint* pD = nullptr;
  {
    auto [f0p0, f0p1, f0p2] = faces[0]->getPoints();
    for (auto* q : {f0p0, f0p1, f0p2})
      if (q && q != pA && q != pB) { pC = q; break; }
    auto [f1p0, f1p1, f1p2] = faces[1]->getPoints();
    for (auto* q : {f1p0, f1p1, f1p2})
      if (q && q != pA && q != pB) { pD = q; break; }
  }
  if (!pC || !pD) return false;
  const double edge_len = Norm(pB->X - pA->X);
  const double opp_dist = Norm(pD->X - pC->X);
  if (!std::isfinite(edge_len) || !std::isfinite(opp_dist)) return false;
  return opp_dist < edge_len;
}

// ============================================================================
// Normal-flip / coplanar-degeneration guards.
// ============================================================================
//
// "reject any edge collapse or edge flip operation that results in a triangle
//  with a normal that is too different from the original triangle normals"
//
// Both helpers below predict the post-op geometry WITHOUT mutating the mesh.
// They return `true` if the operation is safe (cos(old_normal, new_normal) >=
// cos_threshold for every affected face), `false` if it would collapse the
// local surface patch onto itself or produce a degenerate face.
//
// cos_threshold convention: pass a cosine, e.g.
//   0.866 (= cos 30°)  very strict
//   0.5   (= cos 60 deg) strict-but-practical guard
//   0.0   (= cos 90°)  allow up to right-angle rotation (lenient)
//  -0.3   (= ~cos 107°) currently stored in rs.quality_normal_flip_cos
//                       — reject only when normal is flipping toward inversion
// Anything below -0.5 basically only catches full-inversion (180°).

// Predict post-flip geometry of edge l without mutating. Returns the two
// old and two new face normals. Returns false if the edge is not flippable.
inline bool predict_flip_normals(const networkLine* l,
                                  Tddd& n_old_1, Tddd& n_old_2,
                                  Tddd& n_new_1, Tddd& n_new_2) {
  if (!l) return false;
  const auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1]) return false;
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB) return false;

  // Opposite vertices C on face 0, D on face 1.
  networkPoint* pC = nullptr;
  networkPoint* pD = nullptr;
  auto find_opposite = [](networkFace* f, networkPoint* a, networkPoint* b)
                         -> networkPoint* {
    auto [p0, p1, p2] = f->getPoints();
    for (auto* q : {p0, p1, p2})
      if (q && q != a && q != b) return q;
    return nullptr;
  };
  pC = find_opposite(faces[0], pA, pB);
  pD = find_opposite(faces[1], pA, pB);
  if (!pC || !pD) return false;

  n_old_1 = faces[0]->normal;
  n_old_2 = faces[1]->normal;
  if (!isFinite(n_old_1) || !isFinite(n_old_2)) return false;

  // Flip replaces edge (A,B) with edge (C,D); new triangles are (A,D,C) and
  // (B,C,D). Vertex order chosen so that both new normals line up with the
  // hemisphere of the two old normals (orientation is sign-corrected below).
  const Tddd A = pA->X, B = pB->X, C = pC->X, D = pD->X;
  const Tddd raw1 = Cross(D - A, C - A);  // (A,D,C)
  const Tddd raw2 = Cross(C - B, D - B);  // (B,C,D)
  const double n1n = Norm(raw1), n2n = Norm(raw2);
  if (!(n1n > 1e-14) || !(n2n > 1e-14)) return false;   // degenerate triangle
  n_new_1 = raw1 / n1n;
  n_new_2 = raw2 / n2n;
  // Orient each new normal toward the nearest old normal hemisphere.
  const Tddd n_avg_old = Normalize(n_old_1 + n_old_2);
  if (!isFinite(n_avg_old)) return false;
  if (Dot(n_new_1, n_avg_old) < 0.0) n_new_1 = -n_new_1;
  if (Dot(n_new_2, n_avg_old) < 0.0) n_new_2 = -n_new_2;
  return true;
}

// Returns true if edge flip is geometrically safe under cos_threshold and
// the two new triangles have minimum interior angle >= min_angle_rad.
// Pass min_angle_rad = 0 to disable the angle guard.
inline bool flip_preserves_normals(const networkLine* l, double cos_threshold,
                                    double min_angle_rad = 0.0) {
  Tddd n_old_1, n_old_2, n_new_1, n_new_2;
  if (!predict_flip_normals(l, n_old_1, n_old_2, n_new_1, n_new_2))
    return false;
  const Tddd n_avg_old = Normalize(n_old_1 + n_old_2);
  if (Dot(n_new_1, n_avg_old) < cos_threshold) return false;
  if (Dot(n_new_2, n_avg_old) < cos_threshold) return false;
  if (min_angle_rad > 0.0) {
    // Re-derive new-triangle vertex positions (same as predict_flip_normals).
    auto [pA, pB] = l->getPoints();
    const auto faces = l->getBoundaryFaces();
    auto find_opp = [](networkFace* f, networkPoint* a, networkPoint* b)
                       -> networkPoint* {
      auto [p0, p1, p2] = f->getPoints();
      for (auto* q : {p0, p1, p2}) if (q && q != a && q != b) return q;
      return nullptr;
    };
    networkPoint* pC = find_opp(faces[0], pA, pB);
    networkPoint* pD = find_opp(faces[1], pA, pB);
    if (!pC || !pD) return false;
    // New triangles are (A,D,C) and (B,C,D) — same orientation as above.
    if (triangle_min_angle(pA->X, pD->X, pC->X) < min_angle_rad) return false;
    if (triangle_min_angle(pB->X, pC->X, pD->X) < min_angle_rad) return false;
  }
  return true;
}

// Predicts how the area-weighted boundary normal at each affected vertex
// changes when edge (A,B) is flipped to (C,D), without mutating the mesh.
//
// This catches cases where the two flipped faces individually pass the face
// normal/dihedral guards but the local vertex fan normal rotates too much
// because the incident-face set at A/B/C/D changes.
inline bool flip_preserves_vertex_normals(const networkLine* l,
                                          double cos_threshold,
                                          bool include_opposite_vertices = true) {
  if (!l)
    return false;
  const auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1])
    return false;
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB)
    return false;

  auto find_opp = [](networkFace* f, networkPoint* a, networkPoint* b)
                     -> networkPoint* {
    auto [p0, p1, p2] = f->getPoints();
    for (auto* q : {p0, p1, p2})
      if (q && q != a && q != b)
        return q;
    return nullptr;
  };
  networkPoint* pC = find_opp(faces[0], pA, pB);
  networkPoint* pD = find_opp(faces[1], pA, pB);
  if (!pC || !pD || pA == pB || pA == pC || pA == pD ||
      pB == pC || pB == pD || pC == pD)
    return false;

  Tddd n_old_1, n_old_2, n_new_1, n_new_2;
  if (!predict_flip_normals(l, n_old_1, n_old_2, n_new_1, n_new_2))
    return false;

  const double area_new_1 = TriangleArea(pA->X, pD->X, pC->X);
  const double area_new_2 = TriangleArea(pB->X, pC->X, pD->X);
  if (!(area_new_1 > 1e-20) || !(area_new_2 > 1e-20))
    return false;

  auto face_contains = [](networkFace* f, const networkPoint* p) {
    if (!f || !p)
      return false;
    auto [p0, p1, p2] = f->getPoints();
    return p == p0 || p == p1 || p == p2;
  };

  auto boundary_face_contribution = [](networkFace* f) {
    if (!f || !f->BoundaryQ() || !(f->area > 0.0) || !isFinite(f->normal))
      return Tddd{0.0, 0.0, 0.0};
    return f->area * f->normal;
  };

  auto current_vertex_sum = [&](const networkPoint* p) {
    Tddd sum = {0.0, 0.0, 0.0};
    if (!p)
      return sum;
    for (auto* f : p->getBoundaryFaces())
      sum += boundary_face_contribution(f);
    return sum;
  };

  auto predicted_vertex_sum = [&](const networkPoint* p) {
    Tddd sum = current_vertex_sum(p);
    for (auto* f : faces) {
      if (face_contains(f, p))
        sum -= boundary_face_contribution(f);
    }
    if (p == pA || p == pC || p == pD)
      sum += area_new_1 * n_new_1;
    if (p == pB || p == pC || p == pD)
      sum += area_new_2 * n_new_2;
    return sum;
  };

  auto normal_change_ok = [&](const networkPoint* p) {
    const Tddd before = current_vertex_sum(p);
    const Tddd after = predicted_vertex_sum(p);
    const double nb = Norm(before);
    const double na = Norm(after);
    if (!(nb > 1e-14) || !(na > 1e-14))
      return false;
    const double c = Dot(before / nb, after / na);
    return std::isfinite(c) && c >= cos_threshold;
  };

  if (!normal_change_ok(pA) || !normal_change_ok(pB))
    return false;
  if (include_opposite_vertices &&
      (!normal_change_ok(pC) || !normal_change_ok(pD)))
    return false;
  return true;
}

// ============================================================================
// Flip-specific quality predicates (best-first flip strategy)
// ============================================================================

// Returns the minimum interior angle (rad) across the two new triangles
// produced by flipping edge l. Used as the tiebreak key in best-first flip
// selection. Returns 0 if the edge is not flippable.
inline double flip_new_min_angle(const networkLine* l) {
  if (!l) return 0.0;
  const auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1]) return 0.0;
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB) return 0.0;
  auto find_opp = [](networkFace* f, networkPoint* a, networkPoint* b)
                     -> networkPoint* {
    auto [p0, p1, p2] = f->getPoints();
    for (auto* q : {p0, p1, p2}) if (q && q != a && q != b) return q;
    return nullptr;
  };
  networkPoint* pC = find_opp(faces[0], pA, pB);
  networkPoint* pD = find_opp(faces[1], pA, pB);
  if (!pC || !pD) return 0.0;
  return std::min(triangle_min_angle(pA->X, pD->X, pC->X),
                  triangle_min_angle(pB->X, pC->X, pD->X));
}

// Returns true if the dihedral between the two post-flip faces stays close
// to the dihedral of the two pre-flip faces. Specifically, the drop in
// cos(dihedral) must not exceed `max_cos_drop` (default 0.5 ≈ 60° opening).
// Catches flips that would sharpen a nearly-flat patch into a crease.
inline bool flip_dihedral_change_ok(const networkLine* l,
                                     double max_cos_drop = 0.5) {
  Tddd n_old_1, n_old_2, n_new_1, n_new_2;
  if (!predict_flip_normals(l, n_old_1, n_old_2, n_new_1, n_new_2)) return false;
  const double old_dot = Dot(n_old_1, n_old_2);
  const double new_dot = Dot(n_new_1, n_new_2);
  return (old_dot - new_dot) <= max_cos_drop;
}

// Returns true if the naive midpoint of the new flipped edge (C+D)/2 stays
// within `max_ratio` × old_edge_length of the closer old face plane.
// Catches flips whose new midpoint would fall far from the pre-flip surface
// (important because BEM pseudo-quadratic uses (C+D)/2 as the initial X_mid
// for the new edge; large drift produces surface-representation errors).
inline bool flip_midpoint_drift_ok(const networkLine* l,
                                    double max_ratio = 0.1) {
  if (!l) return false;
  const auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1]) return false;
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB) return false;
  auto find_opp = [](networkFace* f, networkPoint* a, networkPoint* b)
                     -> networkPoint* {
    auto [p0, p1, p2] = f->getPoints();
    for (auto* q : {p0, p1, p2}) if (q && q != a && q != b) return q;
    return nullptr;
  };
  networkPoint* pC = find_opp(faces[0], pA, pB);
  networkPoint* pD = find_opp(faces[1], pA, pB);
  if (!pC || !pD) return false;
  const Tddd new_mid = 0.5 * (pC->X + pD->X);
  auto dist_to_face_plane = [&](networkFace* f) {
    auto [p0, p1, p2] = f->getPoints();
    if (!p0) return std::numeric_limits<double>::infinity();
    const Tddd n = f->normal;
    if (!isFinite(n) || Norm(n) < 1e-14)
      return std::numeric_limits<double>::infinity();
    return std::abs(Dot(new_mid - p0->X, n));
  };
  const double drift = std::min(dist_to_face_plane(faces[0]),
                                dist_to_face_plane(faces[1]));
  const double edge_len = Norm(pB->X - pA->X);
  if (edge_len < 1e-20) return true;
  return (drift / edge_len) <= max_ratio;
}

// Returns true if moving `p` to `new_x` would keep every incident face's
// normal within cos_thr AND minimum interior angle >= min_angle_rad.
// Pass min_angle_rad = 0 to disable the angle guard.
inline bool smooth_preserves_normals(const networkPoint* p, const Tddd& new_x,
                                      double cos_thr,
                                      double min_angle_rad = 0.0) {
  if (!p || !isFinite(new_x)) return false;
  for (auto* f : p->getBoundaryFaces()) {
    if (!f) continue;
    auto [p0, p1, p2] = f->getPoints();
    if (!p0 || !p1 || !p2) return false;
    const Tddd v0 = (p0 == p) ? new_x : p0->X;
    const Tddd v1 = (p1 == p) ? new_x : p1->X;
    const Tddd v2 = (p2 == p) ? new_x : p2->X;
    const Tddd raw = Cross(v1 - v0, v2 - v0);
    const double nn = Norm(raw);
    if (!(nn > 1e-14)) return false;
    const Tddd new_n = raw / nn;
    const Tddd old_n = f->normal;
    if (isFinite(old_n) && Dot(old_n, new_n) < cos_thr) return false;
    if (min_angle_rad > 0.0 &&
        triangle_min_angle(v0, v1, v2) < min_angle_rad) return false;
  }
  return true;
}

// Step-size backtracking fused with a normal check. Starts at α=1 and halves
// until `p->X + α·V` keeps every incident
// face's normal within cos_thr. Returns the accepted α, or 0.0 if no step
// in {1, 0.5, 0.25, ..., 0.5^max_halvings} is safe (= skip this vertex).
inline double find_safe_smooth_step(const networkPoint* p, const Tddd& V,
                                     double cos_thr,
                                     double min_angle_rad = 0.0,
                                     int max_halvings = 4) {
  if (!p || !isFinite(V)) return 0.0;
  double alpha = 1.0;
  for (int i = 0; i <= max_halvings; ++i) {
    const Tddd new_x = p->X + alpha * V;
    if (smooth_preserves_normals(p, new_x, cos_thr, min_angle_rad)) return alpha;
    alpha *= 0.5;
  }
  return 0.0;
}

// Predict post-collapse face normals for every face in 1-ring(pA) ∪ 1-ring(pB)
// except the two faces incident to the collapsed edge (which are deleted).
// Returns true if every surviving face keeps Dot(old_n, new_n) >= threshold.
inline bool collapse_preserves_normals(const networkLine* l,
                                        const Tddd& target,
                                        double cos_threshold,
                                        double min_angle_rad = 0.0) {
  if (!l) return false;
  auto [pA, pB] = l->getPoints();
  if (!pA || !pB) return false;
  if (!isFinite(target)) return false;

  std::unordered_set<const networkFace*> doomed;
  for (auto* f : l->getBoundaryFaces())
    if (f) doomed.insert(f);

  for (const auto* p : {pA, pB}) {
    for (auto* f : p->getBoundaryFaces()) {
      if (!f || doomed.count(f)) continue;
      auto [f_p0, f_p1, f_p2] = f->getPoints();
      if (!f_p0 || !f_p1 || !f_p2) return false;
      const Tddd v0 = (f_p0 == p) ? target : f_p0->X;
      const Tddd v1 = (f_p1 == p) ? target : f_p1->X;
      const Tddd v2 = (f_p2 == p) ? target : f_p2->X;
      const Tddd raw = Cross(v1 - v0, v2 - v0);
      const double nn = Norm(raw);
      if (!(nn > 1e-14)) return false;
      const Tddd new_n = raw / nn;
      const Tddd old_n = f->normal;
      if (isFinite(old_n) && Dot(old_n, new_n) < cos_threshold) return false;
      if (min_angle_rad > 0.0 &&
          triangle_min_angle(v0, v1, v2) < min_angle_rad) return false;
    }
  }
  return true;
}

// =============================================================================
// Shared helper: feature-aware smoothing delta
// =============================================================================
// Input:  p          — vertex being smoothed
//         delta_raw  — proposed motion (target − p->X) in 3D
//         feature_angle — crease threshold (rad)
// Output: V_out      — motion to apply (0 if pinned)
// Return: true if vertex may move; false if pinned (corner / indeterminate)
//
// Classification based on incident SharpQ edges:
//   • 0 crease edges → treat as flat: tangent-plane projection
//   • 1 crease line  → slide along crease tangent (averaged edge direction)
//   • 2+ non-collinear creases → corner, pin
//
// Previously this helper used geom_curvature.PD2 as the crease direction without
// counting incident creases. At a box corner PD2 is still defined but motion along it breaks the
// corner. Counting SharpQ edges first is correct and matches the physical
// meaning of a corner. (local_patch also had the PD2 shortcut but its patches
// are small so the effect is smaller.)
inline bool feature_aware_delta_smooth(const networkPoint* p,
                                       const Tddd& delta_raw,
                                       double feature_angle,
                                       Tddd& V_out) {
  if (!p) return false;

  // Gather incident feature edges (BC BCInterface or SharpQ crease) as candidate
  // feature-line directions. Unified with is_feature_edge so smoothing and
  // classify_feature_vertex see exactly the same edges.
  std::vector<Tddd> crease_dirs;
  for (auto* l : p->getBoundaryLines()) {
    if (!is_feature_edge(l, feature_angle)) continue;
    auto [pa, pb] = l->getPoints();
    if (!pa || !pb) continue;
    const networkPoint* other = (pa == p) ? pb : pa;
    if (!other) continue;
    const Tddd d = other->X - p->X;
    const double dn = Norm(d);
    if (dn < 1e-14) continue;
    crease_dirs.push_back(d / dn);
  }

  // Merge collinear directions (both orientations count as the same line).
  const double collinear_cos_thr = std::cos(10.0 * M_PI / 180.0);
  auto is_collinear = [collinear_cos_thr](const Tddd& a, const Tddd& b) {
    return std::abs(Dot(a, b)) > collinear_cos_thr;
  };
  std::vector<Tddd> distinct;
  for (const auto& d : crease_dirs) {
    bool merged = false;
    for (const auto& e : distinct) { if (is_collinear(d, e)) { merged = true; break; } }
    if (!merged) distinct.push_back(d);
  }

  if (distinct.empty()) {
    // Flat — standard tangent-plane projection using area-weighted normal.
    Tddd n = {0., 0., 0.};
    for (auto* f : p->getBoundaryFaces())
      if (f) n += f->area * f->normal;
    const double nn = Norm(n);
    if (!(nn > 0.) || !std::isfinite(nn)) { V_out = delta_raw; return true; }
    const Tddd nhat = n / nn;
    V_out = delta_raw - Dot(delta_raw, nhat) * nhat;
    return true;
  }
  if (distinct.size() == 1) {
    // Single crease line — slide along the averaged incident edge direction.
    Tddd tangent = {0., 0., 0.};
    for (const auto& d : crease_dirs)
      tangent += (Dot(d, distinct[0]) >= 0.0 ? d : -d);
    const double tn = Norm(tangent);
    if (tn < 1e-12) return false;
    const Tddd dir = tangent / tn;
    V_out = Dot(delta_raw, dir) * dir;
    return true;
  }
  // 2+ non-collinear creases → corner, pin.
  return false;
}

// =============================================================================
// Global Trial pass helpers
// =============================================================================
//
// Consolidated here so local trial code and patch scenario code use the same
// boundary, feature, smoothing, and flip guards.

// BC / feature-edge predicates.
//   split    : allowed on any non-null edge — the new midpoint lands on the
//              feature line if the edge itself is a feature, so splitting
//              preserves the crease geometrically.
//   collapse : forbidden on a BCInterface edge (would destroy the
//              Dirichlet/Neumann interface) or on a geometric crease edge
//              (would delete the feature line). Handled uniformly so the
//              global pass and scenario engine agree.
inline bool bc_allows_split(const networkLine* l) {
  return l != nullptr;
}
inline bool bc_allows_collapse(const networkLine* l, double feature_angle) {
  if (!l) return false;
  if (l->BCInterface) return false;
  if (l->SharpQ(feature_angle)) return false;
  return true;
}

// Edge still present in the network after prior mutations in the same pass.
// Lines are removed from water.Lines by Split/Collapse/Flip; snapshot-and-loop
// pass code must guard each edge with this before dereferencing.
inline bool line_alive(const Network& water, const networkLine* l) {
  return l && (water.Lines.find(const_cast<networkLine*>(l))
               != water.Lines.end());
}

// φ / φ_n end-point averages. Lightweight linear DOF interpolation used by
// split (new midpoint) and collapse (kept vertex). The active Trial path
// normally remaps scalars after topology edits, but these helpers remain useful
// for patch-local fallback operations.
inline double avg_phi(const networkPoint* a, const networkPoint* b) {
  return 0.5 * (std::get<0>(a->phiphin) + std::get<0>(b->phiphin));
}
inline double avg_phin(const networkPoint* a, const networkPoint* b) {
  return 0.5 * (std::get<1>(a->phiphin) + std::get<1>(b->phiphin));
}

// Area-weighted face-barycenter Laplacian centroid of p's 1-ring.
//   g = Σ_{f ∈ 1-ring(p)} A_f · b_f / Σ A_f          (b_f = face barycenter)
// Standard face-area Laplacian form. An earlier variant weighted each 1-ring
// NEIGHBOR VERTEX by that neighbor's own total incident-face area, which is
// not a Laplacian and biased the centroid toward dense regions. Face-based
// weighting produces uniform triangles on flat surfaces.
inline Tddd area_weighted_centroid(const networkPoint* p) {
  double wsum = 0.;
  Tddd gsum = {0., 0., 0.};
  for (auto* f : p->getBoundaryFaces()) {
    if (!f) continue;
    const double A = f->area;
    if (!(A > 0.) || !std::isfinite(A)) continue;
    auto [p0, p1, p2] = f->getPoints();
    if (!p0 || !p1 || !p2) continue;
    const Tddd b = (p0->X + p1->X + p2->X) / 3.0;
    wsum += A;
    gsum += A * b;
  }
  if (!(wsum > 0.) || !std::isfinite(wsum)) return p->X;
  return gsum / wsum;
}

// Free-surface target length derived from the bbox horizontal diagonal.
//   L = sqrt(dx² + dy²) / rs.len_target_divisor
inline double compute_target_len(const Network& water,
                                  const SimulationSettings::RemeshingSettings& rs) {
  const auto& bb = water.getBounds();
  const double dx = std::get<1>(std::get<0>(bb)) - std::get<0>(std::get<0>(bb));
  const double dy = std::get<1>(std::get<1>(bb)) - std::get<0>(std::get<1>(bb));
  const double horiz_diag = std::sqrt(dx * dx + dy * dy);
  const double divisor = (rs.len_target_divisor > 0) ? rs.len_target_divisor : 40.0;
  return horiz_diag / divisor;
}

// =============================================================================
// Shared Trial passes: pass_flip / pass_smooth and supporting types
// =============================================================================
//
// `pass_flip` and `pass_smooth` are defined in BEM_remesh_global_passes.cpp.
//
// SizingField is the per-vertex target edge length carrier for pass_smooth.
// Trial currently uses a uniform field via `make_uniform_field`.

struct SizingField {
  bool is_uniform = true;
  double uniform_L = 0.0;
  std::unordered_map<const networkPoint*, double> L_per_vertex;  // unused if uniform
  double L_min = 0.0;
  double L_max = 0.0;

  double at(const networkPoint* p) const {
    if (is_uniform) return uniform_L;
    if (!p) return 0.0;
    auto it = L_per_vertex.find(p);
    return (it != L_per_vertex.end()) ? it->second : 0.0;
  }
  double edge_target(const networkLine* l) const {
    if (!l) return 0.0;
    if (is_uniform) return uniform_L;
    auto [pA, pB] = l->getPoints();
    const double LA = at(pA);
    const double LB = at(pB);
    if (!(LA > 0.) && !(LB > 0.)) return 0.0;
    if (!(LA > 0.)) return LB;
    if (!(LB > 0.)) return LA;
    return std::min(LA, LB);
  }
};

inline SizingField make_uniform_field(double L) {
  SizingField sf;
  sf.is_uniform = true;
  sf.uniform_L = L;
  return sf;
}

inline bool remesh_tri_circumcenter(const Tddd& A, const Tddd& B, const Tddd& C,
                                    Tddd& CC_out) {
  const Tddd u = B - A;
  const Tddd v = C - A;
  const double uu = Dot(u, u);
  const double vv = Dot(v, v);
  const double uv = Dot(u, v);
  const double det = uu * vv - uv * uv;
  if (!(det > 1e-30) || !std::isfinite(det)) return false;
  const double alpha = vv * (uu - uv) / (2.0 * det);
  const double beta  = uu * (vv - uv) / (2.0 * det);
  CC_out = A + alpha * u + beta * v;
  return isFinite(CC_out);
}

inline Tddd remesh_odt_center_for_face(const networkFace* f, double feature_angle) {
  auto [p0, p1, p2] = f->getPoints();
  const Tddd A = p0->X, B = p1->X, C = p2->X;
  const Tddd bary = (A + B + C) / 3.0;
  for (auto* l : f->getLines()) {
    if (!l) continue;
    if (l->BCInterface || l->SharpQ(feature_angle)) return bary;
  }
  Tddd cc;
  if (!remesh_tri_circumcenter(A, B, C, cc)) return bary;
  const double edge_mean = (Norm(B - A) + Norm(C - B) + Norm(A - C)) / 3.0;
  if (!(edge_mean > 0.0) || !std::isfinite(edge_mean)) return bary;
  if (Norm(cc - bary) > 2.0 * edge_mean) return bary;
  return cc;
}

template <class PositionFn>
inline Tddd remesh_quality_smoothing_vector(const networkPoint* p,
                                            const Tddd& current_pX,
                                            PositionFn position) {
  auto faces = p->BCInterface ? p->getFacesDirichlet() : p->getBoundaryFaces();
  std::vector<double> weights;
  std::vector<Tddd> positions;

  for (const auto& f : faces) {
    auto [p0, p1, p2] = f->getPoints(p);
    Tddd X0 = current_pX;
    Tddd X1 = position(p1);
    Tddd X2 = position(p2);
    const double rr = CircumradiusToInradius(X0, X1, X2);
    if (!std::isfinite(rr))
      continue;
    Tddd Xmid = (X2 + X1) * 0.5;
    const auto e = X2 - X1;
    const auto e_norm = Norm(e);
    if (!(e_norm > 0) || !std::isfinite(e_norm))
      continue;
    auto v_raw = Chop(X0 - Xmid, e);
    const auto v_norm = Norm(v_raw);
    if (!(v_norm > 0) || !std::isfinite(v_norm))
      continue;
    Tddd vertical = v_raw / v_norm;
    double height = e_norm * std::sqrt(3.) * 0.5;
    if (!(height > 0) || !std::isfinite(height))
      continue;
    Tddd X_ideal = height * vertical + Xmid;
    Tddd To_ideal = X_ideal - current_pX;
    const double normalized_discrepancy = Norm(To_ideal) / height;

    // Surface-face part of DistorsionMeasureWeightedSmoothingVector_modified.
    // Tetrahedral terms are intentionally ignored in remesh scenario smoothing.
    const double W = 5.0 * (normalized_discrepancy + (rr - 2.0));
    if (!std::isfinite(W))
      continue;
    weights.push_back(W);
    positions.emplace_back(To_ideal);
  }

  if (weights.empty())
    return {0., 0., 0.};

  double sum = 0.0;
  for (const double w : weights)
    sum += w;
  if (!std::isfinite(sum))
    return {0., 0., 0.};
  if (sum < 1e-12) {
    const double inv_n = 1.0 / static_cast<double>(weights.size());
    for (double& w : weights)
      w *= inv_n;
  } else {
    for (double& w : weights)
      w /= sum;
  }

  Tddd V = {0., 0., 0.};
  for (std::size_t i = 0; i < weights.size(); ++i)
    V += weights[i] * positions[i];
  return V;
}

inline bool remesh_smooth_delta(
    const networkPoint* p,
    SimulationSettings::RemeshingSettings::SmoothMode smooth_mode,
    const SizingField& sf,
    double feature_angle,
    Tddd& V_out) {
  using SmoothMode = SimulationSettings::RemeshingSettings::SmoothMode;
  if (!p) return false;

  if (smooth_mode == SmoothMode::CircumradiusToInradius) {
    auto pos_of = [](const networkPoint* q) -> Tddd { return q->X; };
    const Tddd raw = remesh_quality_smoothing_vector(p, p->X, pos_of);
    return feature_aware_delta_smooth(p, raw, feature_angle, V_out);
  }

  double wsum = 0.0;
  Tddd gsum = {0., 0., 0.};
  for (auto* f : p->getBoundaryFaces()) {
    if (!f) continue;
    const double A = f->area;
    if (!(A > 0.0) || !std::isfinite(A)) continue;
    auto [fp0, fp1, fp2] = f->getPoints();
    if (!fp0 || !fp1 || !fp2) continue;

    Tddd center;
    double w = A;
    switch (smooth_mode) {
      case SmoothMode::AreaLaplacian:
        center = (fp0->X + fp1->X + fp2->X) / 3.0;
        break;
      case SmoothMode::OdtLaplacian: {
        center = (fp0->X + fp1->X + fp2->X) / 3.0;
        const double Lb = (sf.at(fp0) + sf.at(fp1) + sf.at(fp2)) / 3.0;
        w = A * Lb;
        break;
      }
      case SmoothMode::OdtCircumcenter:
        center = remesh_odt_center_for_face(f, feature_angle);
        if (!isFinite(center)) continue;
        break;
      case SmoothMode::CircumradiusToInradius:
        continue;
    }
    if (!std::isfinite(w) || !(w > 0.0)) continue;
    wsum += w;
    gsum += w * center;
  }
  if (!(wsum > 0.0)) return false;
  const Tddd raw = (gsum / wsum) - p->X;
  return feature_aware_delta_smooth(p, raw, feature_angle, V_out);
}

inline double remesh_smooth_step_limit_ratio(
    SimulationSettings::RemeshingSettings::SmoothMode smooth_mode) {
  using SmoothMode = SimulationSettings::RemeshingSettings::SmoothMode;
  // The R/r quality smoother is closer to the original node_relocation stage
  // than to ODT/Laplacian smoothing.  That stage effectively moved about
  // 0.05 * (0.3 * local length) per fixed-point iteration.  Keeping the same
  // scale avoids visible step-to-step oscillation when connectivity is almost
  // unchanged and only the quality target alternates.
  if (smooth_mode == SmoothMode::CircumradiusToInradius)
    return 0.015;
  return 0.01;
}

std::pair<int, int> valenceDeviationScore(const networkLine* l, int s_mean = 6);
std::pair<int, int> sharpSectorValenceDeviationScore(
    const networkLine* l, double feature_angle, int s_mean = 6);

// pass_flip / pass_smooth — shared global Trial polish passes.
// Definitions in BEM_remesh_global_passes.cpp (external linkage).
std::size_t pass_flip(
    Network& water,
    SimulationSettings::RemeshingSettings::FlipMode flip_mode,
    double feature_angle, double cos_thr, double min_angle_rad,
    bool use_sharp_sector_valence = false);

std::size_t pass_smooth(
    Network& water,
    SimulationSettings::RemeshingSettings::SmoothMode smooth_mode,
    const SizingField& sf, int sub_iterations,
    double feature_angle, double cos_thr, double min_angle_rad,
    bool skip_sharp_points = false);

// =============================================================================
// Method implementation
// =============================================================================

namespace BEMMeshPipeline {
struct RemeshTrialGeometryInputs;
}

// Local-patch trial-and-error with multi-objective patchQuality score.
// File: BEM_remesh_local_patch.cpp
void remesh_local_patch_trial_multi_objective(
    Network& water, int time_step,
    const SimulationSettings::RemeshingSettings& rs,
    bool skip_post_remesh_quality_rejects,
    const std::string& patch_output_directory,
    double simulation_time,
    PVDWriter* candidate_patches_pvd,
    PVDWriter* remeshed_patches_pvd,
    PVDWriter* trigger_edges_pvd,
    PVDWriter* edges_pvd,
    const BEMMeshPipeline::RemeshTrialGeometryInputs* trial_geometry_inputs = nullptr,
    const std::string& phase_debug_tag = "",
    PVDWriter* remeshed_water_unrepaired_pvd = nullptr,
    int step_retry = 0,
    RemeshWaterSnapshotWriter remeshed_water_snapshot_writer = RemeshWaterSnapshotWriter{});
