#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "BEM_setBoundaryTypes.hpp"
#include "Network.hpp"
#include "basic.hpp"

namespace BEMPreBVP {

struct Options {
  bool print_ok_summary = false;
  bool print_waterline_midpoint_quality = false;
  int max_examples = 6;
  double scalar_overflow_abs = 1e10;
};

struct WaterlineQualityWitness {
  int face_index = -1;
  int subtri_index = -1;
  double value = std::numeric_limits<double>::quiet_NaN();
};

struct WaterlineMidpointStats {
  int waterline_subtri_count = 0;
  double subtri_min_angle_deg = std::numeric_limits<double>::infinity();
  double subtri_max_aspect = 0.0;
  double subtri_min_area = std::numeric_limits<double>::infinity();
  WaterlineQualityWitness min_angle_witness;
  WaterlineQualityWitness max_aspect_witness;
};

struct WaterlineClampedEntity {
  int line_index = -1;
  int endpoint0_index = -1;
  int endpoint1_index = -1;
  int reference_face_index = -1;
  int body_face_index = -1;
  double free_gap = 0.0;
  double body_gap = 0.0;
  double move_ratio = 0.0;
  double local_min_angle_before = std::numeric_limits<double>::quiet_NaN();
  double local_min_angle_after = std::numeric_limits<double>::quiet_NaN();
  double local_max_aspect_before = std::numeric_limits<double>::quiet_NaN();
  double local_max_aspect_after = std::numeric_limits<double>::quiet_NaN();
  Tddd midpoint = {0., 0., 0.};
};

struct WaterlineRefreshSummary {
  bool valid = false;
  int lines_total = 0;
  int ok = 0;
  int ok_via_fallback = 0;
  int no_convergence = 0;
  int quality_damaged_after_repair = 0;
  int quality_rejected_before_repair = 0;
  int raw_contact_rejected_before_repair = 0;
  int skipped = 0;
  int failed_no_dirichlet = 0;
  int failed_no_body = 0;
  int failed_body_gap = 0;
  int failed_no_move = 0;
  int failed_invalid = 0;
  int clamped_count = 0;
  double max_free_gap = 0.0;
  double max_body_gap = 0.0;
  double max_move_ratio = 0.0;
  std::vector<WaterlineClampedEntity> clamped_entities;
};

inline WaterlineRefreshSummary& latest_waterline_refresh_summary() {
  static WaterlineRefreshSummary summary;
  return summary;
}

inline void set_latest_waterline_refresh_summary(const WaterlineRefreshSummary& summary) {
  latest_waterline_refresh_summary() = summary;
}

inline void clear_latest_waterline_refresh_summary() {
  latest_waterline_refresh_summary() = {};
}

struct Stats {
  int water_count = 0;
  int boundary_points = 0;
  int boundary_lines = 0;
  int boundary_faces = 0;
  int unclassified_points = 0;
  int unclassified_lines = 0;
  int bcinterface_points = 0;
  int bcinterface_lines = 0;
  int bcinterface_points_no_contact = 0;
  int bcinterface_lines_no_contact = 0;
  int active_dofs = 0;
  int nonfinite_dofs = 0;
  int overflow_dofs = 0;
  int dirichlet_entities = 0;
  int dirichlet_phi_missing = 0;
  int bcinterface_dirichlet_face_dofs = 0;
  int bcinterface_dirichlet_face_phi_missing = 0;
  int waterline_faces = 0;
  double waterline_edge_min = std::numeric_limits<double>::infinity();
  double waterline_edge_max = 0.0;
  double waterline_face_area_min = std::numeric_limits<double>::infinity();
  double waterline_face_min_angle_deg = std::numeric_limits<double>::infinity();
  double waterline_face_aspect_max = 0.0;
  WaterlineMidpointStats waterline_midpoint_quality;
  WaterlineRefreshSummary waterline_refresh;
  std::vector<std::string> examples;

  bool has_warning() const {
    return unclassified_points > 0 || unclassified_lines > 0 ||
           bcinterface_points_no_contact > 0 || bcinterface_lines_no_contact > 0 ||
           dirichlet_phi_missing > 0 || bcinterface_dirichlet_face_phi_missing > 0 ||
           nonfinite_dofs > 0 || overflow_dofs > 0;
  }

  bool has_numeric_failure() const {
    return nonfinite_dofs > 0 || overflow_dofs > 0;
  }
};

inline double safe_angle_deg(const Tddd& a, const Tddd& b) {
  const double na = Norm(a);
  const double nb = Norm(b);
  if (!(na > 0.0) || !(nb > 0.0) || !std::isfinite(na) || !std::isfinite(nb))
    return std::numeric_limits<double>::infinity();
  const double c = std::clamp(Dot(a, b) / (na * nb), -1.0, 1.0);
  return std::acos(c) * 180.0 / M_PI;
}

inline void accumulate_triangle_quality(WaterlineMidpointStats& out,
                                        const Tddd& a, const Tddd& b, const Tddd& c,
                                        int face_index = -1,
                                        int subtri_index = -1) {
  const double la = Norm(b - a);
  const double lb = Norm(c - b);
  const double lc = Norm(a - c);
  const double area = 0.5 * Norm(Cross(b - a, c - a));
  if (!(std::isfinite(la) && std::isfinite(lb) && std::isfinite(lc) &&
        std::isfinite(area) && area > 0.0))
    return;

  ++out.waterline_subtri_count;
  out.subtri_min_area = std::min(out.subtri_min_area, area);
  const double max_edge = std::max({la, lb, lc});
  const double altitude = 2.0 * area / max_edge;
  if (std::isfinite(altitude) && altitude > 0.0) {
    const double aspect = max_edge / altitude;
    if (std::isfinite(aspect) && aspect > out.subtri_max_aspect) {
      out.subtri_max_aspect = aspect;
      out.max_aspect_witness = {face_index, subtri_index, aspect};
    }
  }
  const double angle0 = safe_angle_deg(b - a, c - a);
  const double angle1 = safe_angle_deg(a - b, c - b);
  const double angle2 = safe_angle_deg(a - c, b - c);
  const double min_angle = std::min({angle0, angle1, angle2});
  if (std::isfinite(min_angle) && min_angle < out.subtri_min_angle_deg) {
    out.subtri_min_angle_deg = min_angle;
    out.min_angle_witness = {face_index, subtri_index, min_angle};
  }
}

inline Tddd midpoint_for_quality(const networkLine* l,
                                 const networkLine* override_line = nullptr,
                                 std::optional<Tddd> override_midpoint = std::nullopt) {
  if (!l)
    return {0., 0., 0.};
  if (override_line && l == override_line && override_midpoint)
    return *override_midpoint;
  if (isFinite(l->X_mid))
    return l->X_mid;
  auto [pA, pB] = l->getPoints();
  if (pA && pB)
    return 0.5 * (pA->X + pB->X);
  return {0., 0., 0.};
}

inline Tddd point_for_quality(const networkPoint* p,
                              const networkPoint* override_point = nullptr,
                              std::optional<Tddd> override_x = std::nullopt) {
  if (!p)
    return {0., 0., 0.};
  if (override_point && p == override_point && override_x)
    return *override_x;
  return p->X;
}

inline bool is_waterline_face(const networkFace* f) {
  if (!f)
    return false;
  auto [p0, p1, p2] = f->getPoints();
  auto [l0, l1, l2] = f->getLines();
  return (p0 && p0->BCInterface) || (p1 && p1->BCInterface) || (p2 && p2->BCInterface) ||
         (l0 && l0->BCInterface) || (l1 && l1->BCInterface) || (l2 && l2->BCInterface);
}

inline WaterlineMidpointStats measure_waterline_midpoint_quality(
    const std::vector<networkFace*>& faces,
    const networkLine* override_line = nullptr,
    std::optional<Tddd> override_midpoint = std::nullopt,
    const networkPoint* override_point = nullptr,
    std::optional<Tddd> override_point_x = std::nullopt) {
  WaterlineMidpointStats stats;
  for (const auto* f : faces) {
    if (!is_waterline_face(f))
      continue;
    auto [p0, p1, p2] = f->getPoints();
    auto [l0, l1, l2] = f->getLines();
    if (!p0 || !p1 || !p2 || !l0 || !l1 || !l2)
      continue;
    const Tddd x0 = point_for_quality(p0, override_point, override_point_x);
    const Tddd x1 = point_for_quality(p1, override_point, override_point_x);
    const Tddd x2 = point_for_quality(p2, override_point, override_point_x);
    const Tddd m0 = midpoint_for_quality(l0, override_line, override_midpoint);
    const Tddd m1 = midpoint_for_quality(l1, override_line, override_midpoint);
    const Tddd m2 = midpoint_for_quality(l2, override_line, override_midpoint);
    accumulate_triangle_quality(stats, x0, m2, m1, f->index, 0);
    accumulate_triangle_quality(stats, x1, m0, m2, f->index, 1);
    accumulate_triangle_quality(stats, x2, m1, m0, f->index, 2);
    accumulate_triangle_quality(stats, m0, m1, m2, f->index, 3);
  }
  return stats;
}

struct WaterlineQualityBaseline {
  std::vector<networkFace*> all_faces;
  WaterlineMidpointStats full_stats;
};

inline void add_unique_waterline_face(std::vector<networkFace*>& faces, networkFace* f) {
  if (!is_waterline_face(f))
    return;
  if (std::find(faces.begin(), faces.end(), f) == faces.end())
    faces.emplace_back(f);
}

inline std::vector<networkFace*> collect_waterline_faces(const std::vector<networkFace*>& faces) {
  std::vector<networkFace*> waterline_faces;
  waterline_faces.reserve(faces.size());
  for (auto* f : faces)
    add_unique_waterline_face(waterline_faces, f);
  return waterline_faces;
}

inline std::vector<networkFace*> affected_waterline_faces(const networkPoint* p) {
  std::vector<networkFace*> faces;
  if (!p)
    return faces;
  for (auto* f : p->getBoundaryFaces())
    add_unique_waterline_face(faces, f);
  return faces;
}

inline std::vector<networkFace*> affected_waterline_faces(const networkLine* l) {
  std::vector<networkFace*> faces;
  if (!l)
    return faces;
  for (auto* f : l->getBoundaryFaces())
    add_unique_waterline_face(faces, f);
  return faces;
}

inline WaterlineQualityBaseline make_waterline_quality_baseline(const std::vector<networkFace*>& faces) {
  WaterlineQualityBaseline baseline;
  baseline.all_faces = faces;
  baseline.full_stats = measure_waterline_midpoint_quality(faces);
  return baseline;
}

inline bool waterline_witness_face_affected(const WaterlineMidpointStats& stats,
                                            const std::vector<networkFace*>& affected_faces) {
  auto contains_face_index = [&](int face_index) {
    if (face_index < 0)
      return false;
    for (const auto* f : affected_faces)
      if (f && f->index == face_index)
        return true;
    return false;
  };
  return contains_face_index(stats.min_angle_witness.face_index) ||
         contains_face_index(stats.max_aspect_witness.face_index);
}

inline WaterlineMidpointStats combine_waterline_quality_baseline(
    const WaterlineMidpointStats& baseline,
    const WaterlineMidpointStats& affected_post) {
  auto combined = baseline;
  if (affected_post.waterline_subtri_count <= 0)
    return combined;
  if (affected_post.subtri_min_angle_deg < combined.subtri_min_angle_deg) {
    combined.subtri_min_angle_deg = affected_post.subtri_min_angle_deg;
    combined.min_angle_witness = affected_post.min_angle_witness;
  }
  if (affected_post.subtri_max_aspect > combined.subtri_max_aspect) {
    combined.subtri_max_aspect = affected_post.subtri_max_aspect;
    combined.max_aspect_witness = affected_post.max_aspect_witness;
  }
  combined.subtri_min_area = std::min(combined.subtri_min_area, affected_post.subtri_min_area);
  return combined;
}

inline WaterlineMidpointStats measure_waterline_midpoint_quality_delta(
    const WaterlineQualityBaseline& baseline,
    const std::vector<networkFace*>& affected_faces,
    const networkLine* override_line = nullptr,
    std::optional<Tddd> override_midpoint = std::nullopt,
    const networkPoint* override_point = nullptr,
    std::optional<Tddd> override_point_x = std::nullopt) {
  if (baseline.full_stats.waterline_subtri_count <= 0)
    return measure_waterline_midpoint_quality(baseline.all_faces, override_line, override_midpoint,
                                              override_point, override_point_x);
  if (affected_faces.empty())
    return baseline.full_stats;
  if (waterline_witness_face_affected(baseline.full_stats, affected_faces))
    return measure_waterline_midpoint_quality(baseline.all_faces, override_line, override_midpoint,
                                              override_point, override_point_x);
  const auto affected_post = measure_waterline_midpoint_quality(affected_faces, override_line, override_midpoint,
                                                                override_point, override_point_x);
  return combine_waterline_quality_baseline(baseline.full_stats, affected_post);
}

inline void merge_waterline_midpoint_quality(WaterlineMidpointStats& into,
                                             const WaterlineMidpointStats& from) {
  into.waterline_subtri_count += from.waterline_subtri_count;
  if (from.subtri_min_angle_deg < into.subtri_min_angle_deg) {
    into.subtri_min_angle_deg = from.subtri_min_angle_deg;
    into.min_angle_witness = from.min_angle_witness;
  }
  if (from.subtri_max_aspect > into.subtri_max_aspect) {
    into.subtri_max_aspect = from.subtri_max_aspect;
    into.max_aspect_witness = from.max_aspect_witness;
  }
  into.subtri_min_area = std::min(into.subtri_min_area, from.subtri_min_area);
}

inline std::string waterline_midpoint_witness_summary(const WaterlineMidpointStats& q) {
  std::ostringstream oss;
  oss << " min_angle_face=" << q.min_angle_witness.face_index
      << " min_angle_subtri=" << q.min_angle_witness.subtri_index
      << " min_angle_value=" << q.min_angle_witness.value
      << " max_aspect_face=" << q.max_aspect_witness.face_index
      << " max_aspect_subtri=" << q.max_aspect_witness.subtri_index
      << " max_aspect_value=" << q.max_aspect_witness.value;
  return oss.str();
}

inline int contact_count(const auto* entity) {
  return static_cast<int>(entity->getContactFaces().size());
}

inline void add_example(Stats& stats, const Options& opt, const std::string& s) {
  if (static_cast<int>(stats.examples.size()) < opt.max_examples)
    stats.examples.push_back(s);
}

inline void check_dof_scalars(Stats& stats, const Options& opt,
                              const char* kind, int entity_index, int face_index,
                              const auto& dof) {
  auto check_one = [&](const char* name, double v, bool overflow_is_failure) {
    if (!std::isfinite(v)) {
      ++stats.nonfinite_dofs;
      std::ostringstream oss;
      oss << kind << "=" << entity_index << " face=" << face_index
          << " active_dof=" << dof.index << " nonfinite " << name << "=" << v;
      add_example(stats, opt, oss.str());
    } else if (overflow_is_failure && std::abs(v) > opt.scalar_overflow_abs) {
      ++stats.overflow_dofs;
      std::ostringstream oss;
      oss << kind << "=" << entity_index << " face=" << face_index
          << " active_dof=" << dof.index << " overflow " << name << "=" << v;
      add_example(stats, opt, oss.str());
    }
  };

  if (dof.index < 0)
    return;
  ++stats.active_dofs;
  check_one("phi", dof.phi, true);
  check_one("phin", dof.phin, true);
  // phi_t/phin_t can legitimately carry a large unset sentinel before the
  // acceleration solve. Non-finite values are still a hard data error.
  check_one("phi_t", dof.phi_t, false);
  check_one("phin_t", dof.phin_t, false);
}

inline void check_dirichlet_phi_defined(Stats& stats, const Options& opt,
                                        const char* kind, int entity_index,
                                        const auto* entity) {
  if (!hasAnyDirichletBoundaryState(entity))
    return;
  ++stats.dirichlet_entities;
  const double phi = std::get<0>(entity->phiphin);
  if (!std::isfinite(phi)) {
    ++stats.dirichlet_phi_missing;
    std::ostringstream oss;
    oss << kind << "=" << entity_index << " Dirichlet phi is not finite: " << phi;
    add_example(stats, opt, oss.str());
  }
}

inline void check_bcinterface_per_face_phi_defined(Stats& stats, const Options& opt,
                                              const char* kind, int entity_index,
                                              const auto* entity) {
  if (!entity || !entity->BCInterface)
    return;
  for (const auto& [f, dof] : entity->dofs) {
    if (!isDirichletBieDofKey(entity, f))
      continue;
    ++stats.bcinterface_dirichlet_face_dofs;
    if (!std::isfinite(dof.phi)) {
      ++stats.bcinterface_dirichlet_face_phi_missing;
      std::ostringstream oss;
      oss << kind << "=" << entity_index << " face=" << (f ? f->index : -1)
          << " BCInterface Dirichlet per-face phi is not finite: " << dof.phi;
      add_example(stats, opt, oss.str());
    }
  }
}

inline void accumulate_face_quality(Stats& stats, const networkFace* f) {
  if (!f)
    return;
  auto [p0, p1, p2] = f->getPoints();
  auto [l0, l1, l2] = f->getLines();
  const bool waterline_face = is_waterline_face(f);
  if (!waterline_face || !p0 || !p1 || !p2)
    return;

  ++stats.waterline_faces;
  const Tddd a = p1->X - p0->X;
  const Tddd b = p2->X - p1->X;
  const Tddd c = p0->X - p2->X;
  const double la = Norm(a);
  const double lb = Norm(b);
  const double lc = Norm(c);
  if (std::isfinite(la) && la > 0.0) {
    stats.waterline_edge_min = std::min(stats.waterline_edge_min, la);
    stats.waterline_edge_max = std::max(stats.waterline_edge_max, la);
  }
  if (std::isfinite(lb) && lb > 0.0) {
    stats.waterline_edge_min = std::min(stats.waterline_edge_min, lb);
    stats.waterline_edge_max = std::max(stats.waterline_edge_max, lb);
  }
  if (std::isfinite(lc) && lc > 0.0) {
    stats.waterline_edge_min = std::min(stats.waterline_edge_min, lc);
    stats.waterline_edge_max = std::max(stats.waterline_edge_max, lc);
  }

  const double area = 0.5 * Norm(Cross(p1->X - p0->X, p2->X - p0->X));
  if (std::isfinite(area) && area > 0.0) {
    stats.waterline_face_area_min = std::min(stats.waterline_face_area_min, area);
    const double max_edge = std::max({la, lb, lc});
    const double altitude = 2.0 * area / max_edge;
    if (std::isfinite(altitude) && altitude > 0.0)
      stats.waterline_face_aspect_max = std::max(stats.waterline_face_aspect_max, max_edge / altitude);
  }

  const double angle0 = safe_angle_deg(p1->X - p0->X, p2->X - p0->X);
  const double angle1 = safe_angle_deg(p0->X - p1->X, p2->X - p1->X);
  const double angle2 = safe_angle_deg(p0->X - p2->X, p1->X - p2->X);
  stats.waterline_face_min_angle_deg = std::min({stats.waterline_face_min_angle_deg, angle0, angle1, angle2});
  merge_waterline_midpoint_quality(stats.waterline_midpoint_quality,
                                   measure_waterline_midpoint_quality({const_cast<networkFace*>(f)}));
}

inline Stats inspect(const std::vector<Network*>& fluid_nets, int time_step, int rk_step,
                     const Options& opt = {}, std::ostream& os = std::cerr) {
  Stats stats;
  stats.water_count = static_cast<int>(fluid_nets.size());
  stats.waterline_refresh = latest_waterline_refresh_summary();

  for (const auto* water : fluid_nets) {
    if (!water)
      continue;
    for (auto* p : water->getBoundaryPoints()) {
      if (!p)
        continue;
      ++stats.boundary_points;
      if (!p->Dirichlet && !p->Neumann && !p->BCInterface) {
        ++stats.unclassified_points;
        std::ostringstream oss;
        oss << "point=" << p->index << " unclassified pos=" << p->X;
        add_example(stats, opt, oss.str());
      }
      if (p->BCInterface) {
        ++stats.bcinterface_points;
        if (contact_count(p) == 0) {
          ++stats.bcinterface_points_no_contact;
          std::ostringstream oss;
          oss << "point=" << p->index << " BCInterface contact_faces=0 pos=" << p->X
              << " adjacent_faces=" << p->getBoundaryFaces().size();
          add_example(stats, opt, oss.str());
        }
      }
      check_dirichlet_phi_defined(stats, opt, "point", p->index, p);
      check_bcinterface_per_face_phi_defined(stats, opt, "point", p->index, p);
      for (const auto& [f, dof] : p->dofs)
        check_dof_scalars(stats, opt, "point", p->index, f ? f->index : -1, dof);
    }

    for (auto* l : water->getBoundaryLines()) {
      if (!l)
        continue;
      ++stats.boundary_lines;
      if (!l->Dirichlet && !l->Neumann && !l->BCInterface) {
        ++stats.unclassified_lines;
        auto [pA, pB] = l->getPoints();
        std::ostringstream oss;
        oss << "line endpoints={" << (pA ? pA->index : -1) << "," << (pB ? pB->index : -1)
            << "} unclassified";
        add_example(stats, opt, oss.str());
      }
      if (l->BCInterface) {
        ++stats.bcinterface_lines;
        if (contact_count(l) == 0) {
          ++stats.bcinterface_lines_no_contact;
          auto [pA, pB] = l->getPoints();
          std::ostringstream oss;
          oss << "line endpoints={" << (pA ? pA->index : -1) << "," << (pB ? pB->index : -1)
              << "} BCInterface contact_faces=0";
          add_example(stats, opt, oss.str());
        }
      }
      check_dirichlet_phi_defined(stats, opt, "line", l->midpoint_index, l);
      check_bcinterface_per_face_phi_defined(stats, opt, "line", l->midpoint_index, l);
      for (const auto& [f, dof] : l->dofs)
        check_dof_scalars(stats, opt, "line", l->midpoint_index, f ? f->index : -1, dof);
    }

    for (const auto* f : water->getBoundaryFaces()) {
      ++stats.boundary_faces;
      accumulate_face_quality(stats, f);
    }
  }

  if (stats.has_warning() || opt.print_ok_summary) {
    const double wl_edge_ratio =
        (std::isfinite(stats.waterline_edge_min) && stats.waterline_edge_min > 0.0)
            ? stats.waterline_edge_max / stats.waterline_edge_min
            : 0.0;
    os << (stats.has_warning() ? Yellow : Green)
       << "[pre_bvp_guard] " << (stats.has_warning() ? "warning" : "ok")
       << " time_step=" << time_step
       << " rk=" << rk_step
       << " active_dofs=" << stats.active_dofs
       << " unclassified={p:" << stats.unclassified_points << ",l:" << stats.unclassified_lines << "}"
       << " BCInterface_no_contact={p:" << stats.bcinterface_points_no_contact
       << "/" << stats.bcinterface_points
       << ",l:" << stats.bcinterface_lines_no_contact
       << "/" << stats.bcinterface_lines << "}"
       << " numeric={nonfinite:" << stats.nonfinite_dofs
       << ",overflow:" << stats.overflow_dofs << "}"
       << " phi_undefined={dirichlet:" << stats.dirichlet_phi_missing
       << "/" << stats.dirichlet_entities
       << ",bcinterface_face_dof:" << stats.bcinterface_dirichlet_face_phi_missing
       << "/" << stats.bcinterface_dirichlet_face_dofs << "}"
       << " waterline_quality={faces:" << stats.waterline_faces
       << ",edge_ratio:" << wl_edge_ratio
       << ",min_area:" << (std::isfinite(stats.waterline_face_area_min) ? stats.waterline_face_area_min : 0.0)
       << ",min_angle_deg:" << (std::isfinite(stats.waterline_face_min_angle_deg) ? stats.waterline_face_min_angle_deg : 0.0)
       << ",max_aspect:" << stats.waterline_face_aspect_max
       << "}";
    if (opt.print_waterline_midpoint_quality) {
      os << " waterline_midpoint_quality={subtri_count:" << stats.waterline_midpoint_quality.waterline_subtri_count
         << ",min_angle_deg:" << (std::isfinite(stats.waterline_midpoint_quality.subtri_min_angle_deg) ? stats.waterline_midpoint_quality.subtri_min_angle_deg : 0.0)
         << ",max_aspect:" << stats.waterline_midpoint_quality.subtri_max_aspect
         << ",min_area:" << (std::isfinite(stats.waterline_midpoint_quality.subtri_min_area) ? stats.waterline_midpoint_quality.subtri_min_area : 0.0)
         << "}";
    }
    if (opt.print_waterline_midpoint_quality && stats.waterline_refresh.valid) {
      os << " waterline_refresh={total:" << stats.waterline_refresh.lines_total
         << ",ok:" << stats.waterline_refresh.ok
         << ",fallback:" << stats.waterline_refresh.ok_via_fallback
         << ",no_convergence:" << stats.waterline_refresh.no_convergence
         << ",quality_damaged_after_repair:" << stats.waterline_refresh.quality_damaged_after_repair
         << ",quality_rejected_before_repair:" << stats.waterline_refresh.quality_rejected_before_repair
         << ",raw_contact_rejected_before_repair:" << stats.waterline_refresh.raw_contact_rejected_before_repair
         << ",skipped:" << stats.waterline_refresh.skipped
         << ",failed:{no_dir:" << stats.waterline_refresh.failed_no_dirichlet
         << ",no_body:" << stats.waterline_refresh.failed_no_body
         << ",body_gap:" << stats.waterline_refresh.failed_body_gap
         << ",no_move:" << stats.waterline_refresh.failed_no_move
         << ",invalid:" << stats.waterline_refresh.failed_invalid
         << "},clamped:" << stats.waterline_refresh.clamped_count
         << ",max_free_gap:" << stats.waterline_refresh.max_free_gap
         << ",max_body_gap:" << stats.waterline_refresh.max_body_gap
         << ",max_move_ratio:" << stats.waterline_refresh.max_move_ratio
         << "}";
    }
    os << colorReset << std::endl;
    for (const auto& ex : stats.examples)
      os << Yellow << "  [pre_bvp_guard] " << ex << colorReset << std::endl;
  }

  return stats;
}

} // namespace BEMPreBVP
