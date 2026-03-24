#pragma once

// Structure penetration detection extracted from main_time_domain.cpp.
// Compiled as a separate translation unit (BEM_penetration.cpp).

#include <algorithm>
#include <iostream>
#include <ranges>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Network.hpp"
#include "BEM_step_failure.hpp"

namespace BEM_Penetration {

struct PenetrationReport {
  std::size_t point_count = 0;
  std::size_t face_count = 0;
  double max_distance = 0.0;
  double worst_tolerance = 0.0;
  double worst_local_length = 0.0;
  std::string worst_body;
};

double penetrationTolerance(const double h_local) {
  constexpr double c_rel = 0.25;
  constexpr double c_abs_eps = 1e-6;
  return (h_local > 0.0) ? c_rel * h_local + c_abs_eps : c_abs_eps;
}

std::vector<Network*> penetrationCandidateSolids(const std::vector<Network*>& fluids,
                                                 const std::vector<Network*>& solids) {
  static std::unordered_map<const Network*, bool> enclosure_cache;
  static std::unordered_set<const Network*> logged;

  std::vector<Network*> out;
  out.reserve(solids.size());

  for (auto* solid : solids) {
    if (!solid)
      continue;

    auto [it, inserted] = enclosure_cache.emplace(solid, false);
    if (inserted) {
      std::size_t inside_count = 0;
      std::size_t total_count = 0;
      for (const auto* fluid : fluids) {
        if (!fluid)
          continue;
        for (const auto* p : fluid->getBoundaryPoints()) {
          ++total_count;
          if (solid->InsideQ(p->X))
            ++inside_count;
        }
      }
      const double inside_fraction =
          (total_count > 0) ? static_cast<double>(inside_count) / static_cast<double>(total_count) : 0.0;
      it->second = inside_fraction > 0.2;

      if (!logged.contains(solid)) {
        std::cout << Yellow << "[structure_contact] body " << solid->getName()
                  << ": inside_fraction=" << inside_fraction
                  << (it->second ? " -> enclosure, excluded from penetration rejection"
                                 : " -> obstacle, monitored for penetration")
                  << colorReset << std::endl;
        logged.insert(solid);
      }
    }

    if (!it->second)
      out.push_back(solid);
  }

  return out;
}

[[nodiscard]] PenetrationReport detectStructurePenetration(const std::vector<Network*>& fluids,
                                                           const std::vector<Network*>& solids,
                                                           const bool check_midpoints) {
  PenetrationReport report;
  for (const auto* water : fluids) {
    for (const auto* p : water->getBoundaryPoints()) {
      const double h_local = localEdgeLength(p);
      const double tolerance = penetrationTolerance(h_local);
      for (const auto* solid : solids) {
        if (!solid->InsideQ(p->X))
          continue;
        auto [near_f, near_x] = solid->Nearest(p->X);
        const double dist = near_f ? Norm(p->X - near_x) : 0.0;
        if (dist <= tolerance)
          continue;
        ++report.point_count;
        if (dist > report.max_distance) {
          report.max_distance = dist;
          report.worst_tolerance = tolerance;
          report.worst_local_length = h_local;
          report.worst_body = solid->getName();
        }
        break;
      }
    }
    if (check_midpoints) {
      // Midpoint DOFs exist only on true-quadratic elements.
      for (const auto* l : water->getBoundaryLines()) {
        const bool has_true_quad = std::ranges::any_of(l->getBoundaryFaces(), [](const auto* f) { return f->isTrueQuadraticElement; });
        if (!has_true_quad)
          continue;
        const double h_local = localEdgeLength(l);
        const double tolerance = penetrationTolerance(h_local);
        for (const auto* solid : solids) {
          if (!solid->InsideQ(l->X_mid))
            continue;
          auto [near_f, near_x] = solid->Nearest(l->X_mid);
          const double dist = near_f ? Norm(l->X_mid - near_x) : 0.0;
          if (dist <= tolerance)
            continue;
          ++report.face_count;
          if (dist > report.max_distance) {
            report.max_distance = dist;
            report.worst_tolerance = tolerance;
            report.worst_local_length = h_local;
            report.worst_body = solid->getName();
          }
          break;
        }
      }
    }
  }
  return report;
}

void throwIfStructurePenetrated(const std::vector<Network*>& fluids,
                                const std::vector<Network*>& solids,
                                const int time_step,
                                const std::string& stage,
                                const bool check_midpoints) {
  const auto active_solids = penetrationCandidateSolids(fluids, solids);
  if (active_solids.empty())
    return;
  const auto report = detectStructurePenetration(fluids, active_solids, check_midpoints);
  if (report.point_count == 0 && report.face_count == 0)
    return;
  throw step_failure("structure penetration detected at " + stage +
                     " on time_step " + std::to_string(time_step) +
                     " (points=" + std::to_string(report.point_count) +
                     ", mids=" + std::to_string(report.face_count) +
                     ", max_distance=" + std::to_string(report.max_distance) +
                     (report.worst_body.empty() ? "" : ", body=" + report.worst_body) +
                     ", tolerance=" + std::to_string(report.worst_tolerance) +
                     ", h_local=" + std::to_string(report.worst_local_length) + ")");
}

} // namespace BEM_Penetration
