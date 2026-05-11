#pragma once

// Debug / diagnostic utilities extracted from main_time_domain.cpp.
// All functions live in an anonymous namespace so they are file-local
// and behave identically to the original code.

#include <algorithm>
#include <array>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <execinfo.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <unordered_set>
#include <vector>

#include "Network.hpp"
#include "BEM_node_face_state.hpp"

namespace {

// ============================================================================
//  Debug corner-point inspection (BEM_DEBUG_CORNER_POINT / BEM_DEBUG_CORNER_TOL)
// ============================================================================

struct DebugCornerTarget {
  bool enabled = false;
  Tddd x = {0., 0., 0.};
  double tol = 5e-2;
};

const DebugCornerTarget& getDebugCornerTarget() {
  static const DebugCornerTarget target = [] {
    DebugCornerTarget t;
    if (const char* env = std::getenv("BEM_DEBUG_CORNER_POINT")) {
      std::stringstream ss(env);
      char comma0 = 0, comma1 = 0;
      if (ss >> t.x[0] >> comma0 >> t.x[1] >> comma1 >> t.x[2]) {
        if (comma0 == ',' && comma1 == ',')
          t.enabled = true;
      }
    }
    if (const char* env_tol = std::getenv("BEM_DEBUG_CORNER_TOL")) {
      try {
        t.tol = std::stod(env_tol);
      } catch (...) {
      }
    }
    return t;
  }();
  return target;
}

networkPoint* findDebugCornerPoint(Network* water, double* nearest_dist = nullptr) {
  const auto& target = getDebugCornerTarget();
  if (!target.enabled || water == nullptr)
    return nullptr;

  networkPoint* best = nullptr;
  double best_dist = std::numeric_limits<double>::max();
  for (auto* p : water->getBoundaryPoints()) {
    const double d = Norm(p->X - target.x);
    if (d < best_dist) {
      best_dist = d;
      best = p;
    }
  }
  if (nearest_dist)
    *nearest_dist = best_dist;
  if (best_dist > target.tol)
    return nullptr;
  return best;
}

std::string faceKeySummary(networkFace* f) {
  if (!f)
    return "nullptr";
  auto [p0, p1, p2] = f->getPoints();
  std::array<int, 3> ids = {p0 ? p0->index : -1, p1 ? p1->index : -1, p2 ? p2->index : -1};
  const bool face_dirichlet = p0 && p1 && p2 &&
                              isDirichletBoundaryState(p0, f) &&
                              isDirichletBoundaryState(p1, f) &&
                              isDirichletBoundaryState(p2, f);
  const bool face_neumann = p0 && p1 && p2 &&
                            isNeumannBoundaryState(p0, f) &&
                            isNeumannBoundaryState(p1, f) &&
                            isNeumannBoundaryState(p2, f);
  std::sort(ids.begin(), ids.end());
  std::ostringstream oss;
  oss << "ptr=" << f
      << " idx=" << f->index
      << " flags={D:" << face_dirichlet << ",N:" << face_neumann << "}"
      << " pts={" << ids[0] << "," << ids[1] << "," << ids[2] << "}"
      << " centroid=" << f->centroid
      << " normal=" << f->normal;
  return oss.str();
}

template <class MapType>
void dumpFaceValueMap(const char* label,
                      const MapType& m,
                      const std::unordered_set<networkFace*>& boundary_face_set) {
  std::vector<std::pair<networkFace*, double>> entries;
  entries.reserve(m.size());
  for (const auto& [f, v] : m)
    entries.emplace_back(f, v);
  std::stable_sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
    return reinterpret_cast<std::uintptr_t>(a.first) < reinterpret_cast<std::uintptr_t>(b.first);
  });
  std::cout << "  " << label << " size=" << entries.size() << std::endl;
  for (const auto& [f, v] : entries) {
    const bool alive = (f == nullptr) || boundary_face_set.count(f);
    std::cout << "    key={" << faceKeySummary(f) << "} alive=" << alive
              << " value=" << std::setprecision(17) << v << std::endl;
  }
}

void dumpDebugCornerPointState(Network* water, const char* stage, int time_step, int rk_step) {
  const auto& target = getDebugCornerTarget();
  if (!target.enabled || water == nullptr)
    return;

  double nearest_dist = std::numeric_limits<double>::max();
  auto* p = findDebugCornerPoint(water, &nearest_dist);
  if (!p) {
    std::cout << Magenta << "[corner_point_debug] stage=" << stage
              << " water=" << water->getName()
              << " time_step=" << time_step
              << " rk_step=" << rk_step
              << " target=" << target.x
              << " no point within tol=" << target.tol
              << " nearest_dist=" << nearest_dist
              << colorReset << std::endl;
    return;
  }

  auto boundary_faces = p->getBoundaryFaces();
  std::unordered_set<networkFace*> boundary_face_set(boundary_faces.begin(), boundary_faces.end());
  std::stable_sort(boundary_faces.begin(), boundary_faces.end(), [](const auto* a, const auto* b) {
    return reinterpret_cast<std::uintptr_t>(a) < reinterpret_cast<std::uintptr_t>(b);
  });

  std::cout << Magenta << "[corner_point_debug] stage=" << stage
            << " water=" << water->getName()
            << " time_step=" << time_step
            << " rk_step=" << rk_step
            << " target=" << target.x
            << " tol=" << target.tol
            << " point_ptr=" << p
            << " point_index=" << p->index
            << " pos=" << p->X
            << " dist=" << nearest_dist
            << " flags={D:" << p->Dirichlet << ",N:" << p->Neumann << ",C:" << p->BCInterface << ",M:" << p->isMultipleNode << "}"
            << " face_count=" << boundary_faces.size()
            << " contact_count=" << p->getContactFaces().size()
            << colorReset << std::endl;

  std::cout << "  boundary_faces:" << std::endl;
  for (auto* f : boundary_faces) {
    auto* cf = p->getNearestContactFace(f);
    Tddd cx = cf ? Nearest(p->getPosition(), ToX(cf)) : Tddd{0., 0., 0.};
    double dist = cf ? Norm(p->getPosition() - cx) : 1E+20;
    std::cout << "    {" << faceKeySummary(f) << "} nearest_contact={"
              << faceKeySummary(cf) << "} nearest_x=" << cx
              << " nearest_dist=" << dist << std::endl;
  }

  std::vector<std::pair<networkFace*, int>> dof_entries;
  for (const auto& [f, d] : p->dofs)
    if (d.index >= 0)
      dof_entries.emplace_back(f, d.index);
  std::stable_sort(dof_entries.begin(), dof_entries.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
  std::cout << "  dofs size=" << dof_entries.size() << std::endl;
  for (const auto& [f, idx] : dof_entries) {
    const bool alive = (f == nullptr) || boundary_face_set.count(f);
    std::cout << "    dof idx=" << idx << " key={" << faceKeySummary(f) << "} alive=" << alive << std::endl;
  }

  std::cout << "  dofs:" << std::endl;
  for (const auto& [f, d] : p->dofs) {
    const bool alive = (f == nullptr) || boundary_face_set.count(f);
    std::cout << "    key={" << faceKeySummary(f) << "} alive=" << alive
              << " phi=" << d.phi << " phin=" << d.phin
              << " phi_t=" << d.phi_t << " phin_t=" << d.phin_t << std::endl;
  }
}

// ============================================================================
//  Crash backtrace handler (BEM_BACKTRACE=1)
// ============================================================================

void crashBacktraceHandler(int sig) {
  const char* sig_name = strsignal(sig);
  if (!sig_name)
    sig_name = "UNKNOWN";

  ::write(STDERR_FILENO, "\n=== crash ===\n", 13);
  ::write(STDERR_FILENO, sig_name, std::strlen(sig_name));
  ::write(STDERR_FILENO, "\n", 1);

  void* frames[128];
  const int n = ::backtrace(frames, 128);
  ::backtrace_symbols_fd(frames, n, STDERR_FILENO);
  ::write(STDERR_FILENO, "=== end ===\n", 12);

  std::_Exit(128 + sig);
}

void installCrashBacktraceIfRequested() {
  const char* env = std::getenv("BEM_BACKTRACE");
  if (!env || std::strcmp(env, "1") != 0)
    return;

  std::signal(SIGSEGV, crashBacktraceHandler);
  std::signal(SIGABRT, crashBacktraceHandler);
  std::signal(SIGBUS, crashBacktraceHandler);
}

// ============================================================================
//  FMM real-field default for time-domain runs
// ============================================================================

void enableRealFieldFmmDefaultForTimeDomain() {
  if (std::getenv("BEM_FMM_REALFIELD_M_CONJ") != nullptr)
    return;
#if defined(_WIN32)
  _putenv_s("BEM_FMM_REALFIELD_M_CONJ", "1");
#else
  ::setenv("BEM_FMM_REALFIELD_M_CONJ", "1", 0 /* overwrite */);
#endif
}

// ============================================================================
//  Corner-connected Neumann line debug log (BEM_CORNER_DEBUG)
//  Extracted from the lambda log_corner_connected_neumann_lines_after_boundary_types
// ============================================================================

void logCornerConnectedNeumannLinesAfterBoundaryTypes(const Network* water, const char* phase) {
  if (!water)
    return;
  auto local_line_length = [](const networkLine* const l) {
    static std::unordered_set<networkPoint*> adjacent_points;
    adjacent_points.clear();
    for (const auto& f : l->getBoundaryFaces()) {
      auto points = f->getPoints();
      adjacent_points.insert(points.begin(), points.end());
    }
    static std::unordered_set<networkLine*> adjacent_lines;
    adjacent_lines.clear();
    for (const auto& p : adjacent_points)
      for (const auto& L : p->getBoundaryLines())
        if (L != l)
          adjacent_lines.insert(L);
    if (adjacent_lines.empty())
      return 0.0;
    auto ret = 0.;
    for (const auto& L : adjacent_lines)
      ret += L->length();
    return ret / adjacent_lines.size();
  };
  static const bool enable_corner_neumann_debug = (std::getenv("BEM_CORNER_DEBUG") && std::string(std::getenv("BEM_CORNER_DEBUG")) != "0");
  if (enable_corner_neumann_debug) {
    auto is_corner_connected_neumann_line = [](const networkLine* l) {
      if (!l || !l->Neumann || l->BCInterface)
        return false;
      auto [p0, p1] = l->getPoints();
      return p0 && p1 && (p0->BCInterface || p1->BCInterface);
    };
    std::vector<networkLine*> candidates;
    for (auto* l : water->getBoundaryLines())
      if (is_corner_connected_neumann_line(l))
        candidates.emplace_back(l);
    std::ranges::sort(candidates, [](const auto* a, const auto* b) { return a->length() < b->length(); });
    if (candidates.size() > 5)
      candidates.resize(5);

    for (auto* l : candidates) {
      auto [p0, p1] = l->getPoints();
      const auto faces = l->getBoundaryFaces();
      const double len = l->length();
      const double local_mean_len = local_line_length(l);
      double alt0 = -1.0, alt1 = -1.0, alt_threshold = -1.0;
      double alt_ratio0 = -1.0, alt_ratio1 = -1.0;
      double aspect_ratio0 = -1.0, aspect_ratio1 = -1.0;
      double min_angle_deg0 = -1.0, min_angle_deg1 = -1.0;
      double area0 = -1.0, area1 = -1.0;
      double mean_area0 = -1.0, mean_area1 = -1.0;
      int common_points = -1;
      int opp0_lines = -1, opp1_lines = -1;
      bool face0_dirichlet = false, face0_neumann = false;
      bool face1_dirichlet = false, face1_neumann = false;
      if (faces.size() == 2 && faces[0] && faces[1] && p0 && p1) {
        auto* f0 = faces[0];
        auto* f1 = faces[1];
        face0_dirichlet = isDirichletBoundaryState(p0, f0) && isDirichletBoundaryState(p1, f0);
        face0_neumann = isNeumannBoundaryState(p0, f0) && isNeumannBoundaryState(p1, f0);
        face1_dirichlet = isDirichletBoundaryState(p0, f1) && isDirichletBoundaryState(p1, f1);
        face1_neumann = isNeumannBoundaryState(p0, f1) && isNeumannBoundaryState(p1, f1);
        auto face_shape_metrics = [](const networkFace* f) {
          std::array<double, 2> out{-1.0, -1.0};
          if (!f)
            return out;
          const auto pts = f->getPoints();
          const Tddd a = pts[0]->X, b = pts[1]->X, c = pts[2]->X;
          const double l01 = Norm(a - b);
          const double l12 = Norm(b - c);
          const double l20 = Norm(c - a);
          const double max_edge = std::max({l01, l12, l20});
          const double area = boundaryFaceArea(const_cast<networkFace*>(f));
          const double altitude = (max_edge > 1e-20) ? 2.0 * area / max_edge : 0.0;
          out[0] = (altitude > 1e-20) ? max_edge / altitude : 1E+100;
          auto clamp_cos = [](double x) { return std::max(-1.0, std::min(1.0, x)); };
          const double ang0 = std::acos(clamp_cos(Dot((b - a) / l01, (c - a) / l20))) * 180.0 / M_PI;
          const double ang1 = std::acos(clamp_cos(Dot((a - b) / l01, (c - b) / l12))) * 180.0 / M_PI;
          const double ang2 = 180.0 - ang0 - ang1;
          out[1] = std::min({ang0, ang1, ang2});
          return out;
        };
        auto [a, this0, b, l1_var, p2_var, l2_var] = f0->getPointsAndLines(l);
        auto [q0, this1, q1, e1, q2, e2] = f1->getPointsAndLines(l);
        if (this0 == l && this1 == l && a == q1 && b == q0) {
          area0 = boundaryFaceArea(f0);
          area1 = boundaryFaceArea(f1);
          const double max_edge0 = std::max({Norm(f0->getPoints()[0]->X - f0->getPoints()[1]->X),
                                             Norm(f0->getPoints()[1]->X - f0->getPoints()[2]->X),
                                             Norm(f0->getPoints()[2]->X - f0->getPoints()[0]->X)});
          const double max_edge1 = std::max({Norm(f1->getPoints()[0]->X - f1->getPoints()[1]->X),
                                             Norm(f1->getPoints()[1]->X - f1->getPoints()[2]->X),
                                             Norm(f1->getPoints()[2]->X - f1->getPoints()[0]->X)});
          alt0 = (max_edge0 > 1e-20) ? 2.0 * area0 / max_edge0 : 0.0;
          alt1 = (max_edge1 > 1e-20) ? 2.0 * area1 / max_edge1 : 0.0;
          alt_threshold = (local_mean_len > 0.0) ? 0.1 * local_mean_len : -1.0;
          alt_ratio0 = (local_mean_len > 0.0) ? alt0 / local_mean_len : -1.0;
          alt_ratio1 = (local_mean_len > 0.0) ? alt1 / local_mean_len : -1.0;
          auto shape0 = face_shape_metrics(f0);
          auto shape1 = face_shape_metrics(f1);
          aspect_ratio0 = shape0[0];
          aspect_ratio1 = shape1[0];
          min_angle_deg0 = shape0[1];
          min_angle_deg1 = shape1[1];
          mean_area0 = localMeanFaceArea(f0);
          mean_area1 = localMeanFaceArea(f1);
          common_points = static_cast<int>(Intersection(p0->getNeighborPointsOnSurfaces(),
                                                        p1->getNeighborPointsOnSurfaces())
                                               .size());
          opp0_lines = p2_var ? static_cast<int>(p2_var->getLines().size()) : -1;
          opp1_lines = q2 ? static_cast<int>(q2->getLines().size()) : -1;
        }
      }
      std::cout << Magenta << "[corner_neumann_debug] " << phase
                << " p0=" << (p0 ? p0->index : -1)
                << " p1=" << (p1 ? p1->index : -1)
                << " len=" << len
                << " local_mean_len=" << local_mean_len
                << " faces=" << faces.size()
                << " alt0=" << alt0
                << " alt1=" << alt1
                << " alt_threshold=" << alt_threshold
                << " alt_ratio0=" << alt_ratio0
                << " alt_ratio1=" << alt_ratio1
                << " aspect_ratio0=" << aspect_ratio0
                << " aspect_ratio1=" << aspect_ratio1
                << " min_angle_deg0=" << min_angle_deg0
                << " min_angle_deg1=" << min_angle_deg1
                << " area0=" << area0
                << " area1=" << area1
                << " local_mean_area0=" << mean_area0
                << " local_mean_area1=" << mean_area1
                << " common_points=" << common_points
                << " line_flags={D:" << l->Dirichlet << ",N:" << l->Neumann << ",C:" << l->BCInterface << "}"
                << " face0_flags={D:" << face0_dirichlet << ",N:" << face0_neumann << "}"
                << " face1_flags={D:" << face1_dirichlet << ",N:" << face1_neumann << "}"
                << " endpoint0_flags={D:" << (p0 ? p0->Dirichlet : false)
                << ",N:" << (p0 ? p0->Neumann : false)
                << ",C:" << (p0 ? p0->BCInterface : false) << "}"
                << " endpoint1_flags={D:" << (p1 ? p1->Dirichlet : false)
                << ",N:" << (p1 ? p1->Neumann : false)
                << ",C:" << (p1 ? p1->BCInterface : false) << "}"
                << " endpoint_corner={" << (p0 ? p0->BCInterface : false) << "," << (p1 ? p1->BCInterface : false) << "}"
                << " endpoint_lines={" << (p0 ? p0->getLines().size() : 0) << "," << (p1 ? p1->getLines().size() : 0) << "}"
                << " opposite_lines={" << opp0_lines << "," << opp1_lines << "}"
                << " x0=" << (p0 ? p0->X : Tddd{0., 0., 0.})
                << " x1=" << (p1 ? p1->X : Tddd{0., 0., 0.})
                << colorReset << std::endl;
    }
  } // if (enable_corner_neumann_debug)
}

} // namespace
