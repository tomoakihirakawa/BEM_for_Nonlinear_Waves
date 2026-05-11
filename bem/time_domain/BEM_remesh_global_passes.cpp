// BEM_remesh_global_passes.cpp
//
// Shared Trial remesh polish passes.

#define BEM
#include "Network.hpp"
#include "BEM_remesh_common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

using FlipMode = SimulationSettings::RemeshingSettings::FlipMode;
using SmoothMode = SimulationSettings::RemeshingSettings::SmoothMode;

namespace {

constexpr int kValenceScoreInf = std::numeric_limits<int>::max();

int sq_int(int x) { return x * x; }

std::pair<int, int> legacyVertexValenceContribution(
    const networkPoint* p, int delta, int s_mean) {
  if (!p)
    return {kValenceScoreInf, kValenceScoreInf};
  const int s = static_cast<int>(p->getLines().size());
  return {sq_int(s - s_mean), sq_int(s + delta - s_mean)};
}

bool hasIncidentSharpEdge(const networkPoint* p, double feature_angle) {
  if (!p)
    return false;
  for (auto* l : p->getBoundaryLines())
    if (l && l->SharpQ(feature_angle))
      return true;
  return false;
}

networkLine* firstIncidentSharpEdge(networkPoint* p, double feature_angle) {
  if (!p)
    return nullptr;
  for (auto* l : p->getBoundaryLines())
    if (l && l->SharpQ(feature_angle))
      return l;
  return nullptr;
}

struct SharpBoundarySector {
  std::vector<networkFace*> faces;
  int target = 1;
  double angle = 0.0;
};

bool closeSharpBoundarySector(const networkPoint* p,
                              SharpBoundarySector& current,
                              std::vector<SharpBoundarySector>& sectors) {
  if (current.faces.empty())
    return true;
  current.target = std::max(1, static_cast<int>(std::lround(current.angle / (M_PI / 3.0))));
  if (!std::isfinite(current.angle) || !(current.angle > 0.0))
    return false;
  sectors.push_back(current);
  current = {};
  (void)p;
  return true;
}

bool buildSharpBoundarySectors(networkPoint* p,
                               double feature_angle,
                               std::vector<SharpBoundarySector>& sectors) {
  sectors.clear();
  auto* start_line = firstIncidentSharpEdge(p, feature_angle);
  if (!start_line)
    return false;

  try {
    const auto lines = p->getBoundaryLinesSort(start_line);
    const auto faces = p->getBoundaryFacesSort(start_line);
    if (lines.empty() || faces.empty() || lines.size() < faces.size())
      return false;

    SharpBoundarySector current;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      auto* f = faces[i];
      if (!f)
        return false;
      current.faces.push_back(f);
      current.angle += f->getAngle(p);

      networkLine* next_line = nullptr;
      if (i + 1 < lines.size())
        next_line = lines[i + 1];
      else if (lines.size() == faces.size())
        next_line = lines.front();

      if (next_line && next_line->SharpQ(feature_angle))
        if (!closeSharpBoundarySector(p, current, sectors))
          return false;
    }
    if (!closeSharpBoundarySector(p, current, sectors))
      return false;
  } catch (...) {
    return false;
  }
  return !sectors.empty();
}

int findSharpBoundarySector(const std::vector<SharpBoundarySector>& sectors,
                            const networkFace* face) {
  if (!face)
    return -1;
  for (std::size_t i = 0; i < sectors.size(); ++i)
    if (std::find(sectors[i].faces.begin(), sectors[i].faces.end(), face) != sectors[i].faces.end())
      return static_cast<int>(i);
  return -1;
}

bool vertexSharpSectorValenceContribution(networkPoint* p,
                                          int delta,
                                          const std::vector<networkFace*>& affected_faces,
                                          double feature_angle,
                                          int s_mean,
                                          int& before_out,
                                          int& after_out,
                                          bool& used_sector_out) {
  used_sector_out = false;
  auto legacy = [&]() {
    auto [b, a] = legacyVertexValenceContribution(p, delta, s_mean);
    before_out = b;
    after_out = a;
    return b != kValenceScoreInf && a != kValenceScoreInf;
  };

  if (!hasIncidentSharpEdge(p, feature_angle))
    return legacy();

  std::vector<SharpBoundarySector> sectors;
  if (!buildSharpBoundarySectors(p, feature_angle, sectors))
    return legacy();

  int affected_sector = -1;
  for (auto* f : affected_faces) {
    const int si = findSharpBoundarySector(sectors, f);
    if (si < 0)
      return legacy();
    if (affected_sector < 0)
      affected_sector = si;
    else if (affected_sector != si)
      return legacy();
  }
  if (affected_sector < 0)
    return legacy();

  int before = 0;
  int after = 0;
  for (std::size_t i = 0; i < sectors.size(); ++i) {
    const int count = static_cast<int>(sectors[i].faces.size());
    const int target = sectors[i].target;
    before += sq_int(count - target);
    after += sq_int(count + (static_cast<int>(i) == affected_sector ? delta : 0) - target);
  }
  before_out = before;
  after_out = after;
  used_sector_out = true;
  return true;
}

} // namespace

// Compare vertex valence balance before and after a flip.
// Return {v_before, v_after}; smaller is better.
std::pair<int, int> valenceDeviationScore(const networkLine* l, int s_mean) {
  if (!l)
    return {kValenceScoreInf, kValenceScoreInf};
  auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1])
    return {kValenceScoreInf, kValenceScoreInf};

  auto [p0, p1] = l->getPoints();
  if (!p0 || !p1)
    return {kValenceScoreInf, kValenceScoreInf};

  auto* p2 = faces[0]->getPointOpposite(l);
  auto* p3 = faces[1]->getPointOpposite(l);
  if (!p2 || !p3 || p2 == p3)
    return {kValenceScoreInf, kValenceScoreInf};

  int s0 = p0->getLines().size();
  int s1 = p1->getLines().size();
  int s2 = p2->getLines().size();
  int s3 = p3->getLines().size();

  // After a flip p0 and p1 lose one incident edge, so degree <= 3 is invalid.
  if (s0 <= 3 || s1 <= 3)
    return {kValenceScoreInf, kValenceScoreInf};

  int v_before = sq_int(s0 - s_mean) + sq_int(s1 - s_mean) + sq_int(s2 - s_mean) + sq_int(s3 - s_mean);
  int v_after = sq_int(s0 - 1 - s_mean) + sq_int(s1 - 1 - s_mean) + sq_int(s2 + 1 - s_mean) + sq_int(s3 + 1 - s_mean);
  return {v_before, v_after};
}

std::pair<int, int> sharpSectorValenceDeviationScore(
    const networkLine* l, double feature_angle, int s_mean) {
  if (!l || l->BCInterface || l->SharpQ(feature_angle))
    return {kValenceScoreInf, kValenceScoreInf};

  auto faces = l->getBoundaryFaces();
  if (faces.size() != 2 || !faces[0] || !faces[1])
    return {kValenceScoreInf, kValenceScoreInf};

  auto [p0, p1] = l->getPoints();
  if (!p0 || !p1)
    return {kValenceScoreInf, kValenceScoreInf};

  auto* p2 = faces[0]->getPointOpposite(l);
  auto* p3 = faces[1]->getPointOpposite(l);
  if (!p2 || !p3 || p2 == p3)
    return {kValenceScoreInf, kValenceScoreInf};

  if (static_cast<int>(p0->getLines().size()) <= 3 ||
      static_cast<int>(p1->getLines().size()) <= 3)
    return {kValenceScoreInf, kValenceScoreInf};

  int before = 0;
  int after = 0;
  bool used_any_sector = false;
  auto add_vertex = [&](networkPoint* p, int delta, std::vector<networkFace*> affected) {
    int b = 0;
    int a = 0;
    bool used_sector = false;
    if (!vertexSharpSectorValenceContribution(
            p, delta, affected, feature_angle, s_mean, b, a, used_sector))
      return false;
    before += b;
    after += a;
    used_any_sector = used_any_sector || used_sector;
    return true;
  };

  if (!add_vertex(p0, -1, {faces[0], faces[1]})) return {kValenceScoreInf, kValenceScoreInf};
  if (!add_vertex(p1, -1, {faces[0], faces[1]})) return {kValenceScoreInf, kValenceScoreInf};
  if (!add_vertex(p2, +1, {faces[0]})) return {kValenceScoreInf, kValenceScoreInf};
  if (!add_vertex(p3, +1, {faces[1]})) return {kValenceScoreInf, kValenceScoreInf};
  if (!used_any_sector)
    return valenceDeviationScore(l, s_mean);
  return {before, after};
}

std::size_t pass_flip(Network& water,
                      SimulationSettings::RemeshingSettings::FlipMode flip_mode,
                      double feature_angle,
                      double cos_thr,
                      double min_angle_rad,
                      bool use_sharp_sector_valence) {
  struct Cand {
    networkLine* l;
    int vg;
    double minang;
  };
  const bool relaxed = (flip_mode == FlipMode::ValenceRelaxed);
  const bool delaunay_required = (flip_mode != FlipMode::ValenceOnly);
  const bool valence_required = (flip_mode != FlipMode::DelaunayOnly);
  const double vertex_cos_thr =
      (cos_thr > std::cos(30.0 * M_PI / 180.0))
          ? cos_thr
          : std::cos(30.0 * M_PI / 180.0);

  auto passes_all_guards = [&](networkLine* l, int& vg_out) -> bool {
    if (!l || !line_alive(water, l))
      return false;
    if (l->BCInterface)
      return false;
    if (l->SharpQ(feature_angle))
      return false;
    if (l->getBoundaryFaces().size() != 2)
      return false;
    if (delaunay_required && !is_non_delaunay(l))
      return false;
    auto [vb, va] = use_sharp_sector_valence
                        ? sharpSectorValenceDeviationScore(l, feature_angle, 6)
                        : valenceDeviationScore(l, 6);
    if (valence_required && (vb == kValenceScoreInf || va == kValenceScoreInf))
      return false;
    vg_out = vb - va;
    if (valence_required) {
      if (relaxed ? (vg_out < 0) : (vg_out <= 0))
        return false;
    }
    if (!flip_preserves_normals(l, cos_thr, min_angle_rad))
      return false;
    if (!flip_preserves_vertex_normals(l, vertex_cos_thr))
      return false;
    if (!flip_dihedral_change_ok(l))
      return false;
    if (!flip_midpoint_drift_ok(l))
      return false;
    return true;
  };

  const int max_sweeps = relaxed ? 5 : 4;
  std::size_t total = 0;
  for (int iter = 0; iter < max_sweeps; ++iter) {
    std::vector<Cand> cands;
    cands.reserve(water.getLines().size() / 4);
    for (auto* l : water.getLines()) {
      int vg = 0;
      if (!passes_all_guards(l, vg))
        continue;
      cands.push_back({l, vg, flip_new_min_angle(l)});
    }
    if (cands.empty())
      break;
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
      if (a.vg != b.vg)
        return a.vg > b.vg;
      return a.minang > b.minang;
    });
    std::size_t round = 0;
    for (const auto& c : cands) {
      int vg2 = 0;
      if (!passes_all_guards(c.l, vg2))
        continue;
      try {
        if (c.l->Flip(true)) {
          ++round;
          ++total;
        }
      } catch (...) {
      }
    }
    if (round == 0)
      break;
  }
  water.setGeometricPropertiesForce();
  return total;
}

std::size_t pass_smooth(Network& water,
                        SimulationSettings::RemeshingSettings::SmoothMode smooth_mode,
                        const SizingField& sf,
                        int sub_iterations,
                        double feature_angle,
                        double cos_thr,
                        double min_angle_rad,
                        bool skip_sharp_points) {
  std::size_t moved = 0;
  auto touches_sharp_line = [&](const networkPoint* p) {
    if (!p)
      return false;
    for (auto* l : p->getBoundaryLines())
      if (l && l->SharpQ(feature_angle))
        return true;
    return false;
  };
  for (int si = 0; si < sub_iterations; ++si) {
    std::vector<std::pair<networkPoint*, Tddd>> updates;
    for (auto* p : water.getBoundaryPoints()) {
      if (!allows_free_surface_smoothing(p))
        continue;
      if (skip_sharp_points && touches_sharp_line(p))
        continue;
      Tddd V;
      if (!remesh_smooth_delta(p, smooth_mode, sf, feature_angle, V))
        continue;
      if (!isFinite(V))
        continue;
      const double disp = Norm(V);
      const double limit = remesh_smooth_step_limit_ratio(smooth_mode) * localEdgeLength(p);
      if (disp > 1e-14) {
        const Tddd V_capped = (disp <= limit) ? V : ((limit / disp) * V);
        const double alpha = find_safe_smooth_step(p, V_capped, cos_thr, min_angle_rad);
        if (alpha <= 0.0)
          continue;
        updates.emplace_back(p, p->X + alpha * V_capped);
      }
    }
    for (auto& [p, x] : updates) {
      p->setXSingle(x);
      ++moved;
    }
    water.setGeometricPropertiesForce();
  }
  return moved;
}
