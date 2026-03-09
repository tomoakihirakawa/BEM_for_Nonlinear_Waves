#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "BEM_freqency_domain.hpp"
#include "basic_vectors.hpp"
#include "rootFinding.hpp"

namespace bem_frequency_domain {

struct FaceField {
  std::array<Complex, 3> phi{};
  std::array<Complex, 3> phin{};
};

struct LinearSolution {
  double omega = 0.0;
  std::unordered_map<const networkFace *, FaceField> face_field;
  std::unordered_map<const networkFace *, std::array<Complex, 3>> grad_phi;
};

struct WaterlineSegment {
  networkPoint *p0 = nullptr;
  networkPoint *p1 = nullptr;
  networkFace *free_face = nullptr;
  Tddd normal{};
  Tddd midpoint{};
  double length = 0.0;
  double nh = 0.0;
};

struct QTFResult {
  std::vector<double> omegas;
  std::vector<std::vector<std::array<Complex, 6>>> qminus;
  std::vector<std::vector<std::array<Complex, 6>>> qplus;
};

struct QTFReport {
  double max_abs_qminus = 0.0;
  double max_abs_qplus = 0.0;
  double max_sym_abs_qminus = 0.0;
  double max_sym_abs_qplus = 0.0;
  double max_sym_rel_qminus = 0.0;
  double max_sym_rel_qplus = 0.0;
};

struct NewmanReport {
  double max_abs = 0.0;
  double max_rel = 0.0;
};

struct QTFPair {
  double omega_m = 0.0;
  double omega_n = 0.0;
  std::array<Complex, 6> qminus{};
  std::array<Complex, 6> qplus{};
};

inline std::array<Complex, 6> conjugate_qtf(const std::array<Complex, 6> &q) {
  std::array<Complex, 6> out = q;
  for (auto &v : out)
    v = std::conj(v);
  return out;
}

inline WaterWaveTheory make_wave_for_omega(const WaterWaveTheory &base, double omega) {
  WaterWaveTheory wave = base;
  if (omega > 0.0)
    wave.set_w_h(omega, base.h);
  return wave;
}

inline Complex incident_phi_hat(const WaterWaveTheory &wave, const Tddd &x) {
  if (!(wave.w > 0.0))
    return Complex{0.0, 0.0};
  const double kx = wave.k * std::cos(wave.theta);
  const double ky = wave.k * std::sin(wave.theta);
  const double z = wave.clampZ(x[2] - wave.h - wave.bottom_z, wave.k);
  const double amp = wave.A * _GRAVITY_ / wave.w;
  const Complex phase = std::exp(Complex{0.0, kx * x[0] + ky * x[1] + wave.phase_shift});
  return Complex{0.0, -1.0} * (amp * cosh_kzh_cosh_kh(wave.k, wave.h, z)) * phase;
}

inline std::array<Complex, 3> incident_grad_phi_hat(const WaterWaveTheory &wave, const Tddd &x) {
  if (!(wave.w > 0.0))
    return {Complex{0.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}};
  const double kx = wave.k * std::cos(wave.theta);
  const double ky = wave.k * std::sin(wave.theta);
  const double z = wave.clampZ(x[2] - wave.h - wave.bottom_z, wave.k);
  const double amp = wave.A * _GRAVITY_ / wave.w;
  const Complex phase = std::exp(Complex{0.0, kx * x[0] + ky * x[1] + wave.phase_shift});
  const Complex common = amp * phase;
  const double cosh_fac = cosh_kzh_cosh_kh(wave.k, wave.h, z);
  const double sinh_fac = sinh_kzh_cosh_kh(wave.k, wave.h, z);
  return {
      common * (kx * cosh_fac),
      common * (ky * cosh_fac),
      Complex{0.0, -1.0} * common * (wave.k * sinh_fac),
  };
}

inline Complex incident_phin_hat(const WaterWaveTheory &wave, const Tddd &x, const Tddd &normal) {
  const auto grad = incident_grad_phi_hat(wave, x);
  return grad[0] * normal[0] + grad[1] * normal[1] + grad[2] * normal[2];
}

inline Complex get_phi_on_face(const networkPoint *p, const networkFace *f) {
  if (!p)
    return Complex{0.0, 0.0};
  if (f) {
    auto it = p->phiOnFace.find(const_cast<networkFace *>(f));
    if (it != p->phiOnFace.end())
      return Complex{it->second, 0.0};
  }
  auto it = p->phiOnFace.find(nullptr);
  if (it != p->phiOnFace.end())
    return Complex{it->second, 0.0};
  return Complex{std::get<0>(p->phiphin), 0.0};
}

inline Complex get_phin_on_face(const networkPoint *p, const networkFace *f) {
  if (!p)
    return Complex{0.0, 0.0};
  if (f) {
    auto it = p->phinOnFace.find(const_cast<networkFace *>(f));
    if (it != p->phinOnFace.end())
      return Complex{it->second, 0.0};
  }
  auto it = p->phinOnFace.find(nullptr);
  if (it != p->phinOnFace.end())
    return Complex{it->second, 0.0};
  return Complex{std::get<1>(p->phiphin), 0.0};
}

inline std::array<Complex, 3> face_gradient(const networkFace *f, const FaceField &field) {
  auto [p0, p1, p2] = f->getPoints();
  const std::array<double, 3> phi_re = {field.phi[0].real(), field.phi[1].real(), field.phi[2].real()};
  const std::array<double, 3> phi_im = {field.phi[0].imag(), field.phi[1].imag(), field.phi[2].imag()};
  const Tddd grad_re = gradP1(T3Tddd{p0->X, p1->X, p2->X}, phi_re);
  const Tddd grad_im = gradP1(T3Tddd{p0->X, p1->X, p2->X}, phi_im);

  std::array<Complex, 3> grad = {
      Complex{grad_re[0], grad_im[0]},
      Complex{grad_re[1], grad_im[1]},
      Complex{grad_re[2], grad_im[2]},
  };

  const Complex phin_avg = (field.phin[0] + field.phin[1] + field.phin[2]) / 3.0;
  grad[0] += phin_avg * f->normal[0];
  grad[1] += phin_avg * f->normal[1];
  grad[2] += phin_avg * f->normal[2];
  return grad;
}

inline void populate_face_gradients(LinearSolution &sol) {
  sol.grad_phi.clear();
  sol.grad_phi.reserve(sol.face_field.size());
  for (const auto &[f, field] : sol.face_field) {
    sol.grad_phi.emplace(f, face_gradient(f, field));
  }
}

inline LinearSolution capture_linear_solution(double omega, const std::unordered_set<networkFace *> &faces) {
  LinearSolution sol;
  sol.omega = omega;
  sol.face_field.reserve(faces.size());
  for (auto *f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    FaceField field;
    field.phi = {get_phi_on_face(p0, f), get_phi_on_face(p1, f), get_phi_on_face(p2, f)};
    field.phin = {get_phin_on_face(p0, f), get_phin_on_face(p1, f), get_phin_on_face(p2, f)};
    sol.face_field.emplace(f, field);
  }
  populate_face_gradients(sol);
  return sol;
}

inline LinearSolution combine_linear_solutions(double omega,
                                               const std::vector<std::pair<const LinearSolution *, Complex>> &parts) {
  LinearSolution out;
  out.omega = omega;
  if (parts.empty() || !parts.front().first)
    return out;

  const auto &base = *parts.front().first;
  out.face_field.reserve(base.face_field.size());
  for (const auto &[f, field] : base.face_field) {
    FaceField sum;
    sum.phi = {Complex{0.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}};
    sum.phin = {Complex{0.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}};
    for (const auto &[part, weight] : parts) {
      if (!part)
        continue;
      const auto it = part->face_field.find(f);
      if (it == part->face_field.end())
        continue;
      for (std::size_t i = 0; i < 3; ++i) {
        sum.phi[i] += weight * it->second.phi[i];
        sum.phin[i] += weight * it->second.phin[i];
      }
    }
    out.face_field.emplace(f, sum);
  }
  populate_face_gradients(out);
  return out;
}

struct IdHash {
  std::size_t operator()(const Id &id) const noexcept {
    const auto *p = std::get<0>(id);
    const auto *f = std::get<1>(id);
    const std::size_t hp = std::hash<const void *>{}(p);
    const std::size_t hf = std::hash<const void *>{}(f);
    return hp ^ (hf << 1);
  }
};

struct IdEq {
  bool operator()(const Id &a, const Id &b) const noexcept {
    return std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b);
  }
};

struct SolutionIndex {
  std::unordered_map<Id, std::size_t, IdHash, IdEq> index;
};

inline SolutionIndex build_solution_index(const Solution &sol) {
  SolutionIndex idx;
  idx.index.reserve(sol.n);
  for (std::size_t i = 0; i < sol.n; ++i)
    idx.index.emplace(sol.id_by_index[i], i);
  return idx;
}

inline Complex lookup_phi(const Solution &sol, const SolutionIndex &idx, const networkPoint *p, const networkFace *f) {
  Id key{const_cast<networkPoint *>(p), const_cast<networkFace *>(f)};
  auto it = idx.index.find(key);
  if (it != idx.index.end())
    return sol.phi[it->second];
  key = {const_cast<networkPoint *>(p), nullptr};
  it = idx.index.find(key);
  if (it != idx.index.end())
    return sol.phi[it->second];
  return Complex{0.0, 0.0};
}

inline LinearSolution capture_linear_solution(double omega, const Solution &sol, const std::unordered_set<networkFace *> &faces) {
  std::unordered_map<Id, std::size_t, IdHash, IdEq> index;
  index.reserve(sol.n);
  for (std::size_t i = 0; i < sol.n; ++i)
    index.emplace(sol.id_by_index[i], i);

  auto lookup = [&](const networkPoint *p, const networkFace *f, const std::vector<Complex> &vals) -> Complex {
    Id key{const_cast<networkPoint *>(p), const_cast<networkFace *>(f)};
    auto it = index.find(key);
    if (it != index.end())
      return vals[it->second];
    key = {const_cast<networkPoint *>(p), nullptr};
    it = index.find(key);
    if (it != index.end())
      return vals[it->second];
    return Complex{0.0, 0.0};
  };

  LinearSolution out;
  out.omega = omega;
  out.face_field.reserve(faces.size());
  for (auto *f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    FaceField field;
    field.phi = {lookup(p0, f, sol.phi), lookup(p1, f, sol.phi), lookup(p2, f, sol.phi)};
    field.phin = {lookup(p0, f, sol.phin), lookup(p1, f, sol.phin), lookup(p2, f, sol.phin)};
    out.face_field.emplace(f, field);
  }
  populate_face_gradients(out);
  return out;
}

inline LinearSolution build_incident_solution(const WaterWaveTheory &wave, const std::unordered_set<networkFace *> &faces) {
  LinearSolution sol;
  sol.omega = wave.w;
  sol.face_field.reserve(faces.size());
  for (auto *f : faces) {
    auto [p0, p1, p2] = f->getPoints();
    FaceField field;
    field.phi = {
        incident_phi_hat(wave, p0->X),
        incident_phi_hat(wave, p1->X),
        incident_phi_hat(wave, p2->X),
    };
    field.phin = {
        incident_phin_hat(wave, p0->X, f->normal),
        incident_phin_hat(wave, p1->X, f->normal),
        incident_phin_hat(wave, p2->X, f->normal),
    };
    sol.face_field.emplace(f, field);
  }
  populate_face_gradients(sol);
  return sol;
}

inline std::vector<WaterlineSegment> collect_waterline_segments(const Network &water,
                                                                const std::unordered_set<networkFace *> &float_faces,
                                                                const std::unordered_set<networkFace *> &free_faces) {
  std::vector<WaterlineSegment> segments;
  for (auto *line : water.getLines()) {
    if (!line)
      continue;
    const auto &fs = line->getBoundaryFaces();
    if (fs.size() < 2)
      continue;
    networkFace *float_face = nullptr;
    networkFace *free_face = nullptr;
    for (auto *f : fs) {
      if (float_faces.contains(f))
        float_face = f;
      if (free_faces.contains(f))
        free_face = f;
    }
    if (!float_face || !free_face)
      continue;
    auto [p0, p1] = line->getPoints();
    if (!p0 || !p1)
      continue;
    const Tddd dx = p1->X - p0->X;
    const double len = Norm(dx);
    if (!(len > 0.0))
      continue;
    WaterlineSegment seg;
    seg.p0 = p0;
    seg.p1 = p1;
    seg.free_face = free_face;
    seg.normal = float_face->normal;
    seg.midpoint = (p0->X + p1->X) * 0.5;
    seg.length = len;
    const double nz = seg.normal[2];
    seg.nh = std::sqrt(std::max(0.0, 1.0 - nz * nz));
    segments.push_back(seg);
  }
  return segments;
}

inline std::array<Complex, 6> integrate_linear_pressure_force(const LinearSolution &sol,
                                                              const std::unordered_set<networkFace *> &faces,
                                                              const Tddd &com,
                                                              double rho) {
  const Complex I(0.0, 1.0);
  const Complex coef = I * sol.omega * rho;

  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});

  // Degree-2 exact quadrature for linear fields on a triangle.
  constexpr double w = 1.0 / 3.0;
  constexpr std::array<std::array<double, 3>, 3> bary = {{
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
  }};

  for (auto *f : faces) {
    const auto it = sol.face_field.find(f);
    if (it == sol.face_field.end())
      continue;
    const auto &field = it->second;
    auto [p0, p1, p2] = f->getPoints();
    const auto &x0 = p0->X;
    const auto &x1 = p1->X;
    const auto &x2 = p2->X;
    const auto &n = f->normal;

    const Complex p_hat0 = coef * field.phi[0];
    const Complex p_hat1 = coef * field.phi[1];
    const Complex p_hat2 = coef * field.phi[2];

    for (const auto &l : bary) {
      const double l0 = l[0], l1 = l[1], l2 = l[2];
      const Complex p_hat = l0 * p_hat0 + l1 * p_hat1 + l2 * p_hat2;
      const Tddd xq = l0 * x0 + l1 * x1 + l2 * x2;
      const std::array<Complex, 3> fq = {p_hat * n[0], p_hat * n[1], p_hat * n[2]};

      out[0] += w * fq[0] * f->area;
      out[1] += w * fq[1] * f->area;
      out[2] += w * fq[2] * f->area;

      const Tddd rq = xq - com;
      const std::array<Complex, 3> tq = {
          rq[1] * fq[2] - rq[2] * fq[1],
          rq[2] * fq[0] - rq[0] * fq[2],
          rq[0] * fq[1] - rq[1] * fq[0],
      };
      out[3] += w * tq[0] * f->area;
      out[4] += w * tq[1] * f->area;
      out[5] += w * tq[2] * f->area;
    }
  }

  return out;
}

inline std::array<Complex, 6> integrate_quadratic_force(const LinearSolution &m,
                                                        const LinearSolution &n,
                                                        const std::unordered_set<networkFace *> &faces,
                                                        const Tddd &com,
                                                        double rho,
                                                        bool use_conjugate) {
  std::array<Complex, 6> out{};
  out.fill(Complex{0.0, 0.0});

  // Quadratic velocity term (Bernoulli). phi^(2) and waterline terms are excluded here.
  // With e^{-i wt}, the spectral coefficient includes a 1/4 factor from time averaging.
  const Complex coef = Complex{-0.25 * rho, 0.0};

  for (auto *f : faces) {
    const auto &gm = m.grad_phi.count(f) ? m.grad_phi.at(f) : face_gradient(f, m.face_field.at(f));
    const auto &gn = n.grad_phi.count(f) ? n.grad_phi.at(f) : face_gradient(f, n.face_field.at(f));

    Complex dot = gm[0] * (use_conjugate ? std::conj(gn[0]) : gn[0]) +
                  gm[1] * (use_conjugate ? std::conj(gn[1]) : gn[1]) +
                  gm[2] * (use_conjugate ? std::conj(gn[2]) : gn[2]);

    const Complex p = coef * dot;
    const std::array<Complex, 3> force = {
        p * f->normal[0] * f->area,
        p * f->normal[1] * f->area,
        p * f->normal[2] * f->area,
    };

    out[0] += force[0];
    out[1] += force[1];
    out[2] += force[2];

    const Tddd r = f->centroid - com;
    const std::array<Complex, 3> torque = {
        r[1] * force[2] - r[2] * force[1],
        r[2] * force[0] - r[0] * force[2],
        r[0] * force[1] - r[1] * force[0],
    };

    out[3] += torque[0];
    out[4] += torque[1];
    out[5] += torque[2];
  }

  return out;
}

inline std::vector<std::vector<Complex>> compute_waterline_eta(
    const std::vector<double> &omegas,
    const std::vector<Solution> &scat_solutions,
    const WaterWaveTheory &wave_base,
    const std::vector<WaterlineSegment> &segments,
    double gravity) {
  const std::size_t nfreq = omegas.size();
  std::vector<std::vector<Complex>> eta(nfreq, std::vector<Complex>(segments.size(), Complex{0.0, 0.0}));
  const Complex I(0.0, 1.0);

  for (std::size_t i = 0; i < nfreq; ++i) {
    const double omega = omegas[i];
    if (!(omega > 0.0))
      continue;
    auto wave = make_wave_for_omega(wave_base, omega);
    wave.A = wave_base.A;
    wave.theta = wave_base.theta;
    wave.phase_shift = wave_base.phase_shift;
    wave.bottom_z = wave_base.bottom_z;

    const auto idx = build_solution_index(scat_solutions[i]);
    for (std::size_t s = 0; s < segments.size(); ++s) {
      const auto &seg = segments[s];
      const Complex phi_inc_0 = incident_phi_hat(wave, seg.p0->X);
      const Complex phi_inc_1 = incident_phi_hat(wave, seg.p1->X);
      const Complex phi_scat_0 = lookup_phi(scat_solutions[i], idx, seg.p0, seg.free_face);
      const Complex phi_scat_1 = lookup_phi(scat_solutions[i], idx, seg.p1, seg.free_face);
      const Complex phi_tot = (phi_inc_0 + phi_inc_1 + phi_scat_0 + phi_scat_1) * 0.5;
      eta[i][s] = (I * omega / gravity) * phi_tot;
    }
  }
  return eta;
}

inline QTFResult compute_waterline_qtf(const std::vector<double> &omegas,
                                       const std::vector<Solution> &scat_solutions,
                                       const WaterWaveTheory &wave_base,
                                       const std::vector<WaterlineSegment> &segments,
                                       const Tddd &com,
                                       double rho,
                                       double gravity,
                                       bool use_symmetry) {
  QTFResult result;
  const std::size_t n = omegas.size();
  result.omegas = omegas;
  const std::array<Complex, 6> zero{};
  result.qminus.assign(n, std::vector<std::array<Complex, 6>>(n, zero));
  result.qplus.assign(n, std::vector<std::array<Complex, 6>>(n, zero));

  if (segments.empty() || n == 0)
    return result;

  const auto eta = compute_waterline_eta(omegas, scat_solutions, wave_base, segments, gravity);
  const Complex coef = Complex{0.25 * rho * gravity, 0.0};

  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t j0 = use_symmetry ? i : 0;
    for (std::size_t j = j0; j < n; ++j) {
      std::array<Complex, 6> fminus{};
      std::array<Complex, 6> fplus{};
      fminus.fill(Complex{0.0, 0.0});
      fplus.fill(Complex{0.0, 0.0});

      for (std::size_t s = 0; s < segments.size(); ++s) {
        const auto &seg = segments[s];
        const Complex eta_i = eta[i][s];
        const Complex eta_j = eta[j][s];
        const Complex w_minus = coef * eta_i * std::conj(eta_j) * seg.nh * seg.length;
        const Complex w_plus = coef * eta_i * eta_j * seg.nh * seg.length;

        const std::array<Complex, 3> df_minus = {w_minus * seg.normal[0], w_minus * seg.normal[1], w_minus * seg.normal[2]};
        const std::array<Complex, 3> df_plus = {w_plus * seg.normal[0], w_plus * seg.normal[1], w_plus * seg.normal[2]};

        const Tddd r = seg.midpoint - com;
        const std::array<Complex, 3> dm_minus = {
            r[1] * df_minus[2] - r[2] * df_minus[1],
            r[2] * df_minus[0] - r[0] * df_minus[2],
            r[0] * df_minus[1] - r[1] * df_minus[0],
        };
        const std::array<Complex, 3> dm_plus = {
            r[1] * df_plus[2] - r[2] * df_plus[1],
            r[2] * df_plus[0] - r[0] * df_plus[2],
            r[0] * df_plus[1] - r[1] * df_plus[0],
        };

        fminus[0] += df_minus[0];
        fminus[1] += df_minus[1];
        fminus[2] += df_minus[2];
        fminus[3] += dm_minus[0];
        fminus[4] += dm_minus[1];
        fminus[5] += dm_minus[2];

        fplus[0] += df_plus[0];
        fplus[1] += df_plus[1];
        fplus[2] += df_plus[2];
        fplus[3] += dm_plus[0];
        fplus[4] += dm_plus[1];
        fplus[5] += dm_plus[2];
      }

      result.qminus[i][j] = fminus;
      result.qplus[i][j] = fplus;
      if (use_symmetry && i != j) {
        result.qminus[j][i] = conjugate_qtf(fminus);
        result.qplus[j][i] = fplus;
      }
    }
  }
  return result;
}

inline QTFResult compute_qtf_result(const std::vector<LinearSolution> &solutions,
                                    const std::unordered_set<networkFace *> &faces,
                                    const Tddd &com,
                                    double rho,
                                    bool use_symmetry) {
  QTFResult result;
  const std::size_t n = solutions.size();
  result.omegas.reserve(n);
  for (const auto &s : solutions)
    result.omegas.push_back(s.omega);
  result.qminus.assign(n, std::vector<std::array<Complex, 6>>(n));
  result.qplus.assign(n, std::vector<std::array<Complex, 6>>(n));

  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t j0 = use_symmetry ? i : 0;
    for (std::size_t j = j0; j < n; ++j) {
      const auto qminus = integrate_quadratic_force(solutions[i], solutions[j], faces, com, rho, true);
      const auto qplus = integrate_quadratic_force(solutions[i], solutions[j], faces, com, rho, false);
      result.qminus[i][j] = qminus;
      result.qplus[i][j] = qplus;
      if (use_symmetry && i != j) {
        result.qminus[j][i] = conjugate_qtf(qminus);
        result.qplus[j][i] = qplus;
      }
    }
  }
  return result;
}

inline QTFReport qtf_symmetry_report(const QTFResult &result) {
  QTFReport report;
  const std::size_t n = result.omegas.size();
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < 6; ++k) {
        const double abs_minus = std::abs(result.qminus[i][j][k]);
        const double abs_plus = std::abs(result.qplus[i][j][k]);
        report.max_abs_qminus = std::max(report.max_abs_qminus, abs_minus);
        report.max_abs_qplus = std::max(report.max_abs_qplus, abs_plus);
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < 6; ++k) {
        const Complex dm = result.qminus[i][j][k] - std::conj(result.qminus[j][i][k]);
        const Complex dp = result.qplus[i][j][k] - result.qplus[j][i][k];
        report.max_sym_abs_qminus = std::max(report.max_sym_abs_qminus, std::abs(dm));
        report.max_sym_abs_qplus = std::max(report.max_sym_abs_qplus, std::abs(dp));
      }
    }
  }
  if (report.max_abs_qminus > 0.0)
    report.max_sym_rel_qminus = report.max_sym_abs_qminus / report.max_abs_qminus;
  if (report.max_abs_qplus > 0.0)
    report.max_sym_rel_qplus = report.max_sym_abs_qplus / report.max_abs_qplus;
  return report;
}

inline NewmanReport qtf_newman_report(const QTFResult &result) {
  NewmanReport report;
  const std::size_t n = result.omegas.size();
  double max_abs_full = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < 6; ++k) {
        max_abs_full = std::max(max_abs_full, std::abs(result.qminus[i][j][k]));
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const auto &qii = result.qminus[i][i];
      const auto &qjj = result.qminus[j][j];
      for (std::size_t k = 0; k < 6; ++k) {
        const Complex qn = 0.5 * (qii[k] + qjj[k]);
        report.max_abs = std::max(report.max_abs, std::abs(result.qminus[i][j][k] - qn));
      }
    }
  }
  if (max_abs_full > 0.0)
    report.max_rel = report.max_abs / max_abs_full;
  return report;
}

inline std::vector<QTFPair> compute_qtf_pairs(const std::vector<LinearSolution> &solutions,
                                              const std::unordered_set<networkFace *> &faces,
                                              const Tddd &com,
                                              double rho,
                                              bool use_symmetry) {
  const auto result = compute_qtf_result(solutions, faces, com, rho, use_symmetry);
  std::vector<QTFPair> pairs;
  const std::size_t n = result.omegas.size();
  pairs.reserve(use_symmetry ? (n * (n + 1)) / 2 : n * n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t j0 = use_symmetry ? i : 0;
    for (std::size_t j = j0; j < n; ++j) {
      QTFPair pair;
      pair.omega_m = result.omegas[i];
      pair.omega_n = result.omegas[j];
      pair.qminus = result.qminus[i][j];
      pair.qplus = result.qplus[i][j];
      pairs.emplace_back(std::move(pair));
    }
  }
  return pairs;
}

inline void write_qtf_csv(const std::filesystem::path &path, const QTFResult &result, bool plus, bool upper_only) {
  std::ofstream ofs(path);
  ofs << std::scientific << std::setprecision(12);
  ofs << "omega_m,omega_n,"
      << "fx_re,fx_im,fy_re,fy_im,fz_re,fz_im,"
      << "mx_re,mx_im,my_re,my_im,mz_re,mz_im\n";

  const std::size_t n = result.omegas.size();
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t j0 = upper_only ? i : 0;
    for (std::size_t j = j0; j < n; ++j) {
      const auto &q = plus ? result.qplus[i][j] : result.qminus[i][j];
      ofs << result.omegas[i] << "," << result.omegas[j];
      for (const auto &v : q)
        ofs << "," << v.real() << "," << v.imag();
      ofs << "\n";
    }
  }
}

inline void write_qtf_csv(const std::filesystem::path &path, const std::vector<QTFPair> &pairs, bool plus) {
  std::ofstream ofs(path);
  ofs << std::scientific << std::setprecision(12);
  ofs << "omega_m,omega_n,"
      << "fx_re,fx_im,fy_re,fy_im,fz_re,fz_im,"
      << "mx_re,mx_im,my_re,my_im,mz_re,mz_im\n";

  for (const auto &pair : pairs) {
    const auto &q = plus ? pair.qplus : pair.qminus;
    ofs << pair.omega_m << "," << pair.omega_n;
    for (const auto &v : q)
      ofs << "," << v.real() << "," << v.imag();
    ofs << "\n";
  }
}

inline void write_qtf_newman_csv(const std::filesystem::path &path, const QTFResult &result, bool upper_only) {
  std::ofstream ofs(path);
  ofs << std::scientific << std::setprecision(12);
  ofs << "omega_m,omega_n,"
      << "fx_re,fx_im,fy_re,fy_im,fz_re,fz_im,"
      << "mx_re,mx_im,my_re,my_im,mz_re,mz_im\n";

  const std::size_t n = result.omegas.size();
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t j0 = upper_only ? i : 0;
    for (std::size_t j = j0; j < n; ++j) {
      ofs << result.omegas[i] << "," << result.omegas[j];
      const auto &qii = result.qminus[i][i];
      const auto &qjj = result.qminus[j][j];
      for (std::size_t k = 0; k < 6; ++k) {
        const Complex qn = 0.5 * (qii[k] + qjj[k]);
        ofs << "," << qn.real() << "," << qn.imag();
      }
      ofs << "\n";
    }
  }
}

} // namespace bem_frequency_domain
