#ifndef SET_DIRECT_INTEGRATION_SIMD_DOUBLE_HPP
#define SET_DIRECT_INTEGRATION_SIMD_DOUBLE_HPP

#include <arm_neon.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "dunavant_rules.hpp"

/*
 * setDirectIntegration の SIMD最適化版（Linear要素専用, double精度）
 *
 * float32x4_t (4レーン) 版と同じアルゴリズムだが、
 * float64x2_t (2レーン) を使用して完全なdouble精度を維持する。
 *
 * 同一セル（バケット）内のターゲットは同じ近傍ソースリストを共有する。
 * この性質を利用して：
 *   1. ターゲットをバケットごとにグループ化
 *   2. ソースフェイスデータを事前にフラット化（LinearFaceDataD）
 *   3. NEON float64x2_t で2ターゲット同時処理（非隣接のみ）
 *   4. 隣接ペアは既存スカラーコールバックで処理
 *   5. 非線形ソースも既存スカラーコールバックで処理
 */

// ---------------------------------------------------------------------------
// ソースフェイスデータ（Linear要素用, double精度）
// ---------------------------------------------------------------------------
struct LinearFaceDataD {
  double X0[3], X1[3], X2[3];      // 頂点座標
  double J_det;                    // ヤコビアン = 2 * area
  double cross[3];                 // 2 * area * normal (= Cross(X1-X0, X2-X0))
  double Xc[3];                    // 重心座標
  int32_t idx[3];                  // グローバルDOFインデックス
  const void* face_vertex_ptrs[3]; // 隣接判定用ポインタ (target4FMM*)
};

// ---------------------------------------------------------------------------
// 事前計算ガウス点テーブル（全フェイス共通の定数, double精度）
// ---------------------------------------------------------------------------
namespace GaussTableSIMD_D {

constexpr std::array<double, 3> modTriShape3(double t0, double t1) {
  return {t0, t1 * (1.0 - t0), (t0 - 1.0) * (t1 - 1.0)};
}

// 25点ガウス求積点 (GW5xGW5)
struct GW5x5TableD {
  double N[25][3];
  double w_eff[25];

  constexpr GW5x5TableD() : N{}, w_eff{} {
    constexpr double pts[25][3] = {
        {0.046910077030668, 0.046910077030668, 0.0140335872156063},
        {0.046910077030668, 0.230765344947158, 0.0283499999999991},
        {0.046910077030668, 0.5, 0.0336962680968792},
        {0.046910077030668, 0.769234655052841, 0.0283499999999991},
        {0.046910077030668, 0.953089922969332, 0.0140335872156063},
        {0.230765344947158, 0.046910077030668, 0.0283499999999991},
        {0.230765344947158, 0.230765344947158, 0.0572713510559978},
        {0.230765344947158, 0.5, 0.0680716331376877},
        {0.230765344947158, 0.769234655052841, 0.0572713510559978},
        {0.230765344947158, 0.953089922969332, 0.0283499999999991},
        {0.5, 0.046910077030668, 0.0336962680968792},
        {0.5, 0.230765344947158, 0.0680716331376877},
        {0.5, 0.5, 0.0809086419753086},
        {0.5, 0.769234655052841, 0.0680716331376877},
        {0.5, 0.953089922969332, 0.0336962680968792},
        {0.769234655052841, 0.046910077030668, 0.0283499999999991},
        {0.769234655052841, 0.230765344947158, 0.0572713510559978},
        {0.769234655052841, 0.5, 0.0680716331376877},
        {0.769234655052841, 0.769234655052841, 0.0572713510559978},
        {0.769234655052841, 0.953089922969332, 0.0283499999999991},
        {0.953089922969332, 0.046910077030668, 0.0140335872156063},
        {0.953089922969332, 0.230765344947158, 0.0283499999999991},
        {0.953089922969332, 0.5, 0.0336962680968792},
        {0.953089922969332, 0.769234655052841, 0.0283499999999991},
        {0.953089922969332, 0.953089922969332, 0.0140335872156063}};

    for (int i = 0; i < 25; ++i) {
      auto ns = modTriShape3(pts[i][0], pts[i][1]);
      N[i][0] = ns[0];
      N[i][1] = ns[1];
      N[i][2] = ns[2];
      w_eff[i] = pts[i][2] * (1.0 - pts[i][0]);
    }
  }
};

// 1点（重心）ガウス求積 — far-field用
struct GW1x1TableD {
  double N[1][3];
  double w_eff[1];

  constexpr GW1x1TableD() : N{}, w_eff{} {
    N[0][0] = 1.0 / 3.0;
    N[0][1] = 1.0 / 3.0;
    N[0][2] = 1.0 / 3.0;
    w_eff[0] = 0.5;
  }
};

// 6点 Dunavant 求積（far non-adjacent 用, degree 4）
struct Dunavant6TableD {
  double N[6][3];
  double w_eff[6];
  constexpr Dunavant6TableD() : N{}, w_eff{} {
    constexpr double pts[6][4] = {
        {0.816847572980459, 0.091576213509771, 0.091576213509771, 0.054975871827661},
        {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.054975871827661},
        {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.054975871827661},
        {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.111690794839006},
        {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.111690794839006},
        {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.111690794839006}};
    for (int i = 0; i < 6; ++i) {
      N[i][0] = pts[i][0];
      N[i][1] = pts[i][1];
      N[i][2] = pts[i][2];
      w_eff[i] = pts[i][3];
    }
  }
};

// 12点 Dunavant 求積（near non-adjacent 用, degree 6）
struct Dunavant12TableD {
  double N[12][3];
  double w_eff[12];
  constexpr Dunavant12TableD() : N{}, w_eff{} {
    constexpr double pts[12][4] = {
        {0.501426509658179, 0.249286745170910, 0.249286745170910, 0.058393137863190},
        {0.249286745170910, 0.501426509658179, 0.249286745170910, 0.058393137863190},
        {0.249286745170910, 0.249286745170910, 0.501426509658179, 0.058393137863190},
        {0.873821971016996, 0.063089014491502, 0.063089014491502, 0.025422453185104},
        {0.063089014491502, 0.873821971016996, 0.063089014491502, 0.025422453185104},
        {0.063089014491502, 0.063089014491502, 0.873821971016996, 0.025422453185104},
        {0.053145049844817, 0.310352451033784, 0.636502499121399, 0.041425537809187},
        {0.053145049844817, 0.636502499121399, 0.310352451033784, 0.041425537809187},
        {0.310352451033784, 0.053145049844817, 0.636502499121399, 0.041425537809187},
        {0.310352451033784, 0.636502499121399, 0.053145049844817, 0.041425537809187},
        {0.636502499121399, 0.053145049844817, 0.310352451033784, 0.041425537809187},
        {0.636502499121399, 0.310352451033784, 0.053145049844817, 0.041425537809187}};
    for (int i = 0; i < 12; ++i) {
      N[i][0] = pts[i][0];
      N[i][1] = pts[i][1];
      N[i][2] = pts[i][2];
      w_eff[i] = pts[i][3];
    }
  }
};

inline constexpr GW5x5TableD gw5x5{};
inline constexpr GW1x1TableD gw1x1{};
inline constexpr Dunavant6TableD dun6{};
inline constexpr Dunavant12TableD dun12{};

} // namespace GaussTableSIMD_D

// ---------------------------------------------------------------------------
// NEON double ヘルパー
// ---------------------------------------------------------------------------
namespace neon_helpers_d {

// Newton-Raphson refined reciprocal square root (3 iterations for double)
// Input x, returns 1/sqrt(x) with ~52-bit precision
inline float64x2_t rsqrt_nr(float64x2_t x) {
  float64x2_t est = vrsqrteq_f64(x);
  est = vmulq_f64(est, vrsqrtsq_f64(vmulq_f64(x, est), est));
  est = vmulq_f64(est, vrsqrtsq_f64(vmulq_f64(x, est), est));
  est = vmulq_f64(est, vrsqrtsq_f64(vmulq_f64(x, est), est));
  return est;
}

} // namespace neon_helpers_d

// ---------------------------------------------------------------------------
// SIMD double版 setDirectIntegration（Linear要素専用）
//
// float SIMD版と同一アルゴリズム。float64x2_t（2レーン）で処理。
// 出力: 各ターゲットの near_indices, near_weights_phi, near_weights_phin,
//       near_run_base_idx, near_run_pos, near_run_len を設定する。
// ---------------------------------------------------------------------------
template <typename BucketType, typename TargetPtr>
void setDirectIntegrationSIMD_double_linear(
    BucketType& B_poles,
    const std::vector<TargetPtr>& targets,
    int near_quadrature_points = 6) {
  using target_t = std::remove_pointer_t<decltype(&*targets[0])>;
  using source_t = typename decltype(B_poles.data1D)::value_type::element_type;

  const auto& dun6 = GaussTableSIMD_D::dun6;
  // near non-adjacent: Dunavant(near_quadrature_points), far: Dunavant 6
  const auto& near_rule = getDunavantP2MRule(near_quadrature_points);
  const int num_near_gp = static_cast<int>(near_rule.size());
  constexpr int MAX_DUN = 25;
  double near_N_buf[MAX_DUN][3];
  double near_w_buf[MAX_DUN];
  for (int i = 0; i < num_near_gp; ++i) {
    near_N_buf[i][0] = near_rule[i][0];
    near_N_buf[i][1] = near_rule[i][1];
    near_N_buf[i][2] = near_rule[i][2];
    near_w_buf[i] = near_rule[i][3];
  }

  // =======================================================================
  // Phase 1: ターゲットをバケットごとにグループ化
  // =======================================================================
  std::unordered_map<const void*, size_t> bucket_to_group;
  struct BucketGroup {
    std::vector<target_t*> tgts;
  };
  std::vector<BucketGroup> groups;

  for (auto& t : targets) {
    auto* b = B_poles.getBucketAtDeepest(t->Xtarget);
    auto it = bucket_to_group.find(static_cast<const void*>(b));
    if (it == bucket_to_group.end()) {
      bucket_to_group[static_cast<const void*>(b)] = groups.size();
      groups.push_back({{&*t}});
    } else {
      groups[it->second].tgts.push_back(&*t);
    }
  }

  // =======================================================================
  // Phase 2+3: バケットごとに近傍ソース収集 → SIMD + スカラー処理
  // =======================================================================
#pragma omp parallel
  {
    // Thread-local accumulators (2 per thread for 2-lane SIMD)
    DirectAccumulator accs[2];

#pragma omp for schedule(dynamic, 1)
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      auto& tgts = groups[gi].tgts;
      const size_t num_targets = tgts.size();

      auto* b_deepest = B_poles.getBucketAtDeepest(tgts[0]->Xtarget);

      // ---- 近傍ソース収集 ----
      std::vector<LinearFaceDataD> linear_faces;
      struct SourceRef {
        source_t* ptr;
        int32_t linear_idx;
      };
      std::vector<SourceRef> all_sources;
      double near_region_sq = 0.;
      size_t cell_count = 0;

      auto collectFromBucket = [&](const auto* B) {
        ++cell_count;
        for (const auto& source : B->data1D_vector) {
          if (!source->fill_direct_entries)
            continue;

          if (source->dof_indices[0] >= 0) {
            const auto& fv = source->face_vertices;
            LinearFaceDataD fd;
            for (int k = 0; k < 3; ++k) {
              fd.X0[k] = fv[0]->Xtarget[k];
              fd.X1[k] = fv[1]->Xtarget[k];
              fd.X2[k] = fv[2]->Xtarget[k];
            }
            double dx1[3], dx2[3];
            for (int k = 0; k < 3; ++k) {
              dx1[k] = fv[1]->Xtarget[k] - fv[0]->Xtarget[k];
              dx2[k] = fv[2]->Xtarget[k] - fv[0]->Xtarget[k];
            }
            double cx = dx1[1] * dx2[2] - dx1[2] * dx2[1];
            double cy = dx1[2] * dx2[0] - dx1[0] * dx2[2];
            double cz = dx1[0] * dx2[1] - dx1[1] * dx2[0];
            fd.J_det = std::sqrt(cx * cx + cy * cy + cz * cz);
            fd.cross[0] = cx;
            fd.cross[1] = cy;
            fd.cross[2] = cz;
            fd.Xc[0] = (fd.X0[0] + fd.X1[0] + fd.X2[0]) / 3.0;
            fd.Xc[1] = (fd.X0[1] + fd.X1[1] + fd.X2[1]) / 3.0;
            fd.Xc[2] = (fd.X0[2] + fd.X1[2] + fd.X2[2]) / 3.0;
            fd.idx[0] = source->dof_indices[0];
            fd.idx[1] = source->dof_indices[1];
            fd.idx[2] = source->dof_indices[2];
            fd.face_vertex_ptrs[0] = static_cast<const void*>(fv[0]);
            fd.face_vertex_ptrs[1] = static_cast<const void*>(fv[1]);
            fd.face_vertex_ptrs[2] = static_cast<const void*>(fv[2]);
            int32_t fidx = static_cast<int32_t>(linear_faces.size());
            linear_faces.push_back(fd);
            all_sources.push_back(SourceRef{source.get(), fidx});
            near_region_sq = source->near_region * source->near_region;
          } else {
            all_sources.push_back(SourceRef{source.get(), -1});
          }
        }
      };

      collectFromBucket(b_deepest);
      for (const auto& b : b_deepest->buckets_near)
        collectFromBucket(b);
      {
        auto* bp = b_deepest->parent;
        while (bp != nullptr) {
          for (auto* B : bp->buckets_near)
            if (!B->hasChildren())
              collectFromBucket(B);
          bp = bp->parent;
        }
      }

      if (const char* env = std::getenv("BEM_NEAR_DUP_DEBUG"); env && std::string(env) != "0") {
        std::unordered_set<const source_t*> seen_sources;
        seen_sources.reserve(all_sources.size() * 2 + 1);
        std::size_t dup_count = 0;
        for (const auto& s : all_sources) {
          if (!seen_sources.emplace(s.ptr).second)
            ++dup_count;
        }
        if (dup_count > 0) {
          static std::atomic<int> dup_report_count{0};
          const int report_id = dup_report_count.fetch_add(1, std::memory_order_relaxed);
          if (report_id < 20) {
            std::cout << Red << "[nearfield:dup] " << Cyan
                      << "dup_sources=" << Yellow << dup_count
                      << Cyan << " total_refs=" << Yellow << all_sources.size()
                      << Cyan << " unique_sources=" << Yellow << seen_sources.size()
                      << Cyan << " near_cells=" << Yellow << cell_count
                      << colorReset << std::endl;
          }
        }
      }

      const size_t num_linear_faces = linear_faces.size();

      // ---- ソースが無い場合：空の結果を設定 ----
      if (all_sources.empty()) {
        for (auto* tgt : tgts) {
          tgt->near_indices.clear();
          tgt->near_weights_phi.clear();
          tgt->near_weights_phin.clear();
          tgt->near_run_base_idx.clear();
          tgt->near_run_pos.clear();
          tgt->near_run_len.clear();
          tgt->near_cell_count = cell_count;
        }
        continue;
      }

      // ---- ターゲットを2個ずつ処理 (float64x2_t = 2レーン) ----
      for (size_t t_base = 0; t_base < num_targets; t_base += 2) {
        const size_t t_end = std::min(t_base + 2, num_targets);
        const size_t t_count = t_end - t_base;

        for (size_t i = 0; i < t_count; ++i)
          accs[i].reset();

        // ターゲット座標をSIMDレジスタにロード (2レーン)
        double tx[2], ty[2], tz[2];
        for (size_t i = 0; i < t_count; ++i) {
          tx[i] = tgts[t_base + i]->Xtarget[0];
          ty[i] = tgts[t_base + i]->Xtarget[1];
          tz[i] = tgts[t_base + i]->Xtarget[2];
        }
        // 端数レーンをパディング
        for (size_t i = t_count; i < 2; ++i) {
          tx[i] = tx[t_count - 1];
          ty[i] = ty[t_count - 1];
          tz[i] = tz[t_count - 1];
        }
        float64x2_t Xt_x = vld1q_f64(tx);
        float64x2_t Xt_y = vld1q_f64(ty);
        float64x2_t Xt_z = vld1q_f64(tz);

        // 隣接判定用のターゲットポインタ
        const void* tgt_ptrs[2];
        for (size_t i = 0; i < t_count; ++i)
          tgt_ptrs[i] = static_cast<const void*>(tgts[t_base + i]);
        for (size_t i = t_count; i < 2; ++i)
          tgt_ptrs[i] = tgt_ptrs[t_count - 1];

        // ==============================================================
        // SIMD pass: 非隣接 Linear フェイス (float64x2_t, 2レーン)
        // ==============================================================
        for (size_t fi = 0; fi < num_linear_faces; ++fi) {
          const auto& face = linear_faces[fi];

          // 隣接判定（各レーン独立）
          bool is_adj[2];
          for (int l = 0; l < 2; ++l) {
            is_adj[l] = (tgt_ptrs[l] == face.face_vertex_ptrs[0] ||
                         tgt_ptrs[l] == face.face_vertex_ptrs[1] ||
                         tgt_ptrs[l] == face.face_vertex_ptrs[2]);
          }

          // near/far判定: 重心とターゲットの距離²
          float64x2_t dx_c = vsubq_f64(vdupq_n_f64(face.Xc[0]), Xt_x);
          float64x2_t dy_c = vsubq_f64(vdupq_n_f64(face.Xc[1]), Xt_y);
          float64x2_t dz_c = vsubq_f64(vdupq_n_f64(face.Xc[2]), Xt_z);
          float64x2_t dist_sq = vfmaq_f64(vfmaq_f64(vmulq_f64(dx_c, dx_c), dy_c, dy_c), dz_c, dz_c);
          uint64x2_t near_mask = vcltq_f64(dist_sq, vdupq_n_f64(near_region_sq));

          uint64_t mask_bits[2];
          vst1q_u64(mask_bits, near_mask);
          const bool all_far = (!mask_bits[0] && !mask_bits[1]);

          // Dunavant テーブル選択: far→D6, near→D(p2m設定)
          int num_gp;
          const double (*gp_N)[3];
          const double* gp_w;
          if (all_far) {
            num_gp = 6;
            gp_N = dun6.N;
            gp_w = dun6.w_eff;
          } else {
            num_gp = num_near_gp;
            gp_N = near_N_buf;
            gp_w = near_w_buf;
          }

          // ソースフェイス定数（broadcast）
          float64x2_t f_Jdet = vdupq_n_f64(face.J_det);
          float64x2_t f_cx = vdupq_n_f64(face.cross[0]);
          float64x2_t f_cy = vdupq_n_f64(face.cross[1]);
          float64x2_t f_cz = vdupq_n_f64(face.cross[2]);

          // Per-face, per-lane, per-DOF accumulator
          float64x2_t WGN0 = vdupq_n_f64(0.), WGN1 = vdupq_n_f64(0.), WGN2 = vdupq_n_f64(0.);
          float64x2_t WGnN0 = vdupq_n_f64(0.), WGnN1 = vdupq_n_f64(0.), WGnN2 = vdupq_n_f64(0.);

          for (int gp = 0; gp < num_gp; ++gp) {
            const double N0 = gp_N[gp][0];
            const double N1v = gp_N[gp][1];
            const double N2v = gp_N[gp][2];
            const double weff = gp_w[gp];

            // ソースガウス点位置（全レーン共通, broadcast）
            const double Xgp_x = N0 * face.X0[0] + N1v * face.X1[0] + N2v * face.X2[0];
            const double Xgp_y = N0 * face.X0[1] + N1v * face.X1[1] + N2v * face.X2[1];
            const double Xgp_z = N0 * face.X0[2] + N1v * face.X1[2] + N2v * face.X2[2];

            // R = Xgp - Xtarget (2 lanes SIMD)
            float64x2_t Rx = vsubq_f64(vdupq_n_f64(Xgp_x), Xt_x);
            float64x2_t Ry = vsubq_f64(vdupq_n_f64(Xgp_y), Xt_y);
            float64x2_t Rz = vsubq_f64(vdupq_n_f64(Xgp_z), Xt_z);

            // |R|² and 1/|R| (Newton-Raphson 3 iterations for double)
            float64x2_t nr2 = vfmaq_f64(vfmaq_f64(vmulq_f64(Rx, Rx), Ry, Ry), Rz, Rz);
            float64x2_t nr_inv = neon_helpers_d::rsqrt_nr(nr2);

            // w_eff / |R|
            float64x2_t w_nr = vmulq_f64(vdupq_n_f64(weff), nr_inv);

            // dot(R, cross) = Rx*cx + Ry*cy + Rz*cz
            float64x2_t dot_R_c = vfmaq_f64(vfmaq_f64(vmulq_f64(Rx, f_cx), Ry, f_cy), Rz, f_cz);

            // WG: J_det * w_eff * nr_inv
            float64x2_t Jdet_w_nr = vmulq_f64(f_Jdet, w_nr);

            // WGn: -dot(R,cross) * w_eff * nr_inv³
            float64x2_t nr_inv2 = vmulq_f64(nr_inv, nr_inv);
            float64x2_t neg_dot_w_nr3 = vnegq_f64(vmulq_f64(vmulq_f64(dot_R_c, w_nr), nr_inv2));

            // Accumulate per DOF
            WGN0 = vfmaq_f64(WGN0, Jdet_w_nr, vdupq_n_f64(N0));
            WGN1 = vfmaq_f64(WGN1, Jdet_w_nr, vdupq_n_f64(N1v));
            WGN2 = vfmaq_f64(WGN2, Jdet_w_nr, vdupq_n_f64(N2v));

            WGnN0 = vfmaq_f64(WGnN0, neg_dot_w_nr3, vdupq_n_f64(N0));
            WGnN1 = vfmaq_f64(WGnN1, neg_dot_w_nr3, vdupq_n_f64(N1v));
            WGnN2 = vfmaq_f64(WGnN2, neg_dot_w_nr3, vdupq_n_f64(N2v));
          } // gp loop

          // Scatter 結果を各レーンの accumulator に書き出し
          double wgn0[2], wgn1[2], wgn2[2];
          double wgnn0[2], wgnn1[2], wgnn2[2];
          vst1q_f64(wgn0, WGN0);
          vst1q_f64(wgn1, WGN1);
          vst1q_f64(wgn2, WGN2);
          vst1q_f64(wgnn0, WGnN0);
          vst1q_f64(wgnn1, WGnN1);
          vst1q_f64(wgnn2, WGnN2);

          for (size_t lane = 0; lane < t_count; ++lane) {
            if (is_adj[lane])
              continue; // 隣接はスカラーで別途処理
            accs[lane].add(face.idx[0], wgn0[lane], wgnn0[lane]);
            accs[lane].add(face.idx[1], wgn1[lane], wgnn1[lane]);
            accs[lane].add(face.idx[2], wgn2[lane], wgnn2[lane]);
          }
        } // fi (linear faces) loop

        // ==============================================================
        // Scalar pass: 隣接 Linear ソース + 全 Non-linear ソース
        // ==============================================================
        for (size_t lane = 0; lane < t_count; ++lane) {
          auto* tgt = tgts[t_base + lane];
          for (const auto& ns : all_sources) {
            if (ns.linear_idx >= 0) {
              if (ns.ptr->isAdjacentTo(tgt)) {
                ns.ptr->fill_direct_entries(tgt, accs[lane]);
              }
            } else {
              if (ns.ptr->fill_direct_entries_nonadj && !ns.ptr->isAdjacentTo(tgt)) {
                ns.ptr->fill_direct_entries_nonadj(tgt, accs[lane]);
              } else {
                ns.ptr->fill_direct_entries(tgt, accs[lane]);
              }
            }
          }
        }

        // ==============================================================
        // Phase 4: 結果抽出 (near_indices / weights / RLE)
        // ==============================================================
        for (size_t lane = 0; lane < t_count; ++lane) {
          auto* tgt = tgts[t_base + lane];
          auto& acc = accs[lane];

          std::sort(acc.touched_list.begin(), acc.touched_list.end());

          tgt->near_indices.clear();
          tgt->near_weights_phi.clear();
          tgt->near_weights_phin.clear();

          const size_t n_touched = acc.touched_list.size();
          tgt->near_indices.reserve(n_touched);
          tgt->near_weights_phi.reserve(n_touched);
          tgt->near_weights_phin.reserve(n_touched);

          for (int32_t idx : acc.touched_list) {
            tgt->near_indices.push_back(idx);
            tgt->near_weights_phi.push_back(acc.phi_acc[idx]);
            tgt->near_weights_phin.push_back(acc.phin_acc[idx]);
            acc.phi_acc[idx] = 0.;
            acc.phin_acc[idx] = 0.;
            acc.touched_flags[idx] = 0;
          }
          acc.update_high_water_mark();

          // Build RLE
          const size_t n = tgt->near_indices.size();
          tgt->near_run_base_idx.clear();
          tgt->near_run_pos.clear();
          tgt->near_run_len.clear();
          if (n > 0) {
            tgt->near_run_base_idx.reserve(n / 4 + 1);
            tgt->near_run_pos.reserve(n / 4 + 1);
            tgt->near_run_len.reserve(n / 4 + 1);
            size_t pos0 = 0;
            int32_t base = tgt->near_indices[0];
            int32_t prev = base;
            for (size_t i = 1; i < n; ++i) {
              const int32_t cur = tgt->near_indices[i];
              if (cur == prev + 1) {
                prev = cur;
                continue;
              }
              tgt->near_run_base_idx.push_back(base);
              tgt->near_run_pos.push_back(static_cast<int32_t>(pos0));
              tgt->near_run_len.push_back(static_cast<int32_t>(i - pos0));
              pos0 = i;
              base = cur;
              prev = cur;
            }
            tgt->near_run_base_idx.push_back(base);
            tgt->near_run_pos.push_back(static_cast<int32_t>(pos0));
            tgt->near_run_len.push_back(static_cast<int32_t>(n - pos0));
          }

          tgt->near_cell_count = cell_count;
        } // lane loop
      } // t_base loop
    } // gi loop
  } // omp parallel
}

// ---------------------------------------------------------------------------
// セルグループ化版 スカラー setDirectIntegration（Linear要素専用, double精度）
//
// SIMD intrinsics を一切使わず、セルグループ化 + LinearFaceDataD のデータ構造
// 最適化のみの効果を測定するためのベンチマーク用。
// ---------------------------------------------------------------------------
template <typename BucketType, typename TargetPtr>
void setDirectIntegrationCellScalar_linear(
    BucketType& B_poles,
    const std::vector<TargetPtr>& targets,
    int near_quadrature_points = 6) {
  using target_t = std::remove_pointer_t<decltype(&*targets[0])>;
  using source_t = typename decltype(B_poles.data1D)::value_type::element_type;

  const auto& dun6 = GaussTableSIMD_D::dun6;
  const auto& near_rule = getDunavantP2MRule(near_quadrature_points);
  const int num_near_gp = static_cast<int>(near_rule.size());
  constexpr int MAX_DUN = 25;
  double near_N_buf[MAX_DUN][3];
  double near_w_buf[MAX_DUN];
  for (int i = 0; i < num_near_gp; ++i) {
    near_N_buf[i][0] = near_rule[i][0];
    near_N_buf[i][1] = near_rule[i][1];
    near_N_buf[i][2] = near_rule[i][2];
    near_w_buf[i] = near_rule[i][3];
  }

  // Phase 1: ターゲットをバケットごとにグループ化
  std::unordered_map<const void*, size_t> bucket_to_group;
  struct BucketGroup {
    std::vector<target_t*> tgts;
  };
  std::vector<BucketGroup> groups;

  for (auto& t : targets) {
    auto* b = B_poles.getBucketAtDeepest(t->Xtarget);
    auto it = bucket_to_group.find(static_cast<const void*>(b));
    if (it == bucket_to_group.end()) {
      bucket_to_group[static_cast<const void*>(b)] = groups.size();
      groups.push_back({{&*t}});
    } else {
      groups[it->second].tgts.push_back(&*t);
    }
  }

  const size_t total_targets = targets.size();
  std::atomic<size_t> processed_targets{0};
  const size_t report_interval = std::max<size_t>(1, total_targets / 20); // ~5%
  std::atomic<size_t> next_report{report_interval};
  const auto progress_t0 = std::chrono::steady_clock::now();

  std::cout << " [nearfield:cell_scalar] begin, targets=" << total_targets
            << ", report_interval=" << report_interval << std::endl;

#pragma omp parallel
  {
    DirectAccumulator acc;

#pragma omp for schedule(dynamic, 1)
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      auto& tgts = groups[gi].tgts;
      const size_t num_targets = tgts.size();

      auto* b_deepest = B_poles.getBucketAtDeepest(tgts[0]->Xtarget);

      // Phase 2: 近傍ソース収集
      std::vector<LinearFaceDataD> linear_faces;
      struct SourceRef {
        source_t* ptr;
        int32_t linear_idx;
      };
      std::vector<SourceRef> all_sources;
      double near_region_sq = 0.;
      size_t cell_count = 0;

      auto collectFromBucket = [&](const auto* B) {
        ++cell_count;
        for (const auto& source : B->data1D_vector) {
          if (!source->fill_direct_entries)
            continue;
          if (source->dof_indices[0] >= 0) {
            const auto& fv = source->face_vertices;
            LinearFaceDataD fd;
            for (int k = 0; k < 3; ++k) {
              fd.X0[k] = fv[0]->Xtarget[k];
              fd.X1[k] = fv[1]->Xtarget[k];
              fd.X2[k] = fv[2]->Xtarget[k];
            }
            double dx1[3], dx2[3];
            for (int k = 0; k < 3; ++k) {
              dx1[k] = fv[1]->Xtarget[k] - fv[0]->Xtarget[k];
              dx2[k] = fv[2]->Xtarget[k] - fv[0]->Xtarget[k];
            }
            double cx = dx1[1] * dx2[2] - dx1[2] * dx2[1];
            double cy = dx1[2] * dx2[0] - dx1[0] * dx2[2];
            double cz = dx1[0] * dx2[1] - dx1[1] * dx2[0];
            fd.J_det = std::sqrt(cx * cx + cy * cy + cz * cz);
            fd.cross[0] = cx;
            fd.cross[1] = cy;
            fd.cross[2] = cz;
            fd.Xc[0] = (fd.X0[0] + fd.X1[0] + fd.X2[0]) / 3.0;
            fd.Xc[1] = (fd.X0[1] + fd.X1[1] + fd.X2[1]) / 3.0;
            fd.Xc[2] = (fd.X0[2] + fd.X1[2] + fd.X2[2]) / 3.0;
            fd.idx[0] = source->dof_indices[0];
            fd.idx[1] = source->dof_indices[1];
            fd.idx[2] = source->dof_indices[2];
            fd.face_vertex_ptrs[0] = static_cast<const void*>(fv[0]);
            fd.face_vertex_ptrs[1] = static_cast<const void*>(fv[1]);
            fd.face_vertex_ptrs[2] = static_cast<const void*>(fv[2]);
            int32_t fidx = static_cast<int32_t>(linear_faces.size());
            linear_faces.push_back(fd);
            all_sources.push_back(SourceRef{source.get(), fidx});
            near_region_sq = source->near_region * source->near_region;
          } else {
            all_sources.push_back(SourceRef{source.get(), -1});
          }
        }
      };

      collectFromBucket(b_deepest);
      for (const auto& b : b_deepest->buckets_near)
        collectFromBucket(b);
      {
        auto* bp = b_deepest->parent;
        while (bp != nullptr) {
          for (auto* B : bp->buckets_near)
            if (!B->hasChildren())
              collectFromBucket(B);
          bp = bp->parent;
        }
      }

      const size_t num_linear_faces = linear_faces.size();

      if (all_sources.empty()) {
        for (auto* tgt : tgts) {
          tgt->near_indices.clear();
          tgt->near_weights_phi.clear();
          tgt->near_weights_phin.clear();
          tgt->near_run_base_idx.clear();
          tgt->near_run_pos.clear();
          tgt->near_run_len.clear();
          tgt->near_cell_count = cell_count;

          const size_t done = processed_targets.fetch_add(1, std::memory_order_relaxed) + 1;
          size_t expected = next_report.load(std::memory_order_relaxed);
          while (done >= expected && expected <= total_targets) {
            if (next_report.compare_exchange_weak(expected,
                                                  expected + report_interval,
                                                  std::memory_order_relaxed)) {
#pragma omp critical(nearfield_progress_log)
              {
                const auto now = std::chrono::steady_clock::now();
                const double elapsed =
                    std::chrono::duration<double>(now - progress_t0).count();
                const int percent = static_cast<int>((100.0 * done) / static_cast<double>(total_targets));
                std::cout << " [nearfield:cell_scalar] progress " << done << "/"
                          << total_targets << " (" << percent << "%), elapsed="
                          << elapsed << " s" << std::endl;
              }
              break;
            }
          }
        }
        continue;
      }

      // Phase 3: ターゲットを1個ずつ処理（SIMDなし）
      for (size_t ti = 0; ti < num_targets; ++ti) {
        auto* tgt = tgts[ti];
        acc.reset();

        const double Xt_x = tgt->Xtarget[0];
        const double Xt_y = tgt->Xtarget[1];
        const double Xt_z = tgt->Xtarget[2];
        const void* tgt_ptr = static_cast<const void*>(tgt);

        // 非隣接 Linear フェイス：スカラー double で積分
        for (size_t fi = 0; fi < num_linear_faces; ++fi) {
          const auto& face = linear_faces[fi];

          // 隣接判定
          if (tgt_ptr == face.face_vertex_ptrs[0] ||
              tgt_ptr == face.face_vertex_ptrs[1] ||
              tgt_ptr == face.face_vertex_ptrs[2])
            continue; // 隣接は後でスカラーコールバック

          // near/far判定
          const double dx_c = face.Xc[0] - Xt_x;
          const double dy_c = face.Xc[1] - Xt_y;
          const double dz_c = face.Xc[2] - Xt_z;
          const double dist_sq = dx_c * dx_c + dy_c * dy_c + dz_c * dz_c;
          const bool is_near = (dist_sq < near_region_sq);

          // Dunavant テーブル選択: far→D6, near→D(p2m設定)
          int num_gp;
          const double (*gp_N)[3];
          const double* gp_w;
          if (!is_near) {
            num_gp = 6;
            gp_N = dun6.N;
            gp_w = dun6.w_eff;
          } else {
            num_gp = num_near_gp;
            gp_N = near_N_buf;
            gp_w = near_w_buf;
          }

          double WGN0 = 0., WGN1 = 0., WGN2 = 0.;
          double WGnN0 = 0., WGnN1 = 0., WGnN2 = 0.;

          for (int gp = 0; gp < num_gp; ++gp) {
            const double N0 = gp_N[gp][0];
            const double N1v = gp_N[gp][1];
            const double N2v = gp_N[gp][2];
            const double weff = gp_w[gp];

            const double Xgp_x = N0 * face.X0[0] + N1v * face.X1[0] + N2v * face.X2[0];
            const double Xgp_y = N0 * face.X0[1] + N1v * face.X1[1] + N2v * face.X2[1];
            const double Xgp_z = N0 * face.X0[2] + N1v * face.X1[2] + N2v * face.X2[2];

            const double Rx = Xgp_x - Xt_x;
            const double Ry = Xgp_y - Xt_y;
            const double Rz = Xgp_z - Xt_z;

            const double nr2 = Rx * Rx + Ry * Ry + Rz * Rz;
            const double nr_inv = 1.0 / std::sqrt(nr2);

            const double w_nr = weff * nr_inv;
            const double dot_R_c = Rx * face.cross[0] + Ry * face.cross[1] + Rz * face.cross[2];

            const double Jdet_w_nr = face.J_det * w_nr;
            const double neg_dot_w_nr3 = -(dot_R_c * w_nr * nr_inv * nr_inv);

            WGN0 += Jdet_w_nr * N0;
            WGN1 += Jdet_w_nr * N1v;
            WGN2 += Jdet_w_nr * N2v;

            WGnN0 += neg_dot_w_nr3 * N0;
            WGnN1 += neg_dot_w_nr3 * N1v;
            WGnN2 += neg_dot_w_nr3 * N2v;
          }

          acc.add(face.idx[0], WGN0, WGnN0);
          acc.add(face.idx[1], WGN1, WGnN1);
          acc.add(face.idx[2], WGN2, WGnN2);
        }

        // 隣接 Linear ソース + 全 Non-linear ソース
        for (const auto& ns : all_sources) {
          if (ns.linear_idx >= 0) {
            if (ns.ptr->isAdjacentTo(tgt)) {
              ns.ptr->fill_direct_entries(tgt, acc);
            }
          } else {
            if (ns.ptr->fill_direct_entries_nonadj && !ns.ptr->isAdjacentTo(tgt)) {
              ns.ptr->fill_direct_entries_nonadj(tgt, acc);
            } else {
              ns.ptr->fill_direct_entries(tgt, acc);
            }
          }
        }

        // Phase 4: 結果抽出
        std::sort(acc.touched_list.begin(), acc.touched_list.end());

        tgt->near_indices.clear();
        tgt->near_weights_phi.clear();
        tgt->near_weights_phin.clear();

        const size_t n_touched = acc.touched_list.size();
        tgt->near_indices.reserve(n_touched);
        tgt->near_weights_phi.reserve(n_touched);
        tgt->near_weights_phin.reserve(n_touched);

        for (int32_t idx : acc.touched_list) {
          tgt->near_indices.push_back(idx);
          tgt->near_weights_phi.push_back(acc.phi_acc[idx]);
          tgt->near_weights_phin.push_back(acc.phin_acc[idx]);
          acc.phi_acc[idx] = 0.;
          acc.phin_acc[idx] = 0.;
          acc.touched_flags[idx] = 0;
        }
        acc.update_high_water_mark();

        const size_t n = tgt->near_indices.size();
        tgt->near_run_base_idx.clear();
        tgt->near_run_pos.clear();
        tgt->near_run_len.clear();
        if (n > 0) {
          tgt->near_run_base_idx.reserve(n / 4 + 1);
          tgt->near_run_pos.reserve(n / 4 + 1);
          tgt->near_run_len.reserve(n / 4 + 1);
          size_t pos0 = 0;
          int32_t base = tgt->near_indices[0];
          int32_t prev = base;
          for (size_t i = 1; i < n; ++i) {
            const int32_t cur = tgt->near_indices[i];
            if (cur == prev + 1) {
              prev = cur;
              continue;
            }
            tgt->near_run_base_idx.push_back(base);
            tgt->near_run_pos.push_back(static_cast<int32_t>(pos0));
            tgt->near_run_len.push_back(static_cast<int32_t>(i - pos0));
            pos0 = i;
            base = cur;
            prev = cur;
          }
          tgt->near_run_base_idx.push_back(base);
          tgt->near_run_pos.push_back(static_cast<int32_t>(pos0));
          tgt->near_run_len.push_back(static_cast<int32_t>(n - pos0));
        }

        tgt->near_cell_count = cell_count;

        const size_t done = processed_targets.fetch_add(1, std::memory_order_relaxed) + 1;
        size_t expected = next_report.load(std::memory_order_relaxed);
        while (done >= expected && expected <= total_targets) {
          if (next_report.compare_exchange_weak(expected,
                                                expected + report_interval,
                                                std::memory_order_relaxed)) {
#pragma omp critical(nearfield_progress_log)
            {
              const auto now = std::chrono::steady_clock::now();
              const double elapsed =
                  std::chrono::duration<double>(now - progress_t0).count();
              const int percent = static_cast<int>((100.0 * done) / static_cast<double>(total_targets));
              std::cout << " [nearfield:cell_scalar] progress " << done << "/"
                        << total_targets << " (" << percent << "%), elapsed="
                        << elapsed << " s" << std::endl;
            }
            break;
          }
        }
      } // ti loop
    } // gi loop
  } // omp parallel
}

// ---------------------------------------------------------------------------
// セルグループ化なし + フラットデータ構造版 スカラー setDirectIntegration
//
// cell_scalar からバケットグループ化を除去。各ターゲットが独立にソースを収集。
// 「フラットデータ構造（LinearFaceDataD）の効果」と
// 「セルグループ化（ソースリスト共有）の効果」を分離するためのベンチマーク用。
// ---------------------------------------------------------------------------
template <typename BucketType, typename TargetPtr>
void setDirectIntegrationFlatScalar_linear(
    BucketType& B_poles,
    const std::vector<TargetPtr>& targets,
    int near_quadrature_points = 6) {
  using target_t = std::remove_pointer_t<decltype(&*targets[0])>;
  using source_t = typename decltype(B_poles.data1D)::value_type::element_type;

  const auto& dun6 = GaussTableSIMD_D::dun6;
  const auto& near_rule = getDunavantP2MRule(near_quadrature_points);
  const int num_near_gp = static_cast<int>(near_rule.size());
  constexpr int MAX_DUN = 25;
  double near_N_buf[MAX_DUN][3];
  double near_w_buf[MAX_DUN];
  for (int i = 0; i < num_near_gp; ++i) {
    near_N_buf[i][0] = near_rule[i][0];
    near_N_buf[i][1] = near_rule[i][1];
    near_N_buf[i][2] = near_rule[i][2];
    near_w_buf[i] = near_rule[i][3];
  }

  const size_t total_targets = targets.size();
  std::atomic<size_t> processed_targets{0};
  const size_t report_interval = std::max<size_t>(1, total_targets / 20); // ~5%
  std::atomic<size_t> next_report{report_interval};
  const auto progress_t0 = std::chrono::steady_clock::now();

  std::cout << " [nearfield:flat_scalar] begin, targets=" << total_targets
            << ", report_interval=" << report_interval << std::endl;

#pragma omp parallel
  {
    DirectAccumulator acc;
    std::vector<LinearFaceDataD> linear_faces;
    struct SourceRef {
      source_t* ptr;
      int32_t linear_idx;
    };
    std::vector<SourceRef> all_sources;

#pragma omp for schedule(dynamic, 4)
    for (size_t ti = 0; ti < targets.size(); ++ti) {
      auto* tgt = &*targets[ti];
      acc.reset();
      linear_faces.clear();
      all_sources.clear();

      auto* b_deepest = B_poles.getBucketAtDeepest(tgt->Xtarget);

      double near_region_sq = 0.;
      size_t cell_count = 0;

      // ---- 各ターゲットが独立にソース収集 ----
      auto collectFromBucket = [&](const auto* B) {
        ++cell_count;
        for (const auto& source : B->data1D_vector) {
          if (!source->fill_direct_entries)
            continue;
          if (source->dof_indices[0] >= 0) {
            const auto& fv = source->face_vertices;
            LinearFaceDataD fd;
            for (int k = 0; k < 3; ++k) {
              fd.X0[k] = fv[0]->Xtarget[k];
              fd.X1[k] = fv[1]->Xtarget[k];
              fd.X2[k] = fv[2]->Xtarget[k];
            }
            double dx1[3], dx2[3];
            for (int k = 0; k < 3; ++k) {
              dx1[k] = fv[1]->Xtarget[k] - fv[0]->Xtarget[k];
              dx2[k] = fv[2]->Xtarget[k] - fv[0]->Xtarget[k];
            }
            double cx = dx1[1] * dx2[2] - dx1[2] * dx2[1];
            double cy = dx1[2] * dx2[0] - dx1[0] * dx2[2];
            double cz = dx1[0] * dx2[1] - dx1[1] * dx2[0];
            fd.J_det = std::sqrt(cx * cx + cy * cy + cz * cz);
            fd.cross[0] = cx;
            fd.cross[1] = cy;
            fd.cross[2] = cz;
            fd.Xc[0] = (fd.X0[0] + fd.X1[0] + fd.X2[0]) / 3.0;
            fd.Xc[1] = (fd.X0[1] + fd.X1[1] + fd.X2[1]) / 3.0;
            fd.Xc[2] = (fd.X0[2] + fd.X1[2] + fd.X2[2]) / 3.0;
            fd.idx[0] = source->dof_indices[0];
            fd.idx[1] = source->dof_indices[1];
            fd.idx[2] = source->dof_indices[2];
            fd.face_vertex_ptrs[0] = static_cast<const void*>(fv[0]);
            fd.face_vertex_ptrs[1] = static_cast<const void*>(fv[1]);
            fd.face_vertex_ptrs[2] = static_cast<const void*>(fv[2]);
            int32_t fidx = static_cast<int32_t>(linear_faces.size());
            linear_faces.push_back(fd);
            all_sources.push_back(SourceRef{source.get(), fidx});
            near_region_sq = source->near_region * source->near_region;
          } else {
            all_sources.push_back(SourceRef{source.get(), -1});
          }
        }
      };

      collectFromBucket(b_deepest);
      for (const auto& b : b_deepest->buckets_near)
        collectFromBucket(b);
      {
        auto* bp = b_deepest->parent;
        while (bp != nullptr) {
          for (auto* B : bp->buckets_near)
            if (!B->hasChildren())
              collectFromBucket(B);
          bp = bp->parent;
        }
      }

      const size_t num_linear_faces = linear_faces.size();

      if (all_sources.empty()) {
        tgt->near_indices.clear();
        tgt->near_weights_phi.clear();
        tgt->near_weights_phin.clear();
        tgt->near_run_base_idx.clear();
        tgt->near_run_pos.clear();
        tgt->near_run_len.clear();
        tgt->near_cell_count = cell_count;

        const size_t done = processed_targets.fetch_add(1, std::memory_order_relaxed) + 1;
        size_t expected = next_report.load(std::memory_order_relaxed);
        while (done >= expected && expected <= total_targets) {
          if (next_report.compare_exchange_weak(expected,
                                                expected + report_interval,
                                                std::memory_order_relaxed)) {
#pragma omp critical(nearfield_progress_log)
            {
              const auto now = std::chrono::steady_clock::now();
              const double elapsed =
                  std::chrono::duration<double>(now - progress_t0).count();
              const int percent = static_cast<int>((100.0 * done) / static_cast<double>(total_targets));
              std::cout << " [nearfield:flat_scalar] progress " << done << "/"
                        << total_targets << " (" << percent << "%), elapsed="
                        << elapsed << " s" << std::endl;
            }
            break;
          }
        }
        continue;
      }

      const double Xt_x = tgt->Xtarget[0];
      const double Xt_y = tgt->Xtarget[1];
      const double Xt_z = tgt->Xtarget[2];
      const void* tgt_ptr = static_cast<const void*>(tgt);

      // 非隣接 Linear フェイス：スカラー double で積分
      for (size_t fi = 0; fi < num_linear_faces; ++fi) {
        const auto& face = linear_faces[fi];

        if (tgt_ptr == face.face_vertex_ptrs[0] ||
            tgt_ptr == face.face_vertex_ptrs[1] ||
            tgt_ptr == face.face_vertex_ptrs[2])
          continue;

        const double dx_c = face.Xc[0] - Xt_x;
        const double dy_c = face.Xc[1] - Xt_y;
        const double dz_c = face.Xc[2] - Xt_z;
        const double dist_sq = dx_c * dx_c + dy_c * dy_c + dz_c * dz_c;
        const bool is_near = (dist_sq < near_region_sq);

        // Dunavant テーブル選択: far→D6, near→D(p2m設定)
        int num_gp;
        const double (*gp_N)[3];
        const double* gp_w;
        if (!is_near) {
          num_gp = 6;
          gp_N = dun6.N;
          gp_w = dun6.w_eff;
        } else {
          num_gp = num_near_gp;
          gp_N = near_N_buf;
          gp_w = near_w_buf;
        }

        double WGN0 = 0., WGN1 = 0., WGN2 = 0.;
        double WGnN0 = 0., WGnN1 = 0., WGnN2 = 0.;

        for (int gp = 0; gp < num_gp; ++gp) {
          const double N0 = gp_N[gp][0];
          const double N1v = gp_N[gp][1];
          const double N2v = gp_N[gp][2];
          const double weff = gp_w[gp];

          const double Xgp_x = N0 * face.X0[0] + N1v * face.X1[0] + N2v * face.X2[0];
          const double Xgp_y = N0 * face.X0[1] + N1v * face.X1[1] + N2v * face.X2[1];
          const double Xgp_z = N0 * face.X0[2] + N1v * face.X1[2] + N2v * face.X2[2];

          const double Rx = Xgp_x - Xt_x;
          const double Ry = Xgp_y - Xt_y;
          const double Rz = Xgp_z - Xt_z;

          const double nr2 = Rx * Rx + Ry * Ry + Rz * Rz;
          const double nr_inv = 1.0 / std::sqrt(nr2);

          const double w_nr = weff * nr_inv;
          const double dot_R_c = Rx * face.cross[0] + Ry * face.cross[1] + Rz * face.cross[2];

          const double Jdet_w_nr = face.J_det * w_nr;
          const double neg_dot_w_nr3 = -(dot_R_c * w_nr * nr_inv * nr_inv);

          WGN0 += Jdet_w_nr * N0;
          WGN1 += Jdet_w_nr * N1v;
          WGN2 += Jdet_w_nr * N2v;

          WGnN0 += neg_dot_w_nr3 * N0;
          WGnN1 += neg_dot_w_nr3 * N1v;
          WGnN2 += neg_dot_w_nr3 * N2v;
        }

        acc.add(face.idx[0], WGN0, WGnN0);
        acc.add(face.idx[1], WGN1, WGnN1);
        acc.add(face.idx[2], WGN2, WGnN2);
      }

      // 隣接 Linear ソース + 全 Non-linear ソース
      for (const auto& ns : all_sources) {
        if (ns.linear_idx >= 0) {
          if (ns.ptr->isAdjacentTo(tgt)) {
            ns.ptr->fill_direct_entries(tgt, acc);
          }
        } else {
          if (ns.ptr->fill_direct_entries_nonadj && !ns.ptr->isAdjacentTo(tgt)) {
            ns.ptr->fill_direct_entries_nonadj(tgt, acc);
          } else {
            ns.ptr->fill_direct_entries(tgt, acc);
          }
        }
      }

      // 結果抽出
      std::sort(acc.touched_list.begin(), acc.touched_list.end());

      tgt->near_indices.clear();
      tgt->near_weights_phi.clear();
      tgt->near_weights_phin.clear();

      const size_t n_touched = acc.touched_list.size();
      tgt->near_indices.reserve(n_touched);
      tgt->near_weights_phi.reserve(n_touched);
      tgt->near_weights_phin.reserve(n_touched);

      for (int32_t idx : acc.touched_list) {
        tgt->near_indices.push_back(idx);
        tgt->near_weights_phi.push_back(acc.phi_acc[idx]);
        tgt->near_weights_phin.push_back(acc.phin_acc[idx]);
        acc.phi_acc[idx] = 0.;
        acc.phin_acc[idx] = 0.;
        acc.touched_flags[idx] = 0;
      }
      acc.update_high_water_mark();

      const size_t n = tgt->near_indices.size();
      tgt->near_run_base_idx.clear();
      tgt->near_run_pos.clear();
      tgt->near_run_len.clear();
      if (n > 0) {
        tgt->near_run_base_idx.reserve(n / 4 + 1);
        tgt->near_run_pos.reserve(n / 4 + 1);
        tgt->near_run_len.reserve(n / 4 + 1);
        size_t pos0 = 0;
        int32_t base = tgt->near_indices[0];
        int32_t prev = base;
        for (size_t i = 1; i < n; ++i) {
          const int32_t cur = tgt->near_indices[i];
          if (cur == prev + 1) {
            prev = cur;
            continue;
          }
          tgt->near_run_base_idx.push_back(base);
          tgt->near_run_pos.push_back(static_cast<int32_t>(pos0));
          tgt->near_run_len.push_back(static_cast<int32_t>(i - pos0));
          pos0 = i;
          base = cur;
          prev = cur;
        }
        tgt->near_run_base_idx.push_back(base);
        tgt->near_run_pos.push_back(static_cast<int32_t>(pos0));
        tgt->near_run_len.push_back(static_cast<int32_t>(n - pos0));
      }

      tgt->near_cell_count = cell_count;

      const size_t done = processed_targets.fetch_add(1, std::memory_order_relaxed) + 1;
      size_t expected = next_report.load(std::memory_order_relaxed);
      while (done >= expected && expected <= total_targets) {
        if (next_report.compare_exchange_weak(expected,
                                              expected + report_interval,
                                              std::memory_order_relaxed)) {
#pragma omp critical(nearfield_progress_log)
          {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed =
                std::chrono::duration<double>(now - progress_t0).count();
            const int percent = static_cast<int>((100.0 * done) / static_cast<double>(total_targets));
            std::cout << " [nearfield:flat_scalar] progress " << done << "/"
                      << total_targets << " (" << percent << "%), elapsed="
                      << elapsed << " s" << std::endl;
          }
          break;
        }
      }
    } // ti loop
  } // omp parallel
}

#endif // SET_DIRECT_INTEGRATION_SIMD_DOUBLE_HPP
