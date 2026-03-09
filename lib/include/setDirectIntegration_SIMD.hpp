#ifndef SET_DIRECT_INTEGRATION_SIMD_HPP
#define SET_DIRECT_INTEGRATION_SIMD_HPP

#include <arm_neon.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

/*
 * setDirectIntegration の SIMD最適化版（Linear要素専用）
 *
 * 同一セル（バケット）内のターゲットは同じ近傍ソースリストを共有する。
 * この性質を利用して：
 *   1. ターゲットをバケットごとにグループ化
 *   2. ソースフェイスデータを事前にフラット化（LinearFaceData）
 *   3. NEON float32x4_t で4ターゲット同時処理（非隣接のみ）
 *   4. 隣接ペアは既存スカラーコールバックで処理
 *   5. 非線形ソースも既存スカラーコールバックで処理
 */

// ---------------------------------------------------------------------------
// ソースフェイスデータ（Linear要素用）
// ---------------------------------------------------------------------------
struct LinearFaceData {
   float X0[3], X1[3], X2[3]; // 頂点座標
   float J_det;                // ヤコビアン = 2 * area
   float cross[3];             // 2 * area * normal (= Cross(X1-X0, X2-X0))
   float Xc[3];                // 重心座標
   int32_t idx[3];             // グローバルDOFインデックス
   const void* face_vertex_ptrs[3]; // 隣接判定用ポインタ (target4FMM*)
};

// ---------------------------------------------------------------------------
// 事前計算ガウス点テーブル（全フェイス共通の定数）
// ---------------------------------------------------------------------------
namespace GaussTableSIMD {

// ModTriShape<3>(t0, t1) = {t0, t1*(1-t0), (t0-1)*(t1-1)}
constexpr std::array<float, 3> modTriShape3(double t0, double t1) {
   return {static_cast<float>(t0),
           static_cast<float>(t1 * (1.0 - t0)),
           static_cast<float>((t0 - 1.0) * (t1 - 1.0))};
}

// 25点ガウス求積点 (GW5xGW5)
// N[gp][j] = ModTriShape<3>(t0, t1)[j] — 形状関数値（ソース点位置 + 重み用）
// w_eff[gp] = ww * (1 - t0) — 有効重み（パラメトリックヤコビアン含む）
struct GW5x5Table {
   float N[25][3];
   float w_eff[25];

   constexpr GW5x5Table() : N{}, w_eff{} {
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
         w_eff[i] = static_cast<float>(pts[i][2] * (1.0 - pts[i][0]));
      }
   }
};

// 1点（重心）ガウス求積 — far-field用
// N = {1/3, 1/3, 1/3} （標準線形形状関数, NOT ModTriShape<3>）
// w_eff = 0.5 （元コードの far-field パスに合わせる）
// N[1][3]/w_eff[1] の配列形状は GW5x5Table とポインタ互換にするため
struct GW1x1Table {
   float N[1][3];
   float w_eff[1];

   constexpr GW1x1Table() : N{}, w_eff{} {
      N[0][0] = 1.0f / 3.0f;
      N[0][1] = 1.0f / 3.0f;
      N[0][2] = 1.0f / 3.0f;
      w_eff[0] = 0.5f;
   }
};

inline constexpr GW5x5Table gw5x5{};
inline constexpr GW1x1Table gw1x1{};

} // namespace GaussTableSIMD

// ---------------------------------------------------------------------------
// NEON ヘルパー
// ---------------------------------------------------------------------------
namespace neon_helpers {

// Newton-Raphson refined reciprocal square root (2 iterations)
// Input x, returns 1/sqrt(x) with ~23-bit precision
inline float32x4_t rsqrt_nr(float32x4_t x) {
   float32x4_t est = vrsqrteq_f32(x);
   est = vmulq_f32(est, vrsqrtsq_f32(vmulq_f32(x, est), est));
   est = vmulq_f32(est, vrsqrtsq_f32(vmulq_f32(x, est), est));
   return est;
}

} // namespace neon_helpers

// ---------------------------------------------------------------------------
// SIMD版 setDirectIntegration（Linear要素専用）
//
// scalar版 target4FMM::setDirectIntegration() の完全互換置換。
// 出力: 各ターゲットの near_indices, near_weights_phi, near_weights_phin,
//       near_run_base_idx, near_run_pos, near_run_len を設定する。
// ---------------------------------------------------------------------------
template <typename BucketType, typename TargetPtr>
void setDirectIntegrationSIMD_linear(
    BucketType& B_poles,
    const std::vector<TargetPtr>& targets) {
   using target_t = std::remove_pointer_t<decltype(&*targets[0])>;
   // Derive source type from bucket's stored element (shared_ptr<source4FMM<target4FMM>>)
   using source_t = typename decltype(B_poles.data1D)::value_type::element_type;

   const auto& gw25 = GaussTableSIMD::gw5x5;
   const auto& gw1 = GaussTableSIMD::gw1x1;

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
      // Thread-local accumulators (4 per thread, reused across batches)
      DirectAccumulator accs[4];

#pragma omp for schedule(dynamic, 1)
      for (size_t gi = 0; gi < groups.size(); ++gi) {
         auto& tgts = groups[gi].tgts;
         const size_t num_targets = tgts.size();

         // Get the deepest bucket (same for all targets in this group)
         auto* b_deepest = B_poles.getBucketAtDeepest(tgts[0]->Xtarget);

         // ---- 近傍ソース収集 ----
         // LinearFaceData: SIMD処理用（Linear要素のみ）
         // source_ptrs: スカラーフォールバック用（全ソース）
         std::vector<LinearFaceData> linear_faces;
         struct SourceRef {
            source_t* ptr;
            int32_t linear_idx; // linear_faces内のインデックス, or -1
         };
         std::vector<SourceRef> all_sources;
         float near_region_sq = 0.f;
         size_t cell_count = 0;

         auto collectFromBucket = [&](const auto* B) {
            ++cell_count;
            for (const auto& source : B->data1D_vector) {
               if (!source->fill_direct_entries)
                  continue;

               if (source->dof_indices[0] >= 0) {
                  // Linear source → LinearFaceData構築
                  const auto& fv = source->face_vertices;
                  LinearFaceData fd;
                  for (int k = 0; k < 3; ++k) {
                     fd.X0[k] = static_cast<float>(fv[0]->Xtarget[k]);
                     fd.X1[k] = static_cast<float>(fv[1]->Xtarget[k]);
                     fd.X2[k] = static_cast<float>(fv[2]->Xtarget[k]);
                  }
                  // cross = (X1-X0) × (X2-X0)
                  double dx1[3], dx2[3];
                  for (int k = 0; k < 3; ++k) {
                     dx1[k] = fv[1]->Xtarget[k] - fv[0]->Xtarget[k];
                     dx2[k] = fv[2]->Xtarget[k] - fv[0]->Xtarget[k];
                  }
                  double cx = dx1[1] * dx2[2] - dx1[2] * dx2[1];
                  double cy = dx1[2] * dx2[0] - dx1[0] * dx2[2];
                  double cz = dx1[0] * dx2[1] - dx1[1] * dx2[0];
                  fd.J_det = static_cast<float>(std::sqrt(cx * cx + cy * cy + cz * cz));
                  fd.cross[0] = static_cast<float>(cx);
                  fd.cross[1] = static_cast<float>(cy);
                  fd.cross[2] = static_cast<float>(cz);
                  fd.Xc[0] = (fd.X0[0] + fd.X1[0] + fd.X2[0]) / 3.0f;
                  fd.Xc[1] = (fd.X0[1] + fd.X1[1] + fd.X2[1]) / 3.0f;
                  fd.Xc[2] = (fd.X0[2] + fd.X1[2] + fd.X2[2]) / 3.0f;
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
                  // Non-linear source → スカラーのみ
                  all_sources.push_back(SourceRef{source.get(), -1});
               }
            }
         };

         // setDirectIntegration と同じバケット走査パターン
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

         // ---- ターゲットを4個ずつ処理 ----
         for (size_t t_base = 0; t_base < num_targets; t_base += 4) {
            const size_t t_end = std::min(t_base + 4, num_targets);
            const size_t t_count = t_end - t_base;

            // Accumulator をリセット
            for (size_t i = 0; i < t_count; ++i)
               accs[i].reset();

            // ターゲット座標をSIMDレジスタにロード
            float tx[4], ty[4], tz[4];
            for (size_t i = 0; i < t_count; ++i) {
               tx[i] = static_cast<float>(tgts[t_base + i]->Xtarget[0]);
               ty[i] = static_cast<float>(tgts[t_base + i]->Xtarget[1]);
               tz[i] = static_cast<float>(tgts[t_base + i]->Xtarget[2]);
            }
            for (size_t i = t_count; i < 4; ++i) {
               tx[i] = tx[t_count - 1];
               ty[i] = ty[t_count - 1];
               tz[i] = tz[t_count - 1];
            }
            float32x4_t Xt_x = vld1q_f32(tx);
            float32x4_t Xt_y = vld1q_f32(ty);
            float32x4_t Xt_z = vld1q_f32(tz);

            // 隣接判定用のターゲットポインタ
            const void* tgt_ptrs[4];
            for (size_t i = 0; i < t_count; ++i)
               tgt_ptrs[i] = static_cast<const void*>(tgts[t_base + i]);
            for (size_t i = t_count; i < 4; ++i)
               tgt_ptrs[i] = tgt_ptrs[t_count - 1];

            // ==============================================================
            // SIMD pass: 非隣接 Linear フェイス
            // ==============================================================
            for (size_t fi = 0; fi < num_linear_faces; ++fi) {
               const auto& face = linear_faces[fi];

               // 隣接判定（各レーン独立）
               bool is_adj[4];
               for (int l = 0; l < 4; ++l) {
                  is_adj[l] = (tgt_ptrs[l] == face.face_vertex_ptrs[0] ||
                               tgt_ptrs[l] == face.face_vertex_ptrs[1] ||
                               tgt_ptrs[l] == face.face_vertex_ptrs[2]);
               }

               // near/far判定: 重心とターゲットの距離²
               float32x4_t dx_c = vsubq_f32(vdupq_n_f32(face.Xc[0]), Xt_x);
               float32x4_t dy_c = vsubq_f32(vdupq_n_f32(face.Xc[1]), Xt_y);
               float32x4_t dz_c = vsubq_f32(vdupq_n_f32(face.Xc[2]), Xt_z);
               float32x4_t dist_sq = vfmaq_f32(vfmaq_f32(vmulq_f32(dx_c, dx_c), dy_c, dy_c), dz_c, dz_c);
               uint32x4_t near_mask = vcltq_f32(dist_sq, vdupq_n_f32(near_region_sq));

               uint32_t mask_bits[4];
               vst1q_u32(mask_bits, near_mask);
               const bool all_far = (!mask_bits[0] && !mask_bits[1] && !mask_bits[2] && !mask_bits[3]);

               // ガウス点テーブル選択
               int num_gp;
               const float(*gp_N)[3];
               const float* gp_w;
               if (all_far) {
                  num_gp = 1;
                  gp_N = gw1.N;
                  gp_w = gw1.w_eff;
               } else {
                  num_gp = 25;
                  gp_N = gw25.N;
                  gp_w = gw25.w_eff;
               }

               // ソースフェイス定数（broadcast）
               float32x4_t f_Jdet = vdupq_n_f32(face.J_det);
               float32x4_t f_cx = vdupq_n_f32(face.cross[0]);
               float32x4_t f_cy = vdupq_n_f32(face.cross[1]);
               float32x4_t f_cz = vdupq_n_f32(face.cross[2]);

               // Per-face, per-lane, per-DOF accumulator
               float32x4_t WGN0 = vdupq_n_f32(0.f), WGN1 = vdupq_n_f32(0.f), WGN2 = vdupq_n_f32(0.f);
               float32x4_t WGnN0 = vdupq_n_f32(0.f), WGnN1 = vdupq_n_f32(0.f), WGnN2 = vdupq_n_f32(0.f);

               for (int gp = 0; gp < num_gp; ++gp) {
                  const float N0 = gp_N[gp][0];
                  const float N1v = gp_N[gp][1];
                  const float N2v = gp_N[gp][2];
                  const float weff = gp_w[gp];

                  // ソースガウス点位置（全レーン共通, broadcast）
                  const float Xgp_x = N0 * face.X0[0] + N1v * face.X1[0] + N2v * face.X2[0];
                  const float Xgp_y = N0 * face.X0[1] + N1v * face.X1[1] + N2v * face.X2[1];
                  const float Xgp_z = N0 * face.X0[2] + N1v * face.X1[2] + N2v * face.X2[2];

                  // R = Xgp - Xtarget (4 lanes SIMD)
                  float32x4_t Rx = vsubq_f32(vdupq_n_f32(Xgp_x), Xt_x);
                  float32x4_t Ry = vsubq_f32(vdupq_n_f32(Xgp_y), Xt_y);
                  float32x4_t Rz = vsubq_f32(vdupq_n_f32(Xgp_z), Xt_z);

                  // |R|² and 1/|R| (Newton-Raphson refined)
                  float32x4_t nr2 = vfmaq_f32(vfmaq_f32(vmulq_f32(Rx, Rx), Ry, Ry), Rz, Rz);
                  float32x4_t nr_inv = neon_helpers::rsqrt_nr(nr2);

                  // w_eff / |R|
                  float32x4_t w_nr = vmulq_f32(vdupq_n_f32(weff), nr_inv);

                  // dot(R, cross) = Rx*cx + Ry*cy + Rz*cz
                  float32x4_t dot_R_c = vfmaq_f32(vfmaq_f32(vmulq_f32(Rx, f_cx), Ry, f_cy), Rz, f_cz);

                  // WG: J_det * w_eff * nr_inv (per Gauss point, broadcast shape fn later)
                  float32x4_t Jdet_w_nr = vmulq_f32(f_Jdet, w_nr);

                  // WGn: -dot(R,cross) * w_eff * nr_inv³
                  float32x4_t nr_inv2 = vmulq_f32(nr_inv, nr_inv);
                  float32x4_t neg_dot_w_nr3 = vnegq_f32(vmulq_f32(vmulq_f32(dot_R_c, w_nr), nr_inv2));

                  // Accumulate per DOF (shape function as scalar broadcast)
                  WGN0 = vfmaq_f32(WGN0, Jdet_w_nr, vdupq_n_f32(N0));
                  WGN1 = vfmaq_f32(WGN1, Jdet_w_nr, vdupq_n_f32(N1v));
                  WGN2 = vfmaq_f32(WGN2, Jdet_w_nr, vdupq_n_f32(N2v));

                  WGnN0 = vfmaq_f32(WGnN0, neg_dot_w_nr3, vdupq_n_f32(N0));
                  WGnN1 = vfmaq_f32(WGnN1, neg_dot_w_nr3, vdupq_n_f32(N1v));
                  WGnN2 = vfmaq_f32(WGnN2, neg_dot_w_nr3, vdupq_n_f32(N2v));
               } // gp loop

               // Scatter 結果を各レーンの accumulator に書き出し
               float wgn0[4], wgn1[4], wgn2[4];
               float wgnn0[4], wgnn1[4], wgnn2[4];
               vst1q_f32(wgn0, WGN0);
               vst1q_f32(wgn1, WGN1);
               vst1q_f32(wgn2, WGN2);
               vst1q_f32(wgnn0, WGnN0);
               vst1q_f32(wgnn1, WGnN1);
               vst1q_f32(wgnn2, WGnN2);

               for (size_t lane = 0; lane < t_count; ++lane) {
                  if (is_adj[lane]) continue; // 隣接はスカラーで別途処理
                  accs[lane].add(face.idx[0], static_cast<double>(wgn0[lane]), static_cast<double>(wgnn0[lane]));
                  accs[lane].add(face.idx[1], static_cast<double>(wgn1[lane]), static_cast<double>(wgnn1[lane]));
                  accs[lane].add(face.idx[2], static_cast<double>(wgn2[lane]), static_cast<double>(wgnn2[lane]));
               }
            } // fi (linear faces) loop

            // ==============================================================
            // Scalar pass: 隣接 Linear ソース + 全 Non-linear ソース
            // ==============================================================
            for (size_t lane = 0; lane < t_count; ++lane) {
               auto* tgt = tgts[t_base + lane];
               for (const auto& ns : all_sources) {
                  if (ns.linear_idx >= 0) {
                     // Linear source: 隣接ペアのみスカラー処理
                     if (ns.ptr->isAdjacentTo(tgt)) {
                        ns.ptr->fill_direct_entries(tgt, accs[lane]);
                     }
                  } else {
                     // Non-linear source: 常にスカラー処理
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

               // Sort by DOF index
               std::sort(acc.touched_list.begin(), acc.touched_list.end());

               // Extract to sparse arrays
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
                  // Reset accumulator entries for reuse
                  acc.phi_acc[idx] = 0.;
                  acc.phin_acc[idx] = 0.;
                  acc.touched_flags[idx] = 0;
               }
               acc.update_high_water_mark();

               // Build RLE (Run-Length Encoding)
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
            } // lane loop (result extraction)
         } // t_base loop
      } // gi loop (bucket groups)
   } // omp parallel
}

#endif // SET_DIRECT_INTEGRATION_SIMD_HPP
