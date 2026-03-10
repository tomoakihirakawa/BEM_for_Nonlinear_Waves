/*
 * metal_nearfield_wrapper.hpp
 *
 * C++ wrapper for Metal-accelerated nearfield direct integration.
 * This header can be included in GCC-compiled code.
 *
 * Provides setDirectIntegrationMetal_linear() as a drop-in replacement
 * for the scalar setDirectIntegration() in the FMM nearfield step.
 *
 * Architecture:
 *   - Targets grouped by FMM bucket (cell group) — same near-field source list
 *   - Local DOF mapping: global DOFs → contiguous local indices per group
 *   - Non-adjacent pairs computed on GPU (float + Kahan summation)
 *   - Adjacent pairs handled on CPU via fill_direct_entries callbacks
 *   - Results mapped back to global DOF indices → near_indices/weights/RLE
 */

#ifndef METAL_NEARFIELD_WRAPPER_HPP
#define METAL_NEARFIELD_WRAPPER_HPP

#include "metal_nearfield.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

/*
 * Metal-accelerated setDirectIntegration for linear elements.
 *
 * Drop-in replacement for:
 *   #pragma omp parallel for
 *   for (auto& t : targets) t->setDirectIntegration(B_poles);
 *
 * Uses Metal GPU for non-adjacent Gauss quadrature integration.
 * Adjacent pairs and non-linear sources fall back to CPU callbacks.
 */
template <typename BucketType, typename TargetPtr>
void setDirectIntegrationMetal_linear(
    BucketType& B_poles,
    const std::vector<TargetPtr>& targets) {

   using target_t = std::remove_pointer_t<decltype(&*targets[0])>;
   using source_t = typename decltype(B_poles.data1D)::value_type::element_type;

   if (targets.empty()) return;

   // =======================================================================
   // Phase 1: Group targets by deepest bucket
   // =======================================================================
   struct CellGroup {
      std::vector<target_t*> tgts;
      const void* bucket_ptr;
   };
   std::unordered_map<const void*, size_t> bucket_to_group;
   std::vector<CellGroup> groups;

   for (auto& t : targets) {
      auto* b = B_poles.getBucketAtDeepest(t->Xtarget);
      auto it = bucket_to_group.find(static_cast<const void*>(b));
      if (it == bucket_to_group.end()) {
         bucket_to_group[static_cast<const void*>(b)] = groups.size();
         groups.push_back({{&*t}, static_cast<const void*>(b)});
      } else {
         groups[it->second].tgts.push_back(&*t);
      }
   }

   const size_t num_groups = groups.size();

   // =======================================================================
   // Phase 2: Collect source faces per group, build local DOF mapping
   // =======================================================================

   // Per-group data
   struct GroupData {
      std::vector<MetalFaceData> faces;         // GPU face data
      std::vector<source_t*> all_source_ptrs;   // all sources (for CPU fallback)
      std::vector<int32_t> linear_source_mask;  // 1 if linear, 0 if not
      std::unordered_map<int32_t, int32_t> global_to_local_dof;
      int32_t num_local_dofs = 0;
      float near_region_sq = 0.f;
      size_t cell_count = 0;
   };
   std::vector<GroupData> group_data(num_groups);

   // Collect sources for each group (sequential — bucket traversal)
   for (size_t gi = 0; gi < num_groups; ++gi) {
      auto& gd = group_data[gi];
      auto* b_deepest = B_poles.getBucketAtDeepest(groups[gi].tgts[0]->Xtarget);

      auto collectFromBucket = [&](const auto* B) {
         ++gd.cell_count;
         for (const auto& source : B->data1D_vector) {
            if (!source->fill_direct_entries)
               continue;

            gd.all_source_ptrs.push_back(source.get());

            if (source->dof_indices[0] >= 0) {
               // Linear source → build MetalFaceData
               gd.linear_source_mask.push_back(1);

               const auto& fv = source->face_vertices;
               MetalFaceData fd;

               for (int k = 0; k < 3; ++k) {
                  fd.X0[k] = static_cast<float>(fv[0]->Xtarget[k]);
                  fd.X1[k] = static_cast<float>(fv[1]->Xtarget[k]);
                  fd.X2[k] = static_cast<float>(fv[2]->Xtarget[k]);
               }

               double dx1[3], dx2[3];
               for (int k = 0; k < 3; ++k) {
                  dx1[k] = fv[1]->Xtarget[k] - fv[0]->Xtarget[k];
                  dx2[k] = fv[2]->Xtarget[k] - fv[0]->Xtarget[k];
               }
               double cx = dx1[1] * dx2[2] - dx1[2] * dx2[1];
               double cy = dx1[2] * dx2[0] - dx1[0] * dx2[2];
               double cz = dx1[0] * dx2[1] - dx1[1] * dx2[0];
               fd.J_det = static_cast<float>(std::sqrt(cx * cx + cy * cy + cz * cz));
               fd.cross_vec[0] = static_cast<float>(cx);
               fd.cross_vec[1] = static_cast<float>(cy);
               fd.cross_vec[2] = static_cast<float>(cz);
               fd.Xc[0] = (fd.X0[0] + fd.X1[0] + fd.X2[0]) / 3.0f;
               fd.Xc[1] = (fd.X0[1] + fd.X1[1] + fd.X2[1]) / 3.0f;
               fd.Xc[2] = (fd.X0[2] + fd.X1[2] + fd.X2[2]) / 3.0f;

               // Map global DOFs to local indices
               for (int j = 0; j < 3; ++j) {
                  int32_t gid = source->dof_indices[j];
                  auto [it, inserted] = gd.global_to_local_dof.try_emplace(gid, gd.num_local_dofs);
                  if (inserted) ++gd.num_local_dofs;
                  fd.local_dof[j] = it->second;
               }

               // Vertex IDs for adjacency check (use pointer as unique ID)
               for (int j = 0; j < 3; ++j)
                  fd.vtx_id[j] = static_cast<int32_t>(reinterpret_cast<uintptr_t>(fv[j]) & 0x7FFFFFFF);

               gd.faces.push_back(fd);
               gd.near_region_sq = source->near_region * source->near_region;
            } else {
               gd.linear_source_mask.push_back(0);
            }
         }
      };

      // Same bucket traversal pattern as setDirectIntegration
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
   }

   // =======================================================================
   // Phase 3: Flatten into contiguous GPU buffers
   // =======================================================================
   const int32_t total_targets = static_cast<int32_t>(targets.size());

   // Count totals
   int32_t total_faces = 0;
   int32_t total_result_slots = 0;
   for (size_t gi = 0; gi < num_groups; ++gi) {
      total_faces += static_cast<int32_t>(group_data[gi].faces.size());
      total_result_slots += static_cast<int32_t>(groups[gi].tgts.size()) * group_data[gi].num_local_dofs;
   }

   if (total_faces == 0) {
      // No linear faces at all — set empty results
      for (auto& t : targets) {
         t->near_indices.clear();
         t->near_weights_phi.clear();
         t->near_weights_phin.clear();
         t->near_run_base_idx.clear();
         t->near_run_pos.clear();
         t->near_run_len.clear();
         t->near_cell_count = 0;
      }
      return;
   }

   // Allocate flat arrays
   std::vector<MetalFaceData> all_faces;
   all_faces.reserve(total_faces);
   std::vector<int32_t> group_face_offset(num_groups + 1);
   std::vector<int32_t> group_dof_count(num_groups);
   std::vector<float> target_pos(total_targets * 3);
   std::vector<int32_t> target_group_id(total_targets);
   std::vector<int32_t> target_vtx_id(total_targets);
   std::vector<int32_t> target_dof_offset(total_targets);

   // Build group_face_offset and flatten faces
   group_face_offset[0] = 0;
   for (size_t gi = 0; gi < num_groups; ++gi) {
      all_faces.insert(all_faces.end(), group_data[gi].faces.begin(), group_data[gi].faces.end());
      group_face_offset[gi + 1] = static_cast<int32_t>(all_faces.size());
      group_dof_count[gi] = group_data[gi].num_local_dofs;
   }

   // Build per-target arrays
   // We need a global target index that matches the GPU buffer layout
   // Order: iterate groups, then targets within each group
   int32_t gpu_tid = 0;
   int32_t dof_offset = 0;
   // Map from target pointer to GPU target index
   std::unordered_map<target_t*, int32_t> target_to_gpu_tid;
   target_to_gpu_tid.reserve(total_targets);

   float max_near_region_sq = 0.f;

   for (size_t gi = 0; gi < num_groups; ++gi) {
      const auto& gd = group_data[gi];
      for (auto* tgt : groups[gi].tgts) {
         target_to_gpu_tid[tgt] = gpu_tid;
         target_pos[gpu_tid * 3 + 0] = static_cast<float>(tgt->Xtarget[0]);
         target_pos[gpu_tid * 3 + 1] = static_cast<float>(tgt->Xtarget[1]);
         target_pos[gpu_tid * 3 + 2] = static_cast<float>(tgt->Xtarget[2]);
         target_group_id[gpu_tid] = static_cast<int32_t>(gi);
         // Vertex ID for adjacency check (use pointer as unique ID)
         target_vtx_id[gpu_tid] = static_cast<int32_t>(reinterpret_cast<uintptr_t>(tgt) & 0x7FFFFFFF);
         target_dof_offset[gpu_tid] = dof_offset;
         dof_offset += gd.num_local_dofs;
         ++gpu_tid;
      }
      if (gd.near_region_sq > max_near_region_sq)
         max_near_region_sq = gd.near_region_sq;
   }

   // =======================================================================
   // Phase 4: GPU setup and compute
   // =======================================================================
   MetalNearfieldContext* ctx = metal_nearfield_init();
   if (ctx == nullptr) {
      std::cerr << "[Metal NF] Failed to initialize. Falling back to scalar." << std::endl;
      for (auto& t : targets)
         t->setDirectIntegration(B_poles);
      return;
   }

   int rc = metal_nearfield_setup(
       ctx,
       total_targets,
       total_faces,
       static_cast<int32_t>(num_groups),
       all_faces.data(),
       group_face_offset.data(),
       group_dof_count.data(),
       target_pos.data(),
       target_group_id.data(),
       target_vtx_id.data(),
       target_dof_offset.data(),
       total_result_slots,
       max_near_region_sq);

   if (rc != 0) {
      std::cerr << "[Metal NF] Setup failed. Falling back to scalar." << std::endl;
      metal_nearfield_destroy(ctx);
      for (auto& t : targets)
         t->setDirectIntegration(B_poles);
      return;
   }

   rc = metal_nearfield_compute(ctx);
   if (rc != 0) {
      std::cerr << "[Metal NF] Compute failed. Falling back to scalar." << std::endl;
      metal_nearfield_destroy(ctx);
      for (auto& t : targets)
         t->setDirectIntegration(B_poles);
      return;
   }

   // =======================================================================
   // Phase 5: Read back GPU results + CPU adjacent handling + result extraction
   // =======================================================================
   const float* phi_hi = metal_nearfield_get_phi_hi(ctx);
   const float* phi_lo = metal_nearfield_get_phi_lo(ctx);
   const float* phin_hi = metal_nearfield_get_phin_hi(ctx);
   const float* phin_lo = metal_nearfield_get_phin_lo(ctx);

   // Build per-group local→global DOF mapping (inverse of global_to_local)
   std::vector<std::vector<int32_t>> local_to_global(num_groups);
   for (size_t gi = 0; gi < num_groups; ++gi) {
      local_to_global[gi].resize(group_data[gi].num_local_dofs);
      for (const auto& [gid, lid] : group_data[gi].global_to_local_dof)
         local_to_global[gi][lid] = gid;
   }

#pragma omp parallel
   {
      DirectAccumulator acc;

#pragma omp for schedule(dynamic, 1)
      for (size_t gi = 0; gi < num_groups; ++gi) {
         const auto& gd = group_data[gi];
         const auto& l2g = local_to_global[gi];
         const int32_t ndofs = gd.num_local_dofs;

         for (auto* tgt : groups[gi].tgts) {
            acc.reset();

            int32_t gpu_id = target_to_gpu_tid.at(tgt);
            int32_t dof_off = target_dof_offset[gpu_id];

            // --- GPU results: non-adjacent linear contributions ---
            for (int32_t d = 0; d < ndofs; ++d) {
               double phi_val = static_cast<double>(phi_hi[dof_off + d])
                              + static_cast<double>(phi_lo[dof_off + d]);
               double phin_val = static_cast<double>(phin_hi[dof_off + d])
                               + static_cast<double>(phin_lo[dof_off + d]);
               if (phi_val != 0.0 || phin_val != 0.0)
                  acc.add(l2g[d], phi_val, phin_val);
            }

            // --- CPU: adjacent linear sources + all non-linear sources ---
            size_t src_idx = 0;
            for (size_t si = 0; si < gd.all_source_ptrs.size(); ++si) {
               auto* src = gd.all_source_ptrs[si];
               if (gd.linear_source_mask[si]) {
                  // Linear source: adjacent only (non-adjacent done on GPU)
                  if (src->isAdjacentTo(tgt))
                     src->fill_direct_entries(tgt, acc);
               } else {
                  // Non-linear source: full processing
                  if (src->fill_direct_entries_nonadj && !src->isAdjacentTo(tgt))
                     src->fill_direct_entries_nonadj(tgt, acc);
                  else
                     src->fill_direct_entries(tgt, acc);
               }
            }

            // --- Extract results: near_indices / weights / RLE ---
            std::sort(acc.touched_list.begin(), acc.touched_list.end());

            const size_t n_touched = acc.touched_list.size();
            tgt->near_indices.clear();
            tgt->near_weights_phi.clear();
            tgt->near_weights_phin.clear();
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

            tgt->near_cell_count = gd.cell_count;
         } // per target
      } // per group
   } // omp parallel

   metal_nearfield_destroy(ctx);
}

#endif /* METAL_NEARFIELD_WRAPPER_HPP */
