#pragma once

// NOTE: This header is included inside `struct BEM_BVP` (see BEM_solveBVP.hpp).
// GMRES + FMM (+ ILU/Jacobi preconditioner) implementation details.

#include <atomic>
#include <deque>
#include <functional>
#include <string>

std::array<double, 2> integrateNearField(const target4FMM& target) const { return target.integrateNearField(cache_phi_val_D_by_index.data(), cache_phin_val_D_by_index.data()); }

inline bool isMidpointDirichletBC(const networkLine* l) const {
  // Keep CORNER classification as mixed boundary.
  // Midpoint Dirichlet-side treatment is applied only to pure Dirichlet edges.
  return l->Dirichlet;
}

// ============================================================================
// compute_Ax_minus_b: BEM行列ベクトル積の中核関数
// ============================================================================
// 役割: FMMを使って A*x - b を計算する
//
// 処理の流れ:
//   1. updateFMM(): 現在の密度値でFMMツリーを更新 (P2M → M2M → M2L → L2L)
//   2. integrateNearField(): 近傍積分（直接計算、Metal GPUまたはCPU）
//   3. integrateFarField(): 遠方積分（FMM多極子展開）
//   4. 結果を結合して A*x - b を返す
//
// 呼び出し元:
//   - return_A_dot_v(): GMRESのKrylov基底構築時 (A*v を計算)
//   - 真の残差計算時: b - A*x = -(A*x - b)
// ============================================================================
V_d compute_Ax_minus_b(bool is_time_derivative = false) {
  const bool do_profile = (profile_fmm_update_time_sum != nullptr && profile_fmm_near_time_sum != nullptr && profile_fmm_far_time_sum != nullptr && profile_fmm_matvec_calls != nullptr);
  auto wall_time = []() -> double {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    using clock = std::chrono::high_resolution_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
#endif
  };

  std::array<double, 6> tmp_elapsed_time{};
  double t_update0 = 0.0;
  if (do_profile)
    t_update0 = wall_time();
  updateFMM(this->B_poles, tmp_elapsed_time);
  if (do_profile) {
    *profile_fmm_update_time_sum += (wall_time() - t_update0);
    *profile_fmm_matvec_calls += 1;
  }
  if (do_profile && profile_fmm_update_step_time_sum != nullptr) {
    for (int i = 0; i < static_cast<int>(tmp_elapsed_time.size()); ++i)
      (*profile_fmm_update_step_time_sum)[i] += tmp_elapsed_time[i];
  }

  V_d Av(matrix_size);
  double near_time_sum = 0.0, far_time_sum = 0.0;

  _Pragma("omp parallel") {
    double near_local = 0.0, far_local = 0.0;

    // Bench (DeepCwind, 3 steps avg; near+far per mat-vec; static=1.00, lower is better):
    // +-------------+--------+
    // | schedule    | ratio  |
    // +-------------+--------+
    // | static      | 1.00   |
    // | dynamic,1   | 1.19   |
    // | dynamic,8   | 1.14   |
    // | guided,8    | 1.09   |
    // +-------------+--------+
    _Pragma("omp for") for (int k = 0; k < points.size(); ++k) {
      auto& a = points[k];
      // Near/Far integrals are functions of the target point a (and the current unknown vector),
      // not of the particular (face,id) row. Compute them once per point and reuse for all IDs.
      bool need_integrals = false;
      for (const auto& [f, i] : a->f2Index) {
        if (!(a->CORNER && isNeumannID_BEM(a, f))) {
          need_integrals = true;
          break;
        }
      }

      std::array<double, 2> near{}, far{};
      if (need_integrals) {
        if (do_profile) {
          const double t0 = wall_time();
          near = integrateNearField(*a);
          near_local += (wall_time() - t0);
          const double t1 = wall_time();
          far = a->integrateFarField();
          far_local += (wall_time() - t1);
        } else {
          near = integrateNearField(*a);
          far = a->integrateFarField();
        }
      }
      const auto [GPhin_point, GnPhi_point] = near + far;

      for (const auto& [f, i] : a->f2Index) {
        if (a->CORNER && isNeumannID_BEM(a, f)) {
          Av[i] = a->phiOnFace_FMM.at(f);
        } else {
          Av[i] = GPhin_point - (GnPhi_point + diag_coeffs[i] * a->phiOnFace_FMM.at(f));
        }
      }
    }

    // Process edge midpoint DOF rows (true quadratic elements)
    if (use_true_quadratic_element) {
      _Pragma("omp for") for (int k = 0; k < static_cast<int>(midpoint_lines.size()); ++k) {
        auto* l = midpoint_lines[k];

        // BIE integration (once per target position, shared by all DOFs of this midpoint)
        auto& mt = static_cast<target4FMM&>(*l);
        std::array<double, 2> near{}, far{};
        if (do_profile) {
          const double t0 = wall_time();
          near = integrateNearField(mt);
          near_local += (wall_time() - t0);
          const double t1 = wall_time();
          far = mt.integrateFarField();
          far_local += (wall_time() - t1);
        } else {
          near = integrateNearField(mt);
          far = mt.integrateFarField();
        }
        const auto [GPhin_mid, GnPhi_mid] = near + far;

        // Write rows for all DOFs of this midpoint (mirrors vertex f2Index pattern)
        for (const auto& [f, i] : l->f2Index) {
          if (i < 0 || i >= static_cast<int>(Av.size()))
            continue;
          if (l->CORNER && isNeumannID_BEM(l, f)) {
            // CORNER constraint row: Av[i] = phi (phi_N = phi_D enforced)
            Av[i] = l->phiOnFace_FMM.at(f);
          } else {
            Av[i] = GPhin_mid - (GnPhi_mid + diag_coeffs[i] * l->phiOnFace_FMM.at(f));
          }
        }
      }
    }

    if (do_profile) {
      _Pragma("omp atomic") near_time_sum += near_local;
      _Pragma("omp atomic") far_time_sum += far_local;
    }
  }

  if (do_profile) {
    *profile_fmm_near_time_sum += near_time_sum;
    *profile_fmm_far_time_sum += far_time_sum;
  }

  return Av;
};

/* -------------------------------------------------------------------------- */

void buildDiagonalPreconditionerForKRing0() {
  TimeWatch watch;
  ilu_preconditioner.reset();
  ilut_preconditioner.reset();
  ilu_faces_to_integrate_by_point.clear();
  ilu_faces_to_integrate_by_midpoint.clear();
  ilu_k_ring_num_cached = 0;
  ilu_topology_hash_cached = 0;

  ilu_diag_inv_preconditioner.assign(matrix_size, 1.0);
  last_A_sparse_nnz = static_cast<double>(matrix_size);
  last_A_sparse_avg_nnz = (matrix_size > 0) ? 1.0 : 0.0;

  _Pragma("omp parallel for") for (int k = 0; k < static_cast<int>(points.size()); ++k) {
    auto* p = points[k];
    for (const auto& [f, i] : p->f2Index) {
      double diag = 1.0;

      // Corner constraint row: A = I on the Neumann ID
      if (p->CORNER && isNeumannID_BEM(p, f)) {
        diag = 1.0;
      } else if (isNeumannID_BEM(p, f)) {
        // Use only the diagonal jump/free-term coefficient (no face integration).
        diag = -diag_coeffs[i];
      } else {
        // Dirichlet IDs: unknown is phin; without faces this term is not defined -> identity.
        diag = 1.0;
      }

      if (!std::isfinite(diag) || std::abs(diag) < 1e-15)
        diag = 1.0;
      ilu_diag_inv_preconditioner[i] = 1.0 / diag;
    }
  }
  last_ilu_build_time = watch()[0];
}

static inline void hash_combine(std::size_t& seed, std::size_t v) noexcept { seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2); }

std::size_t computeILUTopologySignature() const {
  std::size_t h = 0;
  for (const auto* water : WATERS)
    hash_combine(h, water->topologySignature());
  return h;
}

std::size_t computeILUGeometrySignature() const {
  std::size_t h = 0;
  for (const auto* water : WATERS)
    hash_combine(h, water->geometrySignature());
  return h;
}

void buildSparseMatrixForILU() {
  TimeWatch watch_build;
  const bool want_schwarz = (preconditioner_type == "SCHWARZ");
  const bool want_ilut = (preconditioner_type == "ILUT");
  const bool want_milu = (preconditioner_type == "MILU");
  const char* precond_name = want_schwarz ? "Schwarz" : (want_ilut ? "ILUT" : (want_milu ? "MILU" : "ILU"));
  const bool use_k_ring = (ilu_neighborhood_type == "K-RING");
  const int k_ring_num = std::clamp(ilu_kring_num, 0, 20);
  if (use_k_ring && k_ring_num == 0) {
    std::cout << "Building diagonal (Jacobi) preconditioner for " << precond_name << " k-ring num=0..." << std::endl;
    buildDiagonalPreconditionerForKRing0();
    ilu_preconditioner.reset();
    ilut_preconditioner.reset();
    schwarz_preconditioner.reset();
    const auto elapsed = watch_build();
    last_ilu_build_time = elapsed[0];
    std::cout << Magenta << "[BEM] " << Cyan << "build " << precond_name << " preconditioner" << Green << " elapsed=" << elapsed << colorReset << std::endl;
    ilu_geometry_hash_cached = computeILUGeometrySignature();
    preconditioner_type_built = preconditioner_type;
    ilu_neighborhood_type_built = ilu_neighborhood_type;
    ilu_k_ring_num_built = k_ring_num;
    ilu_matrix_size_built = matrix_size;
    ilut_drop_tol_built = ilut_drop_tol;
    ilut_max_entries_per_row_built = ilut_max_entries_per_row;
    ilut_pivot_min_built = ilut_pivot_min;
    schwarz_core_k_built = schwarz_core_k;
    schwarz_overlap_k_built = schwarz_overlap_k;
    schwarz_max_core_size_built = schwarz_max_core_size;
    schwarz_max_block_size_built = schwarz_max_block_size;
    schwarz_pivot_min_built = schwarz_pivot_min;
    schwarz_diag_shift_built = schwarz_diag_shift;
    return;
  }

  const int k_ring_depth = k_ring_num;
  std::cout << "Building sparse matrix for " << precond_name << " preconditioner (" << (use_k_ring ? ("K-RING num=" + std::to_string(k_ring_num)) : "BUCKETS") << ")..." << std::endl;

  // 1) Prepare CRS rows
  std::vector<CRS> crs_nodes(matrix_size);
  std::vector<CRS*> crs_ptrs(matrix_size);
  for (int i = 0; i < matrix_size; ++i) {
    crs_nodes[i].setIndexCRS(i);
    crs_ptrs[i] = &crs_nodes[i];
  }

  // 2) Compute global scale and (optionally) build a spatial bucket of boundary faces
  std::size_t total_faces = 0;
  CoordinateBounds bounds_all;
  double global_scale = 0.;
  for (const auto* water : WATERS) {
    total_faces += water->getBoundaryFaces().size();
    global_scale = std::max(global_scale, water->getScale());
    const auto& b = water->bounds;
    std::get<0>(bounds_all.bounds[0]) = std::min(std::get<0>(bounds_all.bounds[0]), std::get<0>(b[0]));
    std::get<1>(bounds_all.bounds[0]) = std::max(std::get<1>(bounds_all.bounds[0]), std::get<1>(b[0]));
    std::get<0>(bounds_all.bounds[1]) = std::min(std::get<0>(bounds_all.bounds[1]), std::get<0>(b[1]));
    std::get<1>(bounds_all.bounds[1]) = std::max(std::get<1>(bounds_all.bounds[1]), std::get<1>(b[1]));
    std::get<0>(bounds_all.bounds[2]) = std::min(std::get<0>(bounds_all.bounds[2]), std::get<0>(b[2]));
    std::get<1>(bounds_all.bounds[2]) = std::max(std::get<1>(bounds_all.bounds[2]), std::get<1>(b[2]));
  }

  std::unique_ptr<Buckets<networkFace*>> face_bucket;
  double search_range = 0.;
  if (!use_k_ring) {
    // Use a uniform (cube) bounding box and choose dL so that the number of cells stays moderate.
    const auto uniform_bounds = bounds_all.getUniformBounds(1.05);
    const double L = std::max({std::get<1>(uniform_bounds[0]) - std::get<0>(uniform_bounds[0]), std::get<1>(uniform_bounds[1]) - std::get<0>(uniform_bounds[1]), std::get<1>(uniform_bounds[2]) - std::get<0>(uniform_bounds[2])});
    const int cells_per_axis = std::clamp<int>(static_cast<int>(std::llround(std::cbrt(static_cast<double>(std::max<std::size_t>(total_faces, 1))))), 8, 64);
    const double dL = (L > 0.) ? (L / static_cast<double>(cells_per_axis)) : 1.0;

    face_bucket = std::make_unique<Buckets<networkFace*>>(uniform_bounds, dL);
    for (const auto* water : WATERS) {
      for (auto* f : water->getBoundaryFaces()) {
        // NOTE: centroid-only insertion can miss large faces (point may be close to an edge
        // while far from centroid). Insert a few representative points on the face.
        face_bucket->add(f->centroid, f);
        const auto& pts = f->getPoints();
        const auto& x0 = std::get<0>(pts)->X;
        const auto& x1 = std::get<1>(pts)->X;
        const auto& x2 = std::get<2>(pts)->X;
        face_bucket->add(x0, f);
        face_bucket->add(x1, f);
        face_bucket->add(x2, f);
        face_bucket->add((x0 + x1) / 2., f);
        face_bucket->add((x1 + x2) / 2., f);
        face_bucket->add((x2 + x0) / 2., f);
      }
    }
    face_bucket->setVector(); // make neighbor queries faster (read-only afterwards)
    search_range = 0. * dL;   // "幅" (cube range) for collecting near faces
  }

  const double near_switch = (global_scale > 0.) ? (global_scale / 10.) : search_range;

  // 2.5) If using k-ring, precompute candidate faces once per target point (cached on the BVP)
  if (use_k_ring) {
    const std::size_t topology_sig = computeILUTopologySignature();
    const bool cache_hit = (ilu_faces_to_integrate_by_point.size() == points.size() && ilu_k_ring_num_cached == k_ring_num && ilu_topology_hash_cached == topology_sig);
    if (!cache_hit) {
      std::cout << Red << "Precomputing ILU k-ring neighborhoods (num=" << k_ring_num << ")..." << colorReset << std::endl;

      std::unordered_set<networkFace*> boundary_faces_set;
      boundary_faces_set.reserve(std::max<std::size_t>(total_faces * 2, 1024));
      for (const auto* water : WATERS)
        for (auto* f : water->getBoundaryFaces())
          boundary_faces_set.emplace(f);

      ilu_faces_to_integrate_by_point.assign(points.size(), {});

      const bool k_ring_parallel = []() {
        if (const char* env = std::getenv("BEM_KRING_PARALLEL"))
          return (std::string(env) != "0");
        return false; // safer default (avoid thread-unsafe traversal)
      }();

      auto build_k_ring_for_point = [&](int k, std::vector<networkFace*>& local_faces) {
        auto* origin = points[k];
        std::unordered_set<networkFace*> start_faces;
        {
          const auto start_faces_vec = origin->getBoundaryFaces();
          if (!start_faces_vec.empty()) {
            start_faces.insert(start_faces_vec.begin(), start_faces_vec.end());
          } else {
            const auto& fallback_faces = origin->getFaces();
            start_faces.insert(fallback_faces.begin(), fallback_faces.end());
          }
        }

        std::unordered_set<networkFace*> faces_uo;
        if (k_ring_depth == 0) {
          faces_uo = start_faces;
        } else {
          faces_uo = BFS_Flattened(start_faces, k_ring_depth, WATERS);
        }

        local_faces.clear();
        local_faces.reserve(faces_uo.size());
        for (auto* f : faces_uo)
          if (boundary_faces_set.contains(f))
            local_faces.emplace_back(f);

        std::sort(local_faces.begin(), local_faces.end());
        local_faces.erase(std::unique(local_faces.begin(), local_faces.end()), local_faces.end());

        ilu_faces_to_integrate_by_point[k] = local_faces;
      };

      if (k_ring_parallel) {
        _Pragma("omp parallel") {
          thread_local std::vector<networkFace*> local_faces;
          local_faces.clear();
          _Pragma("omp for") for (int k = 0; k < static_cast<int>(points.size()); ++k) { build_k_ring_for_point(k, local_faces); }
        }
      } else {
        thread_local std::vector<networkFace*> local_faces;
        local_faces.clear();
        for (int k = 0; k < static_cast<int>(points.size()); ++k) {
          build_k_ring_for_point(k, local_faces);
        }
      }

      // Build k-ring neighborhoods for midpoint targets (true_quadratic)
      if (use_true_quadratic_element && !midpoint_lines.empty()) {
        ilu_faces_to_integrate_by_midpoint.assign(midpoint_lines.size(), {});
        std::vector<networkFace*> mid_local_faces;
        for (std::size_t mk = 0; mk < midpoint_lines.size(); ++mk) {
          auto [pA, pB] = midpoint_lines[mk]->getPoints();
          std::unordered_set<networkFace*> start_faces;
          for (auto* f : pA->getBoundaryFaces())
            start_faces.insert(f);
          for (auto* f : pB->getBoundaryFaces())
            start_faces.insert(f);

          std::unordered_set<networkFace*> faces_uo;
          if (k_ring_depth == 0)
            faces_uo = start_faces;
          else
            faces_uo = BFS_Flattened(start_faces, k_ring_depth, WATERS);

          mid_local_faces.clear();
          mid_local_faces.reserve(faces_uo.size());
          for (auto* f : faces_uo)
            if (boundary_faces_set.contains(f))
              mid_local_faces.emplace_back(f);
          std::sort(mid_local_faces.begin(), mid_local_faces.end());
          mid_local_faces.erase(std::unique(mid_local_faces.begin(), mid_local_faces.end()), mid_local_faces.end());
          ilu_faces_to_integrate_by_midpoint[mk] = mid_local_faces;
        }
      } else {
        ilu_faces_to_integrate_by_midpoint.clear();
      }

      ilu_k_ring_num_cached = k_ring_num;
      ilu_topology_hash_cached = topology_sig;
    } else {
      std::cout << Green << "ILU k-ring neighborhoods cache hit." << colorReset << std::endl;
    }
  }

  std::cout << Magenta << "[BEM] " << Cyan << "Assembling sparse matrix for " << precond_name << " elapsed=" << watch_build() << colorReset << std::endl;

  // 3) For each target point, collect nearby faces and assemble sparse rows by direct integration
  // NOTE: This loop is typically dominated by allocations & hashing. Reuse thread-local buffers.
  _Pragma("omp parallel") {
    thread_local std::vector<networkFace*> faces_to_integrate;
    faces_to_integrate.clear();
    thread_local std::unordered_map<int, double> coeffs;
    coeffs.clear();

    _Pragma("omp for") for (int k = 0; k < static_cast<int>(points.size()); ++k) {
      auto* origin = points[k];

      const std::vector<networkFace*>* faces_view = nullptr;

      faces_to_integrate.clear();
      coeffs.clear();

      if (use_k_ring) {
        faces_view = &ilu_faces_to_integrate_by_point[k];
      } else {
        // Collect candidate faces once per target point (same for all IDs of the point).
        face_bucket->apply(origin->X, search_range, [&](networkFace* f) { faces_to_integrate.push_back(f); });
        for (auto* f : origin->getBoundaryFaces())
          faces_to_integrate.push_back(f); // guarantee self/adjacent faces included
        // De-duplicate
        std::sort(faces_to_integrate.begin(), faces_to_integrate.end());
        faces_to_integrate.erase(std::unique(faces_to_integrate.begin(), faces_to_integrate.end()), faces_to_integrate.end());
        faces_view = &faces_to_integrate;
      }

      auto accumulate_coeff = [&](networkPoint* p, networkFace* integ_f, const double ig, const double ign) {
        const auto key_f = std::get<1>(pf2ID(p, integ_f));
        const int col_idx = p->f2Index.at(key_f);
        const bool col_is_neumann = isNeumannID_BEM(p, key_f);
        const double val = col_is_neumann ? (-ign) : ig; // G*phin - Gn*phi
        if (val != 0.)
          coeffs[col_idx] += val;
      };

      for (auto* integ_f : *faces_view) {
        if (integ_f == nullptr)
          continue;

        const auto [p0, p1, p2] = integ_f->getPoints(origin);
        const auto* closest_p_to_origin = p0;
        const double dist = Norm((p0->X + p1->X + p2->X) / 3. - origin->X);
        const int how_far = (dist < near_switch) ? 1 : 0;

        if (integ_f->isLinearElement) {
          Tdd ig_ign0{0., 0.}, ig_ign1{0., 0.}, ig_ign2{0., 0.};
          double nr = 0., ig = 0., ign = 0., ww_nr = 0.;
          Tddd R{0., 0., 0.};

          if (p0 == origin) {
            // on-vertex singular case (linear element): use precomputed "near" integration and ignore ign term
            for (const auto& [t0t1, ww, shape3, X, cross, J_det] : integ_f->map_Point_LinearIntegrationInfo_vector[1].at(const_cast<networkPoint*>(closest_p_to_origin))) {
              ig = J_det * (ww / (nr = Norm(R = (X - origin->X))));
              std::get<0>(ig_ign0) += ig * std::get<0>(shape3);
              std::get<0>(ig_ign1) += ig * std::get<1>(shape3);
              std::get<0>(ig_ign2) += ig * std::get<2>(shape3);
            }
          } else {
            for (const auto& [t0t1, ww, shape3, X, cross, J_det] : integ_f->map_Point_LinearIntegrationInfo_vector[how_far].at(const_cast<networkPoint*>(closest_p_to_origin))) {
              ig = J_det * (ww_nr = ww / (nr = Norm(R = (X - origin->X))));
              ign = Dot(R, cross) * ww_nr / (nr * nr);
              std::get<0>(ig_ign0) += ig * std::get<0>(shape3);
              std::get<1>(ig_ign0) -= ign * std::get<0>(shape3);
              std::get<0>(ig_ign1) += ig * std::get<1>(shape3);
              std::get<1>(ig_ign1) -= ign * std::get<1>(shape3);
              std::get<0>(ig_ign2) += ig * std::get<2>(shape3);
              std::get<1>(ig_ign2) -= ign * std::get<2>(shape3);
            }
          }

          accumulate_coeff(p0, integ_f, std::get<0>(ig_ign0), std::get<1>(ig_ign0));
          accumulate_coeff(p1, integ_f, std::get<0>(ig_ign1), std::get<1>(ig_ign1));
          accumulate_coeff(p2, integ_f, std::get<0>(ig_ign2), std::get<1>(ig_ign2));
        } else if (integ_f->isPseudoQuadraticElement) {
          // Use face-precomputed mapping for pseudo-quadratic elements
          auto key_ig_ign = integ_f->map_Point_BEM_IGIGn_info_init.at(const_cast<networkPoint*>(closest_p_to_origin));
          double nr = 0., ig = 0., ign = 0., ww_nr = 0., tmp = 0.;
          Tddd R{0., 0., 0.};

          for (const auto& [t0t1, ww, Nc_N0_N1_N2, X, cross, J_det] : integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector[how_far].at(const_cast<networkPoint*>(closest_p_to_origin))) {
            ig = J_det * (ww_nr = ww / (nr = Norm(R = (X - origin->X))));
            ign = Dot(R, cross) * ww_nr / (nr * nr);
            for (int i = 0; i < 6; ++i) {
              tmp = std::get<0>(Nc_N0_N1_N2)[i];
              std::get<2>(key_ig_ign[i]) += ig * tmp;
              // Rigid-mode technique: for the diagonal (self) term, drop the singular Gn contribution
              // associated with the origin point, and let the jump/free-term coefficient handle it.
              if (std::get<0>(key_ig_ign[i]) != origin)
                std::get<3>(key_ig_ign[i]) -= ign * tmp;

              tmp = std::get<1>(Nc_N0_N1_N2)[i];
              std::get<2>(key_ig_ign[i + 6]) += ig * tmp;
              if (std::get<0>(key_ig_ign[i + 6]) != origin)
                std::get<3>(key_ig_ign[i + 6]) -= ign * tmp;

              tmp = std::get<2>(Nc_N0_N1_N2)[i];
              std::get<2>(key_ig_ign[i + 12]) += ig * tmp;
              if (std::get<0>(key_ig_ign[i + 12]) != origin)
                std::get<3>(key_ig_ign[i + 12]) -= ign * tmp;

              tmp = std::get<3>(Nc_N0_N1_N2)[i];
              std::get<2>(key_ig_ign[i + 18]) += ig * tmp;
              if (std::get<0>(key_ig_ign[i + 18]) != origin)
                std::get<3>(key_ig_ign[i + 18]) -= ign * tmp;
            }
          }

          for (const auto& [p, f, ig_v, ign_v] : key_ig_ign)
            accumulate_coeff(p, f, ig_v, ign_v);
        } else if (integ_f->isTrueQuadraticElement) {
          // True quadratic: raw quadrature (no precomputed tables for this element type)
          const auto [fp0, fp1, fp2] = integ_f->getPoints();
          const auto& [fl0, fl1, fl2] = integ_f->Lines;
          auto dof_idx = getQuadDOFIndices(integ_f);
          const T6Tddd X6_tq = {fp0->X, fp1->X, fp2->X, fl0->X_mid, fl1->X_mid, fl2->X_mid};
          constexpr std::array<bool, 3> all_true_tq{true, true, true};
          std::array<std::array<double, 2>, 6> WGN{}; // [G_weighted, Gn_weighted] per DOF

          int vertex_local = -1;
          if (origin == fp0)
            vertex_local = 0;
          else if (origin == fp1)
            vertex_local = 1;
          else if (origin == fp2)
            vertex_local = 2;
          const bool is_vert_on_elem = (vertex_local >= 0);

          if (is_vert_on_elem) {
            // Duffy transform: permute so origin is at vertex 0
            std::array<int, 3> vperm = {0, 1, 2};
            if (fp1 == origin) {
              vperm = {1, 2, 0};
            } else if (fp2 == origin) {
              vperm = {2, 0, 1};
            }

            for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
              const auto bary_loc = ModTriShape<3>(t0, t1);
              std::array<double, 3> bary{};
              bary[vperm[0]] = bary_loc[0];
              bary[vperm[1]] = bary_loc[1];
              bary[vperm[2]] = bary_loc[2];

              auto N6_geo = TriShape<6>(bary[0], bary[1], all_true_tq);
              auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true_tq);
              auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true_tq);
              auto X_q = Dot(N6_geo, X6_tq);
              auto cross_q = Cross(Dot(dN_dt0, X6_tq), Dot(dN_dt1, X6_tq));
              auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);
              auto R = X_q - origin->X;
              double nr = Norm(R);
              double nr_inv = 1.0 / nr;
              double w_nr = ww * (1. - t0) * nr_inv;
              double WG = Norm(cross_q) * w_nr;
              double WGn = -Dot(R * nr_inv, cross_q) * w_nr * nr_inv;
              for (int kk = 0; kk < 6; ++kk) {
                WGN[kk][0] += WG * N6[kk];
                WGN[kk][1] += (kk == vertex_local) ? 0. : WGn * N6[kk]; // rigid mode: skip Gn for origin
              }
            }
          } else {
            // Non-adjacent: standard GW5xGW5 quadrature
            for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
              auto bary = ModTriShape<3>(t0, t1);
              auto N6_geo = TriShape<6>(bary[0], bary[1], all_true_tq);
              auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true_tq);
              auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true_tq);
              auto X_q = Dot(N6_geo, X6_tq);
              auto cross_q = Cross(Dot(dN_dt0, X6_tq), Dot(dN_dt1, X6_tq));
              auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);
              auto R = X_q - origin->X;
              double nr = Norm(R);
              double nr_inv = 1.0 / nr;
              double w_nr = ww * (1. - t0) * nr_inv;
              double WG = Norm(cross_q) * w_nr;
              double WGn = -Dot(R * nr_inv, cross_q) * w_nr * nr_inv;
              for (int kk = 0; kk < 6; ++kk) {
                WGN[kk][0] += WG * N6[kk];
                WGN[kk][1] += WGn * N6[kk];
              }
            }
          }

          // Accumulate coefficients: vertex DOFs via accumulate_coeff, midpoint DOFs directly
          const std::array<networkPoint*, 3> face_pts_tq = {fp0, fp1, fp2};
          const std::array<networkLine*, 3> face_lines_tq = {fl0, fl1, fl2};
          for (int kk = 0; kk < 6; ++kk) {
            if (std::abs(WGN[kk][0]) < 1e-30 && std::abs(WGN[kk][1]) < 1e-30)
              continue;
            if (kk < 3) {
              accumulate_coeff(face_pts_tq[kk], integ_f, WGN[kk][0], WGN[kk][1]);
            } else {
              int col_idx = dof_idx[kk];
              if (col_idx < 0)
                continue;
              auto lf_key = std::get<1>(lf2ID(face_lines_tq[kk - 3], integ_f));
              bool col_neumann = isNeumannID_BEM(face_lines_tq[kk - 3], lf_key);
              double val = col_neumann ? (-WGN[kk][1]) : WGN[kk][0];
              if (val != 0.)
                coeffs[col_idx] += val;
            }
          }
        }
      }

      // 4) Write rows for all IDs at the target point
      for (const auto& [target_f, row_idx] : origin->f2Index) {
        auto& row = crs_nodes[row_idx];

        // Corner constraint row: A = I on the Neumann ID
        if (origin->CORNER && isNeumannID_BEM(origin, target_f)) {
          row.clearColumnValue();
          row.set(crs_ptrs[row_idx], 1.0);
          continue;
        }

        for (const auto& [col_idx, v] : coeffs)
          row.increment(crs_ptrs[col_idx], v);

        // Diagonal jump term alpha (only multiplies unknown phi)
        if (isNeumannID_BEM(origin, target_f)) {
          const double alpha = diag_coeffs[row_idx];
          row.increment(crs_ptrs[row_idx], -alpha);
        }

        // ILU0_CRS requires an explicit diagonal entry.
        if (!row.contains(crs_ptrs[row_idx]))
          row.set(crs_ptrs[row_idx], 1.0);
      }
    }
  }

  // True quadratic: full BIE integration for midpoint DOF rows
  if (use_true_quadratic_element && !midpoint_lines.empty()) {
    _Pragma("omp parallel") {
      thread_local std::vector<networkFace*> mid_faces_to_integrate;
      mid_faces_to_integrate.clear();
      thread_local std::unordered_map<int, double> mid_coeffs;
      mid_coeffs.clear();

      _Pragma("omp for") for (int mk = 0; mk < static_cast<int>(midpoint_lines.size()); ++mk) {
        auto* l = midpoint_lines[mk];

        auto [pA, pB] = l->getPoints();
        Tddd target_pos = l->X_mid;
        if (!std::isfinite(target_pos[0]) || !std::isfinite(target_pos[1]) || !std::isfinite(target_pos[2]))
          target_pos = 0.5 * (pA->Xtarget + pB->Xtarget);

        mid_coeffs.clear();

        // Get faces to integrate
        const std::vector<networkFace*>* mid_faces_view = nullptr;
        if (use_k_ring) {
          mid_faces_view = &ilu_faces_to_integrate_by_midpoint[mk];
        } else {
          mid_faces_to_integrate.clear();
          face_bucket->apply(target_pos, search_range, [&](networkFace* f) { mid_faces_to_integrate.push_back(f); });
          for (auto* f : pA->getBoundaryFaces())
            mid_faces_to_integrate.push_back(f);
          for (auto* f : pB->getBoundaryFaces())
            mid_faces_to_integrate.push_back(f);
          std::sort(mid_faces_to_integrate.begin(), mid_faces_to_integrate.end());
          mid_faces_to_integrate.erase(std::unique(mid_faces_to_integrate.begin(), mid_faces_to_integrate.end()), mid_faces_to_integrate.end());
          mid_faces_view = &mid_faces_to_integrate;
        }

        for (auto* integ_f : *mid_faces_view) {
          if (integ_f == nullptr || !integ_f->isTrueQuadraticElement)
            continue;

          const auto [fp0, fp1, fp2] = integ_f->getPoints();
          const auto& [fl0, fl1, fl2] = integ_f->Lines;
          auto dof_idx = getQuadDOFIndices(integ_f);
          const T6Tddd X6_m = {fp0->X, fp1->X, fp2->X, fl0->X_mid, fl1->X_mid, fl2->X_mid};
          constexpr std::array<bool, 3> all_true_m{true, true, true};
          std::array<std::array<double, 2>, 6> WGN_m{};

          // Check if this midpoint target belongs to one of the face edges.
          int on_edge = -1;
          if (l == fl0)
            on_edge = 0;
          else if (l == fl1)
            on_edge = 1;
          else if (l == fl2)
            on_edge = 2;

          if (on_edge >= 0) {
            // ===== MIDPOINT ON ELEMENT EDGE: subdivision + Duffy =====
            int va = on_edge, vb = (on_edge + 1) % 3, vc = (on_edge + 2) % 3;
            int mid_local = 3 + on_edge;

            for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
              double omt0 = 1. - t0;
              double s1 = t1 * omt0;
              double s2 = omt0 * (1. - t1);

              // Sub-triangle 1: (M, vb, vc)
              {
                std::array<double, 3> bary{};
                bary[va] = 0.5 * t0;
                bary[vb] = 0.5 * t0 + s1;
                bary[vc] = s2;
                auto N6_geo = TriShape<6>(bary[0], bary[1], all_true_m);
                auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true_m);
                auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true_m);
                auto X_q = Dot(N6_geo, X6_m);
                auto cross_q = Cross(Dot(dN_dt0, X6_m), Dot(dN_dt1, X6_m));
                auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);
                auto R = X_q - target_pos;
                double nr = Norm(R);
                double nr_inv = 1.0 / nr;
                double w_nr = ww * omt0 * 0.5 * nr_inv;
                double WG = Norm(cross_q) * w_nr;
                double WGn = -Dot(R * nr_inv, cross_q) * w_nr * nr_inv;
                for (int kk = 0; kk < 6; ++kk) {
                  WGN_m[kk][0] += WG * N6[kk];
                  WGN_m[kk][1] += (kk == mid_local) ? 0. : WGn * N6[kk]; // rigid mode
                }
              }

              // Sub-triangle 2: (M, vc, va)
              {
                std::array<double, 3> bary{};
                bary[va] = 0.5 * t0 + s2;
                bary[vb] = 0.5 * t0;
                bary[vc] = s1;
                auto N6_geo = TriShape<6>(bary[0], bary[1], all_true_m);
                auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true_m);
                auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true_m);
                auto X_q = Dot(N6_geo, X6_m);
                auto cross_q = Cross(Dot(dN_dt0, X6_m), Dot(dN_dt1, X6_m));
                auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);
                auto R = X_q - target_pos;
                double nr = Norm(R);
                double nr_inv = 1.0 / nr;
                double w_nr = ww * omt0 * 0.5 * nr_inv;
                double WG = Norm(cross_q) * w_nr;
                double WGn = -Dot(R * nr_inv, cross_q) * w_nr * nr_inv;
                for (int kk = 0; kk < 6; ++kk) {
                  WGN_m[kk][0] += WG * N6[kk];
                  WGN_m[kk][1] += (kk == mid_local) ? 0. : WGn * N6[kk]; // rigid mode
                }
              }
            }
          } else {
            // ===== NON-ADJACENT: standard GW5xGW5 =====
            for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
              auto bary = ModTriShape<3>(t0, t1);
              auto N6_geo = TriShape<6>(bary[0], bary[1], all_true_m);
              auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true_m);
              auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true_m);
              auto X_q = Dot(N6_geo, X6_m);
              auto cross_q = Cross(Dot(dN_dt0, X6_m), Dot(dN_dt1, X6_m));
              auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);
              auto R = X_q - target_pos;
              double nr = Norm(R);
              double nr_inv = 1.0 / nr;
              double w_nr = ww * (1. - t0) * nr_inv;
              double WG = Norm(cross_q) * w_nr;
              double WGn = -Dot(R * nr_inv, cross_q) * w_nr * nr_inv;
              for (int kk = 0; kk < 6; ++kk) {
                WGN_m[kk][0] += WG * N6[kk];
                WGN_m[kk][1] += WGn * N6[kk];
              }
            }
          }

          // Accumulate coefficients
          const std::array<networkPoint*, 3> face_pts_m = {fp0, fp1, fp2};
          const std::array<networkLine*, 3> face_lines_m = {fl0, fl1, fl2};
          for (int kk = 0; kk < 6; ++kk) {
            if (std::abs(WGN_m[kk][0]) < 1e-30 && std::abs(WGN_m[kk][1]) < 1e-30)
              continue;
            int col_idx;
            bool col_neumann;
            if (kk < 3) {
              col_idx = pf2Index(face_pts_m[kk], integ_f);
              const auto key_f = std::get<1>(pf2ID(face_pts_m[kk], integ_f));
              col_neumann = isNeumannID_BEM(face_pts_m[kk], key_f);
            } else {
              col_idx = dof_idx[kk];
              if (col_idx < 0)
                continue;
              auto lf_key = std::get<1>(lf2ID(face_lines_m[kk - 3], integ_f));
              col_neumann = isNeumannID_BEM(face_lines_m[kk - 3], lf_key);
            }
            double val = col_neumann ? (-WGN_m[kk][1]) : WGN_m[kk][0];
            if (val != 0.)
              mid_coeffs[col_idx] += val;
          }
        }

        // Write rows for all DOFs of this midpoint (mirrors vertex f2Index pattern)
        for (const auto& [target_f, row_idx] : l->f2Index) {
          auto& row = crs_nodes[row_idx];

          // CORNER constraint row: A = I on the Neumann DOF
          if (l->CORNER && isNeumannID_BEM(l, target_f)) {
            row.clearColumnValue();
            row.set(crs_ptrs[row_idx], 1.0);
            continue;
          }

          for (const auto& [col_idx, v] : mid_coeffs)
            row.increment(crs_ptrs[col_idx], v);

          // Diagonal jump term for Neumann midpoint
          if (isNeumannID_BEM(l, target_f)) {
            const double alpha = diag_coeffs[row_idx];
            row.increment(crs_ptrs[row_idx], -alpha);
          }

          // Ensure diagonal entry exists
          if (!row.contains(crs_ptrs[row_idx]))
            row.set(crs_ptrs[row_idx], 1.0);
        }
      }
    }
  }

  // Sparsity stats (for diagnostics): nnz(A_sparse) and average nnz per row
  {
    std::size_t nnz = 0;
    for (const auto& row : crs_nodes)
      nnz += row.column_value.size();
    last_A_sparse_nnz = static_cast<double>(nnz);
    last_A_sparse_avg_nnz = (matrix_size > 0) ? (static_cast<double>(nnz) / static_cast<double>(matrix_size)) : 0.0;
  }

  if (const char* env = std::getenv("BEM_ILU_MATRIX_DEBUG"); env && std::string(env) != "0") {
    std::vector<char> is_midpoint_row(static_cast<std::size_t>(matrix_size), 0);
    if (use_true_quadratic_element) {
      for (auto* l : midpoint_lines) {
        for (const auto& [f, i] : l->f2Index) {
          if (i >= 0 && i < matrix_size)
            is_midpoint_row[static_cast<std::size_t>(i)] = 1;
        }
      }
    }

    double min_abs_diag = std::numeric_limits<double>::max();
    double max_abs_diag = 0.0;
    int min_diag_row = -1;
    std::size_t tiny_diag_count = 0;
    std::size_t tiny_diag_mid_count = 0;
    std::size_t tiny_diag_vertex_count = 0;
    constexpr double tiny_thr = 1e-12;

    for (int i = 0; i < matrix_size; ++i) {
      const auto& row = crs_nodes[static_cast<std::size_t>(i)];
      const auto it = row.column_value.find(crs_ptrs[static_cast<std::size_t>(i)]);
      const double diag = (it == row.column_value.end()) ? 0.0 : it->second;
      const double ad = std::abs(diag);
      if (ad < min_abs_diag) {
        min_abs_diag = ad;
        min_diag_row = i;
      }
      max_abs_diag = std::max(max_abs_diag, ad);
      if (ad < tiny_thr) {
        ++tiny_diag_count;
        if (is_midpoint_row[static_cast<std::size_t>(i)])
          ++tiny_diag_mid_count;
        else
          ++tiny_diag_vertex_count;
      }
    }

    std::cout << Magenta << "[ILU:matrix] " << Cyan
              << "min|diag|=" << Yellow << min_abs_diag << Cyan
              << " (row=" << Yellow << min_diag_row << Cyan
              << ", kind=" << Yellow << (is_midpoint_row[static_cast<std::size_t>(std::max(0, min_diag_row))] ? "midpoint" : "vertex") << Cyan << ")"
              << "  max|diag|=" << Yellow << max_abs_diag << Cyan
              << "  tiny(|diag|<1e-12)=" << Yellow << tiny_diag_count << Cyan
              << " [vertex=" << Yellow << tiny_diag_vertex_count << Cyan
              << ", midpoint=" << Yellow << tiny_diag_mid_count << Cyan << "]"
              << colorReset << std::endl;
  }

  std::cout << Magenta << "[BEM] " << Cyan << "build sparse matrix (CRS)" << Green << " elapsed=" << watch_build() << ", nnz=" << static_cast<std::size_t>(last_A_sparse_nnz) << ", avg_nnz/row=" << last_A_sparse_avg_nnz << colorReset << std::endl;

  // 5) Create preconditioner (ILU(0) / ILUT / Schwarz)
  ilu_diag_inv_preconditioner.clear();
  if (preconditioner_type == "SCHWARZ") {
    ilu_preconditioner.reset();
    ilut_preconditioner.reset();

    const int core_k = std::clamp(schwarz_core_k, 0, 8);
    const int overlap_k = std::clamp(schwarz_overlap_k, 0, 8);
    const int max_core = std::clamp(schwarz_max_core_size, 1, std::max(1, matrix_size));
    const int max_block = std::clamp(schwarz_max_block_size, max_core, std::max(1, matrix_size));

    // Build adjacency from the sparse matrix pattern (row -> column indices), excluding the diagonal.
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(matrix_size));
    for (int i = 0; i < matrix_size; ++i) {
      const auto& row = crs_nodes[i];
      auto& a = adj[static_cast<std::size_t>(i)];
      a.reserve(row.column_value.size());
      for (const auto& [col_ptr, v] : row.column_value) {
        (void)v;
        if (!col_ptr)
          continue;
        const int j = static_cast<int>(col_ptr->getIndexCRS());
        if (j != i)
          a.push_back(j);
      }
    }

    // Greedy partition into non-overlapping cores; then expand overlap.
    std::vector<char> assigned(static_cast<std::size_t>(matrix_size), 0);
    std::vector<int> mark(static_cast<std::size_t>(matrix_size), -1);
    int stamp = 0;
    std::deque<std::pair<int, int>> q;

    std::vector<std::vector<int>> blocks;
    blocks.reserve(static_cast<std::size_t>((matrix_size + max_core - 1) / max_core));

    for (int seed = 0; seed < matrix_size; ++seed) {
      if (assigned[static_cast<std::size_t>(seed)])
        continue;

      std::vector<int> core;
      core.reserve(static_cast<std::size_t>(max_core));

      // Core BFS on unassigned nodes only.
      ++stamp;
      q.clear();
      q.emplace_back(seed, 0);
      mark[static_cast<std::size_t>(seed)] = stamp;

      while (!q.empty() && static_cast<int>(core.size()) < max_core) {
        const auto [u, d] = q.front();
        q.pop_front();

        if (!assigned[static_cast<std::size_t>(u)]) {
          assigned[static_cast<std::size_t>(u)] = 1;
          core.push_back(u);
          if (static_cast<int>(core.size()) >= max_core)
            break;
        }

        if (d >= core_k)
          continue;
        for (const int v : adj[static_cast<std::size_t>(u)]) {
          if (v < 0 || v >= matrix_size)
            continue;
          auto& mv = mark[static_cast<std::size_t>(v)];
          if (mv == stamp)
            continue;
          mv = stamp;
          if (assigned[static_cast<std::size_t>(v)])
            continue;
          q.emplace_back(v, d + 1);
        }
      }

      if (core.empty()) {
        assigned[static_cast<std::size_t>(seed)] = 1;
        core.push_back(seed);
      }

      // Overlap expansion starting from the core (can include already-assigned indices).
      std::vector<int> block = core;
      block.reserve(static_cast<std::size_t>(max_block));

      ++stamp;
      q.clear();
      for (const int u : core) {
        mark[static_cast<std::size_t>(u)] = stamp;
        q.emplace_back(u, 0);
      }

      while (!q.empty() && static_cast<int>(block.size()) < max_block) {
        const auto [u, d] = q.front();
        q.pop_front();
        if (d >= overlap_k)
          continue;

        for (const int v : adj[static_cast<std::size_t>(u)]) {
          if (v < 0 || v >= matrix_size)
            continue;
          auto& mv = mark[static_cast<std::size_t>(v)];
          if (mv == stamp)
            continue;
          mv = stamp;
          block.push_back(v);
          if (static_cast<int>(block.size()) >= max_block)
            break;
          q.emplace_back(v, d + 1);
        }
      }

      blocks.push_back(std::move(block));
    }

    this->schwarz_preconditioner = std::make_unique<SchwarzAdditive_CRS>(crs_ptrs, blocks, schwarz_pivot_min, schwarz_diag_shift);
    std::cout << "Schwarz preconditioner created (core_k=" << core_k << ", overlap_k=" << overlap_k << ", max_core_size=" << max_core << ", max_block_size=" << max_block << ", pivot_min=" << schwarz_pivot_min << ", diag_shift=" << schwarz_diag_shift << ", blocks=" << schwarz_preconditioner->blocks_count << ", avg_block=" << schwarz_preconditioner->avg_block_size << ")." << std::endl;
  } else if (preconditioner_type == "ILUT") {
    schwarz_preconditioner.reset();
    ilu_preconditioner.reset();
    this->ilut_preconditioner = std::make_unique<ILUT_CRS>(crs_ptrs, ilut_drop_tol, ilut_max_entries_per_row, ilut_pivot_min);
    std::cout << "ILUT preconditioner created (drop_tol=" << ilut_drop_tol << ", max_entries_per_row=" << ilut_max_entries_per_row << ", pivot_min=" << ilut_pivot_min << ")." << std::endl;
  } else if (preconditioner_type == "MILU") {
    schwarz_preconditioner.reset();
    ilut_preconditioner.reset();
    this->ilu_preconditioner = std::make_unique<ILU0_CRS>(crs_ptrs, /*modified=*/true, milu_omega);
    std::cout << "MILU(0) preconditioner created (omega=" << milu_omega << ")." << std::endl;
  } else {
    schwarz_preconditioner.reset();
    ilut_preconditioner.reset();
    this->ilu_preconditioner = std::make_unique<ILU0_CRS>(crs_ptrs);
    std::cout << "ILU(0) preconditioner created." << std::endl;
  }
  const auto elapsed = watch_build();
  last_ilu_build_time = elapsed[0];
  std::cout << Magenta << "[BEM] " << Cyan << "build " << precond_name << " preconditioner" << Green << " elapsed=" << elapsed << colorReset << std::endl;
  ilu_geometry_hash_cached = computeILUGeometrySignature();
  preconditioner_type_built = preconditioner_type;
  ilu_neighborhood_type_built = ilu_neighborhood_type;
  ilu_k_ring_num_built = k_ring_num;
  ilu_matrix_size_built = matrix_size;
  ilut_drop_tol_built = ilut_drop_tol;
  ilut_max_entries_per_row_built = ilut_max_entries_per_row;
  ilut_pivot_min_built = ilut_pivot_min;
  schwarz_core_k_built = schwarz_core_k;
  schwarz_overlap_k_built = schwarz_overlap_k;
  schwarz_max_core_size_built = schwarz_max_core_size;
  schwarz_max_block_size_built = schwarz_max_block_size;
  schwarz_pivot_min_built = schwarz_pivot_min;
  schwarz_diag_shift_built = schwarz_diag_shift;
}

/* -------------------------------------------------------------------------- */

void createCopyMap() {
  copy_map.clear();
  copy_map.reserve(points.size() * 2 * 5);
  copy_map_t.clear();
  copy_map_t.reserve(points.size() * 2 * 5);
  for (const auto& p : points) {
    p->phiOnFace_copy = p->phiOnFace;
    p->phinOnFace_copy = p->phinOnFace;
    p->phitOnFace_copy = p->phitOnFace;
    p->phintOnFace_copy = p->phintOnFace;
    for (const auto& [f, i] : p->f2Index) {
      copy_map.push_back({&p->phiOnFace_copy.at(f), &p->phiOnFace.at(f)});
      copy_map.push_back({&p->phinOnFace_copy.at(f), &p->phinOnFace.at(f)});
      copy_map_t.push_back({&p->phitOnFace_copy.at(f), &p->phitOnFace.at(f)});
      copy_map_t.push_back({&p->phintOnFace_copy.at(f), &p->phintOnFace.at(f)});
    }
  }
}

void copyPhiPhin() {
  // Parallel: each element is independent
  _Pragma("omp parallel for") for (std::size_t i = 0; i < copy_map.size(); ++i) {
    auto& [p_copy, value] = copy_map[i];
    *p_copy = *value;
    *value = 0.;
  }
};

void copyPhiPhin_t() {
  // Parallel: each element is independent
  _Pragma("omp parallel for") for (std::size_t i = 0; i < copy_map_t.size(); ++i) {
    auto& [p_copy, value] = copy_map_t[i];
    *p_copy = *value;
    *value = 0.;
  }
};

void restorePhiPhin() {
  for (auto& [p_copy, value] : copy_map)
    *value = *p_copy;
};

void restorePhiPhin_t() {
  for (auto& [p_copy, value] : copy_map_t)
    *value = *p_copy;
};

/* -------------------------------------------------------------------------- */

std::vector<std::tuple<bool, int, double&, double&>> cache_DorN_phi_phin;
std::vector<double> cache_phi_val_D_by_index;
std::vector<double> cache_phin_val_D_by_index;

void cacheBoundaryValues(std::vector<networkPoint*>& points, int total_unknowns) {
  cache_DorN_phi_phin.clear();
  cache_DorN_phi_phin.reserve(total_unknowns);
  cache_phi_val_D_by_index.assign(total_unknowns, 0.0);
  cache_phin_val_D_by_index.assign(total_unknowns, 0.0);
  for (const auto& p : points) {
    // Ensure pointer/reference stability: avoid unordered_map rehash during operator[] inserts.
    p->phiOnFace_FMM.reserve(p->f2Index.size());
    p->phinOnFace_FMM.reserve(p->f2Index.size());
    for (const auto& [f, i] : p->f2Index) {
      cache_DorN_phi_phin.emplace_back(isDirichletID_BEM(p, f), i, p->phiOnFace_FMM[f], p->phinOnFace_FMM[f]);
      if (i >= 0 && i < total_unknowns) {
        cache_phi_val_D_by_index[i] = p->phiOnFace_FMM[f];
        cache_phin_val_D_by_index[i] = p->phinOnFace_FMM[f];
      }
    }
  }
  // Cache edge midpoint DOFs for true quadratic elements (per-face maps)
  // Note: For phi_t solve, phiOnFace/phinOnFace are swapped with phitOnFace/phintOnFace
  // before calling this function, so we always read from phiOnFace_FMM/phinOnFace_FMM.
  if (use_true_quadratic_element) {
    for (auto* water : WATERS) {
      for (auto* l : water->getBoundaryLines()) {
        l->phiOnFace_FMM.reserve(l->f2Index.size());
        l->phinOnFace_FMM.reserve(l->f2Index.size());
        for (const auto& [f, i] : l->f2Index) {
          if (i >= 0 && i < total_unknowns) {
            cache_DorN_phi_phin.emplace_back(isDirichletID_BEM(l, f), i, l->phiOnFace_FMM[f], l->phinOnFace_FMM[f]);
            cache_phi_val_D_by_index[i] = l->phiOnFace_FMM[f];
            cache_phin_val_D_by_index[i] = l->phinOnFace_FMM[f];
          }
        }
      }
    }
  }
}

void setUnknowns(double value) {
  _Pragma("omp parallel for") for (auto& [isDirichlet, i, phi, phin] : cache_DorN_phi_phin) {
    if (isDirichlet) {
      phin = value;
      cache_phin_val_D_by_index[i] = value;
    } else {
      phi = value;
      cache_phi_val_D_by_index[i] = value;
    }
  }
}

void setUnknowns(const std::vector<double>& values) {
  _Pragma("omp parallel for") for (auto& [isDirichlet, i, phi, phin] : cache_DorN_phi_phin) {
    if (isDirichlet) {
      phin = values[i];
      cache_phin_val_D_by_index[i] = values[i];
    } else {
      phi = values[i];
      cache_phi_val_D_by_index[i] = values[i];
    }
  }
}

void setKnowns(double value) {
  _Pragma("omp parallel for") for (auto& [isDirichlet, i, phi, phin] : cache_DorN_phi_phin) {
    if (isDirichlet) {
      phi = value;
      cache_phi_val_D_by_index[i] = value;
    } else {
      phin = value;
      cache_phin_val_D_by_index[i] = value;
    }
  }
}

void initializeBucket() {
  std::cout << "Setting up buckets..." << std::endl;
  auto obj = WATERS[0];
  auto Center = obj->getCenter();
  auto dL = 1.05 * obj->getScale() / 2.;

  CoordinateBounds bounds(Center[0] - dL, Center[0] + dL, Center[1] - dL, Center[1] + dL, Center[2] - dL, Center[2] + dL);

  // Safe re-init: `Buckets::initialize()` does NOT delete existing child nodes.
  // If we ever re-run initializeBucket() (e.g. when sources escape the root bounds),
  // make sure the old tree is fully destroyed first.
  B_poles.traverseTree([](auto& child) {
    delete child;
    child = nullptr;
  });
  B_poles.children.clear();
  B_poles.level_buckets.clear();
  B_poles.deepest_level_buckets.clear();

  B_poles.initialize(bounds, dL);
  B_poles.setLevel(0, fmm_max_level);
  B_poles.setGrowCondition([](auto bucket) { return bucket->data1D.empty() ? false : (bucket->data1D.size() > static_cast<std::size_t>(fmm_bucket_max_points) && bucket->level < bucket->max_level); });
  std::cout << "Buckets set up complete." << std::endl;
}

void copyToFMM(bool is_time_derivative) {
  // Coordinate scaling math: x' = x/L implies ∂φ/∂n' = L × ∂φ/∂n
  // For Neumann BC input: physical value g → scaled value L×g
  const double phin_scale = (use_coordinate_scaling_ && coordinate_scale_factor_ > 1e-10) ? coordinate_scale_factor_ : 1.0;

  if (is_time_derivative) {
    _Pragma("omp parallel for") for (int i = 0; i < points.size(); ++i) {
      auto& p = points[i];
      if (p->phiOnFace_FMM.empty()) {
        p->phiOnFace_FMM = p->phitOnFace;
        p->phinOnFace_FMM = p->phintOnFace;
        // Apply scaling to Neumann BCs
        if (phin_scale != 1.0) {
          for (auto& [f, val] : p->phinOnFace_FMM) {
            if (isNeumannID_BEM(p, f))
              val *= phin_scale;
          }
        }
      } else {
        for (const auto& [f, val] : p->phitOnFace) {
          p->phiOnFace_FMM[f] = val;
        }
        for (const auto& [f, val] : p->phintOnFace) {
          double scaled_val = val;
          if (phin_scale != 1.0 && isNeumannID_BEM(p, f))
            scaled_val *= phin_scale;
          p->phinOnFace_FMM[f] = scaled_val;
        }
      }
    }
  } else {
    _Pragma("omp parallel for") for (int i = 0; i < points.size(); ++i) {
      auto& p = points[i];
      if (p->phiOnFace_FMM.empty()) {
        p->phiOnFace_FMM = p->phiOnFace;
        p->phinOnFace_FMM = p->phinOnFace;
        // Apply scaling to Neumann BCs
        if (phin_scale != 1.0) {
          for (auto& [f, val] : p->phinOnFace_FMM) {
            if (isNeumannID_BEM(p, f))
              val *= phin_scale;
          }
        }
      } else {
        for (const auto& [f, val] : p->phiOnFace) {
          p->phiOnFace_FMM[f] = val;
        }
        for (const auto& [f, val] : p->phinOnFace) {
          double scaled_val = val;
          if (phin_scale != 1.0 && isNeumannID_BEM(p, f))
            scaled_val *= phin_scale;
          p->phinOnFace_FMM[f] = scaled_val;
        }
      }
    }
  }

  // Copy midpoint per-face maps to FMM copies (with Neumann phin scaling)
  // Note: always use phiOnFace/phinOnFace regardless of is_time_derivative,
  // because the phi_t solve uses a swap mechanism (phi_mid <-> phi_t_mid,
  // phiOnFace <-> phitOnFace) so the "current working values" are always
  // in phiOnFace/phinOnFace.
  if (use_true_quadratic_element) {
    for (auto* water : WATERS) {
      for (auto* l : water->getBoundaryLines()) {
        l->phiOnFace_FMM = l->phiOnFace;
        l->phinOnFace_FMM = l->phinOnFace;
        if (phin_scale != 1.0) {
          for (auto& [f, val] : l->phinOnFace_FMM) {
            if (isNeumannID_BEM(l, f))
              val *= phin_scale;
          }
        }
      }
    }
  }
}

bool createSourcesOnSurfaces() {
  std::cout << "Creating sources on surfaces..." << std::endl;
  auto obj = WATERS[0];

  auto faces = obj->getBoundaryFaces();
  std::atomic<int> reuse_count{0};
  const char* env_line_factor2 = std::getenv("BEM_LINE_FACTOR2_DEBUG");
  const bool debug_line_factor2 = (env_line_factor2 && std::string(env_line_factor2) != "0");
  std::atomic<std::size_t> dbg_phi_checked{0}, dbg_phin_checked{0};
  std::atomic<std::size_t> dbg_phi_ratio2{0}, dbg_phin_ratio2{0};
  std::atomic<int> dbg_examples{0};
  // Parallel source creation: each face's sources are independent (no race condition)
  // Use dynamic scheduling because pseudo-quadratic elements have more sources
  _Pragma("omp parallel for schedule(dynamic, 100)") for (std::size_t face_idx = 0; face_idx < faces.size(); ++face_idx) {
    auto* F = faces[face_idx];
    const auto [p0, p1, p2] = F->Points;

    // 各 face 単位でソース数をチェックし、一致すれば再利用、しなければ新規作成
    // 新設計: 1面1要素（内部に複数P2Mソースを内包）
    const std::size_t expected_sources = 1;
    const bool reuse_sources = (F->sources.size() == expected_sources);
    if (reuse_sources) {
      ++reuse_count;
    } else {
      F->sources.clear();
    }

    //! ソースとして重要なのは，G,Gnを作るための情報のみ. quadpointsには，points_facesが入っている
    //  現在のところ，sourceはupdateで位置などを更新しrebinするのではなく，再度生成して新しくしrebinする方法をとっている
    // create sources
    if (F->isTrueQuadraticElement) {
      /* -------------------------------------------------------------------------- */
      // True quadratic element: 6 DOFs (3 vertices + 3 edge midpoints) with TriShape<6>
      // Geometry: 2次補間（6節点）, Interpolation: ハイブリッド再分配
      /* -------------------------------------------------------------------------- */
      auto dof_indices = getQuadDOFIndices(F);
      // dof_indices を networkFace にも保存
      for (int k = 0; k < 6; ++k)
        F->true_quad_dof_indices[k] = dof_indices[k];
      if (const char* env = std::getenv("BEM_TRUE_QUAD_DOF_DEBUG"); env && std::string(env) != "0") {
        std::array<int, 6> idx = dof_indices;
        std::size_t dup_count = 0;
        for (int a = 0; a < 6; ++a)
          for (int b = a + 1; b < 6; ++b)
            if (idx[a] >= 0 && idx[b] >= 0 && idx[a] == idx[b])
              ++dup_count;
        if (dup_count > 0) {
          static std::atomic<int> warn_dup{0};
          const int n = warn_dup.fetch_add(1, std::memory_order_relaxed);
          if (n < 20) {
            std::cout << Red << "[true-quad:dof-dup] face=" << F
                      << " dup_pairs=" << dup_count
                      << " idx={" << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << "," << idx[4] << "," << idx[5] << "}"
                      << colorReset << std::endl;
          }
        }
      }

      auto X012 = std::array<Tddd, 3>{p0->X, p1->X, p2->X};
      const auto cross_vec = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      const auto normal_vec = Normalize(cross_vec);
      const double J_det = Norm(cross_vec);

      // 6節点座標: 3頂点 + 3辺中点（2次幾何補間用）
      auto [l0_qp, l1_qp, l2_qp] = F->Lines;
      if (const char* env = std::getenv("BEM_TRUE_QUAD_ORDER_DEBUG"); env && std::string(env) != "0") {
        auto edge_matches = [](const networkLine* l, const networkPoint* a, const networkPoint* b) {
          auto [u, v] = l->getPoints();
          return (u == a && v == b) || (u == b && v == a);
        };
        if (!edge_matches(l0_qp, p0, p1) || !edge_matches(l1_qp, p1, p2) || !edge_matches(l2_qp, p2, p0)) {
          static std::atomic<int> warn_count{0};
          const int wid = warn_count.fetch_add(1, std::memory_order_relaxed);
          if (wid < 20) {
            std::cout << Red << "[true-quad:order] mismatch face=" << F
                      << " l0(p0,p1)=" << edge_matches(l0_qp, p0, p1)
                      << " l1(p1,p2)=" << edge_matches(l1_qp, p1, p2)
                      << " l2(p2,p0)=" << edge_matches(l2_qp, p2, p0)
                      << colorReset << std::endl;
          }
        }
      }
      const T6Tddd X6 = {p0->X, p1->X, p2->X, l0_qp->X_mid, l1_qp->X_mid, l2_qp->X_mid};

      // Pointers to DOF values: [phi, phin] for each of the 6 nodes
      // Vertices use phiOnFace_FMM, midpoints use per-face phiOnFace_FMM/phinOnFace_FMM
      auto [l0, l1, l2] = F->Lines;
      auto key0 = std::get<1>(pf2ID(p0, F));
      auto key1 = std::get<1>(pf2ID(p1, F));
      auto key2 = std::get<1>(pf2ID(p2, F));
      auto lf_key0 = std::get<1>(lf2ID(l0, F));
      auto lf_key1 = std::get<1>(lf2ID(l1, F));
      auto lf_key2 = std::get<1>(lf2ID(l2, F));
      std::array<std::array<double*, 2>, 6> pair_pointer_phiphin_6 = {{{&p0->phiOnFace_FMM.at(key0), &p0->phinOnFace_FMM.at(key0)},
                                                                       {&p1->phiOnFace_FMM.at(key1), &p1->phinOnFace_FMM.at(key1)},
                                                                       {&p2->phiOnFace_FMM.at(key2), &p2->phinOnFace_FMM.at(key2)},
                                                                       {&l0->phiOnFace_FMM.at(lf_key0), &l0->phinOnFace_FMM.at(lf_key0)},
                                                                       {&l1->phiOnFace_FMM.at(lf_key1), &l1->phinOnFace_FMM.at(lf_key1)},
                                                                       {&l2->phiOnFace_FMM.at(lf_key2), &l2->phinOnFace_FMM.at(lf_key2)}}};
      if (debug_line_factor2) {
        auto check_ratio = [&](const char* label,
                               int edge_local,
                               double mid_val,
                               double avg_val,
                               std::atomic<std::size_t>& checked,
                               std::atomic<std::size_t>& near2) {
          constexpr double eps = 1e-14;
          if (std::abs(avg_val) < eps)
            return;
          checked.fetch_add(1, std::memory_order_relaxed);
          const double ratio = mid_val / avg_val;
          if (std::abs(std::abs(ratio) - 2.0) < 0.2) {
            near2.fetch_add(1, std::memory_order_relaxed);
            const int n = dbg_examples.fetch_add(1, std::memory_order_relaxed);
            if (n < 20) {
              std::cout << Red << "[line:x2] " << Cyan
                        << "face=" << F
                        << " edge=" << edge_local
                        << " " << label << "_mid=" << mid_val
                        << " avg_end=" << avg_val
                        << " ratio=" << ratio
                        << colorReset << std::endl;
            }
          }
        };
        // edge 0: p0-p1 <-> dof3, edge 1: p1-p2 <-> dof4, edge 2: p2-p0 <-> dof5
        for (int e = 0; e < 3; ++e) {
          const int a = e;
          const int b = (e + 1) % 3;
          const int m = 3 + e;
          const double phi_avg = 0.5 * (*pair_pointer_phiphin_6[a][0] + *pair_pointer_phiphin_6[b][0]);
          const double phin_avg = 0.5 * (*pair_pointer_phiphin_6[a][1] + *pair_pointer_phiphin_6[b][1]);
          const double phi_mid = *pair_pointer_phiphin_6[m][0];
          const double phin_mid = *pair_pointer_phiphin_6[m][1];
          check_ratio("phi", e, phi_mid, phi_avg, dbg_phi_checked, dbg_phi_ratio2);
          check_ratio("phin", e, phin_mid, phin_avg, dbg_phin_checked, dbg_phin_ratio2);
        }
      }

      // P2M rule
      const auto& p2m_rule = getDunavantP2MRule(g_p2m_quadrature_points);
      const int n_p2m = static_cast<int>(p2m_rule.size());

      // Adjacent callback for true quadratic elements (Duffy + rigid-mode).
      // Called when isAdjacentTo(target)=true, i.e., target is one of {p0, p1, p2, l0, l1, l2}.
      // Handles both vertex Duffy and midpoint Duffy(2分割) by on-the-fly quadrature.
      auto fill_direct_true_quad = [F, dof_indices, X6,
                                    points = F->Points](const target4FMM* origin, DirectAccumulator& acc) -> void {
        auto [p0, p1, p2] = points;
        constexpr std::array<bool, 3> all_true{true, true, true};

        std::array<std::array<double, 2>, 6> WGN_WGnN;
        WGN_WGnN.fill({0., 0.});

        auto redistribute_corner_midpoint = [&]() {
          for (int k = 3; k < 6; ++k) {
            if (dof_indices[k] < 0) {
              int v1 = k - 3, v2 = (k - 3 + 1) % 3;
              WGN_WGnN[v1][0] += 0.5 * WGN_WGnN[k][0];
              WGN_WGnN[v1][1] += 0.5 * WGN_WGnN[k][1];
              WGN_WGnN[v2][0] += 0.5 * WGN_WGnN[k][0];
              WGN_WGnN[v2][1] += 0.5 * WGN_WGnN[k][1];
              WGN_WGnN[k] = {0., 0.};
            }
          }
        };

        auto flush_acc = [&]() {
          redistribute_corner_midpoint();
          for (int k = 0; k < 6; ++k)
            if (dof_indices[k] >= 0)
              acc.add(dof_indices[k], WGN_WGnN[k][0], WGN_WGnN[k][1]);
        };

        // --- Vertex Duffy: origin is one of {p0, p1, p2} ---
        if (origin == p0 || origin == p1 || origin == p2) {
          std::array<int, 3> vperm = {0, 1, 2};
          if (p1 == origin)
            vperm = {1, 2, 0};
          else if (p2 == origin)
            vperm = {2, 0, 1};

          const int vertex_local = vperm[0];

          for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
            const auto bary_loc = ModTriShape<3>(t0, t1);
            std::array<double, 3> bary{};
            bary[vperm[0]] = bary_loc[0];
            bary[vperm[1]] = bary_loc[1];
            bary[vperm[2]] = bary_loc[2];

            auto N6_geo = TriShape<6>(bary[0], bary[1], all_true);
            auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true);
            auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true);
            auto X_q = Dot(N6_geo, X6);
            auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
            auto N6 = F->trueQuadN6(bary[0], bary[1]);

            Tddd R;
            double nr = Norm(R = (X_q - origin->Xtarget));
            double tmp = ww * (1. - t0) / nr;
            double WG = Norm(cross_q) * tmp;
            double WGn = -Dot(R / nr, cross_q) * tmp / nr;
            for (int k = 0; k < 6; ++k) {
              WGN_WGnN[k][0] += WG * N6[k];
              WGN_WGnN[k][1] += (k == vertex_local) ? 0. : WGn * N6[k];
            }
          }
          flush_acc();
          return;
        }

        // --- Midpoint Duffy: origin is one of {l0, l1, l2} (edge midpoint) ---
        const auto& [fl0, fl1, fl2] = F->Lines;
        int on_edge = (origin == static_cast<const target4FMM*>(fl0))   ? 0
                      : (origin == static_cast<const target4FMM*>(fl1)) ? 1
                      : (origin == static_cast<const target4FMM*>(fl2)) ? 2
                                                                        : -1;
        if (on_edge < 0)
          return;

        int va = on_edge, vb = (on_edge + 1) % 3, vc = (on_edge + 2) % 3;
        const int mid_local = 3 + on_edge;

        for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
          double omt0 = 1.0 - t0;
          double s1 = t1 * omt0;
          double s2 = omt0 * (1.0 - t1);

          // Sub-triangle 1: (M, vb, vc)
          {
            std::array<double, 3> bary{};
            bary[va] = 0.5 * t0;
            bary[vb] = 0.5 * t0 + s1;
            bary[vc] = s2;
            auto N6_geo = TriShape<6>(bary[0], bary[1], all_true);
            auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true);
            auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true);
            auto X_q = Dot(N6_geo, X6);
            auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
            auto N6 = F->trueQuadN6(bary[0], bary[1]);

            Tddd R;
            double nr = Norm(R = (X_q - origin->Xtarget));
            double tmp = ww * omt0 * 0.5 / nr;
            double WG = Norm(cross_q) * tmp;
            double WGn = -Dot(R / nr, cross_q) * tmp / nr;
            for (int k = 0; k < 6; ++k) {
              WGN_WGnN[k][0] += WG * N6[k];
              WGN_WGnN[k][1] += (k == mid_local) ? 0. : WGn * N6[k];
            }
          }

          // Sub-triangle 2: (M, vc, va)
          {
            std::array<double, 3> bary{};
            bary[va] = 0.5 * t0 + s2;
            bary[vb] = 0.5 * t0;
            bary[vc] = s1;
            auto N6_geo = TriShape<6>(bary[0], bary[1], all_true);
            auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true);
            auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true);
            auto X_q = Dot(N6_geo, X6);
            auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
            auto N6 = F->trueQuadN6(bary[0], bary[1]);

            Tddd R;
            double nr = Norm(R = (X_q - origin->Xtarget));
            double tmp = ww * omt0 * 0.5 / nr;
            double WG = Norm(cross_q) * tmp;
            double WGn = -Dot(R / nr, cross_q) * tmp / nr;
            for (int k = 0; k < 6; ++k) {
              WGN_WGnN[k][0] += WG * N6[k];
              WGN_WGnN[k][1] += (k == mid_local) ? 0. : WGn * N6[k];
            }
          }
        }
        flush_acc();
      };

      // Non-adjacent callback for true_quadratic: pseudo_quadと同様にDunavantを動的選択。
      // near -> g_p2m_quadrature_points, far -> 6点
      auto fill_nonadj_true_quad = [F, near_region = obj->getScale() / 25., X6](const target4FMM* origin, DirectAccumulator& acc) -> void {
        const auto& rule = (Norm(F->centroid - origin->Xtarget) < near_region)
                               ? getDunavantP2MRule(g_p2m_quadrature_points)
                               : getDunavantP2MRule(6);
        constexpr std::array<bool, 3> all_true{true, true, true};

        std::array<std::array<double, 2>, 6> WGN_WGnN;
        WGN_WGnN.fill({0., 0.});
        for (const auto& [b0, b1, b2, ww] : rule) {
          auto N6_geo = TriShape<6>(b0, b1, all_true);
          auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
          auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
          auto X_q = Dot(N6_geo, X6);
          auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
          auto N6 = F->trueQuadN6(b0, b1);

          auto R = X_q - origin->Xtarget;
          double nr = Norm(R);
          double nr_inv = 1.0 / nr;
          double J_qp = Norm(cross_q);
          double WG = J_qp * ww * nr_inv;
          double WGn = -Dot(R * nr_inv, cross_q) * ww * nr_inv * nr_inv;
          for (int k = 0; k < 6; ++k) {
            WGN_WGnN[k][0] += WG * N6[k];
            WGN_WGnN[k][1] += WGn * N6[k];
          }
        }
        for (int k = 0; k < 6; ++k)
          if (F->true_quad_dof_indices[k] >= 0)
            acc.add(F->true_quad_dof_indices[k], WGN_WGnN[k][0], WGN_WGnN[k][1]);
      };

      // Setup source element
      auto setup_true_quad_element = [&](source4FMM<target4FMM>& pole) {
        pole.X = F->centroid;
        pole.normal = normal_vec;
        pole.fill_direct_entries = fill_direct_true_quad;
        pole.fill_direct_entries_nonadj = fill_nonadj_true_quad;
        pole.face_vertices = {p0, p1, p2};
        auto [l0_s, l1_s, l2_s] = F->Lines;
        pole.face_midpoints = {static_cast<const target4FMM*>(l0_s),
                               static_cast<const target4FMM*>(l1_s),
                               static_cast<const target4FMM*>(l2_s)};
        pole.get_X = nullptr;
        pole.get_normal = nullptr;
        pole.get_weighted_source_densities = nullptr;
        // P2M sources using Dunavant quadrature
        pole.p2m_sources.resize(n_p2m);
        constexpr std::array<bool, 3> all_true{true, true, true};
        for (int q = 0; q < n_p2m; ++q) {
          auto& isrc = pole.p2m_sources[q];
          const auto& [b0, b1, b2, ww] = p2m_rule[q];
          const auto N6_geo = TriShape<6>(b0, b1, all_true);
          const auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
          const auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
          const Tddd Xq = Dot(N6_geo, X6);
          const Tddd cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
          const double J_q = Norm(cross_q);
          const Tddd normal_q = (J_q > 1e-30) ? (cross_q / J_q) : normal_vec;

          // Use the same hybrid midpoint redistribution rule as direct integration.
          const auto N6_q = F->trueQuadN6(b0, b1);

          isrc.X = Xq;
          isrc.normal = normal_q;
          isrc.get_weighted_source_densities = [W = J_q * ww, N6_q, pair_pointer_phiphin_6]() -> std::array<double, 2> {
            double phiWN = 0., phinWN = 0.;
            for (int k = 0; k < 6; ++k) {
              const double phi_k = *pair_pointer_phiphin_6[k][0];
              const double phin_k = *pair_pointer_phiphin_6[k][1];
              phiWN += W * N6_q[k] * phi_k;
              phinWN += W * N6_q[k] * phin_k;
            }
            return {phiWN, phinWN};
          };
        }
        pole.max_source_offset = 0.0;
        for (const auto& isrc : pole.p2m_sources)
          pole.max_source_offset = std::max(pole.max_source_offset, Norm(isrc.X - pole.X));
      };

      if (reuse_sources) {
        setup_true_quad_element(*F->sources[0]);
      } else {
        auto pole = std::make_shared<source4FMM<target4FMM>>();
        setup_true_quad_element(*pole);
        F->sources.emplace_back(std::move(pole));
      }

    } else if (F->isPseudoQuadraticElement) {
      /* -------------------------------------------------------------------------- */

      /*
      DodecaPointsのvertexは，quadpoint, quadpoint_l0, quadpoint_l1, quadpoint_l2
      それに対応する形状関数は，Nc_N0_N1_N2
      として保存されている．
      */

      auto insert = [](auto& insert_to, const auto& points_faces) {
        int i = 0;
        for (auto& [p, f] : points_faces) {
          auto key = std::get<1>(pf2ID(p, f));
          insert_to[i++] = {&p->phiOnFace_FMM.at(key), &p->phinOnFace_FMM.at(key)}; //! updateされる値へのポインタ かつ キー．
        }
      };
      std::array<std::array<std::array<double*, 2>, 6>, 4> pair_pointer_phiphin;
      insert(pair_pointer_phiphin[0], F->dodecaPoints[0]->quadpoint.points_faces);
      insert(pair_pointer_phiphin[1], F->dodecaPoints[0]->quadpoint_l0.points_faces);
      insert(pair_pointer_phiphin[2], F->dodecaPoints[0]->quadpoint_l1.points_faces);
      insert(pair_pointer_phiphin[3], F->dodecaPoints[0]->quadpoint_l2.points_faces);

      // b% [1] P2Mソース: Dunavant求積則による多点配置（pseudo-quadratic要素）
      const auto& p2m_rule = getDunavantP2MRule(g_p2m_quadrature_points);
      const int n_p2m = static_cast<int>(p2m_rule.size());
      auto dodecapoint0 = F->dodecaPoints[0]; // P2M用はdodecaPoints[0]を使用

      //^  ----------------------------- set up direct integration -------------------------------------------- */
      //^ 直接積分の callback（要素レベル、1面1回）
      // Vertex-on-element callback for pseudo quadratic elements (Duffy + rigid-mode).
      // Called only when isAdjacentTo(target)=true, i.e., target is one of {p0, p1, p2}.
      // Non-adjacent cases are handled by fill_nonadj_pseudo_quad.
      auto fill_direct_pseudo_quad = [F, points = F->Points](const target4FMM* origin, DirectAccumulator& acc) -> void {
        //! 観測点(target)を特定の番号に位置するようにする．
        auto dodecapoint = F->dodecaPoints[points[1] == origin ? 1 : (points[2] == origin ? 2 : 0)];

        std::array<std::array<double, 2>, 24> ij_WGN_WGnN;
        ij_WGN_WGnN.fill({0., 0.});
        std::array<int32_t, 24> ij_indices;
        std::array<std::array<double, 6>, 4> Nc_N0_N1_N2;
        double nr, tmp, WG, WGn_, WGn;
        Tddd R, cross;

        for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
          auto [xi0, xi1, xi2] = ModTriShape<3>(t0, t1);
          cross = dodecapoint->cross(xi0, xi1);
          Nc_N0_N1_N2 = dodecapoint->N6_new(xi0, xi1);
          nr = Norm(R = (dodecapoint->X(xi0, xi1) - origin->Xtarget));
          tmp = ww * (1. - t0) / nr;
          WG = Norm(cross) * tmp;
          WGn_ = -Dot(R / nr, cross) * tmp / nr;
          for (int i = 0; i < 4; ++i) {
            auto quadpoints = dodecapoint->getPseudoQuadPatch(i - 1);
            for (int j = 0; j < 6; ++j) {
              WGn = WGn_;
              auto [p, _] = quadpoints->points_faces[j];
              ij_WGN_WGnN[6 * i + j] += std::array<double, 2>{WG * Nc_N0_N1_N2[i][j], p == origin ? 0. : WGn * Nc_N0_N1_N2[i][j]};
              auto [p_id, key] = pf2ID(quadpoints->points_faces[j]);
              ij_indices[6 * i + j] = static_cast<int32_t>(p_id->f2Index.at(key));
            }
          }
        }

        for (int k = 0; k < 24; ++k)
          acc.add(ij_indices[k], ij_WGN_WGnN[k][0], ij_WGN_WGnN[k][1]);
      };

      //^ 非隣接面用: dodecaPoints[0]固定・リジッドモード不要・Dunavant求積
      auto fill_nonadj_pseudo_quad = [F, near_region = obj->getScale() / 25.,
                                      dodecapoint = F->dodecaPoints[0]](const target4FMM* origin, DirectAccumulator& acc) -> void {
        std::array<std::array<double, 2>, 24> ij_WGN_WGnN;
        ij_WGN_WGnN.fill({0., 0.});
        std::array<int32_t, 24> ij_indices;

        // near→Dunavant(p2m設定), far→Dunavant 6
        const auto& rule = (Norm(F->centroid - origin->Xtarget) < near_region)
                               ? getDunavantP2MRule(g_p2m_quadrature_points)
                               : getDunavantP2MRule(6);
        for (const auto& [b0, b1, b2, ww] : rule) {
          auto cross = dodecapoint->cross(b0, b1);
          auto Nc_N0_N1_N2 = dodecapoint->N6_new(b0, b1);
          std::array<double, 3> R;
          double nr = Norm(R = (dodecapoint->X(b0, b1) - origin->Xtarget));
          double nr_inv = 1.0 / nr;
          double WG = Norm(cross) * ww * nr_inv;
          double WGn_ = -Dot(R * nr_inv, cross) * ww * nr_inv * nr_inv;
          for (int i = 0; i < 4; ++i) {
            auto quadpoints = dodecapoint->getPseudoQuadPatch(i - 1);
            for (int j = 0; j < 6; ++j) {
              ij_WGN_WGnN[6 * i + j] += std::array<double, 2>{WG * Nc_N0_N1_N2[i][j], WGn_ * Nc_N0_N1_N2[i][j]};
              auto [p_id, key] = pf2ID(quadpoints->points_faces[j]);
              ij_indices[6 * i + j] = static_cast<int32_t>(p_id->f2Index.at(key));
            }
          }
        }

        for (int k = 0; k < 24; ++k)
          acc.add(ij_indices[k], ij_WGN_WGnN[k][0], ij_WGN_WGnN[k][1]);
      };
      //^ -------------------------------------------------------------------------- */

      // 要素の設定（共通部分）
      auto setup_pseudo_quad_element = [&](source4FMM<target4FMM>& pole) {
        pole.X = F->centroid; // 重心（ツリー配置用）
        pole.normal = F->normal;
        pole.fill_direct_entries = fill_direct_pseudo_quad;
        pole.fill_direct_entries_nonadj = fill_nonadj_pseudo_quad;
        pole.face_vertices = {p0, p1, p2};
        // レガシーコールバックをクリア（p2m_sourcesを使用するため）
        pole.get_X = nullptr;
        pole.get_normal = nullptr;
        pole.get_weighted_source_densities = nullptr;
        // P2Mソースの設定（Dunavant求積点 × pseudo-quadratic幾何）
        pole.p2m_sources.resize(n_p2m);
        for (int q = 0; q < n_p2m; ++q) {
          auto& isrc = pole.p2m_sources[q];
          const auto& [b0, b1, b2, ww] = p2m_rule[q];
          // pseudo-quadratic幾何でのX, cross, N6を評価
          auto cross_q = dodecapoint0->cross(b0, b1);
          isrc.X = dodecapoint0->X(b0, b1);
          isrc.normal = Normalize(cross_q);
          double J_det_q = Norm(cross_q);
          auto Nc_N0_N1_N2_q = dodecapoint0->N6_new(b0, b1);
          isrc.get_weighted_source_densities = [W = J_det_q * ww, Nc_N0_N1_N2_q, pair_pointer_phiphin]() -> std::array<double, 2> {
            double phiWN = 0., phinWN = 0.;
            for (int i = 0; i < 4; ++i)
              for (int j = 0; j < 6; ++j) {
                phiWN += W * Nc_N0_N1_N2_q[i][j] * (*pair_pointer_phiphin[i][j][0]);
                phinWN += W * Nc_N0_N1_N2_q[i][j] * (*pair_pointer_phiphin[i][j][1]);
              }
            return std::array<double, 2>{phiWN, phinWN};
          };
        }
        // MAC拡張用: 重心から最も遠い内部ソースまでの距離を計算
        pole.max_source_offset = 0.0;
        for (const auto& isrc : pole.p2m_sources)
          pole.max_source_offset = std::max(pole.max_source_offset, Norm(isrc.X - pole.X));
      };

      if (reuse_sources) {
        setup_pseudo_quad_element(*F->sources[0]);
      } else {
        auto pole = std::make_shared<source4FMM<target4FMM>>();
        setup_pseudo_quad_element(*pole);
        F->sources.emplace_back(std::move(pole));
      }

    } else {
      auto key0 = std::get<1>(pf2ID(p0, F));
      auto key1 = std::get<1>(pf2ID(p1, F));
      auto key2 = std::get<1>(pf2ID(p2, F));
      std::array<double*, 2> pair_pointer_phiphin0 = {&p0->phiOnFace_FMM.at(key0), &p0->phinOnFace_FMM.at(key0)};
      std::array<double*, 2> pair_pointer_phiphin1 = {&p1->phiOnFace_FMM.at(key1), &p1->phinOnFace_FMM.at(key1)};
      std::array<double*, 2> pair_pointer_phiphin2 = {&p2->phiOnFace_FMM.at(key2), &p2->phinOnFace_FMM.at(key2)};
      const int32_t idx0 = static_cast<int32_t>(p0->f2Index.at(key0));
      const int32_t idx1 = static_cast<int32_t>(p1->f2Index.at(key1));
      const int32_t idx2 = static_cast<int32_t>(p2->f2Index.at(key2));

      // b% [1] P2Mソース: Dunavant求積則による多点配置
      auto X012 = ToX(F->getPoints(p0));
      const auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]); // 線形要素では定数
      const auto normal_vec = Normalize(cross);
      const double J_det = Norm(cross);

      const auto& p2m_rule = getDunavantP2MRule(g_p2m_quadrature_points);
      const int n_p2m = static_cast<int>(p2m_rule.size());

      //!  ----------------------------- set up direct integration -------------------------------------------- */

      //$ 直接積分の callback（要素レベル、1面1回）
      const double near_region = obj->getScale() / 25.;
      std::array<std::array<double, 3>, 3> X012_ = {F->Points[0]->Xtarget, F->Points[1]->Xtarget, F->Points[2]->Xtarget};
      auto fill_direct_linear = [F, obj, near_region,
                                 J_det = 2. * F->area,
                                 cross = 2. * F->area * F->normal,
                                 X012_,
                                 points = F->Points, idx0, idx1, idx2](const target4FMM* origin, DirectAccumulator& acc) -> void {
        auto [p0, p1, p2] = points;
        std::array<int, 3> idx = {0, 1, 2};
        std::array<std::array<double, 3>, 3> X012 = X012_;
        if (p1 == origin) {
          p2 = points[0];
          p0 = points[1]; // origin
          p1 = points[2];
          idx = {2, 0, 1};
          X012 = {p0->Xtarget, p1->Xtarget, p2->Xtarget};
        } else if (p2 == origin) {
          p1 = points[0];
          p2 = points[1];
          p0 = points[2]; // origin
          idx = {1, 2, 0};
          X012 = {p0->Xtarget, p1->Xtarget, p2->Xtarget};
        }

        double nr_inv;
        std::array<double, 3> WGN = {0., 0., 0.}, WGnN = {0., 0., 0.}, N012_geometry, R, tmp;
        // b% [2] ソース点を直接積分に置き換えたあとの，積分

        for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
          N012_geometry = ModTriShape<3>(t0, t1);
          nr_inv = 1. / Norm(R = (Dot(N012_geometry, X012) - origin->Xtarget));
          tmp = (ww * (1. - t0) * nr_inv) * N012_geometry;
          WGN += J_det * tmp;
          WGnN += -Dot(R * nr_inv, cross) * tmp * nr_inv;
        }
        acc.add(idx0, WGN[idx[0]], 0.0);
        acc.add(idx1, WGN[idx[1]], WGnN[idx[1]]);
        acc.add(idx2, WGN[idx[2]], WGnN[idx[2]]);
      };

      //^ 非隣接面用: 頂点巡回置換・リジッドモード不要の簡略化パス
      auto fill_nonadj_linear = [near_region,
                                 J_det = 2. * F->area,
                                 cross = 2. * F->area * F->normal,
                                 X012_,
                                 Xc = F->centroid,
                                 idx0, idx1, idx2](const target4FMM* origin, DirectAccumulator& acc) -> void {
        double nr_inv;
        std::array<double, 3> WGN = {0., 0., 0.}, WGnN = {0., 0., 0.}, N012_geometry, R, tmp;
        if (Norm(Xc - origin->Xtarget) < near_region) {
          for (const auto& [t0, t1, ww] : __array_GW5xGW5__) {
            N012_geometry = ModTriShape<3>(t0, t1);
            nr_inv = 1. / Norm(R = (Dot(N012_geometry, X012_) - origin->Xtarget));
            tmp = (ww * (1. - t0) * N012_geometry) * nr_inv;
            WGN += J_det * tmp;
            WGnN -= Dot(R * nr_inv, cross) * (tmp * nr_inv);
          }
        } else {
          for (const auto& [b0, b1, b2, ww] : getDunavantP2MRule(g_p2m_quadrature_points)) {
            N012_geometry = {b0, b1, b2};
            nr_inv = 1. / Norm(R = (Dot(N012_geometry, X012_) - origin->Xtarget));
            tmp = ww * N012_geometry * nr_inv;
            WGN += J_det * tmp;
            WGnN -= Dot(R * nr_inv, cross) * (tmp * nr_inv);
          }
        }
        acc.add(idx0, WGN[0], WGnN[0]);
        acc.add(idx1, WGN[1], WGnN[1]);
        acc.add(idx2, WGN[2], WGnN[2]);
      };

      // 要素の設定（共通部分）
      auto setup_linear_element = [&](source4FMM<target4FMM>& pole) {
        pole.X = Dot(std::array<double, 3>{1. / 3., 1. / 3., 1. / 3.}, X012); // 重心（ツリー配置用）
        pole.normal = normal_vec;
        pole.fill_direct_entries = fill_direct_linear;
        pole.fill_direct_entries_nonadj = fill_nonadj_linear;
        pole.face_vertices = {p0, p1, p2};
        pole.dof_indices[0] = idx0;
        pole.dof_indices[1] = idx1;
        pole.dof_indices[2] = idx2;
        pole.near_region = static_cast<float>(near_region);
        // レガシーコールバックをクリア（p2m_sourcesを使用するため）
        pole.get_X = nullptr;
        pole.get_normal = nullptr;
        pole.get_weighted_source_densities = nullptr;
        // P2Mソースの設定（Dunavant求積点）
        pole.p2m_sources.resize(n_p2m);
        for (int q = 0; q < n_p2m; ++q) {
          auto& isrc = pole.p2m_sources[q];
          const auto& [b0, b1, b2, ww] = p2m_rule[q];
          isrc.X = Dot(std::array<double, 3>{b0, b1, b2}, X012);
          isrc.normal = normal_vec; // 線形要素では法線は面全体で一定
          isrc.get_weighted_source_densities = [W = J_det * ww,
                                                N = std::array<double, 3>{b0, b1, b2},
                                                pp0 = pair_pointer_phiphin0,
                                                pp1 = pair_pointer_phiphin1,
                                                pp2 = pair_pointer_phiphin2]() -> std::array<double, 2> {
            return {W * (N[0] * (*pp0[0]) + N[1] * (*pp1[0]) + N[2] * (*pp2[0])),
                    W * (N[0] * (*pp0[1]) + N[1] * (*pp1[1]) + N[2] * (*pp2[1]))};
          };
        }
        // MAC拡張用: 重心から最も遠い内部ソースまでの距離を計算
        pole.max_source_offset = 0.0;
        for (const auto& isrc : pole.p2m_sources)
          pole.max_source_offset = std::max(pole.max_source_offset, Norm(isrc.X - pole.X));
      };

      if (reuse_sources) {
        setup_linear_element(*F->sources[0]);
      } else {
        auto pole = std::make_shared<source4FMM<target4FMM>>();
        setup_linear_element(*pole);
        F->sources.emplace_back(std::move(pole));
      }
    }
  }
  if (debug_line_factor2) {
    std::cout << Magenta << "[line:x2:summary] " << Cyan
              << "phi_checked=" << Green << dbg_phi_checked.load()
              << Cyan << " phi_near_x2=" << Red << dbg_phi_ratio2.load()
              << Cyan << " phin_checked=" << Green << dbg_phin_checked.load()
              << Cyan << " phin_near_x2=" << Red << dbg_phin_ratio2.load()
              << colorReset << std::endl;
  }
  const bool all_reused = (static_cast<std::size_t>(reuse_count.load()) == faces.size());
  std::cout << "  sources reuse: " << reuse_count.load() << "/" << faces.size() << " faces" << (all_reused ? " (all reused)" : " (some recreated)") << std::endl;
  return all_reused;
};

// -------------------------------------------------------------------------- */
// -------------------------------------------------------------------------- */

bool is_B_pole_initialized = false;
bool is_tree_generated = false;
bool is_fmm_static_initialized = false;
std::size_t fmm_tree_signature_built = 0;

std::size_t computeFMMTreeSignature() const {
  std::size_t h = 0;
  hash_combine(h, static_cast<std::size_t>(B_poles.max_level));
  hash_combine(h, static_cast<std::size_t>(B_poles.level_buckets.size()));
  for (const auto& lv : B_poles.level_buckets) {
    hash_combine(h, lv.size());
    for (const auto* b : lv)
      hash_combine(h, std::hash<const void*>{}(static_cast<const void*>(b)));
  }
  hash_combine(h, B_poles.deepest_level_buckets.size());
  for (const auto* b : B_poles.deepest_level_buckets)
    hash_combine(h, std::hash<const void*>{}(static_cast<const void*>(b)));
  return h;
}

void exportFMMTreeVisualization(const Network* obj) {
  const auto& dir = this->output_directory;
  {
    T8Tddd t8tddd = (CoordinateBounds)(obj->bounds);
    std::ofstream ofs(dir / "obj_bounds.vtp");
    vtkPolygonWrite(ofs, t8tddd);
    ofs.close();
  }

  std::vector<T8Tddd> cube_level_deepest;
  B_poles.forEachAtDeepest([&](auto* B) {
    T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
    cube_level_deepest.push_back(t8tddd); });
  std::ofstream ofs(dir / "cube_level_deepest_points.vtp");
  vtkPolygonWrite(ofs, cube_level_deepest);
  ofs.close();

  for (int i = 0; i <= B_poles.max_level; ++i) {
    std::vector<T8Tddd> cube_level;
    B_poles.forEachAtLevel(i, [&](auto* B) {
      T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
      cube_level.push_back(t8tddd); });
    std::ofstream ofs(dir / ("cube_level_points_" + std::to_string(i) + ".vtp"));
    vtkPolygonWrite(ofs, cube_level);
    ofs.close();
  }
  std::cout << Green << "[FMM] Tree visualization exported to " << dir << colorReset << std::endl;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

V_d solveSystemGMRES(bool is_time_derivative, V_d x0, bool allow_ilu_rebuild_retry = true, const std::function<void()>& postprocess = {}) {
  last_ilu_apply_time_sum = 0.0;
  last_gmres_iter_time_sum = 0.0;

  const V_d x0_initial = x0;
  copyToFMM(is_time_derivative); // phi,phitを選択してFMM用にコピー
  cacheBoundaryValues(points, x0.size());

  // Save midpoint per-face FMM maps before GMRES iteration.
  // The cache references point into phiOnFace_FMM/phinOnFace_FMM,
  // which get zeroed by setKnowns/setUnknowns. We restore after GMRES.
  std::vector<std::pair<std::unordered_map<networkFace*, double>, std::unordered_map<networkFace*, double>>> saved_midpoint_bc;
  if (use_true_quadratic_element) {
    saved_midpoint_bc.reserve(midpoint_lines.size());
    for (auto* l : midpoint_lines)
      saved_midpoint_bc.emplace_back(l->phiOnFace_FMM, l->phinOnFace_FMM);
  }

  this->setUnknowns(0.);
  if (postprocess)
    postprocess();
  auto b = -this->compute_Ax_minus_b(is_time_derivative);

  // コーナー点（DirichletとNeumannの交点）では，phi=phi の条件に変更
  _Pragma("omp parallel for") for (const auto& p : points) {
    if (p->CORNER)
      for (const auto& [f, i] : p->f2Index)
        if (isNeumannID_BEM(p, f))
          b[i] = is_time_derivative ? p->phiphin_t[0] : p->phiphin[0];
  }

  // CORNER midpoint constraint: phi_N = phi_D (same pattern as vertex CORNER above)
  if (use_true_quadratic_element) {
    for (auto* l : midpoint_lines) {
      if (l->CORNER && l->isMultipleNode) {
        for (const auto& [f, idx] : l->f2Index)
          if (isNeumannID_BEM(l, f))
            b[idx] = is_time_derivative ? l->phiphin_t[0] : l->phiphin[0];
      }
    }
  }

  if (!is_time_derivative) {
    for (const auto& p : points)
      for (const auto& [f, i] : p->f2Index)
        p->b_RHS_FMM = b[i];
  }

  // === Check 3: GMRES RHS breakdown (set BEM_BIE_CHECK=1 to enable) ===
  if (use_true_quadratic_element) {
    if (const char* env = std::getenv("BEM_BIE_CHECK"); env && std::string(env) != "0") {
      double max_b_vtx = 0, sum_b_vtx = 0;
      int n_vtx = 0;
      double max_b_mid = 0, sum_b_mid = 0;
      int n_mid = 0;
      for (const auto& p : points) {
        for (const auto& [f, i] : p->f2Index) {
          if (i >= 0 && i < static_cast<int>(b.size())) {
            max_b_vtx = std::max(max_b_vtx, std::abs(b[i]));
            sum_b_vtx += std::abs(b[i]);
            ++n_vtx;
          }
        }
      }
      for (auto* l : midpoint_lines) {
        for (const auto& [f, i] : l->f2Index) {
          if (i >= 0 && i < static_cast<int>(b.size())) {
            max_b_mid = std::max(max_b_mid, std::abs(b[i]));
            sum_b_mid += std::abs(b[i]);
            ++n_mid;
          }
        }
      }
      std::cout << "[BIE-Check3] GMRES RHS b:"
                << " vtx: max=" << max_b_vtx << " mean=" << (n_vtx > 0 ? sum_b_vtx / n_vtx : 0.)
                << " (n=" << n_vtx << ")"
                << " mid: max=" << max_b_mid << " mean=" << (n_mid > 0 ? sum_b_mid / n_mid : 0.)
                << " (n=" << n_mid << ")"
                << std::endl;
    }
  }

  const bool preconditioner_enabled = (preconditioner_type == "ILU" || preconditioner_type == "MILU" || preconditioner_type == "ILUT" || preconditioner_type == "SCHWARZ");
  const bool use_ilu0 = ((preconditioner_type == "ILU" || preconditioner_type == "MILU") && static_cast<bool>(ilu_preconditioner));
  const bool use_ilut = (preconditioner_type == "ILUT" && static_cast<bool>(ilut_preconditioner));
  const bool use_schwarz = (preconditioner_type == "SCHWARZ" && static_cast<bool>(schwarz_preconditioner));
  const bool use_jacobi = ((preconditioner_type == "ILU" || preconditioner_type == "MILU" || preconditioner_type == "ILUT") && !ilu_diag_inv_preconditioner.empty());

  auto apply_jacobi_preconditioner = [&](const V_d& rhs) -> V_d {
    TimeWatch watch;
    V_d out(rhs.size());
    _Pragma("omp parallel for") for (int i = 0; i < static_cast<int>(rhs.size()); ++i) {
      const double d = (i < static_cast<int>(ilu_diag_inv_preconditioner.size())) ? ilu_diag_inv_preconditioner[i] : 1.0;
      out[i] = rhs[i] * d;
    }
    last_ilu_apply_time_sum += watch()[0];
    return out;
  };

  // ==========================================================================
  // return_A_dot_v: GMRESに渡す行列ベクトル積関数
  // ==========================================================================
  // 役割: ベクトル V を受け取り、A*V を返す（GMRESのArnoldi過程で使用）
  //
  // GMRESでの使われ方:
  //   - 初期化時: r₀ = b - A*x₀ を計算
  //   - 各反復: v_{k+1} = A*v_k を計算してKrylov基底を構築
  //   - 真の残差: ||b - A*x|| を計算（リスタート時、収束判定に近い時のみ）
  //
  // 処理の流れ:
  //   1. setUnknowns(V): 未知数（Neumann面のφ、Dirichlet面のφn）を V に設定
  //   2. setKnowns(0.): 既知数を 0 に設定（b の寄与を消す → A*V のみ計算）
  //   3. compute_Ax_minus_b(): A*V - 0 = A*V を計算
  //   4. 前処理適用（オプション）: M⁻¹*A*V を返す
  //
  // 注意: setKnowns(0.) により、compute_Ax_minus_b() は実質 A*V を返す
  // ==========================================================================
  auto return_A_dot_v = [&](const V_d& V) -> V_d {
    this->setUnknowns(V);
    this->setKnowns(0.);
    if (postprocess)
      postprocess();

    // Compute A*V first, then (optionally) apply the preconditioner.
    // Important: keep timing accounting separated (FMM vs preconditioner apply).
    auto Av = this->compute_Ax_minus_b(is_time_derivative);

    if (use_schwarz) {
      TimeWatch watch;
      auto out = schwarz_preconditioner->solve(Av);
      last_ilu_apply_time_sum += watch()[0];
      return out;
    }
    if (use_ilut) {
      TimeWatch watch;
      auto out = ilut_preconditioner->solve(Av);
      last_ilu_apply_time_sum += watch()[0];
      return out;
    }
    if (use_ilu0) {
      TimeWatch watch;
      auto out = ilu_preconditioner->solve(Av);
      last_ilu_apply_time_sum += watch()[0];
      return out;
    }
    if (use_jacobi)
      return apply_jacobi_preconditioner(Av);
    return Av;
  };

  const V_d b_before_precond = b;
  if (use_schwarz) {
    TimeWatch watch;
    b = schwarz_preconditioner->solve(b);
    last_ilu_apply_time_sum += watch()[0];
  } else if (use_ilut) {
    TimeWatch watch;
    b = ilut_preconditioner->solve(b);
    last_ilu_apply_time_sum += watch()[0];
  } else if (use_ilu0) {
    TimeWatch watch;
    b = ilu_preconditioner->solve(b);
    last_ilu_apply_time_sum += watch()[0];
  } else if (use_jacobi)
    b = apply_jacobi_preconditioner(b);

  if (const char* env = std::getenv("BEM_GMRES_DEBUG_PRECOND"); env && std::string(env) != "0") {
    auto vec_stats = [](const V_d& v) {
      double norm2 = 0.0;
      double max_abs = 0.0;
      std::size_t finite_count = 0;
      for (double x : v) {
        if (std::isfinite(x)) {
          ++finite_count;
          norm2 += x * x;
          max_abs = std::max(max_abs, std::abs(x));
        }
      }
      return std::array<double, 3>{std::sqrt(norm2), max_abs, static_cast<double>(finite_count)};
    };
    auto print_top_entries = [&](const V_d& v, const char* label) {
      std::vector<std::pair<double, int>> abs_idx;
      abs_idx.reserve(v.size());
      for (int i = 0; i < static_cast<int>(v.size()); ++i) {
        const double x = v[i];
        if (!std::isfinite(x))
          continue;
        abs_idx.emplace_back(std::abs(x), i);
      }
      std::sort(abs_idx.begin(), abs_idx.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
      const int top_k = std::min<int>(10, static_cast<int>(abs_idx.size()));
      for (int rank = 0; rank < top_k; ++rank) {
        const int row_idx = abs_idx[rank].second;
        const double x = v[row_idx];
        const double alpha = (row_idx >= 0 && row_idx < static_cast<int>(diag_coeffs.size())) ? diag_coeffs[row_idx] : 0.0;
        std::cout << Magenta << "[GMRES:precond] " << Cyan
                  << label << " #" << rank + 1
                  << " row=" << Green << row_idx
                  << Cyan << " val=" << Yellow << x
                  << Cyan << " |val|=" << Yellow << std::abs(x)
                  << Cyan << " diag_coeff=" << Red << alpha
                  << colorReset << std::endl;
      }
    };

    const auto s0 = vec_stats(b_before_precond);
    const auto s1 = vec_stats(b);
    std::cout << Magenta << "[GMRES:precond] " << Cyan
              << "rhs_before: norm=" << Yellow << s0[0] << Cyan << " max=" << Yellow << s0[1] << Cyan << " finite=" << Yellow << static_cast<std::size_t>(s0[2]) << "/" << b_before_precond.size()
              << "  rhs_after: norm=" << Yellow << s1[0] << Cyan << " max=" << Yellow << s1[1] << Cyan << " finite=" << Yellow << static_cast<std::size_t>(s1[2]) << "/" << b.size()
              << colorReset << std::endl;
    print_top_entries(b_before_precond, "rhs_before");
    print_top_entries(b, "rhs_after");
  }

  const int restart_max = std::max(1, solver_restart);

  // FMM breakdown profiling (updateFMM vs integrateNearField vs integrateFarField).
  double fmm_update_time_sum = 0.0;
  double fmm_near_time_sum = 0.0;
  double fmm_far_time_sum = 0.0;
  std::size_t fmm_matvec_calls = 0;
  std::array<double, 6> fmm_update_step_time_sum{};
  profile_fmm_update_time_sum = &fmm_update_time_sum;
  profile_fmm_near_time_sum = &fmm_near_time_sum;
  profile_fmm_far_time_sum = &fmm_far_time_sum;
  profile_fmm_matvec_calls = &fmm_matvec_calls;
  profile_fmm_update_step_time_sum = &fmm_update_step_time_sum;

  struct run_result_t {
    V_d x;
    int iters = 0;
    double estimated_error = 1e300;
    double current_error = 1e300;
  };

  // Spatial mixed precision: self-cell and vertex-adjacent sources use double,
  // other near-field sources use float (GPU). This is set during setDirectIntegration.
  // Switch to full double when:
  //   1. Residual reaches switch_to_double_threshold (default 1e-4), OR
  //   2. Saturation is detected (residual stops decreasing)

  // Threshold-based switching: switch to double when residual reaches this level
  constexpr double switch_to_double_threshold = 1e-4;
  // Saturation detection parameters (fallback if threshold not reached)
  constexpr int saturation_check_window = 3;         // Number of restarts to check
  constexpr double saturation_ratio_threshold = 0.9; // Ratio threshold (1.0 = no improvement)

  // ==========================================================================
  // run_gmres: リスタートGMRESソルバー
  // ==========================================================================
  // 役割: 前処理付きリスタートGMRES法で Ax = b を解く
  //
  // GMRESアルゴリズムの概要:
  //   1. 初期残差 r₀ = b - A*x₀ を計算 → return_A_dot_v(x₀) - b
  //   2. Arnoldi過程でKrylov基底を構築:
  //      - 各反復で v_{k+1} = A*v_k を計算 → return_A_dot_v(v_k)
  //      - Gram-Schmidtで直交化、Hessenberg行列 H を構築
  //   3. 最小二乗問題 min ||β*e₁ - H*y|| を解く
  //   4. 解を更新: x = x₀ + V_k * y
  //   5. リスタート（k回反復後）: 新しい x₀ で 1. に戻る
  //
  // 残差の種類:
  //   - r_est (推定残差): Hessenberg行列から計算（コスト: ほぼゼロ）
  //   - r (真の残差): ||b - A*x|| を実際に計算（コスト: 1回の行列ベクトル積）
  //
  // ==========================================================================

  auto run_gmres = [&](V_d x_init) -> run_result_t {
    run_result_t out;
    out.x = std::move(x_init);

    while (out.iters < solver_max_iter && out.current_error > solver_tol) {
      TimeWatch watch_cycle;
      const int remaining = solver_max_iter - out.iters;
      const int k_cycle_max = std::min(restart_max, remaining);
      const int initial_k = 1; // Start with 1 for stability

      gmres GMRES(return_A_dot_v, b, out.x, initial_k);
      std::cout << Magenta << "  [GMRES] " << Cyan << "Initial error: " << Yellow << GMRES.beta << colorReset << std::endl;
      if (!std::isfinite(GMRES.beta) || GMRES.beta > 1e50) {
        std::cout << Red << "  [GMRES] Warning: very large preconditioned initial residual (beta=" << GMRES.beta
                  << "). Enable BEM_GMRES_DEBUG_PRECOND=1 and BEM_ILU0_DEBUG=1 for diagnostics." << colorReset << std::endl;
      }
      out.iters += initial_k;
      out.estimated_error = GMRES.err;

      int achieved_count = 0;
      int maybe_saturated_restarts = 0;
      while (true) {
        GMRES.Iterate(return_A_dot_v);
#if defined(USE_METAL_M2L)
        std::cout << Magenta << "  [GMRES] " << (isMetalM2LInUse() ? (Red + "Metal") : (Green + "CPU"));
#else
        std::cout << Magenta << "  [GMRES] " << Green << "CPU";
#endif
        std::cout << Cyan << "  Inner iter=" << Green << GMRES.n << Cyan << "  r_est=" << Yellow << std::scientific << std::setprecision(3) << GMRES.err << std::defaultfloat << colorReset << std::endl;
        ++out.iters;
        out.estimated_error = GMRES.err;
#if defined(USE_METAL_M2L)
        if (isMetalM2LInUse() && out.estimated_error <= switch_to_double_threshold) {
          // Switch to full double precision based on threshold
          ++achieved_count;
          if (achieved_count >= 5 || out.estimated_error <= 1E-6) {
            g_metal_m2l_active = false; // switch to full double precision
            GMRES.UpdateSolution();
            out.x = GMRES.x;
            out.current_error = Norm(return_A_dot_v(out.x) - b);
            std::cout << Magenta << "  [GMRES] " << Cyan << "  Switching to full double precision for FMM computations. r_true=" << Red << out.current_error << colorReset << std::endl;
            break; // floatの結果は倍精度のr_trueの結果を保証しなかった．そのため，ここで強制終了して，再度GMRESをリスタートする
          }
        }
#endif
        if (out.estimated_error <= solver_tol * 10.) {
          //check true residual
          GMRES.UpdateSolution();
          out.x = GMRES.x;
          out.current_error = Norm(return_A_dot_v(out.x) - b);
          std::cout << Magenta << "  [GMRES] " << Cyan << "  True residual check: r_true=" << Red << out.current_error << colorReset << std::endl;
          if (out.current_error <= solver_tol)
            break;
          if (maybe_saturated_restarts++ >= 5 /*これは，r_estがsolver_tol*10を超えた回数*/)
            break; // Saturation detected
        }
        if (out.iters >= solver_max_iter && GMRES.n >= static_cast<std::size_t>(k_cycle_max)) {
          // Max iterations reached
          GMRES.UpdateSolution();
          out.x = GMRES.x;
          if (out.estimated_error <= solver_tol * 10.)
            out.current_error = Norm(return_A_dot_v(out.x) - b);
          else
            out.current_error = out.estimated_error;

          break;
        }
      }

      const std::string& res_col = (out.current_error <= solver_tol) ? Green : (out.current_error <= solver_tol * 10.0 ? Yellow : Red);
      std::cout << "  " << Magenta << "[GMRES] " << Cyan << "iter=" << Green << out.iters << Cyan << "  r_est=" << Yellow << std::scientific << std::setprecision(3) << out.estimated_error << Cyan << "  r=" << res_col << out.current_error << std::defaultfloat << colorReset << std::endl;
      last_gmres_iter_time_sum += watch_cycle()[0];
      last_converged_k = GMRES.n;
    }

#if defined(USE_METAL_M2L)
    g_metal_m2l_active = true; // 次回のために復元
#endif

    return out;
  };

  /* --------------------------------------------------------- */

  TimeWatch watch_gmres;
  const auto result = run_gmres(x0);
  const double gmres_time = watch_gmres()[0];

  x0 = result.x;
  const int total_iter = result.iters;
  const double current_error = result.current_error;

  // Print summary
  const std::string bar = "========================================";
  std::cout << "\n"
            << Magenta << "[GMRES] " << Cyan << bar << colorReset << std::endl;

  const std::string& final_col = (current_error <= solver_tol) ? Green : Red;
  const std::string status = (current_error <= solver_tol) ? "CONVERGED" : "NOT CONVERGED";
  const double ms_per_iter = (total_iter > 0) ? (gmres_time * 1000.0 / total_iter) : 0.0;
  std::cout << Magenta << "[GMRES] " << Cyan << "iter=" << Green << total_iter
            << Cyan << "  time=" << Green << gmres_time << "s" << Cyan << " (" << Green << ms_per_iter << "ms/iter" << Cyan << ")"
            << "  r_est=" << final_col << result.estimated_error << Cyan << "  r_true=" << final_col << current_error
            << Cyan << "  [" << final_col << status << Cyan << "]" << colorReset << std::endl;
  std::cout << Magenta << "[GMRES] " << Cyan << bar << colorReset << "\n"
            << std::endl;

  profile_fmm_update_time_sum = nullptr;
  profile_fmm_near_time_sum = nullptr;
  profile_fmm_far_time_sum = nullptr;
  profile_fmm_matvec_calls = nullptr;
  profile_fmm_update_step_time_sum = nullptr;

  // Near/far times are sum of all thread-local times; divide by thread count to get wall-clock estimate
#ifdef _OPENMP
  const int nthreads = omp_get_max_threads();
#else
  const int nthreads = 1;
#endif
  const double fmm_near_time_wall = fmm_near_time_sum / nthreads;
  const double fmm_far_time_wall = fmm_far_time_sum / nthreads;
  const double fmm_total_time = fmm_update_time_sum + fmm_near_time_wall + fmm_far_time_wall;
  if (fmm_matvec_calls > 0 && fmm_total_time > 0.0) {
    const double update_pct = 100.0 * (fmm_update_time_sum / fmm_total_time);
    const double near_pct = 100.0 * (fmm_near_time_wall / fmm_total_time);
    const double far_pct = 100.0 * (fmm_far_time_wall / fmm_total_time);
    const double inv_calls = 1.0 / static_cast<double>(fmm_matvec_calls);
    std::cout << Magenta << "[FMM] " << Cyan << "stage=total"
              << "  calls=" << Green << fmm_matvec_calls << Cyan << "  update=" << Magenta << fmm_update_time_sum << "s" << Cyan << "(" << Magenta << (fmm_update_time_sum * inv_calls) << "s/call" << Cyan << ", " << Magenta << update_pct << "%" << Cyan << ")"
              << "  near=" << Yellow << fmm_near_time_wall << "s" << Cyan << "(" << Yellow << (fmm_near_time_wall * inv_calls) << "s/call" << Cyan << ", " << Yellow << near_pct << "%" << Cyan << ")"
              << "  far=" << Blue << fmm_far_time_wall << "s" << Cyan << "(" << Blue << (fmm_far_time_wall * inv_calls) << "s/call" << Cyan << ", " << Blue << far_pct << "%" << Cyan << ")" << colorReset << std::endl;
  }

  const double fmm_update_step_sum = std::accumulate(fmm_update_step_time_sum.begin(), fmm_update_step_time_sum.end(), 0.0);
  if (fmm_matvec_calls > 0 && fmm_update_step_sum > 0.0) {
    const double denom = (fmm_update_step_sum > 0.0) ? fmm_update_step_sum : 1.0;
    const double pct0 = 100.0 * (fmm_update_step_time_sum[0] / denom);
    const double pct1 = 100.0 * (fmm_update_step_time_sum[1] / denom);
    const double pct2 = 100.0 * (fmm_update_step_time_sum[2] / denom);
    const double pct3 = 100.0 * (fmm_update_step_time_sum[3] / denom);
    const double pct4 = 100.0 * (fmm_update_step_time_sum[4] / denom);
    const double pct5 = 100.0 * (fmm_update_step_time_sum[5] / denom);
    std::cout << Magenta << "[FMM:update] " << Cyan << "stage=total"
              << "  calls=" << Green << fmm_matvec_calls << Cyan << "  updateDensity=" << Green << fmm_update_step_time_sum[0] << "s" << Cyan << "(" << Green << pct0 << "%" << Cyan << ")"
              << "  reset=" << Green << fmm_update_step_time_sum[1] << "s" << Cyan << "(" << Green << pct1 << "%" << Cyan << ")"
              << "  P2M=" << Green << fmm_update_step_time_sum[2] << "s" << Cyan << "(" << Green << pct2 << "%" << Cyan << ")"
              << "  M2M=" << Green << fmm_update_step_time_sum[3] << "s" << Cyan << "(" << Green << pct3 << "%" << Cyan << ")"
              << "  M2L=" << Yellow << fmm_update_step_time_sum[4] << "s" << Cyan << "(" << Yellow << pct4 << "%" << Cyan << ")"
              << "  L2L=" << Green << fmm_update_step_time_sum[5] << "s" << Cyan << "(" << Green << pct5 << "%" << Cyan << ")" << colorReset << std::endl;
  }

  last_gmres_total_iter = total_iter;
  last_gmres_residual_norm = current_error;
  last_gmres_converged = (current_error <= solver_tol);
  last_gmres_hotspots.clear();
  if (!last_gmres_converged) {
    struct row_location_t {
      networkPoint* point = nullptr;
      networkLine* line = nullptr;
      networkFace* face = nullptr;
      bool is_midpoint = false;
    };

    std::vector<row_location_t> row_locations(matrix_size);
    for (auto* p : points) {
      for (const auto& [f, row_idx] : p->f2Index) {
        if (0 <= row_idx && row_idx < static_cast<int>(row_locations.size()))
          row_locations[row_idx] = row_location_t{p, nullptr, f, false};
      }
    }
    if (use_true_quadratic_element) {
      for (auto* l : midpoint_lines) {
        for (const auto& [f, row_idx] : l->f2Index) {
          if (0 <= row_idx && row_idx < static_cast<int>(row_locations.size()))
            row_locations[row_idx] = row_location_t{nullptr, l, f, true};
        }
      }
    }

    const auto residual_vec = return_A_dot_v(x0) - b;
    std::vector<int> rows(residual_vec.size());
    std::iota(rows.begin(), rows.end(), 0);
    auto residual_priority = [&](int row_idx) {
      const double r = residual_vec[row_idx];
      return std::isfinite(r) ? std::abs(r) : std::numeric_limits<double>::infinity();
    };
    std::sort(rows.begin(), rows.end(), [&](int a, int b) {
      return residual_priority(a) > residual_priority(b);
    });

    const int top_k = std::min<int>(8, static_cast<int>(rows.size()));
    last_gmres_hotspots.reserve(top_k);
    const auto face_flag_label = [](const networkFace* f) {
      if (!f)
        return std::string("-");
      std::string s;
      if (f->Dirichlet)
        s += "D";
      if (f->Neumann)
        s += "N";
      if (s.empty())
        s = "-";
      return s;
    };
    const auto make_face_flag_summary = [&](const auto& boundary_faces) {
      std::ostringstream oss;
      oss << "{";
      bool first = true;
      for (const auto* f : boundary_faces) {
        if (!first)
          oss << ",";
        first = false;
        oss << (f ? f->index : -1) << ":";
        oss << face_flag_label(f);
      }
      oss << "}";
      return oss.str();
    };
    const auto make_point_dof_summary = [&](const networkPoint* p) {
      std::ostringstream oss;
      oss << "{";
      bool first = true;
      for (const auto& [f, idx] : p->f2Index) {
        if (!first)
          oss << "; ";
        first = false;
        oss << "f=" << (f ? f->index : -1)
            << ":" << face_flag_label(f)
            << ",idx=" << idx
            << ",phi=" << getPhi(p, f)
            << ",phin=" << getPhin(p, f);
      }
      oss << "}";
      return oss.str();
    };
    const auto make_line_dof_summary = [&](const networkLine* l) {
      std::ostringstream oss;
      oss << "{";
      bool first = true;
      for (const auto& [f, idx] : l->f2Index) {
        if (!first)
          oss << "; ";
        first = false;
        oss << "f=" << (f ? f->index : -1)
            << ":" << face_flag_label(f)
            << ",idx=" << idx
            << ",phi=" << getPhi(l, f)
            << ",phin=" << getPhin(l, f);
      }
      oss << "}";
      return oss.str();
    };
    const auto make_point_contact_summary = [&](const networkPoint* p) {
      std::ostringstream oss;
      oss << "{";
      bool first = true;
      for (const auto* f : p->getBoundaryFaces()) {
        if (!first)
          oss << "; ";
        first = false;
        const auto [cf, Xc, dist] = p->getNearestContactFace_(f);
        oss << "adj=" << (f ? f->index : -1)
            << ":" << face_flag_label(f)
            << "->contact=" << (cf ? cf->index : -1)
            << ":" << face_flag_label(cf)
            << ",dist=" << dist
            << ",X=" << Xc;
        if (cf)
          oss << ",c=" << cf->centroid;
      }
      oss << "}";
      return oss.str();
    };
    const auto make_line_contact_summary = [&](const networkLine* l) {
      std::ostringstream oss;
      oss << "{";
      bool first = true;
      for (const auto* f : l->getBoundaryFaces()) {
        if (!first)
          oss << "; ";
        first = false;
        const auto [cf, Xc, dist] = l->getNearestContactFace_(f);
        oss << "adj=" << (f ? f->index : -1)
            << ":" << face_flag_label(f)
            << "->contact=" << (cf ? cf->index : -1)
            << ":" << face_flag_label(cf)
            << ",dist=" << dist
            << ",X=" << Xc;
        if (cf)
          oss << ",c=" << cf->centroid;
      }
      oss << "}";
      return oss.str();
    };
    for (int rank = 0; rank < top_k; ++rank) {
      const int row_idx = rows[rank];
      GMRESResidualHotspot hotspot;
      hotspot.row = row_idx;
      hotspot.residual = residual_vec[row_idx];
      hotspot.residual_abs = residual_priority(row_idx);
      hotspot.rhs = (0 <= row_idx && row_idx < static_cast<int>(b.size())) ? b[row_idx] : 0.0;
      hotspot.solution = (0 <= row_idx && row_idx < static_cast<int>(x0.size())) ? x0[row_idx] : 0.0;
      hotspot.diag_coeff = (0 <= row_idx && row_idx < static_cast<int>(diag_coeffs.size())) ? diag_coeffs[row_idx] : 0.0;

      const auto& loc = row_locations[row_idx];
      hotspot.is_midpoint = loc.is_midpoint;
      hotspot.face_index = loc.face ? loc.face->index : -1;
      if (loc.is_midpoint && loc.line) {
        hotspot.midpoint_index = loc.line->midpoint_index;
        hotspot.position = loc.line->X_mid;
        const auto [p0, p1] = loc.line->getPoints();
        hotspot.endpoint0 = p0 ? p0->X : Tddd{0., 0., 0.};
        hotspot.endpoint1 = p1 ? p1->X : Tddd{0., 0., 0.};
        hotspot.dirichlet = loc.line->Dirichlet;
        hotspot.neumann = loc.line->Neumann;
        hotspot.corner = loc.line->CORNER;
        hotspot.is_multiple_node = loc.line->isMultipleNode;
        hotspot.contact_faces_count = static_cast<int>(loc.line->ContactFaces.size());
        hotspot.adjacent_face_flags = make_face_flag_summary(loc.line->getBoundaryFaces());
        hotspot.dof_summary = make_line_dof_summary(loc.line);
        hotspot.contact_summary = make_line_contact_summary(loc.line);
      } else if (loc.point) {
        hotspot.point_index = loc.point->index;
        hotspot.position = loc.point->X;
        hotspot.dirichlet = loc.point->Dirichlet;
        hotspot.neumann = loc.point->Neumann;
        hotspot.corner = loc.point->CORNER;
        hotspot.is_multiple_node = loc.point->isMultipleNode;
        hotspot.contact_faces_count = static_cast<int>(loc.point->ContactFaces.size());
        hotspot.adjacent_face_flags = make_face_flag_summary(loc.point->getBoundaryFaces());
        hotspot.dof_summary = make_point_dof_summary(loc.point);
        hotspot.contact_summary = make_point_contact_summary(loc.point);
      }
      last_gmres_hotspots.push_back(hotspot);
    }
  }

  {
    std::string precond_kind = "NONE";
    if (preconditioner_type == "SCHWARZ" && static_cast<bool>(schwarz_preconditioner)) {
      precond_kind = "Schwarz";
    } else if ((preconditioner_type == "ILU" || preconditioner_type == "MILU" || preconditioner_type == "ILUT") && !ilu_diag_inv_preconditioner.empty()) {
      precond_kind = "Jacobi";
    } else if (preconditioner_type == "ILUT" && static_cast<bool>(ilut_preconditioner)) {
      precond_kind = "ILUT";
    } else if (preconditioner_type == "MILU" && static_cast<bool>(ilu_preconditioner)) {
      precond_kind = "MILU(0)";
    } else if (preconditioner_type == "ILU" && static_cast<bool>(ilu_preconditioner)) {
      precond_kind = "ILU(0)";
    }

    std::cout << Magenta << "[Precond] " << Cyan << "type=" << Green << precond_kind << Cyan;
    if (precond_kind == "ILUT") {
      std::cout << "  drop_tol=" << Yellow << ilut_drop_tol << Cyan << "  max_entries/row=" << Yellow << ilut_max_entries_per_row << Cyan << "  pivot_min=" << Yellow << ilut_pivot_min << Cyan;
    } else if (precond_kind == "MILU(0)") {
      std::cout << "  omega=" << Yellow << milu_omega << Cyan;
    } else if (precond_kind == "Schwarz" && static_cast<bool>(schwarz_preconditioner)) {
      std::cout << "  core_k=" << Yellow << schwarz_core_k << Cyan << "  overlap_k=" << Yellow << schwarz_overlap_k << Cyan << "  max_core=" << Yellow << schwarz_max_core_size << Cyan << "  max_block=" << Yellow << schwarz_max_block_size << Cyan << "  pivot_min=" << Yellow << schwarz_pivot_min << Cyan << "  diag_shift=" << Yellow << schwarz_diag_shift << Cyan << "  blocks=" << Yellow << schwarz_preconditioner->blocks_count << Cyan << "  avg_block=" << Yellow
                << schwarz_preconditioner->avg_block_size << Cyan << "  lu_blocks=" << Yellow << schwarz_preconditioner->num_lu_blocks << Cyan << "  fallback=" << Yellow << schwarz_preconditioner->num_fallback_blocks << Cyan;
    }
    std::cout << "  build=" << Magenta << last_ilu_build_time << "s" << Cyan << "  apply=" << Magenta << last_ilu_apply_time_sum << "s" << Cyan << "  gmres_iter_time=" << Magenta << last_gmres_iter_time_sum << "s" << colorReset << std::endl;
  }

  // Restore midpoint per-face FMM maps destroyed by setKnowns(0.) during GMRES.
  if (use_true_quadratic_element && !saved_midpoint_bc.empty()) {
    for (std::size_t k = 0; k < midpoint_lines.size(); ++k) {
      auto* l = midpoint_lines[k];
      l->phiOnFace_FMM = saved_midpoint_bc[k].first;
      l->phinOnFace_FMM = saved_midpoint_bc[k].second;
    }
  }

  copyToFMM(is_time_derivative);
  return x0;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

//   ここにILUを定義しよう

/* -------------------------------------------------------------------------- */

int total_unknowns = 0;
void updateGeometricProperties() {
  // 1. Update diagonal coefficients (geometry dependent)
  diag_coeffs.assign(total_unknowns, 1.0);
  cacheBoundaryValues(points, total_unknowns);

  // Save midpoint per-face FMM maps before rigid-mode test modifies them.
  std::vector<std::pair<std::unordered_map<networkFace*, double>, std::unordered_map<networkFace*, double>>> saved_midpoint_phi_phin;
  if (use_true_quadratic_element) {
    saved_midpoint_phi_phin.reserve(midpoint_lines.size());
    for (auto* l : midpoint_lines)
      saved_midpoint_phi_phin.emplace_back(l->phiOnFace_FMM, l->phinOnFace_FMM);
  }

  if (this->use_rigid_mode) {
    std::cout << "calculate diagonal coefficient (Rigid Mode)" << std::endl;
    _Pragma("omp parallel for") for (const auto& p : points) {
      for (const auto& [f, i] : p->f2Index) {
        p->phinOnFace_FMM[f] = 0.;
        p->phiOnFace_FMM[f] = 1.;
        // Keep dense caches consistent for any subsequent near-field evaluation.
        if (i >= 0 && i < static_cast<int>(cache_phi_val_D_by_index.size())) {
          cache_phi_val_D_by_index[i] = 1.0;
          cache_phin_val_D_by_index[i] = 0.0;
        }
      }
    }
    // Set midpoint phi=1, phin=0 for P2M and near-field consistency
    if (use_true_quadratic_element) {
      for (auto* l : midpoint_lines) {
        for (auto& [f, val] : l->phiOnFace_FMM)
          val = 1.0;
        for (auto& [f, val] : l->phinOnFace_FMM)
          val = 0.0;
        for (const auto& [f, mi] : l->f2Index) {
          if (mi >= 0 && mi < static_cast<int>(cache_phi_val_D_by_index.size())) {
            cache_phi_val_D_by_index[mi] = 1.0;
            cache_phin_val_D_by_index[mi] = 0.0;
          }
        }
      }
    }
    std::cout << Magenta << "[FMM] " << Cyan << "updateFMM (for diagonal coeff)" << colorReset << std::endl;
    std::array<double, 6> tmp_elapsed_time{};
    updateFMM(B_poles, tmp_elapsed_time);
    _Pragma("omp parallel for") for (const auto& a : points) {
      for (const auto& [f, i] : a->f2Index) {
        // This path is not performance critical (only during diagonal-coefficient setup),
        // but keep it consistent with the current dense caches.
        auto GPhin_GnPhi_near = integrateNearField(*a);
        auto GPhin_GnPhi_far = a->integrateFarField();
        auto [GPhin, GnPhi] = GPhin_GnPhi_near + GPhin_GnPhi_far;
        diag_coeffs[i] = -GnPhi;
        if (std::abs(diag_coeffs[i]) < 1e-15)
          diag_coeffs[i] = 1.0;
      }
    }
  } else {
    std::cout << "calculate diagonal coefficient (Solid Angle)" << std::endl;
    _Pragma("omp parallel for") for (const auto& a : points) {
      for (const auto& [f, i] : a->f2Index) {
        diag_coeffs[i] = a->getSolidAngle();
      }
    }
  }

  // Compute diagonal coefficients for midpoint DOFs (true quadratic elements)
  if (use_true_quadratic_element) {
    if (this->use_rigid_mode) {
      // Rigid mode: set phi=1, phin=0 on all DOFs, evaluate BIE at midpoints
      setKnowns(0.0);
      setUnknowns(0.0);
      // Set all phi=1 for rigid body test
      _Pragma("omp parallel for") for (const auto& p : points) {
        for (const auto& [f, i] : p->f2Index) {
          p->phiOnFace_FMM[f] = 1.0;
          p->phinOnFace_FMM[f] = 0.0;
          if (i >= 0 && i < static_cast<int>(cache_phi_val_D_by_index.size())) {
            cache_phi_val_D_by_index[i] = 1.0;
            cache_phin_val_D_by_index[i] = 0.0;
          }
        }
      }
      for (auto* l : midpoint_lines) {
        for (auto& [f, val] : l->phiOnFace_FMM)
          val = 1.0;
        for (auto& [f, val] : l->phinOnFace_FMM)
          val = 0.0;
        for (const auto& [f, mi] : l->f2Index) {
          if (mi >= 0 && mi < static_cast<int>(cache_phi_val_D_by_index.size())) {
            cache_phi_val_D_by_index[mi] = 1.0;
            cache_phin_val_D_by_index[mi] = 0.0;
          }
        }
      }
      std::array<double, 6> tmp_elapsed_time{};
      updateFMM(B_poles, tmp_elapsed_time);
      _Pragma("omp parallel for") for (int k = 0; k < static_cast<int>(midpoint_lines.size()); ++k) {
        auto* l = midpoint_lines[k];
        auto& mt = static_cast<target4FMM&>(*l);
        auto GPhin_GnPhi_near = integrateNearField(mt);
        auto GPhin_GnPhi_far = mt.integrateFarField();
        auto [GPhin, GnPhi] = GPhin_GnPhi_near + GPhin_GnPhi_far;
        // BIE-derived diag_coeff applies to primary DOF; CORNER Neumann DOF gets 1.0 (constraint row)
        for (const auto& [f, i] : l->f2Index) {
          if (l->CORNER && isNeumannID_BEM(l, f)) {
            diag_coeffs[i] = 1.0; // constraint row diagonal
          } else {
            diag_coeffs[i] = -GnPhi;
            if (std::abs(diag_coeffs[i]) < 1e-15)
              diag_coeffs[i] = 0.5;
          }
        }
      }
    } else {
      // Solid angle at midpoint: 0.5 for smooth surface
      for (auto* l : midpoint_lines) {
        for (const auto& [f, i] : l->f2Index) {
          if (i >= 0 && i < static_cast<int>(diag_coeffs.size())) {
            if (l->CORNER && isNeumannID_BEM(l, f))
              diag_coeffs[i] = 1.0; // constraint row
            else
              diag_coeffs[i] = 0.5;
          }
        }
      }
    }
  }

  // === BIE Diagnostic Checks (set BEM_BIE_CHECK=1 to enable) ===
  if (use_true_quadratic_element && this->use_rigid_mode) {
    if (const char* env = std::getenv("BEM_BIE_CHECK"); env && std::string(env) != "0") {
      // --- Check 1: Rigid-mode full residual ---
      // State is phi=1, phin=0 on all DOFs; FMM just updated.
      // GPhin should be ~0 (phin=0 everywhere).
      // Residual = GPhin - (GnPhi + diag*1) should be ~0.
      {
        double max_res_vtx = 0, max_res_mid = 0;
        double max_GPhin_vtx = 0, max_GPhin_mid = 0;
        for (const auto& a : points) {
          for (const auto& [f, i] : a->f2Index) {
            if (a->CORNER && isNeumannID_BEM(a, f))
              continue;
            auto [GPhin, GnPhi] = integrateNearField(*a) + a->integrateFarField();
            double res = GPhin - (GnPhi + diag_coeffs[i] * 1.0);
            max_res_vtx = std::max(max_res_vtx, std::abs(res));
            max_GPhin_vtx = std::max(max_GPhin_vtx, std::abs(GPhin));
          }
        }
        for (int k = 0; k < static_cast<int>(midpoint_lines.size()); ++k) {
          auto* l = midpoint_lines[k];
          auto& mt = static_cast<target4FMM&>(*l);
          int i = l->midpoint_index; // primary DOF index for BIE check
          auto [GPhin, GnPhi] = integrateNearField(mt) + mt.integrateFarField();
          double res = GPhin - (GnPhi + diag_coeffs[i] * 1.0);
          max_res_mid = std::max(max_res_mid, std::abs(res));
          max_GPhin_mid = std::max(max_GPhin_mid, std::abs(GPhin));
        }
        std::cout << "[BIE-Check1] Rigid-mode residual:"
                  << " max|res| vtx=" << max_res_vtx << " mid=" << max_res_mid
                  << " max|GPhin| vtx=" << max_GPhin_vtx << " mid=" << max_GPhin_mid
                  << std::endl;
      }
      // --- Check 2: G*phin pathway (phi=0, phin=1) ---
      // Tests the G kernel path which rigid-mode (phin=0) does NOT test.
      {
        for (auto& p : points) {
          for (auto& [f, val] : p->phiOnFace_FMM)
            val = 0.0;
          for (auto& [f, val] : p->phinOnFace_FMM)
            val = 1.0;
          for (const auto& [f, i] : p->f2Index) {
            if (i >= 0 && i < static_cast<int>(cache_phi_val_D_by_index.size())) {
              cache_phi_val_D_by_index[i] = 0.0;
              cache_phin_val_D_by_index[i] = 1.0;
            }
          }
        }
        for (auto* l : midpoint_lines) {
          for (auto& [f, val] : l->phiOnFace_FMM)
            val = 0.0;
          for (auto& [f, val] : l->phinOnFace_FMM)
            val = 1.0;
          for (const auto& [f, mi] : l->f2Index) {
            if (mi >= 0 && mi < static_cast<int>(cache_phi_val_D_by_index.size())) {
              cache_phi_val_D_by_index[mi] = 0.0;
              cache_phin_val_D_by_index[mi] = 1.0;
            }
          }
        }
        std::array<double, 6> tmp_check2{};
        updateFMM(B_poles, tmp_check2);

        int n_sample = std::min(10, static_cast<int>(midpoint_lines.size()));
        for (int k = 0; k < n_sample; ++k) {
          auto* l = midpoint_lines[k];
          auto& mt = static_cast<target4FMM&>(*l);
          auto [GPhin_near_m, GnPhi_near_m] = integrateNearField(mt);
          auto [GPhin_far_m, GnPhi_far_m] = mt.integrateFarField();
          auto [p0, p1] = l->getPoints();
          auto [GPhin0, GnPhi0] = integrateNearField(*p0) + p0->integrateFarField();
          auto [GPhin1, GnPhi1] = integrateNearField(*p1) + p1->integrateFarField();
          double mid_total = GPhin_near_m + GPhin_far_m;
          double avg_vtx = 0.5 * (GPhin0 + GPhin1);
          std::cout << "[BIE-Check2] mid[" << k << "]"
                    << " near=" << GPhin_near_m << " far=" << GPhin_far_m
                    << " total=" << mid_total
                    << " avg_vtx=" << avg_vtx
                    << " ratio=" << (std::abs(avg_vtx) > 1e-30 ? mid_total / avg_vtx : 0.)
                    << std::endl;
        }
      }
    }
  }

  // Restore midpoint per-face FMM maps after rigid-mode test
  if (use_true_quadratic_element && !saved_midpoint_phi_phin.empty()) {
    for (std::size_t k = 0; k < midpoint_lines.size(); ++k) {
      midpoint_lines[k]->phiOnFace_FMM = saved_midpoint_phi_phin[k].first;
      midpoint_lines[k]->phinOnFace_FMM = saved_midpoint_phi_phin[k].second;
    }
  }

  // Store diagonal coefficient to nodes for existing ParaView output.
  for (const auto& p : points) {
    double alpha_sum = 0.0;
    int alpha_count = 0;
    for (const auto& [f, i] : p->f2Index) {
      (void)f;
      if (i >= 0 && i < static_cast<int>(diag_coeffs.size())) {
        alpha_sum += diag_coeffs[i];
        ++alpha_count;
      }
    }
    p->diag_coeff_BEM = (alpha_count > 0) ? (alpha_sum / static_cast<double>(alpha_count)) : 0.0;
  }
  if (use_true_quadratic_element) {
    for (auto* l : midpoint_lines) {
      // Store primary DOF's diag_coeff for ParaView output
      const int i = l->midpoint_index;
      l->diag_coeff_BEM = (i >= 0 && i < static_cast<int>(diag_coeffs.size())) ? diag_coeffs[i] : 0.0;
    }
  }

  if (const char* env = std::getenv("BEM_DIAG_COEFF_DEBUG"); env && std::string(env) != "0") {
    struct RowMeta {
      bool valid = false;
      bool is_midpoint = false;
      bool is_neumann = false;
      bool is_corner = false;
      Tddd X = {0., 0., 0.};
    };

    std::vector<RowMeta> row_meta(total_unknowns);
    for (const auto& p : points) {
      for (const auto& [f, row_idx] : p->f2Index) {
        if (row_idx < 0 || row_idx >= total_unknowns)
          continue;
        auto& m = row_meta[row_idx];
        if (m.valid)
          continue;
        m.valid = true;
        m.is_midpoint = false;
        m.is_neumann = isNeumannID_BEM(p, f);
        m.is_corner = p->CORNER;
        m.X = p->Xtarget;
      }
    }
    if (use_true_quadratic_element) {
      for (auto* l : midpoint_lines) {
        for (const auto& [f, row_idx] : l->f2Index) {
          if (row_idx < 0 || row_idx >= total_unknowns)
            continue;
          auto& m = row_meta[row_idx];
          m.valid = true;
          m.is_midpoint = true;
          m.is_neumann = isNeumannID_BEM(l, f);
          m.is_corner = l->CORNER;
          m.X = l->Xtarget;
        }
      }
    }

    std::size_t finite_count = 0, non_finite_count = 0;
    double min_abs = std::numeric_limits<double>::infinity();
    double max_abs = 0.0;
    int min_abs_row = -1, max_abs_row = -1;
    std::size_t gt1 = 0, gt10 = 0, gt100 = 0, gt1000 = 0;
    std::vector<std::pair<double, int>> abs_rows;
    abs_rows.reserve(diag_coeffs.size());
    for (int i = 0; i < static_cast<int>(diag_coeffs.size()); ++i) {
      const double a = diag_coeffs[i];
      if (!std::isfinite(a)) {
        ++non_finite_count;
        continue;
      }
      ++finite_count;
      const double aa = std::abs(a);
      abs_rows.emplace_back(aa, i);
      if (aa < min_abs) {
        min_abs = aa;
        min_abs_row = i;
      }
      if (aa > max_abs) {
        max_abs = aa;
        max_abs_row = i;
      }
      if (aa > 1.0)
        ++gt1;
      if (aa > 10.0)
        ++gt10;
      if (aa > 100.0)
        ++gt100;
      if (aa > 1000.0)
        ++gt1000;
    }
    std::sort(abs_rows.begin(), abs_rows.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

    std::cout << Magenta << "[diag_coeffs] " << Cyan
              << "size=" << Green << diag_coeffs.size()
              << Cyan << " finite=" << Green << finite_count
              << Cyan << " non_finite=" << Red << non_finite_count
              << Cyan << " min|a|=" << Yellow << min_abs << Cyan << "(row=" << Yellow << min_abs_row << Cyan << ")"
              << " max|a|=" << Red << max_abs << Cyan << "(row=" << Red << max_abs_row << Cyan << ")"
              << " count(|a|>1)=" << Yellow << gt1
              << Cyan << " >10=" << Yellow << gt10
              << Cyan << " >100=" << Yellow << gt100
              << Cyan << " >1000=" << Red << gt1000
              << colorReset << std::endl;

    const int top_k = std::min<int>(20, static_cast<int>(abs_rows.size()));
    for (int rank = 0; rank < top_k; ++rank) {
      const int row_idx = abs_rows[rank].second;
      const double a = diag_coeffs[row_idx];
      const auto& m = (row_idx >= 0 && row_idx < static_cast<int>(row_meta.size())) ? row_meta[row_idx] : RowMeta{};
      std::cout << Magenta << "[diag_coeffs] " << Cyan
                << "#" << rank + 1
                << " row=" << Green << row_idx
                << Cyan << " alpha=" << Red << a
                << Cyan << " kind=" << (m.is_midpoint ? Yellow + std::string("mid") : Green + std::string("vertex"))
                << Cyan << " bc=" << (m.is_neumann ? Red + std::string("N") : Blue + std::string("D"))
                << Cyan << " corner=" << (m.is_corner ? "1" : "0")
                << Cyan << " X=(" << m.X[0] << "," << m.X[1] << "," << m.X[2] << ")"
                << colorReset << std::endl;
    }
  }

  // 2. Rebuild preconditioner if selected
  if (preconditioner_type == "ILU" || preconditioner_type == "MILU" || preconditioner_type == "ILUT" || preconditioner_type == "SCHWARZ") {
    const int k_ring_num = std::clamp(ilu_kring_num, 0, 20);
    const std::size_t geom_hash = computeILUGeometrySignature();
    const bool have_any_precond = (static_cast<bool>(ilu_preconditioner) || static_cast<bool>(ilut_preconditioner) || static_cast<bool>(schwarz_preconditioner) || !ilu_diag_inv_preconditioner.empty());
    const bool same_settings = (preconditioner_type_built == preconditioner_type && ilu_neighborhood_type_built == ilu_neighborhood_type && ilu_k_ring_num_built == k_ring_num);
    const bool same_ilut_params = (preconditioner_type != "ILUT") || (ilut_drop_tol_built == ilut_drop_tol && ilut_max_entries_per_row_built == ilut_max_entries_per_row && ilut_pivot_min_built == ilut_pivot_min);
    const bool same_schwarz_params = (preconditioner_type != "SCHWARZ") || (schwarz_core_k_built == schwarz_core_k && schwarz_overlap_k_built == schwarz_overlap_k && schwarz_max_core_size_built == schwarz_max_core_size && schwarz_max_block_size_built == schwarz_max_block_size && schwarz_pivot_min_built == schwarz_pivot_min && schwarz_diag_shift_built == schwarz_diag_shift);
    const bool same_matrix = (ilu_matrix_size_built == static_cast<int>(matrix_size));
    const bool same_geometry = (ilu_geometry_hash_cached == geom_hash);

    if (!(have_any_precond && same_settings && same_ilut_params && same_schwarz_params && same_matrix && same_geometry)) {
      ilu_preconditioner.reset();
      ilut_preconditioner.reset();
      schwarz_preconditioner.reset();
      ilu_diag_inv_preconditioner.clear();
      buildSparseMatrixForILU(); // may build ILU(0) or diag preconditioner (k-ring num=0)
    } else {
      last_ilu_build_time = 0.0;
    }
  } else {
    // Do not keep stale preconditioners around when not requested.
    ilu_preconditioner.reset();
    ilut_preconditioner.reset();
    schwarz_preconditioner.reset();
    ilu_diag_inv_preconditioner.clear();
    ilu_faces_to_integrate_by_point.clear();
    ilu_faces_to_integrate_by_midpoint.clear();
    ilu_k_ring_num_cached = 0;
    ilu_topology_hash_cached = 0;
    ilu_geometry_hash_cached = 0;
    preconditioner_type_built.clear();
    ilu_neighborhood_type_built.clear();
    ilu_k_ring_num_built = -1;
    ilu_matrix_size_built = -1;
    ilut_drop_tol_built = -1.0;
    ilut_max_entries_per_row_built = -1;
    ilut_pivot_min_built = -1.0;
    schwarz_core_k_built = -1;
    schwarz_overlap_k_built = -1;
    schwarz_max_core_size_built = -1;
    schwarz_max_block_size_built = -1;
    schwarz_pivot_min_built = -1.0;
    schwarz_diag_shift_built = -1.0;
  }
}

/* -------------------------------------------------------------------------- */

void solveGMRES(TimeWatch& watch, double& time_setup, double& time_solve, const std::function<void()>& postprocess = {}) {
  (void)watch;
  std::cout << Magenta << "[Solver] " << Cyan << "GMRES" << colorReset << std::endl;
  TimeWatch watch_setup;
  TimeWatch log_stage_watch;
  auto log_stage = [&](const std::string& label) {
    std::cout << Magenta << "[BEM] " << Cyan << label << Green << " elapsed=" << log_stage_watch() << colorReset << std::endl;
  };

  //@ -------------------------------------------------------------------------- */
  //@   Geometry update (Network-managed) + topology hash for tree reuse         */
  //@ -------------------------------------------------------------------------- */

  auto obj = WATERS[0];
  this->points = obj->getBoundaryPoints();
  auto faces = obj->getBoundaryFaces();

  //@ -------------------------------------------------------------------------- */
  //@   Coordinate Scaling for FMM numerical stability                           */
  //@ -------------------------------------------------------------------------- */
  use_coordinate_scaling_ = true;
  obj->setGeometricProperties();
  coordinate_scale_factor_ = obj->computeCharacteristicLength();
  if (coordinate_scale_factor_ > 1e-10) {
    std::stringstream ss;
    ss << "Coordinate scaling enabled (L=" << coordinate_scale_factor_ << ")" << colorReset << std::endl;
    log_stage(ss.str());
    for (auto& net : this->WATERS)
      net->applyScaling(coordinate_scale_factor_);
  } else {
    std::stringstream ss;
    ss << "Coordinate scaling skipped (L too small: " << coordinate_scale_factor_ << ")" << colorReset << std::endl;
    log_stage(ss.str());
    use_coordinate_scaling_ = false;
    coordinate_scale_factor_ = 1.0;
  }

  // Midpoint Neumann BC scaling (∂φ/∂n' = L × ∂φ/∂n) is handled by copyToFMM(),
  // which applies phin_scale to the FMM copy. No pre-scaling of phinOnFace needed.

  for (auto& net : this->WATERS)
    net->setGeometricProperties();
  log_stage("setGeometricProperties");

  // Optimized parallel setIntegrationInfo: collect all faces first, then parallelize once
  std::vector<networkFace*> all_faces;
  for (const auto& water : WATERS) {
    const auto& faces = water->getBoundaryFaces();
    all_faces.insert(all_faces.end(), faces.begin(), faces.end());
  }
  _Pragma("omp parallel for schedule(dynamic, 100)") for (std::size_t i = 0; i < all_faces.size(); ++i) all_faces[i]->setIntegrationInfo();
  log_stage("setIntegrationInfo");

  // Full topology hash (faces + connectivity) for tree reuse decisions.
  static std::size_t last_faces_hash = 0;
  static std::size_t last_points_hash = 0;

  std::size_t faces_hash = 0;
  for (const auto& f : faces) {
    faces_hash ^= std::hash<const networkFace*>{}(f);
    auto [p0, p1, p2] = f->getPoints();
    faces_hash ^= std::hash<const networkPoint*>{}(p0);
    faces_hash ^= std::hash<const networkPoint*>{}(p1);
    faces_hash ^= std::hash<const networkPoint*>{}(p2);
  }
  std::size_t points_hash = 0;
  for (const auto& p : this->points)
    points_hash ^= std::hash<const networkPoint*>{}(p);

  const bool topology_same = is_tree_generated && (faces_hash == last_faces_hash) && (points_hash == last_points_hash);
  {
    std::stringstream ss;
    ss << "tree reuse check topology_same=" << topology_same;
    if (!topology_same) {
      ss << " reason:";
      if (!is_tree_generated)
        ss << " tree_not_generated";
      if (faces_hash != last_faces_hash)
        ss << " faces_changed";
      if (points_hash != last_points_hash)
        ss << " points_changed";
    }
    ss << " faces_hash=" << faces_hash << " points_hash=" << points_hash;
    log_stage(ss.str());
  }

  //@ -------------------------------------------------------------------------- */
  //@                       バケットの作成．極の追加 add                         */
  //@ -------------------------------------------------------------------------- */

  if (is_B_pole_initialized == false) {
    this->initializeBucket();
    is_B_pole_initialized = true;
  }
  log_stage("initializeBucket (first time only)");

  // Reorder unknown indices (Morton/Z-order) to improve near-field locality (helps SIMD and cache).
  // Disable with: `export BEM_REINDEX_UNKNOWN_MORTON=0`
  bool enable_reindex_morton = true;
  if (const char* env = std::getenv("BEM_REINDEX_UNKNOWN_MORTON"))
    enable_reindex_morton = (std::string(env) != "0");
  if (enable_reindex_morton) {
    struct dof_entry {
      std::uint64_t key;
      networkPoint* p;
      networkFace* f;
    };

    auto part1by2 = [](std::uint32_t x) -> std::uint64_t {
      std::uint64_t v = x & 0x1fffffU; // 21 bits
      v = (v | (v << 32)) & 0x1f00000000ffffULL;
      v = (v | (v << 16)) & 0x1f0000ff0000ffULL;
      v = (v | (v << 8)) & 0x100f00f00f00f00fULL;
      v = (v | (v << 4)) & 0x10c30c30c30c30c3ULL;
      v = (v | (v << 2)) & 0x1249249249249249ULL;
      return v;
    };

    auto morton3 = [&](std::uint32_t x, std::uint32_t y, std::uint32_t z) -> std::uint64_t { return (part1by2(x) | (part1by2(y) << 1) | (part1by2(z) << 2)); };

    // Bounds for quantization (parallel reduction)
    double xmin = 1e300, xmax = -1e300;
    double ymin = 1e300, ymax = -1e300;
    double zmin = 1e300, zmax = -1e300;

    const auto& pts = this->points;
    _Pragma("omp parallel for reduction(min:xmin,ymin,zmin) reduction(max:xmax,ymax,zmax)") for (std::size_t i = 0; i < pts.size(); ++i) {
      const auto& x = pts[i]->X;
      const double px = std::get<0>(x);
      const double py = std::get<1>(x);
      const double pz = std::get<2>(x);
      xmin = std::min(xmin, px);
      xmax = std::max(xmax, px);
      ymin = std::min(ymin, py);
      ymax = std::max(ymax, py);
      zmin = std::min(zmin, pz);
      zmax = std::max(zmax, pz);
    }

    const int bits = 21;
    const std::uint32_t maxv = (1u << bits) - 1u;
    auto quantize = [&](double v, double v0, double v1) -> std::uint32_t {
      if (!(v1 > v0))
        return 0u;
      double t = (v - v0) / (v1 - v0);
      t = std::clamp(t, 0.0, 1.0);
      return static_cast<std::uint32_t>(t * static_cast<double>(maxv) + 0.5);
    };

    // Parallel DOF array creation (avoiding race on push_back)
    // First, compute offsets for each point
    std::vector<std::size_t> point_dof_counts(this->points.size());
    _Pragma("omp parallel for") for (std::size_t i = 0; i < this->points.size(); ++i) { point_dof_counts[i] = this->points[i]->f2Index.size(); }

    // Prefix sum to get offsets
    std::vector<std::size_t> point_offsets(this->points.size() + 1, 0);
    for (std::size_t i = 0; i < this->points.size(); ++i) {
      point_offsets[i + 1] = point_offsets[i] + point_dof_counts[i];
    }
    std::size_t total = point_offsets.back();

    std::vector<dof_entry> dofs(total);

    // Parallel: compute Morton codes and fill DOF array
    _Pragma("omp parallel for") for (std::size_t i = 0; i < this->points.size(); ++i) {
      auto* p = this->points[i];
      const auto& x = p->X;
      const std::uint32_t qx = quantize(std::get<0>(x), xmin, xmax);
      const std::uint32_t qy = quantize(std::get<1>(x), ymin, ymax);
      const std::uint32_t qz = quantize(std::get<2>(x), zmin, zmax);
      const std::uint64_t key = morton3(qx, qy, qz);

      std::size_t local_idx = point_offsets[i];
      for (const auto& [f, old] : p->f2Index) {
        (void)old;
        dofs[local_idx++] = {key, p, f};
      }
    }

    std::stable_sort(dofs.begin(), dofs.end(), [](const dof_entry& a, const dof_entry& b) {
      if (a.key != b.key)
        return a.key < b.key;
      if (a.p != b.p)
        return a.p < b.p;
      return a.f < b.f;
    });

    int new_idx = 0;
    for (const auto& e : dofs) {
      e.p->f2Index[e.f] = new_idx++;
    }

    // Also reindex midpoint DOFs (true quadratic) after vertex DOFs
    if (use_true_quadratic_element) {
      // Collect all boundary lines with assigned midpoint DOFs
      std::unordered_set<networkLine*> unique_lines;
      for (auto* water : WATERS)
        for (auto* l : water->getBoundaryLines())
          if (!l->f2Index.empty())
            unique_lines.emplace(l);

      if (!unique_lines.empty()) {
        // Extend bounding box to include midpoints
        for (auto* l : unique_lines) {
          auto [pA, pB] = l->getPoints();
          double mx = 0.5 * (std::get<0>(pA->X) + std::get<0>(pB->X));
          double my = 0.5 * (std::get<1>(pA->X) + std::get<1>(pB->X));
          double mz = 0.5 * (std::get<2>(pA->X) + std::get<2>(pB->X));
          xmin = std::min(xmin, mx);
          xmax = std::max(xmax, mx);
          ymin = std::min(ymin, my);
          ymax = std::max(ymax, my);
          zmin = std::min(zmin, mz);
          zmax = std::max(zmax, mz);
        }

        struct mid_dof_entry {
          std::uint64_t key;
          networkLine* l;
        };
        std::vector<mid_dof_entry> mid_dofs;
        mid_dofs.reserve(unique_lines.size());
        for (auto* l : unique_lines) {
          auto [pA, pB] = l->getPoints();
          double mx = 0.5 * (std::get<0>(pA->X) + std::get<0>(pB->X));
          double my = 0.5 * (std::get<1>(pA->X) + std::get<1>(pB->X));
          double mz = 0.5 * (std::get<2>(pA->X) + std::get<2>(pB->X));
          mid_dofs.push_back({morton3(quantize(mx, xmin, xmax), quantize(my, ymin, ymax), quantize(mz, zmin, zmax)), l});
        }
        std::stable_sort(mid_dofs.begin(), mid_dofs.end(), [](const mid_dof_entry& a, const mid_dof_entry& b) {
          return a.key < b.key;
        });
        for (const auto& e : mid_dofs) {
          // Reindex all DOFs of this midpoint (primary DOF first, then per-face DOFs)
          // Preserve the face→index mapping structure but assign new contiguous indices
          e.l->midpoint_index = new_idx; // backward compat: primary DOF index
          // Reindex: nullptr first (Dirichlet/primary DOF), then face-specific DOFs
          if (e.l->f2Index.count(nullptr))
            e.l->f2Index[nullptr] = new_idx++;
          for (auto& [f, idx] : e.l->f2Index) {
            if (f != nullptr)
              idx = new_idx++;
          }
        }
      }
    }

    std::cout << Magenta << "[Solver] " << Cyan << "reindexed unknowns (Morton order), total=" << Green << new_idx << colorReset << std::endl;
  }
  log_stage("reindex unknowns (Morton)");
  copyToFMM(false);

  const bool reused_sources = this->createSourcesOnSurfaces();
  std::cout << Magenta << "[BEM] " << Cyan << "sources reuse" << Green << " all_reused=" << reused_sources << colorReset << std::endl;
  log_stage("createSourcesOnSurfaces");

  if (topology_same && reused_sources) {
    std::cout << "updateTree (rebin)" << std::endl;
    TimeWatch tw;
    B_poles.rebin();
    std::cout << Magenta << "updateTree" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
  } else {
    is_fmm_static_initialized = false;
    // 既存ツリーを保持したまま `add()` すると、各追加でツリー全体を走査してしまい重くなる。
    // ここでは一旦ツリーを破棄してからフラットに追加し、最後に条件に従ってツリーを再生成する。
    B_poles.traverseTree([](auto& child) {
      delete child;
      child = nullptr; });
    B_poles.clearDataKeepTree();

    /* -------------------------------------------------------------------------- */

    std::cout << "バケットの作成．極の追加" << std::endl;
    TimeWatch tw;

    int num_sources_added = 0;
    for (const auto& f : faces)
      for (const auto& s : f->sources)
        num_sources_added += B_poles.add(s);

    std::cout << "num_sources_added :" << num_sources_added << std::endl;
    std::cout << Magenta << "バケットの作成．極の追加" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

    //@ -------------------------------------------------------------------------- */
    //@                                ツリー構造を生成                            */
    //@ -------------------------------------------------------------------------- */

    std::cout << "ツリー構造を生成" << std::endl;
    B_poles.generateTree();
    std::cout << Magenta << "ツリー構造を生成" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
    is_tree_generated = true;
  }
  // Update cached hashes for the next iteration.
  last_faces_hash = faces_hash;
  last_points_hash = points_hash;
  log_stage("updateTree/generateTree");

  static bool tree_viz_exported = false;
  if (!tree_viz_exported) {
    tree_viz_exported = true;
    exportFMMTreeVisualization(obj);
  }

  //@ -------------------------------------------------------------------------- */
  //@                                  FMM                                       */
  //@ -------------------------------------------------------------------------- */

  const std::size_t tree_sig = computeFMMTreeSignature();
  if (is_fmm_static_initialized && tree_sig != fmm_tree_signature_built)
    is_fmm_static_initialized = false;
  bool reuse_static_fmm = is_fmm_static_initialized && topology_same && reused_sources;
  // Disable with: `export BEM_FMM_REUSE_STATIC=0`
  if (const char* env = std::getenv("BEM_FMM_REUSE_STATIC"); env && std::string(env) == "0")
    reuse_static_fmm = false;

  // Sync midpoint Xtarget for true quadratic elements (networkLine inherits target4FMM)
  if (use_true_quadratic_element) {
    midpoint_lines.clear();
    for (auto* water : WATERS) {
      for (auto* l : water->getBoundaryLines()) {
        if (l->f2Index.empty())
          continue;
        auto [pA, pB] = l->getPoints();
        Tddd Xmid = l->X_mid;
        if (!std::isfinite(Xmid[0]) || !std::isfinite(Xmid[1]) || !std::isfinite(Xmid[2]))
          Xmid = 0.5 * (pA->Xtarget + pB->Xtarget);
        l->Xtarget = Xmid;
        midpoint_lines.push_back(l);
      }
    }
  }

  // Build combined target list for FMM (points + midpoint targets)
  std::vector<target4FMM*> all_targets;
  all_targets.reserve(points.size() + midpoint_lines.size());
  for (auto* p : points)
    all_targets.push_back(static_cast<target4FMM*>(p));
  for (auto* l : midpoint_lines)
    all_targets.push_back(static_cast<target4FMM*>(l));

#if defined(USE_METAL_M2L)
  initializeFMM(B_poles, all_targets, reuse_static_fmm, MetalM2LSettings{use_metal_m2l, metal_m2l_threadgroup, metal_m2l_sort_terms}, nearfield_mode);
#else
  initializeFMM(B_poles, all_targets, reuse_static_fmm, {}, nearfield_mode);
#endif
  if (!reuse_static_fmm) {
    is_fmm_static_initialized = true;
    fmm_tree_signature_built = tree_sig;
  }
  log_stage("initializeFMM");

  if (const char* env = std::getenv("BEM_NEAR_RUN_STATS"); env && std::string(env) != "0") {
    // Sequential stats (debug only, not performance critical)
    std::size_t total_terms = 0;
    std::size_t total_runs = 0;
    std::size_t max_run = 0;
    for (const auto* p : points) {
      total_terms += p->near_indices.size();
      total_runs += p->near_run_len.size();
      for (const auto len : p->near_run_len)
        max_run = std::max<std::size_t>(max_run, static_cast<std::size_t>(len));
    }
    const double avg_run = (total_runs > 0) ? (static_cast<double>(total_terms) / static_cast<double>(total_runs)) : 0.0;
    std::cout << Magenta << "[Near] " << Cyan << "terms=" << Green << total_terms << Cyan << "  runs=" << Green << total_runs << Cyan << "  avg_run=" << Yellow << avg_run << Cyan << "  max_run=" << Yellow << max_run << colorReset << std::endl;
  }

  /*
     基本とする形，左辺
     {{G0,G0,G0,G0},{G1,G1,G1,G1}}.{phin,phin,phin,phin} - {{Gn0,Gn0,Gn0,Gn0},{Gn1,Gn1,Gn1,Gn1}}.{phi,phi,phi,phi}
     node1がNeumanの場合，他はDirichletの場合
     左辺，未知変数側：
     {{G0,G0,G0,G0},{G1,G1,G1,G1}}.{phin,0,phin,phin} - {{Gn0,Gn0,Gn0,Gn0},{Gn1,Gn1,Gn1,Gn1}}.{0,phi,0,0}
     右辺，既知変数側：
     - {{G0,G0,G0,G0},{G1,G1,G1,G1}}.{0, phin, 0, 0} + {{Gn0,Gn0,Gn0,Gn0},{Gn1,Gn1,Gn1,Gn1}}.{phi, 0, phi, phi}
  */

  createCopyMap();

  // Sequential sum (small overhead not worth parallelizing)
  total_unknowns = 0;
  for (const auto& p : points)
    total_unknowns += p->f2Index.size();
  // Include midpoint DOFs for true quadratic elements (CORNER lines have 2 DOFs)
  if (use_true_quadratic_element)
    for (auto* l : midpoint_lines)
      total_unknowns += static_cast<int>(l->f2Index.size());

  updateGeometricProperties();
  // restorePhiPhin();
  copyToFMM(false);
  log_stage("createCopyMap + updateGeometricProperties + copyToFMM");

  // Note: Neumann BC scaling is now handled in copyToFMM()

  std::vector<double> x0 = getVectorFromBoundary(WATERS, total_unknowns, [](const networkPoint* p, networkFace* f) { return p->phinOnFace.at(f); }, [](const networkPoint* p, networkFace* f) { return p->phiOnFace.at(f); });

  time_setup = watch_setup()[0];
  TimeWatch watch_solver;
  ans = solveSystemGMRES(false, x0, true, postprocess);
  time_solve = watch_solver()[0];
  log_stage("solveSystemGMRES");

  //@ -------------------------------------------------------------------------- */
  //@   Unscale computed ∂φ/∂n for coordinate scaling                             */
  //@ -------------------------------------------------------------------------- */
  // The computed ∂φ/∂n' (in scaled coords) needs to be DIVIDED by L to get ∂φ/∂n_original
  if (use_coordinate_scaling_ && coordinate_scale_factor_ > 1e-10) {
    std::cout << Magenta << "[BEM] " << Cyan << "Unscaling computed phin by L=" << coordinate_scale_factor_ << colorReset << std::endl;
    for (auto& p : this->points) {
      for (auto& [f, i] : p->f2Index) {
        if (isDirichletID_BEM(p, f)) {
          // For Dirichlet points, the solution ans[i] is ∂φ/∂n' (in scaled coords)
          // Since ∂φ/∂n' = L × ∂φ/∂n, we need ∂φ/∂n_original = ∂φ/∂n' / L
          ans[i] /= coordinate_scale_factor_;
        }
      }
    }
    // Also unscale midpoint Dirichlet unknowns (phin is the unknown for Dirichlet DOFs)
    if (use_true_quadratic_element) {
      for (auto* water : this->WATERS) {
        for (auto* l : water->getBoundaryLines()) {
          for (const auto& [f, idx] : l->f2Index) {
            if (idx >= 0 && idx < static_cast<int>(ans.size()) && isDirichletID_BEM(l, f))
              ans[idx] /= coordinate_scale_factor_;
          }
        }
      }
    }
    // Midpoint Neumann phin unscaling not needed here — phinOnFace was never pre-scaled
    // (scaling is handled only in copyToFMM on the FMM copy).
  }

  std::cout << colorReset << "update p->phiphin and p->phinOnFace for Dirichlet boundary" << colorReset << std::endl;
  storePhiPhin(WATERS, ans);
  log_stage("storePhiPhin");

  //@ -------------------------------------------------------------------------- */
  //@   Remove coordinate scaling                                                 */
  //@ -------------------------------------------------------------------------- */
  if (use_coordinate_scaling_) {
    std::cout << Magenta << "[BEM] " << Cyan << "Removing coordinate scaling" << colorReset << std::endl;
    for (auto& net : this->WATERS)
      net->removeScaling();
    // Reset scaling flags
    use_coordinate_scaling_ = false;
    coordinate_scale_factor_ = 1.0;
  }

  for (auto water : WATERS)
    isSolutionFinite(*water);
}
