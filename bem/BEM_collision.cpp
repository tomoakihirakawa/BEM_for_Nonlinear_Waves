// Separate translation unit for BEM_collision — collision detection/resolution.
// Compiled independently to reduce main_time_domain.cpp build time.

#define BEM
#include "Network.hpp"
#include "BEM_collision.hpp"

// ============================================================================
// Surface collision detection and resolution — implementation
// ============================================================================

// ---------------------------------------------------------------------------
// Problem B: Adjacent face folding detection
// ---------------------------------------------------------------------------
std::unordered_set<networkFace*> detectFoldedFaces(
    const Network& water,
    double normal_reversal_cos,
    double mean_len) {
   if (mean_len <= 0.0)
      mean_len = Mean(extLength(water.getLines()));
   const double sub_tri_area_threshold = mean_len * mean_len * 1e-8;

   std::unordered_set<networkFace*> folded;
   std::size_t check1_count = 0, check2_degen_count = 0, check2_reversal_count = 0;

   // Check 1: Adjacent face normal reversal (original, vertex-based)
   for (auto* l : water.getBoundaryLines()) {
      auto faces = l->getBoundaryFaces();
      if (faces.size() != 2)
         continue;

      auto* f0 = faces[0];
      auto* f1 = faces[1];
      if (!f0 || !f1)
         continue;

      double dot = Dot(f0->normal, f1->normal);
      if (dot < normal_reversal_cos) {
         folded.insert(f0);
         folded.insert(f1);
         ++check1_count;
      }
   }
   const std::size_t check1_faces = folded.size();

   // Check 2: Sub-triangle normal reversal from midpoint deviation
   std::size_t check2_true_quad_faces = 0;
   for (auto* f : water.getBoundaryFaces()) {
      if (folded.count(f))
         continue;
      if (!f || !f->isTrueQuadraticElement)
         continue;
      ++check2_true_quad_faces;
      auto [p0, l0, p1, l1, p2, l2] = f->PLPLPL;
      Tddd m0 = l0->X_mid, m1 = l1->X_mid, m2 = l2->X_mid;
      Tddd fn = f->normal;
      bool any_degen = false, any_reversal = false;
      auto check = [&](const Tddd& a, const Tddd& b, const Tddd& c) {
         Tddd sub_n = Cross(b - a, c - a);
         double len = Norm(sub_n);
         if (len < sub_tri_area_threshold) { any_degen = true; return false; }
         if (Dot(sub_n / len, fn) < normal_reversal_cos) { any_reversal = true; return true; }
         return false;
      };
      if (check(p0->X, m0, m2) || check(m0, p1->X, m1) ||
          check(m2, m1, p2->X) || check(m0, m1, m2)) {
         folded.insert(f);
         if (any_degen) ++check2_degen_count;
         if (any_reversal) ++check2_reversal_count;
      }
   }

   if (!folded.empty()) {
      std::cout << Yellow << "[detectFoldedFaces] total=" << folded.size()
                << " check1_edges=" << check1_count << " check1_faces=" << check1_faces
                << " check2_true_quad_examined=" << check2_true_quad_faces
                << " check2_degen=" << check2_degen_count
                << " check2_reversal=" << check2_reversal_count
                << " mean_len=" << mean_len
                << " sub_tri_area_threshold=" << sub_tri_area_threshold
                << " normal_reversal_cos=" << normal_reversal_cos
                << colorReset << std::endl;
   }

   return folded;
}

// ---------------------------------------------------------------------------
// Topological distance check
// ---------------------------------------------------------------------------
bool isTopologicallyClose(
    const networkFace* f, const networkFace* g, int max_depth) {
   auto f_points = f->getPoints();
   auto g_points = g->getPoints();

   std::unordered_set<networkPoint*> visited;
   std::vector<networkPoint*> current_ring;
   for (auto* p : f_points) {
      current_ring.push_back(p);
      visited.insert(p);
   }

   for (int depth = 0; depth < max_depth; ++depth) {
      std::vector<networkPoint*> next_ring;
      for (auto* p : current_ring) {
         for (auto* neighbor : p->getNeighbors()) {
            if (visited.count(neighbor))
               continue;
            for (auto* gp : g_points) {
               if (neighbor == gp)
                  return true;
            }
            visited.insert(neighbor);
            next_ring.push_back(neighbor);
         }
      }
      current_ring = std::move(next_ring);
   }
   return false;
}

// ---------------------------------------------------------------------------
// Triangle-to-triangle minimum distance
// ---------------------------------------------------------------------------
double triangleTriangleDistance(const networkFace* f, const networkFace* g) {
   auto fp = f->getPoints();
   auto gp = g->getPoints();
   T3Tddd tri_f = {fp[0]->X, fp[1]->X, fp[2]->X};
   T3Tddd tri_g = {gp[0]->X, gp[1]->X, gp[2]->X};

   double min_dist = 1e30;

   for (auto* p : fp) {
      double d = Norm(Nearest(p->X, tri_g) - p->X);
      if (d < min_dist) min_dist = d;
   }
   for (auto* p : gp) {
      double d = Norm(Nearest(p->X, tri_f) - p->X);
      if (d < min_dist) min_dist = d;
   }
   auto [fl0, fl1, fl2] = f->getLines();
   auto [gl0, gl1, gl2] = g->getLines();
   for (auto* l : {fl0, fl1, fl2}) {
      double d = Norm(Nearest(l->X_mid, tri_g) - l->X_mid);
      if (d < min_dist) min_dist = d;
   }
   for (auto* l : {gl0, gl1, gl2}) {
      double d = Norm(Nearest(l->X_mid, tri_f) - l->X_mid);
      if (d < min_dist) min_dist = d;
   }
   {
      Tddd cf = (fp[0]->X + fp[1]->X + fp[2]->X) / 3.0;
      double d = Norm(Nearest(cf, tri_g) - cf);
      if (d < min_dist) min_dist = d;
   }
   {
      Tddd cg = (gp[0]->X + gp[1]->X + gp[2]->X) / 3.0;
      double d = Norm(Nearest(cg, tri_f) - cg);
      if (d < min_dist) min_dist = d;
   }
   return min_dist;
}

// ---------------------------------------------------------------------------
// Problem A: Non-adjacent surface collision detection
// ---------------------------------------------------------------------------
CollisionDetectionResult detectNonAdjacentCollisions(
    Network& water,
    double proximity_threshold) {
   CollisionDetectionResult result;

   water.makeBucketFaces(proximity_threshold * 2.0);

   auto boundary_faces = water.getBoundaryFaces();

   for (auto* f : boundary_faces) {
      auto fp = f->getPoints();
      Tddd centroid_f = (fp[0]->X + fp[1]->X + fp[2]->X) / 3.0;

      auto nearby = water.BucketFaces.getData(centroid_f, proximity_threshold * 2.0);

      for (auto* g : nearby) {
         if (g == f)
            continue;
         if (f > g)
            continue;

         auto gp = g->getPoints();
         bool shares_vertex = false;
         for (auto* p : fp)
            for (auto* q : gp)
               if (p == q) {
                  shares_vertex = true;
                  break;
               }
         if (shares_vertex)
            continue;

         if (isTopologicallyClose(f, g, 2))
            continue;

         double dist = triangleTriangleDistance(f, g);
         if (dist >= proximity_threshold)
            continue;

         Tddd diff = centroid_f - ((gp[0]->X + gp[1]->X + gp[2]->X) / 3.0);
         double d_norm = Norm(diff);
         if (d_norm < 1e-15)
            continue;

         Tddd dir = diff / d_norm;
         if (Dot(f->normal, -dir) > 0 && Dot(g->normal, dir) > 0) {
            result.faces.insert(f);
            result.faces.insert(g);
            result.pairs.push_back({f, g});
         }
      }
   }
   return result;
}

// ---------------------------------------------------------------------------
// Build connected collision zones from problem faces
// ---------------------------------------------------------------------------
std::vector<CollisionZone> buildCollisionZones(
    const std::unordered_set<networkFace*>& problem_faces,
    const std::vector<std::pair<networkFace*, networkFace*>>& collision_pairs,
    int min_zone_faces) {
   std::unordered_map<networkFace*, std::unordered_set<networkFace*>> adj;
   for (auto* f : problem_faces) {
      for (auto* neighbor : f->getNeighbors())
         if (neighbor && problem_faces.count(neighbor))
            adj[f].insert(neighbor);
   }
   for (auto& [f, g] : collision_pairs) {
      if (problem_faces.count(f) && problem_faces.count(g)) {
         adj[f].insert(g);
         adj[g].insert(f);
      }
   }

   std::vector<CollisionZone> zones;
   std::unordered_set<networkFace*> visited;

   for (auto* seed : problem_faces) {
      if (visited.count(seed))
         continue;

      CollisionZone zone;
      std::queue<networkFace*> queue;
      queue.push(seed);
      visited.insert(seed);

      while (!queue.empty()) {
         auto* f = queue.front();
         queue.pop();
         zone.problem_faces.insert(f);

         for (auto* neighbor : adj[f]) {
            if (!visited.count(neighbor)) {
               visited.insert(neighbor);
               queue.push(neighbor);
            }
         }
      }

      if (static_cast<int>(zone.problem_faces.size()) < min_zone_faces) {
         std::cout << Yellow << "  [buildCollisionZones] skipping zone with "
                   << zone.problem_faces.size() << " faces (< min_zone_faces="
                   << min_zone_faces << ")" << colorReset << std::endl;
         continue;
      }

      std::unordered_set<networkPoint*> problem_points;
      for (auto* f : zone.problem_faces)
         for (auto* p : f->getPoints())
            problem_points.insert(p);

      for (auto* f : zone.problem_faces) {
         for (auto* neighbor : f->getNeighbors()) {
            if (neighbor && !zone.problem_faces.count(neighbor))
               zone.boundary_faces.insert(neighbor);
         }
      }

      for (auto* f : zone.boundary_faces)
         for (auto* p : f->getPoints())
            zone.boundary_points.insert(p);

      for (auto* p : problem_points) {
         if (!zone.boundary_points.count(p))
            zone.interior_points.insert(p);
      }

      zone.all_points = zone.interior_points;
      zone.all_points.insert(zone.boundary_points.begin(), zone.boundary_points.end());

      for (auto* f : zone.problem_faces)
         for (auto* l : f->getLines())
            zone.all_lines.insert(l);

      for (auto* p : zone.all_points) {
         zone.saved_phiphin[p] = p->phiphin;
         zone.saved_phiphin_t[p] = p->phiphin_t;
      }

      zones.push_back(std::move(zone));
   }
   return zones;
}

// ===========================================================================
// Phase 2: TetGen-based collision resolution
// ===========================================================================

#ifdef tetgenH

bool resolveCollisionZone(
    Network& water,
    CollisionZone& zone,
    int time_step,
    int zone_idx) {

   if (zone.all_points.size() < 4) {
      std::cout << "  zone " << zone_idx << ": too few points (" << zone.all_points.size() << "), skipping" << std::endl;
      return false;
   }

   // ==================== Phase 0: Compute junction lines ====================
   std::unordered_set<networkLine*> junction_lines;
   {
      std::unordered_set<networkLine*> problem_line_set;
      for (auto* f : zone.problem_faces)
         for (auto* l : f->getLines())
            problem_line_set.insert(l);
      for (auto* f : zone.boundary_faces)
         for (auto* l : f->getLines())
            if (problem_line_set.count(l))
               junction_lines.insert(l);
   }

   // ==================== Phase 1: TetGen tetrahedralization ====================
   std::vector<networkPoint*> points_vec(zone.all_points.begin(), zone.all_points.end());

   // Pre-check: skip if all points are coplanar
   if (points_vec.size() >= 4) {
      auto& x0 = points_vec[0]->X;
      auto& x1 = points_vec[1]->X;
      auto& x2 = points_vec[2]->X;
      Tddd normal = Cross(x1 - x0, x2 - x0);
      double norm_len = Norm(normal);
      if (norm_len > 1e-15) {
         normal /= norm_len;
         bool coplanar = true;
         for (size_t i = 3; i < points_vec.size(); i++) {
            if (std::abs(Dot(points_vec[i]->X - x0, normal)) > 1e-8) {
               coplanar = false;
               break;
            }
         }
         if (coplanar) {
            std::cout << "  zone " << zone_idx << ": all " << points_vec.size()
                      << " points are coplanar, skipping (TetGen requires 3D)" << std::endl;
            return false;
         }
      }
   }

   std::vector<networkFace*> boundary_faces_vec(zone.boundary_faces.begin(), zone.boundary_faces.end());
   tetgenio in_p = generate_tetgenio_input_mixed(points_vec, boundary_faces_vec);
   tetgenio out_p;
   tetgenio in_d = generate_tetgenio_input_from_points(points_vec);
   tetgenio out_d;

   // Try 1: PLC mode with boundary faces as constraints
   {
      tetgenbehavior b;
      b.parse_commandline(const_cast<char*>("pc"));
      try {
         ::tetrahedralize(&b, &in_p, &out_p);
      } catch (...) {
         out_p.numberoftetrahedra = 0;
      }
   }

   std::string tetgen_mode = "PLC";
   tetgenio* p_out = &out_p;

   if (out_p.numberoftetrahedra == 0) {
      tetgen_mode = "Delaunay";
      p_out = &out_d;
      tetgenbehavior b2;
      b2.parse_commandline(const_cast<char*>(""));
      try {
         ::tetrahedralize(&b2, &in_d, &out_d);
      } catch (...) {
         std::cerr << Red << "[collision] zone " << zone_idx
                   << ": TetGen failed (both PLC and Delaunay), skipping" << colorReset << std::endl;
         return false;
      }
   }
   auto& out = *p_out;

   if (out.numberoftetrahedra == 0) {
      std::cout << "  zone " << zone_idx << ": 0 tetrahedra, skipping" << std::endl;
      return false;
   }

   std::cout << "  zone " << zone_idx << ": TetGen produced "
             << out.numberoftetrahedra << " tetras from "
             << out.numberofpoints << " points (" << tetgen_mode << ")" << std::endl;

   // Build index-to-point map
   std::vector<networkPoint*> idx_to_point(out.numberofpoints, nullptr);
   for (auto* p : points_vec)
      idx_to_point[p->index] = p;

   // Handle Steiner points
   std::vector<networkPoint*> steiner_points;
   if (out.numberofpoints > static_cast<int>(points_vec.size())) {
      water.makeBucketPoints(water.getScale() / 50.);
      for (int i = static_cast<int>(points_vec.size()); i < out.numberofpoints; i++) {
         Tddd X = {out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2]};
         auto* p = new networkPoint(&water, X);
         water.BucketPoints.add(X, p);
         idx_to_point[i] = p;
         steiner_points.push_back(p);
      }
      std::cout << "  zone " << zone_idx << ": " << steiner_points.size() << " Steiner points" << std::endl;
   }

   // Create tetrahedra
   std::vector<networkTetra*> zone_tetras;
   zone_tetras.reserve(out.numberoftetrahedra);
   for (int i = 0; i < out.numberoftetrahedra; i++) {
      auto* p0 = idx_to_point[out.tetrahedronlist[i * 4 + 0]];
      auto* p1 = idx_to_point[out.tetrahedronlist[i * 4 + 1]];
      auto* p2 = idx_to_point[out.tetrahedronlist[i * 4 + 2]];
      auto* p3 = idx_to_point[out.tetrahedronlist[i * 4 + 3]];
      if (!p0 || !p1 || !p2 || !p3) continue;
      try {
         auto [success, tet] = genTetra(&water, p0, p1, p2, p3);
         if (tet) zone_tetras.push_back(tet);
      } catch (...) {
      }
   }

   // --- Abort cleanup lambda ---
   auto abort_cleanup = [&]() {
      std::unordered_set<networkFace*> tetra_faces;
      for (auto* t : zone_tetras)
         for (auto* f : t->Faces)
            if (f) tetra_faces.insert(f);
      for (auto* t : zone_tetras) delete t;
      zone_tetras.clear();
      for (auto* f : tetra_faces) {
         if (zone.problem_faces.count(f) || zone.boundary_faces.count(f)) continue;
         if (f->Tetras[0] == nullptr && f->Tetras[1] == nullptr)
            delete f;
      }
      for (auto* p : steiner_points) {
         auto lines = p->getLines();
         for (auto* l : lines)
            if (l->getFaces().empty()) delete l;
      }
      for (auto* p : steiner_points)
         if (p->getLines().empty()) delete p;
   };

   // Classify tetrahedra as water/air
   std::vector<networkFace*> ref_faces(zone.boundary_faces.begin(), zone.boundary_faces.end());
   auto classifyPoint = [&](const Tddd& X) -> bool {
      double min_dist = 1e30;
      networkFace* nearest_face = nullptr;
      Tddd nearest_pt;
      for (auto* f : ref_faces) {
         Tddd pt = Nearest(X, f);
         double d = Norm(X - pt);
         if (d < min_dist) {
            min_dist = d;
            nearest_pt = pt;
            nearest_face = f;
         }
      }
      if (!nearest_face) return false;
      return Dot(X - nearest_pt, nearest_face->normal) < 0;
   };

   std::unordered_set<networkTetra*> water_set;
   for (auto* tet : zone_tetras) {
      auto [p0, p1, p2, p3] = tet->Points;
      Tddd centroid = (p0->X + p1->X + p2->X + p3->X) / 4.0;
      if (classifyPoint(centroid))
         water_set.insert(tet);
   }

   std::cout << "  zone " << zone_idx << ": " << water_set.size() << " water / "
             << (zone_tetras.size() - water_set.size()) << " air" << std::endl;

   if (water_set.empty()) {
      std::cout << Red << "  zone " << zone_idx << ": no water tetras, skipping" << colorReset << std::endl;
      abort_cleanup();
      return false;
   }

   // --- Debug output: tetras before refinement ---
   {
      std::unordered_map<networkPoint*, double> is_water_data;
      for (auto* tet : zone_tetras) {
         double val = water_set.count(tet) ? 1.0 : 0.0;
         auto [a, b, c, d] = tet->Points;
         for (auto* p : {a, b, c, d})
            is_water_data[p] = val;
      }
      std::ofstream ofs("collision_zone" + std::to_string(zone_idx)
                        + "_step" + std::to_string(time_step) + "_tetras_pre.vtu");
      vtkUnstructuredGridWrite(ofs, zone_tetras,
                               {{"is_water", DataVariant3(is_water_data)}});
      std::cout << "  zone " << zone_idx << ": wrote _tetras_pre.vtu" << std::endl;
   }
   // --- Debug output: boundary faces ---
   {
      std::ofstream ofs("collision_zone" + std::to_string(zone_idx)
                        + "_step" + std::to_string(time_step) + "_boundary.vtu");
      vtkUnstructuredGridWrite(ofs, zone.boundary_faces);
      std::cout << "  zone " << zone_idx << ": wrote _boundary.vtu" << std::endl;
   }

   // ==================== Phase 2: Size-based refinement ====================
   double mean_len = Mean(extLength(water.getLines()));
   double min_alt_thresh = mean_len * 0.1;

   auto tetra_min_altitude = [](networkTetra* tet) -> double {
      auto [p0, p1, p2, p3] = tet->Points;
      double vol = TetrahedronVolume(p0->X, p1->X, p2->X, p3->X);
      double max_area = std::max({TriangleArea(p0->X, p1->X, p2->X),
                                  TriangleArea(p0->X, p1->X, p3->X),
                                  TriangleArea(p0->X, p2->X, p3->X),
                                  TriangleArea(p1->X, p2->X, p3->X)});
      return (max_area > 1e-30) ? 3.0 * vol / max_area : 0.0;
   };

   int n_merge = 0, n_trim = 0;
   for (auto* tet : zone_tetras) {
      if (tetra_min_altitude(tet) >= min_alt_thresh) continue;
      if (water_set.count(tet)) {
         water_set.erase(tet);
         n_trim++;
      } else {
         water_set.insert(tet);
         n_merge++;
      }
   }

   if (n_merge || n_trim) {
      std::cout << "  zone " << zone_idx << ": refined +" << n_merge
                << " merged, -" << n_trim << " trimmed -> "
                << water_set.size() << " water" << std::endl;
      {
         std::unordered_map<networkPoint*, double> is_water_data;
         for (auto* tet : zone_tetras) {
            double val = water_set.count(tet) ? 1.0 : 0.0;
            auto [a, b, c, d] = tet->Points;
            for (auto* p : {a, b, c, d})
               is_water_data[p] = val;
         }
         std::ofstream ofs("collision_zone" + std::to_string(zone_idx)
                           + "_step" + std::to_string(time_step) + "_tetras_post.vtu");
         vtkUnstructuredGridWrite(ofs, zone_tetras,
                                  {{"is_water", DataVariant3(is_water_data)}});
         std::cout << "  zone " << zone_idx << ": wrote _tetras_post.vtu" << std::endl;
      }
   }

   if (water_set.empty()) {
      std::cout << Red << "  zone " << zone_idx << ": empty after refinement" << colorReset << std::endl;
      abort_cleanup();
      return false;
   }

   // ==================== Phase 3: Remove isolated water fragments ====================
   {
      std::vector<std::vector<networkTetra*>> comps;
      std::unordered_set<networkTetra*> visited;
      for (auto* seed : water_set) {
         if (visited.count(seed)) continue;
         std::vector<networkTetra*> comp;
         std::queue<networkTetra*> q;
         q.push(seed);
         visited.insert(seed);
         while (!q.empty()) {
            auto* t = q.front();
            q.pop();
            comp.push_back(t);
            for (auto* f : t->Faces) {
               if (!f) continue;
               auto [t0, t1] = f->Tetras;
               auto* nb = (t0 == t) ? t1 : t0;
               if (nb && water_set.count(nb) && !visited.count(nb)) {
                  visited.insert(nb);
                  q.push(nb);
               }
            }
         }
         comps.push_back(std::move(comp));
      }
      if (comps.size() > 1) {
         size_t best = 0;
         for (size_t i = 1; i < comps.size(); i++)
            if (comps[i].size() > comps[best].size()) best = i;
         int removed = 0;
         for (size_t i = 0; i < comps.size(); i++) {
            if (i == best) continue;
            for (auto* t : comps[i]) {
               water_set.erase(t);
               removed++;
            }
         }
         std::cout << "  zone " << zone_idx << ": removed " << removed
                   << " isolated water tetras (" << (comps.size() - 1) << " fragments)" << std::endl;
      }
   }

   if (water_set.empty()) {
      std::cout << Red << "  zone " << zone_idx << ": empty after fragment removal" << colorReset << std::endl;
      abort_cleanup();
      return false;
   }

   // ==================== Phase 4: Surface extraction ====================
   std::unordered_set<networkFace*> zone_tetra_faces;
   for (auto* t : zone_tetras)
      for (auto* f : t->Faces)
         if (f) zone_tetra_faces.insert(f);

   std::unordered_set<networkFace*> new_surface_faces;
   for (auto* f : zone_tetra_faces) {
      auto [t0, t1] = f->Tetras;
      int wc = (t0 && water_set.count(t0) ? 1 : 0)
             + (t1 && water_set.count(t1) ? 1 : 0);
      if (wc == 1)
         new_surface_faces.insert(f);
   }

   std::cout << "  zone " << zone_idx << ": " << new_surface_faces.size() << " new surface faces" << std::endl;

   {
      std::ofstream ofs("collision_zone" + std::to_string(zone_idx)
                        + "_step" + std::to_string(time_step) + "_surface.vtu");
      vtkUnstructuredGridWrite(ofs, new_surface_faces);
      std::cout << "  zone " << zone_idx << ": wrote _surface.vtu" << std::endl;
   }

   if (new_surface_faces.empty()) {
      std::cout << Red << "  zone " << zone_idx << ": no surface faces, aborting" << colorReset << std::endl;
      abort_cleanup();
      return false;
   }

   // ==================== Phase 5: Verify boundary connectivity ====================
   std::unordered_set<networkLine*> new_surface_line_set;
   for (auto* f : new_surface_faces)
      for (auto* l : f->getLines())
         new_surface_line_set.insert(l);

   int n_missing_junctions = 0;
   for (auto* l : junction_lines) {
      if (!new_surface_line_set.count(l))
         n_missing_junctions++;
   }
   if (n_missing_junctions > 0) {
      std::cout << Yellow << "  zone " << zone_idx
                << ": " << n_missing_junctions << "/" << junction_lines.size()
                << " junction lines not in new surface, aborting"
                << colorReset << std::endl;
      abort_cleanup();
      return false;
   }

   // ==================== Phase 6: Cleanup (transactional) ====================
   std::unordered_set<networkFace*> faces_to_keep;
   faces_to_keep.insert(new_surface_faces.begin(), new_surface_faces.end());
   faces_to_keep.insert(zone.boundary_faces.begin(), zone.boundary_faces.end());

   std::unordered_set<networkFace*> faces_to_delete;
   for (auto* f : zone.problem_faces)
      if (!faces_to_keep.count(f))
         faces_to_delete.insert(f);
   for (auto* f : zone_tetra_faces)
      if (!faces_to_keep.count(f))
         faces_to_delete.insert(f);
   for (auto* f : zone.boundary_faces)
      faces_to_delete.erase(f);

   {
      bool pre_check_ok = true;
      for (auto* l : junction_lines) {
         int kept_face_count = 0;
         for (auto* f : l->getFaces()) {
            if (f && !faces_to_delete.count(f))
               kept_face_count++;
         }
         if (kept_face_count != 2) {
            std::cout << Yellow << "  zone " << zone_idx
                      << ": junction line will have " << kept_face_count
                      << " faces after deletion (expected 2), aborting"
                      << colorReset << std::endl;
            pre_check_ok = false;
            break;
         }
      }
      if (pre_check_ok) {
         for (auto* l : new_surface_line_set) {
            if (junction_lines.count(l)) continue;
            int kept_face_count = 0;
            for (auto* f : l->getFaces()) {
               if (f && !faces_to_delete.count(f))
                  kept_face_count++;
            }
            if (kept_face_count != 2) {
               std::cout << Yellow << "  zone " << zone_idx
                         << ": surface line will have " << kept_face_count
                         << " faces after deletion (expected 2), aborting"
                         << colorReset << std::endl;
               pre_check_ok = false;
               break;
            }
         }
      }
      if (!pre_check_ok) {
         abort_cleanup();
         return false;
      }
   }

   // All pre-checks passed. Execute destructive operations.
   for (auto* t : zone_tetras) delete t;
   zone_tetras.clear();

   for (auto* f : faces_to_delete) delete f;

   std::unordered_set<networkLine*> lines_to_check;
   for (auto* l : zone.all_lines)
      lines_to_check.insert(l);
   for (auto* p : zone.interior_points)
      for (auto* l : p->getLines())
         lines_to_check.insert(l);
   for (auto* p : steiner_points)
      for (auto* l : p->getLines())
         lines_to_check.insert(l);

   for (auto* l : lines_to_check)
      if (l->getFaces().empty()) delete l;

   for (auto* p : zone.interior_points)
      if (p->getLines().empty()) delete p;

   std::vector<networkPoint*> surviving_steiner;
   for (auto* p : steiner_points) {
      if (p->getLines().empty())
         delete p;
      else
         surviving_steiner.push_back(p);
   }

   // ==================== Phase 7: Steiner point interpolation ====================
   for (auto* p : surviving_steiner) {
      auto neighbors = p->getNeighbors();
      if (neighbors.empty()) continue;

      double total_w = 0;
      Tdd weighted_phiphin = {0, 0};
      Tdd weighted_phiphin_t = {0, 0};
      for (auto* q : neighbors) {
         double w = 1.0 / std::max(Norm(p->X - q->X), 1e-15);
         auto it = zone.saved_phiphin.find(q);
         if (it != zone.saved_phiphin.end()) {
            weighted_phiphin += w * it->second;
            auto it_t = zone.saved_phiphin_t.find(q);
            weighted_phiphin_t += w * (it_t != zone.saved_phiphin_t.end() ? it_t->second : q->phiphin_t);
         } else {
            weighted_phiphin += w * q->phiphin;
            weighted_phiphin_t += w * q->phiphin_t;
         }
         total_w += w;
      }
      if (total_w > 0) {
         p->phiphin = weighted_phiphin / total_w;
         p->phiphin_t = weighted_phiphin_t / total_w;
      }
   }

   std::cout << Green << "  zone " << zone_idx << ": resolved ("
             << new_surface_faces.size() << " surface faces, "
             << surviving_steiner.size() << " Steiner pts)" << colorReset << std::endl;

   return true;
}

#endif // tetgenH

// ---------------------------------------------------------------------------
// Top-level collision detection and resolution
// ---------------------------------------------------------------------------
bool detectAndResolveCollisions(
    Network& water,
    int time_step,
    const CollisionSettings& settings,
    std::unordered_set<networkLine*>& protected_lines) {
   if (!settings.enabled)
      return true;

   double mean_len = Mean(extLength(water.getLines()));
   double threshold = mean_len * settings.proximity_factor;

   // Phase 1a: Resolve folded faces by edge collapse
   if (settings.detect_folding) {
      auto folded = detectFoldedFaces(water, settings.normal_reversal_cos, mean_len);
      if (!folded.empty()) {
         std::cout << Red << "[collision] time_step " << time_step
                   << ": detected " << folded.size()
                   << " folded faces (adjacent face folding)" << colorReset << std::endl;

         int n_midpoint_reset = 0;
         std::unordered_set<networkLine*> reset_lines;
         for (auto* f : folded) {
            if (!f) continue;
            auto [l0, l1, l2] = f->getLines();
            for (auto* l : {l0, l1, l2}) {
               if (!l || reset_lines.count(l)) continue;
               auto [p0, p1] = l->getPoints();
               if (!p0 || !p1) continue;
               Tddd linear_mid = (p0->X + p1->X) * 0.5;
               double deviation = Norm(l->X_mid - linear_mid);
               if (deviation > mean_len * 1e-2) {
                  l->setXSingle(linear_mid);
                  reset_lines.insert(l);
                  ++n_midpoint_reset;
               }
            }
         }
         if (n_midpoint_reset > 0) {
            water.setGeometricPropertiesForce();
            std::cout << Green << "[collision] time_step " << time_step
                      << ": resolved " << n_midpoint_reset
                      << " folded edge midpoint(s) by linear reset (no face removal)"
                      << colorReset << std::endl;
         }

         constexpr double max_fold_ratio = 0.15;
         auto remaining_folds = (n_midpoint_reset > 0) ? detectFoldedFaces(water, settings.normal_reversal_cos, mean_len) : folded;
         size_t n_boundary_faces = water.getBoundaryFaces().size();
         double fold_ratio = (n_boundary_faces > 0) ? static_cast<double>(remaining_folds.size()) / n_boundary_faces : 0.0;
         if (!remaining_folds.empty()) {
            std::cout << Yellow << "[collision] time_step " << time_step
                      << ": " << remaining_folds.size() << " / " << n_boundary_faces
                      << " folds remain" << (n_midpoint_reset > 0 ? " after midpoint reset" : " (vertex-driven)")
                      << " (ratio=" << fold_ratio << ")"
                      << colorReset << std::endl;
         }
         if (fold_ratio > max_fold_ratio) {
            throw step_failure("fold ratio " + std::to_string(fold_ratio) + " > " + std::to_string(max_fold_ratio)
                               + (n_midpoint_reset > 0 ? " after midpoint reset" : " (vertex-driven, midpoint reset ineffective)")
                               + " at time_step " + std::to_string(time_step));
         }
      }
   }

   // Phase 1b: Detect non-adjacent collisions
   std::unordered_set<networkFace*> problem_faces;
   std::vector<std::pair<networkFace*, networkFace*>> collision_pairs;

   if (settings.detect_non_adjacent) {
      auto result = detectNonAdjacentCollisions(water, threshold);
      if (!result.faces.empty()) {
         std::cout << Red << "[collision] time_step " << time_step
                   << ": detected " << result.faces.size()
                   << " colliding faces in " << result.pairs.size()
                   << " pairs (non-adjacent collision, threshold=" << threshold << ")"
                   << colorReset << std::endl;
         problem_faces.insert(result.faces.begin(), result.faces.end());
         collision_pairs = std::move(result.pairs);
      }
   }

   if (problem_faces.empty()) {
      return true;
   }

   // Phase 2: Build collision zones
   auto zones = buildCollisionZones(problem_faces, collision_pairs, settings.min_zone_faces);
   std::cout << Green << "[collision] time_step " << time_step
             << ": " << zones.size() << " collision zone(s) identified" << colorReset << std::endl;

   for (size_t i = 0; i < zones.size(); ++i) {
      const auto& zone = zones[i];
      std::cout << "  zone " << i << ": "
                << zone.problem_faces.size() << " problem faces, "
                << zone.boundary_faces.size() << " boundary faces, "
                << zone.interior_points.size() << " interior points, "
                << zone.boundary_points.size() << " boundary points"
                << std::endl;
   }

   auto protect_faces = [&](const auto& faces) {
      for (auto* f : faces)
         for (auto* l : f->getLines())
            protected_lines.insert(l);
   };

   if (zones.empty()) {
      std::cout << Yellow << "[collision] time_step " << time_step
                << ": " << problem_faces.size() << " colliding faces but no zones large enough to resolve"
                << colorReset << std::endl;
      protect_faces(problem_faces);
      if (!settings.resolve_with_tetgen) {
         std::cout << Yellow << "[collision] time_step " << time_step
                   << ": monitor-only mode, protected " << protected_lines.size()
                   << " line(s) near detected collisions" << colorReset << std::endl;
         return true;
      }
      return false;
   }

   if (!settings.resolve_with_tetgen) {
      std::cout << Yellow << "[collision] time_step " << time_step
                << ": monitor-only mode, " << zones.size()
                << " collision zone(s) detected; TetGen repair skipped"
                << colorReset << std::endl;
      return true;
   }

   // Phase 3: Resolution (TetGen-based)
#ifdef tetgenH
   int resolved_count = 0;
   std::vector<bool> zone_resolved(zones.size(), false);
   for (size_t i = 0; i < zones.size(); ++i) {
      try {
         if (resolveCollisionZone(water, zones[i], time_step, static_cast<int>(i))) {
            resolved_count++;
            zone_resolved[i] = true;
         }
      } catch (std::exception& e) {
         std::cerr << Red << "[collision] zone " << i
                   << ": resolution failed with exception: " << e.what()
                   << colorReset << std::endl;
      }
   }

   if (resolved_count > 0) {
      water.setGeometricPropertiesForce();
      water.makeBuckets();
      std::cout << Green << "[collision] time_step " << time_step
                << ": resolved " << resolved_count << " / " << zones.size()
                << " collision zone(s)" << colorReset << std::endl;
   }

   bool all_resolved = (resolved_count == static_cast<int>(zones.size()));
   if (!all_resolved) {
      for (size_t i = 0; i < zones.size(); ++i) {
         if (!zone_resolved[i])
            protect_faces(zones[i].problem_faces);
      }
      std::cout << Yellow << "[collision] time_step " << time_step
                << ": " << (zones.size() - resolved_count) << " of " << zones.size()
                << " zone(s) unresolved, protected " << protected_lines.size() << " lines"
                << colorReset << std::endl;
   }
   return all_resolved;
#else
   std::cout << "[collision] TetGen not available, resolution skipped" << std::endl;
   return false;
#endif
}
