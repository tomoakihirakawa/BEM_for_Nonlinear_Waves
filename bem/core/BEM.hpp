#ifndef BEM_H
#define BEM_H

#include "BEM_BoundaryValues.hpp"
#include "BEM_calculateVelocities.hpp"
#include "BEM_setBoundaryTypes.hpp"
#include "BEM_solveBVP.hpp"
#include "Network.hpp"
#include "dunavant_rules.hpp"

// b! ------------------------------------------------------ */
// b!            格子のdivide, merge．それに伴うΦ，Φnの付与          */
// b! ------------------------------------------------------ */

Tdd estimate_phiphin(const networkLine* const l) {

  auto [a, b] = l->getPoints();
  return (a->phiphin + b->phiphin) / 2.;
};

/* ------------------------------------------------------ */
void remesh(Network& water, const Tdd& limit_angle_D, const Tdd& limit_angle_N, bool force = false, int max_count = 100) {
  std::cout << "remeshing" << std::endl;
  water.setGeometricPropertiesForce();
  double mean_length = Mean(extLength(water.getLines()));
  bool isfound = false, ismerged = false;
  int count = 0;

  networkLine* l;
  Tddd X, V;
  networkPoint* q;
  double meanArea;
  V_netFp Fs;
  Tdd phiphin;
  V_netLp lines;
  double local_mean_length;
  do {
    // なくなるまでやるか？
    isfound = false;
    ismerged = false;
    /* ------------------------------------------------------ */
    for (const auto& p : RandomSample(ToVector(water.getPoints()))) {
      /* ------------------------------------------------------ */
      meanArea = Mean(p->getFaceAreas());
      //! ------------------------------------------------------ */
      //!             辺の長さが長すぎるまたは短すぎる場合             */
      //! ------------------------------------------------------ */
      local_mean_length = Mean(extLength(extractLines(Flatten(BFS(p, 2)))));
      lines = p->getLines();
      sortByLength(lines);
      for (const auto& l : Reverse(lines)) {
        auto [p0, p1] = l->getPoints();
        Fs = l->getBoundaryFaces();
        if (l->length() > local_mean_length * 1.5 /*長すぎる*/) {
          phiphin = estimate_phiphin(l);
          /* ------------------------------------------------------ */
          q = l->Split();
          q->phiphin = phiphin;
          isfound = true;
          break;
        }
        //@ ------------------------------------------------------ */
        if (l->length() < local_mean_length / 20.) {
          /*
          b!マージの原則：マージによってノイマン面は変形してはいけない．
          */
          //@ case5 ２点の平均，移動位置は，ノイマンを崩さない方向：ノイマン面の法線方向成分には移動しない．
          auto [a, b] = l->getPoints();
          if (!(a->BCInterface && b->BCInterface)) {
            if (a->BCInterface)
              b = a;
            else if (b->BCInterface)
              a = b;
          }
          phiphin = (a->phiphin + b->phiphin) / 2.;
          V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
          V -= a->getNormalTuple() * Dot(V, a->getNormalTuple());
          V -= b->getNormalTuple() * Dot(V, b->getNormalTuple());
          X = V + a->getXtuple();
          q = l->Collapse();
          q->phiphin = phiphin;
          q->setX(X);
          /* ------------------------------------------------------ */
          ismerged = true;
          break;
        }
        if (ismerged || isfound)
          break;
      }

      if (ismerged || isfound)
        break;
    }
  } while ((ismerged || isfound) && count++ < max_count);

  flipIf(water, limit_angle_D, limit_angle_N, force);
};

/* -------------------------------------------------------------------------- */
/*                                     出力                                    */
/* -------------------------------------------------------------------------- */

template <typename V>
std::unordered_map<networkPoint*, V> init_map(const Network* network, const V& value) {
  std::unordered_map<networkPoint*, V> m;
  for (const auto p : network->getPoints())
    m[p] = value;
  return m;
}

// DOF-keyed init: includes both points and (optionally) boundary line midpoints
template <typename V>
std::unordered_map<BEM_DOF_Base*, V> init_dof_map(const Network* network, const V& value) {
  std::unordered_map<BEM_DOF_Base*, V> m;
  for (auto* p : network->getPoints())
    m[p] = value;
  if (use_true_quadratic_element)
    for (auto* l : network->getBoundaryLines())
      m[l] = value;
  return m;
}

// boundary condition value for ParaView output (point/line 共通)
// 0=BCInterface, 1=multiple Neumann, 2=Neumann, 3=multiple Dirichlet, 4=Dirichlet
inline double boundaryConditionValue(const auto* entity) {
  if (entity->BCInterface) return 0.;
  if (entity->isMultipleNode && entity->Neumann) return 1.;
  if (entity->Neumann) return 2.;
  if (entity->isMultipleNode && entity->Dirichlet) return 3.;
  if (entity->Dirichlet) return 4.;
  return -1.; // undefined
}

VV_VarForOutput dataForOutput(const Network* water, const double dt) {
  try {

    const Tddd tdd0 = {1E+30, 1E+30, 1E+30};
    const double d0 = 1E+30;
    constexpr double output_sharpq_feature_angle = 60.0 * M_PI / 180.0;
    constexpr double output_sharpq_feature_angle_deg = 60.0;

    const uomap_DOF_Tddd dof_tdd0 = init_dof_map(water, tdd0);
    const uomap_DOF_d dof_d0 = init_dof_map(water, d0);

    uomap_DOF_Tddd P_accel_body = dof_tdd0;
    uomap_DOF_Tddd P_velocity_body = dof_tdd0;

    uomap_DOF_Tddd P_accelNeumann = dof_tdd0;
    uomap_DOF_Tddd P_phin_Dirichlet = dof_tdd0;
    uomap_DOF_Tddd P_U_BEM = dof_tdd0;
    uomap_DOF_Tddd P_U_shift_BEM = dof_tdd0;
    uomap_DOF_Tddd P_position = dof_tdd0;
    uomap_DOF_Tddd P_normal_BEM = dof_tdd0;
    uomap_DOF_Tddd P_gradPhi = dof_tdd0;
    uomap_DOF_Tddd P_u_total = dof_tdd0;
    uomap_DOF_Tddd P_u_reloc = dof_tdd0;
    uomap_DOF_Tddd P_u_potential_BEM = dof_tdd0;
    uomap_DOF_Tddd P_u_omega_VPM = dof_tdd0;
    uomap_DOF_Tddd P_velocity_convergence = dof_tdd0;
    uomap_DOF_Tddd P_vecToSurface = dof_tdd0;
    uomap_DOF_Tddd P_clungSurface = dof_tdd0;
    uomap_DOF_Tddd P_uNeumann = dof_tdd0;
    uomap_DOF_Tddd P_V2ContactFaces0 = dof_tdd0, P_V2ContactFaces1 = dof_tdd0, P_V2ContactFaces2 = dof_tdd0, P_V2ContactFaces3 = dof_tdd0, P_V2ContactFaces4 = dof_tdd0, P_V2ContactFaces5 = dof_tdd0, P_u_absorbed = dof_tdd0;

    uomap_DOF_d P_isMultipleNode = dof_d0;
    uomap_DOF_d P_phi = dof_d0;
    uomap_DOF_d P_phi_absorbed = dof_d0;
    uomap_DOF_d P_absorb_gamma = dof_d0;
    uomap_DOF_d P_solidAngle = dof_d0;
    uomap_DOF_d P_b_diff_RHS_FMM = dof_d0;
    uomap_DOF_d P_b_RHS_FMM = dof_d0;
    uomap_DOF_d P_b_RHS_Direct = dof_d0;
    uomap_DOF_d P_phin = dof_d0;
    uomap_DOF_d P_phi_t = dof_d0;
    uomap_DOF_d P_phin_t = dof_d0;
    uomap_DOF_d P_phin_t_from_Hessian = dof_d0;
    uomap_DOF_d P_pressure = dof_d0;
    uomap_DOF_d P_DphiDt = dof_d0;
    uomap_DOF_d P_ContactFaces = dof_d0;
    uomap_DOF_d P_ContactFaces_maxdist = dof_d0;
    uomap_DOF_d P_ContactRange = dof_d0;
    uomap_DOF_d P_facesNeuamnn = dof_d0;
    uomap_DOF_d P_BC = dof_d0;
    uomap_DOF_d P_diag = dof_d0;
    uomap_DOF_d P_isAbsorbed = dof_d0;
    uomap_DOF_d P_isAbsorbed_SDF = dof_d0;
    uomap_DOF_d P_minDepthFromBCInterface = dof_d0;
    uomap_DOF_d P_minDepthFromMultipleNode = dof_d0;
    uomap_DOF_d P_almost_solid_angle = dof_d0;
    uomap_DOF_d P_penetration_dist = dof_d0;
    uomap_DOF_d P_direction_info_count = dof_d0;
    uomap_DOF_d P_contact_faces_count = dof_d0;
    uomap_DOF_d P_body_vertices_count = dof_d0;
    uomap_DOF_d P_isInContact_pass_count = dof_d0;
    uomap_DOF_d P_SharpQ = dof_d0;
    uomap_DOF_d P_SharpQ_direct = dof_d0;
    uomap_DOF_d P_SharpQ_incident_edge_angle_max_deg = dof_d0;
    uomap_DOF_d P_pf_neumann_count = dof_d0;
    uomap_DOF_d P_pf_dirichlet_count = dof_d0;
    uomap_DOF_d P_pf_total_count = dof_d0;
    uomap_DOF_d P_pf_pressure_min = dof_d0;
    uomap_DOF_d P_pf_detach_count_max = dof_d0;
    uomap_DOF_d P_pf_detach_progress_max = dof_d0;
    uomap_DOF_d P_pf_pressure_below_count = dof_d0;
    uomap_DOF_d P_pf_pressure_detachment_eligible_count = dof_d0;
    uomap_DOF_d P_pf_detached_by_pressure_count = dof_d0;
    uomap_DOF_d P_lf_neumann_count = dof_d0;
    uomap_DOF_d P_lf_dirichlet_count = dof_d0;
    uomap_DOF_d P_lf_total_count = dof_d0;
    uomap_DOF_d P_lf_pressure_min = dof_d0;
    uomap_DOF_d P_lf_detach_count_max = dof_d0;
    uomap_DOF_d P_lf_detach_progress_max = dof_d0;
    uomap_DOF_d P_lf_pressure_below_count = dof_d0;
    uomap_DOF_d P_lf_pressure_detachment_eligible_count = dof_d0;
    uomap_DOF_d P_lf_detached_by_pressure_count = dof_d0;

    uomap_F_d F_lagrangian_surface_nearest_mean;
    uomap_F_d F_lagrangian_surface_nearest_max;
    uomap_F_d F_lagrangian_surface_global_fallback;
    uomap_F_d F_SharpQ_face_edge_count;
    uomap_F_d F_SharpQ_face_point_count;
    uomap_F_d F_SharpQ_face_has_any;
    uomap_F_d F_SharpQ_face_direct_edge_count;
    uomap_F_d F_SharpQ_face_max_edge_angle_deg;
    uomap_F_d F_SharpQ_face_direct_has_any;

    auto lag_position = [](const networkPoint* p) -> Tddd {
      // Compare the corrected surface against the same-step pure Lagrangian
      // surface before relocation. After relocation, the current point position
      // already includes vecToSurface, so subtract it back out here.
      return getPosition(p) - p->vecToSurface;
    };

    auto triangle_from_face = [&](const networkFace* f, const auto& position_of) -> T3Tddd {
      auto [a, b, c] = f->getPoints();
      return {position_of(a), position_of(b), position_of(c)};
    };

    auto collect_local_face_candidates = [](const networkFace* f) {
      std::unordered_set<networkFace*> ret;
      if (!f)
        return ret;
      ret.emplace(const_cast<networkFace*>(f));
      for (auto* l : f->getLines())
        for (auto* g : l->getBoundaryFaces())
          if (g)
            ret.emplace(g);
      for (auto* p : f->getPoints())
        for (auto* g : p->getBoundaryFaces())
          if (g)
            ret.emplace(g);
      return ret;
    };

    auto line_max_boundary_normal_angle_deg = [](const networkLine* l) {
      if (!l)
        return 0.0;
      const auto faces = l->getBoundaryFaces();
      double max_angle = 0.0;
      for (std::size_t i = 0; i < faces.size(); ++i) {
        if (!faces[i])
          continue;
        for (std::size_t j = i + 1; j < faces.size(); ++j) {
          if (!faces[j])
            continue;
          const double c = std::clamp(Dot(faces[i]->normal, faces[j]->normal), -1.0, 1.0);
          max_angle = std::max(max_angle, std::acos(c) * 180.0 / M_PI);
        }
      }
      return max_angle;
    };
    auto line_direct_sharp = [&](const networkLine* l) {
      return line_max_boundary_normal_angle_deg(l) > output_sharpq_feature_angle_deg;
    };
    auto point_direct_sharp = [&](const networkPoint* p) {
      if (!p)
        return false;
      for (auto* l : p->getBoundaryLines())
        if (line_direct_sharp(l))
          return true;
      return false;
    };
    auto point_max_incident_edge_angle_deg = [&](const networkPoint* p) {
      double max_angle = 0.0;
      if (!p)
        return max_angle;
      for (auto* l : p->getBoundaryLines())
        max_angle = std::max(max_angle, line_max_boundary_normal_angle_deg(l));
      return max_angle;
    };

    auto nearest_distance_to_lagrangian_faces = [&](const Tddd& X,
                                                    const std::unordered_set<networkFace*>& candidates) {
      double best = std::numeric_limits<double>::infinity();
      for (auto* g : candidates) {
        if (!g)
          continue;
        const auto tri = triangle_from_face(g, lag_position);
        best = std::min(best, Norm(Nearest(X, tri) - X));
      }
      return best;
    };

    {
      const auto dunavant3 = getDunavantP2MRule(3);
      const auto& faces = water->getBoundaryFaces();
      const std::unordered_set<networkFace*> all_faces(faces.begin(), faces.end());
      std::vector<double> lagrangian_surface_nearest_mean(faces.size(), 0.);
      std::vector<double> lagrangian_surface_nearest_max(faces.size(), 0.);
      std::vector<double> lagrangian_surface_global_fallback(faces.size(), 0.);

#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(faces.size()); ++i) {
        const auto& f = faces[i];
        auto [p0, p1, p2] = f->getPoints();
        const Tddd X0_corr = getPosition(p0);
        const Tddd X1_corr = getPosition(p1);
        const Tddd X2_corr = getPosition(p2);
        const Tddd X0_lag = lag_position(p0);
        const Tddd X1_lag = lag_position(p1);
        const Tddd X2_lag = lag_position(p2);

        const double local_scale =
            std::max({Norm(X0_lag - X1_lag), Norm(X1_lag - X2_lag), Norm(X2_lag - X0_lag), 1e-12});
        const auto local_candidates = collect_local_face_candidates(f);

        double weighted_sum = 0.;
        double weight_sum = 0.;
        double max_dist = 0.;
        bool used_global = false;
        for (const auto& q : dunavant3) {
          const double l0 = q[0], l1 = q[1], l2 = q[2], w = q[3];
          const Tddd X_corr = l0 * X0_corr + l1 * X1_corr + l2 * X2_corr;

          double dist = nearest_distance_to_lagrangian_faces(X_corr, local_candidates);
          if (!(dist <= 2.5 * local_scale)) {
            dist = nearest_distance_to_lagrangian_faces(X_corr, all_faces);
            used_global = true;
          }

          weighted_sum += w * dist;
          weight_sum += w;
          max_dist = std::max(max_dist, dist);
        }

        lagrangian_surface_nearest_mean[i] = (weight_sum > 0.) ? (weighted_sum / weight_sum) : 0.;
        lagrangian_surface_nearest_max[i] = max_dist;
        lagrangian_surface_global_fallback[i] = used_global ? 1. : 0.;
      }

      for (std::size_t i = 0; i < faces.size(); ++i) {
        auto* f = faces[i];
        F_lagrangian_surface_nearest_mean[f] = lagrangian_surface_nearest_mean[i];
        F_lagrangian_surface_nearest_max[f] = lagrangian_surface_nearest_max[i];
        F_lagrangian_surface_global_fallback[f] = lagrangian_surface_global_fallback[i];
      }

      for (auto* f : faces) {
        if (!f)
          continue;
        double sharp_edges = 0.;
        double sharp_points = 0.;
        double direct_sharp_edges = 0.;
        double max_edge_angle = 0.;
        for (auto* l : f->getLines())
          if (l && l->SharpQ(output_sharpq_feature_angle))
            sharp_edges += 1.;
        for (auto* l : f->getLines()) {
          const double edge_angle = line_max_boundary_normal_angle_deg(l);
          max_edge_angle = std::max(max_edge_angle, edge_angle);
          if (edge_angle > output_sharpq_feature_angle_deg)
            direct_sharp_edges += 1.;
        }
        for (auto* p : f->getPoints())
          if (p && p->SharpQ(output_sharpq_feature_angle))
            sharp_points += 1.;
        F_SharpQ_face_edge_count[f] = sharp_edges;
        F_SharpQ_face_point_count[f] = sharp_points;
        F_SharpQ_face_has_any[f] = (sharp_edges > 0. || sharp_points > 0.) ? 1. : 0.;
        F_SharpQ_face_direct_edge_count[f] = direct_sharp_edges;
        F_SharpQ_face_max_edge_angle_deg[f] = max_edge_angle;
        F_SharpQ_face_direct_has_any[f] = direct_sharp_edges > 0. ? 1. : 0.;
      }
    }

    std::atomic<int> output_eval_fail_count{0};
    try {
#pragma omp parallel for
      for (const auto& p : ToVector(water->getPoints())) {

        auto [f, _, __] = getEffectiveNearestContactFace(p);

        if (f)
          P_velocity_body[p] = f->getNetwork()->velocityRigidBody(p->X);

        int i = 0;
        for (const auto& [f, d] : p->dofs) {
          if (!f || d.contact_opponent_faces.empty())
            continue;
          auto* cf = d.nearestContactFace();
          auto X = Nearest(p->X, ToX(cf));
          if (i == 0)
            P_V2ContactFaces0[p] = X - p->X;
          if (i == 1)
            P_V2ContactFaces1[p] = X - p->X;
          if (i == 2)
            P_V2ContactFaces2[p] = X - p->X;
          if (i == 3)
            P_V2ContactFaces3[p] = X - p->X;
          if (i == 4)
            P_V2ContactFaces4[p] = X - p->X;
          if (i == 5)
            P_V2ContactFaces5[p] = X - p->X;
          i++;
        }

        try {
          P_accelNeumann[p] = accelNeumann(p);
          P_uNeumann[p] = contactNormalVelocity(p);
          P_phin_Dirichlet[p] = p->getNormalDirichlet_BEM() * p->phin_Dirichlet;
          P_U_shift_BEM[p] = p->vecToSurface;
          P_isMultipleNode[p] = p->isMultipleNode;
          P_isAbsorbed[p] = p->absorbedBy != nullptr;
          P_isAbsorbed_SDF[p] = p->signed_distance;
          P_u_absorbed[p] = p->u_absorbed;
          P_phi_absorbed[p] = p->phi_absorbed;
          P_absorb_gamma[p] = p->absorb_gamma;
          if (p->penetratedBody) {
            auto [near_f, near_X] = p->penetratedBody->Nearest(p->X);
            if (near_f)
              P_penetration_dist[p] = Norm(p->X - near_X);
            else
              P_penetration_dist[p] = 0.;
          } else {
            P_penetration_dist[p] = 0.;
          }

          P_phi[p] = std::get<0>(p->phiphin);
          P_phin[p] = std::get<1>(p->phiphin);
          P_phi_t[p] = std::get<0>(p->phiphin_t);
          P_phin_t[p] = std::get<1>(p->phiphin_t);
          P_phin_t_from_Hessian[p] = phint_Neumann(p);
          P_normal_BEM[p] = p->getNormal_BEM();
          P_ContactFaces[p] = (double)getEffectiveContactFaces(p).size();
          double maxdist = 0;
          // Compute max contact distance from dofs
          for (const auto& [f, d] : p->dofs)
            for (auto* cf : d.contact_opponent_faces) {
              double dist = Norm(p->getPosition() - Nearest(p->getPosition(), ToX(cf)));
              if (dist > maxdist)
                maxdist = dist;
            }
          P_ContactFaces_maxdist[p] = maxdist;
          P_ContactRange[p] = p->contact_range;
          P_facesNeuamnn[p] = (double)p->getFacesNeumann().size();
          {
            double pf_neumann_count = 0.;
            double pf_dirichlet_count = 0.;
            double pf_pressure_min = std::numeric_limits<double>::infinity();
            double pf_detach_count_max = 0.;
            double pf_pressure_below_count = 0.;
            double pf_pressure_detachment_eligible_count = 0.;
            double pf_detached_by_pressure_count = 0.;
            const auto bfs = p->getBoundaryFaces();
            for (const auto& f : bfs) {
              const auto* d_contact = p->findContactState(f);
              if (d_contact && d_contact->nearestContactFace() != nullptr)
                pf_neumann_count += 1.;
              else
                pf_dirichlet_count += 1.;

              if (const auto* d = p->findActiveBieDof(f)) {
                pf_pressure_min = std::min(pf_pressure_min, d->pressure);
                pf_detach_count_max = std::max(pf_detach_count_max, static_cast<double>(d->detach_negative_count));
                if (d->pressure < detachment_pressure_threshold)
                  pf_pressure_below_count += 1.;
                if (d->detached_by_pressure)
                  pf_detached_by_pressure_count += 1.;
              }
              if (pressureDetachmentEligible(p, f))
                pf_pressure_detachment_eligible_count += 1.;
            }
            P_pf_neumann_count[p] = pf_neumann_count;
            P_pf_dirichlet_count[p] = pf_dirichlet_count;
            P_pf_total_count[p] = static_cast<double>(bfs.size());
            P_pf_pressure_min[p] = std::isfinite(pf_pressure_min) ? pf_pressure_min : 0.;
            P_pf_detach_count_max[p] = pf_detach_count_max;
            P_pf_detach_progress_max[p] = (detachment_consecutive_steps > 0)
                                              ? (pf_detach_count_max / static_cast<double>(detachment_consecutive_steps))
                                              : 0.;
            P_pf_pressure_below_count[p] = pf_pressure_below_count;
            P_pf_pressure_detachment_eligible_count[p] = pf_pressure_detachment_eligible_count;
            P_pf_detached_by_pressure_count[p] = pf_detached_by_pressure_count;
          }
          // boundary condition: 0=BCInterface, 1=multiple Neumann, 2=Neumann, 3=multiple Dirichlet, 4=Dirichlet
          P_BC[p] = boundaryConditionValue(p);
          P_diag[p] = p->diag_coeff_BEM;
          P_direction_info_count[p] = p->debug_direction_info_count;
          P_contact_faces_count[p] = p->debug_contact_faces_count;
          P_body_vertices_count[p] = p->debug_body_vertices_count;
          P_isInContact_pass_count[p] = p->debug_isInContact_pass_count;
          P_SharpQ[p] = p->SharpQ(output_sharpq_feature_angle) ? 1. : 0.;
          P_SharpQ_direct[p] = point_direct_sharp(p) ? 1. : 0.;
          P_SharpQ_incident_edge_angle_max_deg[p] = point_max_incident_edge_angle_deg(p);
          P_position[p] = ToX(p);
          P_pressure[p] = p->pressure_BEM;
          P_DphiDt[p] = p->DphiDt(p->u_reloc, 0.);
          P_gradPhi[p] = p->u_potential_BEM;
          P_u_total[p] = p->u_total;
          P_u_reloc[p] = p->u_reloc;
          P_u_potential_BEM[p] = p->u_potential_BEM;
          P_u_omega_VPM[p] = p->u_omega_VPM;
          Tddd convergence_info = {0., 0., 0.};
          gradPhi(p, convergence_info);
          P_velocity_convergence[p] = convergence_info;
          P_vecToSurface[p] = p->vecToSurface;
          P_clungSurface[p] = p->clungSurface;
          P_solidAngle[p] = p->getSolidAngle();
          P_b_diff_RHS_FMM[p] = p->b_diff_RHS_FMM;
          P_b_RHS_FMM[p] = p->b_RHS_FMM;
          P_b_RHS_Direct[p] = p->b_RHS_Direct;
          P_minDepthFromBCInterface[p] = p->minDepthFromBCInterface;
          P_minDepthFromMultipleNode[p] = p->minDepthFromMultipleNode;
          P_almost_solid_angle[p] = p->almost_solid_angle;
        } catch (const std::exception& e) {
          const int fail_idx = output_eval_fail_count.fetch_add(1) + 1;
          P_position[p] = ToX(p);
          P_phi[p] = std::get<0>(p->phiphin);
          P_phin[p] = std::get<1>(p->phiphin);
          P_velocity_convergence[p] = {1E+30, 1E+30, 1E+30};
          P_BC[p] = boundaryConditionValue(p);
          P_SharpQ[p] = p->SharpQ(output_sharpq_feature_angle) ? 1. : 0.;
          P_SharpQ_direct[p] = point_direct_sharp(p) ? 1. : 0.;
          P_SharpQ_incident_edge_angle_max_deg[p] = point_max_incident_edge_angle_deg(p);
          if (fail_idx <= 5) {
#pragma omp critical
            {
              std::cerr << Yellow << "[Output] field eval failed at p->X=" << p->X
                        << " (Dirichlet=" << p->Dirichlet
                        << ", Neumann=" << p->Neumann
                        << ", BCInterface=" << p->BCInterface << "): "
                        << e.what() << colorReset << std::endl;
            }
          }
        }
      }
      if (const int n_fail = output_eval_fail_count.load(); n_fail > 0) {
        std::cerr << Yellow << "[Output] field evaluation failed at " << n_fail
                  << " points in dataForOutput. Stored sentinel values and continued." << colorReset << std::endl;
      }
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };

    // --- Line midpoint data (true quadratic elements) ---
    if (use_true_quadratic_element) {
      for (auto* l : water->getBoundaryLines()) {
        auto [pA, pB] = l->getPoints();
        // Direct midpoint values
        P_phi[l] = l->phiphin[0];
        if (auto* d0 = l->findActiveBieDof(nullptr))
          P_phin[l] = d0->phin;
        else {
          double wa = 0., wp = 0.;
          for (auto* f : l->getBoundaryFaces())
            if (f) {
              if (auto* df = l->findActiveBieDof(f)) {
                wp += df->phin * f->area;
                wa += f->area;
              }
            }
          P_phin[l] = (wa > 0.) ? wp / wa : l->phiphin[1];
        }
        P_phi_t[l] = l->phiphin_t[0];
        P_phin_t[l] = l->phiphin_t[1];
        P_diag[l] = l->diag_coeff_BEM;
        P_isMultipleNode[l] = l->isMultipleNode;
        P_isAbsorbed[l] = l->absorbedBy != nullptr;
        P_isAbsorbed_SDF[l] = l->signed_distance;
        P_ContactRange[l] = l->contact_range;
        P_direction_info_count[l] = l->debug_direction_info_count;
        P_contact_faces_count[l] = l->debug_contact_faces_count;
        P_body_vertices_count[l] = l->debug_body_vertices_count;
        P_isInContact_pass_count[l] = l->debug_isInContact_pass_count;
        P_SharpQ[l] = l->SharpQ(output_sharpq_feature_angle) ? 1. : 0.;
        P_SharpQ_direct[l] = line_direct_sharp(l) ? 1. : 0.;
        P_SharpQ_incident_edge_angle_max_deg[l] = line_max_boundary_normal_angle_deg(l);
        P_position[l] = l->getPosition();
        P_vecToSurface[l] = l->vecToSurface;
        P_clungSurface[l] = l->clungSurface;
        P_u_reloc[l] = l->u_reloc;
        P_u_potential_BEM[l] = l->u_potential_BEM;
        P_u_total[l] = l->u_total;
        P_u_omega_VPM[l] = l->u_omega_VPM;
        P_u_absorbed[l] = l->u_absorbed;
        P_phi_absorbed[l] = l->phi_absorbed;
        P_absorb_gamma[l] = l->absorb_gamma;
        P_U_shift_BEM[l] = l->vecToSurface;
        {
          auto itA = P_pressure.find(pA);
          auto itB = P_pressure.find(pB);
          P_pressure[l] = (itA != P_pressure.end() && itB != P_pressure.end()) ? 0.5 * (itA->second + itB->second) : 0.;
        }
        P_DphiDt[l] = l->DphiDt(l->u_reloc, 0.);
        P_BC[l] = boundaryConditionValue(l);

        // lf_* fields: computed from line's boundary faces
        {
          double lf_neumann = 0., lf_dirichlet = 0.;
          double lf_pressure_min = std::numeric_limits<double>::infinity();
          double lf_detach_count_max = 0., lf_pressure_below = 0.;
          double lf_eligible = 0., lf_detached = 0.;
          for (auto* f : l->getBoundaryFaces()) {
            if (f && isNeumannBoundaryState(l, f))
              lf_neumann += 1.;
            if (f && isDirichletBoundaryState(l, f))
              lf_dirichlet += 1.;
            if (auto* d = l->findActiveBieDof(f)) {
              lf_pressure_min = std::min(lf_pressure_min, d->pressure);
              lf_detach_count_max = std::max(lf_detach_count_max, static_cast<double>(d->detach_negative_count));
              if (d->pressure < detachment_pressure_threshold)
                lf_pressure_below += 1.;
              if (d->detached_by_pressure)
                lf_detached += 1.;
            }
            if (f && pressureDetachmentEligible(l, f))
              lf_eligible += 1.;
          }
          P_lf_neumann_count[l] = lf_neumann;
          P_lf_dirichlet_count[l] = lf_dirichlet;
          P_lf_total_count[l] = static_cast<double>(l->getBoundaryFaces().size());
          P_lf_pressure_min[l] = std::isfinite(lf_pressure_min) ? lf_pressure_min : 0.;
          P_lf_detach_count_max[l] = lf_detach_count_max;
          P_lf_detach_progress_max[l] = (detachment_consecutive_steps > 0) ? (lf_detach_count_max / static_cast<double>(detachment_consecutive_steps)) : 0.;
          P_lf_pressure_below_count[l] = lf_pressure_below;
          P_lf_pressure_detachment_eligible_count[l] = lf_eligible;
          P_lf_detached_by_pressure_count[l] = lf_detached;
        }

        // pf_* are point-only: keep sentinel for lines
        P_pf_neumann_count[l] = 0.;
        P_pf_dirichlet_count[l] = 0.;
        P_pf_total_count[l] = 0.;
        P_pf_pressure_min[l] = 0.;
        P_pf_detach_count_max[l] = 0.;
        P_pf_detach_progress_max[l] = 0.;
        P_pf_pressure_below_count[l] = 0.;
        P_pf_pressure_detachment_eligible_count[l] = 0.;
        P_pf_detached_by_pressure_count[l] = 0.;

        // Fields that need endpoint averaging (gradPhi, velocity_convergence, etc.)
        auto avg_d = [&](uomap_DOF_d& m) {
          auto itA = m.find(pA);
          auto itB = m.find(pB);
          if (itA != m.end() && itB != m.end() && isFinite(itA->second) && isFinite(itB->second))
            m[l] = 0.5 * (itA->second + itB->second);
          else
            m[l] = 0.;
        };
        auto avg_v = [&](uomap_DOF_Tddd& m) {
          auto itA = m.find(pA);
          auto itB = m.find(pB);
          if (itA != m.end() && itB != m.end())
            m[l] = 0.5 * (itA->second + itB->second);
          else
            m[l] = Tddd{0., 0., 0.};
        };
        avg_d(P_solidAngle);
        avg_d(P_b_diff_RHS_FMM);
        avg_d(P_b_RHS_FMM);
        avg_d(P_b_RHS_Direct);
        avg_d(P_ContactFaces);
        avg_d(P_ContactFaces_maxdist);
        avg_d(P_penetration_dist);
        avg_d(P_minDepthFromBCInterface);
        avg_d(P_minDepthFromMultipleNode);
        avg_d(P_almost_solid_angle);
        avg_d(P_phin_t_from_Hessian);
        avg_d(P_facesNeuamnn);
        avg_v(P_accelNeumann);
        avg_v(P_uNeumann);
        avg_v(P_phin_Dirichlet);
        avg_v(P_gradPhi);
        avg_v(P_velocity_convergence);
        avg_v(P_normal_BEM);
      }
    }

    try {
      VV_VarForOutput data = {
          //  {"body accel", P_accel_body},
          //  {"body velocity", P_velocity_body},
          //  {"accelNeumann", P_accelNeumann},
          //  {"uNeumann", P_uNeumann},
          {"isMultipleNode", P_isMultipleNode},
          {"isAbsorbed", P_isAbsorbed},
          {"isAbsorbed_SDF", P_isAbsorbed_SDF},
          {"u_total", P_u_total},
          {"u_reloc", P_u_reloc},
          {"u_potential_BEM", P_u_potential_BEM},
          {"u_omega_VPM", P_u_omega_VPM},
          //  {"u_potential_BEM", P_U_BEM},
          {"u_absorbed", P_u_absorbed},
          {"φ_absorb", P_phi_absorbed},
          {"γ_absorb", P_absorb_gamma},
          {"U_shift_BEM", P_U_shift_BEM},
          {"ContactFaces", P_ContactFaces},
          {"ContactFaces_maxdist", P_ContactFaces_maxdist},
          {"ContactRange", P_ContactRange},
          //  {"faces Neuamnn", P_facesNeuamnn},
          {"grad_phi", P_gradPhi},
          {"velocity_convergence", P_velocity_convergence},
          {"vecToSurface", P_vecToSurface},
          {"clungSurface", P_clungSurface},
          //  {"solidAngle", P_solidAngle},
          {"b_diff_RHS_FMM", P_b_diff_RHS_FMM},
          {"b_RHS_FMM", P_b_RHS_FMM},
          {"b_RHS_Direct", P_b_RHS_Direct},
          {"position", P_position},
          {"φ", P_phi},
          {"φn", P_phin},
          {"φt", P_phi_t},
          {"φnt", P_phin_t},
          {"diag", P_diag},
          {"penetration_distance", P_penetration_dist},
          //  {"φnt hess", P_phin_t_from_Hessian},
          {"boundary condition", P_BC},
          {"pf_neumann_count", P_pf_neumann_count},
          {"pf_dirichlet_count", P_pf_dirichlet_count},
          {"pf_total_count", P_pf_total_count},
          {"pf_pressure_min", P_pf_pressure_min},
          {"pf_detach_count_max", P_pf_detach_count_max},
          {"pf_detach_progress_max", P_pf_detach_progress_max},
          {"pf_pressure_below_count", P_pf_pressure_below_count},
          {"pf_pressure_detachment_eligible_count", P_pf_pressure_detachment_eligible_count},
          {"pf_detached_by_pressure_count", P_pf_detached_by_pressure_count},
          {"lf_neumann_count", P_lf_neumann_count},
          {"lf_dirichlet_count", P_lf_dirichlet_count},
          {"lf_total_count", P_lf_total_count},
          {"lf_pressure_min", P_lf_pressure_min},
          {"lf_detach_count_max", P_lf_detach_count_max},
          {"lf_detach_progress_max", P_lf_detach_progress_max},
          {"lf_pressure_below_count", P_lf_pressure_below_count},
          {"lf_pressure_detachment_eligible_count", P_lf_pressure_detachment_eligible_count},
          {"lf_detached_by_pressure_count", P_lf_detached_by_pressure_count},
          {"lagrangian_surface_nearest_mean", F_lagrangian_surface_nearest_mean},
          {"lagrangian_surface_nearest_max", F_lagrangian_surface_nearest_max},
          {"lagrangian_surface_global_fallback", F_lagrangian_surface_global_fallback},
          {"direction_info_count", P_direction_info_count},
          {"contact_faces_count", P_contact_faces_count},
          {"body_vertices_count", P_body_vertices_count},
          {"isInContact_pass_count", P_isInContact_pass_count},
          {"SharpQ", P_SharpQ},
          {"SharpQ_direct", P_SharpQ_direct},
          {"SharpQ_incident_edge_angle_max_deg", P_SharpQ_incident_edge_angle_max_deg},
          {"SharpQ_face_edge_count", F_SharpQ_face_edge_count},
          {"SharpQ_face_point_count", F_SharpQ_face_point_count},
          {"SharpQ_face_has_any", F_SharpQ_face_has_any},
          {"SharpQ_face_direct_edge_count", F_SharpQ_face_direct_edge_count},
          {"SharpQ_face_max_edge_angle_deg", F_SharpQ_face_max_edge_angle_deg},
          {"SharpQ_face_direct_has_any", F_SharpQ_face_direct_has_any},
          {"pressure", P_pressure},
          //  {"P_V2ContactFaces0", P_V2ContactFaces0},
          //  {"P_V2ContactFaces1", P_V2ContactFaces1},
          //  {"P_V2ContactFaces2", P_V2ContactFaces2},
          //  {"P_V2ContactFaces3", P_V2ContactFaces3},
          //  {"P_V2ContactFaces4", P_V2ContactFaces4},
          //  {"P_V2ContactFaces5", P_V2ContactFaces5},
          //  {"P_minDepthFromBCInterface", P_minDepthFromBCInterface},
          //  {"P_minDepthFromMultipleNode", P_minDepthFromMultipleNode},
          //  {"P_almost_solid_angle", P_almost_solid_angle}
      };
      return data;
    } catch (const error_message&) {
      throw;
    } catch (const std::exception& e) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
    };
    std::cout << __PRETTY_FUNCTION__ << " done" << std::endl;
    return {};
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  }
};

/* ------------------------------------------------------ */

// メッシュ運動制約: 辺方向の相対接近速度で edge crossing を防止
double dt_mesh_motion(const Network& water, double min_dt, const double c = 0.5) {
  for (const auto& p : water.getPoints()) {
    for (const auto& q : p->getNeighbors()) {
      const Tddd dx = q->X - p->X;
      const double dist = Norm(dx);
      if (dist <= 0)
        continue;
      const Tddd edge_dir = dx / dist;
      const Tddd du = p->u_total - q->u_total;
      const double closing = Dot(du, edge_dir);
      if (closing <= 0 || !isFinite(closing))
        continue;
      const double dt = c * dist / closing;
      if (isFinite(dt) && dt < min_dt)
        min_dt = dt;
    }
  }
  return min_dt;
}

// 幾何品質制約: 面の altitude が急激に縮まないように dt を制限
// c = 許容縮み率（1ステップで altitude が最大 c * h_alt だけ縮む）
double dt_geometry(const Network& water, double min_dt, const double c = 0.2) {

  // 既に tiny な sub-triangle の h_alt を基準にする最小閾値
  // これ以下の h_alt は remesh で対処すべきであり、dt 削減では改善しない
  const double h_alt_floor = 1e-4;

  auto check_triangle = [&](const Tddd& X0, const Tddd& X1, const Tddd& X2,
                            const Tddd& U0, const Tddd& U1, const Tddd& U2) {
    const double area = TriangleArea(X0, X1, X2);
    if (area <= 0 || !isFinite(area))
      return;
    const Tddd face_normal = Normalize(Cross(X1 - X0, X2 - X0));

    const std::array<Tddd, 3> X = {X0, X1, X2};
    const std::array<Tddd, 3> U = {U0, U1, U2};

    for (int i = 0; i < 3; i++) {
      const int a = i, b = (i + 1) % 3, opp = (i + 2) % 3;
      const Tddd edge_vec = X[b] - X[a];
      const double edge_len = Norm(edge_vec);
      if (edge_len <= 0)
        continue;

      const double h_alt = 2.0 * area / edge_len;
      if (h_alt < h_alt_floor)
        continue; // 既に tiny — remesh で対処

      const Tddd edge_dir = edge_vec / edge_len;
      const Tddd alt_dir = Normalize(Cross(face_normal, edge_dir));
      const Tddd to_opp = X[opp] - 0.5 * (X[a] + X[b]);
      const double sign = Dot(to_opp, alt_dir);

      const Tddd u_rel = U[opp] - 0.5 * (U[a] + U[b]);
      const double closing_speed = -Dot(u_rel, alt_dir) * (sign > 0 ? 1.0 : -1.0);

      if (closing_speed <= 0 || !isFinite(closing_speed))
        continue;

      const double dt = c * h_alt / closing_speed;
      if (isFinite(dt) && dt < min_dt) {
        min_dt = dt;
        if (bem_verbose())
          std::cout << "[dt_geometry] h_alt=" << h_alt << " closing=" << closing_speed << " dt=" << dt << std::endl;
      }
    }
  };

  for (const auto& f : water.getBoundaryFaces()) {
    if (!f)
      continue;
    auto [p0, p1, p2] = f->getPoints();
    if (!p0 || !p1 || !p2)
      continue;

    if (f->isTrueQuadraticElement) {
      auto [l0, l1, l2] = f->getLines();
      if (!l0 || !l1 || !l2)
        continue;
      check_triangle(p0->X, l0->X_mid, l2->X_mid, p0->u_total, l0->u_total, l2->u_total);
      check_triangle(l0->X_mid, p1->X, l1->X_mid, l0->u_total, p1->u_total, l1->u_total);
      check_triangle(l2->X_mid, l1->X_mid, p2->X, l2->u_total, l1->u_total, p2->u_total);
      check_triangle(l0->X_mid, l1->X_mid, l2->X_mid, l0->u_total, l1->u_total, l2->u_total);
    } else {
      check_triangle(p0->X, p1->X, p2->X, p0->u_total, p1->u_total, p2->u_total);
    }
  }
  return min_dt;
}

void show_info(const Network& net) {
  int total = 0, total_c_face = 0, c = 0, n = 0, d = 0;
  for (const auto& p : net.getPoints()) {
    total++;
    if (p->BCInterface) {
      c++;
      total_c_face += p->getBoundaryFaces().size();
    } else if (p->Neumann)
      n++;
    else if (p->Dirichlet)
      d++;
  }
  std::cout << "net.getPoints() = " << net.getPoints().size() << std::endl;
  std::cout << "Total : " << total << std::endl;
  std::cout << "Total variables: " << d + n + total_c_face << std::endl;
  int doublenode = total - c + total_c_face;
  std::cout << "Total case double-node : " << doublenode << std::endl;
  std::cout << "node reduction : " << (double)(doublenode - total) / (double)doublenode << std::endl;
  std::cout << "BCInterface : " << c << std::endl;
  std::cout << "Total BCInterface faces : " << total_c_face << std::endl;
  std::cout << "Neumann : " << n << std::endl;
  std::cout << "Dirichlet : " << d << std::endl;
};

// b! ------------------------------------------------------ */

#endif
