#pragma once

#include <string>
#include <vector>

#include "BEM_time_domain_types.hpp"

extern bool use_linear_element;
extern bool use_pseudo_quadratic_element;
extern bool use_true_quadratic_element;
extern bool use_quadratic_linear_hybrid;

extern NodeRelocationMethod node_relocation_method;
extern NodeRelocationSurface node_relocation_surface;
extern InterpolationMidpointMode interpolation_midpoint_mode;

extern std::string solver_type;
extern std::string coupling_type;
extern double coupling_tol;
extern std::vector<double> coupling_params;
extern std::string preconditioner_type;
extern std::string ilu_neighborhood_type;
extern int ilu_kring_num;
extern double milu_omega;
extern double ilut_drop_tol;
extern int ilut_max_entries_per_row;
extern double ilut_pivot_min;
extern int schwarz_core_k;
extern int schwarz_overlap_k;
extern int schwarz_max_core_size;
extern int schwarz_max_block_size;
extern double schwarz_pivot_min;
extern double schwarz_diag_shift;
extern double solver_tol;
extern int solver_max_iter;
extern int solver_restart;

extern std::string nearfield_mode;
extern int g_p2m_quadrature_points;
extern double g_mac_theta;
extern int fmm_max_level;
extern int fmm_bucket_max_points;

extern int time_step;
extern double simulation_time;

extern bool enable_pressure_detachment;
extern double detachment_pressure_threshold;
extern int detachment_consecutive_steps;

#if defined(USE_METAL_M2L)
extern bool use_metal_m2l;
extern bool metal_m2l_threadgroup;
extern bool metal_m2l_sort_terms;
#endif
