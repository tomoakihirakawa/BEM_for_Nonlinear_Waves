/*
 * metal_nearfield.h
 * C API for Metal-accelerated nearfield direct integration (setDirectIntegration)
 *
 * This library computes near-field G/Gn kernel weights via Gauss quadrature
 * on the GPU. It replaces the CPU-based setDirectIntegration for linear elements.
 *
 * Architecture:
 *   - Targets are grouped by FMM bucket (cell group). All targets in a group
 *     share the same near-field source face list.
 *   - Each GPU thread processes one target, iterating over all faces in its group.
 *   - Adjacent (target, face) pairs are skipped on GPU and handled on CPU.
 *   - Results are accumulated per local DOF with Kahan summation (float + compensation).
 *   - CPU reads back results and converts to sparse format (near_indices/weights/RLE).
 *
 * Precision: Float + Kahan summation (~7 digits with error compensation)
 *
 * Compiled with Apple Clang, linked as dynamic library.
 * Main BEM code (GCC) calls via C linkage.
 */

#ifndef METAL_NEARFIELD_H
#define METAL_NEARFIELD_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle to the Metal nearfield context */
typedef struct MetalNearfieldContext MetalNearfieldContext;

/*
 * GPU-side face data structure.
 * Must match FaceDataGPU in metal_nearfield.metal.
 * 88 bytes per face, all 4-byte aligned.
 */
typedef struct {
    float X0[3], X1[3], X2[3]; /* vertex coordinates (36 bytes) */
    float J_det;                 /* Jacobian = 2*area (4 bytes) */
    float cross_vec[3];          /* 2*area*normal = (X1-X0)x(X2-X0) (12 bytes) */
    float Xc[3];                 /* centroid (12 bytes) */
    int32_t local_dof[3];        /* local DOF indices within cell group (12 bytes) */
    int32_t vtx_id[3];           /* vertex IDs for adjacency check (12 bytes) */
} MetalFaceData;                 /* Total: 88 bytes */

/*
 * Check if Metal nearfield is available.
 * Returns 1 if available, 0 otherwise.
 */
int metal_nearfield_is_available(void);

/*
 * Initialize the Metal nearfield context.
 * Returns NULL on failure.
 */
MetalNearfieldContext* metal_nearfield_init(void);

/*
 * Setup the nearfield computation.
 *
 * Parameters:
 *   ctx               - Context handle
 *   num_targets       - Total number of targets
 *   num_faces         - Total number of face entries (flattened across all groups)
 *   num_cell_groups   - Number of cell groups
 *   face_data         - MetalFaceData array [num_faces]
 *   group_face_offset - Face range per group [num_cell_groups + 1]
 *   group_dof_count   - Unique DOFs per group [num_cell_groups]
 *   target_pos        - Target positions [num_targets * 3]
 *   target_group_id   - Cell group for each target [num_targets]
 *   target_vtx_id     - Vertex ID per target for adjacency [num_targets]
 *   target_dof_offset - Result buffer offset per target [num_targets]
 *   total_result_slots - Total size of result accumulators
 *   near_region_sq    - Square of near/far threshold
 *
 * Returns 0 on success, non-zero on error.
 */
int metal_nearfield_setup(
    MetalNearfieldContext* ctx,
    int32_t num_targets,
    int32_t num_faces,
    int32_t num_cell_groups,
    const MetalFaceData* face_data,
    const int32_t* group_face_offset,
    const int32_t* group_dof_count,
    const float* target_pos,
    const int32_t* target_group_id,
    const int32_t* target_vtx_id,
    const int32_t* target_dof_offset,
    int32_t total_result_slots,
    float near_region_sq
);

/*
 * Update face geometry (for ALE mesh movement).
 * Only updates geometric fields (X0, X1, X2, J_det, cross_vec, Xc).
 * local_dof and vtx_id remain unchanged.
 */
int metal_nearfield_update_faces(
    MetalNearfieldContext* ctx,
    const MetalFaceData* face_data
);

/*
 * Update target positions (for ALE).
 */
int metal_nearfield_update_targets(
    MetalNearfieldContext* ctx,
    const float* target_pos
);

/*
 * Compute nearfield integration on GPU.
 * Non-adjacent pairs only. Adjacent pairs must be handled on CPU.
 * Results: Kahan-accumulated phi/phin weights per local DOF per target.
 *
 * Returns 0 on success, non-zero on error.
 */
int metal_nearfield_compute(MetalNearfieldContext* ctx);

/*
 * Get result buffer pointers (unified memory, no copy needed).
 * Layout: buffer[target_dof_offset[tid] + local_dof_idx]
 * After compute: phi ≈ phi_hi + phi_lo (Kahan pair)
 */
float* metal_nearfield_get_phi_hi(MetalNearfieldContext* ctx);
float* metal_nearfield_get_phi_lo(MetalNearfieldContext* ctx);
float* metal_nearfield_get_phin_hi(MetalNearfieldContext* ctx);
float* metal_nearfield_get_phin_lo(MetalNearfieldContext* ctx);

/*
 * Release all resources.
 */
void metal_nearfield_destroy(MetalNearfieldContext* ctx);

/*
 * Get Metal device name. Returns static string.
 */
const char* metal_nearfield_device_name(MetalNearfieldContext* ctx);

#ifdef __cplusplus
}
#endif

#endif /* METAL_NEARFIELD_H */
