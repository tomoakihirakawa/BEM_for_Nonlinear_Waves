/*
 * metal_nearfield.metal
 * GPU kernel for nearfield direct integration (linear elements)
 *
 * Each thread processes one target:
 *   1. Iterates over all faces in its cell group
 *   2. Skips adjacent faces (handled on CPU)
 *   3. Computes 25-point (near) or 1-point (far) Gauss quadrature
 *   4. Kahan-accumulates phi and phin weights per local DOF
 *
 * Compile with: metal -c metal_nearfield.metal -o metal_nearfield.air -fno-fast-math
 *   (-fno-fast-math is essential to preserve Kahan summation correctness)
 */

#include <metal_stdlib>
using namespace metal;

/* Face data structure - must match MetalFaceData in metal_nearfield.h (88 bytes) */
struct FaceDataGPU {
    float X0[3];          /* 12 bytes */
    float X1[3];          /* 12 bytes */
    float X2[3];          /* 12 bytes */
    float J_det;          /* 4 bytes */
    float cross_vec[3];   /* 12 bytes */
    float Xc[3];          /* 12 bytes */
    int local_dof[3];     /* 12 bytes */
    int vtx_id[3];        /* 12 bytes */
};

/* Kahan summation: sum += val, compensation stored in comp */
inline void kahan_add(device float& sum, device float& comp, float val) {
    float y = val - comp;
    float t = sum + y;
    comp = (t - sum) - y;
    sum = t;
}

/*
 * Nearfield direct integration kernel (float + Kahan summation)
 *
 * One thread per target. Iterates over all faces in the target's cell group.
 * Skips adjacent faces (target_vtx_id matches any face vertex).
 * Accumulates G/Gn kernel contributions per local DOF with Kahan summation.
 *
 * Near faces (distance < near_region): 25-point Gauss quadrature (ModTriShape<3>)
 * Far faces: 1-point centroid quadrature ({1/3, 1/3, 1/3}, w=0.5)
 */
kernel void nearfield_linear_kahan(
    device const FaceDataGPU* faces         [[buffer(0)]],
    device const int*  group_face_offset    [[buffer(1)]],
    device const int*  group_dof_count      [[buffer(2)]],
    device const float* target_pos          [[buffer(3)]],
    device const int*  target_group_id      [[buffer(4)]],
    device const int*  target_vtx_id        [[buffer(5)]],
    device const int*  target_dof_offset    [[buffer(6)]],
    device float*      result_phi_hi        [[buffer(7)]],
    device float*      result_phi_lo        [[buffer(8)]],
    device float*      result_phin_hi       [[buffer(9)]],
    device float*      result_phin_lo       [[buffer(10)]],
    constant float*    gauss_N              [[buffer(11)]],  /* 25*3 shape functions */
    constant float*    gauss_w              [[buffer(12)]],  /* 25 effective weights */
    constant float&    near_region_sq       [[buffer(13)]],
    constant int&      num_targets          [[buffer(14)]],
    uint tid [[thread_position_in_grid]])
{
    if ((int)tid >= num_targets) return;

    /* Load target data */
    float3 Xt = float3(target_pos[tid * 3],
                        target_pos[tid * 3 + 1],
                        target_pos[tid * 3 + 2]);
    int gid = target_group_id[tid];
    int my_vtx_id = target_vtx_id[tid];
    int face_off = group_face_offset[gid];
    int face_cnt = group_face_offset[gid + 1] - face_off;
    int dof_off = target_dof_offset[tid];
    int dof_cnt = group_dof_count[gid];

    /* Initialize Kahan accumulators */
    for (int d = 0; d < dof_cnt; d++) {
        result_phi_hi[dof_off + d] = 0.0f;
        result_phi_lo[dof_off + d] = 0.0f;
        result_phin_hi[dof_off + d] = 0.0f;
        result_phin_lo[dof_off + d] = 0.0f;
    }

    /* Process each face */
    for (int f = 0; f < face_cnt; f++) {
        FaceDataGPU face = faces[face_off + f];

        /* Adjacency check: skip if target is a vertex of this face */
        if (my_vtx_id == face.vtx_id[0] ||
            my_vtx_id == face.vtx_id[1] ||
            my_vtx_id == face.vtx_id[2]) {
            continue;
        }

        /* Near/far distance check */
        float3 dc = float3(face.Xc[0], face.Xc[1], face.Xc[2]) - Xt;
        bool is_near = dot(dc, dc) < near_region_sq;
        int num_gp = is_near ? 25 : 1;

        /* Face constants */
        float J = face.J_det;
        float3 cr = float3(face.cross_vec[0], face.cross_vec[1], face.cross_vec[2]);

        /* Per-face accumulator */
        float WGN[3] = {0.0f, 0.0f, 0.0f};
        float WGnN[3] = {0.0f, 0.0f, 0.0f};

        for (int gp = 0; gp < num_gp; gp++) {
            float N0, N1, N2, weff;

            if (is_near) {
                N0 = gauss_N[gp * 3];
                N1 = gauss_N[gp * 3 + 1];
                N2 = gauss_N[gp * 3 + 2];
                weff = gauss_w[gp];
            } else {
                /* Far-field: centroid with standard linear shape functions */
                N0 = 1.0f / 3.0f;
                N1 = 1.0f / 3.0f;
                N2 = 1.0f / 3.0f;
                weff = 0.5f;
            }

            /* Source Gauss point position */
            float3 X_gp = float3(
                N0 * face.X0[0] + N1 * face.X1[0] + N2 * face.X2[0],
                N0 * face.X0[1] + N1 * face.X1[1] + N2 * face.X2[1],
                N0 * face.X0[2] + N1 * face.X1[2] + N2 * face.X2[2]
            );

            /* R = X_gp - X_target */
            float3 R = X_gp - Xt;
            float nr2 = dot(R, R);
            float nr_inv = rsqrt(nr2);

            /* G kernel contribution: J_det * w_eff / |R| */
            float w_nr = weff * nr_inv;
            float Jw = J * w_nr;

            /* Gn kernel contribution: -dot(R, cross) * w_eff / |R|^3 */
            float drc = dot(R, cr);
            float neg_dwn3 = -(drc * w_nr * nr_inv * nr_inv);

            /* Accumulate per DOF */
            WGN[0] += Jw * N0;
            WGN[1] += Jw * N1;
            WGN[2] += Jw * N2;

            WGnN[0] += neg_dwn3 * N0;
            WGnN[1] += neg_dwn3 * N1;
            WGnN[2] += neg_dwn3 * N2;
        }

        /* Kahan accumulate to per-target DOF accumulators */
        for (int j = 0; j < 3; j++) {
            int idx = dof_off + face.local_dof[j];
            kahan_add(result_phi_hi[idx], result_phi_lo[idx], WGN[j]);
            kahan_add(result_phin_hi[idx], result_phin_lo[idx], WGnN[j]);
        }
    }
}
