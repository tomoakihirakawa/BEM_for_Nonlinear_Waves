/*
 * metal_m2l.metal
 * Metal compute kernels for M2L (Multipole-to-Local) transformation
 *
 * Precision: Float + Kahan summation (~7 digits with error compensation)
 * For full double precision (~15 digits), use CPU-based M2L instead.
 *
 * M2L computation:
 *   For each row (Local expansion coefficient):
 *     sum0 = Σ AAAY[i] * src_MM[i][0]   (complex multiplication)
 *     sum1 = Σ AAAY[i] * src_MM[i][1]
 *
 * Complex multiplication: (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
 */

#include <metal_stdlib>
using namespace metal;

/* ============================================================================
 * Float precision with Kahan summation
 * ============================================================================ */

/*
 * M2L kernel with Kahan summation for improved precision.
 * Each thread processes one row (one Local expansion coefficient).
 *
 * Input buffers:
 *   src_mm_re0, src_mm_im0: real and imag parts of source MM[0]
 *   src_mm_re1, src_mm_im1: real and imag parts of source MM[1]
 *   coeff_re, coeff_im: AAAY coefficients (real and imag)
 *   term_src_idx: index into src_mm arrays for each term
 *   row_offset: offset into term arrays for each row
 *   row_len: number of terms for each row
 *
 * Output buffer:
 *   results: [row * 4 + 0..3] = (sum_re0, sum_im0, sum_re1, sum_im1)
 */

/*
 * Kahan summation helper using float4 SIMD
 * sum = (sum_re0, sum_im0, sum_re1, sum_im1)
 * c = compensation terms
 */
inline void kahan_add4(thread float4& sum, thread float4& c, float4 val) {
    float4 y = val - c;
    float4 t = sum + y;
    c = (t - sum) - y;
    sum = t;
}

/*
 * Process one term and accumulate with Kahan summation
 * Returns float4(prod_re0, prod_im0, prod_re1, prod_im1)
 */
inline float4 process_term(
    device const float* src_mm_re0,
    device const float* src_mm_im0,
    device const float* src_mm_re1,
    device const float* src_mm_im1,
    float ar, float ai, int32_t src_idx)
{
    float b0r = src_mm_re0[src_idx];
    float b0i = src_mm_im0[src_idx];
    float b1r = src_mm_re1[src_idx];
    float b1i = src_mm_im1[src_idx];

    return float4(
        ar * b0r - ai * b0i,  // prod_re0
        ar * b0i + ai * b0r,  // prod_im0
        ar * b1r - ai * b1i,  // prod_re1
        ar * b1i + ai * b1r   // prod_im1
    );
}

kernel void m2l_float_kahan(
    device const float* src_mm_re0   [[buffer(0)]],
    device const float* src_mm_im0   [[buffer(1)]],
    device const float* src_mm_re1   [[buffer(2)]],
    device const float* src_mm_im1   [[buffer(3)]],
    device const float* coeff_re     [[buffer(4)]],
    device const float* coeff_im     [[buffer(5)]],
    device const int32_t* term_src_idx [[buffer(6)]],
    device const int32_t* row_offset [[buffer(7)]],
    device const int32_t* row_len    [[buffer(8)]],
    device float* results            [[buffer(9)]],
    uint row_id                      [[thread_position_in_grid]])
{
    int32_t offset = row_offset[row_id];
    int32_t len = row_len[row_id];

    // SIMD float4 for Kahan summation: (sum_re0, sum_im0, sum_re1, sum_im1)
    float4 sum = float4(0.0f);
    float4 c = float4(0.0f);  // Kahan compensation

    // 4x unrolled main loop
    int32_t i = 0;
    for (; i + 3 < len; i += 4) {
        int32_t idx0 = offset + i;

        // Prefetch indices
        int32_t src_idx0 = term_src_idx[idx0];
        int32_t src_idx1 = term_src_idx[idx0 + 1];
        int32_t src_idx2 = term_src_idx[idx0 + 2];
        int32_t src_idx3 = term_src_idx[idx0 + 3];

        // Load coefficients
        float ar0 = coeff_re[idx0], ai0 = coeff_im[idx0];
        float ar1 = coeff_re[idx0 + 1], ai1 = coeff_im[idx0 + 1];
        float ar2 = coeff_re[idx0 + 2], ai2 = coeff_im[idx0 + 2];
        float ar3 = coeff_re[idx0 + 3], ai3 = coeff_im[idx0 + 3];

        // Process 4 terms with SIMD accumulation
        float4 prod0 = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1, ar0, ai0, src_idx0);
        kahan_add4(sum, c, prod0);

        float4 prod1 = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1, ar1, ai1, src_idx1);
        kahan_add4(sum, c, prod1);

        float4 prod2 = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1, ar2, ai2, src_idx2);
        kahan_add4(sum, c, prod2);

        float4 prod3 = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1, ar3, ai3, src_idx3);
        kahan_add4(sum, c, prod3);
    }

    // Handle remaining terms (0-3 iterations)
    for (; i < len; ++i) {
        int32_t idx = offset + i;
        int32_t src_idx = term_src_idx[idx];
        float ar = coeff_re[idx];
        float ai = coeff_im[idx];

        float4 prod = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1, ar, ai, src_idx);
        kahan_add4(sum, c, prod);
    }

    // Write SIMD results
    results[row_id * 4 + 0] = sum.x;
    results[row_id * 4 + 1] = sum.y;
    results[row_id * 4 + 2] = sum.z;
    results[row_id * 4 + 3] = sum.w;
}

/* ============================================================================
 * Float precision with threadgroup reduction (row-parallel)
 * ============================================================================ */

#define M2L_TG_MAX 256u

/*
 * One threadgroup handles one row; threads within the group share the work.
 * Partial sums are reduced in threadgroup memory.
 */
kernel void m2l_float_kahan_tg(
    device const float* src_mm_re0   [[buffer(0)]],
    device const float* src_mm_im0   [[buffer(1)]],
    device const float* src_mm_re1   [[buffer(2)]],
    device const float* src_mm_im1   [[buffer(3)]],
    device const float* coeff_re     [[buffer(4)]],
    device const float* coeff_im     [[buffer(5)]],
    device const int32_t* term_src_idx [[buffer(6)]],
    device const int32_t* row_offset [[buffer(7)]],
    device const int32_t* row_len    [[buffer(8)]],
    device float* results            [[buffer(9)]],
    uint tid                         [[thread_index_in_threadgroup]],
    uint tg_id                       [[threadgroup_position_in_grid]],
    uint tg_size                     [[threads_per_threadgroup]])
{
    const uint row_id = tg_id;
    int32_t offset = row_offset[row_id];
    int32_t len = row_len[row_id];

    // SIMD float4 for Kahan summation
    float4 sum = float4(0.0f);
    float4 c = float4(0.0f);

    // Each thread processes elements with stride = tg_size
    // 4x unroll for threads that have enough work
    int32_t stride = (int32_t)tg_size;
    int32_t i = (int32_t)tid;

    // 4x unrolled loop (each iteration jumps by 4*stride)
    for (; i + 3 * stride < len; i += 4 * stride) {
        // Term 0
        {
            int32_t idx = offset + i;
            int32_t src_idx = term_src_idx[idx];
            float4 prod = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1,
                                       coeff_re[idx], coeff_im[idx], src_idx);
            kahan_add4(sum, c, prod);
        }
        // Term 1
        {
            int32_t idx = offset + i + stride;
            int32_t src_idx = term_src_idx[idx];
            float4 prod = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1,
                                       coeff_re[idx], coeff_im[idx], src_idx);
            kahan_add4(sum, c, prod);
        }
        // Term 2
        {
            int32_t idx = offset + i + 2 * stride;
            int32_t src_idx = term_src_idx[idx];
            float4 prod = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1,
                                       coeff_re[idx], coeff_im[idx], src_idx);
            kahan_add4(sum, c, prod);
        }
        // Term 3
        {
            int32_t idx = offset + i + 3 * stride;
            int32_t src_idx = term_src_idx[idx];
            float4 prod = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1,
                                       coeff_re[idx], coeff_im[idx], src_idx);
            kahan_add4(sum, c, prod);
        }
    }

    // Handle remaining iterations
    for (; i < len; i += stride) {
        int32_t idx = offset + i;
        int32_t src_idx = term_src_idx[idx];
        float4 prod = process_term(src_mm_re0, src_mm_im0, src_mm_re1, src_mm_im1,
                                   coeff_re[idx], coeff_im[idx], src_idx);
        kahan_add4(sum, c, prod);
    }

    // Use float4 threadgroup memory for efficient reduction
    // Include Kahan compensation term in final partial sum
    threadgroup float4 tg_sum[M2L_TG_MAX];
    tg_sum[tid] = sum + c;  // Add compensation before reduction

    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Parallel reduction using float4
    for (uint s = tg_size >> 1; s > 0; s >>= 1) {
        if (tid < s) {
            tg_sum[tid] += tg_sum[tid + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    if (tid == 0) {
        float4 final_sum = tg_sum[0];
        results[row_id * 4 + 0] = final_sum.x;
        results[row_id * 4 + 1] = final_sum.y;
        results[row_id * 4 + 2] = final_sum.z;
        results[row_id * 4 + 3] = final_sum.w;
    }
}

