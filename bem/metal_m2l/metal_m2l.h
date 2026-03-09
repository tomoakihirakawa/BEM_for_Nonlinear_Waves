/*
 * metal_m2l.h
 * C API for Metal-accelerated M2L (Multipole-to-Local) transformation
 *
 * This library is compiled with Apple Clang and linked as a dynamic library.
 * The main BEM code (compiled with GCC) calls these functions via C linkage.
 *
 * Precision: Float + Kahan summation (~7 digits with error compensation)
 * For full double precision (~15 digits), use CPU-based M2L instead.
 */

#ifndef METAL_M2L_H
#define METAL_M2L_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle to the Metal M2L context */
typedef struct MetalM2LContext MetalM2LContext;

/*
 * Initialize the Metal M2L context.
 * Parameters:
 *   threadgroup_mode - 0: default kernel, 1: threadgroup parallel kernel
 * Returns NULL on failure (e.g., Metal not available).
 */
MetalM2LContext* metal_m2l_init(int threadgroup_mode);

/*
 * Set up the M2L data structure (float precision with Kahan summation).
 * This should be called once when FMM structure is initialized.
 *
 * Parameters:
 *   ctx             - Context handle
 *   num_buckets     - Total number of buckets (for MM array indexing)
 *   num_coeffs      - Number of coefficients per bucket (typically (N+1)^2)
 *   num_rows        - Total number of M2L rows (all buckets combined)
 *   num_terms       - Total number of M2L terms (all rows combined)
 *
 * Per-row data:
 *   row_offset[]    - Offset into terms array for each row
 *   row_len[]       - Number of terms for each row
 *
 * Per-term data:
 *   term_coeff_re[] - Real part of AAAY coefficient (float)
 *   term_coeff_im[] - Imaginary part of AAAY coefficient (float)
 *   term_src_idx[]  - Source index: bucket_id * num_coeffs + coeff_idx
 *
 * Returns 0 on success, non-zero on error.
 */
int metal_m2l_setup(
    MetalM2LContext* ctx,
    int32_t num_buckets,
    int32_t num_coeffs,
    int32_t num_rows,
    int32_t num_terms,
    const int32_t* row_offset,
    const int32_t* row_len,
    const float* term_coeff_re,
    const float* term_coeff_im,
    const int32_t* term_src_idx
);

/*
 * Get direct pointers to GPU buffers for source MM values.
 * Write source MM data here before calling metal_m2l_compute_no_copy().
 *
 * Layout: src_mm[bucket_id * num_coeffs + coeff_idx] for each component
 *
 *   src_mm_re0[idx] = real(MM[0])
 *   src_mm_im0[idx] = imag(MM[0])
 *   src_mm_re1[idx] = real(MM[1])
 *   src_mm_im1[idx] = imag(MM[1])
 */
float* metal_m2l_get_src_mm_re0_buffer(MetalM2LContext* ctx);
float* metal_m2l_get_src_mm_im0_buffer(MetalM2LContext* ctx);
float* metal_m2l_get_src_mm_re1_buffer(MetalM2LContext* ctx);
float* metal_m2l_get_src_mm_im1_buffer(MetalM2LContext* ctx);

/*
 * Get pointer to result buffer.
 * After compute, results are:
 *   result[row * 4 + 0] = sum_re0 (real part of sum for MM[0])
 *   result[row * 4 + 1] = sum_im0 (imag part of sum for MM[0])
 *   result[row * 4 + 2] = sum_re1 (real part of sum for MM[1])
 *   result[row * 4 + 3] = sum_im1 (imag part of sum for MM[1])
 */
float* metal_m2l_get_result_buffer(MetalM2LContext* ctx);

/*
 * Compute M2L transformation.
 * Assumes source MM data has been written to the buffers.
 * Results are accumulated sums for each row.
 *
 * Returns 0 on success, non-zero on error.
 */
int metal_m2l_compute_no_copy(MetalM2LContext* ctx);

/*
 * Asynchronous compute - submit work to GPU and return immediately.
 * Use metal_m2l_wait() to wait for completion before reading results.
 */
int metal_m2l_compute_async(MetalM2LContext* ctx);

/*
 * Wait for async compute to complete.
 */
void metal_m2l_wait(MetalM2LContext* ctx);

/*
 * Release all resources.
 */
void metal_m2l_destroy(MetalM2LContext* ctx);

/*
 * Check if Metal M2L is available on this system.
 * Returns 1 if available, 0 otherwise.
 */
int metal_m2l_is_available(void);

/*
 * Get the name of the Metal device being used.
 * Returns a static string (do not free).
 */
const char* metal_m2l_device_name(MetalM2LContext* ctx);

#ifdef __cplusplus
}
#endif

#endif /* METAL_M2L_H */
