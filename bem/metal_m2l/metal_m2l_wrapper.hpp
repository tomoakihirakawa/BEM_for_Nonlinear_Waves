/*
 * metal_m2l_wrapper.hpp
 *
 * C++ wrapper for Metal M2L (Multipole-to-Local) transformation library.
 * This header can be included in GCC-compiled code.
 * It provides a RAII wrapper and helper functions.
 *
 * Precision modes:
 *   - GPU: Float + Kahan summation (~7 digits with error compensation)
 *   - CPU: Full double precision (~15 digits)
 */

#ifndef METAL_M2L_WRAPPER_HPP
#define METAL_M2L_WRAPPER_HPP

#include "metal_m2l.h"
#include <vector>
#include <complex>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include <numeric>
#include <atomic>
#if defined(_OPENMP)
#include <omp.h>
#endif

class MetalM2LWrapper {
public:
    MetalM2LWrapper(int threadgroup_mode = 0) : ctx_(nullptr), num_buckets_(0), num_coeffs_(0), num_rows_(0) {
        if (!metal_m2l_is_available()) {
            throw std::runtime_error("Metal M2L is not available on this system");
        }
        ctx_ = metal_m2l_init(threadgroup_mode);
        if (ctx_ == nullptr) {
            throw std::runtime_error("Failed to initialize Metal M2L context");
        }
        std::cout << "[Metal M2L] Initialized on device: " << metal_m2l_device_name(ctx_) << std::endl;
    }

    ~MetalM2LWrapper() {
        if (ctx_ != nullptr) {
            metal_m2l_destroy(ctx_);
            ctx_ = nullptr;
        }
    }

    // Non-copyable
    MetalM2LWrapper(const MetalM2LWrapper&) = delete;
    MetalM2LWrapper& operator=(const MetalM2LWrapper&) = delete;

    // Movable
    MetalM2LWrapper(MetalM2LWrapper&& other) noexcept
        : ctx_(other.ctx_), num_buckets_(other.num_buckets_), num_coeffs_(other.num_coeffs_),
          num_rows_(other.num_rows_),
          bucket_to_idx_(std::move(other.bucket_to_idx_)), row_info_(std::move(other.row_info_)) {
        other.ctx_ = nullptr;
    }

    MetalM2LWrapper& operator=(MetalM2LWrapper&& other) noexcept {
        if (this != &other) {
            if (ctx_ != nullptr) {
                metal_m2l_destroy(ctx_);
            }
            ctx_ = other.ctx_;
            num_buckets_ = other.num_buckets_;
            num_coeffs_ = other.num_coeffs_;
            num_rows_ = other.num_rows_;
            bucket_to_idx_ = std::move(other.bucket_to_idx_);
            row_info_ = std::move(other.row_info_);
            other.ctx_ = nullptr;
        }
        return *this;
    }

    /*
     * Setup M2L data structure from FMM buckets.
     * Uses float + Kahan summation for GPU computation.
     *
     * Template parameters:
     *   Bucket - Bucket type with MomentsLocalExpansion containing m2l_rows and m2l_terms
     *
     * Parameters:
     *   all_buckets - Vector of all bucket pointers (for MM indexing)
     *   sort_terms  - If true, sort terms by source index for improved memory locality
     */
    template<typename Bucket>
    void setup(const std::vector<Bucket*>& all_buckets, bool sort_terms = false) {
        num_buckets_ = static_cast<int32_t>(all_buckets.size());

        if (num_buckets_ == 0) {
            throw std::runtime_error("No buckets provided for Metal M2L setup");
        }

        // Determine num_coeffs from first bucket (MM_ is std::array, always has fixed size)
        num_coeffs_ = static_cast<int32_t>(all_buckets[0]->MomentsLocalExpansion.MM_.size());

        // Create bucket pointer -> index mapping (single map, no duplication)
        bucket_to_idx_.clear();
        bucket_to_idx_.reserve(num_buckets_);
        for (int32_t i = 0; i < num_buckets_; ++i) {
            bucket_to_idx_[all_buckets[i]] = i;
        }

        // Step 1: Count rows and terms per bucket (parallel reduction possible)
        std::vector<int32_t> bucket_row_count(num_buckets_);
        std::vector<int32_t> bucket_term_count(num_buckets_);

#if defined(_OPENMP)
        #pragma omp parallel for schedule(static)
#endif
        for (int32_t b = 0; b < num_buckets_; ++b) {
            bucket_row_count[b] = static_cast<int32_t>(all_buckets[b]->MomentsLocalExpansion.m2l_rows.size());
            bucket_term_count[b] = static_cast<int32_t>(all_buckets[b]->MomentsLocalExpansion.m2l_terms.size());
        }

        // Step 2: Compute prefix sums for offsets
        std::vector<int32_t> bucket_row_offset(num_buckets_ + 1);
        std::vector<int32_t> bucket_term_offset(num_buckets_ + 1);
        bucket_row_offset[0] = 0;
        bucket_term_offset[0] = 0;
        for (int32_t b = 0; b < num_buckets_; ++b) {
            bucket_row_offset[b + 1] = bucket_row_offset[b] + bucket_row_count[b];
            bucket_term_offset[b + 1] = bucket_term_offset[b] + bucket_term_count[b];
        }

        int32_t total_rows = bucket_row_offset[num_buckets_];
        int32_t total_terms = bucket_term_offset[num_buckets_];

        num_rows_ = total_rows;
        row_info_.resize(total_rows);

        // Allocate arrays
        std::vector<int32_t> row_offset(total_rows);
        std::vector<int32_t> row_len(total_rows);
        std::vector<int32_t> term_src_idx(total_terms);

        // Coefficient arrays (float only - Kahan summation handles precision)
        std::vector<float> coeff_re(total_terms);
        std::vector<float> coeff_im(total_terms);
        std::atomic<int32_t> missing_src_bucket{0};

        // Step 3: Collect data from all buckets (parallel)
#if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int32_t b = 0; b < num_buckets_; ++b) {
            const auto* bucket = all_buckets[b];
            const auto& local_exp = bucket->MomentsLocalExpansion;
            int32_t bucket_idx = b;

            int32_t row_base = bucket_row_offset[b];
            int32_t term_base = bucket_term_offset[b];
            int32_t local_term_offset = 0;

            for (size_t r = 0; r < local_exp.m2l_rows.size(); ++r) {
                const auto& row = local_exp.m2l_rows[r];
                int32_t row_idx = row_base + static_cast<int32_t>(r);
                int32_t global_term_offset = term_base + local_term_offset;

                // Store row info
                row_info_[row_idx].dst_mm = row.dst_MM;
                row_info_[row_idx].bucket_idx = bucket_idx;

                row_offset[row_idx] = global_term_offset;
                row_len[row_idx] = static_cast<int32_t>(row.len);

                // Copy terms for this row
                const auto* terms = local_exp.m2l_terms.data() + row.offset;
                for (size_t i = 0; i < row.len; ++i) {
                    const auto& term = terms[i];
                    const std::complex<double>& coeff = term.AAAY;

                    // Use pre-stored bucket pointer and coeff index (O(1) lookup via bucket_to_idx_)
                    auto it = bucket_to_idx_.find(term.src_bucket);
                    int32_t src_bucket_idx = 0;
                    if (it != bucket_to_idx_.end()) {
                        src_bucket_idx = it->second;
                    } else {
                        missing_src_bucket.fetch_add(1, std::memory_order_relaxed);
                    }
                    int32_t src_coeff_idx = term.src_coeff_idx;

                    int32_t term_idx = global_term_offset + static_cast<int32_t>(i);
                    term_src_idx[term_idx] = src_bucket_idx * num_coeffs_ + src_coeff_idx;

                    // Filter out NaN/Inf coefficients to prevent GPU computation from producing NaN
                    float re = static_cast<float>(coeff.real());
                    float im = static_cast<float>(coeff.imag());
                    coeff_re[term_idx] = std::isfinite(re) ? re : 0.0f;
                    coeff_im[term_idx] = std::isfinite(im) ? im : 0.0f;
                }

                local_term_offset += static_cast<int32_t>(row.len);
            }
        }

        if (missing_src_bucket.load(std::memory_order_relaxed) > 0) {
            throw std::runtime_error("Metal M2L setup failed: m2l_terms references a source bucket not present in bucket list");
        }

        // Optional: sort terms within each row by source index to improve memory locality on GPU.
        // This changes summation order; enable explicitly when acceptable.
        if (sort_terms) {
            // Parallel sorting: each row is independent
#if defined(_OPENMP)
            #pragma omp parallel
            {
                // Thread-local temporary buffers
                std::vector<int32_t> order;
                std::vector<int32_t> tmp_src;
                std::vector<float> tmp_re, tmp_im;

                order.reserve(2048);
                tmp_src.reserve(2048);
                tmp_re.reserve(2048);
                tmp_im.reserve(2048);

                #pragma omp for schedule(dynamic, 64)
                for (int32_t row = 0; row < num_rows_; ++row) {
                    const int32_t start = row_offset[row];
                    const int32_t len = row_len[row];
                    if (len < 2) continue;

                    order.assign(len, 0);
                    for (int32_t i = 0; i < len; ++i) order[i] = i;

                    std::stable_sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
                        return term_src_idx[start + a] < term_src_idx[start + b];
                    });

                    tmp_src.resize(len);
                    tmp_re.resize(len);
                    tmp_im.resize(len);
                    for (int32_t i = 0; i < len; ++i) {
                        const int32_t idx = start + order[i];
                        tmp_src[i] = term_src_idx[idx];
                        tmp_re[i] = coeff_re[idx];
                        tmp_im[i] = coeff_im[idx];
                    }
                    for (int32_t i = 0; i < len; ++i) {
                        const int32_t idx = start + i;
                        term_src_idx[idx] = tmp_src[i];
                        coeff_re[idx] = tmp_re[i];
                        coeff_im[idx] = tmp_im[i];
                    }
                }
            }
#else
            // Sequential fallback
            std::vector<int32_t> order;
            std::vector<int32_t> tmp_src;
            std::vector<float> tmp_re, tmp_im;

            order.reserve(2048);
            tmp_src.reserve(2048);
            tmp_re.reserve(2048);
            tmp_im.reserve(2048);

            for (int32_t row = 0; row < num_rows_; ++row) {
                const int32_t start = row_offset[row];
                const int32_t len = row_len[row];
                if (len < 2) continue;

                order.assign(len, 0);
                for (int32_t i = 0; i < len; ++i) order[i] = i;

                std::stable_sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
                    return term_src_idx[start + a] < term_src_idx[start + b];
                });

                tmp_src.resize(len);
                tmp_re.resize(len);
                tmp_im.resize(len);
                for (int32_t i = 0; i < len; ++i) {
                    const int32_t idx = start + order[i];
                    tmp_src[i] = term_src_idx[idx];
                    tmp_re[i] = coeff_re[idx];
                    tmp_im[i] = coeff_im[idx];
                }
                for (int32_t i = 0; i < len; ++i) {
                    const int32_t idx = start + i;
                    term_src_idx[idx] = tmp_src[i];
                    coeff_re[idx] = tmp_re[i];
                    coeff_im[idx] = tmp_im[i];
                }
            }
#endif
            std::cout << "[Metal M2L] Terms sorted by src_idx (metal_m2l_sort_terms=true)" << std::endl;
        }

        // Call C API to setup (float + Kahan mode)
        int result = metal_m2l_setup(
            ctx_, num_buckets_, num_coeffs_, num_rows_, total_terms,
            row_offset.data(), row_len.data(),
            coeff_re.data(), coeff_im.data(),
            term_src_idx.data()
        );

        if (result != 0) {
            throw std::runtime_error("Failed to setup Metal M2L");
        }

        std::cout << "[Metal M2L] Setup complete: " << num_buckets_ << " buckets, "
                  << num_coeffs_ << " coeffs, " << num_rows_ << " rows, " << total_terms << " terms"
                  << " [float+Kahan]" << std::endl;

        // Optional row length statistics (enable with BEM_METAL_M2L_STATS=1)
        const char* stats_env = std::getenv("BEM_METAL_M2L_STATS");
        if (stats_env != nullptr && std::atoi(stats_env) != 0 && !row_len.empty()) {
            int32_t min_len = row_len[0];
            int32_t max_len = row_len[0];
            int64_t sum_len = 0;
            for (int32_t v : row_len) {
                if (v < min_len) min_len = v;
                if (v > max_len) max_len = v;
                sum_len += v;
            }
            double mean_len = static_cast<double>(sum_len) / static_cast<double>(row_len.size());

            std::vector<int32_t> row_len_copy = row_len;
            auto percentile = [&](double p) -> int32_t {
                const size_t n = row_len_copy.size();
                if (n == 0) return 0;
                size_t k = static_cast<size_t>(p * (n - 1));
                std::nth_element(row_len_copy.begin(), row_len_copy.begin() + k, row_len_copy.end());
                return row_len_copy[k];
            };
            int32_t p50 = percentile(0.50);
            int32_t p90 = percentile(0.90);
            int32_t p99 = percentile(0.99);

            std::cout << "[Metal M2L] RowLen stats: min=" << min_len
                      << " mean=" << mean_len
                      << " p50=" << p50
                      << " p90=" << p90
                      << " p99=" << p99
                      << " max=" << max_len
                      << std::endl;
        }
    }

    /*
     * Copy source MM values to GPU buffers.
     *
     * Template parameters:
     *   Bucket - Bucket type with MomentsMultipoleExpansion containing MM_
     */
    template<typename Bucket>
    void copySourceMM(const std::vector<Bucket*>& all_buckets) {
        // Validate bucket count matches setup
        if (static_cast<int32_t>(all_buckets.size()) != num_buckets_) {
            std::cerr << "[Metal M2L] Warning: bucket count mismatch in copySourceMM: "
                      << all_buckets.size() << " vs " << num_buckets_ << std::endl;
        }

        float* re0 = metal_m2l_get_src_mm_re0_buffer(ctx_);
        float* im0 = metal_m2l_get_src_mm_im0_buffer(ctx_);
        float* re1 = metal_m2l_get_src_mm_re1_buffer(ctx_);
        float* im1 = metal_m2l_get_src_mm_im1_buffer(ctx_);

#if defined(_OPENMP)
        #pragma omp parallel for schedule(static)
#endif
        for (int32_t b = 0; b < num_buckets_; ++b) {
            const auto* bucket = all_buckets[b];
            const auto& mm = bucket->MomentsMultipoleExpansion.MM_;
            const int32_t base_idx = b * num_coeffs_;

            for (int32_t c = 0; c < num_coeffs_; ++c) {
                int32_t idx = base_idx + c;
                const std::complex<double>& val0 = mm[c][0];
                const std::complex<double>& val1 = mm[c][1];

                re0[idx] = static_cast<float>(val0.real());
                im0[idx] = static_cast<float>(val0.imag());
                re1[idx] = static_cast<float>(val1.real());
                im1[idx] = static_cast<float>(val1.imag());
            }
        }
    }

    /*
     * Compute M2L transformation on GPU.
     */
    void compute() {
        int result = metal_m2l_compute_no_copy(ctx_);
        if (result != 0) {
            throw std::runtime_error("Metal M2L compute failed");
        }
    }

    /*
     * Submit async compute (returns immediately).
     */
    void computeAsync() {
        int result = metal_m2l_compute_async(ctx_);
        if (result != 0) {
            throw std::runtime_error("Metal M2L async compute failed");
        }
    }

    /*
     * Wait for async compute to complete.
     */
    void wait() {
        metal_m2l_wait(ctx_);
    }

    /*
     * Write results back to bucket MM arrays.
     * Results are accumulated (added to existing values).
     *
     * Note: OpenMP parallelization requires that different rows write to different
     * destination addresses. This is guaranteed by M2L structure (each row corresponds
     * to a unique (bucket, coeff) pair).
     */
    void writeResultsBack() {
        const float* results = metal_m2l_get_result_buffer(ctx_);

#if defined(_OPENMP)
        #pragma omp parallel for schedule(static)
#endif
        for (int32_t row = 0; row < num_rows_; ++row) {
            const RowInfo& info = row_info_[row];
            std::complex<double> delta0(results[row * 4 + 0], results[row * 4 + 1]);
            std::complex<double> delta1(results[row * 4 + 2], results[row * 4 + 3]);

            // Add to destination MM
            (*info.dst_mm)[0] += delta0;
            (*info.dst_mm)[1] += delta1;
        }
    }

    int32_t numBuckets() const { return num_buckets_; }
    int32_t numCoeffs() const { return num_coeffs_; }
    int32_t numRows() const { return num_rows_; }

private:
    MetalM2LContext* ctx_;
    int32_t num_buckets_;
    int32_t num_coeffs_;
    int32_t num_rows_;

    // Bucket pointer to index mapping
    std::unordered_map<const void*, int32_t> bucket_to_idx_;

    // Row information for result writeback
    struct RowInfo {
        std::array<std::complex<double>, 2>* dst_mm;
        int32_t bucket_idx;
    };
    std::vector<RowInfo> row_info_;
};

#endif /* METAL_M2L_WRAPPER_HPP */
