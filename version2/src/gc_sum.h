#pragma once

#include <stdint.h>

class gc_sum_t {
public:
    void split(uint64_t midpoint, gc_sum_t &left, gc_sum_t &right) const;
    void gpu_alloc();
    void gpu_dispose();
    gc_sum_t gpu_copy();

    __device__ __host__ uint64_t get(uint64_t offset) const;
    __device__ __host__ uint64_t get_reverse(uint64_t offset) const;
    __device__ __host__ uint64_t range(uint64_t begin, uint64_t end) const;

    __device__ __host__ uint64_t length() const;

    uint64_t *cumsum;
    uint64_t *gpu_cumsum;
    uint64_t begin;
    uint64_t end;
    uint64_t sum_begin;
    uint64_t sum_end;
};

void cuda_tstudent_split(gc_sum_t gc_sum,
                         uint64_t min_size,
                         uint64_t &best_midpoint,
                         double &best_t);
