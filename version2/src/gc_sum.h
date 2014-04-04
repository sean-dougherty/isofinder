#pragma once

#include <stdint.h>

class gc_sum_t {
public:
    void split(uint64_t midpoint, gc_sum_t &left, gc_sum_t &right) const;
    void gpu_alloc();
    void gpu_dispose();
    gc_sum_t gpu_copy();
    uint64_t range(uint64_t begin, uint64_t end) const;

    __device__ __host__ uint64_t get(uint64_t offset) const {
        return cumsum[begin + offset] - sum_begin;
    }
    __device__ __host__ uint64_t get_reverse(uint64_t offset) const {
        return sum_end - cumsum[begin + offset - 1];
    }
    __device__ __host__ uint64_t length() const {
        return end - begin;
    }

    uint64_t *cumsum;
    uint64_t *gpu_cumsum;
    uint64_t begin;
    uint64_t end;
    uint64_t sum_begin;
    uint64_t sum_end;
};

void create_gc_sum(class seq_t *sequence, class cut_t *cut);
void dispose_gc_sum(class cut_t *cut);
