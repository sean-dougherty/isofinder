#include "gc_sum.h"

#include "sequence.h"
#include "cut.h"

#include <assert.h>
#include <cuda.h>
#include <iostream>

#define xcuda(stmt) {                                                   \
        cudaError_t err = stmt;                                         \
        if (err != cudaSuccess) {                                       \
            std::cerr << __FILE__ << ":" << __LINE__ << ": Failed to run " << #stmt << ". Reason: " << cudaGetErrorString(err) << std::endl; \
            exit(1);                                                    \
        }                                                               \
    }

void create_gc_sum(seq_t *sequence, cut_t *cut) {
    uint64_t *cumsum = new uint64_t[cut->length() + 1];
    cumsum[0] = 0;
    cumsum++;

    uint64_t sum = 0;
    for(uint64_t i = 0; i < cut->length(); i++) {
        if( sequence->is_gc(cut->begin+i) )
            sum++;
        cumsum[i] = sum;
    }

    gc_sum_t &result = cut->gc_sum;
    result.cumsum = cumsum;
    result.begin = 0;
    result.end = cut->length();
    result.sum_begin = 0;
    result.sum_end = sum;

    result.gpu_alloc();
}

void dispose_gc_sum(cut_t *cut) {
    delete [] (cut->gc_sum.cumsum - 1);
    cut->gc_sum.gpu_dispose();
}

void gc_sum_t::gpu_alloc() {
    size_t sizeof_ = sizeof(uint64_t) * (end - begin + 1);
    xcuda( cudaMalloc((void**)&gpu_cumsum, sizeof_) );
    xcuda( cudaMemcpy(gpu_cumsum, cumsum - 1, sizeof_, cudaMemcpyHostToDevice) );
}

void gc_sum_t::gpu_dispose() {
    xcuda( cudaFree(gpu_cumsum) );
}

gc_sum_t gc_sum_t::gpu_copy() {
    gc_sum_t result = *this;
    result.cumsum = gpu_cumsum + 1;
    return result;
}

void gc_sum_t::split(uint64_t midpoint, gc_sum_t &left, gc_sum_t &right) const {
    left.cumsum = this->cumsum;
    left.gpu_cumsum = this->gpu_cumsum;
    left.begin = this->begin;
    left.end = this->begin + midpoint;
    left.sum_begin = this->sum_begin;
    left.sum_end = left.cumsum[left.end - 1];

    right.cumsum = this->cumsum;
    right.gpu_cumsum = this->gpu_cumsum;
    right.begin = left.end;
    right.end = this->end;
    right.sum_begin = left.sum_begin + left.get(midpoint - 1);
    right.sum_end = right.cumsum[right.end - 1];
}

uint64_t gc_sum_t::range(uint64_t begin, uint64_t end) const {
    return cumsum[this->begin + end - 1] - cumsum[this->begin + begin - 1];
}
