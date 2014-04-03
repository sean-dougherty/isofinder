#include "gc_sum.h"

#include <assert.h>
#include <cuda.h>
#include <iostream>

using namespace std;

static const uint Threads_Per_Block = 512;


#define xcuda(stmt) {                                                   \
        cudaError_t err = stmt;                                         \
        if (err != cudaSuccess) {                                       \
            cerr << __FILE__ << ":" << __LINE__ << ": Failed to run " << #stmt << ". Reason: " << cudaGetErrorString(err) << endl; \
            exit(1);                                                    \
        }                                                               \
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

uint64_t gc_sum_t::get(uint64_t offset) const {
    return cumsum[begin + offset] - sum_begin;
}

uint64_t gc_sum_t::get_reverse(uint64_t offset) const {
    return sum_end - cumsum[begin + offset - 1];
}

uint64_t gc_sum_t::length() const {
    return end - begin;
}

__global__ void tstudent(gc_sum_t gc_sum,
                         uint64_t min_size,
                         uint64_t n,
                         double *ts,
                         uint64_t *is) {
    uint i = (blockDim.x * blockIdx.x) + threadIdx.x;
    int tid = threadIdx.x;
    __shared__ double ts_shared[Threads_Per_Block];
    __shared__ int is_shared[Threads_Per_Block];

    ts_shared[tid] = -1.0;
    is_shared[tid] = tid;

    if(i >= n) {
        return;
    }

    uint64_t midpoint = min_size + i;
    double n1 = midpoint;
    double n2 = gc_sum.length() - midpoint;
    double mu1 = gc_sum.get(midpoint - 1) / n1;
    double mu2 = gc_sum.get_reverse(midpoint) / n2;
    double sd = (mu1*(1-mu1)+mu2*(1-mu2))*(1/n1+1/n2)/(n1+n2-2);
    double t;
    if(sd <= 0) {
        t = -1;
    } else {
        sd = sqrt(sd);
        t = abs(mu1 - mu2) / sd;
    }

    ts_shared[tid] = t;

    for(int stride = Threads_Per_Block / 2; stride > 0; stride >>= 1) {
        __syncthreads();

        if(tid < stride) {
            if(ts_shared[tid + stride] > ts_shared[tid]) {
                ts_shared[tid] = ts_shared[tid + stride];
                is_shared[tid] = is_shared[tid + stride];
            }
        }
    }

    if(tid == 0) {
        ts[blockIdx.x] = ts_shared[0];
        is[blockIdx.x] = blockDim.x * blockIdx.x + is_shared[0];
    }

}

void cuda_tstudent_split(gc_sum_t gc_sum,
                         uint64_t min_size,
                         uint64_t &best_midpoint,
                         double &best_t) {
    best_midpoint = 0;
    best_t = -1.0;

    if( gc_sum.length() < 2*min_size )
        return;

    uint64_t n = (gc_sum.length() - min_size) - (min_size) + 1;

    const uint nblocks = (n - 1) / Threads_Per_Block + 1;

    double *ts = new double[nblocks];
    double *gpu_ts;
    uint64_t *is = new uint64_t[nblocks];
    uint64_t *gpu_is;

    xcuda( cudaMalloc((void**)&gpu_ts, nblocks*sizeof(double)) );
    xcuda( cudaMalloc((void**)&gpu_is, nblocks*sizeof(uint64_t)) );
    
    tstudent<<<nblocks, Threads_Per_Block>>>(gc_sum.gpu_copy(),
                                             min_size,
                                             n,
                                             gpu_ts,
                                             gpu_is);
	xcuda( cudaPeekAtLastError() );
	xcuda( cudaThreadSynchronize() );

    xcuda( cudaMemcpy(ts, gpu_ts, nblocks * sizeof(double), cudaMemcpyDeviceToHost) );    
    xcuda( cudaMemcpy(is, gpu_is, nblocks * sizeof(uint64_t), cudaMemcpyDeviceToHost) );    

    for(uint64_t i = 0; i < nblocks; i++) {
        double t = ts[i];
        if(t > best_t) {
            best_t = t;
            best_midpoint = min_size + is[i];
        }
    }

    delete [] ts;
    delete [] is;
    xcuda( cudaFree(gpu_ts) );
    xcuda( cudaFree(gpu_is) );

}
