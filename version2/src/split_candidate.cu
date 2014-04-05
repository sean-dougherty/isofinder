#include "split_candidate.cuh"

#include "gc_sum.h"

#include <cuda.h>
#include <iostream>

#define xcuda(stmt) {                                                   \
        cudaError_t err = stmt;                                         \
        if (err != cudaSuccess) {                                       \
            std::cerr << __FILE__ << ":" << __LINE__ << ": Failed to run " << #stmt << ". Reason: " << cudaGetErrorString(err) << std::endl; \
            exit(1);                                                    \
        }                                                               \
    }

static const uint Threads_Per_Block = 512;

__global__ void tstudent(gc_sum_t gc_sum,
                         uint64_t min_size,
                         uint64_t n,
                         double *ts,
                         uint64_t *is) {
    double best_t = -1.0;
    uint64_t best_i;

    for(uint64_t i = (blockDim.x * blockIdx.x) + threadIdx.x;
        i < n;
        i += gridDim.x * blockDim.x) {

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

            if(t > best_t) {
                best_t = t;
                best_i = i;
            }
        }
    }

    __shared__ double ts_shared[Threads_Per_Block];
    __shared__ uint64_t is_shared[Threads_Per_Block];

    int tid = threadIdx.x;
    ts_shared[tid] = best_t;
    is_shared[tid] = best_i;

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
        is[blockIdx.x] = is_shared[0];
    }

}

void tstudent_split_cuda(gc_sum_t gc_sum,
                         uint64_t min_size,
                         uint64_t &best_midpoint,
                         double &best_t) {
    best_midpoint = 0;
    best_t = -1.0;

    uint64_t n = (gc_sum.length() - min_size) - (min_size) + 1;

    uint64_t nblocks = (n - 1) / Threads_Per_Block + 1;
    if(nblocks > 2048)
        nblocks = 2048;

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
