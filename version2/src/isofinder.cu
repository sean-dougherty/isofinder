#include <assert.h>
#include <math.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <iostream>

#include <cuda.h>

#include "isofinder.cuh"

using namespace std;

#define xcuda(stmt) {							\
	cudaError_t err = stmt;						\
	if (err != cudaSuccess) {					\
	    cerr << __FILE__ << ":" << __LINE__ << ": Failed to run " << #stmt << ". Reason: " << cudaGetErrorString(err) << endl; \
	    exit(1);							\
	}								\
    }

namespace cuda_random {

    struct rng_t {
    private:
	uint state;

    public:
	__host__ __device__ rng_t(uint seed) : state(seed) {}

	__host__ __device__ double next() {
	    uint const a = 16807; //ie 7**5
	    uint const m = UINT32_MAX; //ie 2**32

	    state = max(1, state);
	    state = (ulong(state * a))%m;

	    return double(state) / UINT32_MAX;
	}
    };

    __global__ void experiment(ulong left_window_count,
			       ulong right_window_count,
			       double *gc_windows,
			       ulong experiments_count,
			       double t_filter,
			       uint *result) {
	uint tid = (blockDim.x * blockIdx.x) + threadIdx.x;
	if(tid >= experiments_count)
	    return;
	__shared__ uint t_pass;
	if(threadIdx.x == 0)
	    t_pass = 0;

	__syncthreads();

	struct cut_state_t {
	    ulong window_count;
	    ulong windows_remaining;
	    double x;
	    double x2;

	    __device__ cut_state_t(ulong window_count_) {
		window_count = window_count_;
		windows_remaining = window_count;
		x = 0.0f;
		x2 = 0.0f;
	    }

	    __device__ void add_gc(double gc) {
		windows_remaining--;
		x += gc;
		x2 += gc * gc;
	    }

	    __device__ void complete() {
		x /= window_count;
		x2 -= window_count * (x * x);
	    }
	};

	const double degrees_freedom = (left_window_count + right_window_count) - 2;
	rng_t rng(tid + 1);
	cut_state_t left(left_window_count), right(right_window_count);

	for(ulong i = 0; i < left_window_count + right_window_count; i++) {
	    double p_left = double(left.windows_remaining) / (left.windows_remaining + right.windows_remaining);
	    if(rng.next() < p_left) {
		left.add_gc(gc_windows[i]);
	    } else {
		right.add_gc(gc_windows[i]);
	    }
	}

	left.complete();
	right.complete();

	double t = sqrt((left.x2+right.x2)*(1.0/left.window_count+1.0/right.window_count)/degrees_freedom);
	t = abs(left.x - right.x) / t;

	if(t <= t_filter) {
	    atomicAdd(&t_pass, 1);
	}

	__syncthreads();

	if(threadIdx.x == 0) {
	    result[blockIdx.x] = t_pass;
	}
    }

    TTestState::TTestState() : host_gc_windows(NULL), alloced_window_count(0) {
    }

    double *TTestState::alloc(ulong window_count) {
	return host_gc_windows = new double[window_count];
/*
	if(window_count > alloced_window_count) {
	    if(host_gc_windows) {
		xcuda( cudaFreeHost(host_gc_windows) );
	    }
	    xcuda( cudaMallocHost((void**)&host_gc_windows, window_count * sizeof(double), 0) );
	    alloced_window_count = window_count;
	}
	return host_gc_windows;
*/
    }

    void TTestState::dispose(double *gc_windows) {
	delete [] gc_windows;
/*
	if(alloced_window_count) {
	    xcuda( cudaFreeHost(host_gc_windows) );
	    alloced_window_count = 0;
	}
*/
    }

    GpuState::GpuState(int dev_) : dev(dev_) {
	alloced_window_count = 0;
	alloced_nblocks = 0;
    }

    static const ulong Threads_Per_Block = 512;
    
    void GpuState::ctxt() {
	cudaSetDevice(dev);
    }

    void GpuState::init(TTestState *state, ulong window_count, ulong experiments_count) {
	uint nblocks = (experiments_count - 1) / Threads_Per_Block + 1;

	xcuda( cudaMallocHost((void**)&host_results, nblocks * sizeof(uint), 0) );
	xcuda( cudaMalloc((void **)&device_gc_windows, sizeof(*device_gc_windows) * window_count) );
	xcuda( cudaMalloc((void **)&device_results, sizeof(*device_results) * nblocks) );
	xcuda( cudaMemcpy(device_gc_windows, state->host_gc_windows, sizeof(*device_gc_windows) * window_count, cudaMemcpyHostToDevice) );
    }

    void GpuState::dispose() {
	xcuda( cudaFreeHost(host_results) );
	xcuda( cudaFree(device_results) );	
	xcuda( cudaFree(device_gc_windows) );	
    }

    ulong GpuState::exec(ulong left_window_count,
			 ulong right_window_count,
			 double *host_gc_windows,
			 double t_filter,
			 ulong experiments_count) {

	uint nblocks = (experiments_count - 1) / Threads_Per_Block + 1;

	experiment<<<nblocks, Threads_Per_Block>>>(left_window_count,
						   right_window_count,
						   device_gc_windows,
						   experiments_count,
						   t_filter,
						   device_results);

	xcuda( cudaPeekAtLastError() );
	xcuda( cudaThreadSynchronize() );
	xcuda( cudaMemcpy(host_results, device_results, sizeof(*device_results) * nblocks, cudaMemcpyDeviceToHost) );
	
	ulong npass = 0;
	for(ulong i = 0; i < nblocks; i++) {
	    npass += host_results[i];
	}

	return npass;
    }

}
