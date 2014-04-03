#pragma once

#include <cuda.h>

namespace cuda_random {
    void init(int dev);

    double *alloc_double(ulong count);
    void free_double(double *ptr);

    class TTestState {
    public:
	double *host_gc_windows;
	ulong alloced_window_count;

	TTestState();

	double *alloc(ulong window_count);
	void dispose(double *gc_windows);
    };

    class GpuState {
    public:
	uint *host_results;
	double *device_gc_windows;
	uint *device_results;
	int dev;

	ulong alloced_window_count;
	ulong alloced_nblocks;

	GpuState(int dev);

	void ctxt();
	void init(TTestState *state, ulong window_count, ulong experiments_count);
	void dispose();
	ulong exec(ulong left_window_count,
		   ulong right_window_count,
		   double *host_gc_windows,
		   double t_filter,
		   ulong experiments_count);
    };
}
