#define __device__
#define __host__

#include "split_candidate.h"
#include "split_candidate.cuh"

#include "gc_sum.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

void tstudent_split_cpu(gc_sum_t gc_sum,
                        uint64_t min_size,
                        uint64_t &best_midpoint,
                        double &best_t) {
    best_midpoint = 0;
    best_t = -1;

    for(uint64_t midpoint = min_size; midpoint <= gc_sum.length() - min_size; midpoint++) {
        double n1 = midpoint;
        double n2 = gc_sum.length() - midpoint;
        double mu1 = gc_sum.get(midpoint - 1) / n1;
        double mu2 = gc_sum.get_reverse(midpoint) / n2;
        double sd = (mu1*(1-mu1)+mu2*(1-mu2))*(1/n1+1/n2)/(n1+n2-2);
        if(sd <= 0) {
            continue;
        }
        sd = sqrt(sd);
        double t = fabsl(mu1 - mu2) / sd;
        if( t > best_t ) {
            best_t = t;
            best_midpoint = midpoint;
        }
    }
}

static void compute_window_parameters(uint16_t window_length, split_candidate_t::partition_t &partition) {
    uint64_t window_count = partition.cut.length() / window_length; // todo: ignoring final bases
    assert(window_count > 0);

    double sum = 0.0;
    double sum2 = 0.0;

    for(size_t i = 0; i < window_count; i++) {
        uint16_t window_sum = (uint16_t)partition.cut.gc_sum.range(i*window_length, (i+1)*window_length);
        sum += window_sum;
        sum2 += window_sum * window_sum;
    }

    sum /= window_length;
    sum2 /= window_length * window_length;

    double mu = sum / window_count;

    partition.window_parameters.n = window_count;
    partition.window_parameters.mu = mu;
    partition.window_parameters.sigma2 = sum2 - window_count * (mu * mu);
}
                                      

bool find_best_split_candidate(const cut_t cut,
                               uint16_t window_length,
                               split_candidate_t &result) {

    const uint64_t  Minimum_Cut_Size = 10 * window_length + 2 + 2; // todo: fortran effectively uses 2 more as min

    {
        uint64_t best_midpoint;
        double best_t = -1;

        tstudent_split_cuda(cut.gc_sum, Minimum_Cut_Size, best_midpoint, best_t);

        if(best_t < 0.0) {
            return false;
        }

        cut.split( best_midpoint, result.left.cut, result.right.cut );
    }

    compute_window_parameters(window_length, result.left);
    compute_window_parameters(window_length, result.right);

    split_candidate_t::partition_t::window_parameters_t &left_parms = result.left.window_parameters;
    split_candidate_t::partition_t::window_parameters_t &right_parms = result.right.window_parameters;

    double sd = sqrt( (left_parms.sigma2 + right_parms.sigma2) *
                      (1.0/left_parms.n+1.0/right_parms.n) /
                      (left_parms.n+right_parms.n-2) );
    result.t_filter = fabsl( (left_parms.mu - right_parms.mu) / sd );

    return true;
}
