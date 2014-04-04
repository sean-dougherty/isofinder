#define __device__
#define __host__

#include "split_significance.h"

#include "split_significance.cuh"
#include "sequence.h"
#include "util.h"

#include <boost/math/special_functions/beta.hpp>

using std::shared_ptr;
using std::string;

split_significance_method_t get_split_significance_method(const string &name) {
    if(name == "r")
        return split_significance_method_t::Random;
    else if(name == "p1")
        return split_significance_method_t::Parametric1;
    else if(name == "p2")
        return split_significance_method_t::Parametric2;
    else if(name == "p3")
        return split_significance_method_t::Parametric3;
    else
        err("Invalid significance method %s", name.c_str());
}

string name_of(split_significance_method_t method) {
    switch(method) {
    case split_significance_method_t::Parametric1:
        return "t-student";
    case split_significance_method_t::Parametric2:
        return "t dif.var";
    case split_significance_method_t::Parametric3:
        return "Maximum";
    case split_significance_method_t::Random:
        return "Random";
    default:
        err("Unhandled method");
    }
}

static double ttest0(double t, double df) {
    return boost::math::ibeta(0.5*df, 0.5, df/(df+t*t));
}

static double ttest_difvar(double ave1, double ave2, double var1, double var2, uint64_t n1, uint64_t n2) {
    double t, df, prob;

    t = (ave1-ave2)/sqrt(var1/n1+var2/n2);
    df = square(var1/n1+var2/n2)/(square(var1/n1)/(n1-1)+square(var2/n2)/(n2-1));
    prob = boost::math::ibeta(0.5*df,0.5,df/(df+square(t)));

    return prob;
}

static double neff(double n) {
    return -11.54+4.189*log(n);
}

static double prob(double x, double n) {

    if(fabsl(x)<1.0e-5) {
        return 0.0;
    } else {
        const double expon = neff(n);
        const double f = 0.8;
        const double nu = n - 2;
        return pow( (1.0-boost::math::ibeta(0.5*f*nu,0.5*f,nu/(nu+square(x)))), expon );
    }
}

static double ttest_random(cuda_random::TTestState &ttest,
                           cuda_random::GpuState &gpu,
                           cut_t left_cut,
                           cut_t right_cut,
                           uint16_t window_length,
                           double t_filter) {
    gpu.ctxt();

    ulong left_window_count = left_cut.length() / window_length;
    ulong right_window_count = right_cut.length() / window_length;
    ulong window_count = left_window_count + right_window_count;
    double *gc_windows = ttest.alloc(window_count);

    const ulong Experiments_Count = 50000;

    auto init_gc_windows = [](double *gc_windows, ulong window_count, ulong window_length, cut_t cut) {
        for(ulong i = 0; i < window_count; i++) {
            double gc = cut.gc_sum.range(i*window_length, (i+1)*window_length);
            gc /= window_length;
            gc_windows[i] = gc;
        }
    };

    init_gc_windows( gc_windows,
                     left_window_count,
                     window_length,
                     left_cut);

    init_gc_windows( gc_windows + left_window_count,
                     right_window_count,
                     window_length,
                     right_cut);

    gpu.init(&ttest, window_count, 50000);
    ulong npass = gpu.exec(left_window_count, right_window_count, gc_windows, t_filter, 50000);
    gpu.dispose();

    ttest.dispose(gc_windows);

    return double(npass) / Experiments_Count;
}

bool is_split_significant(const split_candidate_t &candidate,
                          split_significance_method_t method,
                          double significance_level,
                          uint16_t window_length) {
    double p;

    switch(method) {
    case split_significance_method_t::Parametric1: {
        double df = candidate.left.window_parameters.n + candidate.right.window_parameters.n - 2;
        p = 1.0 - ttest0(candidate.t_filter, df);
    } break;
    case split_significance_method_t::Parametric2: {
        uint64_t n1 = candidate.left.window_parameters.n;
        uint64_t n2 = candidate.right.window_parameters.n;
        double mu1 = candidate.left.window_parameters.mu;
        double mu2 = candidate.right.window_parameters.mu;
        double nv1 = candidate.left.window_parameters.sigma2;
        double nv2 = candidate.right.window_parameters.sigma2;

        p = 1.0 - ttest_difvar(mu1, mu2, nv1/(n1-1), nv2/(n2-1), n1, n2);
    } break;
    case split_significance_method_t::Parametric3: {
        p = prob(candidate.t_filter, candidate.left.window_parameters.n + candidate.right.window_parameters.n);
    } break;
    case split_significance_method_t::Random: {
        cuda_random::TTestState ttest;
        cuda_random::GpuState gpu(0);

        p = ttest_random(ttest,
                         gpu,
                         candidate.left.cut,
                         candidate.right.cut,
                         window_length,
                         candidate.t_filter);
    } break;
    default:
        assert(false);
    }
    
    return p >= significance_level;
}
