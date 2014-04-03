#define __host__
#define __device__

#include "isofinder.h"

#include "isofinder.cuh"
#include "queue.h"
#include "unit_tests.h"
#include "util.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <boost/math/special_functions/beta.hpp>

#include <pna.h>

using namespace seqio;
using namespace std;

#define square(x) ((x)*(x))

double absd(double x) {
    return x <= 0.0 ? -x : x;
}

base_frequency_t::base_frequency_t()
    : base_frequency_t(0, 0, 0) {
}

base_frequency_t::base_frequency_t(uint64_t at, uint64_t gc, uint64_t n) {
    frequency[AT] = at;
    frequency[GC] = gc;
    frequency[N] = n;
}

bool base_frequency_t::has_n() const {
    return frequency[N] != 0;
}

uint64_t base_frequency_t::get_at() const {
    return frequency[AT];
}

uint64_t base_frequency_t::get_gc() const {
    return frequency[GC];
}

base_frequency_t base_frequency_t::operator+(const base_frequency_t &other) const {
    base_frequency_t result = *this;
    result += other;
    return result;
}

base_frequency_t &base_frequency_t::operator+=(const base_frequency_t &other) {
    frequency[AT] += other.frequency[AT];
    frequency[GC] += other.frequency[GC];
    frequency[N] += other.frequency[N];

    return *this;
}

base_frequency_t &base_frequency_t::operator-=(const base_frequency_t &other) {
    frequency[AT] -= other.frequency[AT];
    frequency[GC] -= other.frequency[GC];
    frequency[N] -= other.frequency[N];

    return *this;
}

uint64_t base_frequency_t::operator[](category_t category) const {
    return frequency[category];
}

seq_t::seq_t(const string &origin_, uint64_t len_)
    : origin(origin_)
    , len(len_)
    , bases(new char[len_]) {
}

seq_t::~seq_t() {
    delete [] bases;
}

cut_t::cut_t() 
    : cut_t(0, 0, false) {
}

cut_t::cut_t(uint64_t begin_, uint64_t end_, bool is_n_)
    : begin(begin_)
    , end(end_)
    , is_n(is_n_)
    , is_incomplete(false) {
}

void cut_t::split(uint64_t midpoint, cut_t &left, cut_t &right) const {
    assert( (midpoint > 0) && (midpoint < length()) );

    left.begin = this->begin;
    left.end = this->begin + midpoint;
    left.is_n = this->is_n;
    left.is_incomplete = this->is_incomplete;

    right.begin = left.end;
    right.end = this->end;
    right.is_n = this->is_n;
    right.is_incomplete = this->is_incomplete;

    gc_sum.split(midpoint, left.gc_sum, right.gc_sum);
}

uint64_t cut_t::length() const {
    return end - begin;
}

//todo: port ran3()
double rng_t::next() {
    return drand48();
}

SignificanceMethod value_of(const std::string &name) {
    if(name == "r")
        return SignificanceMethod::Random;
    else if(name == "p1")
        return SignificanceMethod::Parametric1;
    else if(name == "p2")
        return SignificanceMethod::Parametric2;
    else if(name == "p3")
        return SignificanceMethod::Parametric3;
    else
        err("Invalid significance method %s", name.c_str());
}

std::string name_of(SignificanceMethod method) {
    switch(method) {
    case SignificanceMethod::Parametric1:
        return "t-student";
    case SignificanceMethod::Parametric2:
        return "t dif.var";
    case SignificanceMethod::Parametric3:
        return "Maximum";
    case SignificanceMethod::Random:
        return "Random";
    default:
        err("Unhandled method");
    }
}

shared_ptr<seq_t> read_seq(const string &path, uint64_t index) {
    PnaReader in(path.c_str());
    shared_ptr<PnaSequenceReader> sin = in.openSequence(index);
    shared_ptr<seq_t> seq = make_shared<seq_t>( path, sin->size() );
    sin->read(seq->bases, seq->len);
    return seq;
}

char random_base(rng_t &rng) {
    static char bases[] = {'A', 'T', 'C', 'G'};

    return bases[int(rng.next() * 4) % 4];
}

base_frequency_t get_base_frequency(char base) {
    static base_frequency_t at(1, 0, 0);
    static base_frequency_t gc(0, 1, 0);
    static base_frequency_t n(0, 0, 1);

    switch(base) {
    case 'A':
    case 'T':
        return at;
    case 'G':
    case 'C':
        return gc;
    case 'N':
        return n;
    default:
        err("%c", base);
    }
}

list<cut_t> find_N_islands(shared_ptr<seq_t> seq) {
    list<cut_t> cuts;
    cut_t cut;

    for(uint64_t i = 0; i < seq->len; i++) {
        base_frequency_t base_frequency = get_base_frequency( seq->bases[i] );

        if( (i == 0) || ((base_frequency.has_n()) != cut.is_n) ) {
            if(i != 0) {
                cut.end = i;
                cuts.push_back(cut);
            }
            cut.begin = i;
            cut.is_n = base_frequency.has_n();
        }
    }

    cut.end = seq->len;
    cuts.push_back(cut);

    return cuts;
}

template<typename Iterator>
Iterator previous(Iterator current) {
    return --current;
}

cut_t merge(const cut_t &a, const cut_t &b, bool is_n = false) {
    assert(a.end == b.begin);

    return cut_t(a.begin, b.end, is_n);
}

void merge_short_N_cuts(shared_ptr<seq_t> sequence, list<cut_t> &cuts) {
    rng_t rng;

    for(auto cut = cuts.begin(); cut != cuts.end(); ++cut) {
        if(cut->is_n && cut->length() < 3000) { // todo: inclusive? use MinimumCutSize?
            cut_t merged_cut = *cut;
            for(uint64_t i = merged_cut.begin; i < merged_cut.end; i++) {
                const char base = random_base(rng);
                sequence->bases[i] = base;
            }

            if(merged_cut.begin != 0) {
                merged_cut = merge( *previous(cut), merged_cut );
                cuts.erase( previous(cut) );
            }
            if(merged_cut.end != sequence->len) {
                merged_cut = merge( merged_cut, *next(cut) );
                cuts.erase( next(cut) );
            }

            *cut = merged_cut;
        }
    }
}

template<typename cuts_container>
void dump(cuts_container &cuts, bool onlyN) {
    for(auto &cut: cuts) {
        if(onlyN == cut.is_n) {
            printf("%9lu --> %9lu %9lu\n", cut.begin + 1, cut.end, cut.length());
        }
    }
}

void verify_bounds(shared_ptr<seq_t> sequence, list<cut_t> &cuts) {
    assert(cuts.front().begin == 0);
    assert(cuts.back().end == sequence->len);

    for(auto cut = cuts.begin(); next(cut) != cuts.end(); ++cut) {
        assert(cut->end == next(cut)->begin);
    }
}

list<cut_t> initialize_cuts(shared_ptr<seq_t> sequence) {
    list<cut_t> cuts = find_N_islands(sequence);

    printf("--------------\n" );
    printf("--- N Cuts ---\n" );
    printf("--------------\n" );
    dump(cuts, true);

    verify_bounds(sequence, cuts);

    printf("---------------------\n");
    printf("--- Merged N Cuts ---\n");
    printf("---------------------\n");
    merge_short_N_cuts(sequence, cuts);
    dump(cuts, true);

    verify_bounds(sequence, cuts);

    return cuts;
}

void mark_incomplete_cuts(shared_ptr<seq_t> sequence, vector<cut_t> &cuts) {
    for(auto cut = cuts.begin(); cut != cuts.end(); ++cut) {
        if(cut->is_n) {
            if(cut->begin != 0) {
                previous(cut)->is_incomplete = true;
            }
            if(cut->end != sequence->len) {
                next(cut)->is_incomplete = true;
            }
        }
    }

    cuts.front().is_incomplete = true;
    cuts.back().is_incomplete = true;
}

void dump_analysis_parameters(shared_ptr<seq_t> sequence,
                              const list<cut_t> &cuts,
                              double significance_level,
                              SignificanceMethod significance_method,
                              uint16_t window_length) {
    cout << "==================================================================" << endl;
    cout << "DNA file:                " << sequence->origin << endl;
    cout << "Length:                  " << sequence->len << endl;
    //cout << "G+C (%):                 " << float(base_frequency[base_frequency_t::GC]) / (sequence->len - n_frequency) * 100 << endl;
    //cout << "Undefined:               " << n_frequency << " (" << float(n_frequency) / sequence->len * 100 << "%)" << endl;
    cout << "Sig. Level:              " << significance_level << endl;
    cout << "Method:                  " << name_of(significance_method) << endl;
    cout << "Coarse-graining level:   " << window_length << endl;
    cout << "==================================================================" << endl;
}

struct split_candidate_t {
    double t = -1.0;
    double t_filter;
    struct partition_t {
        cut_t cut;
        struct window_parameters_t {
            uint64_t n;
            double mu;
            double sigma2;
        } window_parameters;
    } left, right;

    bool is_valid() {
        return t > 0.0;
    }
};

void create_gc_sum(shared_ptr<seq_t> sequence, cut_t &cut) {
    uint64_t *cumsum = new uint64_t[cut.length() + 1];
    cumsum[0] = 0;
    cumsum++;

    uint64_t sum = 0;
    for(uint64_t i = 0; i < cut.length(); i++) {
        sum += get_base_frequency(sequence->bases[cut.begin+i]).get_gc();
        cumsum[i] = sum;
    }

    gc_sum_t &result = cut.gc_sum;
    result.cumsum = cumsum;
    result.begin = 0;
    result.end = cut.length();
    result.sum_begin = 0;
    result.sum_end = sum;

    result.gpu_alloc();
}

void dispose_gc_sum(cut_t cut) {
    delete [] (cut.gc_sum.cumsum - 1);
    cut.gc_sum.gpu_dispose();
}

void __tstudent_split(gc_sum_t gc_sum,
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
        double t = absd(mu1 - mu2) / sd;
        if( t > best_t ) {
            best_t = t;
            best_midpoint = midpoint;
        }
    }
}

split_candidate_t find_maximal_tstudent_split(shared_ptr<seq_t> sequence,
                                              const cut_t cut,
                                              uint64_t min_size,
                                              uint16_t window_length) {

    split_candidate_t result;

    {
        uint64_t best_midpoint;
        double best_t = -1;

        cuda_tstudent_split(cut.gc_sum, min_size, best_midpoint, best_t);

        result.t = best_t;
        if(best_t > -1) {
            cut.split( best_midpoint, result.left.cut, result.right.cut );
        }
    }

    if(result.is_valid()) {
        auto compute_window_parameters =
            [sequence, window_length] (split_candidate_t::partition_t &partition) {

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
        };

        compute_window_parameters(result.left);
        compute_window_parameters(result.right);

        split_candidate_t::partition_t::window_parameters_t &left_parms = result.left.window_parameters;
        split_candidate_t::partition_t::window_parameters_t &right_parms = result.right.window_parameters;

        double sd = sqrt( (left_parms.sigma2 + right_parms.sigma2) *
                          (1.0/left_parms.n+1.0/right_parms.n) /
                          (left_parms.n+right_parms.n-2) );
        result.t_filter = absd( (left_parms.mu - right_parms.mu) / sd );

    }

    return result;
}

double ttest0(double t, double df) {
    return boost::math::ibeta(0.5*df, 0.5, df/(df+t*t));
}

double ttest_difvar(double ave1, double ave2, double var1, double var2, uint64_t n1, uint64_t n2) {
    double t, df, prob;

    t = (ave1-ave2)/sqrt(var1/n1+var2/n2);
    df = square(var1/n1+var2/n2)/(square(var1/n1)/(n1-1)+square(var2/n2)/(n2-1));
    prob = boost::math::ibeta(0.5*df,0.5,df/(df+square(t)));

    return prob;
}

double neff(double n) {
    return -11.54+4.189*log(n);
}

double prob(double x, double n) {

    if(abs(x)<1.0e-5) {
        return 0.0;
    } else {
        const double expon = neff(n);
        const double f = 0.8;
        const double nu = n - 2;
        return pow( (1.0-boost::math::ibeta(0.5*f*nu,0.5*f,nu/(nu+square(x)))), expon );
    }
}

double random_ttest(cuda_random::TTestState &ttest,
                    cuda_random::GpuState &gpu,
                    shared_ptr<seq_t> sequence,
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

    auto init_gc_windows = [](double *gc_windows, ulong window_count, char *bases, ulong window_length) {
        for(ulong i = 0; i < window_count; i++) {
            double gc = 0;
            for(ushort j = 0; j < window_length; j++) {
                base_frequency_t base_frequency = get_base_frequency(bases[i * window_length + j]);
                gc += base_frequency.get_gc();
            }
            gc /= window_length;
            gc_windows[i] = gc;
        }
    };

    init_gc_windows( gc_windows,
                     left_window_count,
                     sequence->bases + left_cut.begin,
                     window_length );

    init_gc_windows( gc_windows + left_window_count,
                     right_window_count,
                     sequence->bases + right_cut.begin,
                     window_length );

    gpu.init(&ttest, window_count, 50000);
    ulong npass = gpu.exec(left_window_count, right_window_count, gc_windows, t_filter, 50000);
    gpu.dispose();

    ttest.dispose(gc_windows);

    return double(npass) / Experiments_Count;
}

class find_split_context_t {
    function<bool(cut_t, pair<cut_t, cut_t> &)> impl;
public:
    find_split_context_t(function<bool(cut_t, pair<cut_t, cut_t> &)> impl_) : impl(impl_) {}

    bool find_split(cut_t cut, pair<cut_t, cut_t> &result) {
        return impl(cut, result);
    }
};

vector<cut_t> __find_isochores(list<cut_t> initial_cuts,
                               vector<find_split_context_t> &split_contexts) {
    vector<cut_t> processed_cuts;
    queue_t<cut_t> working_cuts;

    for(auto cut: initial_cuts) {
        if(cut.is_n) {
            processed_cuts.push_back(cut);
        } else {
            working_cuts.push(cut);
        }
    }

    atomic_uint working_cuts_remaining(working_cuts.size());
    mutex processed_mutex;

    thread threads[split_contexts.size()];
    for(size_t i = 0; i < split_contexts.size(); i++) {
        find_split_context_t &split_context = split_contexts[i];

        threads[i] = thread([&split_context, &processed_cuts, &working_cuts, &working_cuts_remaining, &processed_mutex] {
                cut_t cut;

                while(working_cuts.pop(cut)) {
                    pair<cut_t, cut_t> split;
                    if( split_context.find_split(cut, split) ) {
                        working_cuts.push(split.first);
                        working_cuts.push(split.second);
                        working_cuts_remaining++;
                    } else {
                        {
                            unique_lock<mutex> lock(processed_mutex);
                            processed_cuts.push_back(cut);
                        }
                        if(0 == --working_cuts_remaining) {
                            working_cuts.close();
                        }
                    }
                }
            });
    }

    for(auto &thread: threads) {
        thread.join();
    }

    sort(processed_cuts.begin(), processed_cuts.end(),
         [](const cut_t &a, const cut_t &b) {
             return a.begin < b.begin;
         });

    return processed_cuts;
}

vector<cut_t> find_isochores(shared_ptr<seq_t> sequence,
                             double significance_level,
                             SignificanceMethod significance_method,
                             uint16_t window_length) {
    list<cut_t> initial_cuts = initialize_cuts(sequence);
    dump_analysis_parameters(sequence, initial_cuts, significance_level, significance_method, window_length);

    for(auto &cut: initial_cuts) {
        create_gc_sum(sequence, cut);
    }

    vector<find_split_context_t> contexts;
    const int nthreads = 1;
    for(int i = 0; i < nthreads; i++) {
        auto find_split = [i, sequence, significance_level, significance_method, window_length]
            (cut_t cut, pair<cut_t, cut_t> &result) {

            //cuda_random::TTestState ttest;
            //cuda_random::GpuState gpu(i);

            const uint64_t  Minimum_Cut_Size = 10 * window_length + 2 + 2; // todo: fortran effectively uses 2 more as min

            split_candidate_t candidate = find_maximal_tstudent_split(sequence,
                                                                      cut,
                                                                      Minimum_Cut_Size,
                                                                      window_length);
            if(!candidate.is_valid()) {
                return false;
            }

/*
            cout << cut.begin << "  " << cut.end << "  " << candidate.left.cut.end << "   "
            << candidate.left.window_parameters.sigma2 << "  " << candidate.right.window_parameters.sigma2 << "  "
            << candidate.t_filter << "   " << candidate.left.window_parameters.n << "  " << candidate.right.window_parameters.n;
*/
            double p;

            switch(significance_method) {
            case SignificanceMethod::Parametric1: {
                double df = candidate.left.window_parameters.n + candidate.right.window_parameters.n - 2;
                p = 1.0 - ttest0(candidate.t_filter, df);
            } break;
            case SignificanceMethod::Parametric2: {
                uint64_t n1 = candidate.left.window_parameters.n;
                uint64_t n2 = candidate.right.window_parameters.n;
                double mu1 = candidate.left.window_parameters.mu;
                double mu2 = candidate.right.window_parameters.mu;
                double nv1 = candidate.left.window_parameters.sigma2;
                double nv2 = candidate.right.window_parameters.sigma2;

                p = 1.0 - ttest_difvar(mu1, mu2, nv1/(n1-1), nv2/(n2-1), n1, n2);
            } break;
            case SignificanceMethod::Parametric3: {
                p = prob(candidate.t_filter, candidate.left.window_parameters.n + candidate.right.window_parameters.n);
            } break;
            case SignificanceMethod::Random: {
/*
                p = random_ttest(ttest,
                                 gpu,
                                 sequence,
                                 candidate.left.cut,
                                 candidate.right.cut,
                                 window_length,
                                 candidate.t_filter);
*/
                assert(false);
            } break;
            default:
                assert(false);
            }

            //cout << "  " << p << endl;

            if(p >= significance_level) {
                result.first = candidate.left.cut;
                result.second = candidate.right.cut;
                return true;
            } else {
                return false;
            }
        };
        contexts.push_back(find_split_context_t(find_split));
    }

    vector<cut_t> final_cuts = __find_isochores(initial_cuts, contexts);

    mark_incomplete_cuts(sequence, final_cuts);

/*
    printf("--------------\n");
    printf("--- Result ---\n");
    printf("--------------\n");
    dump(final_cuts, false);
*/

    for(auto &cut: initial_cuts) {
        dispose_gc_sum(cut);
    }

    return final_cuts;
}

void write_results(vector<cut_t> &cuts, const string &path) {
    ofstream out(path.c_str());

    for(auto cut: cuts) {
        if(!cut.is_n && !cut.is_incomplete) {
            out << cut.begin << '\t' << cut.end << '\t' << cut.length() << endl;
        }
    }
}

int main(int argc, const char **argv) {
    if( (argc > 1) && (0 == strcmp(argv[1], "--ut")) ) {
        unit_tests::run();
        return 0;
    }

    // todo: arg error handling
    int argi = 1;
    const char *sequence_path = argv[argi++];
    uint64_t sequence_index = atoi(argv[argi++]);
    double significance = atof(argv[argi++]);
    SignificanceMethod significance_method = value_of(argv[argi++]);
    uint16_t window_length = atoi(argv[argi++]);
    const char *output_path = argv[argi++];
    

    vector<cut_t> cuts = find_isochores( read_seq(sequence_path, sequence_index),
                                         significance,
                                         significance_method,
                                         window_length);

    write_results(cuts, output_path);

    return 0;
}
