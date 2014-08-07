#define __host__
#define __device__

#include "cut.h"
#include "queue.h"
#include "split_candidate.h"
#include "split_significance.cuh"
#include "split_significance.h"
#include "unit_tests.h"
#include "util.h"

#include <assert.h>
#include <stdio.h>

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <seqio.h>


using namespace std;

void dump_analysis_parameters(seq_t &sequence,
                              const list<cut_t> &cuts,
                              double significance_level,
                              split_significance_method_t significance_method,
                              uint16_t window_length) {
    cout << "==================================================================" << endl;
    cout << "DNA file:                " << sequence.origin << endl;
    cout << "Length:                  " << sequence.len << endl;
    cout << "Sig. Level:              " << significance_level << endl;
    cout << "Method:                  " << name_of(significance_method) << endl;
    cout << "Coarse-graining level:   " << window_length << endl;
    cout << "==================================================================" << endl;
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

vector<cut_t> find_isochores(seq_t &sequence,
                             double significance_level,
                             split_significance_method_t significance_method,
                             uint16_t window_length) {
    list<cut_t> initial_cuts = initialize_cuts(sequence);
    dump_analysis_parameters(sequence, initial_cuts, significance_level, significance_method, window_length);

    for(auto &cut: initial_cuts) {
        create_gc_sum(&sequence, &cut);
    }

    vector<find_split_context_t> contexts;
    const int nthreads = 1;
    for(int i = 0; i < nthreads; i++) {
        auto find_split = [i, sequence, significance_level, significance_method, window_length]
            (cut_t cut, pair<cut_t, cut_t> &result) {

            split_candidate_t candidate;
            if(!find_best_split_candidate(cut,
                                          window_length,
                                          candidate)) {
                return false;
            }

            if( is_split_significant(candidate, significance_method, significance_level, window_length) ) {
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

    for(auto &cut: initial_cuts) {
        dispose_gc_sum(&cut);
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

    out.close();
}

int main(int argc, const char **argv) {
    int argi = 1;
    const char *sequence_path = argv[argi++];
    double significance = atof(argv[argi++]);
    split_significance_method_t significance_method = get_split_significance_method(argv[argi++]);
    uint16_t window_length = atoi(argv[argi++]);
    const char *output_path = argv[argi++];

    vector<cut_t> cuts;
    {
        char *sequence_buffer = nullptr;
        uint64_t sequence_buffer_length;
        uint64_t sequence_length;

        seqio_sequence_iterator iterator;
        seqio_sequence_options sequence_options = SEQIO_DEFAULT_SEQUENCE_OPTIONS;
        sequence_options.base_transform = SEQIO_BASE_TRANSFORM_CAPS_GATCN;
        seqio_create_sequence_iterator(sequence_path,
                                       sequence_options,
                                       &iterator);
        
        seqio_sequence sequence;
        while( (SEQIO_SUCCESS == seqio_next_sequence(iterator, &sequence))
               && (sequence != nullptr) ) {

            char const *name, *comment;
            seqio_const_dictionary metadata;
            seqio_get_metadata(sequence, &metadata);
            seqio_get_value(metadata, SEQIO_KEY_NAME, &name);
            seqio_get_value(metadata, SEQIO_KEY_COMMENT, &comment);

            seqio_read_all(sequence, &sequence_buffer, &sequence_buffer_length, &sequence_length);

            seq_t seq(sequence_path, name, sequence_buffer, sequence_length);
            vector<cut_t> sequence_cuts
                = find_isochores(seq,
                                 significance,
                                 significance_method,
                                 window_length);

            cuts.insert(cuts.end(), sequence_cuts.begin(), sequence_cuts.end());

            seqio_dispose_sequence(&sequence);
        }

        seqio_dispose_sequence_iterator(&iterator);
        seqio_dispose_buffer(&sequence_buffer);
    }

    write_results(cuts, output_path);

    return 0;
}
