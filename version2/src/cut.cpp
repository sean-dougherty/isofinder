#define __device__
#define __host__

#include "cut.h"

#include <assert.h>

using std::list;

template<typename cuts_container>
void dump(cuts_container &cuts, bool onlyN) {
    for(auto &cut: cuts) {
        if(onlyN == cut.is_n) {
            printf("%9lu --> %9lu %9lu\n", cut.begin + 1, cut.end, cut.length());
        }
    }
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

template<typename Iterator>
Iterator previous(Iterator current) {
    return --current;
}

static cut_t merge(const cut_t &a, const cut_t &b, bool is_n = false) {
    assert(a.end == b.begin);

    return cut_t(a.begin, b.end, is_n);
}

void merge_short_N_cuts(seq_t &sequence, list<cut_t> &cuts) {
    for(auto cut = cuts.begin(); cut != cuts.end(); ++cut) {
        if(cut->is_n && cut->length() < 3000) { // todo: inclusive? use MinimumCutSize?
            cut_t merged_cut = *cut;
            for(uint64_t i = merged_cut.begin; i < merged_cut.end; i++) {
                sequence.set_random_base(i);
            }

            if(merged_cut.begin != 0) {
                merged_cut = merge( *previous(cut), merged_cut );
                cuts.erase( previous(cut) );
            }
            if(merged_cut.end != sequence.len) {
                merged_cut = merge( merged_cut, *next(cut) );
                cuts.erase( next(cut) );
            }

            *cut = merged_cut;
        }
    }
}

void verify_bounds(seq_t &sequence, list<cut_t> &cuts) {
    assert(cuts.front().begin == 0);
    assert(cuts.back().end == sequence.len);

    for(auto cut = cuts.begin(); next(cut) != cuts.end(); ++cut) {
        assert(cut->end == next(cut)->begin);
    }
}

list<cut_t> find_N_islands(seq_t &seq) {
    list<cut_t> cuts;
    cut_t cut;

    for(uint64_t i = 0; i < seq.len; i++) {
        bool is_n = seq.is_n(i);

        if( (i == 0) || (is_n != cut.is_n) ) {
            if(i != 0) {
                cut.end = i;
                cuts.push_back(cut);
            }
            cut.begin = i;
            cut.is_n = is_n;
        }
    }

    cut.end = seq.len;
    cuts.push_back(cut);

    return cuts;
}

list<cut_t> initialize_cuts(seq_t &sequence) {
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

void mark_incomplete_cuts(seq_t &sequence, list<cut_t> &cuts) {
    for(auto cut = cuts.begin(); cut != cuts.end(); ++cut) {
        if(cut->is_n) {
            if(cut->begin != 0) {
                previous(cut)->is_incomplete = true;
            }
            if(cut->end != sequence.len) {
                next(cut)->is_incomplete = true;
            }
        }
    }

    cuts.front().is_incomplete = true;
    cuts.back().is_incomplete = true;
}
