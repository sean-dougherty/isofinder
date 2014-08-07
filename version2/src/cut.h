#pragma once

#include "gc_sum.h"
#include "sequence.h"

#include <list>
#include <stdint.h>

class cut_t {
public:
    cut_t();
    cut_t(uint64_t begin_, uint64_t end_, bool is_n_);
    
    void split(uint64_t midpoint, cut_t &left, cut_t &right) const;

    uint64_t length() const;

    gc_sum_t gc_sum;
    uint64_t begin; // inclusive
    uint64_t end;   // exclusive
    bool is_n;
    bool is_incomplete; // is it immediately before or after N island?
};

std::list<cut_t> initialize_cuts(seq_t &sequence);
void mark_incomplete_cuts(seq_t &sequence, std::list<cut_t> &cuts);
void merge_short_N_cuts(seq_t &sequence, std::list<cut_t> &cuts);
void verify_bounds(seq_t &sequence, std::list<cut_t> &cuts);
