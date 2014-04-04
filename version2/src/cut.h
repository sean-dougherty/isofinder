#pragma once

#include "gc_sum.h"
#include "sequence.h"

#include <list>
#include <memory>
#include <stdint.h>
#include <vector>

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

std::list<cut_t> initialize_cuts(std::shared_ptr<seq_t> sequence);
void mark_incomplete_cuts(std::shared_ptr<seq_t> sequence, std::vector<cut_t> &cuts);
void merge_short_N_cuts(std::shared_ptr<seq_t> sequence, std::list<cut_t> &cuts);
void verify_bounds(std::shared_ptr<seq_t> sequence, std::list<cut_t> &cuts);

template<typename cuts_container>
void dump(cuts_container &cuts, bool onlyN) {
    for(auto &cut: cuts) {
        if(onlyN == cut.is_n) {
            printf("%9lu --> %9lu %9lu\n", cut.begin + 1, cut.end, cut.length());
        }
    }
}
