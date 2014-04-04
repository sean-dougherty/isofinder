#pragma once

#include "cut.h"

#include <stdint.h>

struct split_candidate_t {
    double t_filter;

    struct partition_t {
        cut_t cut;
        struct window_parameters_t {
            uint64_t n;
            double mu;
            double sigma2;
        } window_parameters;
    } left, right;
};

bool find_best_split_candidate(const cut_t cut,
                               uint16_t window_length,
                               split_candidate_t &result);
