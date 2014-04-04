#pragma once

#include "split_candidate.h"

#include <string>
#include <stdint.h>

enum class split_significance_method_t {
    Parametric1,
    Parametric2,
    Parametric3,
    Random
};

split_significance_method_t get_split_significance_method(const std::string &name);
std::string name_of(split_significance_method_t method);

bool is_split_significant(const split_candidate_t &candidate,
                          split_significance_method_t method,
                          double significance_level,
                          uint16_t window_length);


