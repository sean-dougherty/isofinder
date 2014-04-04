#pragma once

#include "gc_sum.h"

void tstudent_split_cuda(gc_sum_t gc_sum,
                         uint64_t min_size,
                         uint64_t &best_midpoint,
                         double &best_t);
