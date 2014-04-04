#define __host__
#define __device__

#include "gc_sum.h"

#include <assert.h>

#include <iostream>

using namespace std;

namespace unit_tests {

    void gc_sum() {
        auto validate = [] (uint64_t *vals, uint64_t n, gc_sum_t gc) {
            assert(gc.begin >= 0);
            assert(gc.end <= n);

            uint64_t sum = 0;
            for(uint64_t i = gc.begin; i < gc.end; i++) {
                sum += vals[i];
                uint64_t result = gc.get(i - gc.begin);
                assert(sum == result);
            }

            uint64_t rsum = 0;
            for(int i = (int)gc.end - 1; i >= (int)gc.begin; i--) {
                rsum += vals[i];
                int offset = i - gc.begin;
                uint64_t result = gc.get_reverse(offset);
                assert(rsum == result);
            }
        };

        uint64_t vals[] =   {   1, 2, 3,  4,  5,  6,  7,  8};
        uint64_t cumsum[] = {0, 1, 3, 6, 10, 15, 21, 28, 36};

        gc_sum_t gc;
        gc.cumsum = cumsum + 1;
        gc.begin = 0;
        gc.end = 8;
        gc.sum_begin = 0;
        gc.sum_end = 36;

        validate(vals, 8, gc);

        gc_sum_t a, b;
        gc.split(4, a, b);
        validate(vals, 8, a);
        validate(vals, 8, b);

        gc_sum_t aa, ab;
        a.split( 2, aa, ab );
        validate(vals, 8, aa);
        validate(vals, 8, ab);

        gc_sum_t ba, bb;
        b.split( 2, ba, bb );
        validate(vals, 8, ba);
        validate(vals, 8, bb);
    }

    void run() {
        gc_sum();
        cout << "Unit tests passed." << endl;
    }
}
