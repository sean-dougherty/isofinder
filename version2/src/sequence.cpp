#include "sequence.h"

using std::string;

seq_t::seq_t(const string &origin_,
             const string &name_,
             char *bases_,
             uint64_t len_)
    : origin(origin_)
    , name(name_)
    , len(len_)
    , bases(bases_) {
}

bool seq_t::is_gc(uint64_t index) {
    char base = bases[index];

    return (base == 'G') || (base == 'C');
}

bool seq_t::is_n(uint64_t index) {
    char base = bases[index];

    return (base == 'N');
}

void seq_t::set_random_base(uint64_t index) {
    static char Bases[] = {'A', 'T', 'C', 'G'};

    bases[index] = Bases[int(drand48() * 4) % 4];
}
