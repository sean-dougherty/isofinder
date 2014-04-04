#include "sequence.h"

#include <pna.h>

using namespace seqio;
using std::make_shared;
using std::shared_ptr;
using std::string;

seq_t seq_t::read(const string &path, uint64_t index) {
    PnaReader in(path.c_str());
    shared_ptr<PnaSequenceReader> sin = in.openSequence(index);
    seq_t seq( path, sin->size() );
    sin->read(seq.bases, seq.len);
    return seq;
}

seq_t::seq_t(const string &origin_, uint64_t len_)
    : origin(origin_)
    , len(len_)
    , bases(new char[len_]) {
}

void seq_t::dispose() {
    delete [] bases;
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
