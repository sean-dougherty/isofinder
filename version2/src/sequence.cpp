#include "sequence.h"

#include <pna.h>

using namespace seqio;
using std::make_shared;
using std::shared_ptr;
using std::string;

seq_t::seq_t(const string &origin_, uint64_t len_)
    : origin(origin_)
    , len(len_)
    , bases(new char[len_]) {
}

seq_t::~seq_t() {
    delete [] bases;
}

shared_ptr<seq_t> read_seq(const string &path, uint64_t index) {
    PnaReader in(path.c_str());
    shared_ptr<PnaSequenceReader> sin = in.openSequence(index);
    shared_ptr<seq_t> seq = make_shared<seq_t>( path, sin->size() );
    sin->read(seq->bases, seq->len);
    return seq;
}

char random_base() {
    static char bases[] = {'A', 'T', 'C', 'G'};

    return bases[int(drand48() * 4) % 4];
}
