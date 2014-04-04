#pragma once

#include <memory>
#include <stdint.h>
#include <string>

class seq_t { // todo: rename to sequence_t
public:
    seq_t(const std::string &origin_, uint64_t len_);
    ~seq_t();

    const std::string origin;
    uint64_t len; // todo: rename to length
    char *bases;
};

std::shared_ptr<seq_t> read_seq(const std::string &path, uint64_t index);
char random_base();
