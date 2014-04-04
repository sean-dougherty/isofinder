#pragma once

#include <stdint.h>
#include <string>

class seq_t { // todo: rename to sequence_t
public:
    static seq_t read(const std::string &path, uint64_t index);

    seq_t(const std::string &origin_, uint64_t len_);

    void dispose();

    bool is_gc(uint64_t index);
    bool is_n(uint64_t index);
    void set_random_base(uint64_t index);

    const std::string origin;
    uint64_t len; // todo: rename to length

private:
    char *bases;
};
