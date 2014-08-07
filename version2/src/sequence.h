#pragma once

#include <stdint.h>
#include <string>

class seq_t {
public:
    seq_t(const std::string &origin_,
          const std::string &name_,
          char *bases_,
          uint64_t len_);

    bool is_gc(uint64_t index);
    bool is_n(uint64_t index);
    void set_random_base(uint64_t index);

    const std::string origin;
    const std::string name;
    uint64_t len;

private:
    char *bases;
};
