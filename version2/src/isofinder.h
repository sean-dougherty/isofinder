#pragma once

#include <stddef.h>
#include <stdint.h>

#include <memory>
#include <string>

#include "gc_sum.h"

class base_frequency_t {
public:
    enum category_t {
	AT = 0, GC = 1, N = 2
    };

    base_frequency_t();
    base_frequency_t(uint64_t at, uint64_t gc, uint64_t n);

    bool has_n() const;
    uint64_t get_at() const;
    uint64_t get_gc() const;

    base_frequency_t operator+(const base_frequency_t &other) const;
    base_frequency_t &operator+=(const base_frequency_t &other);
    base_frequency_t &operator-=(const base_frequency_t &other);
    uint64_t operator[](category_t category) const;

private:
    uint64_t frequency[3];
};

class seq_t { // todo: rename to sequence_t
public:
    seq_t(const std::string &origin_, uint64_t len_);
    ~seq_t();

    const std::string origin;
    uint64_t len; // todo: rename to length
    char *bases;
};

enum class SignificanceMethod {
    Parametric1,
    Parametric2,
    Parametric3,
    Random
};

SignificanceMethod value_of(const std::string &name);
std::string name_of(SignificanceMethod method);

class cut_t {
public:
    cut_t();
    cut_t(uint64_t begin_, uint64_t end_, bool is_n_);
    
    void split(uint64_t midpoint, cut_t &left, cut_t &right) const;

    uint64_t length() const;

    gc_sum_t gc_sum;
    uint64_t begin; // inclusive
    uint64_t end;   // exclusive
    bool is_n;
    bool is_incomplete; // is it immediately before or after N island?
};

class rng_t {
public:
    double next();
};
