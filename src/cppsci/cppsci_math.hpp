/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(jianopt@163.com)"
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#pragma once

#include "cppsci_time.hpp"

namespace cppsci {

/**
 * Square
 */
double square(double n);

/**
 * Return a random floating number between 0 and 1.
 */
double rand();

/*
 * Set the seed of the random engine
 */
void seed(unsigned t);

struct Rand {
    Rand() {}

    Rand(int seed_) : seed(seed_) {}

    double operator ()() {
        double r;
        assert(seed >= 0);
        seed = ((long long)a * seed + c) % modules;
        r = double(seed)/modules;
        return r;
    }

    double operator ()(int s) {
        assert(s >= 0);
        s = ((long long)a * s + c) % modules;
        return double(s)/modules;
    }

protected:
    int seed = 1;
    int modules = 2147483647;
    int a = 1103515245;
    int c = 12345;

};

} // namespace jncpp

