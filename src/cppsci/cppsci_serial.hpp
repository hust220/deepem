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

#include "cppsci_traits.hpp"
#include "cppsci_array.hpp"

namespace cppsci {

template<typename T>
inline void serialize(STD_ ostream &stream, const T *v) {
    stream.write(reinterpret_cast<const char *>(&v), sizeof(T*));
}

template<typename T>
inline void parse(STD_ istream &stream, T *v) {
    stream.read(reinterpret_cast<char *>(&v), sizeof(T*));
}

inline void serialize(STD_ ostream &stream, const int &v) {
    stream.write(reinterpret_cast<const char *>(&v), sizeof v);
}

inline void parse(STD_ istream &stream, int &v) {
    stream.read(reinterpret_cast<char *>(&v), sizeof v);
}

inline void serialize(STD_ ostream &stream, const float &v) {
    stream.write(reinterpret_cast<const char *>(&v), sizeof v);
}

inline void parse(STD_ istream &stream, float &v) {
    stream.read(reinterpret_cast<char *>(&v), sizeof v);
}

inline void serialize(STD_ ostream &stream, const double &v) {
    stream.write(reinterpret_cast<const char *>(&v), sizeof v);
}

inline void parse(STD_ istream &stream, double &v) {
    stream.read(reinterpret_cast<char *>(&v), sizeof v);
}

inline void serialize(STD_ ostream &stream, const char &v) {
    stream.write(reinterpret_cast<const char *>(&v), sizeof v);
}

inline void parse(STD_ istream &stream, char &v) {
    stream.read(reinterpret_cast<char *>(&v), sizeof v);
}

inline void serialize(STD_ ostream &stream, const Str &v) {
    stream.write(reinterpret_cast<const char *>(v.c_str()), v.size()+1);
}

inline void parse(STD_ istream &stream, Str &v) {
    char c[2];
    STD_ ostringstream ostream;

    while (stream.read(c, 1)) {
        if (c[0] == '\0') break;
        ostream << c[0];
    }
    v = ostream.str();
}

template<typename T>
inline void serialize(STD_ ostream &stream, const V<T> &v) {
    int n = v.size();
    serialize(stream, n);

    for (auto && i : v) {
        serialize(stream, i);
    }
}

template<typename T>
inline void parse(STD_ istream &stream, V<T> &v) {
    int n;
    parse(stream, n);

    v.resize(n);
    for (auto && i : v) {
        parse(stream, i);
    }
}

template<typename T>
inline void serialize(STD_ ostream &stream, const L<T> &v) {
    int n = v.size();
    serialize(stream, n);

    for (auto && i : v) {
        serialize(stream, i);
    }
}

template<typename T>
inline void parse(STD_ istream &stream, L<T> &v) {
    int n;
    parse(stream, n);

    v.resize(n);
    for (auto && i : v) {
        parse(stream, i);
    }
}

template<typename T>
inline void serialize(STD_ ostream &stream, const Q<T> &v) {
    int n = v.size();
    serialize(stream, n);

    for (auto && i : v) {
        serialize(stream, i);
    }
}

template<typename T>
inline void parse(STD_ istream &stream, Q<T> &v) {
    int n;
    parse(stream, n);

    v.resize(n);
    for (auto && i : v) {
        parse(stream, i);
    }
}

template<typename T, int N>
inline void serialize(STD_ ostream &stream, const A<T, N> &v) {
    for (auto && i : v) {
        serialize(stream, i);
    }
}

template<typename T, int N>
inline void parse(STD_ istream &stream, A<T, N> &v) {
    for (auto && i : v) {
        parse(stream, i);
    }
}

template<typename T>
inline void serialize(STD_ ostream &stream, const Array<T> &v) {
    Shape shape = v.shape();
    serialize(stream, shape);
    stream.write(reinterpret_cast<const char *>(v.data()), sizeof(T) * v.size());
}

template<typename T>
inline void parse(STD_ istream &stream, Array<T> &v) {
    Shape shape;
    parse(stream, shape);
    v = Array<T>(shape);
    stream.read(reinterpret_cast<char *>(v.data()), sizeof(T) * v.size());
}

} // namespace cppsci

