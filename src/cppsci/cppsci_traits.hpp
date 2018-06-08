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

#include <type_traits>
#include <functional>
#include <vector>
#include <list>
#include <deque>
#include <map>
#include <utility>
#include <memory>
#include <set>
#include <tuple>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <chrono>
#include <cassert>
#include <regex>
#include <iterator>
#include <cfloat>
#include <sstream>
#include <regex>
#include <cctype>
#include <cmath>
#include <random>
#include <climits>
#include <array>
#include <complex>
#include <iomanip>
#include <numeric>
#include <ios>

#if defined(CPPSCI_OS_WIN)
#include <windows.h>
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include <thread>
#include "omp.h"

#include "cppsci_macros.hpp"
#include "cppsci_platform.hpp"

#define JN_PP_CAT(a, b)  a##b
#define JN_PP_CAT2(a, b) JN_PP_CAT(a, b)
#define JN_PP_CAT3(a, b) JN_PP_CAT2(a, b)

#define JN_PP_STR(a)  #a
#define JN_PP_STR2(a) JN_PP_STR(a)
#define JN_PP_STR3(a) JN_PP_STR2(a)

namespace cppsci {

#define JN_ENABLE(_cond) typename ::std::enable_if<_cond, int>::type = 1

template<bool Condition, typename Type>
using jn_enable_t = typename ::std::enable_if<Condition, Type>::type;

template<bool A, bool B> struct jn_static_or                 { enum { value = true  }; };
template<>               struct jn_static_or<false, false>   { enum { value = false }; };
#define JN_STATIC_OR(A, B)   jn_static_or<A, B>::value

template<bool A, bool B> struct jn_static_and                { enum { value = false }; };
template<>               struct jn_static_and<true, true>    { enum { value = true  }; }; 
#define JN_STATIC_AND(A, B)  jn_static_and<A, B>::value

template<bool A>         struct jn_static_not                { enum { value = false }; };
template<>               struct jn_static_not<false>         { enum { value = true  }; }; 
#define JN_STATIC_NOT(A)     jn_static_not<A>::value

/**
 * Calculate the number of types.
 */
template<typename _T, typename... _Types>
struct jn_num_types {
    enum { N = 1 + jn_num_types<_Types...>::N };
};

template<typename _T>
struct jn_num_types<_T> { enum { N = 1 }; };

/**
 * Retrieve the first type.
 */
template<typename _T, typename... _Types>
struct jn_first {
    using type = _T;
};

template<typename... _Types>
using jn_first_t = typename jn_first<_Types...>::type;

#ifdef JN_PRECISION
using Num = JN_PRECISION;
#else
using Num = double;
#endif

/**
 * Definition of jn_if.
 */
template<bool condition, typename T, typename F> struct jn_if {
    using type = T;
};

template<typename T, typename F>
struct jn_if<false, T, F> {
    using type = F;
};

template<bool condition, typename T, typename F>
using jn_if_t = typename jn_if<condition, T, F>::type;

template<typename _CharType, typename _CharTraits = STD_ char_traits<_CharType>>
using BStr = STD_ basic_string<_CharType, _CharTraits>;
using Str = STD_ string;
using WStr = STD_ wstring;

template<typename _Type>
using Ptr = _Type *;

template<typename _Type>
using SPtr = STD_ shared_ptr<_Type>;
template<typename _Type>
using SP = SPtr<_Type>;

template<typename _Type>
using UPtr = STD_ unique_ptr<_Type>;
template<typename _Type>
using UP = UPtr<_Type>;

template<typename _Fty>
using Fn = STD_ function<_Fty>;
template<typename _Fty>
using Function = STD_ function<_Fty>;

template<typename _Type, int _N>
using A = STD_ array<_Type, _N>;
template<int _N> using An = A<Num, _N>;
template<int _N> using Ai = A<int, _N>;
template<int _N> using Ad = A<double, _N>;
template<int _N> using Af = A<float, _N>;
template<int _N> using Ac = A<char, _N>;
template<int _N> using As = A<Str, _N>;
template<int _N> using Ab = A<bool, _N>;

template<typename _Type>
using V = STD_ vector<_Type>;
using Vn = V<Num>;
using Vi = V<int>;
using Vd = V<double>;
using Vf = V<float>;
using Vc = V<char>;
using Vs = V<Str>;
using Vb = V<bool>;

template<typename _Type>
using L = STD_ list<_Type>;
using Ln = L<Num>;
using Li = L<int>;
using Ld = L<double>;
using Lf = L<float>;
using Lc = L<char>;
using Ls = L<Str>;
using Lb = L<bool>;

template<typename _Type>
using Q = STD_ deque<_Type>;
using Qn = Q<Num>;
using Qi = Q<int>;
using Qd = Q<double>;
using Qf = Q<float>;
using Qc = Q<char>;
using Qs = Q<Str>;
using Qb = Q<bool>;

template<typename _Type>
using T = STD_ set<_Type>;
using Tn = T<Num>;
using Ti = T<int>;
using Td = T<double>;
using Tf = T<float>;
using Tc = T<char>;
using Ts = T<Str>;
using Tb = T<bool>;

template<typename _Type1, typename _Type2>
using Pr = STD_ pair<_Type1, _Type2>;

template<typename... _Types>
using Tp = STD_ tuple<_Types...>;

template<typename _KeyType, typename _Type>
using M = STD_ map<_KeyType, _Type>;
template<typename _Type>
using Mn = M<Num, _Type>;
template<typename _Type>
using Mi = M<int, _Type>;
template<typename _Type>
using Md = M<double, _Type>;
template<typename _Type>
using Mf = M<float, _Type>;
template<typename _Type>
using Mc = M<char, _Type>;
template<typename _Type>
using Ms = M<Str, _Type>;
template<typename _Type>
using Mb = M<bool, _Type>;

template<typename T, typename U>
T lexical_cast(U && u) {
    ::std::stringstream stream;
    T t;

    stream << u;
    stream >> t;
    return t;
}

#define CPPSCI_INT(a) ::cppsci::lexical_cast<int>(a)
#define CPPSCI_DBL(a) ::cppsci::lexical_cast<double>(a)
#define CPPSCI_FLT(a) ::cppsci::lexical_cast<float>(a)
#define CPPSCI_STR(a) ::cppsci::lexical_cast<Str>(a)

template<typename Stream_>
void stream_push(Stream_ && stream) {
}

template<typename Stream_, typename Value_, typename... Values_>
void stream_push(Stream_ && stream, Value_ && value, Values_ && ...values) {
    stream << value;
    stream_push(stream, values...);
}

template<typename List1, typename List2>
bool eles_equal(const List1 &l1, const List2 &l2) {
    auto it1 = l1.begin();
    auto it2 = l2.begin();
    for (; it1 != l1.end() && it2 != l2.end(); it1++, it2++) {
        if (*it1 != *it2) return false;
    }
    if (it1 != l1.end() || it2 != l2.end()) return false;
    return true;
}

} // namespace cppsci
