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

#include "cppsci_math.hpp"
#include "cppsci_traits.hpp"

namespace cppsci {

template<typename _Val, typename Derived_> class BasicArray;
template<typename _Val> class Array;
template<typename _Val> class MapArray;
template<typename _Val> class SubArray;

#define JN_IS_COMPLEX(_type) jn_is_complex<typename jn_array_val<_type>::type>::value

template<typename T>
struct jn_is_complex {
private:
    using F = typename STD_ decay<T>::type;
    template<typename _Val> static STD_ true_type check(STD_ complex<_Val>);
    static STD_ false_type check(...);
public:
    enum { value = STD_ is_same<decltype(check(STD_ declval<F>())), STD_ true_type>::value };
};

template<typename T>
struct jn_complex_val {
private:
    using F = typename STD_ decay<T>::type;
    template<typename _Val> static _Val check(STD_ complex<_Val>);
    static T check(...);
public:
    using type = decltype(check(STD_ declval<F>()));
};
template<typename T>
using jn_complex_val_t = typename jn_complex_val<T>::type;

template<typename T>
struct jn_array_val {
private:
    using F = typename STD_ decay<T>::type;
    template<typename _Val, typename Derived_> static _Val check(BasicArray<_Val, Derived_>);
                template<typename _Val                   > static _Val check(Array<_Val>);
                template<typename _Val                   > static _Val check(SubArray<_Val>);
                template<typename _Val                   > static _Val check(MapArray<_Val>);
    static T check(...);
public:
    using type = decltype(check(STD_ declval<F>()));
};
template<typename T>
using jn_array_val_t = typename jn_array_val<T>::type;

/**
 * Judge whether a type is an array.
 */
template<typename T>
struct jn_is_array {
private:
    using F = typename ::std::decay<T>::type;
    template<typename _Val, typename Derived_> static ::std::true_type check(const BasicArray<_Val, Derived_> *);
//    template<typename _Val> static ::std::true_type check(Array<_Val>);
//    template<typename _Val> static ::std::true_type check(MapArray<_Val>);
//    template<typename _Val> static ::std::true_type check(SubArray<_Val>);
    static ::std::false_type check(...);
public:
    enum { value = ::std::is_same<decltype(check(::std::declval<F *>())), ::std::true_type>::value };
};

#define JN_IS_ARRAY(_type) jn_is_array<_type>::value

#define ARRAY_EACH(shape, ind) \
    for (Shape ind = array_helper::ind_begin(shape); array_helper::ind_in(ind, shape); array_helper::ind_next(ind, shape))

#define ARRAY_EACH_R(shape, ind) \
    for (Shape ind = array_helper::ind_end(shape); array_helper::ind_in(ind, shape); array_helper::ind_prev(ind, shape))

#define ARRAY3_EACH(shape, i, a, b, c) \
    for (int i = 0, a = 0; a < shape[0]; a++) \
    for (int b = 0; b < shape[1]; b++) \
    for (int c = 0; c < shape[2]; c++, i++)

#define ARRAY3_EACH_R(shape, i, a, b, c) \
    for (int i = shape[0] * shape[1] * shape[2] - 1, a = shape[0] - 1; a >= 0; a--) \
    for (int b = shape[1]-1; b >= 0; b--) \
    for (int c = shape[2]-1; c >= 0; c--, i--)

} // namespace cppsci


