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

#include "cppsci_array_core.hpp"
#include "cppsci_array_traits.hpp"

namespace cppsci {

/**
 * Create an indices array.
 */
template<typename type>
inline type indices(double a, double b, double step = 1, int n = 0) {
    if (step != 0 && n == 0) n = int(::std::ceil((b-a)/step));
    else if (step == 0 && n != 0) step = (b-a)/(n-1);
    else throw "Array::range error!";

    type m(n);
    for (int i = 0; i < n; i++) m[i] = (typename type::value_type)(a + i * step);
    return ::std::move(m);
}

/**
 * Create an indices array.
 */
template<typename type>
inline type indices(int a) {
    type m({1});
    m[0] = a;
    return ::std::move(m);
}

/**
 * Create an indices array.
 */
template<typename type, typename T>
inline type indices_v(const T &t) {
    int n = STD_ distance(STD_ begin(t), STD_ end(t));
    if (n == 0) return type{};
    else if (n == 1) return indices<type>(t[0]);
    else if (n == 2) return indices<type>(t[0], t[1]);
    else if (n == 3) return indices<type>(t[0], t[1], t[2]);
    else return indices<type>(t[0], t[1], t[2], t[3]);
}

/**
 * Sub array
 */
template<typename NumType_>
class SubArray : public BasicArray<NumType_, SubArray<NumType_>> {
public:
    using value_type = NumType_;
    using type = SubArray<value_type>;
    using self_type = SubArray<value_type>;
    using base_type = BasicArray<value_type, type>;

    STD_ vector<value_type *> m_data;

    SubArray() {
        base_type::set_shape({0});
    }

    template<typename _A, JN_ENABLE(JN_IS_ARRAY(_A))>
    self_type &operator =(const _A &a) {
        for (int i = 0; i < base_type::size(); i++) at(i) = a(i);
        return *this;
    }

    template<typename _V, JN_ENABLE(JN_STATIC_NOT(JN_IS_ARRAY(_V)))>
    self_type &operator =(const _V &v) {
        for (int i = 0; i < base_type::size(); i++) at(i) = v;
        return *this;
    }

    virtual value_type &at(int a) { return *(m_data[a]); }
    virtual const value_type &at(int a) const { return *(m_data[a]); }

};

template<typename T>
struct ArrayValue_M {
    using type = typename STD_ remove_reference<decltype(STD_ declval<T>()[0])>::type;
};
template<typename T>
using ArrayValue_T = typename ArrayValue_M<T>::type;

template<typename T>
struct SubArray_M {
    using value_type = ArrayValue_T<T>;
    using type = SubArray<jn_if_t<STD_ is_const<T>::value, const value_type, value_type>>;
};
template<typename T>
using SubArray_T = typename SubArray_M<T>::type;

template<typename _Array, typename _Ls>
SubArray_T<_Array> array_subi(_Array &&array, _Ls &&range) {
    SubArray_T<_Array> sub;
    sub.set_shape(range.shape());
    sub.m_data.resize(sub.size());
    for (int i = 0; i < sub.size(); i++) sub.m_data[i] = &array(range(i));
    return STD_ move(sub);
}

template<typename _Array, typename _Ls>
SubArray_T<_Array> array_subv(_Array &&array, _Ls &&v) {
    SubArray_T<_Array> sub;

    int nd = array.dim();
    if (v.size() < nd) for (int i = 0; i < nd-v.size(); i++) v.push_back(Arrayi{});

    Shape shape(nd);
    for (int i = 0; i < nd; i++) {
        if (v[i].size() == 0) shape[i] = array.shape(i);
        else shape[i] = v[i].size();
    }
    sub.set_shape(shape);

    sub.m_data.resize(sub.size());

    auto in_range = [&nd, &v](const Shape &ind) {
        for (int i = 0; i < nd; i++) {
            if (v[i].begin() != v[i].end() && STD_ find(v[i].begin(), v[i].end(), ind[i]) == v[i].end()) return false;
        }
        return true;
    };

    int m = 0;
    int n = 0;
    ARRAY_EACH(array.shape(), ind) {
        if (in_range(ind)) {
            sub.m_data[n] = &(array[m]);
            n++;
        }
        m++;
    }

    return STD_ move(sub);
}

template<typename _Array>
SubArray_T<_Array> array_sub(_Array &&array, STD_ initializer_list<Vd> range) {
    STD_ vector<Arrayi> v;
    for (auto && r : range) v.push_back(indices_v<Arrayi>(r));
    return array_subv(array, v);

}

template<typename _Array, typename _Ls>
SubArray_T<_Array> array_sub(_Array &&array, _Ls &&range) {
    STD_ vector<Arrayi> v;
    for (auto && r : range) v.push_back(indices_v<Arrayi>(r));
    return array_subv(array, v);

}

template<typename _Array>
SubArray_T<_Array> array_take(_Array && array, int val, int axis = 0, int mode = ARRAY_TAKE_MODE_CLIP) {
    ERROR_CHECK(axis < 0 || axis >= array.dim(), "ILLEGAL axis!!!");

    SubArray_T<_Array> sub;
    int nd = array.dim();

    // Set shape
    Shape shape;
    for (int i = 0; i < nd; i++) if (i != axis) shape.push_back(array.shape(i));
    sub.set_shape(shape);

    // Set m_data
    sub.m_data.resize(array.size()/array.shape(axis));
    int i = 0;
    int j = 0;
    ARRAY_EACH(array.shape(), ind) {
        if (ind[axis] == val) {
            sub.m_data[j] = &(array[i]);
            j++;
        }
        i++;
    }

    return STD_ move(sub);
}

template<typename _Array, typename _LS>
SubArray_T<_Array> array_takes(_Array && array, _LS && ls, int axis = 0, int mode = ARRAY_TAKE_MODE_CLIP) {
    ERROR_CHECK(axis < 0 || axis >= array.dim(), "ILLEGAL axis!!!");

    SubArray_T<_Array> sub;
    int nd = array.dim();

    // Set shape
    Shape shape = array.shape();
    shape[axis] = ls.size();
    sub.set_shape(shape);

    // Set m_data
    sub.m_data.resize(array.size()/array.shape(axis)*ls.size());
    int i = 0;
    ARRAY_EACH(array.shape(), ind) {
        if (ind[axis] == array.shape(axis) - 1) {
            int d = 0;
            for (auto && n : ls) {
                if (n >= array.shape(axis)-1) {
                    if (n > array.shape(axis)-1 && mode != ARRAY_TAKE_MODE_CLIP) ERROR_REPORT("Indices error!");
                    Shape v = ind;
                    v[axis] = d;
                    int sum = 0; for (int k = 0; k < nd; k++) sum += sub.coeff(k) * v[k];
                    sub.m_data[sum] = &(array[i]);
                }
                d++;
            }
        }
        else if (ind[axis] == 0) {
            int d = 0;
            for (auto && n : ls) {
                if (n <= 0) {
                    if (n < 0 && mode != ARRAY_TAKE_MODE_CLIP) ERROR_REPORT("Indices error!");
                    Shape v = ind;
                    v[axis] = d;
                    int sum = 0; for (int k = 0; k < nd; k++) sum += sub.coeff(k) * v[k];
                    sub.m_data[sum] = &(array[i]);
                }
                d++;
            }
        }
        else {
            auto it = STD_ find(ls.begin(), ls.end(), ind[axis]);
            if (it != ls.end()) {
                int d = STD_ distance(ls.begin(), it);
                Shape v = ind;
                v[axis] = d;
                int sum = 0; for (int k = 0; k < nd; k++) sum += sub.coeff(k) * v[k];
                sub.m_data[sum] = &(array[i]);
            }
        }
        i++;
    }

    return STD_ move(sub);
}


} // namespace cppsci

