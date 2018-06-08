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

namespace cppsci {

#define JN_DEF_BINARY_OP(T, V, W) \
    static W check(T, V); \
    static W check(V, T);

    template<typename T, typename V>
    struct jn_simple_op_type {
    private:
        using F = typename STD_ decay<T>::type;
        using G = typename STD_ decay<V>::type;
        template<typename K> static K check(K, K);
        //static STD_ false check(...);
        JN_DEF_BINARY_OP(bool,               unsigned int,       unsigned int)
            JN_DEF_BINARY_OP(bool,               unsigned long,      unsigned long)
            JN_DEF_BINARY_OP(bool,               unsigned long long, unsigned long long)
            JN_DEF_BINARY_OP(bool,               short,              short)
            JN_DEF_BINARY_OP(bool,               int,                int)
            JN_DEF_BINARY_OP(bool,               long,               long)
            JN_DEF_BINARY_OP(bool,               long long,          long long)
            JN_DEF_BINARY_OP(bool,               float,              float)
            JN_DEF_BINARY_OP(bool,               double,             double)
            JN_DEF_BINARY_OP(bool,               long double,        long double)

            JN_DEF_BINARY_OP(unsigned short,     unsigned int,       unsigned int)
            JN_DEF_BINARY_OP(unsigned short,     unsigned long,      unsigned long)
            JN_DEF_BINARY_OP(unsigned short,     unsigned long long, unsigned long long)
            JN_DEF_BINARY_OP(unsigned short,     short,              short)
            JN_DEF_BINARY_OP(unsigned short,     int,                int)
            JN_DEF_BINARY_OP(unsigned short,     long,               long)
            JN_DEF_BINARY_OP(unsigned short,     long long,          long long)
            JN_DEF_BINARY_OP(unsigned short,     float,              float)
            JN_DEF_BINARY_OP(unsigned short,     double,             double)
            JN_DEF_BINARY_OP(unsigned short,     long double,        long double)

            JN_DEF_BINARY_OP(unsigned int,       unsigned long,      unsigned long)
            JN_DEF_BINARY_OP(unsigned int,       unsigned long long, unsigned long long)
            JN_DEF_BINARY_OP(unsigned int,       int,                int)
            JN_DEF_BINARY_OP(unsigned int,       long,               long)
            JN_DEF_BINARY_OP(unsigned int,       long long,          long long)
            JN_DEF_BINARY_OP(unsigned int,       float,              float)
            JN_DEF_BINARY_OP(unsigned int,       double,             double)
            JN_DEF_BINARY_OP(unsigned int,       long double,        long double)

            JN_DEF_BINARY_OP(unsigned long,      unsigned long long, unsigned long long)
            JN_DEF_BINARY_OP(unsigned long,      long,               long)
            JN_DEF_BINARY_OP(unsigned long,      long long,          long long)
            JN_DEF_BINARY_OP(unsigned long,      float,              float)
            JN_DEF_BINARY_OP(unsigned long,      double,             double)
            JN_DEF_BINARY_OP(unsigned long,      long double,        long double)

            JN_DEF_BINARY_OP(unsigned long long, long long,          long long)
            JN_DEF_BINARY_OP(unsigned long long, float,              float)
            JN_DEF_BINARY_OP(unsigned long long, double,             double)
            JN_DEF_BINARY_OP(unsigned long long, long double,        long double)

            JN_DEF_BINARY_OP(short,              int,                int)
            JN_DEF_BINARY_OP(short,              long,               long)
            JN_DEF_BINARY_OP(short,              long long,          long long)
            JN_DEF_BINARY_OP(short,              float,              float)
            JN_DEF_BINARY_OP(short,              double,             double)
            JN_DEF_BINARY_OP(short,              long double,        long double)

            JN_DEF_BINARY_OP(int,                long,               long)
            JN_DEF_BINARY_OP(int,                long long,          long long)
            JN_DEF_BINARY_OP(int,                float,              float)
            JN_DEF_BINARY_OP(int,                double,             double)
            JN_DEF_BINARY_OP(int,                long double,        long double)

            JN_DEF_BINARY_OP(long,               long long,          long long)
            JN_DEF_BINARY_OP(long,               float,              float)
            JN_DEF_BINARY_OP(long,               double,             double)
            JN_DEF_BINARY_OP(long,               long double,        long double)

            JN_DEF_BINARY_OP(long long,          float,              float)
            JN_DEF_BINARY_OP(long long,          double,             double)
            JN_DEF_BINARY_OP(long long,          long double,        long double)

            JN_DEF_BINARY_OP(float,              double,             double)
            JN_DEF_BINARY_OP(float,              long double,        long double)

            JN_DEF_BINARY_OP(double,             long double,        long double)
    public:
            using type = decltype(check(STD_ declval<F>(), STD_ declval<G>()));
    };

    template<typename T, typename V>
    struct jn_complex_op_type {
    private:
        using F = typename STD_ decay<T>::type;
        using G = typename STD_ decay<V>::type;
        template<typename _V1, typename _V2, JN_ENABLE(JN_STATIC_OR(JN_IS_COMPLEX(_V1), JN_IS_COMPLEX(_V2)))>
        static auto check(_V1, _V2)
        -> typename jn_simple_op_type<typename jn_complex_val<_V1>::type, typename jn_complex_val<_V2>::type>::type;
        static auto check(...)
            -> typename jn_simple_op_type<F, G>::type;
    public:
        using type = decltype(check(STD_ declval<F>(), STD_ declval<G>()));
    };

#define JN_BINARY_OP_TYPE(T, V) \
    typename jn_complex_op_type<typename STD_ decay<T>::type, typename STD_ decay<V>::type>::type

#define JN_ARR_BIN_RT(type1, type2) \
    Array<JN_BINARY_OP_TYPE(typename jn_array_val<type1>::type, typename jn_array_val<type2>::type)>

#define JN_DEF_ARRAY_BINARY_OP(name, name2) \
    template<typename T, typename V, JN_ENABLE(JN_STATIC_OR(JN_IS_ARRAY(T), JN_IS_ARRAY(V)))> \
    auto name(const T &a1, const V &a2) -> JN_ARR_BIN_RT(T, V) \
    { \
        using rt = JN_ARR_BIN_RT(T, V); \
        return rt::name2(a1, a2); \
    }

    JN_DEF_ARRAY_BINARY_OP(operator +,  plus         )
        JN_DEF_ARRAY_BINARY_OP(operator -,  minus        )
        JN_DEF_ARRAY_BINARY_OP(operator *,  times        )
        JN_DEF_ARRAY_BINARY_OP(operator /,  divide       )
        JN_DEF_ARRAY_BINARY_OP(operator >,  gt           )
        JN_DEF_ARRAY_BINARY_OP(operator <,  lt           )
        JN_DEF_ARRAY_BINARY_OP(operator >=, ge           )
        JN_DEF_ARRAY_BINARY_OP(operator <=, le           )
        JN_DEF_ARRAY_BINARY_OP(operator ==, eq           )

        JN_DEF_ARRAY_BINARY_OP(operator &,  bitAnd       )
        JN_DEF_ARRAY_BINARY_OP(operator |,  bitOr        )
        JN_DEF_ARRAY_BINARY_OP(operator &&, logicalAnd   )
        JN_DEF_ARRAY_BINARY_OP(operator ||, logicalOr    )

} // namespace cppsci

namespace std {

using namespace cppsci;

JN_DEF_ARRAY_BINARY_OP(pow, pow);

template<typename T, JN_ENABLE(JN_IS_ARRAY(T))>                                                                   
T log(const T &a1)                    
{                                                                                      
    return T::log(a1);                                                           
}

template<typename T, JN_ENABLE(JN_IS_ARRAY(T))>                                                                   
T exp(const T &a1)                    
{                                                                                      
    return T::exp(a1);                                                           
}

template<typename T, JN_ENABLE(JN_STATIC_AND(JN_IS_ARRAY(T), JN_IS_COMPLEX(T)))>                                                                   
T conj(const T &a1)                    
{                                                                                      
    return T::conj(a1);                                                           
}

template<typename T, JN_ENABLE(JN_STATIC_AND(JN_IS_ARRAY(T), JN_IS_COMPLEX(T)))>                                                                   
Array<typename jn_complex_val<typename jn_array_val<T>::type>::type> real(const T &a1)                    
{                                                                                      
    return T::real(a1);                                                           
}

template<typename T, JN_ENABLE(JN_STATIC_AND(JN_IS_ARRAY(T), JN_IS_COMPLEX(T)))>                                                                   
Array<typename jn_complex_val<typename jn_array_val<T>::type>::type> imag(const T &a1)                    
{                                                                                      
    return T::imag(a1);                                                           
}

template<typename T, JN_ENABLE(JN_IS_ARRAY(T))>
Array<typename jn_complex_val<typename jn_array_val<T>::type>::type> abs(const T &a1)                    
{                                                                                      
    return T::abs(a1);                                                           
}

template<typename T, JN_ENABLE(JN_IS_ARRAY(T))>                                                                   
T diag(const T &a1)                    
{                                                                                      
    return T::diag(a1);                                                           
}

} // namespace std

namespace cppsci {

inline void extend_line(double *l, int length, int size1, int size2) {
    if (size1 != 0) {
        double *ll = l + size1-1;
        double *rr = ll+1;
        int n = size1 / length;
        int r = size1 % length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < length; j++) *ll-- = *rr++;
            rr -= 2*length;
        }
        for (int j = 0; j < r; j++) *ll-- = *rr++;
    }

    if (size2 != 0) {
        double *ll = l + size1 + length -1;
        double *rr = ll+1;
        int n = size2 / length;
        int r = size2 % length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < length; j++) *rr++ = *ll--;
            ll += 2*length;
        }
        for (int j = 0; j < r; j++) *rr++ = *ll--;
    }
}

class ArrayLines : public Arrayd {
public:
    using base_type = Arrayd;

    Shape array_shape;
    Shape array_coeff;

    int axis;
    int size1;
    int size2;
    int nd;
    int length;

    ArrayLines(const Shape &_shape, int _axis, int _size1 = 0, int _size2 = 0) :
        array_shape(_shape), axis(_axis), size1(_size1), size2(_size2), 
        Arrayd({STD_ accumulate(_shape.begin(), _shape.end(), 1, [](int a, int b){return a *b;}) / _shape[_axis], _shape[_axis]+_size1+_size2})
        {
            nd = _shape.size();
            length = _shape[_axis];
            array_coeff.resize(nd);
            int l = length+size1+size2;
            int k = l;
            for (int j = nd-1; j >= 0; j--) {
                if (j == axis) array_coeff[j] = 1;
                else {
                    array_coeff[j] = k;
                    k *= _shape[j];
                }
            }
            //MPI_LOG << "array_coeff: ";
            //for (auto && n : array_coeff) MPI_LOG << n << ' '; MPI_LOG << STD_ endl;
        }

    void read_array(const Arrayd &array) {
        //        MPI_LOG << __FUNCTION__ << STD_ endl;
        //        MPI_LOG << "lss: " << length << ' ' << size1 << ' ' << size2 << STD_ endl;
        //        MPI_LOG << "shape: " << base_type::shape(0) << ' ' << base_type::shape(1) << STD_ endl;
        int m = 0;
        //        MPI_LOG << "nd: " << nd << STD_ endl;
        ARRAY_EACH(array.shape(), ind) {
            //            for (int i = 0;i < nd; i++) MPI_LOG << ind[i] << ' '; MPI_LOG << STD_ endl;
            //            MPI_LOG << "coeff1: "; for (int i = 0;i < nd; i++) MPI_LOG << array_coeff[i] << ' '; MPI_LOG << STD_ endl;
            int n = 0;
            for (int j = 0; j < nd; j++) n += array_coeff[j]*ind[j];
            //            MPI_LOG << "coeff2: "; for (int i = 0;i < nd; i++) MPI_LOG << array_coeff[i] << ' '; MPI_LOG << STD_ endl;
            //            MPI_LOG << n+size1 << ' ' << m << STD_ endl;
            //            MPI_LOG << "coeff3: "; for (int i = 0;i < nd; i++) MPI_LOG << array_coeff[i] << ' '; MPI_LOG << STD_ endl;
            //            MPI_LOG << base_type::size() << STD_ endl;
            base_type::at(n+size1) = array[m];
            //            MPI_LOG << "coeff4: "; for (int i = 0;i < nd; i++) MPI_LOG << array_coeff[i] << ' '; MPI_LOG << STD_ endl;
            m++;
        }
        double *p = base_type::data();
        for (int j = 0; j < base_type::shape(0); j++) {
            extend_line(p, length, size1, size2);
            p += length+size1+size2;
        }
    }

    void write_array(Arrayd &array) {
        int m = 0;
        ARRAY_EACH(array.shape(), ind) {
            int n = 0;
            for (int j = 0; j < nd; j++) n += array_coeff[j]*ind[j];
            array[m] = base_type::at(n+size1);
            m++;
        }
    }

};

inline Arrayi randperm(int n) {
    Arrayi outp({n});
    STD_ iota(outp.begin(), outp.end(), 0);
//    STD_ shuffle(outp.begin(), outp.end(), STD_ mt19937(STD_ random_device{}()));
    STD_ shuffle(outp.begin(), outp.end(), STD_ mt19937(::std::rand()));
    return STD_ move(outp);
}

inline Arrayi randperm(int n, int seed) {
    Rand lcg(seed);
    Arrayi outp({n});
    STD_ iota(outp.begin(), outp.end(), 0);
    for (int i = n-1; i > 0; i--) {
        int j = int(lcg()*(i+1));
        std::swap(outp[i], outp[j]);
    }
    return std::move(outp);
}

template<typename _V>
Array<_V> repmat(const Array<_V> &inp, Shape shape) {
    int nd = inp.dim();
    if (shape.size() < nd) {
        for (int i = 0; i < nd - shape.size(); i++) {
            shape.push_back(1);
        }
    }
    Shape fullshape = shape;
    for (int i = 0; i < nd; i++) fullshape[i] *= inp.shape(i);
    Array<_V> outp(fullshape);
    int n = 0;
    ARRAY_EACH(fullshape, ind) {
        Shape ind2 = ind;
        for (int i = 0; i < nd; i++) ind2[i] = ind2[i] % inp.shape(i);
        outp[n] = inp(ind2);
        n++;
    }
    return STD_ move(outp);
}

template<typename _V>
Array<_V> expand(const Array<_V> &inp, Shape shape) {
    int nd = inp.dim();
    if (shape.size() < nd) {
        for (int i = 0; i < nd - shape.size(); i++) {
            shape.push_back(1);
        }
    }
    Shape fullshape = shape;
    for (int i = 0; i < nd; i++) fullshape[i] *= inp.shape(i);
    Array<_V> outp(fullshape);
    int n = 0;
    ARRAY_EACH(fullshape, ind) {
        Shape ind2 = ind;
        for (int i = 0; i < nd; i++) ind2[i] = ind2[i] / shape[i];
        outp[n] = inp(ind2);
        n++;
    }
    return STD_ move(outp);
}

template<typename _V>
Array<_V> flipdim(const Array<_V> &a, int n) {
    Shape shape = a.shape();
    int nd = shape.size();
    int sz = a.size();
    Array<_V> b(shape);

    int area = STD_ accumulate(STD_ next(shape.begin(), n), shape.end(), 1, [](int n, int i){return n*i;});
    int slis = sz/area;
    int cols = STD_ accumulate(STD_ next(shape.begin(), n+1), shape.end(), 1, [](int n, int i){return n*i;});
    int rows = area/cols;

    const _V *a_it = a.data();
    _V *b_it = b.data();
    for (int i = 0; i < slis; i++) {
        for (int j = 0; j < rows; j++) {
            const _V *a_it2 = a_it + (rows-1-j)*cols;
            _V *b_it2 = b_it + j*cols;
            for (int k = 0; k < cols; k++) {
                *b_it2 = *a_it2;
                a_it2++;
                b_it2++;
            }
        }
        a_it += area;
        b_it += area;
    }

//    int i = 0;
//    ARRAY_EACH(shape, ind) {
//        Shape ind2 = ind;
//        ind2[n] = shape[n] - ind2[n] - 1;
//        outp[i] = inp(ind2);
//        i++;
//    }
    return STD_ move(b);
}

template<typename _V>
Array<_V> flipud(const Array<_V> &inp) {
    return flipdim(inp, 0);
}

template<typename _V>
Array<_V> fliplr(const Array<_V> &inp) {
    return flipdim(inp, 1);
}

template<typename _V>
Array<_V> flipall(const Array<_V> &inp) {
    Shape shape = inp.shape();
    int nd = inp.dim();

    Array<_V> outp(shape);
    int i = 0;
    ARRAY_EACH(shape, ind) {
        Shape ind2 = ind;
        for (int j = 0; j < nd; j++) ind2[j] = shape[j] - ind2[j] - 1;
        outp[i] = inp(ind2);
        i++;
    }
    return STD_ move(outp);
}

template<typename _V>
Array<_V> meanv(const Array<_V> &inp, int axis = 0) {
    int i, j;
    int nd = inp.dim();
    Shape ishape = inp.shape();
    Shape oshape = ishape;
    oshape[axis] = 1;

    Array<_V> outp(oshape);
    i = 0;
    Shape ii(nd);
    ARRAY_EACH(oshape, oi) {
        ii = oi;
        ii[axis] = 0;
        const _V *p = &(inp(ii));
        _V sum = 0;
        for (j = 0; j < ishape[axis]; j++) {
            sum += *p;
            p += inp.coeff(axis);
        }
        outp[i] = _V(sum / double(ishape[axis]));
        i++;
    }
    return STD_ move(outp);
}

template<typename _V>
Array<_V> sumv(const Array<_V> &inp, int axis = 0) {
    int i, j;
    int nd = inp.dim();
    Shape ishape = inp.shape();
    Shape oshape = ishape;
    oshape[axis] = 1;

    Array<_V> outp(oshape);
    i = 0;
    Shape ii(nd);
    ARRAY_EACH(oshape, oi) {
        ii = oi;
        ii[axis] = 0;
        const _V *p = &(inp(ii));
        _V sum = 0;
        for (j = 0; j < ishape[axis]; j++) {
            sum += *p;
            p += inp.coeff(axis);
        }
        outp[i] = sum;
        i++;
    }
    return STD_ move(outp);
}

template<typename _V>
Arrayi argmaxv(const Array<_V> &inp, int axis = 0) {
    int i, j;
    int nd = inp.dim();
    Shape ishape = inp.shape();
    Shape oshape = ishape;
    oshape[axis] = 1;

    Arrayi outp(oshape);
    i = 0;
    Shape ii(nd);
    ARRAY_EACH(oshape, oi) {
        ii = oi;
        ii[axis] = 0;
        const _V *p = &(inp(ii));
        _V max = *p;
        int arg = 0;
        for (j = 0; j < ishape[axis]; j++) {
            if (max < *p) {
                max = *p;
                arg = j;
            }
            p += inp.coeff(axis);
        }
        outp[i] = arg;
        i++;
    }
    return STD_ move(outp);
}

template<typename _V>
Array<_V> maxv(const Array<_V> &inp, int axis = 0) {
    int i, j;
    int nd = inp.dim();
    Shape ishape = inp.shape();
    Shape oshape = ishape;
    oshape[axis] = 1;

    Array<_V> outp(oshape);
    i = 0;
    Shape ii(nd);
    ARRAY_EACH(oshape, oi) {
        ii = oi;
        ii[axis] = 0;
        const _V *p = &(inp(ii));
        _V max = *p;
        for (j = 0; j < ishape[axis]; j++) {
            if (max < *p) max = *p;
            p += inp.coeff(axis);
        }
        outp[i] = max;
        i++;
    }
    return STD_ move(outp);
}

template<typename _V>
Arrayi argminv(const Array<_V> &inp, int axis = 0) {
    int i, j;
    int nd = inp.dim();
    Shape ishape = inp.shape();
    Shape oshape = ishape;
    oshape[axis] = 1;

    Arrayi outp(oshape);
    i = 0;
    Shape ii(nd);
    ARRAY_EACH(oshape, oi) {
        ii = oi;
        ii[axis] = 0;
        const _V *p = &(inp(ii));
        _V min = *p;
        int arg = 0;
        for (j = 0; j < ishape[axis]; j++) {
            if (min > *p) {
                min = *p;
                arg = j;
            }
            p += inp.coeff(axis);
        }
        outp[i] = arg;
        i++;
    }
    return STD_ move(outp);
}

template<typename _V>
Array<_V> minv(const Array<_V> &inp, int axis = 0) {
    int i, j;
    int nd = inp.dim();
    Shape ishape = inp.shape();
    Shape oshape = ishape;
    oshape[axis] = 1;

    Array<_V> outp(oshape);
    i = 0;
    Shape ii(nd);
    ARRAY_EACH(oshape, oi) {
        ii = oi;
        ii[axis] = 0;
        const _V *p = &(inp(ii));
        _V min = *p;
        for (j = 0; j < ishape[axis]; j++) {
            if (min > *p) min = *p;
            p += inp.coeff(axis);
        }
        outp[i] = min;
        i++;
    }
    return STD_ move(outp);
}

template<typename _A1, typename _A2>
bool approx_eq(const _A1 &a1, const _A2 &a2) {
    return Arrayd::minus(a1,a2).all([](double n){return n<1e-5;});
}

/**
 * see Matlab's mapstd function
 */
inline Arrayd mapstd(Arrayd a) {
    Shape shape = a.shape();
    int cols = shape.back();
    int n = a.size() / cols;
    Arrayd b(shape);
    double *a_it = a.data();
    double *b_it = b.data();
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < cols; j++) sum += a_it[j];
        double avg = sum / cols;
        sum = 0;
        for (int j = 0; j < cols; j++) {
            double x = a_it[j]-avg;
            sum += x*x;
        }
        double std = STD_ sqrt(sum / (cols - 1));
        for (int j = 0; j < cols; j++) {
            *b_it = (a_it[j] - avg) / std;
            b_it++;
        }
        a_it += cols;
    }
    return STD_ move(b);
}

/**
 * TODO
 */
//inline Arrayd std(Arrayd a) {}

/**
 * Return a random array with a specified shape.
 */
template<typename _A = Arrayd>
_A rand(Shape shape) {
    using value_type = typename _A::value_type;

    _A a(shape);
    for (auto && n : a) {
        n = value_type(::cppsci::rand());
    }
    return STD_ move(a);
}

template<typename RT_, typename Array_>
RT_ array_to(Array_ && inp) {
    RT_ rt(inp.size());
    STD_ copy(inp.begin(), inp.end(), rt.begin());
    return STD_ move(rt);
}

template<typename Array_>
using jn_array_take_t = SubArray<jn_if_t<STD_ is_const<Array_>::value, const jn_array_val_t<Array_>, jn_array_val_t<Array_>>>;

template<typename Array_>
using jn_array_sub_t = jn_array_take_t<Array_>;

/**
 * Return a part of the array.
 *
 * TODO
 */
template<typename Array_>
jn_array_sub_t<Array_> array_sub(Array_ && array, const STD_ vector<STD_ vector<int>> & indices) {
    jn_array_sub_t<Array_> sub;
    return STD_ move(sub);
}

/**
 * Return a part of the array.
 *
 * TODO
 */
template<typename Array_>
jn_array_take_t<Array_> array_take(Array_ && array, const STD_ vector<STD_ vector<int>> & indices) {
    jn_array_take_t<Array_> sub;
    return STD_ move(sub);
}

/**
 * Convert a list to an Array.
 */
template<typename T>
using to_array_t = decltype(*(STD_ begin(STD_ declval<T>())));

template<typename V, typename T>
Array<V> to_array(const T &t) {
    int size = STD_ distance(STD_ begin(t), STD_ end(t));
    Array<V> v({size});
    STD_ copy(STD_ begin(t), STD_ end(t), STD_ begin(v));
    return std::move(v);
}

template<typename T>
to_array_t<T> to_array(const T &t) {
    using RT = to_array_t<T>;
    return to_array<RT>(t);
}

/**
 * Read an array from file.
 */
template<typename T>
Array<T> array_read(STD_ istream &stream) {
    int ndim;
    stream >> ndim;

    Shape shape(ndim);
    for (int i = 0; i < ndim; i++) stream >> shape[i];

    Array<T> arr(shape);
    for (int i = 0; i < arr.size(); i++) stream >> arr[i];

    return STD_ move(arr);
}

template<typename T>
Array<T> array_read(STD_ string filename) {
    STD_ ifstream ifile(filename.c_str());
    auto && arr = array_read<T>(ifile);
    ifile.close();
    return STD_ move(arr);
}

/**
 * Write an array from file.
 */
template<typename _Arr>
void array_write(_Arr && arr, STD_ ostream &stream) {
    int ndim = arr.dim();
    stream << ndim << ' ';

    for (int i = 0; i < ndim; i++) stream << arr.shape(i) << ' ';

    for (int i = 0; i < arr.size(); i++) stream << arr[i] << ' ';
}

template<typename _Arr>
void array_write(_Arr && arr, STD_ string fname) {
    STD_ ofstream ofile(fname.c_str());
    array_write(arr, ofile);
    ofile.close();
}

/**
 * This is only for 2-D array!!!
 * To be generated to n-D array.
 */
template<typename T, typename _Arr1, typename _Arr2>
Array<T> vstack(const _Arr1 &a1, const _Arr2 &a2) {
    Array<T> arr({a1.shape(0)+a2.shape(0), a1.shape(1)});
    int i = 0;
    for (int j = 0; j < a1.size(); i++, j++) arr[i] = a1[j];
    for (int j = 0; j < a2.size(); i++, j++) arr[i] = a2[j];
    return STD_ move(arr);
}

/**
 * This is only for 2-D array!!!
 * To be generated to n-D array.
 */
template<typename T, typename _Arr1, typename _Arr2>
Array<T> hstack(const _Arr1 &a1, const _Arr2 &a2) {
    if (a1.size() == 0) return a2;
    else if (a2.size() == 0) return a1;

    int rows = a1.shape(0);
    int cols = a1.shape(1) + a2.shape(1);
    Array<T> arr({rows, cols});
    for (int i = 0; i < rows; i++) {
        int j = 0;
        for (int k = 0; k < a1.shape(1); k++, j++) arr(i, j) = a1(i, k);
        for (int k = 0; k < a2.shape(1); k++, j++) arr(i, j) = a2(i, k);
    }
    return STD_ move(arr);
}

/**
 * Returns an array with another shape.
 */
template<typename _T, typename _D>
Array<_T> reshape(const BasicArray<_T, _D> & a, const Shape & shape) {
    Array<_T> arr = a;
    arr.reshape(shape);
    return STD_ move(arr);
}

/**
 * Multiplies A and B.
 */
template<typename T, typename _A, typename _B>
Array<T> multiply(const Array<_A> & a, const Array<_B> & b, CBLAS_TRANSPOSE transa = CblasNoTrans, CBLAS_TRANSPOSE transb = CblasNoTrans) {

    int m, n, k;
    if (transa == CblasNoTrans) m = a.shape(0); else m = a.shape(1);
    if (transa == CblasNoTrans) k = a.shape(1); else k = a.shape(0);
    if (transb == CblasNoTrans) n = b.shape(1); else n = b.shape(0);

//    ERROR_CHECK(a.shape(1) != b.shape(0), "nCols of A should be equal to nRows of B!");

    Array<T> c({m, n});

    double alpha = 1.0;
    double beta = 0.0;

    int lda, ldb, ldc;
    if (transa == CblasNoTrans) lda = k; else lda = m;
    if (transb == CblasNoTrans) ldb = n; else ldb = k;
    ldc = n;

    //std::cout << string::format("transa:%c transb:%c m:%d n:%d k:%d lda:%d ldb:%d\n", transa, transb, m, n, k, lda, ldb);

    cblas_dgemm(CblasRowMajor, transa, transb, m, n, k, alpha, a.data(), lda, b.data(), ldb, beta, c.data(), ldc);

    return STD_ move(c);
}

} // namespace cppsci


