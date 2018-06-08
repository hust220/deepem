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

#include "cppsci_macros.hpp"
#include "cppsci_traits.hpp"
#include "cppsci_math.hpp"
#include "cppsci_string.hpp"
#include "cppsci_fft.hpp"
#include "cppsci_mpi.hpp"
#include "cppsci_array_traits.hpp"

namespace cppsci {

using Shape = ::std::vector<int>;

enum { ARRAY_TAKE_MODE_CLIP };

namespace array_helper {

    inline Shape ind_begin(const Shape &shape) {
        return Shape(shape.size(), 0);
    }

    inline Shape ind_end(const Shape &shape) {
        Shape ind = shape;
        for (auto && n : ind) n--;
        return ::std::move(ind);
    }

    inline bool ind_in(const Shape &ind, const Shape &shape) {
        int dim = int(ind.size());
        for (int i = 0; i < dim; i++) if (ind[i] < 0 || ind[i] >= shape[i]) return false;
        return true;
    }

    inline void ind_next(Shape &ind, const Shape &shape) {
        int dim = int(ind.size());
        ind[dim-1]++;
        bool flag;
        do {
            flag = false;
            for (int i = dim-1; i > 0; i--) {
                if (ind[i] >= shape[i]) {
                    ind[i] = 0;
                    ind[i-1]++;
                    flag = true;
                }
            }
        } while (flag);
    }

    inline void ind_prev(Shape &ind, const Shape &shape) {
        int dim = int(ind.size());
        ind[dim-1]--;
        bool flag;
        do {
            flag = false;
            for (int i = dim-1; i > 0; i--) {
                if (ind[i] < 0) {
                    ind[i] = shape[i]-1;
                    ind[i-1]--;
                    flag = true;
                }
            }
        } while (flag);
    }

    template<typename... Int_>
    inline void set_shape_helper(int n_ind, Shape & shape, int a, Int_ ...inds) {
        int dim = int(shape.size());
        shape[dim-n_ind] = a;
        set_shape_helper(n_ind-1, shape, inds...);
    }

    template<>
    inline void set_shape_helper<>(int n_ind, Shape & shape, int a) {
        if (n_ind != 1) DIE("Fatal error");
        int dim = int(shape.size());
        shape[dim-n_ind] = a;
    }

    template<typename... Int_>
    inline int index_helper(int n_ind, const Shape & shape, const Shape & coeff, int a, Int_ ...inds) {
        int dim = int(shape.size());
        int n = coeff[dim-n_ind];
        //        for (int i = int(dim - n_ind + 1); i < dim; i++) n *= shape[i];
        return a * n + index_helper(n_ind - 1, shape, coeff, inds...);
    }

    template<>
    inline int index_helper<>(int n_ind, const Shape & shape, const Shape & coeff, int a) {
        if (n_ind != 1) DIE("Fatal error");
        return a;
    }

    template<typename... Int_>
    inline int index(const Shape & shape, const Shape &coeff, Int_ ...inds) {
        return index_helper(int(shape.size()), shape, coeff, inds...);
    }

    inline int index(const Shape &shape, const Shape &coeff, const Shape &ind) {
        int dim = int(shape.size());
        int sum = 0;
        for (int i = 0; i < dim; i++) sum += coeff[i] * ind[i];
        return sum;
    }

} // namespace array_helper

template<typename NumType_> class SubArray;

/**
 * The basic array class.
 */
template<typename NumType_, typename Derived_>
class BasicArray {
public:
    using value_type = NumType_;
    using derived_type = Derived_;
    using type = BasicArray<value_type, derived_type>;
    using self_type = type;

    /**
     * Get ith elements of the shape of the array.
     */
    int shape(int i) const {
        return m_shape[i];
    }

    /**
     * Get the shape of the array.
     */
    const Shape &shape() const {
        return m_shape;
    }

    /**
     * Get ith elements of the coefficients of the array.
     */
    int coeff(int i) const {
        return m_coeff[i];
    }

    /**
     * Get the coefficients of the array.
     */
    const Shape &coeff() const {
        return m_coeff;
    }

    /**
     * Get the ith element of the array.
     */
    virtual value_type &at(int ind) = 0;

    /**
     * Get the ith element of a constant array.
     */
    virtual const value_type &at(int ind) const = 0;

    /** Get the element of a constant array through indices.
     */
    value_type &at(const Shape &ind) {
        return at(array_helper::index(m_shape, m_coeff, ind));
    }

    /**
     * Get the element of a constant array through indices.
     */
    const value_type &at(const Shape &ind) const {
        return at(array_helper::index(m_shape, m_coeff, ind));
    }

    /**
     * Get the element of an array through indices.
     */
    template<typename... Int_>
    value_type &at(Int_ ...inds) {
        assert(jn_num_types<Int_...>::N == dim());
        return at(array_helper::index_helper(int(m_shape.size()), m_shape, m_coeff, inds...));
    }

    /**
     * Get the element of a constant array through indices.
     */
    template<typename... Int_>
    const value_type &at(Int_ ...inds) const {
        assert(jn_num_types<Int_...>::N == dim());
        return at(array_helper::index_helper(int(m_shape.size()), m_shape, m_coeff, inds...));
    }

    value_type &operator ()(const Shape &ind) { return at(ind); }
    const value_type &operator ()(const Shape &ind) const { return at(ind); }
    template<typename... Int_> value_type &operator ()(Int_ ...inds) { return at(inds...); }
    template<typename... Int_> const value_type &operator ()(Int_ ...inds) const { return at(inds...); }

    int size() const { return m_size; }
    int dim() const { return int(m_shape.size()); }

    /**
     * Set the shape of the Array by integers.
     */
    template<typename... Int_>
    void set_shape(Int_ ...inds) {
        int dim = jn_num_types<Int_...>::N;
        m_shape.resize(dim);
        array_helper::set_shape_helper(dim, m_shape, inds...);
        set_size();
        set_coeff();
    }

    /**
     * Set the shape of the Array.
     */
    void set_shape(const Shape & _shape) {
        m_shape = _shape;
        set_size();
        set_coeff();
    }

    /**
     * Reshape the Array.
     */
    derived_type &reshape(const Shape & _shape) {
        if (m_size != ::std::accumulate(_shape.begin(), _shape.end(), 1, [](int a, int b){return a*b;})) {
            DIE("Reshape error");
        }
        set_shape(_shape);
        return *(derived_type *)this;
    }

    /**
     * Return the flattened array.
     *
     * This method is equal to reshape({size()}).
     */
    derived_type &flatten() {
        set_shape({size()});
        return *(derived_type *)this;
    }

    /**
     * Set all values.
     *
     * Set the values of all the elements of the array to a same value.
     */
    template<typename _N>
    void set_all(const _N & n) {
        for (int i = 0; i < size(); i++) at(i) = n;
    }

    /**
     * Iterator class of array.
     */
    template<typename _Arr, typename _V>
    class iterator : public ::std::iterator<::std::random_access_iterator_tag, _V> {
    public:
        using array_type = _Arr;
        using array_value_type = _V;
        using difference_type = int;
        using self = iterator<array_type, array_value_type>;

        array_type *array;
        int n;

        //iterator() : array(nullptr), n(-1) {}
        //iterator(array_type *_array) : array(_array), n(0) {}
        iterator(array_type *_array, int _n) : array(_array) { setn(_n); }
        iterator(const self &rhs) : array(rhs.array) { setn(rhs.n); }

        // Operators : misc
    public:
        inline self& operator=(const self &rhs) {setn(rhs.n); array = rhs.array; return *this;}
        inline self& operator+=(int rhs) { setn(n+rhs); return *this;}
        inline self& operator-=(int rhs) { setn(n-rhs); return *this;}
        inline array_value_type& operator*() {return array->at(n); }
        inline array_value_type* operator->() {return &(array->at(n)); }

        // Operators : arithmetic
        inline self& operator++() { setn(n+1); return *this;}
        inline self& operator--() { setn(n-1); return *this;}
        inline self operator++(int) {self tmp(*this); setn(n+1); return tmp;}
        inline self operator--(int) {self tmp(*this); setn(n-1); return tmp;}
        inline int operator-(const self& rhs) const {return n-rhs.n;}
        inline self operator+(int rhs) {return self(array, n+rhs);}
        inline self operator-(int rhs) {return self(array, n-rhs);}
        friend inline self operator+(int lhs, const self& rhs) {return self(rhs.array, lhs+rhs.n);}
        friend inline self operator-(int lhs, const self& rhs) {return self(rhs.array, lhs-rhs.n);}

        // Operators : comparison
        inline bool operator==(const self& rhs) {return n == rhs.n && array == rhs.array;}
        inline bool operator!=(const self& rhs) {return n != rhs.n || array != rhs.array;}
        inline bool operator>=(const self& rhs) {return n >= rhs.n && array == rhs.array;}
        inline bool operator<=(const self& rhs) {return n <= rhs.n && array == rhs.array;}
        inline bool operator> (const self& rhs) {return n >  rhs.n && array == rhs.array;}
        inline bool operator< (const self& rhs) {return n <  rhs.n && array == rhs.array;}

        void setn(int _n) { if (_n < 0 || _n >= array->size()) n = array->size(); else n = _n; }
    };

    /**
     * Return the begin iterator of the array.
     */
    iterator<self_type, value_type> begin() {
        return iterator<self_type, value_type>(this, 0);
    }

    /**
     * Return the begin iterator of the array.
     */
    iterator<const self_type, const value_type> begin() const {
        return iterator<const self_type, const value_type>(this, 0);
    }

    /**
     * Return the end iterator of the array.
     */
    iterator<self_type, value_type> end() {
        return iterator<self_type, value_type>(this, -1);
    }

    /**
     * Return the end iterator of the array.
     */
    iterator<const self_type, const value_type> end() const {
        return iterator<const self_type, const value_type>(this, -1);
    }

    /**
     * Return the position of the minimum element.
     */
    int argmin() const {
        int ind = 0;
        value_type min = at(0);
        for (int i = 1; i < size(); i++) if (min > at(i)) {
            min = at(i);
            ind = i;
        }
        return ind;
    }

    /**
     * Return the position of the maximum element.
     */
    int argmax() const {
        int ind = 0;
        value_type max = at(0);
        for (int i = 1; i < size(); i++) if (max < at(i)) {
            max = at(i);
            ind=i;
        }
        return ind;
    }

    /**
     * Return the number of true elements.
     */
    bool count() const {
        int n = 0;
        for (int i = 0; i < size(); i++) {
            if ((at(i))) n++;
        }
        return n;
    }

    /**
     * Return the number elements satisfying a specific condition.
     */
    template<typename _F>
    bool count(_F && f) const {
        int n = 0;
        for (int i = 0; i < size(); i++) {
            if (f(at(i))) n++;
        }
        return n;
    }

    /**
     * Return an array representing if there exists a true element along an axis.
     */
    template<typename _Val = bool>
    Array<_Val> axis_any(int axis) const {
        return axis_any(axis, [](const value_type &v)->bool{return v;});
    }

    /**
     * Return an array representing if there exists
     * an element satisfying a specific condition along an axis.
     */
    template<typename _Val = bool, typename _F>
    Array<_Val> axis_any(int axis, _F && f) const {
        int nd = dim();

        // Set shape
        Shape shape_;
        for (auto && n : shape()) if (n != axis) shape_.push_back(n);
       
        // Set coeff
        Shape coeff_(nd);
        int k = 1;
        for (int i = nd-1; i>=0; i--) {
            if (i == axis) {
                coeff_[i] = 0;
            }
            else {
                coeff_[i] = k;
                k *= shape(i);
            }
        }

        // Set value
        Array<_Val> output(shape_, false);
        int n = 0;
        ARRAY_EACH(shape(), ind) {
            if (f(at(n))) {
                int sum = 0;
                for (int i = 0; i < nd; i++) sum += coeff_[i]*ind[i];
                output[sum] = true;
            }
            n++;
        }
        return ::std::move(output);
    }

    /** Return if there exists a true element.  */
    bool any() const {
        for (int i = 0; i < size(); i++) {
            if ((at(i))) return true;
        }
        return false;
    }

    /** Return if there exists an element satisfying a specific condition.  */
    template<typename _F>
    bool any(_F && f) const {
        for (int i = 0; i < size(); i++) {
            if (f(at(i))) return true;
        }
        return false;
    }

    /** Return if all the elements are true.  */
    bool all() const {
        for (int i = 0; i < size(); i++) {
            if (!(at(i))) return false;
        }
        return true;
    }

    /** Return if all the elements satisfy a specific condition.  */
    template<typename _F>
    bool all(_F && f) const {
        for (int i = 0; i < size(); i++) {
            if (!f(at(i))) return false;
        }
        return true;
    }

    /** Return the minimum element.  */
    template<typename T = value_type>
    T min() const {
        T min  = at(0);
        for (int i = 1; i < size(); i++) if (min > at(i)) min = at(i);
        return min;
    }

    /** Return the maximum element.  */
    template<typename T = value_type>
    T max() const {
        T max  = at(0);
        for (int i = 1; i < size(); i++) if (max < at(i)) max = at(i);
        return max;
    }

    /** Return the summation.  */
    template<typename T = value_type>
    T sum() const {
        T sum  = 0;
        for (int i = 0; i < size(); i++) sum  += at(i);
        return sum;
    }

    /** Return the summation along an axis.  */
    template<typename T = value_type>
    Array<T> sum(int axis) const {
        int nd = dim();
        Shape shape2(nd-1);
        Shape coeff2(nd);
        int k = 1;
        for (int i = nd-1, j = nd-2; i >= 0; i--) {
            if (i != axis) {
                shape2[j] = shape(i);
                coeff2[i] = k;
                k *= shape(i);
                j--;
            }
            else {
                coeff2[i] = 0;
            }
        }
        Array<T> output(shape2, 0);
        k = 0;
        ARRAY_EACH(shape(), ind) {
            int n = 0;
            for (int i = 0; i < nd; i++) n += ind[i] * coeff2[i];
            output[n] += at(k);
            k++;
        }
        //output.selfDivide(shape(axis));
        return ::std::move(output);
    }

    /** Return the mean.  */
    template<typename T = value_type>
    T mean() const {
        T mean = 0;
        for (int i = 0; i < size(); i++) mean += at(i);
        return T(mean / (double)size());
    }

    /** Return the variance.  */
    template<typename T = value_type>
    T var() const {
        T m = mean();
        T sum = 0;
        for (int i = 0; i < size(); i++) sum += ::std::pow(at(i)-m, 2);
        return T(sum / (double)size());
    }

    /** Return the norm.  */
    template<typename T = value_type>
    T norm()  const {
        T norm = 0;
        for (int i = 0; i < size(); i++) norm += square(at(i));
        return ::std::sqrt(norm);
    }

    /** Return the trace.  */
    template<typename T = value_type>
    T trace() const {
        int sum = 0;
        for (auto && n : coeff()) sum += n;
        T trace = 0;
        for (int i = 0; i < size(); i += sum) trace += at(i);
        return trace;
    }

    /** Return the transposed array.  */
    template<typename _RT = Array<value_type>>
    _RT transpose() {
        int nd = dim();
        int a = shape(0);
        int b = shape(1);
        ERROR_CHECK(nd!=2, "Transpose operation is only supported by 2-dimensional matrix now!");
        _RT n({b,a});
        for (int i = 0; i < b; i++) for (int j = 0; j < a; j++) {
            n(i, j) = at(j, i);
        }
        return ::std::move(n);
    }

    /** Transpose the array itself.  */
    derived_type &selfTranspose() {
        int nd = dim();
        int a = shape(0);
        int b = shape(1);
        ERROR_CHECK(a!=b || nd!=2, "Transpose operation is only supported by 2-dimensional square matrix now!");
        for (int i = 0; i < a; i++) for (int j = i+1; j < b; j++) {
            ::std::swap(at(i,j), at(j,i));
        }
        return *(derived_type *)this;
    }

    /** Reverse the array itself.  */
    derived_type &selfReverse() {
        int l = size();
        for (int i = 0; i < l/2; i++) ::std::swap(at(i), at(l-1-i)); 
        return *(derived_type*)this;
    }

    /** Apply a function to every element of the array itself.  */
    template<typename _F>
    derived_type &selfApply(_F && f) {
        for (int i = 0; i < size(); i++) at(i) = f(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's conjugation.  */
    derived_type &selfConj() {
        for (int i = 0; i < size(); i++) at(i) = ::std::conj(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's square.  */
    derived_type &selfSquare() {
        for (int i = 0; i < size(); i++) at(i) *= at(i);
        return *(derived_type*)this;
    }

    /** Change the array itself to it's exponation.  */
    derived_type &selfExp   () {
        for (int i = 0; i < size(); i++) at(i) = ::std::exp(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's logarithm.  */
    derived_type &selfLog   () {
        for (int i = 0; i < size(); i++) at(i) = ::std::log(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's squared root.  */
    derived_type &selfSqrt  () {
        for (int i = 0; i < size(); i++) at(i) = ::std::sqrt(at(i));
        return *(derived_type*)this;
    }

/**
 * @name Arithmetic methods
 * @{
 */
#define ARRAY_DEF_METHOD(name) \
    template<typename _N> auto name(const _N & n) \
        -> jn_enable_t<JN_IS_ARRAY(_N), derived_type &> \
    { \
        for (int i = 0; i < size(); i++) at(i) = ARRAY_METHOD(at(i), n(i)); \
        return *(derived_type*)this; \
    } \
    template<typename _N> auto name(const _N & n) \
        -> jn_enable_t<JN_STATIC_NOT(JN_IS_ARRAY(_N)), derived_type &> \
    { \
        for (int i = 0; i < size(); i++) at(i) = ARRAY_METHOD(at(i), n); \
        return *(derived_type*)this; \
    }

#define ARRAY_METHOD(m, n) m + n
        ARRAY_DEF_METHOD(selfPlus)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m + n
        ARRAY_DEF_METHOD(operator +=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m - n
        ARRAY_DEF_METHOD(selfMinus)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m - n
        ARRAY_DEF_METHOD(operator -=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) n - m
        ARRAY_DEF_METHOD(selfMinusedBy)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m * n
        ARRAY_DEF_METHOD(selfTimes)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m * n
        ARRAY_DEF_METHOD(operator *=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m / n
        ARRAY_DEF_METHOD(selfDivide)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m / n
        ARRAY_DEF_METHOD(operator /=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (n / m)
        ARRAY_DEF_METHOD(selfDividedBy)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m && n)
        ARRAY_DEF_METHOD(selfLogicalAnd)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m || n)
        ARRAY_DEF_METHOD(selfLogicalOr)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m & n)
        ARRAY_DEF_METHOD(selfBitAnd)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m | n)
        ARRAY_DEF_METHOD(selfBitOr)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) ::std::pow(m, n)
        ARRAY_DEF_METHOD(selfPow)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) ::std::pow(n, m)
        ARRAY_DEF_METHOD(selfPowedBy)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m > n)
        ARRAY_DEF_METHOD(selfGt)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m < n)
        ARRAY_DEF_METHOD(selfLt)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m == n)
        ARRAY_DEF_METHOD(selfEq)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m >= n)
        ARRAY_DEF_METHOD(selfGe)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m <= n)
        ARRAY_DEF_METHOD(selfLe)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) n
        ARRAY_DEF_METHOD(operator =)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) n
        ARRAY_DEF_METHOD(assign)
#undef ARRAY_METHOD

#undef ARRAY_DEF_METHOD
/**
 * @}
 */

    /** Print the identification.  */
    void identify(::std::string name) const {
#ifdef USEMPI
        if (!MPI_IS_ROOT) return;
#endif
        ::std::cout << name << " (";
        for (auto && n : shape()) ::std::cout << n << ' ';
        ::std::cout << ')'
            << " num:" << size()
            << " min:" << min()
            << " max:" << max()
            << " sum:" << sum<double>()
            << " mean:" << mean<double>()
            << ::std::endl;
    }

    /** Pretty print the array.  */
    void print(::std::ostream &stream = ::std::cout) const {
#ifdef USEMPI
        if (!MPI_IS_ROOT) return;
#endif
        stream << ::std::right;
        int dim = int(shape().size());
        const auto & w = shape();
        ::std::vector<int> v(dim);
        int n = 1;
        for (int i = dim - 1; i >= 0; i--) {
            n *= shape(i);
            v[i] = n;
        }
        for (int i = 0; i < size(); i++) {
            for (int j = 0; j < dim; j++) {
                if (i % v[j] == 0) {
                    stream << '[';
                }
                else if (i % v[dim-1] == 0) {
                    stream << ' ';
                    int d = (i % v[j]) / v[j+1];
                    if (d >= 3 && d < w[j]-3) {
                        i += v[j+1]*(w[j]-6)-1;
                        stream << "...\n";
                        goto array_print_for;
                    }
                }
            }
            //if (i % v[dim-1] != v[dim-1]-1) stream << ::std::setw(5);
            stream << ::std::setw(5);
            if (i % v[dim-1] >= 3 && i % v[dim-1] < w[dim-1]-3) {
                stream << "...";
                i += w[dim-1]-7;
                continue;
            }
            else {
                stream << at(i);
            }
            for (int j = dim - 1; j >= 0; j--) {
                if (i % v[j] == v[j] - 1) stream << ']';
            }
            if (i % v[dim-1] == v[dim-1]-1) stream << '\n';
            else stream << ' ';
array_print_for:;
        }
    }

    /** Print all the elements of the array.  */
    void printFull(::std::ostream &stream = ::std::cout) const {
#ifdef USEMPI
        if (!MPI_IS_ROOT) return;
#endif
        int k = m_shape[dim()-1];
        for (int i = 0; i < size(); i++) {
            stream << at(i);
            if (i % k == k-1) stream << "\n";
            else stream << "\t";
        }
    }

protected:
    int m_size;
    Shape m_shape;
    Shape m_coeff;

    /**
     * Set the size of the Array.
     *
     * This meethod is equal to reshape
     */
    void set_size() {
        m_size = 1;
        for (auto &&n : m_shape) {
            m_size *= n;
        }
    }

    /**
     * Set the coefficients of the Array.
     */
    void set_coeff() {
        int k = 1;
        int d = dim();
        m_coeff.resize(d);
        for (int i = d - 1; i >= 0; i--) {
            m_coeff[i] = k;
            k *= m_shape[i];
        }
    }

    template<typename _T, typename _V>
    int get_ind(_T && _inds, _V && _shape) {
        int dim = _inds.size();
        int n = 0;
        int m = 1;
        for (int i = dim - 1; i >= 0; i--) {
            n += _inds[i] * m;
            m *= _shape[i];
        }
        return n;
    }

    template<typename _T, typename _V>
    bool next_ind(_T && _inds, _V && _range) {
        _inds.back()++;
        int dim = _inds.size();
        bool flag;
        do {
            flag = false;
            for (int i = dim - 1; i > 0; i--) {
                if (_inds[i] >= _range[i][1]) {
                    _inds[i] = _range[i][0];
                    _inds[i-1]++;
                    flag = true;
                    break;
                }
            }
        } while (flag);
        return _inds[0] < _range[0][1];
    }

};

/** Operator << reload for array.
 */
template<typename Derived_, typename _Val>
::std::ostream &operator <<(::std::ostream &stream, const BasicArray<_Val, Derived_> &array) {
    array.print(stream);
    return stream;
}

/** Operator == reload for array.
 */
template<typename _V, typename _D1, typename _D2>
bool operator ==(const BasicArray<_V, _D1> &a1, const BasicArray<_V, _D2> &a2) {
    if (a1.shape() != a2.shape()) return false;
    for (int i = 0; i < a1.size(); i++) {
        if (a1(i) != a2(i)) return false;
    }
    return true;
}

/** Class Array.
 */
template<typename NumType_>
class Array : public BasicArray<NumType_, Array<NumType_>> {
public:
    using value_type = NumType_;
    using type = Array<value_type>;
    using self_type = Array<value_type>;
    using base_type = BasicArray<value_type, self_type>;

    virtual value_type &at(int a) { return m_data[a]; }
    virtual const value_type &at(int a) const { return m_data[a]; }

    value_type *data() { return m_data; }
    const value_type *data() const { return m_data; }

    value_type &operator [](int i) { return m_data[i]; }
    const value_type &operator [](int i) const { return m_data[i]; }

    /**
     * @name Array constructors and destructors
     * @{
     */

    /**
     * Default constructor.
     */
    Array() {
        init();
    }

    /**
     * Copy constructor.
     */
    Array(const self_type & array) {
        init();
        copy(array);
    }

    /**
     * Move constructor.
     */
    Array(self_type && array) {
        init();
        move(array);
    }

    /**
     * Extended copy constructor.
     */
    template<typename T_, typename DT_>
    Array(const BasicArray<T_, DT_> &array) {
        init();
        copy(array);
    }

    /**
     * Copy assignment.
     */
    self_type &operator =(const self_type &array) {
        copy(array);
        return *this;
    }

    /**
     * Move assignment.
     */
    self_type &operator =(self_type &&array) {
        move(array);
        return *this;
    }

    /**
     * Extended copy assignment.
     */
    template<typename _T, typename _DT>
    self_type &operator =(const BasicArray<_T, _DT> &array) {
        copy(array);
        return *this;
    }

    /**
     * Constructor from a shape with the type of initializer_list.
     *
     * @code
     * Array array({2,2});
     * @endcode
     */
    Array(const ::std::initializer_list<int> & _shape) {
        base_type::set_shape(Shape(_shape));
        m_data = new value_type[base_type::size()];
    }

    /**
     * Constructor from a shape.
     *
     * @code
     * Shape shape{2,2};
     * Array array(shape);
     * @endcode
     */
    Array(const Shape & _shape) {
        base_type::set_shape(_shape);
        m_data = new value_type[base_type::size()];
    }

    /**
     * Constructor from multiple integer parameters.
     *
     * @code
     * Array array(2,2);
     * @endcode
     */
    template<typename... Int_, typename = jn_enable_t<STD_ is_same<jn_first_t<Int_...>, int>::value, int>>
    Array(Int_ ...inds) {
        base_type::set_shape(inds...);
        m_data = new value_type[base_type::size()];
    }

    /**
     * Array constructor from shape and data.
     *
     * @code
     * Array array({2,2},{1,2,3,4});
     * @endcode
     */
    Array(const Shape & _shape, const ::std::initializer_list<value_type> & _data) {
        base_type::set_shape(_shape);
        m_data = new value_type[base_type::size()];
        int j = 0; for (auto && n : _data) { m_data[j] = n; j++; }
    }

    /**
     * Array constructor from shape and value.
     *
     * @code
     * Array array({2,2},0);
     * @endcode
     */
    Array(const Shape & shape_, const value_type & value_) {
        base_type::set_shape(shape_);
        int n = base_type::size();
        m_data = new value_type[n];
        for (int i = 0; i < n; i++) m_data[i] = value_;
    }

    /**
     * Array constructor from range with the type of ::std::vector<::std::vector<int>>.
     *
     * @code
     * Array a({3,3},{1,2,3,4,5,6,7,8,9});
     * Array b(a, {{1,2},{2,3}});
     * @endcode
     */
    template<typename _V, typename _D>
    Array(const BasicArray<_V, _D> & array, const ::std::vector<::std::vector<int>> &range) {
        int dim = array.dim();
        if (range.size() != dim) DIE("Array constructor error");

        Shape shape(dim);
        for (int i = 0; i < dim; i++) shape[i] = range[i][1] - range[i][0];
        base_type::set_shape(shape);
        m_data = new value_type[base_type::size()];

        ::std::vector<int> inds(dim);
        for (int i = 0; i < dim; i++) inds[i] = range[i][0];
        int n = 0;
        do {
            m_data[n] = array(ind(inds, array.shape()));
            n++;
        } while (base_type::next_ind(inds, range));
    }

    /**
     * Array constructor from range with the type of ::std::vector<::std::vector<int>>.
     *
     * @code
     * Array a({3,3},{1,2,3,4,5,6,7,8,9});
     * Arrayi range({3}, {1,3,4,6});
     * Array b(a, range);
     * @endcode
     */
    template<typename _V, typename _D>
    Array(const BasicArray<_V, _D> & array, const Array<int> &range) {
        base_type::set_shape(range.shape());
        m_data = new value_type[base_type::size()];
        for (int i = 0; i < base_type::size(); i++) m_data[i] = array(range(i));
    }

    ~Array() {
        clear();
    }
    /**
     * @}
     */

    /**
     * @name FFT (fast fourior transformation) methods
     * @{
     */
    self_type &selfFFT()       {       ::cppsci::fft(data(), data(),                   base_type::shape( )); return *this; }
    self_type &selfIFFT()      {      ::cppsci::ifft(data(), data(),                   base_type::shape( )); return *this; }
    self_type &selfFFTshift()  {  ::cppsci::fftshift(data(), data(), base_type::dim(), base_type::shape(0)); return *this; }
    self_type &selfFFT1()      {      ::cppsci::fft1(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfIFFT1()     {     ::cppsci::ifft1(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfFFTshift1() { ::cppsci::fftshift1(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfFFT2()      {      ::cppsci::fft2(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfIFFT2()     {     ::cppsci::ifft2(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfFFTshift2() { ::cppsci::fftshift2(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfFFT3()      {      ::cppsci::fft3(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfIFFT3()     {     ::cppsci::ifft3(data(), data(),                   base_type::shape(0)); return *this; }
    self_type &selfFFTshift3() { ::cppsci::fftshift3(data(), data(),                   base_type::shape(0)); return *this; }
    template<typename T> static self_type       fft(T && m) { self_type n(m.shape());       ::cppsci::fft(m.data(), n.data(),          m.shape( )); return ::std::move(n); }
    template<typename T> static self_type      ifft(T && m) { self_type n(m.shape());      ::cppsci::ifft(m.data(), n.data(),          m.shape( )); return ::std::move(n); }
    template<typename T> static self_type  fftshift(T && m) { self_type n(m.shape());  ::cppsci::fftshift(m.data(), n.data(), m.dim(), m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type      fft1(T && m) { self_type n(m.shape());      ::cppsci::fft1(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type     ifft1(T && m) { self_type n(m.shape());     ::cppsci::ifft1(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type fftshift1(T && m) { self_type n(m.shape()); ::cppsci::fftshift1(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type      fft2(T && m) { self_type n(m.shape());      ::cppsci::fft2(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type     ifft2(T && m) { self_type n(m.shape());     ::cppsci::ifft2(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type fftshift2(T && m) { self_type n(m.shape()); ::cppsci::fftshift2(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type      fft3(T && m) { self_type n(m.shape());      ::cppsci::fft3(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type     ifft3(T && m) { self_type n(m.shape());     ::cppsci::ifft3(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    template<typename T> static self_type fftshift3(T && m) { self_type n(m.shape()); ::cppsci::fftshift3(m.data(), n.data(),          m.shape(0)); return ::std::move(n); }
    /**
     * @}
     */

    /**
     * @name Array arithmetic methods
     *
     * @code
     * Array a({3},{1,2,3});
     * Array b({3},{3,4,5});
     * ::std::cout << Array::plus(a, b) << ::std::endl;
     * // [4, 6, 8]
     * @endcode
     *
     * @{
     */
#define ARRAY_DEF_METHOD(name) \
    template<typename _M, typename _N> \
    static auto name(const _M & m, const _N & n) \
        -> jn_enable_t<JN_STATIC_AND(JN_IS_ARRAY(_M), JN_IS_ARRAY(_N)), type> \
    { \
        type r(m.shape()); \
        for (int i = 0; i < r.size(); i++) r(i) = ARRAY_METHOD(m(i), n(i)); \
        return ::std::move(r);\
    } \
\
    template<typename _M, typename _N> \
    static auto name(const _M & m, const _N & n) \
        -> jn_enable_t<JN_STATIC_AND(JN_IS_ARRAY(_M), JN_STATIC_NOT(JN_IS_ARRAY(_N))), type> \
    { \
        type r(m.shape()); \
        for (int i = 0; i < r.size(); i++) r(i) = ARRAY_METHOD(m(i), n); \
        return ::std::move(r);\
    } \
\
    template<typename _M, typename _N> \
    static auto name(const _M & m, const _N & n) \
        -> jn_enable_t<JN_STATIC_AND(JN_STATIC_NOT(JN_IS_ARRAY(_M)), JN_IS_ARRAY(_N)), type> \
    { \
        type r(n.shape()); \
        for (int i = 0; i < r.size(); i++) r(i) = ARRAY_METHOD(m, n(i)); \
        return ::std::move(r);\
    }

#define ARRAY_METHOD(m, n) m + n
    ARRAY_DEF_METHOD(plus)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m - n
    ARRAY_DEF_METHOD(minus)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m * n
    ARRAY_DEF_METHOD(times)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m / n
    ARRAY_DEF_METHOD(divide)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m > n
    ARRAY_DEF_METHOD(gt)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m >= n
    ARRAY_DEF_METHOD(ge)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m < n
    ARRAY_DEF_METHOD(lt)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m <= n
    ARRAY_DEF_METHOD(le)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m == n
    ARRAY_DEF_METHOD(eq)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m & n
    ARRAY_DEF_METHOD(bitAnd)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m | n
    ARRAY_DEF_METHOD(bitOr)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m && n
    ARRAY_DEF_METHOD(logicalAnd)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m || n
    ARRAY_DEF_METHOD(logicalOr)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) ::std::pow(m, n)
    ARRAY_DEF_METHOD(pow)
#undef ARRAY_METHOD

#undef ARRAY_DEF_METHOD

    /**
     * @}
     */

    /**
     * Return an array with all the elements being mapped by the function.
     */
    template<typename _F, typename _M1>
    static type map(_F && f, _M1 && m1) {
        type r(m1.shape()); 
        for (int i = 0; i < r.size(); i++) r[i] = f(m1[i]); 
        return ::std::move(r);
    }

    /**
     * Return an array with all the elements being mapped by the function.
     */
    template<typename _F, typename _M1, typename _M2>
    static type map(_F && f, _M1 && m1, _M2 && m2) {
        ERROR_CHECK(m1.shape() != m2.shape(), "The shape of the two arrays should be same!");
        type r(m1.shape()); 
        for (int i = 0; i < r.size(); i++) r[i] = f(m1(i), m2(i)); 
        return ::std::move(r);
    }

    /**
     * Same as map(f, m1)
     */
    template<typename _M, typename _F>
    static type apply(_M && m1, _F && f) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = f(m1(i)); 
        return ::std::move(m2);
    }

    /**
     * Return a new array with all the elements being reversed
     */
    template<typename _M, typename _F> static type reverse(_M && m1) {
        type m2(m1.shape()); 
        int l = m1.size();
        for (int i = 0; i < l; i++) m2(i) = m1(l-1-i); 
        return ::std::move(m2);
    }

    /**
     * Return a conjugate array.
     */
    template<typename _M> static type conj   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::conj(m1(i)); 
        return ::std::move(m2);
    }

    /**
     * Return the real part of the array.
     */
    template<typename _M> static type real   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::real(m1(i)); 
        return ::std::move(m2);
    }

    /**
     * Return the imaginary part of the array.
     */
    template<typename _M> static type imag   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::imag(m1(i)); 
        return ::std::move(m2);
    }

    /**
     * Return the square of the array.
     */
    template<typename _M> static type square (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = m1(i)*m1(i); 
        return ::std::move(m2);
    }

    /**
     * Return the exponation of the array.
     */
    template<typename _M> static type exp    (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::exp(m1(i)); 
        return ::std::move(m2);
    }

    /**
     * Return the logarithm of the array.
     */
    template<typename _M> static type log    (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::log(m1(i)); 
        return ::std::move(m2);
    }

    /**
     * Return the squared root of the array.
     */
    template<typename _M> static type sqrt   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::sqrt(m1(i)); 
        return ::std::move(m2);
    }

    /**
     * Return the not of the array.
     */
    template<typename _M> static type logicalNot(_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = !m1(i); 
        return ::std::move(m2);
    }

    /**
     * Return the absolution of the array.
     */
    template<typename _M> static type abs(const _M & m1) {
        type m2(m1.shape());
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::abs(m1(i));
        return ::std::move(m2);
    }

    /**
     * Return the difference of the array.
     */
    template<typename _M> static type diff(const _M & m1) {
        type m2({m1.size()-1});
        for (int i = 0; i < m2.size(); i++) m2(i) = m1(i+1)-m1(i);
        return ::std::move(m2);
    }

    /**
     * Return the diagonized array.
     */
    template<typename _M> static type diag(const _M & m1) {
        if (m1.dim() > 1) {
            int l = *::std::min_element(m1.shape().begin(), m1.shape().end());
            type m2({l});
            int k = 1; for (int i = 1; i < l; i++) k *= m1.shape(i);
            for (int i = 0, j = 0; i < l; i++, j += k) {
                m2[i] = m1[k];
            }
            return ::std::move(m2);
        }
        else if (m1.dim() == 1) {
            int n = m1.size();
            type m2({n, n}, 0);
            for (int i = 0; i < n; i++) m2(i, i) = m1[i];
            return ::std::move(m2);
        }
        else {
            ERROR_REPORT("It's an empty array!");
        }
    }

    /**
     * Extract array satisfying the condition.
     */
    template<typename _M>
    static type extract(const Array<bool> &condition, _M && m) {
        int n = 0;
        for (int i = 0; i < condition.size(); i++) if (condition[i]) n++;
        type r({n});
        for (int i = 0, j = 0; i < condition.size(); i++) if (condition[i]) {
            r[j] = m[i];
            j++;
        }
        return ::std::move(r);
    }

    /**
     * Create an array satisfying the range.
     */
    static type range(double a, double b, double step = 1, int n = 0) {
        if (step != 0 && n == 0) n = int(::std::ceil((b-a)/step));
        else if (step == 0 && n != 0) step = (b-a)/(n-1);
        else throw "Array::range error!";

        type m(n);
        for (int i = 0; i < n; i++) m[i] = value_type(a + i * step);
        return ::std::move(m);
    }

    /**
     * Create an array satisfying the range.
     */
    static type range(int a) {
        type m({a});
        for (int i = 0; i < a; i++) m[i] = value_type(i);
        return ::std::move(m);
    }

    /**
     * Create an array satisfying the range.
     */
    template<typename T>
    static type rangev(const T &t) {
        int n = STD_ distance(STD_ begin(t), STD_ end(t));
        if (n == 0) return range(0);
        else if (n == 1) return range(t[0]);
        else if (n == 2) return range(t[0], t[1]);
        else if (n == 3) return range(t[0], t[1], t[2]);
        else return range(t[0], t[1], t[2], t[3]);
    }

    /**
     * Create an array satisfying the space.
     */
    static type linspace(double start, double stop, int num) {
        type m(num);
        double d = (stop - start) / (num - 1);
        for (int i = 0; i < num; i++) m(i) = start + i * d;
        return ::std::move(m);
    }

    /**
     * Create indices satisfying the condition.
     */
    static type where(const Array<bool> &condition) {
        ::std::vector<Shape> ls;
        int i = 0;
        ARRAY_EACH(condition.shape(), ind) {
            if (condition[i]) ls.push_back(ind);
            i++;
        }
        type m({condition.dim(), int(ls.size())});
        i = 0;
        ARRAY_EACH(m.shape(), ind) {
            m[i] = ls[ind[1]][ind[0]];
            i++;
        }
        return ::std::move(m);
    }

    /**
     * Return the convolve of two arrays.
     */
    template<typename _M, typename _N>
    static type convolve(const _M &m, const _N &n) {
        Array<::std::complex<double>> cm = m, cn = n;
        cm.selfFFT();
        return type::abs(cn.selfFFT().selfConj().selfTimes(cm).selfIFFT());
    }

protected:

    value_type *m_data = nullptr;

    /**
     * @name Don't use theses methods!!!
     * @{
     */
    void init() {
        base_type::set_shape({0});
        m_data = nullptr;
    }

    template<typename... Int_>
    void realloc(Int_ ..._inds) {
        clear();
        base_type::set_shape(_inds...);
        m_data = new value_type[base_type::size()];
    }

    void realloc(const Shape &_shape) {
        clear();
        base_type::set_shape(_shape);
        m_data = new value_type[base_type::size()];
    }

    void clear() {
        if (m_data != nullptr) delete [] m_data;
        init();
    }

    template<typename _T>
    void copy(const _T &array) {
        if (base_type::size() != array.size()) realloc(array.shape());
        for (int i = 0; i < base_type::size(); i++) m_data[i] = array(i);
    }

    template<typename _T>
    void move(_T &array) {
        std::swap(base_type::m_shape, array.m_shape);
        std::swap(base_type::m_coeff, array.m_coeff);
        std::swap(base_type::m_size, array.m_size);
        std::swap(m_data, array.m_data);
    }

    /**
     * @}
     */

};

/**
 * Map a pointer to an array
 */
template<typename NumType_>
class MapArray : public BasicArray<NumType_, MapArray<NumType_>> {
public:
    using value_type = NumType_;
    using type = MapArray<value_type>;
    using self_type = MapArray<value_type>;
    using base_type = BasicArray<value_type, type>;

    virtual value_type &at(int a) { return m_data[a]; }
    virtual const value_type &at(int a) const { return m_data[a]; }

    value_type *data() { return m_data; }
    const value_type *data() const { return m_data; }

    MapArray(value_type * data_, const Shape & shape_) {
        m_data = data_;
        base_type::set_shape(shape_);
    }

    template<typename... Int_, typename = jn_enable_t<STD_ is_same<jn_first_t<Int_...>, int>::value, int>>
    MapArray(value_type * data_, Int_ ...inds_) {
        m_data = data_;
        base_type::set_shape(inds_...);
    }

protected:
    value_type *m_data;

};

/**
 * Reshape an array.
 */
template<typename NumType_, typename _OrigArray>
class ReshapeArray : public BasicArray<NumType_, ReshapeArray<NumType_, _OrigArray>> {
public:
    using value_type = NumType_;
    using orig_array_type = _OrigArray;

    using type = ReshapeArray<value_type, orig_array_type>;
    using self_type = ReshapeArray<value_type, orig_array_type>;
    using base_type = BasicArray<value_type, type>;

    _OrigArray *orig_array;

    ReshapeArray(_OrigArray * _orig_array, const Shape & _shape) {
        orig_array = _orig_array;
        base_type::set_shape(_shape);
    }

    virtual value_type &at(int a) { return orig_array->at(a); }
    virtual const value_type &at(int a) const { return orig_array->at(a); }
};

template<typename T>
struct array_reshape_m {
    using value_type = typename STD_ remove_reference<decltype(STD_ declval<T>().at(0))>::type;
    using type = ReshapeArray<jn_if_t<STD_ is_const<T>::value, const value_type, value_type>, typename STD_ remove_reference<T>::type>;
};
template<typename T>
using array_reshape_t = typename array_reshape_m<T>::type;

template<typename T>
array_reshape_t<T> array_reshape(T &&t, Shape shape) {
    return array_reshape_t<T>(&t, shape);
}

/**
 * @name Class Redefinitions
 * @{
 */
using Arrayb  = Array<bool>;
using Arrayc  = Array<char>;
using Arrays  = Array<::std::string>;
using Arrayi  = Array<int>;
using Arrayf  = Array<float>;
using Arrayd  = Array<double>;
using Arrayci = Array<::std::complex<int>>;
using Arraycf = Array<::std::complex<float>>;
using Arraycd = Array<::std::complex<double>>;

using MapArrayb  = MapArray<bool>;
using MapArrayc  = MapArray<char>;
using MapArrays  = MapArray<::std::string>;
using MapArrayi  = MapArray<int>;
using MapArrayf  = MapArray<float>;
using MapArrayd  = MapArray<double>;
using MapArrayci = MapArray<::std::complex<int>>;
using MapArraycf = MapArray<::std::complex<float>>;
using MapArraycd = MapArray<::std::complex<double>>;

using SubArrayb  = SubArray<bool>;
using SubArrayc  = SubArray<char>;
using SubArrays  = SubArray<::std::string>;
using SubArrayi  = SubArray<int>;
using SubArrayf  = SubArray<float>;
using SubArrayd  = SubArray<double>;
using SubArrayci = SubArray<::std::complex<int>>;
using SubArraycf = SubArray<::std::complex<float>>;
using SubArraycd = SubArray<::std::complex<double>>;

/**
 * @}
 */

} // namespace cppsci

