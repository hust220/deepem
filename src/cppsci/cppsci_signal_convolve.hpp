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

#include "cppsci_array.hpp"
#include "cppsci_fft.hpp"
#include "cppsci_time.hpp"

namespace cppsci { namespace signal {

    enum {
        CONVOLVE_PLAIN,
        CONVOLVE_FFT
    };

    enum {
        CONVOLVE_SAME,
        CONVOLVE_VALID,
        CONVOLVE_FULL
    };

    template<typename _A, typename _B, typename _C>
    struct Convolve {
        Array<_C> c;

        Convolve(const Array<_A> &a, const Array<_B> &b, int result_mode, int compute_mode) {
            if (compute_mode == CONVOLVE_PLAIN) run_plain(a, b, result_mode);
            else if (compute_mode == CONVOLVE_FFT) run_fft(a, b, result_mode);
            else ERROR_REPORT("Convolve compute mode error!");
        }

        void run_plain(const Array<_A> &a, const Array<_B> &b, int result_mode) {
            set_shape(a.shape(), b.shape());

            if (result_mode == CONVOLVE_VALID) {
                run_plain_valid(a, b);
            }
            else if (result_mode == CONVOLVE_SAME) {
                for (int i = 0; i < nd; i++) shape_c.push_back(shape_a[i]);
                ERROR_REPORT("Have not been implemented!");
            }
            else if (result_mode == CONVOLVE_FULL) {
                run_plain_full(a, b);
            }
        }

        void run_plain_valid(const Array<_A> &a, const Array<_B> &b) {
            if (!shape_c.empty()) shape_c.clear();
            for (int i = 0; i < nd; i++) shape_c.push_back(shape_a[i]-shape_b[i]+1);
            c = Array<_C>(shape_c, 0);

            if (nd == 2) {
                int a_rows = shape_a[0];
                int a_cols = shape_a[1];
                int b_rows = shape_b[0];
                int b_cols = shape_b[1];
                int c_rows = shape_c[0];
                int c_cols = shape_c[1];

                int a_size = a.size();
                int b_size = b.size();
                int c_size = c.size();

                const _B *b_end = b.data() + b_size - 1;

                const _A *a_it = a.data();
                _C *c_it = c.data();

                int incx = 1;
                int incy = 1;

                // Using daxpy is faster!!!
                for (int c_row = 0; c_row != c_rows; c_row++) {
                    const _A *a_it2 = a_it;
                    const _B *b_it = b_end;
                    for (int b_row = 0; b_row != b_rows; b_row++) {
                        for (int b_col = 0; b_col != b_cols; b_col++) {
                            daxpy(&c_cols, b_it, a_it2, &incx, c_it, &incy);
                            b_it--;
                            a_it2++;
                        }
                        a_it2 += a_cols - b_cols;
                    }
                    c_it += c_cols;
                    a_it += a_cols;
                }

                // Below codes are slower than that using daxpy!
//                    for (int c_row = 0; c_row != c_rows; c_row++) {
//                        for (int c_col = 0; c_col != c_cols; c_col++) {
//                            _C sum = 0;
//                            const _A *a_it2 = a_it;
//                            const _B *b_it = b_end;
//                            for (int b_row = 0; b_row != b_rows; b_row++) {
//                                for (int b_col = 0; b_col != b_cols; b_col++) {
//                                    sum += (*a_it2) * (*b_it);
//                                    b_it--;
//                                    a_it2++;
//                                }
//                                a_it2 += a_cols - b_cols;
//                            }
//                            (*c_it) = sum;
//                            a_it++;
//                            c_it++;
//                        }
//                        a_it += b_cols - 2;
//                    }
            }
            else if (nd == 3) {
                int a_slis = shape_a[0];
                int a_rows = shape_a[1];
                int a_cols = shape_a[2];

                int b_slis = shape_b[0];
                int b_rows = shape_b[1];
                int b_cols = shape_b[2];

                int c_slis = shape_c[0];
                int c_rows = shape_c[1];
                int c_cols = shape_c[2];

                int a_size = a.size();
                int b_size = b.size();
                int c_size = c.size();

                const _B *b_end = b.data() + b_size - 1;

                const _A *a_it = a.data();
                _C *c_it = c.data();

                int incx = 1;
                int incy = 1;

                // Using daxpy is faster!!!
                for (int c_sli = 0; c_sli != c_slis; c_sli++) {
                    for (int c_row = 0; c_row != c_rows; c_row++) {
                        const _A *a_it2 = a_it;
                        const _B *b_it = b_end;
                        for (int b_sli = 0; b_sli != b_slis; b_sli++) {
                            for (int b_row = 0; b_row != b_rows; b_row++) {
                                for (int b_col = 0; b_col != b_cols; b_col++) {
                                    daxpy(&c_cols, b_it, a_it2, &incx, c_it, &incy);
                                    b_it--;
                                    a_it2++;
                                }
                                a_it2 += a_cols - b_cols;
                            }
                            a_it2 += (a_rows - b_rows) * a_cols;
                        }
                        c_it += c_cols;
                        a_it += a_cols;
                    }
                    a_it += (a_rows - c_rows) * a_cols;
                }

            }
            else {
                Array<_B> b_ = b;
                reverse(b_);

                int m = 0;
                ARRAY_EACH(shape_c, ind) {
                    V<Vi> r;
                    for (int i = 0; i < nd; i++) r.push_back({ind[i], ind[i]+shape_b[i]});
                    c[m] = Array<_C>::times(array_sub(a, r), b_).sum();
                    m++;
                }
            }
        }

        void run_plain_full(const Array<_A> &a, const Array<_B> &b) {
            if (!shape_c.empty()) shape_c.clear();
            for (int i = 0; i < nd; i++) shape_c.push_back(shape_a[i]+shape_b[i]-1);
            c = Array<_C>(shape_c, 0);

            if (nd == 2) {
                ERROR_REPORT("Not implemented yet!");
            }
            else if (nd == 3) {
                int a_slis = shape_a[0];
                int a_rows = shape_a[1];
                int a_cols = shape_a[2];

                int b_slis = shape_b[0];
                int b_rows = shape_b[1];
                int b_cols = shape_b[2];

                int c_slis = shape_c[0];
                int c_rows = shape_c[1];
                int c_cols = shape_c[2];

                int a_size = a.size();
                int b_size = b.size();
                int c_size = c.size();

                const _B *b_end = b.data() + b_size - 1;

                const _A *a_it = a.data();
                _C *c_it = c.data();

                int incx = 1;
                int incy = 1;

                // Using daxpy is faster!!!
                for (int a_sli = 0; a_sli != a_slis; a_sli++) {
                    for (int a_row = 0; a_row != a_rows; a_row++) {
                        _C *c_it2 = c_it;
                        const _B *b_it = b_end;
                        for (int b_sli = 0; b_sli != b_slis; b_sli++) {
                            for (int b_row = 0; b_row != b_rows; b_row++) {
                                for (int b_col = 0; b_col != b_cols; b_col++) {
                                    daxpy(&a_cols, b_it, a_it, &incx, c_it2, &incy);
                                    b_it--;
                                    c_it2++;
                                }
                                c_it2 += c_cols - b_cols;
                            }
                            c_it2 += (c_rows - b_rows) * c_cols;
                        }
                        c_it += c_cols;
                        a_it += a_cols;
                    }
                    c_it += (c_rows - a_rows) * c_cols;
                }
            }
            else {
                ERROR_REPORT("Not implemented yet!");
            }
        }

        template<typename _DA, typename _DB>
        void run_fft(const BasicArray<_A, _DA> &a, const BasicArray<_B, _DB> &b, int result_mode) {
            set_shape(a.shape(), b.shape());

            if (result_mode == CONVOLVE_VALID) run_fft_valid(a, b);
            else if (result_mode == CONVOLVE_SAME) run_fft_same(a, b);
            else if (result_mode == CONVOLVE_FULL) run_fft_full(a, b);
        }

        template<typename _DA, typename _DB>
        void run_fft_valid(const BasicArray<_A, _DA> &a, const BasicArray<_B, _DB> &b) {
            shape_c = shape_a;

            auto && aa = array_reshape(a, shape_a);
            auto && bb = indent_b(b);
            auto && cc = Array<_C>::convolve(aa, bb);

            if (nd == 2) {
                int rows = shape_a[0]-shape_b[0]+1;
                int cols = shape_a[1]-shape_b[1]+1;
                c = Array<_C>({rows, cols});
                int m = 0;
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        c[m] = cc(i, j);
                        m++;
                    }
                }
            }
            else {
                STD_ vector<STD_ vector<int>> range;
                for (int i = 0; i < nd; i++) {
                    range.push_back({0, shape_c[i]-shape_b[i]+1});
                }
                c = array_sub(cc, range);
            }
        }

        template<typename _DA, typename _DB>
        void run_fft_same(const BasicArray<_A, _DA> &a, const BasicArray<_B, _DB> &b) {
            if (!shape_c.empty()) shape_c.clear();
            for (int i = 0; i < nd; i++) shape_c.push_back(shape_a[i]+shape_b[i]-1);

            auto && aa = indent_a_full(a);
            auto && bb = indent_b(b);
            auto && cc = Array<_C>::convolve(aa, bb);

            STD_ vector<STD_ vector<int>> range;
            for (int i = 0; i < nd; i++) {
                range.push_back({shape_b[i]-1, shape_c[i]});
            }
            c = array_sub(cc, range);
        }

        template<typename _DA, typename _DB>
        void run_fft_full(const BasicArray<_A, _DA> &a, const BasicArray<_B, _DB> &b) {
            if (!shape_c.empty()) shape_c.clear();
            for (int i = 0; i < nd; i++) shape_c.push_back(shape_a[i]+shape_b[i]-1);

            auto && aa = indent_a_full(a);
            auto && bb = indent_b(b);
            auto && cc = Array<_C>::convolve(aa, bb);

            c = STD_ move(cc);

        }

    private:
        int nd;
        Shape shape_a, shape_b;
        Shape shape_c;

        void set_shape(const Shape &sa, const Shape &sb) {
            int nda = sa.size();
            int ndb = sb.size();
            nd = STD_ max(nda, ndb);

            if (!shape_a.empty()) shape_a.clear();
            for (int i = 0; i < nd-nda; i++) shape_a.push_back(1);
            for (auto && n : sa) shape_a.push_back(n);

            if (!shape_b.empty()) shape_b.clear();
            for (int i = 0; i < nd-ndb; i++) shape_b.push_back(1);
            for (auto && n : sb) shape_b.push_back(n);

        }

        template<typename T>
        T &reverse(T &in) {
            Shape shape = in.shape();
            int nd = in.dim();
            if (nd == 2) {
                // Accelerate when nd == 2
                int rows = in.shape(0);
                int cols = in.shape(1);
                int m = 0;
                for (int i = 0; i < rows; i++) {
                    if (i >= rows/2) break;
                    for (int j = 0; j < cols; j++) {
                        STD_ swap(in[m], in(rows-1-i, cols-1-j));
                        m++;
                    }
                }
            }
            else {
                Shape v(in.dim(), 0);
                auto symm = [&shape](Shape v) -> Shape {
                    for (int i = 0; i < v.size(); i++) {
                        v[i] = shape[i]-1-v[i];
                    }
                    return ::std::move(v);
                };
                ARRAY_EACH(shape, ind) {
                    if (ind[0] >= shape[0]/2) break;
                    ::std::swap(in(ind), in(symm(ind)));
                }
            }
            return in;
        }

        template<typename T>
        Array<_C> indent_a_full(const T &a) {
            Array<_C> aa(shape_c, 0);
            if (nd == 2) {
                // Accelerate when nd == 2
                int rows = shape_a[0];
                int cols = shape_a[1];
                int m = 0;
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        aa(shape_b[0]-1+i, shape_b[1]-1+j) = a(m);
                        m++;
                    }
                }
            }
            else {
                STD_ vector<STD_ vector<int>> range;
                for (int i = 0; i < nd; i++) {
                    range.push_back({shape_b[i]-1, shape_c[i]});
                }
                array_sub(aa, range) = a;
            }
            return STD_ move(aa);
        }

        template<typename T>
        Array<_C> indent_b(const T &b) {
            Array<_C> bb(shape_c, 0);
            if (nd == 2) {
                // Accelerate when nd == 2
                int rows = shape_b[0];
                int cols = shape_b[1];
                int m = 0;
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        bb(rows-1-i, cols-1-j) = b(m);
                        m++;
                    }
                }
            }
            else {
                Array<typename T::value_type> b_ = b;
                reverse(b_);
                STD_ vector<STD_ vector<int>> range;
                for (int i = 0; i < nd; i++) {
                    range.push_back({0, shape_b[i]});
                }
                array_sub(bb, range) = b_;
            }
            return STD_ move(bb);
        }

    };

    template<typename _C, typename _A, typename _B, typename _DA, typename _DB>
    inline Array<_C> convolve(const BasicArray<_A, _DA> &a, const BasicArray<_B, _DB> &b,
                              int result_mode=CONVOLVE_FULL, int compute_mode = CONVOLVE_FFT) {
        return STD_ move(Convolve<_A, _B, _C>(a, b, result_mode, compute_mode).c);
    }

    template<typename _A, typename _B, typename _DA, typename _DB>
    inline Arrayd convolve(const BasicArray<_A, _DA> &a, const BasicArray<_B, _DB> &b,
                           int result_mode=CONVOLVE_FULL, int compute_mode = CONVOLVE_FFT) {
        return STD_ move(Convolve<_A, _B, double>(a, b, result_mode, compute_mode).c);
    }

}} // namespace cppsci::signal

