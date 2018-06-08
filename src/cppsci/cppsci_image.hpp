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

namespace cppsci { namespace image {

    enum {
        ROTATE_NEAREST,
        ROTATE_BILINEAR,
        ROTATE_BICUBIC
    };

    /**
     * Refer to: http://blog.csdn.net/fengbingchun/article/details/17713429
     */
    template<typename _A>
    Arrayd rotate_nearest(const Array<_A> &a, double degree) {
        using namespace ::std;

        if (degree == 0) return a;

        double alpha = degree / 180.0 * 3.1415927;
        int rows = a.shape(0);
        int cols = a.shape(1);

        double x = (cols-1)/2.0;
        double y = (rows-1)/2.0;

        double max_x = fabs(cos(alpha)*x) + fabs(sin(alpha)*y);
        double max_y = fabs(sin(alpha)*x) + fabs(cos(alpha)*y);

        int brows = int(ceil(max_y*2));
        int bcols = int(ceil(max_x*2));
        int bsz = brows*bcols;
        double bx = (bcols-1)/2.0;
        double by = (brows-1)/2.0;

//        cout << "bx:" << bx << " by:" << by << endl;

        Arrayd b({brows, bcols});
        Arrayd c({bsz, 2});
        double *c_it = c.data();
        for (int i = 0; i < brows; i++) {
            for (int j = 0; j < bcols; j++) {
                c_it[0] = -bx + j;
                c_it[1] =  by - i;
                c_it += 2;
            }
        }

        c_it = c.data();
        for (int i = 0; i < bsz; i++) {
            double cx = cos(-alpha)*c_it[0]-sin(-alpha)*c_it[1];
            double cy = sin(-alpha)*c_it[0]+cos(-alpha)*c_it[1];
            c_it[0] = cx;
            c_it[1] = cy;
            c_it += 2;
        }

        double *b_it = b.data();
        c_it = c.data();
        for (int i = 0; i < brows; i++) {
            for (int j = 0; j < bcols; j++) {
                int jj = int(c_it[0] - (-x+0.5) + 1);
                int ii = int((y-0.5) -  c_it[1] + 1);
                if (ii >= 0 && ii < rows && jj >= 0 && jj < cols) *b_it = a[ii*cols+jj];
                else *b_it = 0;
                c_it += 2;
                b_it++;
            }
        }

        return move(b);
    }

}} // end namespace cppsci::ndimage

