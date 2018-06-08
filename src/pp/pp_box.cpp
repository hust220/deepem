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

#include "pp_box.hpp"

namespace deepem {

using namespace ::cppsci;

Arrayi box_read(Str filename) {
    L<A<int, 4>> ls;
    Str line;
    Vs v;

    STD_ ifstream ifile(filename.c_str());
    while (ifile) {
        STD_ getline(ifile, line);
        string::tokenize(line, v, "\t ");
        if (v.size() == 4) {
            ls.push_back({CPPSCI_INT(v[0]), CPPSCI_INT(v[1]), CPPSCI_INT(v[2]), CPPSCI_INT(v[3])});
        }
    }
    ifile.close();

    int n = ls.size();
    Arrayi a({n, 4});
    int *it = a.data();
    for (auto && p : ls) {
        for (int i = 0; i < 4; i++) {
            *it = p[i];
            it++;
        }
    }

    return STD_ move(a);
}

} // end namespace deepem

