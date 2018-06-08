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

#include "pp_mrc.hpp"

namespace deepem {

using namespace ::cppsci;

STD_ ostream &operator <<(STD_ ostream& os, const MrcHead &head) {
    os << "Mrc Head {"
#define ELT(T,N,V,I,C) <<" "<< #N <<":"<<head.N
#define SEP
    MRCHEAD_DATA
#undef ELT
#undef SEP
    << " }" << std::endl;
	return os;
}

MrcHead mrc_head_read(STD_ istream &stream) {
    MrcHead head;
    stream.read((char*)(&head.NC),256*sizeof(int));
    if (head.MODE > 100) {
        head.is_ieee_le = false;
        int *p = &head.NC;
        for (int i = 0; i < 256; i++) {
            p[i] = SWAP32(p[i]);
        }
    }
    JN_INFO(head);
    return STD_ move(head);
}

} // end namespace deepem


