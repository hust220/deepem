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

#include "../cppsci/cppsci.hpp"

namespace deepem {

using namespace ::cppsci;

/**
 * Head of mrc file.
 */
struct MrcHead {
#define MRCHEAD_DATA \
        ELT(int,    NC,       0,  1, "of Columns    (fastest changing in map)"                             ) SEP \
        ELT(int,    NR,       0,  2, "of Rows"                                                             ) SEP \
        ELT(int,    NS,       0,  3, "of Sections   (slowest changing in map)"                             ) SEP \
        ELT(int,    MODE,     2,  4, "Data type,0,1,2,3,4,5 we use 2(float)"                               ) SEP \
        ELT(int,    NCSTART,  0,  5, "Number of first COLUMN  in map"                                      ) SEP \
        ELT(int,    NRSTART,  0,  6, "Number of first ROW     in map"                                      ) SEP \
        ELT(int,    NSSTART,  0,  7, "Number of first SECTION in map"                                      ) SEP \
        ELT(int,    NX,       0,  8, "Number of intervals along X"                                         ) SEP \
        ELT(int,    NY,       0,  9, "Number of intervals along Y"                                         ) SEP \
        ELT(int,    NZ,       0,  10, "Number of intervals along Z"                                        ) SEP \
        ELT(float,  X_length, 0,  11, "Cell Dimensions (Angstroms)"                                        ) SEP \
        ELT(float,  Y_length, 0,  12, "\""                                                                 ) SEP \
        ELT(float,  Z_length, 0,  13, "\""                                                                 ) SEP \
        ELT(float,  Alpha,    0,  14, "Cell Angles     (Degrees)"                                          ) SEP \
        ELT(float,  Beta,     0,  15, "\""                                                                 ) SEP \
        ELT(float,  Gamma,    0,  16, "\""                                                                 ) SEP \
        ELT(int,    MAPC,     0,  17, "Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)"                 ) SEP \
        ELT(int,    MAPR,     0,  18, "Which axis corresponds to Rows   (1,2,3 for X,Y,Z)"                 ) SEP \
        ELT(int,    MAPS,     0,  19, "Which axis corresponds to Sects. (1,2,3 for X,Y,Z)"                 ) SEP \
        ELT(float,  AMIN,     0,  20, "Minimum density value"                                              ) SEP \
        ELT(float,  AMAX,     0,  21, "Maximum density value"                                              ) SEP \
        ELT(float,  AMEAN,    0,  22, "Mean    density value    (Average)"                                 ) SEP \
        ELT(int,    ISPG,     0,  23, "Space group number"                                                 ) SEP \
        ELT(int,    NSYMBT,   0,  24, "Number of bytes used for storing symmetry operators"                ) SEP \
        ELT(int,    LSKFLG,   0,  25, "Flag for skew transformation, =0 none, =1 if foll"                  ) SEP \
        ELT(int,    SKWMAT11, 0,  26, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT12, 0,  27, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT13, 0,  28, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT21, 0,  29, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT22, 0,  30, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT23, 0,  31, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT31, 0,  32, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT32, 0,  33, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT33, 0,  34, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWTRN_X, 0,  35, "Skew translation t if LSKFLG .ne. 0.Skew transformation is from"    ) SEP \
        ELT(int,    SKWTRN_Y, 0,  36, "standard orthogonal coordinate frame (as used for atoms) to"        ) SEP \
        ELT(int,    SKWTRN_Z, 0,  37, "orthogonal map frame, as Xo(map) = S * (Xo(atoms) - t)"             ) SEP \
        ELT(int,    NONE1,    0,  38, "(some of these are used by the MSUBSX routines"                     ) SEP \
        ELT(int,    NONE2,    0,  39, "in MAPBRICK, MAPCONT and FRODO)"                                    ) SEP \
        ELT(int,    NONE3,    0,  40, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE4,    0,  41, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE5,    0,  42, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE6,    0,  43, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE7,    0,  44, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE8,    0,  45, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE9,    0,  46, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE10,   0,  47, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE11,   0,  48, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE12,   0,  49, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE13,   0,  50, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE14,   0,  51, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE15,   0,  52, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    MAP,      0,  53, "Character string 'MAP ' to identify file type"                      ) SEP \
        ELT(int,    MACHST,   0,  54, "Machine stamp indicating the machine type which wrote file"         ) SEP \
        ELT(int,    ARMS,     0,  55, "Rms deviation of map from mean density"                             ) SEP \
        ELT(int,    NLABL,    0,  56, "Number of labels being used"                                        ) // end of macro
        
#define ELT(T,N,V,I,C) T N = V;
#define SEP
        MRCHEAD_DATA
#undef ELT
#undef SEP
        // 10  80 character text labels (ie. A4 format)
        int LABEL[200] = {0};
        bool is_ieee_le = true;

 };

/**
 * Push the MrcHead object to stream.
 */
STD_ ostream &operator <<(STD_ ostream& os, const MrcHead &head);

/**
 * Read mrc head from the stream.
 */
MrcHead mrc_head_read(STD_ istream &stream);

/**
 * Read data from the mrc file.
 */
template<typename _V = double>
Array<_V> mrc_read(Str filename) {
    STD_ ifstream ifile(filename.c_str());
    MrcHead &&head = mrc_head_read(ifile);
    int nx = head.NC;
    int ny = head.NR;
    int nz = head.NS;
    int size = nx*ny*nz;
    Array<_V> a({nz, ny, nx});
    ifile.seekg(1024);
    if (head.MODE == 1) {
        int *buf = new int[size];
        ifile.read((char *)buf, size*sizeof(int));
        for (int i = 0; i < size; i++) a[i] = _V(buf[i]);
        delete [] buf;
    }
    else if (head.MODE == 2) {
        float *buf = new float[size];
        ifile.read((char *)buf, size*sizeof(float));
        for (int i = 0; i < size; i++) a[i] = _V(buf[i]);
        delete [] buf;
    }
    ifile.close();
    return STD_ move(a);
}

/**
 * Read the sli-th slide from the mrc file.
 */
template<typename _V = double>
Array<_V> mrc_read(Str filename, int sli) {
    STD_ ifstream ifile(filename.c_str());
    MrcHead &&head = mrc_head_read(ifile);
    int nx = head.NC;
    int ny = head.NR;
    int nz = head.NS;
    int size = nx*ny;
    float temp;
    JN_INFO(nx);
    JN_INFO(ny);
    JN_INFO(nz);
    Array<_V> a({ny, nx});
    if (head.MODE == 1) {
        int *buf = new int[size];
        ifile.seekg(1024+sizeof(int)*size*sli);
        ifile.read((char *)buf, size*sizeof(int));
        if (head.is_ieee_le) for (int i = 0; i < size; i++) a[i] = _V(buf[i]);
        else for (int i = 0; i < size; i++) {
            ((int *)&temp)[0] = SWAP32(((int *)buf)[i]);
            a[i] = _V(temp);
        }
        delete [] buf;
    }
    else if (head.MODE == 2) {
        float *buf = new float[size];
        ifile.seekg(1024+sizeof(float)*size*sli);
        ifile.read((char *)buf, size*sizeof(float));
        if (head.is_ieee_le) for (int i = 0; i < size; i++) a[i] = _V(buf[i]);
        else for (int i = 0; i < size; i++) {
            ((int *)&temp)[0] = SWAP32(((int *)buf)[i]);
            a[i] = _V(temp);
        }
        delete [] buf;
    }
    ifile.close();
    return STD_ move(a);
}

} // end namespace deepem

