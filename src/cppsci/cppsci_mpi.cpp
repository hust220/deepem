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

#include "cppsci_mpi.hpp"

#ifdef USEMPI

int mpi_rank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int mpi_size() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

MpiAlloc mpi_alloc(int size) {
    int mpirank = mpi_rank();
    int mpisize = mpi_size();

    MpiAlloc alloc;
    alloc.size = size;
    alloc.step = 1;
    alloc.bin = int(STD_ ceil(size/double(mpisize)));
    alloc.beg = mpirank*alloc.bin;
    int end_ = alloc.beg+alloc.bin;
    alloc.end = (end_>size?size:end_);

    return STD_ move(alloc);
}

MpiAlloc mpi_alloc(int size, int step) {
    int mpirank = mpi_rank();
    int mpisize = mpi_size();

    ::std::vector<int> v;
    for (int i = 0; i < size; i+=step) {
        v.push_back(i);
    }

    auto &&alloc = mpi_alloc(v.size());
    alloc.size = size;
    alloc.step = step;
    alloc.beg = v[alloc.beg];
    alloc.end = v[alloc.end];
    alloc.bin = alloc.bin*step;

    return STD_ move(alloc);
}

#endif // USEMPI


