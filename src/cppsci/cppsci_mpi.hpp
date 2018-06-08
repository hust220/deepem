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

#include "cppsci_mkl.hpp"
#include "cppsci_macros.hpp"

//#define USEMPI

#ifdef USEMPI

#include <vector>
#include <mpi.h>
#define MPI_INITIALIZE(argc, argv) MPI_Init(argc, argv)
#define MPI_FINALIZE MPI_Finalize()
#define MPI_IS_ROOT (mpi_rank() == 0)

#ifdef MPI_INFO_ALL
# define MPI_LOG ::std::cout << "Node " << mpi_rank() << ": "
#else
# define MPI_LOG if (MPI_IS_ROOT) ::std::cout
#endif

/**
 * Get the rank of the process.
 */
int mpi_rank();

/**
 * Get the number of all the processes.
 */
int mpi_size();

/**
 * MPI allocation
 */
struct MpiAlloc {
    int beg;
    int end;
    int bin;
    int size;
    int step;
};

/**
 * Allocate 0-size to processes.
 */
MpiAlloc mpi_alloc(int size);

/**
 * Allocate 0-size to processes with step.
 */
MpiAlloc mpi_alloc(int size, int step);

#else // USEMPI

#define MPI_INITIALIZE(argc, argv)
#define MPI_FINALIZE
#define MPI_LOG ::std::cout

#endif // USEMPI


