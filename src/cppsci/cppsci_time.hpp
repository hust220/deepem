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

#include "cppsci_traits.hpp"
#include "cppsci_mpi.hpp"

namespace cppsci {

inline double dtime() {
    return std::chrono::duration_cast<std::chrono::duration<double>>
        (std::chrono::system_clock::now().time_since_epoch()).count();
}

inline double &last_time() {
    //thread_local static double s_last_time;
    //return s_last_time;
    assert(omp_get_max_threads()<=256);
    static double s_last_time[256];
    return s_last_time[omp_get_thread_num()];
}

inline void tic() {
    last_time() = dtime();
}

inline void toc(std::ostream &stream = std::cout) {
    MPI_LOG << "Elapsed time is " << dtime() - last_time() << " seconds." << std::endl;
}

struct Timer {
    double last_time;

    Timer() {
        last_time = dtime();
    }

    double toc() {
        double t = last_time;
        last_time = dtime();
        return last_time - t;
    }
};


}

