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

#ifndef PI
#  define PI 3.14159265358979323846
#endif // PI

#ifdef JN_SHOW_INFO

#  define JN_INFO(x) MPI_LOG << #x << ": " << (x) << ::std::endl
#  define JN_INFOS(x) MPI_LOG << x << ::std::endl
#  define JN_INFOV(x) MPI_LOG << #x << ": "; for (auto && n : (x)) MPI_LOG << n << ' '; MPI_LOG << ::std::endl
#  define JN_INFOA(x) (x).identify(#x);
#  define JN_RUN(x) x;

#  define JN_INFOV2(x, name) MPI_LOG << name << ": " << (x) << ::std::endl
#  define JN_INFOA2(x, name) (x).identify(name);

#else

#  define JN_INFO(x)
#  define JN_INFOS(x)
#  define JN_INFOV(x)
#  define JN_INFOA(x)
#  define JN_RUN(x)

#  define JN_INFOV2(x, name)
#  define JN_INFOA2(x, name)

#endif // JN_SHOW_INFO

#  define STD_ ::std::
#  define CPPSCI_ ::cppsci::

#define SWAP32(x) ( (((x)&0x000000FF)<<24) | (((x)&0x0000FF00)<<8) | (((x)&0x00FF0000)>>8) | (((x)&0xFF000000)>>24) )


