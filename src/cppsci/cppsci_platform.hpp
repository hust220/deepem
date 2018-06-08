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

/* Clang/LLVM. ---------------------------------------------- */
#if defined(__clang__)
#  define CPPSCI_CLANG

/* Intel ICC/ICPC. ------------------------------------------ */
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#  define CPPSCI_ICC

/* GNU GCC/G++. --------------------------------------------- */
#elif defined(__GNUC__) || defined(__GNUG__)
#  define CPPSCI_GCC

/* Hewlett-Packard C/aC++. ---------------------------------- */
#elif defined(__HP_cc) || defined(__HP_aCC)
#  define CPPSCI_HPCC

/* IBM XL C/C++. -------------------------------------------- */
#elif defined(__IBMC__) || defined(__IBMCPP__)
#  define CPPSCI_IBMC

/* Microsoft Visual Studio. --------------------------------- */
#elif defined(_MSC_VER)
#  define CPPSCI_MSC

/* Portland Group PGCC/PGCPP. ------------------------------- */
#elif defined(__PGI)

/* Oracle Solaris Studio. ----------------------------------- */
#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)

#endif

#if defined(__APPLE__) && defined(__GNUC__)  
#  define CPPSCI_OS_MACX  
#elif defined(__MACOSX__)  
#  define CPPSCI_OS_MACX  
#elif defined(macintosh)  
#  define CPPSCI_OS_MAC9  
#elif defined(__CYGWIN__)  
#  define CPPSCI_OS_CYGWIN  
#elif defined(MSDOS) || defined(_MSDOS)  
#  define CPPSCI_OS_MSDOS  
#elif defined(__OS2__)  
#  if defined(__EMX__)  
#    define CPPSCI_OS_OS2EMX  
#  else  
#    define CPPSCI_OS_OS2  
#  endif  
#elif !defined(SAG_COM) && (defined(WIN64) || defined(_WIN64) || defined(__WIN64__))  
#  define CPPSCI_OS_WIN32  
#  define CPPSCI_OS_WIN64  
#  define CPPSCI_OS_WIN
#elif !defined(SAG_COM) && (defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))  
#  define CPPSCI_OS_WIN32  
#  define CPPSCI_OS_WIN
#elif defined(__MWERKS__) && defined(__INTEL__)  
#  define CPPSCI_OS_WIN32  
#  define CPPSCI_OS_WIN
#elif defined(__sun) || defined(sun)  
#  define CPPSCI_OS_SOLARIS  
#elif defined(hpux) || defined(__hpux)  
#  define CPPSCI_OS_HPUX  
#elif defined(__ultrix) || defined(ultrix)  
#  define CPPSCI_OS_ULTRIX  
#elif defined(sinix)  
#  define CPPSCI_OS_RELIANT  
#elif defined(__linux__) || defined(__linux)  
#  define CPPSCI_OS_LINUX  
#elif defined(__FreeBSD__)  
#  define CPPSCI_OS_FREEBSD  
#  define CPPSCI_OS_BSD4  
#elif defined(__NetBSD__)  
#  define CPPSCI_OS_NETBSD  
#  define CPPSCI_OS_BSD4  
#elif defined(__OpenBSD__)  
#  define CPPSCI_OS_OPENBSD  
#  define CPPSCI_OS_BSD4  
#elif defined(__bsdi__)  
#  define CPPSCI_OS_BSDI  
#  define CPPSCI_OS_BSD4  
#elif defined(__sgi)  
#  define CPPSCI_OS_IRIX  
#elif defined(__osf__)  
#  define CPPSCI_OS_OSF  
#elif defined(_AIX)  
#  define CPPSCI_OS_AIX  
#elif defined(__Lynx__)  
#  define CPPSCI_OS_LYNX  
#elif defined(__GNU_HURD__)  
#  define CPPSCI_OS_HURD  
#elif defined(__DGUX__)  
#  define CPPSCI_OS_DGUX  
#elif defined(__QNXNTO__)  
#  define CPPSCI_OS_QNX6  
#elif defined(__QNX__)  
#  define CPPSCI_OS_QNX  
#elif defined(_SEQUENT_)  
#  define CPPSCI_OS_DYNIX  
#elif defined(_SCO_DS)                   /* SCO OpenServer 5 + GCC */  
#  define CPPSCI_OS_SCO  
#elif defined(__USLC__)                  /* all SCO platforms + UDK or OUDK */  
#  define CPPSCI_OS_UNIXWARE  
#  define CPPSCI_OS_UNIXWARE7  
#elif defined(__svr4__) && defined(i386) /* Open UNIX 8 + GCC */  
#  define CPPSCI_OS_UNIXWARE  
#  define CPPSCI_OS_UNIXWARE7  
#else  
#  error "Qt has not been ported to this OS - talk to qt-bugs@trolltech.com"  
#endif  

#if defined(CPPSCI_OS_MAC9) || defined(CPPSCI_OS_MACX)  
#  define CPPSCI_OS_MAC  
#endif  

#if defined(CPPSCI_OS_MAC9) || defined(CPPSCI_OS_MSDOS) || defined(CPPSCI_OS_OS2) || defined(CPPSCI_OS_WIN32) || defined(CPPSCI_OS_WIN64)  
#  undef CPPSCI_OS_UNIX  
#elif !defined(CPPSCI_OS_UNIX)  
#  define CPPSCI_OS_UNIX  
#endif
