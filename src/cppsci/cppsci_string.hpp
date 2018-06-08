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

namespace cppsci {
namespace string {

void tokenize(const Str &str, Vs &tokens, const Str &delimiters = " ");
void tokenize(const Str &str, Vs &tokens, const Str &delimiters, const Str &temp);

template<typename... Strs_>
::std::string merge(Strs_ && ...strs) {
    ::std::stringstream stream;
    stream_push(stream, strs...);
    return stream.str();
}

template<typename... Pars_>
void disp(Pars_ && ...pars) {
    stream_push(::std::cout, pars...);
}

template<typename... Pars_>
::std::string format(std::string fmt, Pars_ && ...pars) {
    int count = snprintf(NULL, 0, fmt.c_str(), pars...);
    ::std::string buf;
    buf.resize(count);
    sprintf(&(buf[0]), fmt.c_str(), pars...);
    return std::move(buf);
}

template<typename _Interval, typename _Ls>
static Str join(_Interval && interval, _Ls && ls) {
    STD_ stringstream stream;
    int i = 0;
    for (const auto & s : ls) {
        if (i != 0) stream << interval;
        stream << s;
        i++;
    }
    return stream.str();
}

} // namespace string
} // namespace cppsci

#ifndef DIE
#  define DIE(...) do {\
    ::std::cerr << ::cppsci::string::merge(__VA_ARGS__, " in ", __FILE__, ":", __LINE__) << ::std::endl;\
    ::std::abort();\
} while(0)
#endif

#ifndef ERROR_REPORT
#  define ERROR_REPORT(what) DIE(what)
#endif

#ifndef ERROR_CHECK
#  define ERROR_CHECK(condition,what) if(condition){ERROR_REPORT(what);}
#endif


