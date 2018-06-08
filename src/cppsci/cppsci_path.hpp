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

#include "cppsci_string.hpp"
#include "cppsci_platform.hpp"

namespace cppsci {

inline STD_ pair<Str, Str> path_splitext(Str name) {
    STD_ smatch match;
    STD_ regex_match(name, match, STD_ regex("^(.*)(\\.[^.]+)$"));
    return {match[1], match[2]};
}

inline Str path_basename(Str name) {
    STD_ smatch match;
    STD_ regex_match(name, match, STD_ regex("^.*/([^/]+)$"));
    return match[1];
}

inline bool path_exists(const Str& path) {
#if defined(CPPSCI_OS_WIN)
    struct _stat info;
    return _stat(path.c_str(), &info) == 0;
#else 
    struct stat info;
    return stat(path.c_str(), &info) == 0;
#endif
}

inline bool path_make(const Str& path) {
#if defined(CPPSCI_OS_WIN)
    int ret = _mkdir(path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
#endif
    if (ret == 0)
        return true;

    switch (ret) {
        case ENOENT:
            // parent didn't exist, try to create it
            {
                int pos = path.find_last_of('/');
                if (pos == Str::npos)
#if defined(CPPSCI_OS_WIN)
                    pos = path.find_last_of('\\');
                if (pos == Str::npos)
#endif
                    return false;
                if (!path_make( path.substr(0, pos) ))
                    return false;
            }
            // now, try to create again
#if defined(CPPSCI_OS_WIN)
            return 0 == _mkdir(path.c_str());
#else 
            return 0 == mkdir(path.c_str(), mode);
#endif

        case EEXIST:
            // done!
            return path_exists(path);

        default:
            return false;
    }
}

} // namespace cppsci


