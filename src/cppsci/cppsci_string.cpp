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

#include "cppsci_string.hpp"

namespace cppsci {
namespace string {

void tokenize(const Str &str, Vs &tokens, const Str &delimiters) {
    tokens.clear();
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void tokenize(const Str &str, Vs &tokens, const Str &delimiters, const Str &temp) {
    tokens.clear();
    using pair_t = ::std::pair<Str::size_type, Str::size_type>;
    ::std::vector<pair_t> vec;
    Str::size_type first_i, first_j, second_i, second_j;
    int expected = 0;
    for (Str::size_type i = 0; i < str.size(); i++) {
        int flag = 0;
        Str::size_type j;
        for (j = 0; j < temp.size(); j++) {
            if (str[i] == temp[j]) {
                if (j % 2 == 0 && expected == 0) { flag = 1; break;
                } else if (j % 2 == 1 && expected == 1) { flag = 2; break; }
            }
        }
        if (flag == 1) {
            first_i = i; first_j = j; expected = 1;
        } else if (flag == 2 && j - first_j == 1) {
            second_i = i; second_j = j; expected = 0;
            vec.push_back(::std::make_pair(first_i, second_i));
        }
    }
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (::std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p){
        return pos != Str::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
        while (::std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p){
            return pos != Str::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
}

} // namespace string
} // namespace cppsci

