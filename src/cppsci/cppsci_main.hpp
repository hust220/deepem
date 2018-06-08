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
#include "cppsci_string.hpp"
#include "cppsci_par.hpp"

// #define CPPSCI_REG_COMPONENT(i, par)\
// int JN_PP_CAT(main_, i)(const ::cppsci::Par &);\
// ::cppsci::MainComponent JN_PP_CAT(component_, i)(JN_PP_STR(i), JN_PP_CAT(main_, i));\
// int JN_PP_CAT(main_, i)(const ::cppsci::Par &par)

namespace cppsci {

class RomePickerMain {
public:
    using func_t = STD_ function<int(const Par &)>;
    using func_map_t = M<Str, func_t>;

    func_map_t methods;

    func_t main_help = [](const Par &par) -> int {
        STD_ cout
        << "Wrong command, "<< string::join(' ', par.globalPar) << std::endl
        << "Possible commands : " << std::endl;
        for (auto && p : RomePickerMain::instance().methods) {
            // Vs v;
            // string::tokenize(p.first, v, "_");
            // STD_ cout << "    " << string::join(' ', v) << '\n';
            STD_ cout << "rome_picker " << p.first << std::endl;
        }
		return 1;
    };

    static RomePickerMain &instance() {
        static RomePickerMain main;
        return main;
    }

    static void run(int argc, char **argv) {
        Par par(argc, argv);
        instance().get_func(par.globalPar)(par);
    }

    int get_score(Str s, Qs p) {
        Vs v;
        string::tokenize(s, v, "_");
        if (v.size() != p.size()) return 0;
        auto it1 = v.begin();
        auto it2 = p.begin();
        int score = 0;
        for (; it1 != v.end() && it2 != p.end(); it1++, it2++) {
            if (*it1 == "*") continue;
            else if (*it1 == *it2) score++;
            else return 0;
        }
        return score;
    }

    func_t &get_func(Qs path) {
        using pair_t = STD_ pair<Str, int>;
        M<Str, int> scores;
        for (auto && p : methods) {
            scores[p.first] = get_score(p.first, path);
        }
        auto it = STD_ max_element(scores.begin(), scores.end(), [](const pair_t &p1, const pair_t &p2){
            return p1.second <= p2.second;
        });
        if (it->second == 0) return main_help;
        else {
            return methods[it->first];
        }
    }

};

class RomePickerComponent {
public:
    template<typename F>
    static void addComponent(const ::cppsci::Str &name, F &&f) {
        RomePickerMain::instance().methods[name] = f;
    }
    static void runComponent(int argc, char **argv) {
        RomePickerMain::run(argc, argv);
    }
};

} // namespace cppsci

