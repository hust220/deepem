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

class Par {
    bool debug_Par = false;
public:
    using Qs = Qs;

    M<Str, Qs> pars;
    M<Str, Str> intros;
    Qs globalPar;

    int argc;
    char **argv;

    Par() = default;

    Par(const Par &par) = default;

    Par(int argc, char **argv);

    Par(Str str);

    void read(int argc, char **argv);

    void read(Str par_file);

    friend STD_ ostream &operator <<(STD_ ostream &out, const Par &par);

    template<typename K>
    bool has(K &&k) const {
        return pars.find(k) != pars.end();
    }

    template<typename K, typename U, typename... V>
    bool has(K &&k, U &&u, V && ...rest) const {
        return has(k) || has(u, rest...);
    }

    template<typename T>
    void set(T &&v) const {}

    template<typename T, typename K, typename... V>
    void set(T &&v, K &&s, V && ...rest) const {
        if (s != "" && pars.find(s) != pars.end()) {
            if(debug_Par) std::cout<<__FILE__<<" "<<__LINE__<<" pars.at(s).size() = "<<pars.at(s).size()<<std::endl;
            v = parse<typename STD_ decay<T>::type>(pars.at(s)[0]);
        }
        else set(v, rest...);
    }

    template<typename T>
    void setv(T &&v) const {}

    template<typename T, typename _First, typename... _Rest>
    void setv(T &&t, _First &&first, _Rest && ...rest) const {
        if (first != "" && pars.find(first) != pars.end()) {
            auto && l = pars.at(first);
            t.resize(l.size());
            STD_ copy(l.begin(), l.end(), t.begin());
        }
        else {
            setv(t, rest...);
        }
    }

    STD_ list<Str> &keys_chain() const {
        static STD_ list<Str> chain;
        return chain;
    }

    Str get() const {
        STD_ ostringstream stream;
        stream << "jian::Par::get error! Didn't found parameters for keys:";
        for (auto && key : keys_chain()) stream << " " << key;
        throw stream.str();
    }

    template<typename K, typename... V>
    Str get(K &&s, V && ..._pars) const {
        if (pars.count(s)) {
            keys_chain().clear();
            return pars.at(s)[0];
        } else {
            keys_chain().push_back(s);
            return get(_pars...);
        }
    }

    Qs getv() const {
        STD_ ostringstream stream;
        stream << "jian::Par::getv error! Didn't found parameters for keys:";
        for (auto && key : keys_chain()) stream << " " << key;
        throw stream.str();
    }

    template<typename K, typename... V>
    Qs getv(K &&s, V && ..._pars) const {
        if (pars.count(s)) {
            keys_chain().clear();
            return pars.at(s);
        }
        else {
            keys_chain().push_back(s);
            return getv(_pars...);
        }
    }

    template<typename T>
    T parse(const Str &s) const {
        return lexical_cast<T>(s);
    }

};

} // namespace cppsci

