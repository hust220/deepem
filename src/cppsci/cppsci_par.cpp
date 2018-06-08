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
#include "cppsci_par.hpp"

namespace cppsci {

Par::Par(int argc, char **argv) {
    read(argc, argv);
}

Par::Par(Str str) {
    read(str);
}

void Par::read(Str par_file) {
    STD_ ifstream ifile(par_file.c_str());
	if (!ifile) {
		std::cout << "Note: The parameter file '" << par_file << "' is not found!!!" << std::endl;
	}
    Str line;
    Vs tokens;
    while (ifile) {
        STD_ getline(ifile, line);
        string::tokenize(line, tokens, " ");
        if (!tokens.empty()) {
            pars[tokens[0]] = Qs{};
            for (auto it2 = STD_ next(tokens.begin()); it2 != tokens.end(); it2++) {
                pars[tokens[0]].push_back(*it2);
                if(debug_Par) std::cout<<__FILE__<<" "<<__LINE__<<" "<<tokens[0]<<" "<<*it2<<std::endl;
            }
            if (tokens.size()==1) { // no value for this option
                pars[tokens[0]].push_back(std::string(""));
            }
            if (tokens[0] == "par") {
                for (auto && p : pars["par"]) read(p);
            }
        }
    }
    ifile.close();

}

void Par::read(int argc, char **argv) {
    this->argc = argc;
    this->argv = argv;
    Str key;
    Qs values;
    int n = 0;
    for (int i = 1; i < argc; i++)
    {
        if (Str(argv[i]).size() >=2 && argv[i][0] == '-' && STD_ isalpha(argv[i][1]))
        {
            if (n != 0) {
                pars[key] = values;
                if (key == "par") {
                    for (auto && p : values) {
                        read(p);
                    }
                }
            } else {
                globalPar = values;
                // pars["global"] = values;
            }
            Str str(argv[i]);
            key = str.substr(1, str.size() - 1);
            if(debug_Par) std::cout<<__FILE__<<" "<<__LINE__<<" key : "<<key<<std::endl;
            values.clear();
            n++;
        } else {
            if(debug_Par) std::cout<<__FILE__<<" "<<__LINE__<<" "<<argv[i]<<std::endl;
            values.push_back(argv[i]);
        }
    }
    if (n != 0) {
        pars[key] = values;
        if (key == "par") {
            for (auto && p : values) read(p);
        }
    } else {
        globalPar = values;
        // pars["global"] = values;
    }
}

STD_ ostream &operator <<(STD_ ostream &out, const Par &par) {
    for (auto && i : par.pars) {
        out << i.first;
        for (auto && j : i.second) {
            out << ' ' << j;
        }
        out << STD_ endl;
    }
    return out;
}

} // namespace cppsci

