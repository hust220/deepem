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

#include <cstdio>
#include <iostream>
#include <cstdlib>

#include "../src/cppsci/cppsci.hpp"
#include "../src/pp/pp_train.hpp"
#include "../src/pp/pp_test.hpp"

int rome_picker_predict(const ::cppsci::Par &par)
{
    if (par.has("h", "help")) {
        ::deepem::pp_test_help();
    }
    else {
        double t1 = ::cppsci::dtime();

        ::deepem::pp_test(par);

        double t2 = ::cppsci::dtime();

        MPI_LOG << "Partical picking test costs : " << t2-t1 << " seconds." << STD_ endl;
    }

    return 0;
}

int rome_picker_train(const ::cppsci::Par &par)
{
    if (par.has("h", "help")) {
        ::deepem::pp_train_help();
    }
    else {
        double t1 = ::cppsci::dtime();
        
        ::deepem::pp_train(par);
        
        double t2 = ::cppsci::dtime();
        
        MPI_LOG << "Partical picking training costs : " << t2-t1 << " seconds." << STD_ endl;
    }
    
    return 0;
}

int rome_picker_test(const ::cppsci::Par &par)
{
    using namespace cppsci;
    auto && a = array_read<double>(par.get("a"));
    auto && b = array_read<double>(par.get("b"));
    auto &&c = multiply<double>(a, b);
    std::cout << c << std::endl;
    std::cout << cppsci::string::format("%.10f", std::pow(0.123456789, 3)) << std::endl;
    return 0;
}

int main(int argc, char **argv)
{
    MPI_INITIALIZE(&argc, &argv);
    
    cppsci::RomePickerComponent::addComponent("train", rome_picker_train);
    cppsci::RomePickerComponent::addComponent("predict", rome_picker_predict);
    cppsci::RomePickerComponent::addComponent("test", rome_picker_test);
    
    cppsci::RomePickerComponent::runComponent(argc, argv);
    
    MPI_FINALIZE;
}
