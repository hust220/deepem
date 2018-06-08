#include "unit_test.h"
#include "cppsci_optimize.hpp"

TEST_CASE(test_optimize)
{
    using namespace cppsci;

    auto f = [](double x) ->double {
        return (x-2)*x*(x+2)*(x+2);
    };
    auto r = minimize_scalar(f);
    ::std::cout << "xmin: " << r.xmin << ::std::endl;
    ::std::cout << "iter: " << r.iter << ::std::endl;
    ::std::cout << "fval: " << r.fval << ::std::endl;
    ::std::cout << "funcalls: " << r.funcalls << ::std::endl;
}



