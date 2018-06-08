#include "unit_test.h"
#include "cppsci_stats.hpp"

TEST_CASE(test_stats)
{
    using namespace cppsci;

    ::std::cout << __FUNCTION__ << ::std::endl;
    for (int i = 0; i < 10; i++)
        ::std::cout << ::std::fixed << ::std::setprecision(12) << norm_cdf(i) << ::std::endl;
}



