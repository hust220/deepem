#include <fstream>
#include "unit_test.h"
#include "cppsci_traits.hpp"

TEST_CASE(test_traits)
{
    using namespace cppsci;

    {
        Vi a{1,2,3,4};
        Lf b{1,2,3,4};
        TEST_CHECK(eles_equal(a,b)==true);
    }

    {
        Vi a{1,2,3};
        Lf b{1,2,3,4};
        TEST_CHECK(eles_equal(a,b)==false);
    }

}



