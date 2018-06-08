#include "unit_test.h"
#include "cppsci_array.hpp"
#include "cppsci_utils.hpp"

TEST_CASE(test_utils)
{
    using namespace cppsci;

    auto arr = Arrayd::linspace(1,4,4);
    Arrayd arr2({4},{1,2,3,4});
    TEST_CHECK((arr-arr2).all([](double d){return d<1e-5;}));

    Arrayd a2({3,3,3},{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27});
    ArrayWindowRoller roller(a2.shape(), {2, 2, 2}, {1, 1, 1});
    Arrayd sum({2,2,2}, 0);
    do {
        sum += array_sub(a2, roller.range);
    } while (roller.roll());
    Arrayd a3({2,2,2},{60,68,84,92,132,140,156,164});
    TEST_CHECK((sum-a3).all([](double d){return d<1e-5;}));
}



