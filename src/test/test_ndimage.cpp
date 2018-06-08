#include "unit_test.h"
#include "cppsci_ndimage.hpp"

TEST_CASE(test_ndimage)
{
    using namespace cppsci;

    Arrayd a = Arrayd::range(12);
    a.reshape({4,3});
    auto && b = ndimage::map_coordinates(a, Arrayd({2,2},{0.5,2,0.5,1}),1);
    Arrayd c({2},{2,7});
    TEST_CHECK((b-c).all([](double n){return n<1e-5;}));
}



