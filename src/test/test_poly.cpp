#include "unit_test.h"
#include "cppsci_array.hpp"
#include "cppsci_poly.hpp"

TEST_CASE(test_poly)
{
    using namespace cppsci;

    Arrayd a1({3},{1,2,3});
    Arrayd b1 = polyvander(a1, 2);
    Arrayd c1({3,3},{1,1,1,1,2,4,1,3,9});
    TEST_CHECK(Arrayb::eq(b1,c1).all());

//    b1.printFull();
    Arrayd d1 = b1.sum(1);
//    d1.printFull();
    Arrayd e1({3}, {3,7,13});
    TEST_CHECK(Arrayb::eq(d1,e1).all());

    Arrayd a2({5},{1,2,3,4,5});
    Arrayd b2({5, 1}, {3,5,7,6,5});
    auto && rt = polyfit(a2, b2, 3);

    Arrayd a3({4},{1,2,3,4});
    Arrayd b3({4},{1,2});
    auto && c3 = polyval<double>(a3, b3);
}



