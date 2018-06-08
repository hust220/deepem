#include "unit_test.h"
#include "cppsci_array.hpp"
#include "cppsci_filter.hpp"
#include "cppsci_morph.hpp"

TEST_CASE(test_morph)
{
    using namespace cppsci;

    Arrayi a1({3,5},{0,1,2,2,3,1,5,2,3,4,2,3,3,4,0});
    Arrayi b1({3,5},{0,1,2,3,2,1,2,3,2,1,2,3,2,1,0});
    auto && c1 = distance_transform_cdt(a1);
    TEST_CHECK((c1-b1).all([](int n){return n==0;}));

    Arrayi a2({3,3, 5},{
            1,1,2,2,3,
            1,5,2,3,4,
            2,3,3,4,0,

            0,1,2,2,3,
            1,5,2,3,4,
            2,3,3,4,0, 

            1,1,2,2,3,
            1,5,2,3,4,
            2,3,3,4,0});
    Arrayi b2({3,3, 5},{
            1,2,3,3,2,
            2,3,3,2,1,
            3,3,2,1,0,

            0,1,2,3,2,
            1,2,3,2,1,
            2,3,2,1,0, 

            1,2,3,3,2,
            2,3,3,2,1,
            3,3,2,1,0});
    auto && c2 = distance_transform_cdt(a2);
    ::std::cout << "c2\n" << c2 << ::std::endl;
    TEST_CHECK((c2-b2).all([](int n){return n==0;}));
}



