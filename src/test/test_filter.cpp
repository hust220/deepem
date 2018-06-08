#include "unit_test.h"
#include "cppsci_array.hpp"
#include "cppsci_filter.hpp"

TEST_CASE(test_filter)
{
    using namespace cppsci;

    Arrayd a1({4}, {1,2,3,4});
    ArrayLines lines(a1.shape(), 0, 3, 3);
    lines.read_array(a1);
    Arrayd a2({10}, {3,2,1,1,2,3,4,4,3,2});
    TEST_CHECK((lines-a2).all([](double n){return n<1e-5;}));

    Arrayd b1(a1.shape());
    lines.write_array(b1);
    TEST_CHECK((b1-a1).all([](double n){return n<1e-5;}));

    Arrayd c1 = cppsci::gaussian_filter(a1, 0.5);

    Arrayd d1(a1.shape());
    correlate1d(a1, d1, Arrayd({5},{0.1,0.2,0.3,0.2,0.1}));

}



