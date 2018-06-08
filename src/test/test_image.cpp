#include <fstream>
#include "unit_test.h"
#include "cppsci_image.hpp"

TEST_CASE(test_image)
{
    using namespace cppsci;

    STD_ cout << "test_image" << STD_ endl;

    Arrayd a1({1, 10}, {1,2,3,4,5,6,7,8,9,10});
    image::rotate_nearest(a1, 20).print();
    //TEST_CHECK(approx_eq(multiply<double>(a13, b13), c13));

}



