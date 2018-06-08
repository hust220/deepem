#include "unit_test.h"
#include "cppsci_array.hpp"

TEST_CASE(test_fft)
{
    using namespace cppsci;

    Arrayd data({2,2,2},{0,1,2,3,4,5,6,7});
    Arrayd f({2,2,2},{7,6,5,4,3,2,1,0});
    fftshift3(data.data(), data.data(), 2);
    TEST_CHECK(Arrayd::minus(data,f).all([](double n){return n==0;}));
    f.selfFFTshift3();
    TEST_CHECK(Arrayd::plus(data,f).all([](double n){return n==7;}));
}



