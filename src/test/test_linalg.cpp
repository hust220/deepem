#include "unit_test.h"
#include "cppsci_array.hpp"
#include "cppsci_linalg.hpp"

TEST_CASE(test_linalg)
{
    using namespace cppsci;

    ::std::cout << "linalg" << ::std::endl;
    Arrayd a1({2,3},{1,2,3,4,5,6});
    linalg::svd rt(a1);
//    ::std::cout << rt.s << '\n' << rt.u << '\n' << rt.vt << ::std::endl;

    auto && b1 = linalg::pinv(a1);
    auto && c1 = linalg::dot(a1, linalg::dot(b1, a1));
    TEST_CHECK((b1-c1).all([](double d){return d<1e-5;}));

    Arrayd a2({2,2},{1,2,3,4});
    linalg::eig eigrt(a2);
    ::std::cout << eigrt.a << ::std::endl;
    ::std::cout << eigrt.wr << ::std::endl;
    ::std::cout << eigrt.wi << ::std::endl;
    ::std::cout << eigrt.vl << ::std::endl;
    ::std::cout << eigrt.vr << ::std::endl;
//    ::std::cout << "a1\n" << a1 << ::std::endl;
//    ::std::cout << "b1\n" << b1 << ::std::endl;
//    ::std::cout << "b1*a1*b1\n" << linalg::dot(b1, linalg::dot(a1, b1)) << ::std::endl;
//    ::std::cout << "a1*b1*a1\n" << linalg::dot(a1, linalg::dot(b1, a1)) << ::std::endl;

}



