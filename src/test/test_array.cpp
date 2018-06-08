#include <fstream>
#include "unit_test.h"
#include "cppsci_array.hpp"

TEST_CASE(test_array)
{
    using namespace cppsci;

    Arrayd arr({2, 2}, {0, 1, 2, 3});
    TEST_CHECK(arr[0] == 0);
    TEST_CHECK(arr(0) == 0);
    TEST_CHECK(arr(1, 1) == 3);
    TEST_CHECK(arr({1, 1}) == 3);
    TEST_CHECK(arr.trace() == 3);

    Arrayd arr2({2, 2, 3}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11});
    TEST_CHECK(arr2[1] == 1);
    TEST_CHECK(arr2(6) == 6);
    TEST_CHECK(arr2(1, 1, 1) == 10);

    auto && sub2 = array_sub(arr2, V<Vi>{{0,2},{0,2},{1,3}});
    TEST_CHECK(sub2(1, 1, 1) == 11);
    TEST_CHECK(sub2({1, 1, 1}) == 11);
    TEST_CHECK(sub2.trace() == 12);

    sub2.selfTimes(2).selfPlus(5).selfMinus(5);
    TEST_CHECK(sub2(1, 1, 1) == 22);
    TEST_CHECK(arr2(1, 1, 2) == 22);
    sub2 *= 0.5;
    TEST_CHECK(sub2(1, 1, 1) == 11);
    sub2.selfDivide(0.1);
    TEST_CHECK(sub2(1, 1, 1) == 110);
    sub2 /= 10;
    TEST_CHECK(sub2(1, 1, 1) == 11);
    TEST_CHECK(sub2.argmax() == 7);
    TEST_CHECK(sub2.max() == 11);
    TEST_CHECK(sub2.argmin() == 0);
    TEST_CHECK(sub2.min() == 1);
    TEST_CHECK(sub2.sum() == 48);

    Arrayd arr3({2, 2}, {5, 4, 3, 2});
    arr3.selfMinus(arr);
    TEST_CHECK(arr3(1, 1) == -1);
    TEST_CHECK(arr3.sum() == 8);
    arr3.selfTimes(arr);
    TEST_CHECK(arr3(1, 1) == -3);
    TEST_CHECK(arr3.max() == 3);

    auto && arr4 = Arrayd::plus(arr3, arr);
    TEST_CHECK(arr4.trace() == 0);

    Arrayd arr5({79}, 0);
    TEST_CHECK(arr5.sum() == 0);

    //------------------------------------------------------
    // Test sort
    Arrayd arr6({15}, {3.1, 2.2,-1.8, -5.3, 7.543, 8.7, 3, 2, 9, 0, -1, -5, 4, -2, 3});
    ::std::sort(arr6.begin(), arr6.end());

    //------------------------------------------------------
    // Test complex times
    Arraycd a7({2}, {::std::complex<double>(1,2), ::std::complex<double>(2,3)});
    Arraycd b7({2}, {::std::complex<double>(3,4), ::std::complex<double>(5,6)});
//    ::std::cout << Arraycd::times(a7, b7) << ::std::endl;
//    ::std::cout << a7 * b7 << ::std::endl;

    //----------------------------------------------------
    // Test meanv, sumv, argminv, minv, argmaxv, maxv
    Arrayd a8({2, 3},{1,2,3,4,5,6});

    Arrayd a8_mean({1,3},{2.5,3.5,4.5});
    TEST_CHECK((meanv(a8)-a8_mean).all([](double a){return a<1e-5;}));

    Arrayd a8_sum({1,3},{5,7,9});
    TEST_CHECK((sumv(a8)-a8_sum).all([](double a){return a<1e-5;}));

    Arrayd a8_argmin({1,3},{0,0,0});
    TEST_CHECK((argminv(a8)-a8_argmin).all([](double a){return a<1e-5;}));

    Arrayd a8_min({1,3},{1,2,3});
    TEST_CHECK((minv(a8)-a8_min).all([](double a){return a<1e-5;}));

    Arrayd a8_argmax({1,3},{1,1,1});
    TEST_CHECK((argmaxv(a8)-a8_argmax).all([](double a){return a<1e-5;}));

    Arrayd a8_max({1,3},{4,6,6});
    TEST_CHECK((maxv(a8)-a8_max).all([](double a){return a<1e-5;}));

    //-----------------------------------------------
    // Test repmat, expand, flipdim, flipud, fliplr, flipall
    Arrayd a9({2,2,2},{1,2,3,4,5,6,7,8});

    auto && b9 = repmat(a9, {2,2});
    TEST_CHECK(b9(3,2,0) == 5);

    auto && b92 = expand(a9, {2,2});
    TEST_CHECK(b92(3,2,0) == 7);

    auto && c9 = flipdim(a9, 2);
    Arrayd c9_flip({2,2,2},{2,1,4,3,6,5,8,7});
    TEST_CHECK(approx_eq(c9,c9_flip));

    auto && d9 = flipud(a9);
    Arrayd d9_flip({2,2,2},{5,6,7,8,1,2,3,4});
    TEST_CHECK(approx_eq(d9,d9_flip));

    auto && e9 = fliplr(a9);
    Arrayd e9_flip({2,2,2},{3,4,1,2,7,8,5,6});
    TEST_CHECK(approx_eq(e9,e9_flip));

    auto && f9 = flipall(a9);
    Arrayd f9_flip({2,2,2},{8,7,6,5,4,3,2,1});
    TEST_CHECK(approx_eq(f9,f9_flip));

    //-----------------------------------------------
    // Test subarray
    Arrayd a10({3,3},{1,2,3,4,5,6,7,8,9});
    Arrayd b10({3,1},{3,6,9});
    array_sub(a10, V<Vi>{{}, {2}}).identify("b10");
    TEST_CHECK(approx_eq(array_sub(a10, V<Vi>{{},{2}}),b10));
    Arrayd c10({2,3},{1,2,3,4,5,6});
    array_sub(a10, V<Vi>{{0,2}, {}}).identify("c10");
    TEST_CHECK(approx_eq(array_sub(a10, V<Vi>{{0,2},{}}),c10));
    Arrayd d10({1,3},{7,8,9});
    array_sub(a10, V<Vi>{{2}}).identify("d10");
    TEST_CHECK(approx_eq(array_sub(a10, V<Vi>{{2}}),d10));

    //----------------------------------------------------
    // Test operators
    Arrayd a11({2,2},{1,2,3,4});
    std::cout << "operators:" << std::endl;
    std::cout << a11 - 0.5 << std::endl;

    //----------------------------------------------------
    //Test read and write
    Arrayd a12({2,2},{1,3,5,9});
    STD_ string fname = STD_ tmpnam(nullptr);
    STD_ ofstream ofile(fname);
    array_write(a12, ofile);
    ofile.close();

    STD_ ifstream ifile(fname);
    auto && b12 = array_read<int>(ifile);
    ifile.close();
    TEST_CHECK(approx_eq(a12, b12));

    Arrayd a13({3,2}, {1,2,3,4,5,6});
    Arrayd b13({2,2}, {1,2,3,4});
    Arrayd c13({3,2}, {7, 10, 15, 22, 23, 34});
    Arrayd d13({2,3}, {1,3,5,2,4,6});
    TEST_CHECK(approx_eq(multiply<double>(a13, b13), c13));
    TEST_CHECK(approx_eq(multiply<double>(d13, b13, CblasTrans, CblasNoTrans), c13));

    // -----------------------------------------------------
    // Test ReshapeArray
    Arrayd a14({4}, {1,2,3,4});
    auto && b14 = array_reshape(a14, {2,2});
    Arrayd c14({2,2}, {1,2,3,4});
    TEST_CHECK(approx_eq(b14, c14));

    // -----------------------------------------------------
    // Test mapstd
    Arrayd a15({2, 4}, {1,2,3,4,5,6,7,8});
    Arrayd b15({2, 4}, {-1.1619, -0.3873, 0.3873, 1.1619, -1.1619, -0.3873, 0.3873, 1.1619}); 
    TEST_CHECK(approx_eq(mapstd(a15), b15));

    // -----------------------------------------------------
    // Test flipdim
    Arrayd a16({2, 4}, {1,2,3,4,5,6,7,8});
    Arrayd b16({2, 4}, {5,6,7,8,1,2,3,4}); 
    TEST_CHECK(approx_eq(flipdim(a16, 0), b16));

}



