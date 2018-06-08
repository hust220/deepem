#include <sstream>
#include "unit_test.h"
#include "cppsci_serial.hpp"

TEST_CASE(test_serial)
{
    using namespace cppsci;

    STD_ cout << "test_serial" << STD_ endl;

    {
        STD_ stringstream stream;
        int a = 5;
        int b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(a == b);
    }

    {
        STD_ stringstream stream;
        float a = 1.5;
        float b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(a == b);
    }

    {
        STD_ stringstream stream;
        double a = 1.5;
        double b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(a == b);
    }

    {
        STD_ stringstream stream;
        char a = '#';
        char b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(a == b);
    }

    {
        STD_ stringstream stream;
        Str a = "hello";
        Str b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(a == b);
    }

    {
        STD_ stringstream stream;
        Vf a = {1.1, 1.3, 5.2, 4.7};
        Vf b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(eles_equal(a,b));
    }

    {
        STD_ stringstream stream;
        Str a = "hello";
        Str b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(a == b);
    }

    {
        STD_ stringstream stream;
        Arrayd a({2,3},{1,2,3,4,5,6});
        Arrayd b;
        serialize(stream, a);
        parse(stream, b);
        TEST_CHECK(approx_eq(a, b));
    }

}



