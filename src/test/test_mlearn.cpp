#include "unit_test.h"
#include "cppsci_mlearn.hpp"

TEST_CASE(test_mlearn)
{
    using namespace cppsci;
    using namespace cppsci::mlearn;

    std::cout << "test_mlearn" << std::endl;

    Arrayd a({2,2,2},{1,2,3,4,5,6,7,8});
    TEST_CHECK(a(1,1,0) == 7);

    auto && b = repmat(a, {2,2});
    auto && c = sigm(a);
    TEST_CHECK(b(3,2,0) == 5);
    TEST_CHECK(c(0,1,0) == 1.0/(1.0+std::exp(-3.0)));

    CNN cnn;
    cnn.layers = {
        cnn_input_layer(),
        cnn_conv_layer(6, 5),
        cnn_samp_layer(2),
        cnn_conv_layer(12, 5),
        cnn_samp_layer(2)
    };

    auto && train_x = array_read<double>("train_x.txt");
    train_x.selfDivide(255);
    auto && train_y = array_read<double>("train_y.txt");
    auto && test_x = array_read<double>("test_x.txt");
    test_x.selfDivide(255);
    auto && test_y = array_read<double>("test_y.txt");

    JN_INFOA(train_x);
    JN_INFOA(train_y);
    cnn_setup(cnn, train_x, train_y);
    CNNOpts opts{1, 50, 1};
    cnn_train(cnn, train_x, train_y, opts);

    auto r = cnn_test(cnn, test_x, test_y);

    STD_ cout << "err: " << r.first << STD_ endl;

}



