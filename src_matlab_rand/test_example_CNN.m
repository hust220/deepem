function a=test_example_CNN()

load mnist_uint8;

addpath(genpath('DeepLearnToolbox'));

train_x = double(reshape(train_x',28,28,60000))/255;
train_x = train_x(:,:,1:1000);

test_x = double(reshape(test_x',28,28,10000))/255;
test_x = test_x(:,:,1:1000);

train_y = double(train_y');
train_y = train_y(:,1:1000);

test_y = double(test_y');
test_y = test_y(:,1:1000);

%% ex1 Train a 6c-2s-12c-2s Convolutional neural network 
%will run 1 epoch in about 200 second and get around 11% error. 
%With 100 epochs you'll get around 1.2% error
rand('state',0);
cnn.layers = {
    struct('type', 'i') %input layer
    struct('type', 'c', 'outputmaps', 6, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %sub sampling layer
    struct('type', 'c', 'outputmaps', 12, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %subsampling layer
};

disp('cnnsetup');
cnn = cnnsetup(cnn, train_x, train_y);

opts.alpha = 1;
opts.batchsize = 50;
opts.numepochs = 1;

disp('cnntrain');
cnn = cnntrain(cnn, train_x, train_y, opts);

disp('cnntest');
[er, bad] = cnntest(cnn, test_x, test_y);

%plot mean squared error
disp('plot');
figure; plot(cnn.rL);

disp(['error: ' er])
% assert(er<0.12, 'Too big error');