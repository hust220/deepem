addpath(genpath('DeepLearnToolbox'));

positive1=read_train_data(name_prefix,name_length,positive1_mic_start_num,positive1_mic_end_num,mic_path,positive1_box_path,boxsize,dim_x,dim_y);
negative1=read_train_data(name_prefix,name_length,negative1_mic_start_num,negative1_mic_end_num,mic_path,negative1_box_path,boxsize,dim_x,dim_y);

if(do_train_again)
    positive2=read_train_data(name_prefix,name_length,positive2_mic_start_num,positive2_mic_end_num,mic_path,positive2_box_path,boxsize,dim_x,dim_y);
    negative2=read_train_data(name_prefix,name_length,negative2_mic_start_num,negative2_mic_end_num,mic_path,negative2_box_path,boxsize,dim_x,dim_y);
end

jn_infoa(positive1, 'positive1');
jn_infoa(negative1, 'negative1');

ran1=randperm(size(positive1,3));
for i=1:num_positive1
    for j=1:rotation_n
        a=imrotate(positive1(:,:,ran1(i)),rotation_angel*(j-1),'nearest');
        train_x(rotation_n*(i-1)+j,:)=a(:);
    end
end

if(do_train_again)
    ran2=randperm(size(positive2,3));
    for i=1:num_positive2
        for j=1:rotation_n
            a=imrotate(positive2(:,:,ran2(i)),rotation_angel*(j-1),'nearest');
            train_x(rotation_n*(i-1)+j+num_positive1*rotation_n,:)=a(:);
        end
    end
end

ran3=randperm(size(negative1,3));
for i=1:num_negative1
    for j=1:rotation_n
        a=imrotate(negative1(:,:,ran3(i)),rotation_angel*(j-1),'nearest');
        train_x(rotation_n*(i-1)+j+(num_positive2*do_train_again+num_positive1)*rotation_n,:)=a(:);
    end
end

if(do_train_again)
    ran4=randperm(size(negative2,3));
    for i=1:num_negative2
        for j=1:rotation_n
            a=imrotate(negative2(:,:,ran4(i)),rotation_angel*(j-1),'nearest');
            train_x(rotation_n*(i-1)+j+(num_positive2*do_train_again+num_positive1+num_negative1)*rotation_n,:)=a(:);
        end
    end
end
jn_infoa(train_x, 'train_x');
train_x=mapstd(train_x);

train_y(1:(rotation_n*(num_positive1+num_positive2*do_train_again)),1)=1;
train_y(rotation_n*(num_positive2*do_train_again+num_positive1)+1:rotation_n*(num_positive1+num_negative1+(num_positive2+num_negative2)*do_train_again),1)=0;

for i=1:num_p_test
    a=positive1(:,:,ran1(i+num_positive1));
    test_x(i,:)=a(:);
    test_y(i,1)=1;
end

for i=1:num_n_test
    a=negative1(:,:,ran3(i+num_negative1));
    test_x(i+num_p_test,:)=a(:);
    test_y(i+num_p_test,1)=0;
end

test_x=mapstd(test_x);

train_x = double(reshape(train_x',boxsize,boxsize,rotation_n*(num_positive1+num_negative1+do_train_again*(num_positive2+num_negative2))));
test_x = double(reshape(test_x',boxsize,boxsize,num_p_test+num_n_test));
train_y = double(train_y');
test_y = double(test_y');

rand('state',0)

cnn.layers = {
    struct('type', 'i')
    struct('type', 'c', 'outputmaps', FL_feature_map, 'kernelsize', FL_kernelsize)
    struct('type', 's', 'scale',SL_poolingsize)
    struct('type', 'c', 'outputmaps', TL_feature_map, 'kernelsize',TL_kernalsize)
    struct('type', 's', 'scale',FOL_poolingsize)
    struct('type', 'c', 'outputmaps', FIL_feature_map,'kernelsize',FIL_kernalsize)
    struct('type', 's', 'scale',SIL_poolingsize)
    };

opts.alpha = 1;
opts.batchsize = 50;
opts.numepochs = 20;

jn_infoa(train_x, 'train_x');
jn_infoa(train_y, 'train_y');

cnn = cnnsetup(cnn, train_x, train_y);
cnn = cnntrain(cnn, train_x, train_y, opts);

save(OUTPUT_CNN_NAME,'cnn');
% result=cnntest_m(cnn,test_x); % testing dataset
% 
% ri=0;
% for i=1:num_p_test
%     if(result(i)>0.5) % set the threhold as 0.5
%         ri=ri+1;
%     end
% end
% 
% for i=1:num_n_test
%     if(result(i+num_p_test)<=0.5)
%         ri=ri+1;
%     end
% end
% 
% p=ri/(num_p_test+num_n_test);
