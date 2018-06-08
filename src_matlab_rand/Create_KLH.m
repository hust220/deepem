clear;
OUTPUT_CNN_NAME='KLH.rand';
boxsize=272;

dim_x=2048;
dim_y=2048;

name_length=2;
name_prefix='';
mic_path='../test/KLHdata/mic/';

num_positive1=800;
num_negative1=800;
positive1_box_path='../test/KLHdata/positive272/';
negative1_box_path='../test/KLHdata/negative272/';
positive1_mic_start_num=1;
positive1_mic_end_num=50;
negative1_mic_start_num=1;
negative1_mic_end_num=40;

do_train_again=0;
num_positive2=0;
num_negative2=0;
positive2_box_path='';
negative2_box_path='';
positive2_mic_start_num=0;
positive2_mic_end_num=0;
negative2_mic_start_num=0;
negative2_mic_end_num=0;

rotation_angel=90;
rotation_n=360/rotation_angel;

num_p_test=300;
num_n_test=300;

%FL_kernelsize=51;
FL_kernelsize=51;
%TL_kernalsize=21;
TL_kernalsize=21;
%FIL_kernalsize=10;
FIL_kernalsize=10;

SL_poolingsize=3;
FOL_poolingsize=2;
SIL_poolingsize=2;

FL_feature_map=6;
TL_feature_map=12;
FIL_feature_map=12;

create_CNN;
