clear;
addpath(genpath('DeepLearnToolbox'));
load('KLH.acc.mat');
data_path='../test/KLHdata/';

boxsize=272;

start_mic_num=57;
end_mic_num=57;

dim_x=2048;
dim_y=2048;
scan_step=20;
range1=70;
range2=40;

min_std=22;
max_std=34;

name_length=2;%for the other dataset length=5
name_prefix='';%for the othe dataset 'image_'

rotation_angel=90;
rotation_n=360/rotation_angel;

select_CNN;
