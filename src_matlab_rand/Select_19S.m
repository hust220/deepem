clear;
addpath(genpath('DeepLearnToolbox'));
load('fcnn32_rot160.mat');
data_path='19Sdata/';

boxsize=160;

start_mic_num=30180;
end_mic_num=30180;

dim_x=1855;
dim_y=1919;
scan_step=20;
range1=100;
range2=40;

min_std=3.4;
max_std=3.7;

name_length=5;
name_prefix='image_';

rotation_angel=90;
rotation_n=360/rotation_angel;

select_CNN;