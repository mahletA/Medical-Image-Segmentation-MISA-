clc; clear;
%%  MICO2D run on split slices %%

img = load_nii('D:\MAIA_Masters\Semester3_Spain\MISA\l1_preprocessing\braindata\t1_icbm_normal_1mm_pn0_rf0.nii');
shape = size(img.img);

