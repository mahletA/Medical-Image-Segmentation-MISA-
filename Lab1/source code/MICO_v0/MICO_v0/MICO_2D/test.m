clc; close

fixed1 = imread('./Results/n0_rf20_slice90.png');
im1 = imread('./Results/noise0_bias20_q1.png');

disp('Noise 0 Bias 20')
SD(fixed1, im1)
MI_GG(fixed1, im1)

fixed2 = imread('./Results/n0_rf40_slice90.png');
im2 = imread('./Results/noise0_bias40_q3.png');

disp('Noise 0 Bias 40')
SD(fixed1, im1)
MI_GG(fixed1, im1)
