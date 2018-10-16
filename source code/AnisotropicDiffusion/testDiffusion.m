function testDiffusion ()
num_iter = 100;
delta_t = 0.25;
kappa = 5;
option = 1;

im1 = rgb2gray(imread('D:\MAIA_Masters\Semester3_Spain\MISA\l1_preprocessing\source code\MICO_v0\MICO_v0\MICO_2D\brainweb64.tif'));
%im1 = (dicomread('slice_3D.dcm'));

%blurring
H = fspecial('disk',3);
blurred = imfilter(im1,H,'replicate');

ad = anisodiff(im1,num_iter,kappa,delta_t,option);
figure, subplot 131, imshow(im1), subplot 132, imshow(ad,[]), subplot 133, imshow(blurred,[])



