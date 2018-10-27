clc; close
% Script to calculate similarity measurements between bias corrected and
% original image

ref = imread('./Results/n0_rf0_slice90.png');
im1 = imread('./Results/n0_rf20_slice90.png');
im2 = imread('./Results/n0_rf40_slice90.png');

q = 5;

%Dice = [];
MI1 = [];
MI2 = [];
SD1 = [];
SD2 = [];
for i=1:q
    filename1 = strcat('./Results/noise0_bias20_q', int2str(i),'.png');
    img1 = imread(filename1);
    filename2 = strcat('./Results/noise0_bias40_q', int2str(i),'.png');
    img2 = imread(filename2);
    MI1 = [MI1, MI_GG(ref,img1)];
    MI2 = [MI2, MI_GG(ref,img2)];
    SD1 = [SD1, SD(ref, img1)];
    SD2 = [SD2, SD(ref, img2)];
end

plot(1:q, MI1,1:q, MI2,'LineWidth', 2);
xlabel('q')
ylabel('MI')
set(gca,'XTick',(1:1:5))
grid on;
legend('High bias field', 'Low bias field');
saveas(gcf,'./Results/MI.png');


plot(1:q, SD1,1:q, SD2,'LineWidth', 2);
xlabel('q')
ylabel('SSD')
grid on
set(gca,'XTick',(1:1:5))
legend('High bias field', 'Low bias field');
saveas(gcf,'./Results/SSD.png');
