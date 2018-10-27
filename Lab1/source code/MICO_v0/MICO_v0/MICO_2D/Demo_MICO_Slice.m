% This Matlab file demomstrates the method for simultaneous segmentation and bias field correction
% in Chunming Li et al's paper:
%    "Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation",
%     Magnetic Resonance Imaging, vol. 32 (7), pp. 913-923, 2014
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://imagecomputing.org/~cmli/

clc;close all;clear all;
iterNum = 200;
N_region=3;  q=1;


% Load 3D data and select only one slice
img = load_nii('D:\MAIA_Masters\Semester3_Spain\MISA\l1_preprocessing\Lab1\braindata\t1_icbm_normal_1mm_pn0_rf20.nii');
shape = size(img.img);
slice = 90;
Img = double(img.img(:,:,slice))';

%load ROI
A=255;
Img_original = Img;
[nrow,ncol] = size(Img);n = nrow*ncol;

ROI = (Img>20); ROI = double(ROI);

tic

Bas=getBasisOrder3(nrow,ncol);
N_bas=size(Bas,3);
for ii=1:N_bas
    ImgG{ii} = Img.*Bas(:,:,ii).*ROI;
    for jj=ii:N_bas
        GGT{ii,jj} = Bas(:,:,ii).*Bas(:,:,jj).*ROI;
        GGT{jj,ii} = GGT{ii,jj} ;
    end
end


energy_MICO = zeros(3,iterNum);

b=ones(size(Img));
for ini_num = 1:1
    C=rand(3,1);
    C=C*A;
    M=rand(nrow,ncol,3);
    a=sum(M,3);
    for k = 1 : N_region
        M(:,:,k)=M(:,:,k)./a;
    end
    
    [e_max,N_max] = max(M,[], 3);
    for kk=1:size(M,3)
        M(:,:,kk) = (N_max == kk);
    end
    
    M_old = M; chg=10000;
    energy_MICO(ini_num,1) = get_energy(Img,b,C,M,ROI,q);
    
    
    for n = 2:iterNum
        pause(0.1)
        
        [M, b, C]=  MICO(Img,q,ROI,M,C,b,Bas,GGT,ImgG,1, 1);
        energy_MICO(ini_num,n) = get_energy(Img,b,C,M,ROI,q);
        
        figure(2),
        if(mod(n,1) == 0)
            PC=zeros(size(Img));
            for k = 1 : N_region
                PC=PC+C(k)*M(:,:,k);
            end
            subplot(241),imshow(uint8(Img)),title('original')
            subplot(242),imshow(PC.*ROI,[]); colormap(gray);
            iterNums=['segmentation: ',num2str(n), ' iterations'];
            title(iterNums);
            subplot(243),imshow(b.*ROI,[]),title('bias field')
            img_bc = Img./b;  % bias field corrected image
            subplot(244),imshow(uint8(img_bc.*ROI),[]),title('bias corrected')
            subplot(2,4,[5 6 7 8]),plot(energy_MICO(ini_num,:))
            xlabel('iteration number');
            ylabel('energy');
            pause(0.1)
        end
    end
end

[M,C]=sortMemC(M,C);
seg=zeros(size(Img));
for k = 1 : N_region
    seg=seg+k*M(:,:,k);   % label the k-th region 
end
figure;
subplot(141),imshow(Img,[]),title('Original image');
subplot(142),imshow(seg.*ROI,[]),title('Segmentation result');
subplot(143),imshow(b.*ROI,[]),title('bias field')
subplot(144),imshow(uint8(img_bc.*ROI),[]),title('bias corrected')

img_bc = (img_bc - min(img_bc(:)))/(max(img_bc(:))- min(img_bc(:)));
imwrite(img_bc.*ROI, './Results/noise0_bias20_q1.png');

% seg = (seg - min(seg(:)))/(max(seg(:)) - min(seg(:)));
% imwrite(seg.*ROI, './Results/seg_bias40_q5.png');
