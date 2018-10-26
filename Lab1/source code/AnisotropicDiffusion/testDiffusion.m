function testDiffusion ()
%% 2D for 1 slice
clear all, close all
num_iter = 100;
delta_t = 0.25;
kappa = 6.2;% 6.2 option 2 27 option 1
option = 2;

img_ref = load_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_1\Medical-Image-Segmentation-MISA-\braindata\t1_icbm_normal_1mm_pn0_rf0.nii');
img = load_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_1\Medical-Image-Segmentation-MISA-\braindata\t1_icbm_normal_1mm_pn5_rf0.nii');
shape = size(img.img);
slice = 76; % a sample slice

img_slice = double(img.img(:,:,slice))';
img_ref_slice = double(img_ref.img(:,:,slice))';
    if mean2(img_slice) == 0
        fprintf('error');
    else

        %blurring
        H = fspecial('disk',3);
        blurred = imfilter(img_slice,H,'replicate');

        %anisotropic diffusion
        ad = anisodiff(img_slice,num_iter,kappa,delta_t,option);
        ad1 = (ad-min(ad(:)))/(max(ad(:))-min(ad(:)));
        ad1 = uint8(round(ad1*255));
        figure, subplot 131, imshow(img_slice,[]), subplot 132, imshow(ad,[]), subplot 133, imshow(blurred,[])
%       imwrite(ad1,"C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_1\iso_smoothing_01.png")
            
        %quality of original image 
           
        [peaksnr_ori,snr_ori] = psnr(uint8(img_slice),uint8(img_ref_slice));

        [peaksnr,snr] = psnr(uint8(ad),uint8(img_ref_slice));
        mserror = immse(ad,img_ref_slice);
        fprintf('\n The Peak-SNR value of original image is %0.4f', peaksnr_ori);
        fprintf('\n The Peak-SNR value of filtered image is %0.4f', peaksnr);
        fprintf('\n The MSE value is %0.4f \n', mserror);

    end
    
%% 2D PSNR analysis 
clear all, close all
num_iter = 100;
delta_t = 0.25;

i = 1;

img_ref = load_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_1\Medical-Image-Segmentation-MISA-\braindata\t1_icbm_normal_1mm_pn0_rf0.nii');
img = load_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_1\Medical-Image-Segmentation-MISA-\braindata\t1_icbm_normal_1mm_pn5_rf0.nii');

slice = 76; % a sample slice

img_slice = double(img.img(:,:,slice))';
img_ref_slice = double(img_ref.img(:,:,slice))';
    if mean2(img_slice) == 0
        fprintf('error');
    else
        for kappa= 1:0.2:50
            %anisotropic diffusion
            ad1 = anisodiff(img_slice,num_iter,kappa,delta_t,1);
            ad2 = anisodiff(img_slice,num_iter,kappa,delta_t,2);    
            %quality
            if i ==1
                % quality of the original image 
                [peaksnr_ori,snr_ori] = psnr(uint8(img_slice),uint8(img_ref_slice));
            end
            
            [peaksnr1(i),snr] = psnr(uint8(ad1),uint8(img_ref_slice));
            [peaksnr2(i),snr] = psnr(uint8(ad2),uint8(img_ref_slice));
            i = i+1;
%           
        end
    end
 
    figure,p=plot(1:0.2:50,peaksnr1,1:0.2:50,peaksnr2)
    p(1).LineWidth = 1
    p(2).LineWidth = 1
    xlabel("k"), ylabel("PSNR (dB)")
    legend("Leclerc","Lorentz")
    grid

%% 3D

clear all, close all
num_iter = 100;
delta_t = 0.25;
kappa = 6.2; %0.9 for pn 1, 6.2 for pn 5
option = 2;

img = load_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_1\Medical-Image-Segmentation-MISA-\braindata\t1_icbm_normal_1mm_pn5_rf20.nii');
shape = size(img.img);

% for all volume
for slice = 1:shape(3)
    Img = double(img.img(:,:,slice));
    if mean2(Img) == 0
         
        continue;
        
    else
        ad = anisodiff(Img,num_iter,kappa,delta_t,option);
        ad = (ad-min(ad(:)))/(max(ad(:))-min(ad(:)));
        ad = uint8(round(ad*255));
        
        img.img(:,:,slice) = ad(:,:);
        save_nii(img, 'C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_1\Medical-Image-Segmentation-MISA-\results\pn5_rf20_filtered.nii');
       
    end
end


