% MISA Lab 2 preprocessing
clear all, close all
T1_img = load_untouch_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_2\5\T1.nii');
Flair_img = load_untouch_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_2\5\T2_FLAIR.nii');
GT_img = load_untouch_nii('C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_2\5\LabelsForTesting.nii');

shape = size(T1_img.img);
%slice = 30;
for slice = 1:shape(3)
    
    T1_slice = double(T1_img.img(:,:,slice));
    Flair_slice = double(Flair_img.img(:,:,slice));
    GT_slice = double(GT_img.img(:,:,slice));
    
    % thresholding
    mask = GT_slice > 0;
    
    T1_slice1 = T1_slice .* mask;
    T1_slice = (T1_slice1-min(T1_slice(:)))/(max(T1_slice(:))-min(T1_slice(:)));
    T1_slice = uint8(round(T1_slice*255));
    Flair_slice1 = Flair_slice .* mask;
    Flair_slice = (Flair_slice1-min(Flair_slice(:)))/(max(Flair_slice(:))-min(Flair_slice(:)));
    Flair_slice = uint8(round(Flair_slice*255));
    
    %figure, imshow(T1_slice,[]);
    
    
    %figure, imshow(T1_slice);
    T1_img.img(:,:,slice) = T1_slice(:,:);
    save_untouch_nii(T1_img, 'C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_2\5\preprocessedT1.nii');
    Flair_img.img(:,:,slice) = Flair_slice(:,:);
    save_untouch_nii(Flair_img, 'C:\Users\dono_\Documents\1.MAIA Girona\MISA\Lab_2\5\preprocessedFlair.nii');
end

