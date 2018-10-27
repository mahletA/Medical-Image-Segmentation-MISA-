% Comparison of Observed image, filtered and Bias corrected 
orig = load_nii('D:\MAIA_Masters\Semester3_Spain\MISA\l1_preprocessing\Lab1\braindata\t1_icbm_normal_1mm_pn5_rf20.nii');
slice = 90;
Orig = double(orig.img(:,:,slice))';
figure;imshow(Orig,[]);
Orig= (Orig- min(Orig(:)))/ (max(Orig(:))- min(Orig(:)) );
imwrite(Orig,'Orig.png')

filtered = load_nii('D:\MAIA_Masters\Semester3_Spain\MISA\l1_preprocessing\Lab1\Results\pn5_rf20_filtered.nii');
Filtered = double(filtered.img(:,:,slice))';
figure;imshow(Filtered,[]);
Filtered= (Filtered- min(Filtered(:)))/ (max(Filtered(:))- min(Filtered(:)) );
imwrite(Filtered,'Filtered.png')

bc_fil = load_nii('D:\MAIA_Masters\Semester3_Spain\MISA\l1_preprocessing\Lab1\Results\pn5_rf20_filtered_bc.nii');
Bc_fil = double(bc_fil.img(:,:,slice))';
figure;imshow(Bc_fil,[]);
Bc_fil= (Bc_fil- min(Bc_fil(:)))/ (max(Bc_fil(:))- min(Bc_fil(:)) );
imwrite(Bc_fil,'Bc_fil.png')
