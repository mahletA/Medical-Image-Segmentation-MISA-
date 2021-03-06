
%% Loading data 

clear all;
clc;
profile on;
T1 = load_untouch_nii('../Data/5/T1.nii');
T2 = load_untouch_nii('../Data/5/T2_Flair.nii');

ground_truth = load_untouch_nii('../Data/5/LabelsForTesting.nii');

x1 = [];
x2 = [];

x1_original = [];
x2_original = [];
GT = [];

N = size(T1.img,1) * size(T1.img,2);
for i = 1: size(T1.img,3)
    x1_original = vertcat(x1_original,reshape(double(T1.img(:,:,i)),N,1));
    x2_original = vertcat(x2_original,reshape(double(T2.img(:,:,i)),N,1));
    
    GT = vertcat(GT,reshape(double(ground_truth.img(:,:,i)),N,1));

end

%Remove the backgroud pixels
mask = GT > 0;
mask2 = (GT == 0) * -50;
x1 = (x1_original .* mask) + mask2;
x2 = (x2_original .* mask) + mask2;
x1(x1 == -50 ) = [];
x2(x2 == -50 ) = [];

X = horzcat(x1, x2);


%% Initialize 
K=3; % for 3 regions, i.e. CSF, GM, WM
[alpha, Mu, Sigma1] = Initialize (X,K);

Sigma = [];
for i=1:K
    Sigma = vertcat(Sigma,Sigma1); % Initialize all sigmas the same way (vertically stack)
end

iteration = 0; 

N = size(X,1);
f = size(X,2);

sprintf('Iteration ------- Log Likelihood diff.')
diff = 1000;
logL = log_likelihood(X,alpha,Mu,Sigma);

while diff > 0.1
    iteration = iteration + 1;
    
    %[alpha_new, Mu_new, Sigma_new] = M_Step (X, W);
    %[W] = E_Step(X, alpha_new, Mu_new, Sigma_new);
    %logL_new = log_likelihood(X,alpha_new,Mu_new,Sigma_new);
    
    
    %% EXPECTATION %%

    W = zeros(N,K);

    for k =1:K
        if k == 1
            W(:,k) = mvnpdf(X ,Mu(k,:), Sigma(1:f,:))*alpha(k);
        else
            W(:,k)= mvnpdf(X, Mu(k,:), Sigma((f*(k-1))+1: f*k,:))*alpha(k);
        end
    end

    for i=1:N
        W(i,:) = W(i,:)/ sum(W(i,:));
    end 
    
    

    
     %% MAXIMIZATION %%
    
    Net = sum(W);

    % UPDATE alpha
    alpha_new = Net / N ;

    % UPDATE MU
    Mu_new = zeros(K, f);
    for k = 1:K
        Mu_new(k,:) = sum(W(:,k).* X)/Net(k);
    end
    
    % UPDATE SIGMA
    Sigma_new = zeros(K*f, f);
    for k = 1:K
        
         temp =  (W(:,k) .* (X - Mu_new(k,:)))' * (X - Mu_new(k,:));
         
         if k ==1
            Sigma_new(1:f,:) =  temp/Net(k); 
         else
            Sigma_new(f*(k-1)+1:f*k,:) = temp/Net(k);
         end
    end
    %Sigma_new
    
    %% Log likelihood
    temp = zeros(N,1);
    for k =1:K
        if k == 1
            temp =  temp + mvnpdf(X ,Mu_new(k,:), Sigma_new(1:f,:))*alpha_new(k);
        else
            temp = temp + mvnpdf(X, Mu_new(k,:), Sigma_new((f*(k-1))+1: f*k,:))*alpha_new(k);
        end
    end

    logL_new = sum(log(temp));
    
    diff = abs(logL-logL_new);
    
    % Update variables
    logL = logL_new;
    Sigma = Sigma_new;
    alpha = alpha_new;
    Mu = Mu_new;
    
    sprintf('%i    -------    %0.4f',iteration,diff)
    
    
end

%% TESTING
% Reconstruct segmented results
[~, I]= (max(W'));
final_seg = zeros(size(x1_original));
N_orig =  size(x1_original,1);

%Reintroduce the background
index = 1;
for i=1:N_orig
    if mask(i) == 0
        final_seg(i) = 0;
    else
        final_seg(i) = I(index);
        index = index + 1;
    end
end

%test the dice of a slice with ground truth
final_seg = reshape(final_seg, size(T1.img));

seg = final_seg(:,:,25)';

% Tweak segmentation labels to match ground truth
seg2 = seg; seg2(seg == 2) = 3; seg2(seg == 3) = 2;
%seg2(seg ==3) =1; seg2(seg == 2) = 3; seg2(seg ==1) = 2;

truth = double(ground_truth.img(:,:,25))';
% T1 = load_untouch_nii('../Data/5/preprocessedT1.nii');
% T2 = load_untouch_nii('../Data/5/preprocessedFlair.nii');
% t1_im = T1.img(:,:,25)';
% t2_im = T2.img(:,:,25)';
figure;imshow(truth,[])
figure; imshow(seg2,[])
figure;imshow(seg,[])
%figure; imshow(t1_im,[])
%figure; imshow(t2_im,[])


dice(truth, seg2)


%% Dice for the whole volume
dice_coeff = [];
for slice=1:size(final_seg,3)
    seg = final_seg(:,:,slice)';

    % Tweak segmentation labels to match ground truth
    seg2 = seg; seg2(seg == 2) = 3; seg2(seg == 3) = 2;
    truth = double(ground_truth.img(:,:,slice))';
    if size(dice(truth, seg2),1) == 3
        dice_coeff = horzcat(dice_coeff,dice(truth, seg2));
    end
    

end
mean_dice = mean(dice_coeff,2)

