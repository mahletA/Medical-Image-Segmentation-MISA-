%% Loading data 

clear;
clc;
profile on;
T1 = load_untouch_nii('../Data/4/preprocessedT1.nii');
T2 = load_untouch_nii('../Data/4/preprocessedFlair.nii');

x1 = [];
x2 = [];

x1_original = [];
x2_original = [];

N = size(T1.img,1) * size(T1.img,2);
for i = 1: size(T1.img,3)
    x1_original = vertcat(x1_original,reshape(double(T1.img(:,:,i)),N,1));
    x2_original = vertcat(x2_original,reshape(double(T2.img(:,:,i)),N,1));
end
x1 = x1_original;
x2 = x2_original;
x1(x1 <= 6) = [];
x2(x2 <= 6) = [];

X = horzcat(x1, x2);


%% Initialize 
k=3; % for 3 regions, i.e. CSF, GM, WM
[alpha, Mu, Sigma1] = Initialize (X,k);

Sigma = [];
for i=1:k
    Sigma = vertcat(Sigma,Sigma1); % Initialize all sigmas the same way (vertically stack)
end

%% Epectation
[W] = E_Step(X, alpha, Mu, Sigma);

%% Maximization
logL = log_likelihood(X,alpha,Mu,Sigma);
iteration = 0; 

sprintf('Iteration ------- Log Likelihood diff.');
diff = 1000;

N = size(X,1);
f = size(X,2);
while diff > 100
    iteration = iteration + 1;
    
    %[alpha_new, Mu_new, Sigma_new] = M_Step (X, W);
    %[W] = E_Step(X, alpha_new, Mu_new, Sigma_new);
    %logL_new = log_likelihood(X,alpha_new,Mu_new,Sigma_new);
    %% MAXIMIZATION %%
    K = size(W,2);
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

%        temp = zeros(1,f);
%         for i= 1:N
%             temp = temp + W(i,k) .* ((X(i,:) - Mu_new(k,:))' * (X(i,:) - Mu_new(k,:)));
%         end
         temp =  W(i,k) .* ((X(:,:) - Mu_new(k,:))' * (X(:,:) - Mu_new(k,:)));
         if k ==1
            Sigma_new(1:f,:) =  temp/Net(k); 
         else
            Sigma_new(f*(k-1)+1:f*k,:) = temp/Net(k);
         end
    end
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
    
    
    %% Log likelihood
    temp = zeros(N,1);
    for k =1:K
        if k == 1
            temp =  temp +mvnpdf(X ,Mu(k,:), Sigma(1:f,:))*alpha(k);
        else
            temp = temp + mvnpdf(X, Mu(k,:), Sigma((f*(k-1))+1: f*k,:))*alpha(k);
        end
    end

    logL_new = sum(log(temp));
    
    diff = abs(logL-logL_new);
    logL = logL_new;
    sprintf('%i    -------    %0.4f',iteration,diff);
    
    
end


% Reconstruct segmented results
[~, I]= (max(W'));
final_seg = zeros(size(x1_original))
N_orig =  size(x1_original,1);

index = 1;
for i=1:N_orig
    if x1_original(i) == 0
        final_seg(i) = 0
    else
        final_seg(i) = I(index)
        index = index + 1
    end
end
    
final_seg = reshape(final_seg, size(T1.img));
imshow(final_seg(:,:,25)',[]);

