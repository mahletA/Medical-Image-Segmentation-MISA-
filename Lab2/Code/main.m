%% Loading data 

clear;
clc;

T1 = load_untouch_nii('../Data/1/T1.nii');
T2 = load_untouch_nii('../Data/1/T2_FLAIR.nii');

x1 = [];
x2 = [];
N = size(T1.img,1) * size(T1.img,2);
for i = 1: size(T1.img,3)
    x1 = vertcat(x1,reshape(double(T1.img(:,:,i)),N,1));
    x2 = vertcat(x2,reshape(double(T2.img(:,:,i)),N,1));
end
  
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
while diff > 100
    iteration = iteration + 1;
    
    [alpha_new, Mu_new, Sigma_new] = M_Step (X, W);
    [W] = E_Step(X, alpha_new, Mu_new, Sigma_new);
    
    logL_new = log_likelihood(X,alpha_new,Mu_new,Sigma_new);
    diff = abs(logL-logL_new);
    logL = logL_new;
    sprintf('%i    -------    %0.4f',iteration,diff);
end

