%% INITIALIZATION %%

function [alpha, Mu, Sigma] = Initialize (X,k)

% X is an N x f matrix
N = size(X,1);
alpha = 1/k* ones(1,k);

[~,Mu] = kmeans(X,k);
Sigma = cov(X);

disp ("******* INITIALIZATION DONE! *******")
end
