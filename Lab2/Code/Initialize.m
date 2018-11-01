%% INITIALIZATION %%

function [alpha, Mu, Sigma] = Initialize (X,K)

% X is an N x f matrix
N = size(X,1);
f = size(X,2);
alpha = 1/K* ones(1,K);

[idx,Mu] = kmeans(X,K);

Sigma = cov(X);
% Sigma = zeros(K*f,f);
% for k = 1:K
%      if k ==1
%         Sigma(1:f,:) =  cov(X(idx==k)); 
%      else
%         Sigma(f*(k-1)+1:f*k,:) = cov(X(idx==k)); 
%      end   
% end

disp ("******* INITIALIZATION DONE! *******")
end
