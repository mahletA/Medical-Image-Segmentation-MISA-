%% EXPECTATION %%

function [W] = E_Step(X, alpha, Mu, Sigma)
N = size(X,1);
f = size(X,2);
K = size(alpha,2);

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

%disp ("******* EXPECTATION DONE! *******")
end