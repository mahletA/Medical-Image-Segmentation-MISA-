%% MAXIMIZATION %%

function [alpha_new, Mu_new, Sigma_new] = M_Step (X, W)
N = size(X,1);
f = size(X,2);
K = size(W,2);

Net = sum(W);

%% UPDATE alpha
alpha_new = Net / N ;

%% UPDATE MU
Mu_new = zeros(K, f);
for k = 1:K
    Mu_new(k,:) = sum(W(:,k).* X)/Net(k);
end
    
%% UPDATE SIGMA
Sigma_new = zeros(K*f, f);
for k = 1:K
    
    temp = zeros(1,f);
    for i= 1:N
        temp = temp + W(i,k) .* ((X(i,:) - Mu_new(k,:))' * (X(i,:) - Mu_new(k,:)));
    end
   
     if k ==1
        Sigma_new(1:f,:) =  temp/Net(k); 
     else
        Sigma_new(f*(k-1)+1:f*k,:) = temp/Net(k);
     end
end

%disp ("******* MAXIMIZATION DONE! *******")




end