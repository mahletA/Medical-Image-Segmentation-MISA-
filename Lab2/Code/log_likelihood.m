function [l] = log_likelihood(X,alpha,Mu,Sigma)
N = size(X,1);
f = size(X,2);
K = size(alpha,2);

temp = zeros(N,1);
for k =1:K
    if k == 1
        temp =  temp +mvnpdf(X ,Mu(k,:), Sigma(1:f,:))*alpha(k);
    else
        temp = temp + mvnpdf(X, Mu(k,:), Sigma((f*(k-1))+1: f*k,:))*alpha(k);
    end
end

l = sum(log(temp));

end