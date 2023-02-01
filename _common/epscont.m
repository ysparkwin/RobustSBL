function x = epscont(n,m,sigma,epsilon,mu,lambda)
% Generate an array of n x m random deviates from an epsilon-contaminated 
% complex Gaussian distribution:
%     f(x)= (1-epsilon)*CN(mu,sigma^2) + epsilon*CN(mu,(lambda sigma)^2)
% The variance is varx= (1-epsilon)*sigma^2 + epsilon * (lambda*sigma)^2
%
% INPUT:
% n   = row length (positive integer) of the array
% m   = column length (positive integer) of the array
% mu  = complex number ( mean of random variable x)
% lambda = multiplier for epsilon contaminated variance (positive real) 
% sigma = variance of true normal  (positive real) 
%
% USAGE: 
% n = 500; epsilon = 0.2; mu = 1+1i; lambda=2; sigma=1;
% x = epscont(n,epsilon,mu,lambda,sigma);
%   =_d (1-eps)*n_1 + eps*n_2 
% var(x)
% varx = (1-epsilon)*sigma^2 + epsilon * (lambda*sigma)^2

if nargin < 4
    mu = 0;
    lambda = 1;
    epsilon = 0;
end
    
if epsilon~=0 
    b = logical(binornd(1,epsilon,n,m));
    cnt = sum(b(:));
    n2 = (lambda*sigma/sqrt(2))*(randn(cnt,1)+ 1i*randn(cnt,1))+mu;
    n1 = (sigma/sqrt(2))*(randn(n*m-cnt,1)+ 1i*randn(n*m-cnt,1)) + mu;
    x = zeros(n,m); 
    x(b) = n2;
    x(~b) = n1;
else
    x= sigma * complex(randn(n,m),randn(n,m))/sqrt(2);
end