function [Lhk,fhk] = a08e02getPDE(h,k,N)
% Assignment 8, Programming exercise 2b
a = 1;
alpha = 1;
beta = 1;
lambda = k/h^2;
Lhk = gallery('tridiag', N, -a*lambda, 1+2*a*lambda, -a*lambda);
%% RHS
fhk = zeros(N,1);
% BC terms
fhk(1) = a*lambda*alpha;
fhk(N) = a*lambda*beta;
end
