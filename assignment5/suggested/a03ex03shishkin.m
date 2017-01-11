function xh = a03ex03shishkin(N, sigma)
% Assignment 5, Programming exercise 4c
% generates a column vector of size 2N+1 describing a "Shishkin" grid xh

% set the maximum value for sigma to use Shishkin grid
tol = 0.5;

% two different step sizes
H = (1-sigma)/N;
h = sigma/N;

xh = (0:1:(2*N))';

if sigma>tol
    xh = linspace(0,1,2*N)';
else 
    xh(1:N)=xh(1:N,:)*H;
    xh(N+1:end)=(1-sigma)+(xh(N+1:end)-N)*h;
end


