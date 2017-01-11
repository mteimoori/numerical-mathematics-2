function u = a02e05HeatEqGen(t,x,a)
% Assignment 02, Exercise 05
%
% computes the exact solution for a given right hand side vector of Fourier
% coefficients

% test run a02e05HeatEqGen(0,pi,[0.5,0,0.5])

% Fourier expansion with a vector of coeffitients a
u_func = @(x,t,n) exp(-(pi*n)^2*t).*cos(n*pi*x); 
u = 0;
for i = 1:length(a);
    u = u + a(i).*u_func(x,t,i-1);
end      
