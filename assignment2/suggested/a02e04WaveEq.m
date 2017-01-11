function u = a02e04WaveEq(t,x,a,b)
% Assignment 02, Exercise 04
%
% computes the exact solution of the wave equation 
% for given t \in (0,\infty), x,a,b \in R
% using D'Alemberts formular:
%   u(x,t) = 1/2*(\phi(x-c*t)+g(x+c*t)) 
%          + 1/2c * int_{x-c*t}^{x+c*t}(\psi(y) d{y}
% where \phi and \psi are the initial data. Here, c=1.

% \phi is given by
phi     = @(x) a./((x-3).^2+1)+a./((x+3).^2+1); 
% the integral of \psi is
psi_int = @(x) -b*exp(-0.5*x.^2);

u = 0.5*(phi(x-t)+phi(x+t)) + 0.5 * (psi_int(x+t)-psi_int(x-t));