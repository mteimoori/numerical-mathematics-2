function u = a02e05HeatEq(t,x)
% Assignment 02, Exercise 05
%
% computes the exact solution of the heat equation for given t \in (0,\infty), x \in R.
% Solution by separation of variables

u = 0.5+0.5*exp(-4*t).*cos(2*x); 