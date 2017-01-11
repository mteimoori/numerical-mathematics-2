function [err, uex] = a05ex04error(eps,xh,uh)
% Assignment 5, Programming exercise 4b, 
% Returns the error err between uh and the restricted exact solution uex


% exact solution
u = @(x) x - (exp(-(1-x)/eps) - exp(-1/eps))./(1-exp(-1/eps));

uex = u(xh);

err = norm(uex - uh, Inf);
