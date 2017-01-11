function [uh] = a06ex03sylsolver(p,i)
% Assignment 6, Programming exercise 4
% Solves the BVPs from A04Ex03 with the Sylvester equation


N    = 2^p - 1;    % number of interior points in one direction
h    = 1/(N + 1);  % lattice spacing

% pick r.h.s. according to i
if i == 1
    f = @(x,y) - 6*x.^3.*y - 6*x.*y.^3 + 12*x.*y; 
else
    f = @(x,y) 10*pi^2*sin(3*pi*x).*sin(pi*y);
end


[xh,yh]=meshgrid(linspace(h,1-h,N)); % create mesh using MATLABs meshgrid

% evaluate the inhomogenity f
fh = f(xh,yh);
% transform into lexicographically ordered vector
Fh = reshape(fh',N,N);

% assemble Lh using the Kronecker product

I = speye(N,N);
E = sparse(2:N,1:N-1,1,N,N);
S = E + E' - 2*I;
% full(S)

% Uh = sylvester(full(S),full(S),-h^2*Fh);
Uh = lyap((S),h^2*Fh);

uh = reshape(Uh, N^2,1);
end
