function [Lh,fh] = a04ex03getPDE(p,i)
% Assignment 4, Programming exercise 3, Sample solution
% Compute the operator Lh and the indices of boundary 
% points for a Poisson problem on the domain (0,1)^2 with N = 2^p - 1 
% interior points in each direction using a five point stencil.
%
%       [    -1    ]
% 1/h^2 [ -1 +4 -1 ]
%       [    -1    ]
%
% Note:
%   i determines the problem to solve
%
%   Lh is a N^2 x N^2 sparse matrix
%   fh is a N^2 vector representing the inhomogenity at the grid points
%
%
% Example:
%       [Lh,fh] = a04ex03getPDE(p,i);


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
fh = reshape(fh',(N)^2,1);

% assemble Lh using the Kronecker product

I = speye(N,N);
E = sparse(2:N,1:N-1,1,N,N);
T = E + E' - 2*I;
Lh = -1/h^2*(kron(T,I) + kron(I,T));
end
