function [xh, Lh, fh] = a03ex04getBVP(p)
%% function [xh, Lh, fh] = a03ex04getBVP(p)
%% Return the inner grid discretization xh, the matrix Lh and the right-hand
%% side fh corresponding to the problem -u''-4u'+u = f

%% Grid discretization
a=0; b=1 ; N = 2^p-1 ; h = 1. / (2.^p);
xh = linspace(a+h,b-h,N)';

%% Matrix
% -D^- D^+ uh(ih) = 1/h^2 ( -uh((i-1)h) + 2uh(ih) - uh((i+1)h) )
minusDminusDplush = spdiags(ones(N,1)*[-1,2,-1],[-1,0,1],N,N)/h^2;
% -4D^0 uh(ih) = 4/(2h) ( uh((i-1)h) - uh((i+1)h) ) 
minusfourDzero = spdiags(ones(N,1)*[1,-1],[-1,1],N,N)*2/h;
% Lh = -D^- D^+ uh(ih) - 4D^0 uh(ih) + I uh(ih)
Lh = minusDminusDplush + minusfourDzero + speye(N,N);

%% Right-hand side

% Part coming from the computation of f(ih)
fh = -3. * xh.^3 + 40 * xh.^2 -14 * xh - 7;
uh_a = 1;
uh_b = 2;

% Part coming from the transmission of the boundary conditions
fh(1) = fh(1) + (1/h^2-2/h) * uh_a;
fh(N) = fh(N) + (1/h^2+2/h) * uh_b;