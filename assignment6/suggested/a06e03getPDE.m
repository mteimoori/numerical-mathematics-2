function [Lh,fh] = a06e03getPDE(p,beta)
% Assignment 6, Programming exercise 3a
% Solves the problem the BVP
% (2-x^2) y'' - xy' + 16y = 0, x\in (-1,1)
% y(-1) = -1/2, y'(1) = \beta
% Return the matrix Lh and the right-hand
% side fh corresponding to the problem
% 

%% Grid discretization
N = 2^p;
h = 2. / (2.^p);
xh = linspace(-1+h,1,N)';

% coefficiant (2-x^2)
ax = 2-xh.^2;
% coefficiant (-x)
bx = -xh;
% coefficiant 162
c = 16;

% D+-
DplusDminus = gallery('tridiag', ax(2:N)./h^2, -2*ax(1:N)./h^2, ax(1:N-1)./h^2);
DplusDminus(N,N-1) = 2/h^2;
% full(DplusDminus*h^2)

% D0
Dfirst = gallery('tridiag',-bx(2:N)./(2*h), zeros(N,1), bx(1:N-1)./(2*h));
Dfirst(N,N-1) = 0.0;
% full(Dfirst)


% form of Lh follows from given problem
Lh = DplusDminus + Dfirst + c*speye(N,N);

%% Right-hand side
% Part coming from the computation of f(ih)
fh = zeros(N,1);

% BC terms
fh(1) = 0.5*(1/h^2*ax(1) - 1/(2*h)*bx(1));
fh(N) = beta - 2*beta/h;



