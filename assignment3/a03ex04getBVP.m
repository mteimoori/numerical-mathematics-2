% Assignment 03, Exercise 04a/b) by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543)
function [ xh, Lh, fh ] = a03ex04getBVP(p)

% part a): analytical solution for f(x)
syms x;
u(x) = 1 + 4 * x^2 - 3 * x^3;
du = diff(u,x);
ddu = diff(du,x);
f_analytical = -ddu - 4 * du + u;

% part b):
N = 2.^p - 1;
h = 1 / (N + 1);     % mesh size

a = 1;
b = 4;
c = 1;


Sparse_a = sparse(gallery('tridiag',N,1,-2,1));
Sparse_b = sparse(gallery('tridiag',N,-1,0,1));
Sparse_c = speye(N,N);

Lh = -a / h.^2 * Sparse_a - b / (2 * h) * Sparse_b + c * Sparse_c;

xh = zeros(N,1); %creates zero vector for the values of xh
f_analytical = zeros(N,1); %creates zero vector for the analytical solution of f

    for k = 1:N %this for loop computes the values for xh and the analytical solution of f for h1,h2,...,hN
        f_analytical(k) = -7-14*(k * h)+40*(k * h).^2-3*(k * h).^3;
        xh(k) = 1 + 4 * (k * h).^2 - 3 * (k * h).^3;
    end

fh = Lh * xh; %numerical solution of the BVP

figure(p)
x = 2:N-1;
plot(x,fh(2:N-1),x,f_analytical(2:N-1)); %this plot compares the values between the analytical and the numerical solution
end

