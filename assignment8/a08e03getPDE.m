% Assignment 08, Exercise 03a)/b) by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543)

%% 3a) Determine f such that u(x,y,z) = x * (1 - x) * sin(pi * y) * (1 - cos(2 * pi * z))

% syms x y z;
% u = @(x,y,z) (x * (1 - x) * sin(pi * y) * (1 - cos(2 * pi * z)));

%f = diff(u,2,x) + diff(u,2,y) + diff(u,2,z) = 2*sin(pi*y)*(cos(2*pi*z) - 1) - 4*pi^2*x*cos(2*pi*z)*sin(pi*y)*(x - 1) - pi^2*x*sin(pi*y)*(x - 1)*(cos(2*pi*z) - 1)


%% 3b)
function [ Lh,fh ] = a08e03getPDE( p,f )


N = 2^p-1;
h = 1 / (N + 1);

[xh,yh,zh]=meshgrid(linspace(h,1-h,N));

fh = f(xh,yh,zh);
fh_temp = fh(:,:,1);
for i = 2:N
fh_temp = vertcat(fh_temp,fh(:,:,i));
end
fh = reshape(fh_temp',N^3,1);

I = speye(N,N);
E = sparse(2:N,1:N-1,1,N,N);
D = E+E'-3*I;
A = kron(D,I)+kron(I,D);

Lh = kron(A,I)+kron(I,A);

Lh(Lh==-12) = -6;
Lh(Lh==2) = 1;

Lh = 1 / h^2 * Lh;

end

