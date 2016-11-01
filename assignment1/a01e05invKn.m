function [x] = a01e05invKn(b)
N = length(b);
Kn = a01e04sparse(N);
[L,U]=lu(Kn);
% Perform the forward substitution
y(1) = b(1);
for k = 2:N
y(k) = b(k) - L(k,k-1)*y(k-1) ;
end
% Perform the backward substitution
x(N) = y(N)/U(N,N);
for k=(N-1):-1:1
x(k) =(y(k)- U(k,k+1)*x(k+1))/U(k,k) ;
end
x = x';
end