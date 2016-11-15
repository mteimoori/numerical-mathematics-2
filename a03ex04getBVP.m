function [ xh, Lh, fh ] = a03ex04getBVP(p)

N = 2.^p - 1;
h = 1 / (N + 1);     % mesh size

a = 1;
b = 4;
c = 1;

Sparse_a = zeros(N,N); % It ist better to later rewrite values in an existing, fixed-size matrix than to change the size of the matrix every time 
Sparse_b = zeros(N,N);
Sparse_c = speye(N,N); % creates the (0,1,0) matrix

for i = 1:N % these two for-loops add the correct values to the other two matrices, (1,-2,1) and (-1,0,1)
    for j = 1:N
        if i == j
            Sparse_a(i,j) = -2; 
        end
        if j+1 == i
                Sparse_a(i,j) = 1;
                Sparse_b(i,j) = -1;
        end
        if i+1 == j
                    Sparse_a(i,j) = 1;
                    Sparse_b(i,j) = 1;
        end

    end
end

Sparse_a = sparse(Sparse_a);

Sparse_b = sparse(Sparse_b);

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

