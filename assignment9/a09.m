% Assignment 09, Exercise 02, by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543) 
for NN=10:10:300
N=NN;
coeff = a09e02getpoly(N);
exact = @(x) exp(x)-(exp(1)+1)*x-1;
h= 1/N;
%initial vectors%
approximate=zeros(N,1);
exact_solution=zeros(N,1);
for i=1:N
    s = 0;
    for j=1:N
        s = s + coeff(j,1)*(i*h)^j;
    end
    approximate(i,1) = s;
    exact_solution(i,1) = exact(i*h);
end
steps(N,1) = h;
error(N,1) = norm((exact_solution - approximate),Inf);
end
figure(1)
steps
error
size1 = size(steps)
size2=size(error)
loglog (steps,error)
xlabel('h')
ylabel('error norm')
figure(2)
plot ((h:h:1)',exact_solution)
hold on
plot ((h:h:1)',approximate)

% part f%
%because Galerkin approximation calculates the inverse of a full matrix which takes a lot og memory and computation time%
