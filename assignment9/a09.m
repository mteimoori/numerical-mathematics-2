N=5000;
coeff = a09e02getpoly(N);

exact = @(x) exp(x)+(-exp(1)-1)*x+1;
x=0;
h= 1/N; 
approximate=zeros(N-1,1);
exact_solution=zeros(N-1,1);
for i=1:N-1
    s = 0;
    for j=1:N
        s = s + coeff(j,1)*((exp(i*h)-1)^j);
    end
    approximate(i,1) = s+2;
    exact_solution(i,1) = exact(i*h);
end
exact_solution;
approximate = approximate ;
error = exact_solution - approximate

figure(1)

plot ((1:1:N-1)',exact_solution)
hold on
plot ((1:1:N-1)',approximate)