function [Uhk] = a08e02solver(p,lam,T,v, h,k)
% Assignment 8, Programming exercise 2c
N = 2^p;
m = T/k;
[Lhk,fhk] = a08e02getPDE(h,k,N);
% Setup mesh in x-direction:
x1 = 0;
x2 = pi;
x = linspace(x1+h,x2, N)';
% Setup vector U^0 at time l=0: U(x,0) = v(x)
Ul = v(x);
sizeu = size(Ul);
sizeu = size(fhk);
% plot inital condition at tk =0 %
plot (x,Ul)
title ('Initial condition at tk=0')
xlabel ('x')
ylabel ('U')
%initialization
t=0;
U=zeros(N,m);
time_vector=zeros(m,1);
size(Ul)
U(:,1) = Ul;
%iterate on time
for tstep=1:m
    time_vector(tstep)=t;
    t = t+k;
    Ul_1 = Lhk\(Ul+fhk); 
    U(:,tstep) = Ul_1;
    Ul = Ul_1;
end
figure(2)
mesh (time_vector,x,U)
title ('Backward Euler method')
xlabel ('t')
ylabel ('x')
zlabel ('U')