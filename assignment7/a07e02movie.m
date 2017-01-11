clear;
clc;

v=@(x)(exp(-x.^2));
f=@(x)(x.^2);
p=4;

for t = 0:0.05:1.5



xmin=-3;
xmax=4;
N=2^p;
h=1/N;
k=h/2;
l=t/k;
xh=(xmin-l:h:xmax+l+1)';


U = a07e02LaxFrie( t,xmin,xmax,v,f,h,k );


figure(1)
plot(xh,U);
xlim([-3 4]);
ylim([0 1]);

pause(1/25);



end



