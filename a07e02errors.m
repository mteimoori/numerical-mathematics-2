clear;
clc;

v=@(x)(exp(-x.^2));
f=@(x)(-x);
for p = 3:10

xmin=-1;
xmax=1;
t=1.5;
N=2^p;
h=1/N;
k=h/2;
l=t/k;
xh=(xmin-l:h:xmax+l+1)';


U = a07e02LaxFrie( t,xmin,xmax,v,f,h,k );

u = exp(-(xh+t).^2);

error(p) = max(abs(U-u));

end


