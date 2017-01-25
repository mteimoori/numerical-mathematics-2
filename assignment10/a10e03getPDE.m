function [x,U] = a10e03getPDE(a,b,alpha,beta,f,N)
syms x;
phi=(beta-alpha)/(b-a)*(x-a)+alpha;
phider1=diff(phi);
phider2=diff(diff(phi));
for i=1:N, 
    psi(1,i)=x^(i-1); 
end
psi=(x-a)*(x-b)*psi;
psider1=diff(psi);
psider2=diff(diff(psi));
psitrans=transpose(psi);
E1=-psider2+psi;
E2=-phider2+phi-f;
c=-(int(psitrans*E1,x,a,b))\int(psitrans*E2,x,a,b);
e=E1*c+E2;

approx=phi+psi*c;
x=(a-10):.2:(b+10);
y_approx=subs(approx);
approx=phi+psi*c;
x=(a-10):.2:(b+10);
y_approx=subs(approx);

clear x;
syms x;
exac=dsolve('a2*D2y+a1*Dy+a0*y=f','y(a)=phi1','y(b)=phi2','x');
exac=subs(exac);
x=(a-10):.2:(b+10);
y_exac=subs(exac);

e_norm2=sqrt(int(e^2,a,b));
f_norm2=sqrt(int(approx^2,a,b));
err=abs(vpa(e_norm2/f_norm2*100,2));

plot(x,y_approx,'-o','Color','r'),grid
axis([a-pi b+pi -100 100])
hold on
title('\fontsize{14} Solutions of the Ordinary Differential Equation through Galerkin Method')
xlabel('x')
ylabel('\Phi(x)')
plot(x,y_exac),grid
legend('Approximative solution','Exact solution')
hold off
grid on
end