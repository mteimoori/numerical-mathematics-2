% The following code is testing the a04ex02getPDE for different
% approximations for u' on the example from Assignment 3 Ex. 4

N = 100;    % set the number steps
h = 1.0./(N-1);
xh = linspace(0,1,N)';
consts = [1, -4, 1, 1, 2];
f = @(x) -3. * x.^3 + 40 * x.^2 -14 * x - 7;

[Lh,fh] = a04ex02getPDE(xh,f,consts,'0');
u = Lh\fh;
u_0 = [consts(1); u; consts(5)];

[Lh,fh] = a04ex02getPDE(xh,f,consts,'+');
u = Lh\fh;
u_p = [consts(1); u; consts(5)];

[Lh,fh] = a04ex02getPDE(xh,f,consts,'-');
u = Lh\fh;
u_m = [consts(1); u; consts(5)];

%%
plot(xh, u_0,'r')
hold on
plot(xh, u_p,'b')
plot(xh, u_m,'g')
hold off
legend('D^0 approximation', 'D^+ approximation', 'D^- approximation','Location','southeast')
title(['Numerical solution for N=', num2str(N),' inner steps'])
