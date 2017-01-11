% Assignment 6, Programming exercise 3b, 
% Returns the error err between yh and the restricted exact solution yex


% exact solution
y = @(x) x.^4 - 2*x.^2 + 1/2; 

maxp = 15;
for p = 2:maxp
    N = 2^p+1; 
    h = 1. / (2.^p);
    xh = linspace(-1,1,N)';
    steps(p-1) = h;
    
    yex = y(xh);

    beta = 0;
    [Lh,fh] = a06e03getPDE(p,beta);
    yh = Lh\fh;
    yh = [-0.5; yh];

    errors(p-1) = norm(yex - yh, Inf);
end

figure(1)
plot(xh, yh, 'b--')
hold on
plot(xh, yex, 'r')
hold off
xlabel('x')
ylabel('y(x)')
title('Solution of the BVP')

% behavior at infinity
S = diff(log(errors))./diff(log(steps));
s = mean(S(2:end));

figure(2)
loglog(steps,errors)
xlabel('h')
ylabel('error norm')
title('Convergance plot of the numerical solution')
text(steps(6),errors(7),strcat('O(N^',num2str(s,0),')'))

% test different beta values
Betas = [1,-1,5,-5];
figure(3)
plot(xh, yh, '--')
hold on

for m = 1:length(Betas)
    [Lh,fh] = a06e03getPDE(p,Betas(m));
    yh = Lh\fh;
    yh = [-0.5; yh];
    plot(xh, yh)
end
hold off
xlabel('x')
ylabel('y(x)')
legend('\beta =  0', '\beta =  1', '\beta = -1', '\beta =  5', '\beta = -5')
title('Solution of the BVP for different values of \beta')
