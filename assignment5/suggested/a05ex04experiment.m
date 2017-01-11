% Assignment 05, Programming exercise 4d
% plots the exact and approximated solution
close all
% set the values 
eps = 0.001;
N   = [5,50,500,5000];
flag={'+','0','-','0'};
k   = length(N);

% compute the error and present the results for the four cases

for m=1:4
    figure(m)
    hold on
    for j=1:k
        if m==4
            sigma=4*eps*log(2*N(j));
            xh = a03ex03shishkin(N(j),sigma);
        else
            xh = linspace(0,1,2*N(j)+1)';
        end
        
        uh = a05ex04solve(eps,xh,flag{m});
        [err(j),uex] = a05ex04error(eps,xh,uh);
        plot(xh,uh)
        drawnow
    end

    plot(xh,uex,'r-')
    if m~=4
        string=sprintf('Operator D%s, uniform grid\nError Inf-norm: %.3e, %.3e, %.3e, %.3e',flag{m},err);
    else
        string=sprintf('Operator D%s, shishkin grid\nError Inf-norm: %.3e, %.3e, %.3e, %.3e',flag{m},err);
    end
    title(string)
    hold off
    xlabel('x')
    ylabel('Solution u')
    legend('N=5','N=50','N=500','N=5000', 'exact','Location','Northwest')
    drawnow
end

