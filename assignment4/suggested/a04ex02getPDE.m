function [Lh,fh] = a04ex02getPDE(xh,f,consts,flag)
% Assignment 4, Programming exercise 2, Two point BVP
% Compute the operator Lh according to the flag
% the right-hand side vector fh 
% on the inner grid xh
% corresponding to the problem 
%   -a*u''+b*u'+c*u = f
%   u(o) = alpha,  u(1) = beta
% where the model parameter are contained in
% consts = [a, b, c, alpha, beta]

% Number of internal points 
N = length(xh)-2; 
% mesh steps
h(1:N+1) = xh(2:N+2) - xh(1:N+1); 

% the right hand side from the equation
fh = f(xh(2:N+1));

if flag == '+'  % D+
    Dfirst = spdiags([-1./h(2:N+1)' 1./h(2:N+1)'],0:1,N,N);
    fh(1) = fh(1) + consts(4)*(consts(1)*2/(h(1)*(h(1)+h(2))));
    fh(N) = fh(N) + consts(5)*(consts(1)*2/(h(N)*(h(N)+h(N+1))) -...
            consts(2)*(1/h(N+1)));
elseif flag == '-' % D-
    Dfirst = spdiags([-1./h(1:N)' 1./h(1:N)'],-1:0,N,N);
    fh(1) = fh(1) + consts(4)*(consts(1)*2/(h(1)*(h(1)+h(2))) +...
            consts(2)/h(1));
    fh(N) = fh(N) + consts(5)*(consts(1)*2/(h(N)*(h(N)+h(N+1))));
elseif flag == '0' %D0
    Dfirst = spdiags([-1./(h(1:N)+h(2:N+1))' 1./(h(1:N)+h(2:N+1))'],[-1 1],N,N);
    fh(1) = fh(1) + consts(4)*(consts(1)*2/(h(1)*(h(1)+h(2))) + ...
            consts(2)/(h(1)+h(2)));
    fh(N) = fh(N) + consts(5)*(consts(1)*2/(h(N)*(h(N)+h(N+1))) -...
            consts(2)*(1/(h(N)+h(N+1))));
else error(['wrong flag (options: ''', '+''', ' ''', '-''', ' ''', '0'')'])
    
end

%D+-
DplusDminus = spdiags([2./(h(1:N).*(h(1:N)+h(2:N+1)))' -2./(h(1:N).*h(2:N+1))'...
                       2./(h(2:N+1).*(h(1:N)+h(2:N+1)))'],-1:1,N,N);

% form of Lh follows from given problem
Lh = - consts(1)*DplusDminus + consts(2)*Dfirst + consts(3)*speye(N,N);

% fh = f(xh(2:N+1));

% Part coming from the transmission of the boundary conditions
% fh(1) = fh(1) + (1/h(1)^2-2/h(1)) * consts(4);
% fh(N) = fh(N) + (1/h(N+1)^2+2/h(N+1)) * consts(5);
