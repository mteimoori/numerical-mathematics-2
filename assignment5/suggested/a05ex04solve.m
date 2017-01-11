function uh = a05ex04solve(eps,xh,flag)
% Assignment 5, Programming exercise 4a, Sample solution
% Solves the problem -eps*u''+u' = 1
% for for given eps, grid xh, and approximation
% for the first derivative selected by flag i
% xh - is the full grid (with boundary elements)
% uh - is the solution on the full grid

% Number of internal points 
N = length(xh)-2; 
% mesh steps
h(1:N+1) = xh(2:N+2) - xh(1:N+1); 
fh=ones(N,1);
% flag
% if you use D+/- instead of D0 for the approximation of the first order
% derivative, the elements of the sub- and superdiagonal will change.
% The term under the fraction bar will either be h(i) (in case of D-) or h(i+1) 
% (in case of D+) instead of the sum (h(i)+h(i+1)) - see definition of D+,D-
if flag == '-'  % D-
    Dfirst = gallery('tridiag', -1./h(2:N), 1./h(1:N), zeros(N-1,1));
elseif flag == '+' % D+
    Dfirst = gallery('tridiag',zeros(N-1,1), -1./h(1:N), 1./h(1:N-1));
elseif flag == '0' %D0
    Dfirst = gallery('tridiag',-1./(h(2:N)+h(3:N+1)), zeros(N,1), 1./(h(1:N-1)+h(2:N)));
%     full(Dfirst)
else error(['wrong flag (options: ''', '+''', ' ''', '-''', ' ''', '0'')'])
    
end

%D+-
DplusDminus = gallery('tridiag',...
                      2./(h(2:N).*(h(2:N)+h(3:N+1))),...
                      -2./(h(1:N).*h(2:N+1)),...
                      2./(h(2:N).*(h(1:N-1)+h(2:N))));

% form of Lh follows from given problem
Lh = - eps*DplusDminus + Dfirst;

uh_in = Lh\fh;

uh = [0; uh_in; 0];
 



