% Assignment 10, by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543) 
function [x,U] = a10e03getPDE(a,b,alpha,beta,f,N)
syms x;
n = N+2; % number of points
nel = n-1; % number of elements
h = (b-a)/(nel); % Element length
x_vec = a:h:b;

% inital shape functions
for i=1:n-1
    shape_funcs(i,1) = hat(x_vec(i),i,h);
end
%construct LHS
A = zeros(n,n);
for i=1:2:n-1
    A(i,i) = int(diff(shape_funcs(1), x)*diff(shape_funcs(1), x) - shape_funcs(1)*shape_funcs(1), x_vec(i),x_vec(i+1));
    A(i,i+1) = int(diff(shape_funcs(1), x)*diff(shape_funcs(2), x) - shape_funcs(1)*shape_funcs(2), x_vec(i),x_vec(i+1));
    A(i+1,i) = int(diff(shape_funcs(2), x)*diff(shape_funcs(1), x) - shape_funcs(2)*shape_funcs(1), x_vec(i),x_vec(i+1));
    A(i+1,i+1) = int(diff(shape_funcs(2), x)*diff(shape_funcs(2), x) - shape_funcs(2)*shape_funcs(2), x_vec(i),x_vec(i+1));
end
%construct RHS: missing
%soling system: missing
end
% shape function on elements
function [phi] = hat(xi,i,h)
    syms x    
    if(mod(i,2)==1)
        phi = (x - xi)/h;
    else
        phi = (xi - x)/h;
    end
end