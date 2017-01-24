function [coeff] = a09e02getpoly(N)
A=zeros(N,N);
b=zeros(N,1);
    for k=1:N
        for p=1:N
            A(k,p)= (k*p)/(k+p-1);
        end
    end
    for e=1:N
        fun = @(x) -exp(x).*(x.^e);
        b(e,1) =  integral(fun,0,1)-1;
    end
    coeff = A\b;
end