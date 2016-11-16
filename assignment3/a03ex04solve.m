% Assignment 03, Exercise 04c) by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543)
function error = a03ex04solve()

a = 1; %these parameters are independant of p
b = 4;
c = 1;

for p = 1:15
    
    N = 2.^p - 1;
    h = 1 / (N + 1);
    
    Sparse_a = sparse(gallery('tridiag',N,1,-2,1));
    Sparse_b = sparse(gallery('tridiag',N,-1,0,1));
    Sparse_c = speye(N,N);

    Lh = -a / h.^2 * Sparse_a - b / (2 * h) * Sparse_b + c * Sparse_c;

    xh = zeros(N,1);
    f_analytical = zeros(N,1);
    
    if p == 1 %this if condition makes sure that the arrays don't get overwritten on every increase of p
    error_abs = zeros(N,1); 
    error = zeros(15,1);
    h_array = zeros(15,1);
    end
    

        for k = 1:N
            f_analytical(k) = -7-14*(k * h)+40*(k * h).^2-3*(k * h).^3;
            xh(k) = 1 + 4 * (k * h).^2 - 3 * (k * h).^3;
        end

    fh = Lh * xh;
    
    
    for l=2:N-1 %this excludes the values for the first and last entry of f
    error_abs(l) = fh(l) - f_analytical(l);   
    end
    error(p)=norm(error_abs,inf);
    h_array(p)=h; %stores the values of h in an array which can be plotted



end

loglog(h_array,error); 
end


