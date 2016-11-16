function error = a03ex04solve()

a = 1; %these parameters are independant of p
b = 4;
c = 1;

for p = 1:13 %for p > 13, error message: "Out of memory"
    
    N = 2.^p - 1;
    h = 1 / (N + 1);
    
    Sparse_a = zeros(N,N);
    Sparse_b = zeros(N,N);
    Sparse_c = speye(N,N);

    for i = 1:N
        for j = 1:N
            if i == j
                Sparse_a(i,j) = -2;
            end
            if j+1 == i
                    Sparse_a(i,j) = 1;
                    Sparse_b(i,j) = -1;
            end
            if i+1 == j
                        Sparse_a(i,j) = 1;
                        Sparse_b(i,j) = 1;
            end

        end
    end

    Sparse_a = sparse(Sparse_a);

    Sparse_b = sparse(Sparse_b);

    Lh = -a / h.^2 * Sparse_a - b / (2 * h) * Sparse_b + c * Sparse_c;

    xh = zeros(N,1);
    f_analytical = zeros(N,1);
    
    if p == 1 %this if condition makes sure that the arrays don't get overwritten on every increase of p
    error_abs = zeros(N,1); 
    error = zeros(13,1);
    error_log = zeros(13,1);
    h_array = zeros(13,1);
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


