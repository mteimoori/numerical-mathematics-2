function x = a01e05invKn( b )
    % Assignment 01, Exercise 05
    % solves K_n*x = b with the LU method
    n = length(b);
    % Initialize vectors x and y.
    y = zeros(n,1);
    x = y;

    % The general LU decomposition for Kn matrix
    % L_{i,i}   = 1 
    % L_{i,i-1} = -(i-1)/i
    % L_{i,j}   = 0 if i<j and |i-j|>1
    % 
    % U_{i,i}   = (i+1)/i 
    % U_{i-1,i} = -1
    % U_{i,j}   =0 if i>j and |i-j|>1
    %       

    % Forward substition with elements of the lower triangular matrix L
    y(1) = b(1);
    for i = 2:n
        y(i) = b(i) + (i-1)*y(i-1)/i;
    end 

    % Backward substitution with the upper triangular matrix U
    x(n) = y(n)/((n+1)/n);
    for i = n-1:-1:1
        x(i) = (y(i)+x(i+1))/(2-(i-1)/i);
    end
end

