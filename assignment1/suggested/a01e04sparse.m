function Kn = a01e04sparse( n )
    % Assignment 01, Exercise 04
    % create n-dim sparse matrix with structure from exercise 1
    Kn = spdiags(ones(n,1)*[-1 2 -1],[-1 0 1],n,n); 
end

