function [ s ] = a01e04sparse( n )
%returns sparse matrix of Kn
s = sparse(gallery('tridiag', n, -1, 2, -1));
end

