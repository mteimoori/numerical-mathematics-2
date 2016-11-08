% Assignment 01, Exercise 04,05, by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543) 
function [ s ] = a01e04sparse( n )
%returns sparse matrix of Kn
s = sparse(gallery('tridiag', n, -1, 2, -1));
end

