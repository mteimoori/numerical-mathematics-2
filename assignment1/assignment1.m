% Assignment 01, Exercise 04,05, by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543) 
% exercise 4, part a
Kn = a01e04sparse(1000);
% exercise 4, part b
    %{it stores only the nonzero elements of the matrix with their indices and reduces computation time by eliminating operations on zero elements.this leads to saving computaional time and memory usage.%}
% exercise 5, part a
[L,U]=lu(Kn);
b = rand(1000,1);
x = a01e05invKn(b);
x_check = Kn\b;
isEqual = isequal(x,x_check); %check if our function works correctly
isEqual 
% exercise 5, part b
for dim = [10, 100, 1000]
    Kn = a01e04sparse(dim);
    b = rand(dim,1);
    tic
        Knb = Kn\b;
   tictoced = toc;
        fprintf('timing for Kn\\b when size = %d is %f \n', dim,tictoced);
    tic    
        fullKnb = full(Kn)\b;
    tictoced = toc;
        fprintf('timing for full(Kn)\\b when size = %d is %f \n', dim,tictoced);
   tic
        ours = a01e05invKn(b);
   tictoced = toc;
        fprintf('timing for a01e05invKn when size = %d is %f \n', dim,tictoced);
    tic
        invFullKnb = inv(full(Kn))*b;
    tictoced = toc;
        fprintf('timing for inv(full(Kn))*b when size = %d is %f \n', dim,tictoced);
    
end
%time increases respectively 
% exercise 5, part c
points = 10;
t = zeros(points);
dim = 10;
for n = 1:points
    dim = dim * 2;
    Kn = a01e04sparse(dim);
    b = rand(dim,1);
    tic;
    x = a01e05invKn(b);
    t(n) = toc;
end
plot(t)
% as you see it grows exponentialy when n goes toward infinity
