% Assignment 02, Exercise 03 by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543) 
function [ ] = a02e04Surf()
x = -15:.2:15;
t = (0:.2:10);
[X,T] = meshgrid(x,t);
count=0;
for dim = [0 1;1 0;1 1]'
    count=count+1;
    U = a02e04WaveEq(T, X, dim(1),dim(2));
    subplot(3,1,count)
    surf(X, T, U);
end

end

