% Assignment 08, Exercise 02, by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543) 
p = 10;
h = 0.02;
k = 0.05;
lam = k/(h^2);
T = 2;
v = @(x) (cos(2.*x)+x);
a08e02solver(p,lam,T,v,h,k)