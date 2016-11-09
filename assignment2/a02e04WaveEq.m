% Assignment 02, Exercise 03 by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543) 
function [ u ] = a02e04WaveEq(t, x, a, b)
fun = @(l) b.*l.*exp(-0.5*(l.^2));
Qarray = @(a,b)arrayfun(@(ak,bk)integral(fun,ak,bk),a,b);
 u = 0.5.*((a./((x-t-3).^2+1))+(a./((x-t+3).^2+1))+(a./((x+t-3).^2+1))+(a./((x+t+3).^2+1)))+0.5.*Qarray(x-t,x+t);
 
end

