function u = a02e04WaveEq( t,x,a,b )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

prompt = 'What is the value of t?';
t = input(prompt);
prompt = 'What is the value of x?';
x = input(prompt);
prompt = 'What is the value of a?';
a = input(prompt);
prompt = 'What is the value of b?';
b = input(prompt);


phi_xplusct = a/(((x+t)-3)^2+1) + a/(((x+t)+3)^2+1);

phi_xminusct = a/(((x-t)-3)^2+1) + a/(((x-t)+3)^2+1);

psi = @(x_int) b.*x_int.*exp(-x_int.^2/2);

u = 1/2 * (phi_xplusct + phi_xminusct) + 1/2 * integral(psi,x-t,x+t);

end

