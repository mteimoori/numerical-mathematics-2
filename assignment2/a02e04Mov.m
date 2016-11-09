% Assignment 02, Exercise 03 by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543)
function u = a02e04Mov(a,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
prompt = 'What is the value of a?';
a = input(prompt);
if ~isa(a, 'double')
   error('Bad type!');
end

prompt = 'What is the value of b?';
b = input(prompt);
if ~isa(b, 'double')
   error('Bad type!');
end

t_range = 0:0.1:10;
x_range = -15:0.5:15;
t_range_size = size(t_range);
x_range_size = size(x_range);
z = t_range_size(2) * x_range_size(2);

u = zeros(t_range_size(2),x_range_size(2));

set_y_axis = max([abs(a) abs(b)]) + (max([abs(a) abs(b)])) / 5;

curve = animatedline;
set(gca,'XLim',[-15 15],'YLim',[-set_y_axis set_y_axis]);

psi = @(x_int) b.*x_int.*exp(-x_int.^2/2);

for i = 1:t_range_size(2)
    delete(curve);
    curve = animatedline;
    for j = 1:x_range_size(2)
        
        phi_xplusct = a/(((x_range(j)+t_range(i))-3)^2+1) + a/(((x_range(j)+t_range(i))+3)^2+1);

        phi_xminusct = a/(((x_range(j)-t_range(i))-3)^2+1) + a/(((x_range(j)-t_range(i))+3)^2+1);

        
    
        u(i,j) = 1/2 * (phi_xplusct + phi_xminusct) + 1/2 * integral(psi,x_range(j)-t_range(i),x_range(j)+t_range(i));
        
       
    end
    
    addpoints(curve,x_range,u(i,:));
       drawnow
end

end
