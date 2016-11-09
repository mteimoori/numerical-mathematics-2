function a02e04Surf()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

t_range = 0:0.2:10;
x_range = -15:0.2:15;
t_range_size = size(t_range);
x_range_size = size(x_range);
z = t_range_size(2) * x_range_size(2);

U = zeros(t_range_size(2),x_range_size(2));


a=0;
b=1;

psi = @(x_int) b.*x_int.*exp(-x_int.^2/2);

for i = 1:t_range_size(2)
    
    for j = 1:x_range_size(2)
        
        phi_xplusct = a/(((x_range(j)+t_range(i))-3)^2+1) + a/(((x_range(j)+t_range(i))+3)^2+1);

        phi_xminusct = a/(((x_range(j)-t_range(i))-3)^2+1) + a/(((x_range(j)-t_range(i))+3)^2+1);

        
        U(i,j) = 1/2 * (phi_xplusct + phi_xminusct) + 1/2 * integral(psi,x_range(j)-t_range(i),x_range(j)+t_range(i));

    end
    
end
name_figure1 = 'a=0, b=1';
figure('Name',name_figure1,'NumberTitle','off');
surf(x_range,t_range,U)

a=1;
b=0;

psi = @(x_int) b.*x_int.*exp(-x_int.^2/2);

for i = 1:t_range_size(2)
    
    for j = 1:x_range_size(2)
        
        phi_xplusct = a/(((x_range(j)+t_range(i))-3)^2+1) + a/(((x_range(j)+t_range(i))+3)^2+1);

        phi_xminusct = a/(((x_range(j)-t_range(i))-3)^2+1) + a/(((x_range(j)-t_range(i))+3)^2+1);

        
        U(i,j) = 1/2 * (phi_xplusct + phi_xminusct) + 1/2 * integral(psi,x_range(j)-t_range(i),x_range(j)+t_range(i));

    end
    
end
name_figure2 = 'a=1, b=0';
figure('Name',name_figure2,'NumberTitle','off');
surf(x_range,t_range,U)


a=1;
b=1;

psi = @(x_int) b.*x_int.*exp(-x_int.^2/2);

for i = 1:t_range_size(2)
    
    for j = 1:x_range_size(2)
        
        phi_xplusct = a/(((x_range(j)+t_range(i))-3)^2+1) + a/(((x_range(j)+t_range(i))+3)^2+1);

        phi_xminusct = a/(((x_range(j)-t_range(i))-3)^2+1) + a/(((x_range(j)-t_range(i))+3)^2+1);

       
        U(i,j) = 1/2 * (phi_xplusct + phi_xminusct) + 1/2 * integral(psi,x_range(j)-t_range(i),x_range(j)+t_range(i));

    end
    
end
name_figure3 = 'a=1, b=1';
figure('Name',name_figure3,'NumberTitle','off');
surf(x_range,t_range,U)
end