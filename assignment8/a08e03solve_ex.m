% Assignment 08, Exercise 03c) by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543)
function errors = a08e03solve_ex

u = @(x,y,z) (x .* (1 - x) .* sin(pi .* y) .* (1 - cos(2 .* pi .* z)));

f = @(x,y,z) (2.*sin(pi.*y).*(cos(2.*pi.*z) - 1) - 4.*pi^2.*x.*cos(2.*pi.*z).*sin(pi.*y).*(x - 1) - pi^2.*x.*sin(pi.*y).*(x - 1).*(cos(2.*pi.*z) - 1));

Nsteps = 6;
steps = zeros(Nsteps);
errors  = zeros(Nsteps);

for p = 1:Nsteps
    
N = 2^p-1;
h = 1 / (N + 1);
steps(p) = h;
p_array(p)= p;

[xh,yh,zh]=meshgrid(linspace(h,1-h,N));
u_eval = u(xh,yh,zh);
u_eval_temp = u_eval(:,:,1);
for i = 2:N
u_eval_temp = vertcat(u_eval_temp,u_eval(:,:,i));
end
u_eval = reshape(u_eval_temp',N^3,1);



    [Lh, fh] = a08e03getPDE(p,f); 
    
    uh = Lh\fh;
        
    errors(p) = norm(u_eval - uh, Inf);
    

end
figure(1)
loglog(steps,errors) 
xlabel('h')
ylabel('error norm')
legend('errors(p)','Location','northeast')
title('Error between the computed approximation and the restricted exact solution in the maximum norm for p(1...7)')


EOCp = diff(log(errors))./diff(log(steps));

refinement = p_array(2:end)';
steps = steps(2:end,1);
errors = errors(2:end,1);
EOCp = EOCp(:,1);
Table = table(refinement,steps,errors,EOCp)
end


