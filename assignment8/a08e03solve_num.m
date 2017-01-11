% Assignment 08, Exercise 03d) by Georgia Markouleki(387232), Manuel Widdel(379704), Mohammad Teimoori(370543)
function errors_num = a08e03solve_num

f = @(x,y,z) (exp( -10 .* ((x - 0.5).^2 + (y - 0.5).^2 + (z - 0.5).^2)));

Nsteps = 4;
% steps = zeros(Nsteps);
% errors  = zeros(Nsteps);

for p = 1:Nsteps
    

%% set up coarse grid
N_coarse = 2^p - 1;
h_coarse = 1 / (N_coarse + 1);
steps_coarse(p) = h_coarse;


[xh,yh,zh]=meshgrid(linspace(h_coarse,1-h_coarse,N_coarse));

fh_coarse = f(xh,yh,zh);
fh_coarse_temp = fh_coarse(:,:,1);
for i = 2:N_coarse
fh_coarse_temp = vertcat(fh_coarse_temp,fh_coarse(:,:,i));
end
fh_coarse = reshape(fh_coarse_temp',N_coarse^3,1);

I_coarse = speye(N_coarse,N_coarse);
E_coarse = sparse(2:N_coarse,1:N_coarse-1,1,N_coarse,N_coarse);
D_coarse = E_coarse+E_coarse'-3*I_coarse;
A_coarse = kron(D_coarse,I_coarse)+kron(I_coarse,D_coarse);

Lh_coarse = kron(A_coarse,I_coarse)+kron(I_coarse,A_coarse);

Lh_coarse(Lh_coarse==-12) = -6;
Lh_coarse(Lh_coarse==2) = 1;

Lh_coarse = 1 / h_coarse^2 * Lh_coarse;

%% set up fine grid

N_fine = 2^(p + 1) - 1;
h_fine = 1 / (N_fine + 1);
steps_fine(p) = h_fine;

[xh,yh,zh]=meshgrid(linspace(h_fine,1-h_fine,N_fine));

fh_fine = f(xh,yh,zh);
fh_fine_temp = fh_fine(:,:,1);
for i = 2:N_fine
fh_fine_temp = vertcat(fh_fine_temp,fh_fine(:,:,i));
end
fh_fine = reshape(fh_fine_temp',N_fine^3,1);

I_fine = speye(N_fine,N_fine);
E_fine = sparse(2:N_fine,1:N_fine-1,1,N_fine,N_fine);
D_fine = E_fine+E_fine'-3*I_fine;
A_fine = kron(D_fine,I_fine)+kron(I_fine,D_fine);

Lh_fine = kron(A_fine,I_fine)+kron(I_fine,A_fine);

Lh_fine(Lh_fine==-12) = -6;
Lh_fine(Lh_fine==2) = 1;

Lh_fine = 1 / h_fine^2 * Lh_fine;
        
%% compute uh_coarse and uh_fine and the error

uh_coarse = Lh_coarse\fh_coarse;
uh_fine = Lh_fine\fh_fine;


%% cut off the last value of uh_fine until the ratio of the sizes is an integer
    size_coarse = size(uh_coarse);
    size_fine = size(uh_fine);
    
    inter = size_fine(:,1) / size_coarse(:,1);
    reminder = rem(inter,1);
    while reminder ~= 0
    
    uh_fine = uh_fine(1:end-1);
    size_fine = size(uh_fine);
    inter = size_fine(:,1) / size_coarse(:,1);
    reminder = rem(inter,1);
    end  
    
    uh_fine = uh_fine(1:inter:end);

    errors(p) = norm(uh_coarse - uh_fine, Inf);
    
    EOCp = diff(log(errors))./diff(log(steps_coarse));
    


    
    
    p_array(p)= p; 
end
 
refinement = (2:1:p)'
steps = (steps_coarse(2:end))'
errors = (errors(2:end))'
EOCp = EOCp'
Table = table(refinement,steps,errors,EOCp)



end

