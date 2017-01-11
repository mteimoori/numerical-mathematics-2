clear

% Maximum mesh size (logscale)
pmax=15;
% Mesh size
h=1. ./ (2.^(1:pmax));
% Allocation of the error vector
error = zeros(size(h));

for p=1:pmax;
  % Get back the discrete problem datas
  [xh, Lh, fh] = a03ex04getBVP(p);
  % Solve the problem
  uh = Lh \ fh;
  % Compute the exact solution
  uexh = 1+4*xh.^2-3*xh.^3;
  % Compute the error
  error(p) = max(abs(uh-uexh));
end
%%
S = diff(log(error))./diff(log(h));
s = mean(S);

loglog(h,error);
% plot(h,error);
title('Numerical error vs. mesh size')
ylabel('error')
xlabel('h')
text(h(8),error(6),'O(h^2)')