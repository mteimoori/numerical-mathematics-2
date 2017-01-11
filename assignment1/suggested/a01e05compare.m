%% Assignment 01, Exercise 05 b-c
% script that compares the runtime and the residual of different methods

%% compare different methods runtime and the residual
Ns = [10, 100, 1000];
runtime = zeros(3,4);

residuum_LU = zeros(1,3);
residuum_backslash = zeros(1,3);
residuum_backslash_full = zeros(1,3);
residuum_inverse_full = zeros(1,3);

% set the value how may time to run computation to compute the runtime 
T = 100; 

for j=1:3
    n = Ns(j);
    b = rand(n,1);
    Kn = a01e04sparse(n);

    % how many time to run 
    % t_LU = zeros(T);
    % t_backslash = zeros(T);
    % t_backslash_full = zeros(T);
    % t_inverse_full = zeros(T);
    t_LU = 0;
    t_backslash = 0;
    t_backslash_full = 0;
    t_inverse_full = 0;


    % cumulate runtime
    for i=1:T
        % LU method
        tic
        x_LU = a01e05invKn(b);
        t_LU = t_LU + toc;
        % matlab backslash method
        tic
        x_bs = Kn\b;
        t_backslash = t_backslash + toc;
        % matlab backslash full method
        tic
        x_bs_full = full(Kn)\b;
        t_backslash_full = t_backslash_full + toc;
        % matlab inverse full method
        tic
        x_inv = inv(full(Kn))*b;
        t_inverse_full = t_inverse_full + toc;
    end

    runtime(j,:) = [t_LU, t_backslash, t_backslash_full, t_inverse_full]/T;

    residuum_LU(j) = norm(Kn*x_LU-b);
    residuum_backslash(j) = norm(Kn*x_bs-b);
    residuum_backslash_full(j) = norm(Kn*x_bs_full-b);
    residuum_inverse_full(j) = norm(Kn*x_inv-b);

end

% graphical output
figure(1)
subplot(1,2,1)
loglog(Ns, runtime(:,1),'*')
hold on
loglog(Ns, runtime(:,2),'o')
loglog(Ns, runtime(:,3),'d')
loglog(Ns, runtime(:,4),'x')
hold off
legend('LU','matlab backslash','matlab backslash full','matlab inverse full', 'Location','northwest')
title('runtime of different methods')
ylabel('runtime [s]')
xlabel('N')

subplot(1,2,2)
loglog(Ns, residuum_LU,'*')
hold on
loglog(Ns, residuum_backslash,'o')
loglog(Ns, residuum_backslash_full,'d')
loglog(Ns, residuum_inverse_full,'x')
hold off
legend('LU','matlab backslash','matlab backslash full','matlab inverse full', 'Location','northwest')
title('residual of the solution with different methods')
ylabel('norm')
xlabel('N')

%% runtime for different N (size of the matrix Kn
% initialisation of necessary variables
n = round(logspace(1,3),0);
n_size = length(n);

t_LU = zeros(1,n_size);
t_backslash = zeros(1,n_size);
t_backslash_full = zeros(1,n_size);
t_inverse_full = zeros(1,n_size);

% doing the comparison for different order of magnitude of n
for i=1:n_size
    Kn = a01e04sparse(n(i));
    b = rand(n(i),1);
    % cumulate runtime
    for j=1:T
        % LU method
        tic
        x = a01e05invKn(b);
        t_LU(i) = t_LU(i) + toc;
        % matlab backslash method
        tic
        x = Kn\b;
        t_backslash(i) = t_backslash(i) + toc;
        % matlab backslash full method
        tic
        x = full(Kn)\b;
        t_backslash_full(i) = t_backslash_full(i) + toc;
        % matlab inverse full method
        tic
        x = inv(full(Kn))*b;
        t_inverse_full(i) = t_inverse_full(i) + toc;
    end
end
% averaging
t_LU = t_LU/T;
t_backslash = t_backslash/T;
t_backslash_full = t_backslash_full/T;
t_inverse_full = t_inverse_full/T;
% graphical output
figure(2)
loglog(n,t_LU)
hold on
loglog(n,t_backslash)
loglog(n,t_backslash_full)
loglog(n,t_inverse_full)
hold off
legend('LU','matlab backslash','matlab backslash full','matlab inverse full', 'Location','northwest')
ylabel('runtime [s]')
xlabel('n')
title('Cumulative runtime of different methods')

%% runtime behavior at infinity
% initialisation of necessary variables
n = round(logspace(1,5),0);
n_size = length(n);

t_LU = zeros(1,n_size);
t_backslash = zeros(1,n_size);


% doing the comparison for different order of magnitude of n
for i=1:n_size
    Kn = a01e04sparse(n(i));
    b = rand(n(i),1);
    % cumulate runtime
    for j=1:T
        % LU method
        tic
        x_LU = a01e05invKn(b);
        t_LU(i) = t_LU(i) + toc;
        % matlab backslash method
        tic
        x_bs = Kn\b;
        t_backslash(i) = t_backslash(i) + toc;
    end
end
%averaging
t_LU = t_LU/T;
% behavior at infinity
S = diff(log(t_LU))./diff(log(n));
s_LU = mean(S(35:n_size-1));
c_LU = sum((t_LU(10:n_size)).*((n(10:n_size)).^s_LU))/...
    sum((n(10:n_size).^2).^s_LU);
Y_LU = (c_LU*n.^s_LU);


% graphical output
figure(3)
loglog(n,t_LU)
hold on
% loglog(n,t_backslash)
loglog(n(10:n_size),Y_LU(10:n_size))
% loglog(n(10:n_size),Y_bs(10:n_size))
hold off
% legend('LU method')%,'matlab backslash')
title('runtime for the LU method')
ylabel('runtime [s]')
xlabel('N')
text(n(39),t_LU(44),'O(N)')
