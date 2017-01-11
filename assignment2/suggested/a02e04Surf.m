function a02e04Surf(h,dt)
% Assignment 02, Exercise 04
%
% generates three surface plots of the exact solution u 
% with t \in [0, 10] and x \in [-15, 15].
%
% dt is time stepsize in Time
% h  is mesh stepsize for the x coordinate
% Run, for example: a02e04Surf(0.1,0.2)

% Fill in unset optional values for meshgrid steps
switch nargin
    case 0
        h = 0.1;
        dt = 0.2;
    case 1
        dt = 0.2;
end

% Create the meshgrid
[T,X] = meshgrid(0:dt:10,-15:h:15); 

% For the plots use the following three sets of parameter values 
%   (a, b) \in {(0, 1), (1, 0), (1, 1)}.
AB = [0,1;1,0;1,1];
for i = 1:3
    % Set parameters
    a = AB(i,1);
    b = AB(i,2);
    % graphical output
    figure(i)
    U = a02e04WaveEq(T,X,a,b);
    surf(T,X,U);
    xlabel('t')
    ylabel('x')
    title(['Solution of the wave equation for (a,b)=(',num2str(a),',',num2str(b),')'])
end