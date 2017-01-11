function a02e05Surf(h,dt)
% Assignment 02, Exercise 05
%
% generates a surface plot of the exact solution u 
% with t \in [0, 10] and x \in [-\pi, \pi].
% 
% dt is time stepsize in Time
% h  is mesh stepsize for the x coordinate
% Run, for example: a02e05Surf(0.1,0.2)

% Fill in unset optional values for meshgrid steps
switch nargin
    case 0
        h = 0.1;
        dt = 0.2;
    case 1
        dt = 0.2;
end

% Create the meshgrid
[T,X] = meshgrid(0:dt:10,-pi:h:pi); 

U = a02e05HeatEq(T,X);
surf(T,X,U);
title('Solution of the heat equation with initial condition')
xlabel('t')
ylabel('x')
zlabel('u(t,x)')

