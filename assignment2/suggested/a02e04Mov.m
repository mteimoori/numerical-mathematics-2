function a02e04Mov(a,b,steps,h)
% Assignment 02, Exercise 04
%
% generates for given parameter values a, b \in R 
% an animation of the exact solution.
% Each frame should display a line plot of the exact solution 
% for a fixed t and for x \in [-15, 15]. 
% Your animation should show the evolution of the exact solution 
% over t \in [0, 10] with at least 100 frames.

% default number of time frames
switch nargin
    case 2
        steps = 100;
        h = 100;
    case 3
        h = 100; 
end

t     = linspace(0,10,steps);
x     = linspace(-15,15,h);
for i = 1:steps
    u = a02e04WaveEq(t(i),x,a,b);
    plot(x,u,'Color','r','LineWidth',3)
    xlabel('x')
    title('Exact solution of the heat eaquation in time')
    ylim([-1.5 1.5])
    pause(1/24)
end
