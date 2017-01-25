f = @(x) x^3-6*x+((pi^2)/4+1)*cos(pi^2*x);
[x,U] = a10e03getPDE(-1,1,-1,1,f,10)