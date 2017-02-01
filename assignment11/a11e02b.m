h = 0.25;
x=[0:h:1];
y=[0:h:1];

[X,Y] = meshgrid(x,y);
tri = delaunay(X,Y);
triplot(tri,X,Y);
z = zeros(size(X));
z(1,1) = 1;
trisurf(tri,X,Y,z)