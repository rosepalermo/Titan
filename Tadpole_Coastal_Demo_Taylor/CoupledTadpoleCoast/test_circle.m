function [init,p] = test_circle(p);

theta = 0 : 0.1 : 2*pi;
radius = 5;
lakex = radius * cos(theta);
lakey = radius * sin(theta);

eps = 50;
dx = 0.05; dy = 0.05;
x = (min(lakex)-eps):dx:(max(lakex)+eps);
y = (min(lakey)-eps):dy:(max(lakey)+eps);
[X,Y] = meshgrid(x,y);
Xinon = reshape(X,[],1);
Yinon = reshape(Y,[],1);

[in, on] = inpoly([Xinon,Yinon]',[lakex;lakey]);
lake = in + on;
init = ~reshape(lake,length(y),length(x));

p.dx = dx; p.dy = dy; p.Nx = length(x); p.Ny = length(y);