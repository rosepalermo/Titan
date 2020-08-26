function [init,p] = test_square(p);

side_length = 10;
lakex = [linspace(0,side_length,100) side_length*ones(1,100) linspace(side_length,0,100) zeros(1,100)];
lakey = [zeros(1,100) linspace(0,side_length,100) side_length*ones(1,100) linspace(side_length,0,100)];

eps = 50;
dx = 0.05; dy = 0.05;
x = (min(lakex)-eps):dx:(max(lakex)+eps);
y = (min(lakey)-eps):dy:(max(lakey)+eps);
[X,Y] = meshgrid(x,y);
Xinon = reshape(X,[],1);
Yinon = reshape(Y,[],1);

[in, on] = inpoly([Xinon,Yinon]',[lakex;lakey]);
lake = in + on;
lake = reshape(lake,length(y),length(x));
init = double(~lake);

p.dx = dx; p.dy = dy; p.Nx = length(x); p.Ny = length(y);


    p.strength = 1;
    p.So = 1;
    p.dxo =0.05;
    p.Ao = 100;
    p.Ao_cells = length(find(lake));
end