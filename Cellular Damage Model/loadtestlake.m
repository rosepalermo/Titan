load('xycontours.mat')
lakex = x_1m_t1;
lakey = y_1m_t1;
eps = 200;
dx = 5; dy = 5;
x = (min(lakex)-eps):dx:(max(lakex)+eps);
y = (min(lakey)-eps):dy:(max(lakey)+eps);
[X,Y] = meshgrid(x,y);
Xinon = reshape(X,[],1);
Yinon = reshape(Y,[],1);

[in, on] = inpolygon(Xinon,Yinon,lakex,lakey);
lake = in + on;
lake = reshape(lake,length(y),length(x));
land = ~lake;
[shoreline] = addidshoreline_cardonly(lake,land);
