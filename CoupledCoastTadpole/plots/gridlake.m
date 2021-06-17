function [lake,X,Y] = gridlake(lakex,lakey,dx,dy,eps)

%make a grid larger than lake by eps
% eps = 5;
% dx = 0.05; dy = 0.05;
x = (min(lakex)-eps):dx:(max(lakex)+eps);
y = (min(lakey)-eps):dy:(max(lakey)+eps);
% minxy = min([min(lakex) min(lakey)]);
% maxxy = max([max(lakex) max(lakey)]);
% x = (minxy-eps):dx:(maxxy+eps);
% y = (minxy-eps):dy:(maxxy+eps);
[X,Y] = meshgrid(x,y);
Xinon = reshape(X,[],1);
Yinon = reshape(Y,[],1);


%points in and on the polygon are the lake
[in, on] = inpoly([Xinon,Yinon]',[lakex;lakey]);
lake = in + on;
lake(lake>1) =1;
lake = reshape(lake,length(y),length(x));
% figure()
% imagesc(x,y,lake)
% shading interp