function [strengthloss,totalSloss_frac,test] = calcstrengthloss(lake,dt,tend,p,dx,fetch_on);
input = lake;
t = dt:dt:tend;
p.dt = dt;
p.dx = dx;
p.dy = dx;
Sinit = p.So*(p.dx^2)./(p.dxo^2);
strength_init = Sinit*ones(size(lake));
p.strength_init = Sinit;
strength_in = strength_init;
for i = 1:length(t)    
[lake,strength_out,~,dam_matrix,~] = coastal_erosion(input,fetch_on,strength_in,p,[],[],[]);
test(i) = sum(strength_out-strength_in,'all');%sum(dam_matrix,'all');
% figure()
% imagesc(dam_matrix); colormap gray
strength_in = strength_out;
input = lake;
end

strengthloss = strength_init-strength_out;
totalSloss_frac = sum(strengthloss,'all')./sum(Sinit.*ones(size(lake)),'all');
end