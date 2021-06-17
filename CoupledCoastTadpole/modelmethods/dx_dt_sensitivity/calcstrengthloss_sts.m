function [strengthloss,totalSloss_frac] = calcstrengthloss_sts(lake,tend,p,sts,fetch_on);
strength_in = ones(size(lake));
input = lake;
t = 0:dt:tend;
p.dt = dt;
p.dx = dx;
p.dy = dx;
Sinit = p.So*p.dx./p.dx;
p.sts = sts;
for i = 1:length(t)
    
[lake,strength,~,~,~] = coastal_erosion(input,fetch_on,strength_in,p,[],[],[]);
input = lake;
strength_in = strength;
end

strengthloss = ones(size(lake))-strength;
totalSloss_frac = sum(strengthloss,'all')./sum(Sinit.*ones(size(lake)),'all');
end