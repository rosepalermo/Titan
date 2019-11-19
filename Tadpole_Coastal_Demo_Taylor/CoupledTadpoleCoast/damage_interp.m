close all; clear all

load('initialtopo.mat')

figure()
x_all = [NaN];
y_all = [NaN];
dam_all = [NaN];
% iterate over 0.1 m and calculate damage for uniform erosion
for SL = 10
    % define the shoreline
    lake = init<SL; 
    [damage] = addidshoreline(lake,~lake); % corners and edges
    indshoreline = find(damage);
    dam = damage(indshoreline);
    [y_ind,x_ind] = ind2sub(size(lake),indshoreline);
    scatter3(x_ind,y_ind,dam,[],dam)
    x_all = [x_all;x_ind];
    y_all = [y_all;y_ind];
    dam_all = [dam_all;dam];
    view(2)
    hold on
end
caxis([0 10])
axis equal
x_all = x_all(2:end); y_all = y_all(2:end); 


%%
SL1 = 1; SL2 = 10;

lake1 = init<SL1;
lake2 = init<SL2;

[damage1] = addidshoreline(lake1,~lake1); % corners and edges
[damage2] = addidshoreline(lake2,~lake2); % corners and edges
indshoreline1 = find(damage1);
dam1 = damage1(indshoreline1);
indshoreline2 = find(damage2);
dam2 = damage2(indshoreline2);
[y_ind1,x_ind1] = ind2sub(size(lake1),indshoreline1);
[y_ind2,x_ind2] = ind2sub(size(lake2),indshoreline2);

% figure
% scatter3(x_ind1,y_ind1,dam1,[],dam1)
% view(2)
% hold on
% scatter3(x_ind2,y_ind2,dam2,[],dam2)

sl2 = damage2;
sl2(sl2>0) = 1;
intermediate = lake2 - lake1 + sl2;
% figure()
% imagesc(intermediate)
[y_indint,x_indint] = ind2sub(size(intermediate),find(intermediate));

daminterp = damage1 + damage2;
% daminterp(daminterp=0)=NaN;
[y_daminterp,x_daminterp] = ind2sub(size(daminterp),find(daminterp));

x = 1:400; y = x;
[X,Y] = meshgrid(x,y);
Vq = griddata(x_daminterp,y_daminterp,daminterp(find(daminterp)),x_indint,y_indint);

figure
scatter3(x_indint,y_indint,Vq,[],Vq)
view(2)
caxis([0 10])
axis equal

figure(); imagesc(init); caxis([SL1 SL2])
axis equal
set(gca,'Ydir','normal')