% color = {'k';[0.3 0.3 0.3];[0.7 0.7 0.7]};
% color = {'k';'r';'b'};
addpath('C:\Users\Rose Palermo\Documents\Titan\Modeling\11_18_riverandwave_higherstrength')
addpath('C:\Users\Rose Palermo\Documents\GitHub\Titan2\Wavelets\Current Working folder')


hh = figure();
% hh.Position = [440,378,974,420];
load('wavet2v1.mat')
runsplot = [1;3;5];
runsplot = [1:5];

for i = 1:length(runsplot)
x = ordered_sl_save{runsplot(i),1}{1,1}(:,1);
y = ordered_sl_save{runsplot(i),1}{1,1}(:,2);
subaxis(2,2,1,'Spacing',0.05,'Margin',0.05)
hold on
plot(x,y)%,'Color',color{i})
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
axis equal

subaxis(2,2,2,'Spacing',0.05,'Margin',0.05)
hold on
plot(x,y)%,'Color',color{i})
axis equal
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'XLim',([2100 2400])); set(gca,'YLim',([1800 2000]))

% [period{i},global_Save{i}] = dowave(theta,deltad,n,x,y,savename,save_on,fetch,i);

end
%%
addpath('D:\Titan\Modeling\AGU final folder')
load('uniformt2v1.mat')
% runsplot = [1; 10; 20];
runsplot = [1:2:20];

% h = figure();
% h.Position = [440,378,974,420];

for ii = 1:length(runsplot)

shoreline = addidshoreline_cardonly(lake_save{runsplot(ii),1},~lake_save{runsplot(ii),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(ii),1},shoreline);
xuniform = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
yuniform = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
subaxis(2,2,3,'Spacing',0.05,'Margin',0.05)
hold on
plot(xuniform,yuniform)%,'Color',color{ii})
axis equal
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])

subaxis(2,2,4,'Spacing',0.05,'Margin',0.05)
hold on
plot(xuniform,yuniform)%,'Color',color{ii})
axis equal
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'XLim',([2100 2400])); set(gca,'YLim',([1800 2000]))

end