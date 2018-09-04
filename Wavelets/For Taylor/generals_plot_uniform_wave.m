load('uniform_rednoise.mat')
runsplot = [1:5:35];
runsplot = [1;30;60];
h = figure();
h.Position = [440,378,974,420];
color = {'k';[0.3 0.3 0.3];[0.7 0.7 0.7]};

for i = 1:3

shoreline = addidshoreline_cardonly(lake_save{runsplot(i),1},~lake_save{runsplot(i),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(i),1},shoreline);
xuniform = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
yuniform = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
subaxis(1,2,1,'Spacing',0.05,'Margin',0.05)
hold on
plot(xuniform,yuniform,'Color',color{i})
axis equal tight
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
subaxis(1,2,2,'Spacing',0.05,'Margin',0.05)
hold on
plot(xuniform,yuniform,'Color',color{i})
axis equal
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'XLim',([450 750])); set(gca,'YLim',([600 850]))

end



hh = figure();
hh.Position = [440,378,974,420];
load('wave_rednoise.mat')
runsplot = [1;7;14];
for i = 1:3
x = ordered_sl_save{runsplot(i),1}{1,1}(:,1);
y = ordered_sl_save{runsplot(i),1}{1,1}(:,2);
subaxis(1,2,1,'Spacing',0.05,'Margin',0.05)
hold on
plot(x,y,'Color',color{i})
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
axis equal tight

subaxis(1,2,2,'Spacing',0.05,'Margin',0.05)
hold on
plot(x,y,'Color',color{i})
axis equal
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'XLim',([450 750])); set(gca,'YLim',([600 850]))
end
