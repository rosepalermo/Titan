load('uniform_rednoise.mat')
runsplot = [1:5:35];
runsplot = [1;30;60];
figure()
color = {'k';[0.3 0.3 0.3];[0.7 0.7 0.7]};

for i = 1:3

shoreline = addidshoreline_cardonly(lake_save{runsplot(i),1},~lake_save{runsplot(i),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(i),1},shoreline);
xuniform = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
yuniform = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
subplot(2,1,1)
hold on
plot(xuniform./1000,yuniform./1000,'Color',color{i})
axis equal tight
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
subplot(2,1,2)
hold on
plot(xuniform./1000,yuniform./1000,'Color',color{i})
axis equal tight
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
set(gca,'XLim',([0.45 0.75])); set(gca,'YLim',([0.6 0.9]))

end




% 
% subplot(1,2,2)
% hold on
% load('wave_rednoise.mat')
% runsplot = [1;3;5;7;9;11;13];
% runsplot = (1:7);
% for i = 1:7
% x = ordered_sl_save{runsplot(i),1}{1,1}(:,1);
% y = ordered_sl_save{runsplot(i),1}{1,1}(:,2);
% plot(x,y)
% end
% axis equal tight
