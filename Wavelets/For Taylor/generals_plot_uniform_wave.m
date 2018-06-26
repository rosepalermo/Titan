load('uniform_rednoise.mat')
runsplot = [1:5:35];
runsplot = [1:10:70];
figure()
% subplot(1,2,1)
hold on
for i = 1:7
shoreline = addidshoreline_cardonly(lake_save{runsplot(i),1},~lake_save{runsplot(i),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(i),1},shoreline);
xuniform = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
yuniform = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
plot(xuniform,yuniform)
end
axis equal tight
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
