
% the original shorelines I gave you were time steps 1, 3, 5 for wave and
% 1, 10, 20 for uniform. 1 is Rednoise.

load('wavet1v1.mat')
runsplot = [1:5];

% t1v1
for i = 1:length(runsplot)
    SL_matrix(i).cord{1,1}(:,1) = ordered_sl_save{runsplot(i),1}{1,1}(:,1);
    SL_matrix(i).cord{1,1}(:,2) = ordered_sl_save{runsplot(i),1}{1,1}(:,2);
end

save('wave_shorelines_all.mat')


%% t1v1 UNIFORM
clear all

load('uniformt1v1.mat')
runsplot = [1:20];

for i = 1:length(runsplot)
shoreline = addidshoreline_cardonly(lake_save{runsplot(i),1},~lake_save{runsplot(i),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(i),1},shoreline);
SL_matrix(i).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(i).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
end

save('uniform_shorelines_all.mat')
