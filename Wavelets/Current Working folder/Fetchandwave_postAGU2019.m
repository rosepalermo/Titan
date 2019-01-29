%this code organizes the shoreline coordinates for main and calculates fetch/wave weighting for each shoreline.
% Shorelines are first eroded by rivers, then waves
% Rose Palermo
% Jan 2019

%% load shoreline coordinates
addpath('D:\Titan\Modeling\river_and_wave_1_2019')
addpath('C:\Users\Rose Palermo\Documents\GitHub\Titan2\Wavelets\Current Working folder')

%% t2v1 WAVE 

load('wavet2v1.mat')
runsplot = [1;30;50];

% t1v1 _ wavet30
SL_matrix(1).cord{1,1}(:,1) = ordered_sl_save{runsplot(2),1}{1,1}(:,1);
SL_matrix(1).cord{1,1}(:,2) = ordered_sl_save{runsplot(2),1}{1,1}(:,2);

% t1v1 _ wavet50
SL_matrix(2).cord{1,1}(:,1) = ordered_sl_save{runsplot(3),1}{1,1}(:,1);
SL_matrix(2).cord{1,1}(:,2) = ordered_sl_save{runsplot(3),1}{1,1}(:,2);



%% t1v1 UNIFORM

load('uniformt1v1.mat')
runsplot = [1; 10; 20];

shoreline = addidshoreline_cardonly(lake_save{runsplot(2),1},~lake_save{runsplot(2),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(2),1},shoreline);
SL_matrix(3).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(3).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

shoreline = addidshoreline_cardonly(lake_save{runsplot(3),1},~lake_save{runsplot(3),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(3),1},shoreline);
SL_matrix(4).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(4).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));


%% t1v1 WAVE 

load('wavet1v1.mat')
runsplot = [1;30;50];

% t1v1 _ wavet30
SL_matrix(5).cord{1,1}(:,1) = ordered_sl_save{runsplot(2),1}{1,1}(:,1);
SL_matrix(5).cord{1,1}(:,2) = ordered_sl_save{runsplot(2),1}{1,1}(:,2);

% t1v1 _ wavet50
SL_matrix(6).cord{1,1}(:,1) = ordered_sl_save{runsplot(3),1}{1,1}(:,1);
SL_matrix(6).cord{1,1}(:,2) = ordered_sl_save{runsplot(3),1}{1,1}(:,2);



%% t2v1 UNIFORM

load('uniformt2v1.mat')
runsplot = [1; 10; 20];

shoreline = addidshoreline_cardonly(lake_save{runsplot(2),1},~lake_save{runsplot(2),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(2),1},shoreline);
SL_matrix(7).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(7).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

shoreline = addidshoreline_cardonly(lake_save{runsplot(3),1},~lake_save{runsplot(3),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(3),1},shoreline);
SL_matrix(8).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(8).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
%% t2v2 WAVE 

load('wavet2v2.mat')
runsplot = [1;30;50];

% t1v2 _ wavet30
SL_matrix(9).cord{1,1}(:,1) = ordered_sl_save{runsplot(2),1}{1,1}(:,1);
SL_matrix(9).cord{1,1}(:,2) = ordered_sl_save{runsplot(2),1}{1,1}(:,2);

% t1v2 _ wavet50
SL_matrix(10).cord{1,1}(:,1) = ordered_sl_save{runsplot(3),1}{1,1}(:,1);
SL_matrix(10).cord{1,1}(:,2) = ordered_sl_save{runsplot(3),1}{1,1}(:,2);



%% t1v2 UNIFORM

load('uniformt1v2.mat')
runsplot = [1; 10; 20];

shoreline = addidshoreline_cardonly(lake_save{runsplot(2),1},~lake_save{runsplot(2),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(2),1},shoreline);
SL_matrix(11).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(11).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

shoreline = addidshoreline_cardonly(lake_save{runsplot(3),1},~lake_save{runsplot(3),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(3),1},shoreline);
SL_matrix(12).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(12).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));


%% t1v2 WAVE 

load('wavet1v2.mat')
runsplot = [1;30;50];

% t1v2 _ wavet30
SL_matrix(13).cord{1,1}(:,1) = ordered_sl_save{runsplot(2),1}{1,1}(:,1);
SL_matrix(13).cord{1,1}(:,2) = ordered_sl_save{runsplot(2),1}{1,1}(:,2);

% t1v2 _ wavet50
SL_matrix(14).cord{1,1}(:,1) = ordered_sl_save{runsplot(3),1}{1,1}(:,1);
SL_matrix(14).cord{1,1}(:,2) = ordered_sl_save{runsplot(3),1}{1,1}(:,2);



%% t2v2 UNIFORM

load('uniformt2v2.mat')
runsplot = [1; 10; 20];

shoreline = addidshoreline_cardonly(lake_save{runsplot(2),1},~lake_save{runsplot(2),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(2),1},shoreline);
SL_matrix(15).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(15).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

shoreline = addidshoreline_cardonly(lake_save{runsplot(3),1},~lake_save{runsplot(3),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(3),1},shoreline);
SL_matrix(16).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(16).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

%% CALCULATE FETCH AND WAVE WEIGHTING FOR EACH SHORELINE AND SAVE
% 
% 
% for i = 1:length(SL_matrix)
%     i
%     [WaveArea,FetchArea] = fetch_wavefield_cell(SL_matrix(i).cord);
%     WaveArea_save{i} = WaveArea;
%     FetchArea_save{i} = FetchArea;
%     figure();scatter3(SL_matrix(i).cord{1,1}(:,1),SL_matrix(i).cord{1,1}(:,2),WaveArea{1,1},[],WaveArea{1,1});view(2);colorbar
% end

figure()
for i = 1:8
    plot(SL_matrix(i).cord{1,1}(:,1),SL_matrix(i).cord{1,1}(:,2))
    hold on
    axis equal
end
figure()
for i = 9:16
    plot(SL_matrix(i).cord{1,1}(:,1),SL_matrix(i).cord{1,1}(:,2))
    hold on
    axis equal
end
% save('waveandfetch_4AGU18.mat','SL_matrix','WaveArea_save','FetchArea_save')
% save('shorelines_4AGU18.mat','SL_matrix')