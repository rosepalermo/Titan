%this code organizes the shoreline coordinates for main and calculates fetch/wave weighting for each shoreline.
% Rose Palermo
% Dec 2018

%% load shoreline coordinates
addpath('D:\Titan\Modeling\AGU final folder')
addpath('C:\Users\Rose Palermo\Documents\GitHub\Titan2\Wavelets\Current Working folder')

%% t1v1 RED NOISE
load('wavet1v1.mat')

% t1v1 _ REDNOISE
SL_matrix(1).cord{1,1}(:,1) = ordered_sl_save{1,1}{1,1}(:,1);
SL_matrix(1).cord{1,1}(:,2) = ordered_sl_save{1,1}{1,1}(:,2);

%% t1v1 WAVE 

load('wavet1v1.mat')
runsplot = [1;3;5];

% t1v1 _ wavet30
SL_matrix(2).cord{1,1}(:,1) = ordered_sl_save{runsplot(2),1}{1,1}(:,1);
SL_matrix(2).cord{1,1}(:,2) = ordered_sl_save{runsplot(2),1}{1,1}(:,2);

% t1v1 _ wavet50
SL_matrix(3).cord{1,1}(:,1) = ordered_sl_save{runsplot(3),1}{1,1}(:,1);
SL_matrix(3).cord{1,1}(:,2) = ordered_sl_save{runsplot(3),1}{1,1}(:,2);



%% t1v1 UNIFORM

load('uniformt1v1.mat')
runsplot = [1; 10; 20];

shoreline = addidshoreline_cardonly(lake_save{runsplot(2),1},~lake_save{runsplot(2),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(2),1},shoreline);
SL_matrix(4).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(4).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

shoreline = addidshoreline_cardonly(lake_save{runsplot(3),1},~lake_save{runsplot(3),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(3),1},shoreline);
SL_matrix(5).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(5).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));


%% t3v1 INITIAL AFTER RIVERS

load('xy_beta1_6_v1.mat')

% t1v1 _ REDNOISE
SL_matrix(6).cord{1,1}(:,1) = x_0_5m_t2v1;
SL_matrix(6).cord{1,1}(:,2) = y_0_5m_t2v1;

%% t1v1 WAVE 

%% LAKE POWELL
filename = 'LakePowell_gp.csv';
M = csvread(filename);
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);
lat = 110.978;
lon = 89.012;
x0 = (x0-min(x0)) * lon;%changed from lon/lat*lon to just lon
y0 = (y0-min(y0)) * lat;

dist = sqrt((x0(2:end) - x0(1:end-1)).^2 +  (y0(2:end) - y0(1:end-1)).^2);
dist_off = find(dist>1);
x0_fixed = [x0(1:dist_off(1)-1);x0(dist_off(2)+1:dist_off(3)-1);x0(dist_off(1)+1:dist_off(2)-1);x0(dist_off(12)+1:dist_off(13)-1);x0(dist_off(13)+1:end);x0(dist_off(11)+1:dist_off(12)-1);x0(dist_off(10)+1:dist_off(11)-1);x0(dist_off(6)+1:dist_off(7)-1);x0(dist_off(9)+1:dist_off(10)-1);x0(dist_off(3)+1:dist_off(4)-1);x0(dist_off(7)+1:dist_off(8)-1);x0(dist_off(5)+1:dist_off(6)-1);x0(dist_off(4)+1:dist_off(5)-1)];
y0_fixed = [y0(1:dist_off(1)-1);y0(dist_off(2)+1:dist_off(3)-1);y0(dist_off(1)+1:dist_off(2)-1);y0(dist_off(12)+1:dist_off(13)-1);y0(dist_off(13)+1:end);y0(dist_off(11)+1:dist_off(12)-1);y0(dist_off(10)+1:dist_off(11)-1);y0(dist_off(6)+1:dist_off(7)-1);y0(dist_off(9)+1:dist_off(10)-1);y0(dist_off(3)+1:dist_off(4)-1);y0(dist_off(7)+1:dist_off(8)-1);y0(dist_off(5)+1:dist_off(6)-1);y0(dist_off(4)+1:dist_off(5)-1)];
x0_fixed(1245) = []; y0_fixed(1245) = [];
x0 = x0_fixed*1000; y0 = y0_fixed*1000; % convert to meters

SL_matrix(7).cord{1,1}(:,1) = x0;
SL_matrix(7).cord{1,1}(:,2) = y0;

%% SCOTLAND

filename = 'Scotland_gp.csv';
M = csvread(filename);
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);
lat = 111360;
lon = 60772;
x0 = (x0-min(x0)) * lon;%changed from lon/lat*lon to just lon
y0 = (y0-min(y0)) * lat;
dist = sqrt((x0(2:end) - x0(1:end-1)).^2 +  (y0(2:end) - y0(1:end-1)).^2);
dist_off = find(dist>1000);
x0_fixed = [flipud(x0(924:2345));x0(1:923);x0(2346:end)];
y0_fixed = [flipud(y0(924:2345));y0(1:923);y0(2346:end)];
plot(x0_fixed,y0_fixed)
x0 = x0_fixed; y0 = y0_fixed;

SL_matrix(8).cord{1,1}(:,1) = x0;
SL_matrix(8).cord{1,1}(:,2) = y0;


%% LIGEIA MARE
filename = 'lg_4wavelets.xls';
M = xlsread(filename);

% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/lg_4wavelets_wavelet_updated'; 
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);

SL_matrix(9).cord{1,1}(:,1) = x0;
SL_matrix(9).cord{1,1}(:,2) = y0;

%% LAKE SEBAGO
M = xlsread('sebago_2_xy_UTM.xls');
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);
x0 = [x0(3319+25:3459);x0(2223:2787);x0(3460:4155);x0(5442:7123);x0(87:1664);x0(1941:2222);x0(4156:5044);x0(8329:8703);x0(7124:8328);x0(5381:5440);x0(5045:5379)];
y0 = [y0(3319+25:3459);y0(2223:2787);y0(3460:4155);y0(5442:7123);y0(87:1664);y0(1941:2222);y0(4156:5044);y0(8329:8703);y0(7124:8328);y0(5381:5440);y0(5045:5379)];
x_island = [x0(1665:1940);x0(2788:3318)];
y_island = [y0(1665:1940);y0(2788:3318)];
SL_matrix(10).cord{1,1}(:,1) = x0;
SL_matrix(10).cord{1,1}(:,2) = y0;

%% t2v1 WAVE 

load('wavet2v1.mat')
runsplot = [1;30;50];

% t1v1 _ wavet30
SL_matrix(11).cord{1,1}(:,1) = ordered_sl_save{runsplot(2),1}{1,1}(:,1);
SL_matrix(11).cord{1,1}(:,2) = ordered_sl_save{runsplot(2),1}{1,1}(:,2);

% t1v1 _ wavet50
SL_matrix(12).cord{1,1}(:,1) = ordered_sl_save{runsplot(3),1}{1,1}(:,1);
SL_matrix(12).cord{1,1}(:,2) = ordered_sl_save{runsplot(3),1}{1,1}(:,2);



%% t2v1 UNIFORM

load('uniformt2v1.mat')
runsplot = [1; 10; 20];

shoreline = addidshoreline_cardonly(lake_save{runsplot(2),1},~lake_save{runsplot(2),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(2),1},shoreline);
SL_matrix(13).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(13).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

shoreline = addidshoreline_cardonly(lake_save{runsplot(3),1},~lake_save{runsplot(3),1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(3),1},shoreline);
SL_matrix(14).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
SL_matrix(14).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

%% CALCULATE FETCH AND WAVE WEIGHTING FOR EACH SHORELINE AND SAVE
% 
% 
for i = 1:length(SL_matrix)
    i
    [WaveArea,FetchArea] = fetch_wavefield_cell(SL_matrix(i).cord);
    WaveArea_save{i} = WaveArea;
    FetchArea_save{i} = FetchArea;
    figure();scatter3(SL_matrix(i).cord{1,1}(:,1),SL_matrix(i).cord{1,1}(:,2),WaveArea{1,1},[],WaveArea{1,1});view(2);colorbar
end

% 
% for i = 1:length(SL_matrix)
%     figure()
%     plot(SL_matrix(i).cord{1,1}(:,1),SL_matrix(i).cord{1,1}(:,2))
% end
% save('waveandfetch_4AGU18.mat','SL_matrix','WaveArea_save','FetchArea_save')
save('shorelines_4AGU18.mat','SL_matrix')