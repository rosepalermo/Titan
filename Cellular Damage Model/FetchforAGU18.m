% % load shoreline coordinates
addpath('D:\Titan\Modeling\AGU final folder')
%Ligeia Mare
% filename = 'lg_4wavelets.xls';
% M = xlsread(filename);
% 
% %get rid of duplicate points (when it goes exactly around a pixel)
% % indices to unique values in column 3
% [~, ind] = unique(M(:, 4:5), 'rows');
% % duplicate indices
% duplicate_ind = setdiff(1:size(M, 1), ind);
% % duplicate values
% duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
% M(duplicate_ind,:)=[];
% 
% x0 = M(:,4);
% y0 = M(:,5);
% x0 = x0(1:2:end);
% y0 = y0(1:2:end);
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

% %% t2v1 WAVE 
% 
% load('wavet2v1.mat')
% runsplot = [1;30;50];
% 
% % t1v1 _ wavet30
% SL_matrix(7).cord{1,1}(:,1) = ordered_sl_save{runsplot(2),1}{1,1}(:,1);
% SL_matrix(7).cord{1,1}(:,2) = ordered_sl_save{runsplot(2),1}{1,1}(:,2);
% 
% % t1v1 _ wavet50
% SL_matrix(8).cord{1,1}(:,1) = ordered_sl_save{runsplot(3),1}{1,1}(:,1);
% SL_matrix(8).cord{1,1}(:,2) = ordered_sl_save{runsplot(3),1}{1,1}(:,2);
% 
% 
% 
% %% t2v1 UNIFORM
% 
% load('uniformt2v1.mat')
% runsplot = [1; 10; 20];
% 
% shoreline = addidshoreline_cardonly(lake_save{runsplot(2),1},~lake_save{runsplot(2),1});
% [sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(2),1},shoreline);
% SL_matrix(9).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
% SL_matrix(9).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
% 
% shoreline = addidshoreline_cardonly(lake_save{runsplot(3),1},~lake_save{runsplot(3),1});
% [sl_cell,cells2trash] = order_cw_lastpoint(lake_save{runsplot(3),1},shoreline);
% SL_matrix(10).cord{1,1}(:,1) = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
% SL_matrix(10).cord{1,1}(:,2) = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));

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