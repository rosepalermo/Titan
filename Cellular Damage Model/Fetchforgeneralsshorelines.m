% load shoreline coordinates
filename = 'lg_4wavelets.xls';
M = xlsread(filename);

%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0 = M(:,4);
y0 = M(:,5);
x0 = x0(1:2:end);
y0 = y0(1:2:end);

A(1).cord{1,1}(:,1)=x0;
A(1).cord{1,1}(:,2)=y0;

% REDNOISE
load('wave_rednoise.mat')
A(2).cord{1,1}(:,1) = ordered_sl_save{1,1}{1,1}(:,1);
A(2).cord{1,1}(:,2) = ordered_sl_save{1,1}{1,1}(:,2);

% WAVE MODIFIED
A(3).cord{1,1}(:,1) = ordered_sl_save{7,1}{1,1}(:,1);
A(3).cord{1,1}(:,2) = ordered_sl_save{7,1}{1,1}(:,2);


% RIVER MODIFIED
load('xycontours.mat')
A(4).cord{1,1}(:,1) = x_1m_t2(1:2:end)';
A(4).cord{1,1}(:,2) = y_1m_t2(1:2:end)';

% UNIFORM MODIFIED
load('uniform_rednoise.mat')
shoreline = addidshoreline_cardonly(lake_save{75,1},~lake_save{75,1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{75,1},shoreline);
A(5).cord{1,1}(:,1) = sl_cell{1,1}(:,1);
A(5).cord{1,1}(:,2) = sl_cell{1,1}(:,2);


for i = 1:length(A)
    i
    [WaveArea,FetchArea] = fetch_wavefield_cell(A(i).cord);
    WaveArea_save{i} = WaveArea;
    FetchArea_save{i} = FetchArea;
end