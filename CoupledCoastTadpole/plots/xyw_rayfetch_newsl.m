% IC 5
file_name = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_5_v1_uniform_Kc0_1.mat';
[x,y,wave_weight] = calc_xyw(file_name,1,true);


% uniform 5
file_name = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_5_v1_uniform_Kc0_1.mat';
[x,y,wave_weight] = calc_xyw(file_name,'end',true);


% wave 5
file_name = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_5_v1_wave_Kc0_1.mat';
[x,y,wave_weight] = calc_xyw(file_name,'end',true);

% IC 9
file_name = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_9_v1_uniform_Kc0_1.mat';
[x,y,wave_weight] = calc_xyw(file_name,1,true);


% uniform 9
file_name = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_9_v1_uniform_Kc0_1.mat';
[x,y,wave_weight] = calc_xyw(file_name,'end',true);


% wave 9
file_name = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_9_v1_wave_Kc0_1.mat';
[x,y,wave_weight] = calc_xyw(file_name,'end',true);
% % wave and uniform
% file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/112920/river_mwave_muniform.mat');
% [x,y,wave_weight] = calc_xyw(file_name,14,true);