% IC
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_15.mat');
[x,y,wave_weight] = calc_xyw(file_name,1,true);


% uniform t1
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1.mat');
[x,y,wave_weight] = calc_xyw(file_name,7,true);

% uniform t2
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1.mat');
[x,y,wave_weight] = calc_xyw(file_name,14,true);

% wave t1
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1.mat');
[x,y,wave_weight] = calc_xyw(file_name,166,true);
% wave t2
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1.mat');
[x,y,wave_weight] = calc_xyw(file_name,248,true);

% wave and uniform
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_mwave_muniform.mat');
[x,y,wave_weight] = calc_xyw(file_name,'end',true);