% IC
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_uniform_Kc0_1.mat');
[x,y,wave_weight] = calc_xyw(file_name,1,true);

% uniform
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_uniform_Kc0_1.mat');
[x,y,wave_weight] = calc_xyw(file_name,'end',true);

% wave
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_wave_Kc0_1.mat');
[x,y,wave_weight] = calc_xyw(file_name,'end',true);
