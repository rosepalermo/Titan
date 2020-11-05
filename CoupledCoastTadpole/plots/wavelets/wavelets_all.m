% wavelets 101620

pmin = 16;
pmax = 128;

% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_uniform_Kc0_1_1_xyw.mat')
idx = sub2ind(size(fetch_matrix),y,x);
wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_uniform_Kc0_1_1')

% uniform
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_uniform_Kc0_1_end_xyw.mat')
wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_uniform_Kc0_1_end')

% wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_wave_Kc0_1_end_xyw.mat')
wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/102920/river_wave_Kc0_1_end')