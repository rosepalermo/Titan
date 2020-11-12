% wavelets 101620

pmin = 16;
pmax = 32;

% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_uniform_Kc0_1_1_xyw.mat')
idx = sub2ind(size(fetch_matrix),y,x);
wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_uniform_Kc0_1_1')

% % uniform t1
% load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_uniform_Kc0_1_2_xyw.mat')
% wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_uniform_Kc0_1_2')
% uniform t2
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_uniform_Kc0_1_3_xyw.mat')
wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_uniform_Kc0_1_3')

% % wave
% load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_wave_Kc0_15_12_xyw.mat')
% wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_wave_Kc0_15_12')
% wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_wave_Kc0_15_end_xyw.mat')
wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/110920/river_wave_Kc0_15_end')