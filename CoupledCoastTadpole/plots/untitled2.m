% wavelets 2_17_21

pmin = 512;
pmax = 1024;
dx = 62.5;

%% IC 5
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_5_v1_uniform_Kc0_1_1_xyw.mat')
% idx = sub2ind(size(fetch_matrix),y,x);
[~,eq_14_1]= wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_1',1);

%% uniform 5
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_5_v1_uniform_Kc0_1_end_xyw.mat')
[~,eq_14_u2]=wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_7',1);
%% wave 5
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_5_v1_wave_Kc0_1_end_xyw.mat')
[~,eq_14_u3]=wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_14',1);

%% IC 9
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_9_v1_uniform_Kc0_1_1_xyw.mat')
% idx = sub2ind(size(fetch_matrix),y,x);
[~,eq_14_1]= wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_1',1);
%% uniform 9
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_9_v1_uniform_Kc0_1_end_xyw.mat')
[~,eq_14_w2]=wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_166',1);
%% wave 9
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/ray_fetch_test/river_idx_9_v1_wave_Kc0_1_end_xyw.mat')
[~,eq_14_w3]=wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_248',1);

%%
% ligeia mare
% load('Ligeia.mat')
% wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_248',1)


%% wave and uniform
% load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_mwave_muniform_end_xyw.mat')
% wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_mwave_muniform_end',1)
