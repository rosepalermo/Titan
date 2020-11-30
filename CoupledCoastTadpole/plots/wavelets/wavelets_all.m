% wavelets 101620

pmin = 512;
pmax = 1024;
dx = 62.5;

%% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_1_xyw.mat')
% idx = sub2ind(size(fetch_matrix),y,x);
wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_1',1)

%% uniform t1
% load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_7_xyw.mat')
% wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_7',1)
%% uniform t2
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_14_xyw.mat')
wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1_14',1)

%% wave
% load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_166_xyw.mat')
% wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_166',1)
%% wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_248_xyw.mat')
wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_248',1)

%%
% ligeia mare
load('Ligeia.mat')
% wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_248',1)


%% wave and uniform
% load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_mwave_muniform_end_xyw.mat')
% wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_mwave_muniform_end',1)
