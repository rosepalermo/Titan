% calc fetch for all sl_river IC
%% uniform
% small
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_uniform_Kc0_1_sl.mat')
p.Kwave = p.Kuniform;
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_uniform_Kc0_1_sl_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;

%% wave
% small
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc1_sl.mat');
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc1_sl_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;

%% wave and uniform
% small

%% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/riverIC_sl.mat')
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/riverIC_sl_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')