% calc fetch for all sl
%% uniform
% small
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015_sl_small.mat')
p.Kwave = p.Kuniform;
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015_sl_small_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;

% large
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015_sl.mat')
p.Kwave = p.Kuniform;
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015_sl_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;
%% wave
% small
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15_sl_small.mat')
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15_sl_small_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;
% large
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15_sl.mat')
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15_sl_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;
%% wave and uniform
% small
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform_sl_small.mat')
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform_sl_small_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;
% large
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform_sl.mat')
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform_sl_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')
clear all;
%% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_IC_sl.mat')
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_IC_sl_fetch.mat','x','y','lake','dam_matrix','wave_weight_matrix','fetch_matrix','indshoreline_ordered','cells2trash','p')