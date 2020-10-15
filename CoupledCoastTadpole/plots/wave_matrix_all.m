% shorelines for wavelets
close all; clear all
%% uniform only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015_sl_small_fetch.mat')

figure()
imagesc(fetch_matrix)
% set(gca,'Clim',[0 0.3])
colorbar
% hold on; scatter(x,y)
title('Uniform')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015_sl_small.mat','x','y','lake','p')
%% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15_sl_small_fetch.mat')

figure()
imagesc(fetch_matrix)
% set(gca,'Clim',[0 0.3])
colorbar
% hold on; scatter(x,y)
title('Wave')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15_sl_small.mat','x','y','lake','p')
%% uniform and wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform_sl_small.mat')

figure()
imagesc(fetch_matrix)
% set(gca,'Clim',[0 0.3])
colorbar
% hold on; scatter(x,y)
title('W & U')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform_sl_small_fetch.mat','x','y','lake','p')
%% initial condition
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_IC_sl.mat')

figure()
imagesc(fetch_matrix)
% set(gca,'Clim',[0 0.3])
colorbar
% hold on; scatter(x,y)
title('IC')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_IC_sl_fetch.mat','x','y','lake','p')
%% figure

