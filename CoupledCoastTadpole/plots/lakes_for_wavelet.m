% shorelines for wavelets
close all; clear all
%% uniform only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015.mat')

% lake = g.output(:,:,166)<=p.sealevel_init; % larger size
lake = g.output(:,:,105)<=p.sealevel_init; % smaller size
figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{2}; % 2 smaller size, 1 larger size
p.lake_only = 1;
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
hold on; scatter(x,y)
title('Uniform')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015_sl_small.mat','x','y','lake','p')
%% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15.mat')

% lake = g.output(:,:,139)<=p.sealevel_init;% larger size
lake = g.output(:,:,99)<=p.sealevel_init;% smaller size

figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{1};
p.lake_only = 1;
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
hold on; scatter(x,y)
title('Wave')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15_sl_small.mat','x','y','lake','p')
%% uniform and wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform.mat')

% lake = g.output(:,:,85)<=p.sealevel_init; %larger size
lake = g.output(:,:,59)<=p.sealevel_init; %smaller size

figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{2};
p.lake_only = 1;
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
hold on; scatter(x,y)
title('W & U')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform_sl_small.mat','x','y','lake','p')
%% initial condition
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15.mat')

lake = g.output(:,:,1)<=p.sealevel_init;
figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{2};
p.lake_only = 1;
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
hold on; scatter(x,y)
title('IC')
save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_IC_sl.mat','x','y','lake','p')
%% figure

