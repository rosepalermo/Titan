% shorelines for wavelets
close all; clear all
%% uniform only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015.mat')

lake = g.output(:,:,166)<=p.sealevel_init;
figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{1};
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
hold on; scatter(x,y)
%% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15.mat')

lake = g.output(:,:,139)<=p.sealevel_init;
figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{1};
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
%% uniform and wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_mwave_muniform.mat')

lake = g.output(:,:,85)<=p.sealevel_init;
figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{2};
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
%% initial condition
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15.mat')

lake = g.output(:,:,1)<=p.sealevel_init;
figure()
imagesc(lake)
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
lake = F_lake_all{2};
[indshoreline_ordered] = get_ordered_sl(lake,p);
[y,x] = ind2sub(size(lake),indshoreline_ordered);
%% figure

