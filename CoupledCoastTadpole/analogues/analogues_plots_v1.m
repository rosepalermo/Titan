% wavelets 101620

pmin = 1040;
pmax = 2048;
dx = 30;

%% Croatia1
load('Croatia1v2_sl_1_xyw.mat')
% idx = sub2ind(size(fetch_matrix),y,x);
dx=1;
[~,eq_14_cro1]= wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/Croatia1_',1);

%% Crooked Lake
load('FLCrookedLakev2_sl_1_xyw.mat')
% idx = sub2ind(size(fetch_matrix),y,x);
dx=1;
[~,eq_14_cro1]= wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/FLCrookedLakev2_',1);

%% Fort Peck Lake
dx = 30;
load('FortPeckLake_sl_1_xyw.mat')
[~,eq_14_fpl]= wavelets(x*dx,y*dx,wave_weight,pmin,pmax,['/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/FortPeckLake_'],1);
%% Sebago
dx = 30;
load('Sebago_sl_1_xyw.mat')
[~,eq_14_seb]=wavelets(x*dx,y*dx,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/Sebago_',1);

%%
% ligeia mare
% load('Ligeia.mat')
% wavelets(x,y,wave_weight,pmin,pmax,'/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1_248',1)
