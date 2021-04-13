% Butler
file_name = ['/Users/rosepalermo/Documents/GitHub/Titan/CoupledCoastTadpole/analogues/ButlerLake_sl.mat'];
load(file_name)
pmin = 128; pmax = 256;
p.Ao = 1; p.con8 = 1; p.Kwave = 1; p.dx = 30; p.dy = 30;p.Nx = size(lake,2); p.Ny = size(lake,1);p.nrays = 180; p.delta = 0.05;p.dt = 1; p.dxo = 1;p.So = 1;
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,1,true,lake,p);

[~,eq_14]= wavelets(x*p.dx,y*p.dx,wave_weight,pmin,pmax,['/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/Butler_'],1);

save_name = ['/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/Butler_results.mat'];
save(save_name,'x','y','wave_weight','eq_14')

%% Fort Peck Lake
file_name = ['/Users/rosepalermo/Documents/GitHub/Titan/CoupledCoastTadpole/analogues/FortPeckLake_sl.mat'];
load(file_name)
pmin = 128; pmax = 256;
p.Ao = 1; p.con8 = 1; p.Kwave = 1; p.dx = 30; p.dy = 30;p.Nx = size(lake,2); p.Ny = size(lake,1);p.nrays = 180; p.delta = 0.05;p.dt = 1; p.dxo = 1;p.So = 1;
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,1,true,lake,p);

[~,eq_14]= wavelets(x*p.dx,y*p.dx,wave_weight,pmin,pmax,['/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/FortPeckLake_'],1);

save_name = ['/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/FortPeckLake_results.mat'];
save(save_name,'x','y','wave_weight','eq_14')

%% Sebago
file_name = ['/Users/rosepalermo/Documents/GitHub/Titan/CoupledCoastTadpole/analogues/Sebago_sl.mat'];
load(file_name)
pmin = 128; pmax = 256;
p.Ao = 1; p.con8 = 1; p.Kwave = 1; p.dx = 30; p.dy = 30;p.Nx = size(lake,2); p.Ny = size(lake,1);p.nrays = 180; p.delta = 0.05;p.dt = 1; p.dxo = 1;p.So = 1;
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,1,true,lake,p);

[~,eq_14]= wavelets(x*p.dx,y*p.dx,wave_weight,pmin,pmax,['/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/Sebago_'],1);

save_name = ['/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/results/Sebago_results.mat'];
save(save_name,'x','y','wave_weight','eq_14')
