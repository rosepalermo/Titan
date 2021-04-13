% analogue analysis
%% Croatia1
p.Ao = 1; p.con8 = 1; p.Kwave = 1; p.dx = 1; p.dy = 1;p.Nx = size(lake,2); p.Ny = size(lake,1);p.nrays = 180; p.delta = 0.05;p.dt = 1; p.dxo = 1;p.So = 1;
file_name = ['/Users/rosepalermo/Documents/GitHub/Titan/CoupledCoastTadpole/analogues/Croatia1v2_sl.mat'];
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,1,true,lake,p);

%% FlCrookedLAke_v1
p.Ao = 1; p.con8 = 1; p.Kwave = 1; p.dx = 1; p.dy = 1;p.Nx = size(lake,2); p.Ny = size(lake,1);p.nrays = 180; p.delta = 0.05;p.dt = 1; p.dxo = 1;p.So = 1;
file_name = ['/Users/rosepalermo/Documents/GitHub/Titan/CoupledCoastTadpole/analogues/FLCrookedLakev2_sl.mat'];
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,1,true,lake,p);