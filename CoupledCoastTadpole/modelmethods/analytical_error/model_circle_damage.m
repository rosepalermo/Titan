% analytical error test

%% calculate strength loss for each
% wave
load('circle_init.mat')
load('p_circle_init.mat')
p.Kwave = 0.001;
p.Ao = 1;%25*pi;
p.runname = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/wave_sf1p1');
p.doWaveErosion = 1;
p.size_final = 1.1;
lake = ~init;
p.dt = 1;
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);

% solution = Tadpole(init*10,p);

%% uniform
load('circle_init.mat')
load('p_circle_init.mat')
p.con8 = 0;
p.saveint = 100;
p.Kuniform = 0.04;
p.Ao = 25*pi;
p.runname = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/uniform_sf1p1');
p.doUniformErosion = 1;
p.size_final = 2;
solution = Tadpole(init*10,p);
figure()
imagesc(solution(:,:,end)); title('con8 = 0')

%% analyze output