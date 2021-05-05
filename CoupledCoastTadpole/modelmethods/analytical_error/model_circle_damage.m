% models for analytical error test

%% run models
% wave
load('circle_init.mat')
load('p_circle_init.mat')
p.Kwave = 0.001;
p.Ao = 25*pi;
p.runname = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/wave_sf1p1');
p.doWaveErosion = 1;
p.size_final = 1.1;
solution = Tadpole(init*10,p);

% uniform
load('circle_init.mat')
load('p_circle_init.mat')
p.Kuniform = 0.001;
p.Ao = 25*pi;
p.runname = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/circle_analytical_error/uniform_sf1p1');
p.doUniformErosion = 1;
p.size_final = 1.1;
solution = Tadpole(init*10,p);

%% analyze output