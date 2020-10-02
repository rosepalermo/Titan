% make plots fetch dist through time for wave and uniform

%% uniform only 0.5 size final
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_uniform_Kc0_15.mat');

figure()
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_uniform_Kc0_15_fetch');
make_video_fetch_hist(v,p,g)

%% wave only 0.5 size final
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_wave_Kc0_15.mat');

figure()
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_wave_Kc0_15_fetch');
make_video_fetch_hist(v,p,g)

%% med wave med uniform 0.5 size final
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_mwave_muniform.mat')

figure()
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_mwave_muniform_fetch');
make_video_fetch_hist(v,p,g)