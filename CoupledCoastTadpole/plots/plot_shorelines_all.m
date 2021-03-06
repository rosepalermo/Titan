% plot shoreline of wave, uniform, and mix on top of each other
figure()

% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_uniform_Kc0_15.mat');
plotICsl(p,g)
hold on

% uniform only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_uniform_Kc0_15.mat');
plotendsl(p,g)
hold on

% high uniform low wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_lwave_huniform.mat')
plotendsl(p,g)

% med uniform med wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_mwave_muniform.mat')
plotendsl(p,g)

% low uniform high wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_hwave_luniform_2_.mat')
plotendsl(p,g)

% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_wave_Kc0_15.mat');
plotendsl(p,g)

legend('IC','uniform','uniform>wave','uniform=wave','uniform<wave','wave')
set(gca,'FontSize',14)
xlabel('xpts')
ylabel('ypts')
axis equal