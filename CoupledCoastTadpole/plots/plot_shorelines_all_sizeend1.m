% plot shoreline of wave, uniform, and mix on top of each other
figure()

% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_uniform_Kc0_15.mat');
plotICsl(p,g)
hold on

% uniform only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/size_final_1/rand_uniform_Kc0_015.mat')
plotisl(p,g,)
hold on

% high uniform low wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/size_final_1/size_final_1rand_lwave_huniform.mat')
plotisl(p,g,)

% med uniform med wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/size_final_1/size_final_1rand_mwave_muniform.mat')
plotisl(p,g,)

% low uniform high wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/size_final_1/size_final_1rand_hwave_luniform_2_.mat')
plotisl(p,g,)

% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/size_final_1/size_final_1rand_wave_Kc0_015.mat')
plotisl(p,g,)

legend('IC','uniform','uniform>wave','uniform=wave','uniform<wave','wave')
set(gca,'FontSize',14)
xlabel('xpts')
ylabel('ypts')
axis equal