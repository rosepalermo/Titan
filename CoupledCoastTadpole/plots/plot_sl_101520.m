% plot all shorelines 10_15_20
colors_ = {'k','b','r'};
% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_uniform_Kc0_1.mat')
plotICsl(p,g,'k')
hold on

% uniform only
plotendsl(p,g,'b')
hold on

% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc0_15.mat')
plotendsl(p,g,'r')
hold on

legend('River IC','uniform','wave')