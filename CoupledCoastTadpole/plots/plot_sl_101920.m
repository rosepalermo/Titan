% plot all shorelines 10_19_20
colors_ = {'k','b','r'};
% IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc1.mat')
plotICsl(p,g,'k')
hold on
lake = g.output(:,:,1)<=p.sealevel_init;
% save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/riverIC_sl.mat','lake','g','x','y','p')

% uniform only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_uniform_Kc0_1.mat')
plotendsl(p,g,'b')
hold on
% save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_uniform_Kc0_1_sl.mat','lake','g','x','y','p')


% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc1.mat')
plotendsl(p,g,'r')
hold on
% save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc1_sl.mat','lake','g','x','y','p')

legend('River IC','uniform','wave')