% plot eq14 vs fetch

% load river IC
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/riverIC_sl_fetch_eq14.mat')
ind = sub2ind(size(lake),y,x);
y_limits = [0 nanmax(eq14holder)];
x_limits = [0 3.5e7];
figure()
scatter(fetch_matrix(ind),eq14holder,'k.')
hold on
ylabel('Roughness')
xlabel('Fetch')
set(gca,'Ylim',y_limits)
set(gca,'Xlim',x_limits)
set(gca,'FontSize',14)


% load wave
figure()
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc1_sl_fetch_eq14.mat')
ind = sub2ind(size(lake),y,x);
scatter(fetch_matrix(ind),eq14holder,'r.')
ylabel('Roughness')
xlabel('Fetch')
set(gca,'Ylim',y_limits)
set(gca,'Xlim',x_limits)
set(gca,'FontSize',14)

% load uniform
figure()
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_uniform_Kc0_1_sl_fetch_eq14.mat')
ind = sub2ind(size(lake),y,x);
scatter(fetch_matrix(ind),eq14holder,'b.')
set(gca,'Ylim',y_limits)
set(gca,'Xlim',x_limits)
set(gca,'FontSize',14)

ylabel('Roughness')
xlabel('Fetch')
% legend('river','wave','uniform')

%%
% UIST
load('uist.mat')
scatter(Fetch_matrix{1},eq_14,'k.')
% set(gca,'Ylim',y_limits)
% set(gca,'Xlim',x_limits)
set(gca,'FontSize',14)

ylabel('Roughness')
xlabel('Fetch')

%%
% LIGEIA
load('Ligeia.mat')
scatter(Fetch_matrix{1},eq_14,'k.')
% set(gca,'Ylim',y_limits)
% set(gca,'Xlim',x_limits)
set(gca,'FontSize',14)

ylabel('Roughness')
xlabel('Fetch')

%%
% Sebago
load('Sebago.mat')
scatter(Fetch_matrix{1},eq_14,'k.')
% set(gca,'Ylim',y_limits)
% set(gca,'Xlim',x_limits)
set(gca,'FontSize',14)

ylabel('Roughness')
xlabel('Fetch')