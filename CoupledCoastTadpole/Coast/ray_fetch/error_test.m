% test error and time for rayfetch vs fetch_vis_approx

% load 
load('fetch_input8con_lake','F_lake'); % lake is 1 for water, 0 for land
load('fetch_input8con','fetch_sl_cells'); % x and y are ordered clockwise if i increases downward, first point != last point

% visilibity run
tic
eps = 0.1;
epsilon = 1e-4;
snap_distance = 0.05; % Looks like we get some invalid (empty array) visibility polygons if snap_distance >= eps.

[WA_vislib,FA_vislib] = fetch_VisiLibity(fetch_sl_cells,eps,epsilon,snap_distance);
time_vis = toc;
WA_vislib_ = cell2mat(WA_vislib);

cellsize = 62.5;
nrays = 36:360; % number of rays around the circle. Recommend nrays >= 36
delta = 0.002:0.001:0.5; % distance increment along each ray, in cells. Recommend delta <= 0.5
nrays_test = 180;
delta_test = 0.05;
time_rays = zeros(size(nrays));
time_delta = zeros(size(delta));

% % loop over rays
for n = 1:length(nrays)
    tic
    [warea,~,~,~,~,~,~,~] = GetFetchArea(fetch_sl_cells,F_lake,nrays(n),delta_test,cellsize);
    time_rays(n) = toc;
    warea_rays_save{n} = warea;
    diff_rays{n} = WA_vislib_ - cell2mat(warea_rays_save{n});
    sum_diff_rays(n) = sum(diff_rays{n});
    sum_abs_diff_rays(n) = sum(abs(diff_rays{n}));
    sum_warea(n) = sum(cell2mat(warea_rays_save{n}),'all');
    diff_vis_rays(n) = sum_warea(n)-sum(WA_vislib_);
end

% loop over delta
for n = 1:length(delta)
    tic
    [warea,~,~,~,~,~,~,~] = GetFetchArea(fetch_sl_cells,F_lake,nrays_test,delta(n),cellsize);
    time_delta(n) = toc;
    warea_delta_save{n} = warea;
    diff_delta{n} = WA_vislib_ - cell2mat(warea_delta_save{n});
    sum_diff_delta(n) = sum(diff_delta{n});
    sum_abs_diff_delta(n) = sum(abs(diff_delta{n}));
    sum_warea_delta(n) = sum(cell2mat(warea_delta_save{n}),'all');
    diff_vis_delta(n) = sum_warea_delta(n)-sum(WA_vislib_);
end

%% plots
% plot time
figure()
subplot(2,1,1)
plot(delta,time_delta,'k','LineWidth',2)
xlabel('delta')
ylabel('time (s)')
set(gca,'FontSize',12)

subplot(2,1,2)
plot(nrays,time_rays,'k','LineWidth',2)
xlabel('number of rays')
ylabel('time (s)')
set(gca,'FontSize',12)

% plot difference
figure()
subplot(2,1,1)
plot(delta,diff_vis_delta,'k','LineWidth',2)
xlabel('delta')
ylabel('sum of wave area difference between vis and ray')
set(gca,'FontSize',12)

subplot(2,1,2)
plot(nrays,diff_vis_rays,'k','LineWidth',2)
xlabel('number of rays')
ylabel('sum of wave area difference between vis and ray')
set(gca,'FontSize',12)

% % plot absolute difference
% figure()
% subplot(2,1,1)
% plot(delta,sum_abs_diff_delta,'k','LineWidth',2)
% xlabel('delta')
% ylabel('sum of abs (difference between vis and ray)')
% set(gca,'FontSize',12)
% 
% subplot(2,1,2)
% plot(nrays,sum_abs_diff_rays,'k','LineWidth',2)
% xlabel('number of rays')
% ylabel('sum of abs (difference between vis and ray)')
% set(gca,'FontSize',12)
%%
% plot warea
h = figure();
% subplot(2,1,1)
ax1 = axes(h);
plot(ax1,delta,sum_warea_delta./sum_warea_delta(1),'k','LineWidth',2)
xlabel('delta')
% ylabel('sum of wave weighting')
set(gca,'FontSize',16)
ylim([0.84 1])
hold on
ax2 = axes(h);
% subplot(2,1,2)
plot(ax2, nrays,sum_warea./sum_warea(end),'r','LineWidth',2)
xlabel('number of rays')
ylabel('relative error')
set(gca,'FontSize',16)
ylim([0.84 1])
ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
ax1.Box = 'off';
ax2.Box = 'off';
ax2.Color = 'none';
linkaxes([ax1 ax2],'y')
ax2.XColor = 'r';
% ax2.YColor = 'r';
%%
% plot warea slope
figure()
subplot(2,1,1)
plot(delta(1:end-1),diff(sum_warea_delta)/0.001,'k','LineWidth',2)
xlabel('delta')
ylabel('slope of wave weighting')
set(gca,'FontSize',12)

subplot(2,1,2)
plot(nrays(1:end-1),diff(sum_warea),'k','LineWidth',2)
xlabel('number of rays')
ylabel('slope of wave weighting')
set(gca,'FontSize',12)

warea = warea_delta_save{49};
warea = cell2mat(warea);
figure()
histogram(warea)
hold on
histogram(WA_vislib_)
xlabel('wave area')
ylabel('#')
legend('rayfetch','VisiLibity')
set(gca,'FontSize',16)
title('180 rays; 0.05 delta')

warea = warea_rays_save{1};
warea = cell2mat(warea);
figure()
histogram(warea)
hold on
histogram(WA_vislib_)
ylabel('#')
xlabel('wave area')
legend('rayfetch','VisiLibity')
set(gca,'FontSize',16)
title('36 rays; 0.05 delta')

warea = warea_rays_save{end};
warea = cell2mat(warea);
figure()
histogram(warea)
hold on
histogram(WA_vislib_)
xlabel('wave area')
ylabel('#')
legend('rayfetch','VisiLibity')
set(gca,'FontSize',16)
title('360 rays; 0.05 delta')
