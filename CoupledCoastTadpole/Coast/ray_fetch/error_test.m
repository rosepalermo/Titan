% test error and time for rayfetch vs fetch_vis_approx

% load 
load('fetch_input8con_lake','F_lake'); % lake is 1 for water, 0 for land
load('fetch_input8con','fetch_sl_cells'); % x and y are ordered clockwise if i increases downward, first point != last point

% visilibity run
% tic
% [WA_vislib,FA_vislib] = fetch_vis_approx(fetch_sl_cells);% first is wave, second is fetch!!
% time_vis = toc;
WA_vislib_ = cell2mat(WA_vislib);

cellsize = 62.5;
nrays = 36:360; % number of rays around the circle. Recommend nrays >= 36
delta = 0.01:0.01:0.5; % distance increment along each ray, in cells. Recommend delta <= 0.5
nrays_test = 360/5;
delta_test = 0.5;
time_rays = zeros(size(nrays));
time_delta = zeros(size(delta));

% loop over rays
for n = 1:length(nrays)
    tic
    [warea,~,~,~,~,~,~,~] = GetFetchArea(fetch_sl_cells,F_lake,nrays(n),delta_test,cellsize);
    time_rays(n) = toc;
    warea_rays_save{n} = warea;
    diff_rays{n} = WA_vislib_ - cell2mat(warea_rays_save{n});
    sum_diff_rays(n) = sum(diff_rays{n});
    sum_abs_diff_rays(n) = sum(abs(diff_rays{n}));
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
plot(delta,sum_diff_delta,'k','LineWidth',2)
xlabel('delta')
ylabel('sum of fetch difference between vis and ray')
set(gca,'FontSize',12)

subplot(2,1,2)
plot(nrays,sum_diff_rays,'k','LineWidth',2)
xlabel('number of rays')
ylabel('sum of fetch difference between vis and ray')
set(gca,'FontSize',12)

% plot absolute difference
figure()
subplot(2,1,1)
plot(delta,sum_abs_diff_delta,'k','LineWidth',2)
xlabel('delta')
ylabel('sum of abs (difference between vis and ray)')
set(gca,'FontSize',12)

subplot(2,1,2)
plot(nrays,sum_abs_diff_rays,'k','LineWidth',2)
xlabel('number of rays')
ylabel('sum of abs (difference between vis and ray)')
set(gca,'FontSize',12)
