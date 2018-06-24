% plot the shorelines


% for scotland
% convert long and lat to x and y (in km) -- NEED TO DO FOR EACH
%the below conversion factors are calculated for 57 degrees latitude, a
%central point along the cape, by a National Geospatial-Intelligence Agency
%calculator:
%http://msi.nga.mil/MSISiteContent/StaticFiles/Calculators/degree.html
lat = 111.360;
lon = 60.772;
SLdata(:,1) = (SLdata(:,1)-min(SLdata(:,1))) * lon;%changed from lon/lat*lon to just lon
SLdata(:,2) = (SLdata(:,2)-min(SLdata(:,2))) * lat;




% load the data
% load Ligeia Mare
filename = 'lg_4wavelets.xls';
M = xlsread(filename);
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];
xlgm=M(:,4);
ylgm=M(:,5);

% figure(); plot(xlgm/1000,ylgm/1000)
% xlabel('km')
% ylabel('km')
% set(gca,'FontSize',14)
% axis equal tight

%  Red noise t1
load('wave_rednoise.mat')
x0 = ordered_sl_save{1,1}{1,1}(:,1);
y0 = ordered_sl_save{1,1}{1,1}(:,2);
M = [x0 y0];
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 1) M(duplicate_ind, 2)];
M(duplicate_ind,:)=[];
xt1 = M(:,1);
yt1 = M(:,2);

%  Red noise eroded by waves
x0 = ordered_sl_save{7,1}{1,1}(:,1);
y0 = ordered_sl_save{7,1}{1,1}(:,2);
M = [x0 y0];
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 1) M(duplicate_ind, 2)];
M(duplicate_ind,:)=[];
xwaves = M(:,1);
ywaves = M(:,2);

% red noise eroded by rivers at t = 2
load('xycontours.mat')
xrivers = x_1m_t2(1:4:end)';
yrivers = y_1m_t2(1:4:end)';


% Red noise eroded by uniform model
load('uniform_rednoise.mat')
shoreline = addidshoreline_cardonly(lake_save{30,1},~lake_save{30,1});
[sl_cell,cells2trash] = order_cw_lastpoint(lake_save{30,1},shoreline);
xuniform = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
yuniform = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
%%
fig1 = figure();
plot(xt1,yt1,'k')
hold on
plot(xrivers,yrivers,'Color',[0.7 0.7 0.7])
plot(xuniform,yuniform,'b')
plot(xwaves,ywaves,'r')
xlabel('km')
ylabel('km')
set(gca,'FontSize',14)
legend('initial conditions','rivers','uniform','waves','Location','southeast')
fig1.Position = ([604,101,775,639])

figure()
subplot(2,2,1)
plot(xt1,yt1,'k')
hold on
plot(xrivers,yrivers,'Color',[0.7 0.7 0.7])
plot(xuniform,yuniform,'b')
plot(xwaves,ywaves,'r')
set(gca,'XLim',([500 800])); set(gca,'YLim',([650 800]));
xlabel('km')
ylabel('km')
set(gca,'FontSize',14)
set(gca,'XLim',([500 800])); set(gca,'YLim',([650 800]));
legend('initial conditions','rivers','uniform','waves','Location','southeast')


load('x_1m_t2_wavelet.mat')
subplot(2,2,2)
scatter3(xx/1e3,yy/1e3,norm_rness_unsmoothed,[],norm_rness_unsmoothed,'.')
colormap jet
colorbar
view(2)
axis equal tight
xlabel('km')
ylabel('km')
title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
set(gca,'XLim',([500 800])); set(gca,'YLim',([650 800]));
set(gca,'Clim',[0 1])
set(gca,'FontSize',14)

load('uniform_rednoise_wavelet.mat')
subplot(2,2,3)
scatter3(xx/1e3,yy/1e3,norm_rness_unsmoothed,[],norm_rness_unsmoothed,'.')
colormap jet
colorbar
view(2)
axis equal tight
xlabel('km')
ylabel('km')
title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
set(gca,'XLim',([500 800])); set(gca,'YLim',([650 800]));
set(gca,'Clim',[0 1])
set(gca,'FontSize',14)

load('wave_rednoise_wavelet.mat')
subplot(2,2,3)
scatter3(xx/1e3,yy/1e3,norm_rness_unsmoothed,[],norm_rness_unsmoothed,'.')
colormap jet
colorbar
view(2)
axis equal tight
xlabel('km')
ylabel('km')
title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
set(gca,'XLim',([500 800])); set(gca,'YLim',([650 800]));
set(gca,'Clim',[0 1])
set(gca,'FontSize',14)