
save_on = false;

addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Cellular Damage Model')
for ii = 8
    i = ii;

clearvars -except period global_Save i save_on
load('shorelines_4AGU18.mat')

% 2 = REDNOISE
savename{1} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/rednoise_v1_'; 
% 3 = WAVE t1v1_30
savename{2} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/wave_t1v1_30_'; 
% 3 = WAVE t1v1_50
savename{3} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/wave_t1v1_50_';
% 5 = UNIFORM t1v1_30
savename{4} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/uniform_t1v1_30_'; 
% 5 = UNIFORM t1v1_50
savename{5} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/uniform_t1v1_50_'; 
% 4 = RIVERS t3v1
savename{6} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/river_t3v1_'; 
% 6 = Lake Powell
savename{7} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/LakePowell_'; 
% 7 = Scotland
savename{8} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/scotland_'; 
% 9 = LGM
savename{9} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/lgm'; 
% 9 = LGM
savename{10} = '/Volumes/GoogleDrive/My Drive/RoseAndrewShare/AGU 2018/figures from rose/sebago'; 

fetch = [];
if save_on
savename = savename{i};
end
if i < 7 % one of the models
%     xx = A(5).cord{1,1}(:,2); yy = A(5).cord{1,1}(:,1);
%     A(5).cord{1,1}(:,1) = yy;
%     A(5).cord{1,1}(:,2) = xx;
x0 = SL_matrix(i).cord{1,1}(:,1);
y0 = SL_matrix(i).cord{1,1}(:,2);

% fetch = WaveArea_save{1, i}{1, 1};
M = [x0 y0];
% indices to unique values in column 3
[~, ind] = unique(M, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 1) M(duplicate_ind, 2)];
M(duplicate_ind,:)=[];
% fetch(duplicate_ind) = [];


x0=M(:,1);
y0=M(:,2);
end

if i == 7
% %Lake Powell
filename = 'LakePowell_gp.csv';
M = csvread(filename);
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);
lat = 110.978;
lon = 89.012;
x0 = (x0-min(x0)) * lon;%changed from lon/lat*lon to just lon
y0 = (y0-min(y0)) * lat;

dist = sqrt((x0(2:end) - x0(1:end-1)).^2 +  (y0(2:end) - y0(1:end-1)).^2);
dist_off = find(dist>1);
x0_fixed = [x0(1:dist_off(1)-1);x0(dist_off(2)+1:dist_off(3)-1);x0(dist_off(1)+1:dist_off(2)-1);x0(dist_off(12)+1:dist_off(13)-1);x0(dist_off(13)+1:end);x0(dist_off(11)+1:dist_off(12)-1);x0(dist_off(10)+1:dist_off(11)-1);x0(dist_off(6)+1:dist_off(7)-1);x0(dist_off(9)+1:dist_off(10)-1);x0(dist_off(3)+1:dist_off(4)-1);x0(dist_off(7)+1:dist_off(8)-1);x0(dist_off(5)+1:dist_off(6)-1);x0(dist_off(4)+1:dist_off(5)-1)];
y0_fixed = [y0(1:dist_off(1)-1);y0(dist_off(2)+1:dist_off(3)-1);y0(dist_off(1)+1:dist_off(2)-1);y0(dist_off(12)+1:dist_off(13)-1);y0(dist_off(13)+1:end);y0(dist_off(11)+1:dist_off(12)-1);y0(dist_off(10)+1:dist_off(11)-1);y0(dist_off(6)+1:dist_off(7)-1);y0(dist_off(9)+1:dist_off(10)-1);y0(dist_off(3)+1:dist_off(4)-1);y0(dist_off(7)+1:dist_off(8)-1);y0(dist_off(5)+1:dist_off(6)-1);y0(dist_off(4)+1:dist_off(5)-1)];
x0_fixed(1245) = []; y0_fixed(1245) = [];
x0 = x0_fixed*1000; y0 = y0_fixed*1000; % convert to meters

fetch = [];
end
% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/LakePowell_gp_wavelet'; 
% % 


%% Scotland
if i == 8
filename = 'Scotland_gp.csv';
M = csvread(filename);
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);
lat = 111360;
lon = 60772;
x0 = (x0-min(x0)) * lon;%changed from lon/lat*lon to just lon
y0 = (y0-min(y0)) * lat;
dist = sqrt((x0(2:end) - x0(1:end-1)).^2 +  (y0(2:end) - y0(1:end-1)).^2);
dist_off = find(dist>1000);
x0_fixed = [flipud(x0(924:2345));x0(1:923);x0(2346:end)];
y0_fixed = [flipud(y0(924:2345));y0(1:923);y0(2346:end)];
% plot(x0_fixed,y0_fixed)
x0 = x0_fixed; y0 = y0_fixed;
fetch = [];
end

%% load Ligeia Mare
% clear
if i == 9
filename = 'lg_4wavelets.xls';
M = xlsread(filename);

% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/lg_4wavelets_wavelet_updated'; 
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);
end

%% Sebago
if i == 10
M = xlsread('sebago_2_xy_UTM.xls');
%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);
x0 = [x0(3319+25:3459);x0(2223:2787);x0(3460:4155);x0(5442:7123);x0(87:1664);x0(1941:2222);x0(4156:5044);x0(8329:8703);x0(7124:8328);x0(5381:5440);x0(5045:5379)];
y0 = [y0(3319+25:3459);y0(2223:2787);y0(3460:4155);y0(5442:7123);y0(87:1664);y0(1941:2222);y0(4156:5044);y0(8329:8703);y0(7124:8328);y0(5381:5440);y0(5045:5379)];
x_island = [x0(1665:1940);x0(2788:3318)];
y_island = [y0(1665:1940);y0(2788:3318)];

lat = 111103;
lon = 80837;
x0 = (x0-min(x0)) * lon;%changed from lon/lat*lon to just lon
y0 = (y0-min(y0)) * lat;
end
%%  Red noise eroded by waves
% clear
% load('wave_rednoise_wrong.mat')
% x0 = ordered_sl_save{7,1}{1,1}(:,1);
% y0 = ordered_sl_save{7,1}{1,1}(:,2);
% M = [x0 y0];
% %get rid of duplicate points (when it goes exactly around a pixel)
% % indices to unique values in column 3
% [~, ind] = unique(M, 'rows');
% % duplicate indices
% duplicate_ind = setdiff(1:size(M, 1), ind);
% % duplicate values
% duplicate_val = [M(duplicate_ind, 1) M(duplicate_ind, 2)];
% M(duplicate_ind,:)=[];
% x0 = M(:,1);
% y0 = M(:,2);
% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/Figures in prep/wave_rednoise_wavelet'; 
% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/wave_rednoise_wavelet_updated'; 

%%  Red noise t1
% clear
% load('wave_rednoise.mat')
% x0 = ordered_sl_save{1,1}{1,1}(:,1);
% y0 = ordered_sl_save{1,1}{1,1}(:,2);
% M = [x0 y0];
% %get rid of duplicate points (when it goes exactly around a pixel)
% % indices to unique values in column 3
% [~, ind] = unique(M, 'rows');
% % duplicate indices
% duplicate_ind = setdiff(1:size(M, 1), ind);
% % duplicate values
% duplicate_val = [M(duplicate_ind, 1) M(duplicate_ind, 2)];
% M(duplicate_ind,:)=[];
% x0 = M(:,1);
% y0 = M(:,2);
% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/rednoiset1_wavelet_updated';


%% red noise eroded by rivers at t = 2
% clear
% load('xycontours.mat')
% x0 = x_1m_t2(1:4:end)';
% y0 = y_1m_t2(1:4:end)';
% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/x_1m_t2_wavelet_updated'; 


%% Red noise eroded by uniform model
% clear
% load('uniform_rednoise.mat')
% shoreline = addidshoreline_cardonly(lake_save{30,1},~lake_save{30,1});
% [sl_cell,cells2trash] = order_cw_lastpoint(lake_save{30,1},shoreline);
% x0 = X(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
% y0 = Y(sub2ind(size(X),sl_cell{1,1}(:,1),sl_cell{1,1}(:,2)));
% savename = '/Users/rosepalermo/Documents/Research/Titan/Notes/Generals/uniform_rednoise_wavelet_updated';
%%
% CONSIDER NOT INTERPOLATING HERE. WE DON'T NEED EVENLY SPACED POINTS TO
% COMPUTE AZIMUTHS, AND IT CAN CREATE STEPS IN THE AZIMUTH VS. DISTANCE
% SERIES IF INTERPOLATING CREATES MULTIPLE SUCCESSIVE POINTS WITH THE SAME
% AZIMUTH.

% % % interpolate
% dist=sqrt((x0(2:end)-x0(1:end-1)).^2+(y0(2:end)-y0(1:end-1)).^2); % distance increment
% dist_cum = [0;cumsum(dist)];
% dist_total = sum(dist);
% dist_interp = linspace(0,dist_total,length(x0)/8);
% x_int = interp1(dist_cum,x0,dist_interp);
% y_int = interp1(dist_cum,y0,dist_interp);
% figure()
% scatter(x0,y0,'r')
% hold on
% scatter(x_int,y_int,'k')
% 
% x0 = x_int';
% y0 = y_int';


% % rearrange ligeia for debugging
% x0 = [x0(23186:end);x0(1:23185)];
% y0 = [y0(23186:end);y0(1:23185)];



% add first AND SECOND points to the end for meander
x = [x0;x0(1:2)];
y = [y0;y0(1:2)];



% transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% streamwise distance
[theta, dtheta, deltad] = meander_titan(x,y);
x = x(1:end-2);
y = y(1:end-2);


%% Taylor's synthetic coast tests
% % TRY SMOOTHING AZIMUTH
% Lsm = 7; % smoothing window length in # points
% thetasm = movmean([theta; theta; theta],Lsm); thetasm = thetasm(length(x)+1:2*length(x)); % note that we took advantage of the periodicity of the spectrum in t

% % HERE IS A SYNTHETIC SIGNAL WITH THE SAME LENGTH AS theta
% 
% t=deltad*(0:length(theta)-1);
% t = t(:);
% T = t(end)+deltad;
% hwin = 0.5*(1-cos(2*pi*t/T));
% synth = 4*sin(2*pi*t/(T/4)) + hwin.*( 1*sin(2*pi*t/(T/64)) );
% 
% n=2; i = 2;
% dowave(synth,deltad,n,x,y,savename,save_on,[],i);
% 
% 
% th = linspace(0,2*pi,60);
% r = 2 + rand(size(th))-0.5;
% x = r.*cos(th);
% y = r.*sin(th);
% N=100; %1000 data points
% th = linspace(0,2*pi,N)'; %theta from 0 to 2pi
% amp = 2*ones(N,1);
% for i=1:5 %higher numbers here make more higher frequency fluctuations
%     a = rand()-0.5; %random number from -0.5 to 0.5
%     b = rand()-0.5;
%     amp = amp + a/i*cos(i*th) + b/i*sin(i*th); %change circle by a fourier function over it
% end
% x=(amp.*cos(th))'; y = (amp.*sin(th))';
% % add first AND SECOND points to the end for meander
% x = [x0;x0(1:2)];
% y = [y0;y0(1:2)];
% 
% % transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% % streamwise distance
% [theta, dtheta, deltad] = meander_titan(x,y);
% x = x(1:end-2);
% y = y(1:end-2);
% 
% i=2;
% n=2;
% dowave(theta,deltad,n,x,y,savename,save_on,[],i);
% 
% % NOW SHUFFLE ORDER AS A TEST:
% ihalf = round(length(t)/2);
% x = [x(ihalf+1:end);x(1:ihalf)];
% y = [y(ihalf+1:end);y(1:ihalf)];
% synth = [synth(ihalf+1:end); synth(1:ihalf)];
% 
% n=2;
% dowave(synth,deltad,n,x,y,savename); % NO DIFFERENCE IN ROUGHNESS MAP!
%%
% do the wavelet transforms, comparing the spectra to those of an 
% autoregressive process of order n (a.k.a. an AR(n) process)
n=2;
% dowave_duplicate(theta,deltad,n,x,y,savename);

% theta = [theta(length(theta)*2/3:end);theta;theta(1:length(theta)*1/3)];
% x = [x(length(x)*2/3:end);x;x(1:length(x)*1/3)];
% y = [y(length(y)*2/3:end);y;y(1:length(y)*1/3)];
[period{i},global_Save{i}] = dowave(theta,deltad,n,x,y,savename,save_on,fetch,i);

% % plot where the first point is
% figure
% % plot(x,y,'k','LineWidth',2)
% hold on
% scatter(x,y,'.','k')
% scatter(x(1),y(1),'*','r')
% scatter(x(30),y(30),'*','b')
% scatter(x(ceil(length(x)/3)),y(ceil(length(y)/3)),'*','m')
% xlabel('X')
% ylabel('Y')
% axis tight
% axis equal

end
