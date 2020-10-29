function [period, eq14] = dowave_greece(y,dt,ord,xx,yy,savename,save_on,fetch,i)
%%
% Computes the continuous wavelet transform of a function y with point
% spacing dt, compares it with the spectrum for a random process with
% the same autocorrelation structure, and plots the results. The optional
% input argument ord allows the user to specify the autocorrelation of the
% random ivprocess used to estimate the background spectrum: for an AR(1),
% set ord=1; for an AR(2), set ord=2, etc.e
%
%
% Modified from Torrence & Compo's routine.

% remove the mean and normalize to unit variance. Could also do first-order
% detrending if desired.
% variance = std(y)^2;
% y = (y - mean(y))/sqrt(variance);

% close the loop
figure(); scatter3(xx,yy,[(0:length(xx)-1)*dt],100,[(0:length(xx)-1)*dt],'.');view(2);axis equal tight; colormap jet; colorbar

n = length(y);
t = (0:length(y)-1)*dt;  % construct time array
xlim = [min(t),max(t)];  % plotting range
pad = 0;      % pad the time series with zeroes (recommended - not doing this for Titan because closed loop)
dj = 0.25;    % this will do 4 sub-octaves per octave
s0 = 2*dt; % figure()
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%
% % contour(t,log2(period),log2(power),log2(levels));  %*** or use 'contourf'
% % imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% colormap jet
% set(gca,'FontSize',12)
% set(gca,'Clim',[-20 15])
% this says start at a scale of 2*dt (analogous to the Nyquist wavelength of a time series)
% npow = 10;   % number of powers of two to do.
% j1 = npow/dj;    % this says do 7 powers-of-two with dj sub-octaves each
j1=-1; % This will use the default number of powers of two
mother = 'Morlet'; % Mother wavelet

% Wavelet transform:
[wave,period,scale,coi] = wavelet(y,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;  % compute wavelet power spectrum

% Background spectrum and significance levels:

% Determine the order of the AR(n) if it was not specified
if nargin < 3
    maxlag=20; % number of lags up to which we'll calculate partial autocorr.
    
    % determine the order and noise amplitude of the autocorrelated process
    pacf = parcorfn(y,maxlag);
    siglvl = 2/sqrt(n);
    j=1;
    while abs(pacf(j))>siglvl && j<=maxlag
        j=j+1;
    end
    ord=j-2;
end


[a,e] = aryule(y,ord); % Estimate autocorrelation coefficients and
% noise amplitude using Yule-Walker method. More
% specifically, by solving the Yule-Walker equations,
% we can determine the parameters for an all-pole
% filter that when excited with white noise will
% produce an AR signal whose statistics match those
% of the real data. a contains the coefficients for
% the AR(ord) process, and e is the variance of the
% Gaussian noise.


% Find the spectrum of a random signal with the same autocorrelation and
% variance characteristics of the data:

% Method 1: a Monte Carlo approach
%
% iter=100; % number of Monte Carlo iterations to do
% fft_theor=zeros(size(wave)); % this will hold the background spectrum
%
% h = waitbar(0,'Performing Monte Carlo iterations...');
%
% for i=1:iter
%     % generate a random, autocorrelated time series with the same length,
%     % order, and noise amplitude as the data
%     AR=filter(1,a,sqrt(e)*randn(n,1));
%
%     % find the wavelet power spectrum
%     AR = wavelet(AR,dt,pad,dj,s0,j1,mother);
%     AR = (abs(AR)).^2;
%
%     % add it to the array
%     fft_theor = fft_theor + AR;
%
%     waitbar(i/iter)
% end
%
% close(h)
%
% fft_theor = fft_theor/iter; % find the mean background spectrum
% global_fft_theor = (sum(fft_theor')/n); % average over position to get
%                                         % global theoretical spectrum

% Method 2: a filter approach
%
% Calculate the frequency response of the filter used to construct the
% AR(ord) process. The frequency response of the filter is an estimate of
% the fourier transform of the AR(ord) process, and agrees very well with the
% Monte Carlo result.
fft_theor = freqz(sqrt(e),a,1./period,1/dt);
fft_theor = fft_theor.*conj(fft_theor); % the power spectrum

% [signif,fft_theor] = wave_signif(std(y)^2,dt,scale,0,lag1,-1,-1,mother);
signif = wave_signif_ARn(std(y)^2,dt,scale,0,fft_theor,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:
% global_ws = variance*(sum(power')/n);   % time-average over all times
global_ws = (sum(power')/n);   % time-average over all times
% while retaining unit variance
dof = n - scale;  % the -scale corrects for padding at edges
% global_signif = wave_signif(std(y)^2,dt,scale,1,lag1,-1,dof,mother);
global_signif = wave_signif_ARn(std(y)^2,dt,scale,1,fft_theor,-1,dof,mother);

eq14 = dj.*dt./0.776./length(y)*sum(sum(power./scale'));

% PLOTTING
figure

%--- Plot time series
subplot('position',[0.1 0.75 0.65 0.2])
scatter3(t,y,[(0:length(xx)-1)*dt],100,[(0:length(xx)-1)*dt],'.');view(2)
hold on;
plot(t,y,'k')
set(gca,'XLim',xlim(:))
xlabel('Alongshore location')
ylabel('Azimuth')
title('data')
hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.37 0.65 0.28])
% levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));

% contour(t,log2(period),log2(power),log2(levels));  %*** or use 'contourf'
imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
colormap jet

xlabel('t')
ylabel('Period (units of Alongshore location)')
title('Wavelet Power Spectrum')
load('climits_wps.mat')
set(gca,'Clim',climits)
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% hold on
% contour(t,log2(period),sig95,[-99,1],'r','linewidth',1);
% hold on
% % cone-of-influence, anything "below" is dubious
% plot(t,log2(coi),'k')
% hold off


%--- Plot global wavelet spectrum
subplot('position',[0.77 0.37 0.2 0.28])
semilogx(fft_theor,log2(period),'k')
hold on
semilogx(global_signif,log2(period),'r')
semilogx(global_ws,log2(period))
hold off
xlabel('Power (amplitude^2)')
title('Global Wavelet Spectrum')
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel','')
set(gca,'XLim',[0,1.25*max(global_ws)])


%Model lakes
pmin1 = 2^2;
pmax1 = 2^3.8;


h = figure();
imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
colormap jet
xlabel('Alongshore location')
ylabel('Period (units of Alongshore location)')
title('Wavelet Power Spectrum')
load('climits_wps.mat')
set(gca,'Clim',climits)
set(gca,'XLim',xlim(:))
    set(gca,'YLim',log2([pmin1,pmax1]), ...
    	'YDir','reverse', ...
    	'YTick',log2(Yticks(:)), ...
    	'YTickLabel',Yticks)
% set(gca,'YLim',log2([min(period),max(period)]), ...
%     'YDir','reverse', ...
%     'YTick',log2(Yticks(:)), ...
%     'YTickLabel',Yticks)
set(gca,'FontSize',14)
load('position_wps.mat')
h.Position = position_wps;

fracloc = 0.3;
% don't like doing bespoke data editing but here we go anyway
%big kluge
% this was for scotland - not using but I renamed the values
absolutepower = size(power);
realabsolutepower = absolutepower(2);


goodi = 1:realabsolutepower;


% this computation is i specific and comes after plotting
pband1 = period >= pmin1 & period <= pmax1;

powernorm_sub = power(pband1,goodi);
eq14 = dj.*dt./0.776./length(y)*sum(powernorm_sub./scale(pband1)');




%% shoreline w/ eq4
    figure()
    scatter3(xx,yy,eq14',[],eq14','filled')
    %     scatter3(xx/1e3+0.75,yy/1e3+1.65,eq14,30,eq14,'filled')
    
    view(2)
    axis equal tight
    colorbar
    %  set(gca,'XLim',([2.100 2.400])); set(gca,'YLim',([1.800 2.000])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
    %   set(gca,'XLim',([2.0500 2.4500])); set(gca,'YLim',([1.7500 2.0500])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
    
    
    load('clim_eq14.mat')
    set(gca,'Clim',clim_eq14)
    set(gca,'FontSize',14)
%     set(gca,'xtick',[],'ytick',[])
%     set(gca,'xticklabel',[],'yticklabel',[])
    %text(xx(1:length(t)/20:end-20)/1e3,yy(1:length(t)/20:end-20)/1e3,plotidx,'FontSize',14);
    if save_on
        fig = '.eps'; rn3z ='eq14zoom'; figname = strcat(savename,rn3z);
        print(figname,'-depsc')
        fig = '.fig'; rn3z ='eq14zoom'; figname = strcat(savename,rn3z);
        saveas(gcf,figname)
    end
    %  colorbar
    set(gca,'fontsize',18)
    %     set(gca,'fontweight','bold')
    grid off
    % Full figure
    
    %     figure()
    %     scatter3(xx,yy,eq14,[],eq14,'.')
    %    % title('sum of eq14')
    %     view(2)
    %     colorbar
    %     set(gca,'Clim',[0 0.00005])
    %     axis equal tight
    %     %text(xx(1:length(t)/20:end-20),yy(1:length(t)/20:end-20),plotidx,'FontSize',14);

if save_on
    fig = '.eps'; EQ14_ ='EQ14'; figname = strcat(savename,EQ14_);
    print(figname,'-depsc')
    fig = '.fig'; EQ14_ ='EQ14'; figname = strcat(savename,EQ14_);
    saveas(gcf,figname)
end

% Plot fetch vs roughness
figure()
scatter(log(fetch),eq14')
xlabel('log fetch')
ylabel('wavelet variance')
set(gca,'FontSize',14)

