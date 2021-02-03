function [period, eq14] = dowave_greece(y,dt,ord,xx,yy,savename,save_on,fetch,pmin,pmax)
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

% % close the loop
% figure(); scatter3(xx,yy,[(0:length(xx)-1)*dt],100,[(0:length(xx)-1)*dt],'.');view(2);grid off;axis equal tight; colormap parula; h = colorbar; xlabel('meters'); ylabel('meters'); ylabel(h,'# points'); set(gca,'FontSize',16)
% % 
%     if save_on
%         fig = '.eps'; fig_suf ='_pts'; figname = strcat(savename,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end

n = length(y);
t = (0:length(y)-1)*dt;  % construct time array
xlimits = [min(t),max(t)];  % plotting range
pad = 0;      % pad the time series with zeroes (recommended - not doing this for Titan because closed loop)
dj = 0.25;    % this will do 4 sub-octaves per octave
s0 = 2*dt; % figure()
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%
% % contour(t,log2(period),log2(power),log2(levels));  %*** or use 'contourf'
% % imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% colormap parula
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
% h = figure();
% 
% %--- Plot time series
% subplot('position',[0.1 0.75 0.65 0.2]);
% scatter3(t,y,[(0:length(xx)-1)*dt],100,[(0:length(xx)-1)*dt],'.');view(2);
% hold on;
% plot(t,y,'k')
% set(gca,'XLim',xlimits(:))
% % xlabel('Alongshore location')
% ylabel('Azimuth')
% title('data')
% set(gca,'FontSize',16)
% hold off
% 
% %--- Contour plot wavelet power spectrum
% subplot('position',[0.1 0.37 0.65 0.28]);
% % levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% 
% % contour(t,log2(period),log2(power),log2(levels));  %*** or use 'contourf'
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% colormap hsv
% 
% xlabel('alongshore position (meters)')
% ylabel('Period (units of Alongshore location)')
% title('Wavelet Power Spectrum')
% load('climits_wps.mat')
% set(gca,'Clim',climits)
% set(gca,'XLim',xlimits(:))
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks)
% set(gca,'FontSize',16)
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% hold on
% contour(t,log2(period),sig95,[-99,1],'r','linewidth',1);
% hold on
% % cone-of-influence, anything "below" is dubious
% plot(t,log2(coi),'k')
% hold off


%--- Plot global wavelet spectrum
% subplot('position',[0.77 0.37 0.2 0.28]);
% semilogx(fft_theor,log2(period),'k')
% hold on
% semilogx(global_signif,log2(period),'r')
% semilogx(global_ws,log2(period))
% legend('fft estimate','significance','spatial average')
% hold off
% xlabel('Power (amplitude^2)')
% title('Global Wavelet Spectrum')
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel','')
% set(gca,'XLim',[0,1.25*max(global_ws)])
% h.Position = [569  -217   950   699];
%     if save_on
%         fig = '.eps'; fig_suf ='_wps_all'; figname = strcat(savename,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
% 
%     figure()
% semilogx(fft_theor,log2(period),'k')
% hold on
% semilogx(global_signif,log2(period),'r')
% semilogx(global_ws,log2(period))
% legend('fft estimate','significance','spatial average')
% hold off
% xlabel('Power (amplitude^2)')
% title('Global Wavelet Spectrum')
% set(gca,'YLim',log2([4,256]), ...
% 	'YDir','reverse')
% set(gca,'XLim',[0,1.25*max(global_ws)])
%         if save_on
%         fig = '.eps'; fig_suf ='global'; figname = strcat(savename,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
%Model lakes
% pmin = 2^2;
% pmax = 2^3.8;


% h = figure();
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% colormap parula
% xlabel('Alongshore location')
% ylabel('Period (units of Alongshore location)')
% title('Wavelet Power Spectrum')
% load('climits_wps.mat')
% set(gca,'Clim',climits)
% set(gca,'XLim',xlimits(:))
%     set(gca,'YLim',log2([min(period),max(period)]), ...
%     	'YDir','reverse', ...
%     	'YTick',log2(Yticks(:)), ...
%     	'YTickLabel',Yticks)
% set(gca,'YLim',log2([min(period),max(period)]), ...
%     'YDir','reverse', ...
%     'YTick',log2(Yticks(:)), ...
%     'YTickLabel',Yticks)
% set(gca,'FontSize',16)
% load('position_wps.mat')
% h.Position = position_wps;
%     if save_on
%         fig = '.eps'; fig_suf ='_wps'; figname = strcat(savename,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end


% don't like doing bespoke data editing but here we go anyway
%big kluge
% this was for scotland - not using but I renamed the values
absolutepower = size(power);
realabsolutepower = absolutepower(2);


goodi = 1:realabsolutepower;


% this computation is i specific and comes after plotting
pband1 = period >= pmin & period <= pmax;

powernorm_sub = power(pband1,goodi);
eq14 = dj.*dt./0.776./length(y)*sum(powernorm_sub./scale(pband1)');




%% shoreline w/ eq4
%     figure()
%     scatter3(xx,yy,eq14',[],eq14','filled')
%     view(2)
%     grid off
%     axis tight
%     axis equal
% %     xlim([0 200])
% %     ylim([20 180])
% %         set(gca,'Ydir','reverse')
%     h = colorbar;
%     ylabel(h, 'azimuthal variance')
%     set(gca,'Clim',[0 3e-4])
%     xlabel('meters');ylabel('meters');
% %     load('clim_eq14.mat')
% %     set(gca,'Clim',clim_eq14)
%     set(gca,'FontSize',16)
% %     set(gca,'ColorScale','log')
% %     set(gca,'YTickLabel',[])
% %     set(gca,'XTickLabel',[])
% %     set(gca,'fontsize',18)
%     grid off
%     if save_on
%         fig = '.eps'; fig_suf ='_eq14_sl'; 
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
%% Plot fetch on shoreline
%         figure()
%     scatter3(xx,yy,fetch,[],fetch,'filled')
%     view(2)
%     grid off
%     axis tight
%     axis equal
% %     xlim([0 200])
% %     ylim([20 180])
% %         set(gca,'Ydir','reverse')
%     h = colorbar;
%     ylabel(h, 'fetch area')
% %     set(gca,'Clim',[1e9 1e10])
%     xlabel('meters');ylabel('meters');
% %     load('clim_eq14.mat')
% %     set(gca,'Clim',clim_eq14)
%     set(gca,'FontSize',16)
% %     set(gca,'ColorScale','log')
% %     set(gca,'YTickLabel',[])
% %     set(gca,'XTickLabel',[])
% %     set(gca,'fontsize',18)
%     if save_on
%         fig = '.eps'; fig_suf ='_ww_sl'; 
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
    
%% Plot fetch vs roughness
% figure()
% % scatter(log(fetch),eq14')
% X = [log(fetch),eq14'];
% X(X==inf) = NaN;
% X(X==-inf) = NaN;
% hist3(X,'CdataMode','auto'); view(2)
% xlabel('log fetch')
% ylabel('azimuthal variance (radians^2)')
% set(gca,'FontSize',16)
%     if save_on
%         fig = '.eps'; fig_suf ='fvr'; 
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
    
    %% Plot fetch vs roughness
    %% linear
% figure()
% % scatter(log(fetch),eq14')
% [B,~,idx] = histcounts(fetch);
% % plot(B);
% % xlabel('bin')
% % ylabel('N')
%     if save_on
%         fig = '.eps'; fig_suf ='hist_lin'; 
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
% figure()
% idx = idx+1; % this is because it starts at 0
% meaneq14 = accumarray(idx(:),eq14,[],@mean);
% meaneq14(meaneq14==0)=NaN;
% meanfetch = accumarray(idx(:),fetch,[],@mean);
% meanfetch(meanfetch==0)=NaN;
% medianeq14 = accumarray(idx(:),eq14,[],@median);
% medianeq14(medianeq14==0)=NaN;
% medianfetch = accumarray(idx(:),fetch,[],@median);
% medianfetch(medianfetch==0)=NaN;
% stdeq14 = accumarray(idx(:),eq14,[],@std);
% stdeq14(stdeq14==0)=NaN;
% B = [1 B]';
% SEM = stdeq14./sqrt(B);                         % Standard Error Of The Mean
% CI95 = SEM .* tinv(0.975, B-1);              % 95% Confidence Intervals
% stdfetch = accumarray(idx(:),fetch,[],@std);
% stdfetch(stdfetch==0)=NaN;
% 
% plot(meanfetch,meaneq14,'k','LineWidth',2)
% hold on
% plot(fetch,eq14,'.','Color',[0.8 0.8 0.8])
% % scatter(medianfetch,medianeq14,'k*')
% errorbar(meanfetch,meaneq14,CI95,'k')
% % legend('mean','median')
% xlabel('weighted fetch area')
% ylabel('azimuthal variance (radians^2)')
% set(gca,'FontSize',16)
%     if save_on
%         fig = '.eps'; fig_suf ='fvr_mean'; 
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
    %% log
%     figure()
% scatter(log(fetch),eq14')
[B,~,idx] = loghistcounts(fetch);
% plot(B);
% xlabel('bin')
% ylabel('N')
%     if save_on
%         fig = '.eps'; fig_suf ='hist_log'; 
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
% figure()
idx = idx+1; % this is because it starts at 0
meaneq14 = accumarray(idx(:),eq14,[],@mean);
meaneq14(meaneq14==0)=NaN;
meanfetch = accumarray(idx(:),fetch,[],@mean);
meanfetch(meanfetch==0)=NaN;
medianeq14 = accumarray(idx(:),eq14,[],@median);
medianeq14(medianeq14==0)=NaN;
medianfetch = accumarray(idx(:),fetch,[],@median);
medianfetch(medianfetch==0)=NaN;
stdeq14 = accumarray(idx(:),eq14,[],@std);
stdeq14(stdeq14==0)=NaN;
B = [1 B]';
SEM = stdeq14./sqrt(B);                         % Standard Error Of The Mean
CI95 = SEM .* tinv(0.975, B-1);              % 95% Confidence Intervals
stdfetch = accumarray(idx(:),fetch,[],@std);
stdfetch(stdfetch==0)=NaN;

p_2 = semilogx(meanfetch(B>1),meaneq14(B>1),'k','LineWidth',2);
hold on
p2_2 = plot(fetch,eq14,'.','Color',[0.8 0.8 0.8]);
% plot(fetch,eq14,'.','Color','g')
% scatter(medianfetch,medianeq14,'k*')
p3_2 = errorbar(meanfetch(B>1),meaneq14(B>1),CI95(B>1),'k');
% legend('mean','median')
ylim([0 5e-4])
xlabel('weighted fetch area')
ylabel('azimuthal variance (radians^2)')
% legend('mean','data','95% CI')
set(gca,'FontSize',16)
    if save_on
        fig = '.eps'; fig_suf ='fvr_mean_log'; 
        range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
        figname = strcat(savename,range,fig_suf,'.jpg');
        saveas(gcf,figname)
    end

