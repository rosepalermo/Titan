function dowave(y,dt,ord,xx,yy,savename,save_on)

% Computes the continuous wavelet transform of a function y with point
% spacing dt, compares it with the spectrum for a random process with 
% the same autocorrelation structure, and plots the results. The optional
% input argument ord allows the user to specify the autocorrelation of the
% random process used to estimate the background spectrum: for an AR(1),
% set ord=1; for an AR(2), set ord=2, etc.
%
%
% Modified from Torrence & Compo's routine.

% remove the mean and normalize to unit variance. Could also do first-order
% detrending if desired.
variance = std(y)^2;
y = (y - mean(y))/sqrt(variance);

n = length(y);
t = (0:length(y)-1)*dt;  % construct time array
xlim = [min(t),max(t)];  % plotting range
pad = 0;      % pad the time series with zeroes (recommended)
dj = 0.25;    % this will do 4 sub-octaves per octave
s0 = 2*dt;    % this says start at a scale of 2*dt (analogous to the Nyquist wavelength of a time series)
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


% PLOTTING
figure

%--- Plot time series
subplot('position',[0.1 0.75 0.65 0.2])
plot(t,y)
set(gca,'XLim',xlim(:))
xlabel('t')
ylabel('y')
title('data')
set(gca,'FontSize',12)
hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.37 0.65 0.28])
% levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));

% contour(t,log2(period),log2(power),log2(levels));  %*** or use 'contourf'
imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
colormap gray
set(gca,'FontSize',12)


xlabel('t')
ylabel('Period (units of t)')
title('Wavelet Power Spectrum')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(t,log2(period),sig95,[-99,1],'r','linewidth',1);
set(gca,'FontSize',12)
% hold on
% cone-of-influence, anything "below" is dubious
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
set(gca,'FontSize',12)
fig1 = gcf;
fig1.Position = ([0,0,750,400]);
if save_on
    fig = '.png'; pws ='pws'; figname = strcat(savename,pws,fig);
    saveas(gca,figname)
end

% TAYLOR 9 MAY 2018

% normalize the wavelet spectrum by the global wavelet spectrum and display
% it

powernorm = power./repmat(global_ws(:),[1 n]);
% powernorm = power; % should be normalized by global otherwise as the power increases with period, higher periods will dominate roughness
% figure
% imagesc(t,log2(period),log2(powernorm));
% colorbar
% set(gca,'clim',[-5 5])
% title('powernorm')

% note the differences in power along the shoreline at wavelengths
% (periods) of about 2^11 to 2^14 m! Let's use the normalized power summed
% over that range of periods as a measure of relative roughness within that
% wavelength range.

% Ligeia Mare
pmin1 = 2^11;
pmax1 = 2^17; 
pmin2 = 2^10;
pmax2 = 2^15;


% Model lakes 
pmin1 = 8;
pmax1 = 256; 
pmin2 = 8;
pmax2 = 16;

% % test new
% pmin1 = 8;
% pmax1 = 16;



%Lake Powell
% pmin1 = 2^9;
% pmax1 = 2^14; 
% pmin2 = 2^9;
% pmax2 = 2^14;

pband1 = period >= pmin1 & period <= pmax1;
powernorm_sub = powernorm(pband1,:);

rness = sum(powernorm_sub);
norm_rness_unsmoothed = rness./mean(global_ws);
% norm_rness_unsmoothed = (rness-min(rness))./max(rness);
% norm_rness_unsmoothed = (rness)./max(max(power));
% norm_rness_unsmoothed = (rness)./mean(power(pband1,:));
figure
plot(t/1000,norm_rness_unsmoothed,'-b')
xlabel('distance along coast (km)')
ylabel('roughness')
set(gca,'FontSize',14)
if save_on
    fig = '.png'; rn ='rn'; figname = strcat(savename,rn,fig);
    saveas(gca,figname)
end

rms(norm_rness_unsmoothed)
var(norm_rness_unsmoothed)



figure
scatter3(xx/1e3,yy/1e3,norm_rness_unsmoothed,[],norm_rness_unsmoothed,'.')
colormap jet
colorbar
view(2)
axis equal tight
xlabel('km')
ylabel('km')
title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
%set(gca,'Clim',[0 1])
set(gca,'FontSize',14)
if save_on
    fig = '.png'; rn3 ='rn3'; figname = strcat(savename,rn3,fig);
    saveas(gca,figname)
end

% norm_rness zoomed
figure
scatter3(xx/1e3,yy/1e3,norm_rness_unsmoothed,[],norm_rness_unsmoothed,'.')
colormap jet
colorbar
view(2)
axis equal tight
xlabel('km')
ylabel('km')
title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
set(gca,'XLim',([0.5 0.8])); set(gca,'YLim',([0.5 0.8])); %set(gca,'Clim',[0 0.75])
set(gca,'FontSize',14)
if save_on
    fig = '.png'; rn3z ='rn3z'; figname = strcat(savename,rn3z,fig);
    saveas(gca,figname)
end

if save_on
    save(savename,'xx','yy','norm_rness_unsmoothed','pmin1','pmax1','period','power')
end


% there's a signal there, but it's noisy, so let's smooth it 
Lsm = 5*pmax1; % smoothing window length in meters
rnesssm = movmean([rness rness rness],round(Lsm/dt)); rnesssm = rnesssm(n+1:2*n); % note that we took advantage of the periodicity of the spectrum in t
% hold on
% plot(t/1000,rnesssm,'-r') % better!

% % normalize it to range from 0 to 1
% rnesssm = rnesssm-min(rnesssm);
% rnesssm = rnesssm/max(rnesssm);

% % plot it along the coast
% figure
% scatter3(xx/1e3,yy/1e3,rnesssm,[],rnesssm,'.')
% colormap jet
% colorbar
% view(2)
% axis equal tight
% xlabel('km')
% ylabel('km')
% title(['roughness in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])

% % Ratio of roughness in smaller band to larger band
% pband2 = period >= pmin2 & period <= pmax2;
% powerratio_sub = powernorm(pband2,:);
% powerratio_ness = sum(powerratio_sub);
% rationess = rness./powerratio_ness;
% rationess = rationess./max(rationess);
% figure
% plot(t/1000,rationess,'-b')
% xlabel('distance along coast (km)')
% ylabel('roughness ratio of smaller pmin to pmax')

% % plot it along the coast
% figure
% scatter3(xx/1e3,yy/1e3,rationess,[],rationess,'.')
% colormap jet
% colorbar
% view(2)
% axis equal tight
% xlabel('km')
% ylabel('km')
% title(['roughness ratio in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) 'compared to the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])


% m=zeros(length(power(1,:)),1);
% figure
% for i=1:length(power(1,:))
%     p=polyfit(log2(period'),log2(power(:,i)),1);
%     m(i)=p(1);
%     v(i) = var(log2(power(:,i)));
%     stdev(i) = std(log2(power(:,i)));
%     
%     %     figure
%     %     plot(log2(period'),log2(power(:,i)),'k')
%     scatter((period'),(power(:,i)),'k','MarkerFaceAlpha',0.2)
%     hold on
%     %     yslope=polyval(p,log2(period'));
%     %     plot(log2(period'),yslope,'r')
%     
% end
% 
% figure
% for i=1:length(power(1,:))
%     p=polyfit(log2(period'),log2(power(:,i)),1);
%     m(i)=p(1);
%     v(i) = var(log2(power(:,i)));
%     stdev(i) = std(log2(power(:,i)));
%     
%     %     figure
%         scatter(log2(period'),log2(power(:,i)),'k','MarkerFaceAlpha',0.2)
% %     plot((period'),(power(:,i)),'k')
%     hold on
%     %     yslope=polyval(p,log2(period'));
%     %     plot(log2(period'),yslope,'r')
%     
% end


%pause
% 
% %plot the slope of the power spectrum at each point
% % xx=xx(1:end-2);
% % yy=yy(1:end-2);
% figure
% scatter3(xx,yy,m,[],m,'.')
% colormap(parula)
% colorbar
% set(gca,'Clim',[mean(m)-2*std(m),mean(m)+2*std(m)])
% view(2)
% title('m')
% axis equal
% axis equal
% 
% figure()
% plot(t,m)
% title('m')
% 
% figure
% scatter3(xx,yy,v,[],v,'.')
% colormap(parula)
% colorbar
% set(gca,'Clim',[mean(v)-2*std(v),mean(v)+2*std(v)])
% view(2)
% title('var')
% axis equal
% axis equal
% 
% figure()
% plot(t,v)
% title('var')
% 
% figure
% scatter3(xx,yy,stdev,[],stdev,'.')
% colormap(parula)
% colorbar
% set(gca,'Clim',[mean(stdev)-2*std(stdev),mean(stdev)+2*std(stdev)])
% view(2)
% title('stdev')
% axis equal
% axis equal
% 
% figure()
% plot(t,stdev)
% title('stdev')
% 
% % 
% % %moving average--need the financial toolbox
% % floor(0.1*length(power(1,:)))
% % %mper=[m(length(m)/2:end);m;m(1:length(m)/2)]; %make think it's periodic (add either end to the other side--result is 2x length of the data set)
% % mper=m;
% % mtest=tsmovavg(mper,'t',floor(0.01*length(mper)),1); %window=1/10 of lake points
% % %movav=mtest((length(m)/2)+2:(length(m))+(length(m)/2)+1);
% % movav=mtest;
% % figure
% % scatter3(xx,yy,movav,[],movav,'.')
% % colormap(flipud(jet))
% % colorbar
% % %set(gca,'Clim',[0.8,2.2])
% % view(2)
% % axis equal
% % axis equal
% 
% normm=m./max(m);
% figure
% scatter3(xx,yy,normm,[],normm,'.')
% colormap(parula)
% colorbar
% title('normalized spectral slope')
% %set(gca,'Clim',[0.8,2.2])
% view(2)
% axis equal
% axis equal
% 
% 
% %plot just time series
% figure
% ax1 = subplot(2,1,1)
% plot(t,y,'k','LineWidth',1.5)
% set(gca,'XLim',xlim(:))
% xlabel('alongshore distance (m)')
% ylabel('azimuth (radians)')
% set(gca,'Fontsize', 16)
% 
% 
% %plot wavelet power spectrum
% ax2 = subplot(2,1,2)
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% colormap parula
% hold on
% contour(t,log2(period),sig95,[-99,1],'r','linewidth',1);
% xlabel('alongshore distance (m)')
% ylabel('Period (units of t)')
% set(gca,'XLim',xlim(:))
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks)
% set(gca,'Fontsize', 16)
% linkaxes([ax1,ax2],'x')
% 
% %plot just time series
% figure
% ax1 = subplot(2,1,1)
% scatter3(t,y,t,[],t,'.')
% set(gca,'XLim',xlim(:))
% xlabel('alongshore distance (m)')
% ylabel('azimuth (radians)')
% set(gca,'Fontsize', 16)
% view(2)
% title('time series with alongshore coloring')
% 
% 
% %plot wavelet power spectrum
% ax2 = subplot(2,1,2)
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% colormap parula
% hold on
% contour(t,log2(period),sig95,[-99,1],'r','linewidth',1);
% xlabel('alongshore distance (m)')
% ylabel('Period (units of t)')
% set(gca,'XLim',xlim(:))
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks)
% set(gca,'Fontsize', 16)
% linkaxes([ax1,ax2],'x')
% 
% figure()
% title('lake with alongshore coloring')
% scatter3(xx,yy,t,[],t,'.')
% view(2)
% axis equal
% 
% 
% 
% 
% % 
% % imshow(A)
% % hold on
% % scatter3(xx./10,-yy./10,movav,[],movav,'.')
% % scatter(xx./100,yy./100,'*')
% % colormap(flipud(jet))
% % colorbar
% 
% % name='spectralslope';
% % name=repmat(name,length(xmid),1);
% % S = geoshape(xx,yy,name,movav);
% % shapewrite(S,shapefilename)
% 
% 
% 
% save=[xx yy normm];
% csvwrite(savename,save)
% 

