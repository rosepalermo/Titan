function [period, global_Save] = dowave_sat(y,dt,ord,xx,yy,savename,save_on,fetch,i)
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

n = length(y);
t = (0:length(y)-1)*dt;  % construct time array
xlim = [min(t),max(t)];  % plotting range
pad = 0;      % pad the time series with zeroes (recommended - not doing this for Titan because closed loop)
dj = 0.25;    % this will do 4 sub-octaves per octave
s0 = 2*dt;    % this says start at a scale of 2*dt (analogous to the Nyquist wavelength of a time series)
% npow = 10;   % number of powers of two to do.
% j1 = npow/dj;    % this says do 7 powers-of-two with dj sub-octaves each
j1=-1; % This will use the default number of powers of two
mother = 'Morlet'; % Mother wavelet

% Wavelet transform:
[wave,period,scale,coi] = wavelet(y,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;  % compute wavelet power spectrum
figure()
plot(power,period)

energy = power*dt;
% figure()
% plot(energy,period)

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
plot(t,y)
set(gca,'XLim',xlim(:))
xlabel('alongshore distance (m)')
ylabel('y')
% title('data')
set(gca,'FontSize',12)
hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.37 0.65 0.28])
% levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));

% contour(t,log2(period),log2(power),log2(levels));  %*** or use 'contourf'
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
colormap gray
set(gca,'FontSize',12)


xlabel('alongshore distance (m)')
ylabel('Period')
% title('Wavelet Power Spectrum')
set(gca,'XLim',xlim(:))
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks)
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% hold on
% contour(t,log2(period),sig95,[-99,1],'r','linewidth',1);
% set(gca,'FontSize',12)
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
global_Save = global_ws;
hold off
xlabel('Power (amplitude^2)')
% title('Global Wavelet Spectrum')
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

% powernorm = power./repmat(global_ws(:),[1 n]);
powernorm = power; % should be normalized by global otherwise as the power increases with period, higher periods will dominate roughness
% figure
% imagesc(t,log2(period),log2(powernorm));
% colorbar
% set(gca,'clim',[-5 5])
% title('powernorm')

% note the differences in power along the shoreline at wavelengths
% (periods) of about 2^11 to 2^14 m! Let's use the normalized power summed
% over that range of periods as a measure of relative roughness within that
% wavelength range.

if i ==1
    % Ligeia Mare
    pmin1 = 2^10;
    pmax1 = 2^15;
    pmin2 = 2^10;
    pmax2 = 2^15;
end


if i == 2|i == 3|i == 4|i == 5
    %Model lakes
    pmin1 = 2^2;
    pmax1 = 2^4;
    pmin2 = 2^3;
    pmax2 = 2^4;
end

if i == 6| i == 7
    %Lake Powell or scotland
    pmin1 = 2^9;
    pmax1 = 2^14;
    pmin2 = 2^9;
    pmax2 = 2^14;
end

if i == 12
    pmin1 = 0;
    pmax1 = 2^5;
end

if i == 13
    pmin1 = 2^1.8045;
    pmax1 = 2^3.3045;
    pmin2 = 2^4.55;
    pmax2 = 2^5.805;
end

pband1 = period >= pmin1 & period <= pmax1;
pband2 = period >= pmin2 & period <= pmax2;

powernorm_sub = powernorm(pband1,:);
powernorm_sub2 = powernorm(pband2,:);

eq14_1 = sum(dj.*dt./0.776./length(y)*sum(powernorm_sub./scale(pband1)'))
eq14_2 = sum(dj.*dt./0.776./length(y)*sum(powernorm_sub2./scale(pband2)'))


sumeq14 = eq14_1 + eq14_2
