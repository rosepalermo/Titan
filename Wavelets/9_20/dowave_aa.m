function [period, eq14] = dowave(y,dt,ord,xx,yy,savename,save_on,fetch,i)
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
imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
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
% figure()
% plot(power,period)

% energy = power*dt;
% figure()
% plot(energy,period)
if i == 10
    plotidx = cellstr(num2str(t(1:length(t)/20:end-20)'));
else
    plotidx = cellstr(num2str(ceil(t(1:length(t)/20:end-20)')));
end

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

% %% PLOTTING the Wavelet
% figure
% 
% %--- Plot time series
% subplot('position',[0.1 0.75 0.65 0.2])
% plot(t,y)
% set(gca,'XLim',xlim(:))
% xlabel('alongshore distance (m)')
% ylabel('y')
% % title('data')
% set(gca,'FontSize',12)
% hold off
% %text(t(1:length(t)/20:end-20),y(1:length(t)/20:end-20),plotidx,'FontSize',14);
% 
% %--- Contour plot wavelet power spectrum
% subplot('position',[0.1 0.37 0.65 0.28])
% % levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% 
% % contour(t,log2(period),log2(power),log2(levels));  %*** or use 'contourf'
% % imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% colormap jet
% set(gca,'FontSize',12)
% 
% 
% xlabel('alongshore distance (m)')
% ylabel('Period')
% % title('Wavelet Power Spectrum')
% set(gca,'XLim',xlim(:))
% % % set(gca,'YLim',log2([min(period),max(period)]), ...
% % 	'YDir','reverse', ...
% % 	'YTick',log2(Yticks(:)), ...
% % 	'YTickLabel',Yticks)
% % % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% % hold on
% % contour(t,log2(period),sig95,[-99,1],'r','linewidth',1);
% % set(gca,'FontSize',12)
% % hold on
% % cone-of-influence, anything "below" is dubious
% % plot(t,log2(coi),'k')
% % hold off
% 
% 
% %--- Plot global wavelet spectrum
% subplot('position',[0.77 0.37 0.2 0.28])
semilogx(fft_theor,log2(period),'k')
hold on
semilogx(global_signif,log2(period),'r')
semilogx(global_ws,log2(period))
global_Save = global_ws;
hold off
% xlabel('Power (amplitude^2)')
% % title('Global Wavelet Spectrum')
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel','')
set(gca,'XLim',[0,1.25*max(global_ws)])
% set(gca,'FontSize',12)
% fig1 = gcf;
% fig1.Position = ([0,0,750,400]);
% if save_on
%     fig = '.eps'; pws ='pws'; figname = strcat(savename,pws);
%     print(figname,'-depsc')
%     fig = '.fig'; pws ='pws'; figname = strcat(savename,pws);
%     saveas(gcf,figname)
%     
% end



if i ==9
    % Ligeia Mare
    pmin1 = 2^10;
    pmax1 = 2^14;
    pmin2 = 2^10;
    pmax2 = 2^15;
end


if i == 1|i == 2|i == 3|i == 4|i == 5|i == 6
    %Model lakes
    pmin1 = 2^2;
    pmax1 = 2^4;
    pmin2 = 2^3;
    pmax2 = 2^4;
end

if i == 7| i == 8 | i == 11 | i ==12
    %Lake Powell or SCOTLAND
    pmin1 = 2^8;
    pmax1 = 2^10.5;
    pmin2 = 2^9;
    pmax2 = 2^14;
end

if i == 10
    % Sebago
    pmin1 = 2^7;
    pmax1 = 2^10;
end

size(power)
fracloc = 0.3;
% don't like doing bespoke data editing but here we go anyway
%big kluge
% this was for scotland - not using but I renamed the values
absolutepower = size(power);
realabsolutepower = absolutepower(2);


if i == 11
    %eastern Scotland isle
    
    goodi = 1:floor(realabsolutepower*fracloc);

elseif i == 12  
   goodi = ceil(realabsolutepower*fracloc):realabsolutepower;
 
else
    goodi = 1:realabsolutepower;
   
end

% this computation is i specific and comes after plotting 
pband1 = period >= pmin1 & period <= pmax1;

powernorm_sub = power(pband1,goodi);
eq14 = dj.*dt./0.776./length(y)*sum(powernorm_sub./scale(pband1)');


%% new figure
% figure() 
% plot(t/100,eq14)
% xlabel('distance along coast (km)')
% title('sum of eq14')
% if save_on
%     fig = '.eps'; EQ14vt ='EQ14vt'; figname = strcat(savename,EQ14vt);
%     print(figname,'-depsc')
%     fig = '.fig'; EQ14vt ='EQ14vt'; figname = strcat(savename,EQ14vt);
%     saveas(gcf,figname)
% end
% text(t(1:length(t)/20:end-20),eq14(1:length(t)/20:end-20),plotidx,'FontSize',14);

% figure()
% h = histogram(eq14,10,'Normalization','probability');
% 
% ppp = size(eq14)
% % h = findobj(gca,'Type','patch');
% h.FaceColor = 'k';
% % h.FaceColor = [0.6 0.6 0.6];
% h.EdgeColor = 'w';
% title('sum of eq14')
% if save_on
%     fig = '.eps'; EQ14hist ='EQ14hist'; figname = strcat(savename,EQ14hist);
%     print(figname,'-depsc')
%     fig = '.fig'; EQ14hist ='EQ14hist'; figname = strcat(savename,EQ14hist);
%     saveas(gcf,figname)
% end


%% shoreline w/ eq4
if i < 100%i == 1|i == 2|i == 3|i == 4|i == 5 |i == 6 
    figure()
    offset = 10000;
    scatter3(xx/1e3+offset,yy/1e3+offset,eq14,20,eq14,'filled')
    %     scatter3(xx/1e3+0.75,yy/1e3+1.65,eq14,30,eq14,'filled')

    view(2)
    axis equal tight
    
  %  set(gca,'XLim',([2.100 2.400])); set(gca,'YLim',([1.800 2.000])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
 %   set(gca,'XLim',([2.0500 2.4500])); set(gca,'YLim',([1.7500 2.0500])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])

 

  set(gca,'Clim',[0 0.00005])
    set(gca,'FontSize',14)
    set(gca,'xtick',[],'ytick',[])
    set(gca,'xticklabel',[],'yticklabel',[])
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
    [max(xx) max(yy)]
    % Full figure
    
%     figure()
%     scatter3(xx,yy,eq14,[],eq14,'.')
%    % title('sum of eq14')
%     view(2)
%     colorbar
%     set(gca,'Clim',[0 0.00005])
%     axis equal tight
%     %text(xx(1:length(t)/20:end-20),yy(1:length(t)/20:end-20),plotidx,'FontSize',14);
else
    figure()
    scatter3(xx,yy,eq14,[],eq14,'.')
    title('sum of eq14')
    view(2)
    colorbar
    set(gca,'Clim',[0 0.00005])
    axis equal tight
    %text(xx(1:length(t)/20:end-20),yy(1:length(t)/20:end-20),plotidx,'FontSize',14);
end
if save_on
    fig = '.eps'; EQ14_ ='EQ14'; figname = strcat(savename,EQ14_);
    print(figname,'-depsc')
    fig = '.fig'; EQ14_ ='EQ14'; figname = strcat(savename,EQ14_);
    saveas(gcf,figname)
end


% %% wavelet power zoom
% figure()
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% imagesc(t,log2(period(pband1)),log2(power(pband1,:)));  %*** uncomment for 'image' plot
% colormap jet
% set(gca,'FontSize',12)
% xlabel('alongshore distance (m)')
% ylabel('Period')
% title('Wavelet Power Spectrum')
% set(gca,'XLim',xlim(:))
% set(gca,'CLim',[-12 8])
% colorbar
% if save_on
%     fig = '.eps'; wpszoom_ ='wps_zoom'; figname = strcat(savename,wpszoom_);
%     print(figname,'-depsc')
%     fig = '.fig'; wpszoom_ ='wps_zoom'; figname = strcat(savename,wpszoom_);
%     saveas(gcf,figname)
% end


%% test before AGU 18
% rness = sum(powernorm_sub);
% figure() 
% plot(t/100,rness)
% xlabel('distance along coast (km)')
% title('sum of wavelet power')
% if save_on
%     fig = '.eps'; rnessvt ='rnessvt'; figname = strcat(savename,rnessvt);
%     print(figname,'-depsc')
% end
% if i == 12
%     figure()
%     scatter3(t,y,rness,[],rness,'.')
%     title('sum of wavelet power')
%     view(2)
%     set(gca,'Clim',[0 30])
%     colorbar
%     axis equal tight
% elseif i == 2|i == 3|i == 4|i == 5
%     figure()
%     scatter3(xx/1e3,yy/1e3,rness,[],rness,'.')
%     view(2)
%     axis equal tight
%     set(gca,'XLim',([0.5 0.8])); set(gca,'YLim',([0.5 0.8])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
%     set(gca,'Clim',[0 30])
%     set(gca,'FontSize',14)
%     set(gca,'xtick',[],'ytick',[])
%     set(gca,'xticklabel',[],'yticklabel',[])
%     if save_on
%         fig = '.eps'; rn3z ='rnesszoom'; figname = strcat(savename,rn3z);
%         print(figname,'-depsc')
%     end
%     figure()
%     scatter3(xx,yy,rness,[],rness,'.')
%     title('sum of wavelet power')
%     view(2)
%     colorbar
%     set(gca,'Clim',[0 30])
%     axis equal tight
% else
%     figure()
%     scatter3(xx,yy,rness,[],rness,'.')
%     title('sum of wavelet power')
% %     set(gca,'Clim',[0 30])
%     view(2)
%     colorbar
%     axis equal tight
% end
% if save_on
%     fig = '.eps'; RNESS_ ='RNESS'; figname = strcat(savename,RNESS_);
%     print(figname,'-depsc')
% end
% 
% figure()
% h = histogram(rness,10,'Normalization','probability')
% % h = findobj(gca,'Type','patch');
% h.FaceColor = 'k';
% % h.FaceColor = [0.6 0.6 0.6];
% h.EdgeColor = 'w';
% title('sum of wavelet power')
% if save_on
%     fig = '.eps'; rness_hist ='RNESShist'; figname = strcat(savename,rness_hist);
%     print(figname,'-depsc')
% end
% 
% rness_energy = sum(energy(pband1,:));
% 
% figure() 
% plot(t/100,rness_energy)
% title('sum of wavelet power*dt')
% xlabel('distance along coast (km)')
% if save_on
%     fig = '.eps'; rness_energyvt ='rness_energyvt'; figname = strcat(savename,rness_energyvt);
%     print(figname,'-depsc')
% end
% 
% figure()
% h = histogram(rness_energy,10,'Normalization','probability')
% % h = findobj(gca,'Type','patch');
% h.FaceColor = 'k';
% % h.FaceColor = [0.6 0.6 0.6];
% h.EdgeColor = 'w';
% title('sum of wavelet power*dt')
% if save_on
%     fig = '.eps'; rness_energy_hist ='rness_energyhist'; figname = strcat(savename,rness_energy_hist);
%     print(figname,'-depsc')
% end
% 
% if i == 12
%     figure()
%     scatter3(t,y,rness_energy,[],rness_energy,'.')
%     title('sum of wavelet power*dt')
%     view(2)
%     colorbar
%     set(gca,'Clim',[0 30])
%     axis equal tight
%     elseif i == 2|i == 3|i == 4|i == 5
%     figure()
%     scatter3(xx/1e3,yy/1e3,rness_energy,[],rness_energy,'.')
%     view(2)
%     axis equal tight
%     set(gca,'XLim',([0.5 0.8])); set(gca,'YLim',([0.5 0.8])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
%     set(gca,'CLim',[0 30])
%     set(gca,'FontSize',14)
%     set(gca,'xtick',[],'ytick',[])
%     set(gca,'xticklabel',[],'yticklabel',[])
%     if save_on
%         fig = '.eps'; rn3z ='Energyzoom'; figname = strcat(savename,rn3z);
%         print(figname,'-depsc')
%     end
%     figure()
%     scatter3(xx,yy,rness_energy,[],rness_energy,'.')
%     title('sum of wavelet power*dt')
%     view(2)
%     colorbar
%     set(gca,'CLim',[0 30])
%     axis equal tight
% else
%     figure()
%     scatter3(xx,yy,rness_energy,[],rness_energy,'.')
%     title('sum of wavelet power*dt')
%     view(2)
%     colorbar
% %     set(gca,'CLim',[0 50])
%     axis equal tight
% end
% if save_on
%     fig = '.eps'; RENERGY ='RENERGY'; figname = strcat(savename,RENERGY);
%     print(figname,'-depsc')
% end
%% test before generals
% norm_rness_unsmoothed = rness./sum(global_ws);
% norm_rness_unsmoothed = rness./sum(global_ws(pband1));
% norm_rness_unsmoothed = rness./mean(global_ws);
% norm_rness_unsmoothed = (rness-min(rness))./max(rness);
% norm_rness_unsmoothed = rness./mean(rness);
% norm_rness_unsmoothed = (rness)./max(rness);
% norm_rness_unsmoothed = (rness)./max(max(power));
% norm_rness_unsmoothed = (rness)./mean(power(pband1,:));
% figure
% plot(t/1000,rness,'g')
% xlabel('distance along coast (km)')
% ylabel('roughness')
% title('roughness')
% set(gca,'FontSize',14)
% ylim([0 90])
% if save_on
%     fig = '.eps'; rn ='rn'; figname = strcat(savename,rn);
%     print(figname,'-depsc')
% end

% figure()
% % h = histogram(rness(1:length(y)/2),10,'Normalization','probability')
% % h = histogram(rness(length(y)/2:end),10,'Normalization','probability')
% h = histogram(rness,10,'Normalization','probability')
% % h = findobj(gca,'Type','patch');
% h.FaceColor = 'k';
% % h.FaceColor = [0.6 0.6 0.6];
% h.EdgeColor = 'w';
% xlabel('roughness')
% ylabel('frequency')
% set(gca,'FontSize',20)
% ylim([0 0.5])
% title('roughness')
% % set(gca,'xlim',[0 90])
% % yticks([0 0.25 0.5])
% % xticks([0 20 40 60 80])
% % if save_on
% %     fig = '.eps'; his ='his'; figname = strcat(savename,his);
% %     print(figname,'-depsc')
% % end


% rms_ness = rms(rness);
% var_ness = var(rness);
% mean_ness = mean(rness);
% median_ness = median(rness);
% skewness_ness = skewness(rness);

% if ~isempty(fetch)
%     figure()
%     scatter(fetch,norm_rness_unsmoothed,'.','k')
%     xlabel('Wave Weighting')
%     ylabel('Roughness')
% end

% % figure
% % % if ~isempty(fetch)
% % %     subplot(1,2,1)
% % % end
% % scatter3(xx/1e3,yy/1e3,norm_rness_unsmoothed,[],norm_rness_unsmoothed,'.')
% % colormap jet
% % colorbar
% % view(2)
% % axis equal tight
% % % xlabel('km')
% % % ylabel('km')
% % % title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
% % % set(gca,'Clim',[0 mean(rness)+2*std(rness)])
% % % set(gca,'Clim',[0 90])
% % % set(gca,'Clim',[0 9])
% % set(gca,'FontSize',14)
% % set(gca,'xtick',[],'ytick',[])
% % set(gca,'xticklabel',[],'yticklabel',[])
% % 
% % % if ~isempty(fetch)
% % %     subplot(1,2,2)
% % %     scatter3(xx/1e3,yy/1e3,fetch,[],fetch,'.')
% % %     colormap jet
% % %     colorbar
% % %     view(2)
% % %     axis equal tight
% % %     xlabel('km')
% % %     ylabel('km')
% % %     % title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
% % %     % set(gca,'Clim',[0 0.75])
% % %     set(gca,'FontSize',14)
% % % end
% % if save_on
% %     fig = '.eps'; rn3 ='rn3'; figname = strcat(savename,rn3);
% %     print(figname,'-depsc')
% % end
% % 
% % % norm_rness zoomed
% % figure
% % scatter3(xx/1e3,yy/1e3,norm_rness_unsmoothed,[],norm_rness_unsmoothed,'.')
% % colormap jet
% % colorbar
% % view(2)
% % axis equal tight
% % % xlabel('km')
% % % ylabel('km')
% % % title(['normalized sum of the power spectrum in the ' num2str(pmin1/1e3) '-' num2str(pmax1/1e3) ' km band'])
% % set(gca,'XLim',([0.5 0.8])); set(gca,'YLim',([0.5 0.8])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
% % set(gca,'FontSize',14)
% % % set(gca,'Clim',[0 90])
% % set(gca,'Clim',[0 9])
% % set(gca,'xtick',[],'ytick',[])
% % set(gca,'xticklabel',[],'yticklabel',[])
% % if save_on
% %     fig = '.eps'; rn3z ='rn3z'; figname = strcat(savename,rn3z);
% %     print(figname,'-depsc')
% % end
% % 
% % if save_on
% %     save(savename,'xx','yy','norm_rness_unsmoothed','pmin1','pmax1','period','power')
% % end


% there's a signal there, but it's noisy, so let's smooth it 
% Lsm = 5*pmax1; % smoothing window length in meters
% rnesssm = movmean([rness rness rness],round(Lsm/dt)); rnesssm = rnesssm(n+1:2*n); % note that we took advantage of the periodicity of the spectrum in t
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
% % % 
% % ff = figure()
% % ff.Position = [793 267 560 538];
% % %plot wavelet power spectrum
% % ax2 = subplot(2,1,2)
% % % figure()
% % Yticks = (floor(log2(min(period)))):ceil((log2(max(period))));
% % if i ==1
% %     imagesc(t./1e6,log2(period),log2(power));  %*** uncomment for 'image' plot
% % else
% %     imagesc(t,log2(period),log2(power));  %*** uncomment for 'image' plot
% % end
% % colormap gray
% % colorbar
% % if i == 1
% %     xlabel('alongshore distance (1000s of km)')
% % else
% %     xlabel('alongshore distance (m)')
% % end
% % ylabel('wavelength')
% % if i == 1
% %     set(gca,'XLim',xlim(:)./1e6)
% % else
% %     set(gca,'XLim',xlim(:))
% % end
% % set(gca,'cLim',[-25 10])
% % set(gca,'YLim',log2([min(period),max(period)]), ...
% % 	'YDir','reverse', ...
% % 	'YTick',(Yticks(:)), ...
% % 	'YTickLabel',Yticks)
% % set(gca,'Fontsize', 16)
% % % linkaxes([ax1,ax2],'x')
% % % 
% % % %plot just time series
% % fff = figure()
% % fff.Position = [793 267 560 538];ax1 = subplot(2,1,1)
% % z = zeros(size(t));
% % yplot = y';
% % if i == 1
% % surface([t./1e6;t./1e6],[yplot;yplot],[z;z],[t;t],...
% %         'facecol','no',...
% %         'edgecol','interp',...
% %         'linew',2);
% % else
% % surface([t;t],[yplot;yplot],[z;z],[t;t],...
% %         'facecol','no',...
% %         'edgecol','interp',...
% %         'linew',2);
% % end
% % % set(gca,'XLim',xlim(:)./1e6)
% % if i == 1
% %     xlabel('alongshore distance (1000s of km)')
% % else
% %     xlabel('alongshore distance (m)')
% % end
% % ylabel('azimuth (radians)')
% % colormap parula
% % set(gca,'Fontsize', 16)
% % view(2)
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

