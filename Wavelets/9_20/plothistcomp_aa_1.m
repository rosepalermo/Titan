%% file to look at histograms from all shorelines and compute some
% statistics

clear
close all

% load the file with all of the small-scale variance data
load alleq14_n6

% some general controls

fs = 10;    % font size
titleboost = 1.2 %
edgespace = .05 ;
edgemax = 1.5 ;
radscale = 10^-4;
ylimits = [0 0.2]
xlimits = [0 edgemax];

% this is the shoreline that you get. in short:
% 1 rednoise_v1_';
datatitle(1,:) = string('Red Noise Initial Conditions');
datatitle(2,:) = 'Early Waves';
datatitle(3,:) = 'Wave Erosion';
datatitle(4,:) = 'Early Uniform Erosion';
datatitle(5,:) = 'Uniform Erosion';
datatitle(6,:) = 'River Incisision';
datatitle(7,:) = 'Lake Powell';
datatitle(8,:) = 'Uist';
datatitle(9,:) = 'Ligea Mare';
datatitle(10,:) = 'Lake Sebago';
datatitle(11,:) = 'Western Uist';
datatitle(12,:) = 'Eastern Uist';

%% main loop





for jj = 1:12
    hh(jj) = subplot(3,4, jj)
end

for ii = 1:12
    
    % load a specific set of values
    
    subplot(hh(ii))
    
    if ii < 11
        vars = cell2mat(eq14save(ii))/radscale;
        
    elseif ii == 11
        vars =  cell2mat(eq14save(8))/radscale;
        vars = vars(1:floor(size(vars)/3));
    elseif ii == 12
        vars =  cell2mat(eq14save(8))/radscale;
        vars = vars(floor(size(vars)/3):end);
    end
    % else if ii = 12
    
    % make a histogram with these controls
    meanvar(ii) = mean(vars);
    skewvar(ii) = skewness(vars);
    
    edges = (0:edgespace:edgemax);
    
    h = histogram(vars,edges,'Normalization','probability');
    
    histvalues(ii,:) = h.Values;
    histbins(ii,:) = h.BinEdges;
    
    meantext= ['Mean = ' num2str(meanvar(ii))];
    skewtext= ['Skewness = ' num2str(skewvar(ii))];
    thetext = {meantext, skewtext};
    % text(edgemax * 0.65, ylimits * 0.8, thetext)
    
    ylim(ylimits)
    xlim(xlimits)
    
    h.FaceColor = 'k';
    h.EdgeColor = 'w';
    xlabel('Wavelet Variance (10^4) (^c)^2')
    ylabel('Density')
    set(gca,'fontweight','bold')
    set(gca,'fontsize',fs)
    title(datatitle(ii,:),'fontsize',fs*titleboost)
    hold all
    
end

figure

for jj = 1:3
    hh2(jj) = subplot(2,3,jj);
end

subplot(hh2(1))

% waves figure
for j = 1:3
    plot(histbins(j,1:(end-1)),histvalues(j,:), 'linewidth',2)
    hold on
    
    legend('Initial','Time 1','Time 2')
    
    set(gca,'fontweight','bold')
    set(gca,'fontsize',fs)
    
    ylim(ylimits)
    xlim(xlimits)
    xlabel('Wavelet Variance (10^4) (^c)^2')
    ylabel('Density')
    title('Wave Erosion','fontsize',fs*titleboost)
end

subplot(hh2(2))
for j = [1 4 5]
    plot(histbins(j,1:(end-1)),histvalues(j,:), 'linewidth',2)
    hold on
    
    legend('Initial','Time 1','Time 2')
    
    set(gca,'fontweight','bold')
    set(gca,'fontsize',fs)
    ylim(ylimits)
    xlim(xlimits)
    xlabel('Wavelet Variance (10^4) (^c)^2')
    ylabel('Density')
    title('Uniform Erosion','fontsize',fs*titleboost)
end


subplot(hh2(3))
for j = [1 6]
    plot(histbins(j,1:(end-1)),histvalues(j,:), 'linewidth',2)
    hold on
    
    %ylim([0 ylimits])
    legend('Initial','Time 1','Time 2')
    
    set(gca,'fontweight','bold')
    set(gca,'fontsize',fs)
    ylim(ylimits)
    xlim(xlimits)
    xlabel('Wavelet Variance (10^4) (^c)^2')
    ylabel('Density')
    title('River Incision','fontsize',fs*titleboost)
end


