%% file to look at histograms from all shorelines and compute some
% statistics

% clear
% close all
% 
% % load the file with all of the small-scale variance data
load histwaveuniform.mat

% some general controls

fs = 28;    % font size
titleboost = 1.2 %
edgespace = .005 ;
edgemax = .15 ;
radscale = 10^-4;
% ylimits = [0 0.2]
xlimits = [0 edgemax];


%% main loop


for ii = 1:2
    
    % load a specific set of values
    figure
    %vars = cell2mat(eq14save(ii))/radscale;
    
    
    if ii < 11
        vars = cell2mat(eq14save(ii))/radscale;
        
    elseif ii == 11
        vars =  cell2mat(eq14save(8))/radscale;
        vars = vars(1:floor(size(vars)/3));
    elseif ii == 12
        vars =  cell2mat(eq14save(8))/radscale;
        vars = vars(floor(size(vars)/3):end);
    end
    
    
    
    
    % make a histogram with these controls
    
    meanvar(ii) = mean(vars);
    skewvar(ii) = skewness(vars);
    
    edges = (0:edgespace:edgemax);
    
    
    
    h = histogram(vars,edges,'Normalization','probability');
    
    
    meantext= ['Mean = ' num2str(meanvar(ii))];
    skewtext= ['Skewness = ' num2str(skewvar(ii))];
    thetext = {meantext, skewtext};
   % text(edgemax * 0.65, max(h.Values) * 0.8, thetext)
    
%        ylim(ylimits)
    xlim(xlimits)
   
    set(gca,'ytick',0:0.1:0.2)
set(gca,'xtick',0:0.15:0.15)
    
    h.FaceColor = 'k';
    h.EdgeColor = 'w';
     xlabel('Wavelet Variance (10^4) (^c)^2')
    ylabel('Density')
    set(gca,'fontweight','bold')
    set(gca,'fontsize',fs)
    %title(datatitle(ii,:),'fontsize',fs*titleboost)
    hold all
    
end


% % waves figure
% for i = 1:3
% figure
%
%
% plot




