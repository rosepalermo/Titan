%% file to look at histograms from all shorelines and compute some
% statistics

clear
close all

% load the file with all of the small-scale variance data
load alleq14_2

% some general controls

fs = 20;    % font size
edgespace = .15 ;
edgemax = 4 ;
radscale = 10^-4;

% this is the shoreline that you get. in short:
% 1 rednoise_v1_';
datatitle(1,:) = string('Red Noise Initial Conditions');
datatitle(2,:) = 'Early Waves';
datatitle(3,:) ='Later Waves';
datatitle(4,:) = 'Early Uniform Erosion';
datatitle(5,:) = 'Later Uniform Erosion';
datatitle(6,:) = 'Oops Just Rivers';
% 7 LakePowell_';
datatitle(7,:) = 'Ricky Powell';
% 8 scotland_';
datatitle(8,:) = 'Highlands';
% 9 lgm'; - Ligea Mare
datatitle(9,:) = 'Ligea Mare';
% 10 sebago';
datatitle(10,:) = 'Sebago';
datatitle(11,:) = 'West Scotland';
datatitle(12,:) ='East Scotland';

%% main loop


for ii = 1:12
    
    % load a specific set of values
    figure
    vars = cell2mat(eq14save(ii))/radscale; 
    
    % make a histogram with these controls
    
    meanvar(ii) = mean(vars);
    skewvar(ii) = skewness(vars);
    
    edges = (0:edgespace:edgemax);
    
    
    
    h = histogram(vars,edges,'Normalization','probability');
    
    
    
    meantext= ['Mean = ' num2str(meanvar(ii))];
    skewtext= ['Skewness = ' num2str(skewvar(ii))]; 
    thetext = {meantext, skewtext};
    text(edgemax * 0.65, max(h.Values) * 0.8, thetext)
    
   ylim([0 0.5])
    
    
    h.FaceColor = 'k';
    h.EdgeColor = 'w';
    xlabel('Small-scale Wavelet Power Spectrum Variance (10^4) (^c)^2')
    ylabel('Density')
    set(gca,'fontweight','bold')
    set(gca,'fontsize',fs)
    title(datatitle(ii,:),'fontsize',fs*1.5)
    hold all
    
end


% % waves figure
% for i = 1:3
% figure
% 
% 
% plot




