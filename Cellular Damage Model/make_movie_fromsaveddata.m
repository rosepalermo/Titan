% addpath('D:\Titan\Modeling\AGU final folder')
addpath('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/River_and_wave_9_19/')
load('wavet2v110.mat')
% make a movie from the figures
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/River_and_wave_9_19/wavet2v110_zoom');
open(v);
for k = 1:10
%     lake = zeros(size(X));
%     lake(shoreline_save{k,1}) = 1;
%     Xinon = reshape(X,[],1);
%     Yinon = reshape(Y,[],1);
%     % [lakex,lakey] = ind2sub(size(X),);
%     [in, on] = inpoly([Xinon,Yinon]',[ordered_sl_save{k,1}{1,1}]');
%     lake = in + on;
%     lake = ~reshape(lake,size(X));
    
    imagesc(lake_save{k})
    colormap gray
    set(gca,'YDir','Normal')
    set(gca,'XLim',([2500 3700])); set(gca,'YLim',([300 1300])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);
%%
% addpath('D:\Titan\Modeling\AGU final folder')
% load('uniformt1v1.mat')
% % make a movie from the figures
% v = VideoWriter('uniform_t1v1_zoom');
% open(v);
% for k = 1:20
%     imagesc(~lake_save{k,1})
%     colormap gray
%     set(gca,'YDir','Normal')
%     set(gca,'XLim',([675 825])); set(gca,'YLim',([180 300])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
%     frame = getframe(gcf);
%     writeVideo(v,frame);
%     
% end
% close(v);
% 
% 
% %% Rivers
% 
% addpath('D:\Titan\Modeling\AGU final folder')
% load('wavet1v1.mat')
% % make a movie from the figures
% v = VideoWriter('rivers_v1_zoom');
% open(v);
% lake = zeros(size(X));
% lake(shoreline_save{1,1}) = 1;
% Xinon = reshape(X,[],1);
% Yinon = reshape(Y,[],1);
% % [lakex,lakey] = ind2sub(size(X),);
% [in, on] = inpoly([Xinon,Yinon]',[ordered_sl_save{1,1}{1,1}]');
% lake = in + on;
% lake = ~reshape(lake,size(X));
% 
% imagesc(lake)
% colormap gray
% set(gca,'YDir','Normal')
% set(gca,'XLim',([675 825])); set(gca,'YLim',([180 300])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
% frame = getframe(gcf);
% writeVideo(v,frame);
% % clearvars -except v
% 
% load('uniformt2v1.mat')
% imagesc(~lake_save{1,1})
% colormap gray
% set(gca,'YDir','Normal')
% set(gca,'XLim',([840 1010])); set(gca,'YLim',([270 415])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
% frame = getframe(gcf);
% 
% 
% writeVideo(v,frame);
% close(v);