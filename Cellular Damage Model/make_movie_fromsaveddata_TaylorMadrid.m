% addpath('D:\Titan\Modeling\AGU final folder')
addpath('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/River_and_wave_9_19/')
% load('wavet2v110.mat')
load('makemovies_9_19_w18_u15.mat')
% make a movie from the figures
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/River_and_wave_9_19/wavet2v1_more20');
open(v);
k = 0;
for k=1:18
    for kk=1:10
    
    imagesc(wave{k})
    colormap gray
    set(gca,'YDir','Normal')
    frame = getframe(gcf);
    writeVideo(v,frame);
    axis equal
    end
end
close(v);

% make a movie from the figures
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/River_and_wave_9_19/uniformt2v1_more20');
open(v);
% k = 0;
for k=1:15
    for kk=1:10
    imagesc(uniform{k})
    colormap gray
    set(gca,'YDir','Normal')
    frame = getframe(gcf);
    writeVideo(v,frame);
    axis equal
    end
end
close(v);
