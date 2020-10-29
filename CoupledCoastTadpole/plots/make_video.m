% make video of lake through time

figure()

v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_wave_Kc1_zoom_video');
open(v);
k = 0;
for k=1:size(g.output,3)
    for ii = 1:10
    imagesc(g.output(:,:,k)>p.sealevel_init)
    set(gca,'YDir','Normal')
    set(gca,'Ylim',[45  110])
    set(gca,'Xlim',[60  140])
    axis equal
        colormap gray
    frame = getframe(gcf);
    writeVideo(v,frame);

%     axis equal
    end
end
close(v);