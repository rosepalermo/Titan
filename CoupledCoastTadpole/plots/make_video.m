% make video of lake through time

figure()

v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/test1');
open(v);
k = 0;
for k=1:length(output)
    
    imagesc(output(:,:,k)>p.sealevel_init)
    set(gca,'YDir','Normal')
    frame = getframe(gcf);
    writeVideo(v,frame);
    axis equal
end
close(v);