% make video of lake through time

figure()

v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/save_more/rand_uniform_Kc0_0001');
open(v);
k = 0;
for k=1:size(g.output,3)
    
    imagesc(g.output(:,:,k)>p.sealevel_init)
    set(gca,'YDir','Normal')
    frame = getframe(gcf);
    writeVideo(v,frame);
    axis equal
end
close(v);