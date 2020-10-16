% make video of lake through time

figure()

v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/101520/river_uniform_Kc0_1_video');
open(v);
k = 0;
for k=1:size(g.output,3)
    
    imagesc(g.output(:,:,k)>p.sealevel_init)
    set(gca,'YDir','Normal')
%     set(gca,'Xlim',[95.7765  152.9194])
%     set(gca,'Ylim',[76.5584  147.7993])
    axis equal
    
    frame = getframe(gcf);
    writeVideo(v,frame);
%     axis equal
end
close(v);