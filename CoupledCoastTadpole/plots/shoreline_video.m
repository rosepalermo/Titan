% make a shoreline movie

figure()

v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/091820/rand_hwave_luniform_sl');

open(v);
k = 0;
for k=1:length(sl_save)
    
    plot(sl_save{k}(:,1),sl_save{k}(:,2),'k','LineWidth',1.5)
    set(gca,'YDir','Normal')
    frame = getframe(gcf);
    writeVideo(v,frame);
    axis equal
end
close(v);