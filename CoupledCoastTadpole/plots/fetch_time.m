% make plots fetch dist through time for wave and uniform
time = g.t;

figure()
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/091820/rand_uniform_Kc0_0015');
open(v);
if p.doUniformErosion && ~p.doWaveErosion % need to calculate wave matrix
    wave_weight_matrix = zeros(size(g.output));
    p.Kwave = p.Kuniform;
    ind = 1:size(g.output,3);
for i = 1:length(ind)
    [~,wave_weight_matrix(:,:,i),indshoreline_ordered,~,~] = get_dam_wave((g.output(:,:,i)<=p.sealevel_init),p);
    [N,edges] = histcounts(wave_weight_matrix(:,:,i));
    subplot(1,2,1)
    hold on
    plot(edges(2:end),N)
    xlim([0 1])
    
    subplot(1,2,2)
    hold on
    imagesc(wave_weight_matrix(:,:,i))
    colorbar
    caxis([0 0.6])
    axis equal
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end

else
    ind = 1:size(g.output,3);
for i = 1:length(ind)
    subplot(1,2,1)
    hold on
    [N,edges] = histcounts(g.wave_matrix_save(:,:,i));
    plot(edges(2:end),N)
    xlim([0 1])
    
    subplot(1,2,2)
    hold on
    imagesc(g.wave_matrix_save(:,:,i))
    colorbar
    caxis([0 0.6])
    axis equal
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
end

