% make plots fetch dist through time for wave and uniform
if p.doUniformErosion
    time = g.t;
end
if p.doWaveErosion
    time = g.t/10;
end
figure()
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/save_more/rand_uniform_Kc0_0001_fetch');
open(v);
if p.doUniformErosion % need to calculate wave matrix
    wave_weight_matrix = zeros(size(g.output));
    p.Kwave = p.Kuniform;
    ind = 1:10:size(g.output,3);
for i = 1:length(ind)
    [~,wave_weight_matrix(:,:,i),indshoreline_ordered,~,~] = get_dam_wave((g.output(:,:,i)<=p.sealevel_init),p);
    [N,edges] = histcounts(wave_weight_matrix(:,:,i));
    plot(edges(2:end),N)
    hold on
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
end


figure()
v = VideoWriter('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/save_more/rand_wave_Kc0_0001_fetch');
open(v);
    wave_weight_matrix = zeros(size(g.output));
if p.doWaveErosion % need to calculate wave matrix
    ind = 1:size(g.output,3);
for i = 1:length(ind)
    [~,wave_weight_matrix(:,:,i),indshoreline_ordered,~,~] = get_dam_wave((g.output(:,:,i)<=p.sealevel_init),p);
    [N,edges] = histcounts(wave_weight_matrix(:,:,i));
    plot(edges(2:end),N)
    hold on
    time_n = strrep(num2str(time(ind(i))),'.','_');
    title(string(strcat('time',' ',string(time_n))))
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
end