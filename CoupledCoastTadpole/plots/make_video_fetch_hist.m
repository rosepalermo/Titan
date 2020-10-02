function make_video_fetch_hist(v,p,g)
folder = fileparts(which('getpath_CCT.m'));
addpath(genpath(folder));
open(v);
    wave_weight_matrix = zeros(size(g.output));
    ind = 1:size(g.output,3);
if ~isfield(g,'wave_matrix_save') % need to calculate wave matrix
    p.Kwave = p.Kuniform;
for i = 1:length(ind)
    [~,wave_weight_matrix(:,:,i),indshoreline_ordered,~,~] = get_dam_wave((g.output(:,:,i)<=p.sealevel_init),p);
    histogram(wave_weight_matrix(:,:,i)./max(wave_weight_matrix(:)),[0:0.001:1],'DisplayStyle','stairs');
    hold on
    time_n = strrep(num2str(g.t(ind(i))),'.','_');
    title(string(strcat('time',' ',string(time_n))))
    frame = getframe(gcf);
    writeVideo(v,frame);
end
else
    for i = 1:length(ind)
    histogram(g.wave_matrix_save(:,:,i)./max(max(g.wave_matrix_save(:,:,i))),[0:0.001:1],'DisplayStyle','stairs');
    hold on
    time_n = strrep(num2str(g.t(ind(i))),'.','_');
    title(string(strcat('time',' ',string(time_n))))
        xlim([0 1])
    frame = getframe(gcf);
    writeVideo(v,frame);
    end


close(v);
end
end