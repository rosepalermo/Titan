
% old wave run
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_wave_Kc0_15.mat');
figure()
for i = 1:size(g.output,3);
    imagesc(g.wave_matrix_save(:,:,i))
    caxis([0 1])
    title(i)
    drawnow
end

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015.mat')
figure()
for i = 1:size(g.output,3);
    imagesc(g.dam_uniform_save(:,:,i))
    caxis([0 0.03])
    title(i)
    drawnow
end

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15.mat')
figure()
for i = 1:size(g.output,3);
    imagesc(g.output(:,:,i)<=p.sealevel_init)
%     caxis([0 1])
% plotisl(p,g,i)
    title(i)
    drawnow
end

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/093020/rand_mwave_muniform.mat')
