% shorelines for wavelets
close all; clear all
%% uniform only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_uniform_Kc0_015.mat')
for i = 1:size(g.output,3)
    length_uniform(i) = length(find(g.output(:,:,i)<=p.sealevel_init));
end
lake = g.output(:,:,166);



%% wave only
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15.mat')
for i = 1:size(g.output,3)
    length_wave(i) = length(find(g.output(:,:,i)<=p.sealevel_init));
end
lake = g.output(:,:,139);


%% uniform and wave
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_hwave_luniform_2_.mat')
for i = 1:size(g.output,3)
    length_uandw(i) = length(find(g.output(:,:,i)<=p.sealevel_init));
end
lake = g.output(:,:,85);


%% initial condition
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/100220/rand_wave_Kc0_15.mat')

length_ic = ones(1,size(g.output,3))*length(find(g.output(:,:,1)<=p.sealevel_init));
lake = g.output(:,:,1);


%% figure
figure(); hold on; plot(length_ic); plot(length_uandw); plot(length_uniform); plot(length_wave); legend('IC','U = W','U','W')
xlabel('iteration')
ylabel('# lake cells')
