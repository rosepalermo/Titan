% uniform 8con - 11
% uniform 4con - 11
% wave 8con - 224
% wave 4con - 274


% uniform 8con - 11
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/modelexamples/river_uniform_con81_Kc0_1.mat')

ind = 1:size(g.output,3);
nLakeCells_u8 = zeros(1,length(ind));
for i = 1:length(ind)
nLakeCells_u8(i) = length(find(g.output(:,:,ind(i))<=p.sealevel_init));
end

% uniform 4con - 11
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/modelexamples/river_uniform_Kc0_1.mat')

ind = 1:size(g.output,3);
nLakeCells_u4 = zeros(1,length(ind));
for i = 1:length(ind)
nLakeCells_u4(i) = length(find(g.output(:,:,ind(i))<=p.sealevel_init));
end

%wave - 8con
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/modelmethods/modelexamples/river_wave_con81_Kc0_1.mat')


ind = 1:size(g.output,3);
nLakeCells = zeros(1,length(ind));
for i = 1:length(ind)
nLakeCells(i) = length(find(g.output(:,:,ind(i))<=p.sealevel_init));
end

pct_Ao = nLakeCells./nLakeCells(1);
time = 1:length(g.t);
figure(1)
hold on
if p.doWaveErosion
    plot(time,pct_Ao,'r--','LineWidth',2)
elseif p.doUniformErosion
    plot(time,pct_Ao,'k--','LineWidth',2)
end
set(gca,'FontSize',14)
xlabel('time')
ylabel('Lake area relative to init')