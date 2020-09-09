% make plots of things changing through time

% area through time
ind = 1:size(g.output,3);
nLakeCells = zeros(1,length(ind));
for i = 1:length(ind)
nLakeCells(i) = length(find(g.output(:,:,ind(i))<=p.sealevel_init));
end

pct_Ao = nLakeCells./p.Ao_cells;
time = g.t;

% % figure()
% hold on
% plot(time,pct_Ao,'r','LineWidth',2)
% set(gca,'FontSize',14)
% xlabel('time')
% ylabel('Lake area/A_o area')

neroded = zeros(size(pct_Ao));
for i = 2:length(ind)
    neroded(i) = nLakeCells(i) - nLakeCells(i-1);
end
figure(1)
hold on
if p.doWaveErosion
    plot(time,neroded,'r','LineWidth',2)
elseif p.doUniformErosion
    plot(time,neroded,'k','LineWidth',2)
end
set(gca,'FontSize',14)
xlabel('time')
ylabel('Lake area/A_o area')


%calculate mean damage through time
if p.doWaveErosion
    for i = 1:length(ind)-1
        temp = g.dam_wave_save(:,:,i);
        temp(temp==0) =nan;
        meandam(i) = nanmean(temp,'all');
    end
end
if p.doUniformErosion
    for i = 1:length(ind)-1
        temp = g.dam_uniform_save(:,:,i);
        temp(temp==0) =nan;
        meandam(i) = nanmean(temp,'all');
    end
end
figure(2)
hold on
if p.doWaveErosion
    plot(time(1:end-1), meandam,'r','LineWidth',2)
elseif p.doUniformErosion
    plot(time(1:end-1), meandam,'k','LineWidth',2)
end
set(gca,'FontSize',14)
xlabel('time')
ylabel('mean damage')

%% calculate mean wave weighting through time
if p.doWaveErosion
    for i = 1:length(ind)-1
        temp = g.wave_matrix_save(:,:,i);
        temp(temp==0) =nan;
        meandam(i) = nanmean(temp,'all');
    end
end
figure(3)
hold on
if p.doWaveErosion
    plot(time(1:end-1), meandam,'r','LineWidth',2)
elseif p.doUniformErosion
    plot(time(1:end-1), meandam,'k','LineWidth',2)
end
set(gca,'FontSize',14)
xlabel('time')
ylabel('mean wave weight')

% figure()
% for i = 1:10:(length(ind))
% subplot(2,ceil(ind(end)/10),i)
% temp = g.dam_wave_save(:,:,i);
% temp(temp==0) =nan;
% histogram(temp)
% end