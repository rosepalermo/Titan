% make plots of things changing through time

% area through time
ind = 1:length(output);
for i = 1:32
nLakeCells(i) = length(find(output(:,:,ind(i))<=p.sealevel_init));
end

pct_Ao = nLakeCells./p.Ao_cells;
time = ind*p.dtmax;

figure()
plot(time,pct_Ao,'k','LineWidth',2)
set(gca,'FontSize',14)
xlabel('time')
ylabel('Lake area/A_o area')


figure()
for i = 1:length(ind)
subplot(2,5,i)
temp = g.dam_uniform_save(:,:,i);
temp(temp==0) =nan;
histogram(temp)
end