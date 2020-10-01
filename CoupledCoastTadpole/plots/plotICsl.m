function plotICsl(p,g)
lake = g.output(:,:,1)<=p.sealevel_init;
shoreline = addidshoreline(lake,~lake);
indsl = find(shoreline);
[y,x] = ind2sub(size(shoreline),indsl);
scatter(x,y,'.')
set(gca,'Ydir','reverse')
xlim([0 p.Nx])
ylim([0 p.Ny])
end