% plot end shoreline from data
function plotendsl(p,g,i)
lake = g.output(:,:,i)<=p.sealevel_init;
shoreline = addidshoreline(lake,~lake);
indsl = find(shoreline);
[y,x] = ind2sub(size(shoreline),indsl);
scatter(x,y,'.')
set(gca,'Ydir','reverse')
xlim([0 p.Nx])
ylim([0 p.Ny])
end