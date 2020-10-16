function plotICsl(p,g,colors_)
lake = g.output(:,:,1)<=p.sealevel_init;
[F_lake_all,~,~,~] = find_first_order_lakes(lake);
for i = 1:length(F_lake_all)
lengthF(i) = length(find(F_lake_all{i}));
end
[~,indmaxL] = max(lengthF);
lake = F_lake_all{indmaxL};
[sl_cell,~,~,~] = order_shoreline_bwbound(lake,p);
for i = 1:length(sl_cell)
    y = sl_cell{i}(:,1);
    x = sl_cell{i}(:,2);
    plot(x,y,colors_)
    hold on
end
% shoreline = addidshoreline(lake,~lake);
% indsl = find(shoreline);
% [y,x] = ind2sub(size(shoreline),indsl);
% scatter(x,y,'.')
set(gca,'Ydir','reverse')
xlim([0 p.Nx])
ylim([0 p.Ny])
end