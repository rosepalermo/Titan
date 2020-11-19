% plot end shoreline from data
function [x,y] = plotisl(p,g,i)
if i == 'end'
    lake = g.output(:,:,end)<=p.sealevel_init;
else
    lake = g.output(:,:,i)<=p.sealevel_init;
end
[F_lake_all,~,~,~] = find_first_order_lakes(lake,p);
for i = 1:length(F_lake_all)
lengthF(i) = length(find(F_lake_all{i}));
end
[~,indmaxL] = max(lengthF);
lake = F_lake_all{indmaxL};
p.sl_analysis = 1;
[sl_cell,~,~,~] = order_shoreline_bwbound(lake,p);
for i = 1:length(sl_cell)
    sl_cell_length = length(sl_cell{i});
    y = sl_cell{i}(:,1);
    x = sl_cell{i}(:,2);
%     plot(x,y)
%     hold on
end
[~,indmaxslcell] = max(sl_cell_length);
y = sl_cell{indmaxslcell}(:,1);
x = sl_cell{indmaxslcell}(:,2);
% set(gca,'Ydir','reverse')
% xlim([0 p.Nx])
% ylim([0 p.Ny])
end