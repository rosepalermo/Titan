uiopen('/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/Butler_chain.tif',1)
butler = double(Butler_chain);
butler_lake_tile = butler>20;
p.con8 = 1;
[F_lake_all,~,~,~] = find_first_order_lakes(butler_lake_tile,p);
% for ff = 1:length(F_lake_all) %% breaks after 200..
% imagesc(F_lake_all{ff})
% ff
% pause
% end
%%
butler_lake = F_lake_all{339};
[indshoreline_ocw,~,~,~] = order_shoreline_bwbound(butler_lake,p);
imagesc(butler_lake)
hold on
for i = 1:length(indshoreline_ocw)
    plot(indshoreline_ocw{i}(:,2),indshoreline_ocw{i}(:,1),'r','LineWidth',2)
end