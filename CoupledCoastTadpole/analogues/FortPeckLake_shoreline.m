uiopen('/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/FortPeckLake.tif',1)
FortPeckLake = double(FortPeckLake);
FortPeckLake_lake_tile = FortPeckLake>75;
p.con8 = 1;
[F_lake_all,~,~,~] = find_first_order_lakes(FortPeckLake_lake_tile,p);
% parfor ff = 1:171%length(F_lake_all) %% breaks after 200..
%     nlakecells(ff) = length(find(F_lake_all{ff}));
% end
%%
FortPeckLake_lake = F_lake_all{171};
[indshoreline_ocw,~,~,~] = order_shoreline_bwbound(FortPeckLake_lake,p);

%%
imagesc(lake)
hold on
for i = 1:length(indshoreline_ocw)
    plot(indshoreline_ocw{i}(:,2),indshoreline_ocw{i}(:,1),'r','LineWidth',2)
end