uiopen('/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/Sebago.tif',1)
sebago = double(Sebago);
sebago_lake_tile = sebago>20;
p.con8 = 1;
[F_lake_all,~,~,~] = find_first_order_lakes(sebago_lake_tile,p);
sebago_lake = F_lake_all{31};
[indshoreline_ocw,~,~,~] = order_shoreline_bwbound(sebago_lake,p);
imagesc(sebago_lake)
hold on
for i = 1:length(indshoreline_ocw)
    plot(indshoreline_ocw{i}(:,2),indshoreline_ocw{i}(:,1),'r','LineWidth',2)
end
axis equal