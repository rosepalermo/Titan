% upload Croatia1

C = kml2struct('/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/FLCrookedLake_v2Kml.kml');

% for i  = 1:length(C)
%     plot(C(i).Lat,C(i).Lon)
%     hold on
% end

% figure()
[x,y,utmzone] = deg2utm(C(1).Lat,C(1).Lon);
[lake,X,Y] = gridlake(x',y',1,1,10);

figure()
imagesc(X(1,:)',Y(:,1),lake)
hold on
plot(x,y,'r','LineWidth',2)
p.con8 = 1;
[indshoreline_ocw,~,~,~] = order_shoreline_bwbound(lake,p);
save('FLCrookedLake_v2_sl','indshoreline_ocw','lake')
