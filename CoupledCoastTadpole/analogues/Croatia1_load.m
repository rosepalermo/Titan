% upload Croatia1

C = kml2struct('/Users/rosepalermo/Documents/Research/Titan/SurfaceWaterTiff/Croatia1_v1.kml');

for i  = 1:length(C)
    plot(C(i).Lat,C(i).Lon)
    hold on
end

legend('1','2','3','4','5','6')
plot(C(1).Lat(1),C(1).Lon(1),'*')
plot(C(2).Lat(1),C(2).Lon(1),'*')
plot(C(3).Lat(1),C(3).Lon(1),'*')
plot(C(4).Lat(1),C(4).Lon(1),'*')
plot(C(5).Lat(1),C(5).Lon(1),'*')
plot(C(6).Lat(1),C(6).Lon(1),'*')


CroatiaLat = [C(1).Lat;C(3).Lat;C(5).Lat;flipud(C(4).Lat);flipud(C(6).Lat);C(2).Lat;C(1).Lat(1)];
CroatiaLon = [C(1).Lon;C(3).Lon;C(5).Lon;flipud(C(4).Lon);flipud(C(6).Lon);C(2).Lon;C(1).Lon(1)];

% figure()
% plot(CroatiaLat,CroatiaLon)
[x,y,utmzone] = deg2utm(CroatiaLat,CroatiaLon);
[lake,X,Y] = gridlake(x',y',1,1,90);

% imagesc(X(1,:)',Y(:,1),lake)
% hold on
% plot(x,y,'r','LineWidth',2)
p.con8 = 1;
[indshoreline_ocw,~,~,~] = order_shoreline_bwbound(lake,p);
save('Croatia1_sl','indshoreline_ocw','lake')
