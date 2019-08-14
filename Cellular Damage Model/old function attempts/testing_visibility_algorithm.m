% testing visibility algorithm

addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Cellular Damage Model')
addpath('/Users/rosepalermo/Documents/GitHub/VisiLibity1/src')
load wavet2v2_lake.mat
load test_vis_polygon_updated.mat
% [observer_y,observer_x] = ind2sub(size(lake),find(lake));
x = fetch_sl_cells{1,1}(:,1);
y = fetch_sl_cells{1,1}(:,2);
% I need to find a point eps inside of the lake at every lake point..
% dx_eps = (xend - x(iv))*1e-8;
% dy_eps = (yend - y(iv))*1e-8;
% xseg = [ones(size(theta))*x(iv)+dx_eps; xend; nan(size(theta))];
% yseg = [ones(size(theta))*y(iv)+dy_eps; yend; nan(size(theta))];

epsilon = 1e-5;
snap_distance = 0.05;

environment = {flipud([fetch_sl_cells{1,1};fetch_sl_cells{1,1}(1,:)])};

x_wrapped = [x; x(1,:)];
y_wrapped = [y; y(1,:)];
diff = [x_wrapped(2:end), y_wrapped(2:end)] - [x, y];
perp = [diff(:,2), -diff(:,1)];
norm = sqrt(diff(:,1).^2 + diff(:,2).^2);
eps = .001;
perp = eps*perp./norm;
x_inside = x + perp(:,1);
y_inside = y + perp(:,2);



tic
for l = 1:length(x_inside)
V{l} = visibility_polygon( [x_inside(l) y_inside(l)] , environment , epsilon , snap_distance );
end
toc
%%

% load('test_polygon_island.mat')
% lake = repmat(polygon,2);

% land = ~lake;
% x = 1:length(lake); y = 1:length(lake);
% [X,Y] = meshgrid(x,y);
% [shoreline] = addidshoreline(lake,land); % corners are part of the shoreline!
% i =1;
% [L,fol] = find_first_order_lakes(lake);
% for ff = 2:2%length(fol)
%     
%     F_lake = (L == fol(ff));
%     if (exist('erodedind','var')) | (i == 1) | (ff>1)
%         disp('fetch')
%         clearvars fetch_sl_cells indshoreline WaveArea_cell
%         
%         %order the shoreline and islands
%         %             shoreline = addidshoreline_cardonly(lake,land); %rewrite shoreline to be card only for fetch.. will damage corners separately later
%         shoreline = addidshoreline_cardonly(F_lake,~F_lake); %rewrite shoreline to be card only for fetch.. will damage corners separately later
%         disp('ordering')
%         tic
%         [indshoreline_ocw,cells2trash] = order_cw_lastpoint(F_lake,shoreline); % ccw ordered ind = indshoreline
%         disp('ordered')
%         toc
%         for l = 1: length(indshoreline_ocw)
%             indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{l,1}(:,1),indshoreline_ocw{l,1}(:,2));
%             fetch_sl_cells{1,l}(:,1) = X(indshoreline{l,1});
%             fetch_sl_cells{1,l}(:,2) = Y(indshoreline{l,1});
%         end
%     end
% end
% 
% figure();
% imagesc(lake)
% 
% %Robustness constant
% epsilon = 1e-5;
% 
% 
% %Snap distance (distance within which an observer location will be snapped to the
% %boundary before the visibility polygon is computed)
% snap_distance = 0.05;
% 
% environment = fetch_sl_cells(1,1);
% observer_x = 300; observer_y = 350;
% addpath('/Users/rosepalermo/Documents/GitHub/VisiLibity1/src')
% figure()
% %Plot Environment
% patch( environment{1}(:,1) , environment{1}(:,2) , 0.1*ones(1,length(environment{1}(:,1)) ) , ...
%        'w' , 'linewidth' , 1.5 );
% for i = 2 : size(environment,2)
%     patch( environment{i}(:,1) , environment{i}(:,2) , 0.1*ones(1,length(environment{i}(:,1)) ) , ...
%            'k' , 'EdgeColor' , [0 0 0] , 'FaceColor' , [0.8 0.8 0.8] , 'linewidth' , 1.5 );
% end
% hold on
% %Plot observer
%     plot3( x_inside(100) , y_inside(100) , 0.3 , ...
%            'o' , 'Markersize' , 9 , 'MarkerEdgeColor' , 'y' , 'MarkerFaceColor' , 'k' );
% 
% %Compute and plot visibility polygon
% V{1} = visibility_polygon( [x_inside(100) y_inside(100)] , environment , epsilon , snap_distance );
%     patch( V{1}(:,1) , V{1}(:,2) , 0.1*ones( size(V{1},1) , 1 ) , ...
%            'r' , 'linewidth' , 1.5 );
%     plot3( V{1}(:,1) , V{1}(:,2) , 0.1*ones( size(V{1},1) , 1 ) , ...
%            'b*' , 'Markersize' , 5 );
%        
%        axis equal