% function [FetchArea,normFetchArea] = fetch_cell(shoreline_fetch);

%%
%ROSE FIGURE OUT WHAT TO DO TO CLOSE THE POLYGONS BECAUSE YOU DELETED x =
%[x;x(1)]; and y = [y;y(1)];



%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%

x_tot = NaN;
y_tot = NaN;
for ii = 1:length(shoreline_fetch) %add up all of the shoreline pieces (islands too)
    x_tot = [x_tot; shoreline_fetch{ii,1}(:,1)];
    y_tot = [y_tot; shoreline_fetch{ii,1}(:,2)];
end
x_tot = x_tot(2:end); y_tot = y_tot(2:end); % get rid of the NaN
x_ = x_tot - mean(x);
y_ = y_tot - mean(y);

length_cell=cellfun(@length,shoreline_fetch);
bounding_idx = length_cell == max(length_cell);
slbounding = shoreline_fetch(length_cell == max(length_cell));
x_bounding = slbounding{1,1}(:,1);
y_bounding = slbounding{1,1}(:,2);

D = 2*max(sqrt(x_.^2 + y_.^2)) + 10;

%find 1/2 the smallest dist from one point to the next
eps = 0.5*sqrt(min((x_tot(1:end-1)-x_tot(2:end)).^2 + (y_tot(1:end-1)-y_tot(2:end)).^2));

% angles around from point
theta = linspace(0,2*pi,1000);%1000 steps looked good

%find lake area - island area to get total possible fetch area
FetchArea = zeros(size(x_tot));
islandarea = 0;
for ii = find(~bounding_idx)
    islandarea = islandarea + polyarea(shoreline_fetch{ii,1}(:,1),shoreline_fetch{ii,1}(:,2));
end
LakeArea = polyarea(x_bounding,y_bounding)-islandarea;

x = x_tot; y = y_tot;

% % plot lake area, rays that extend from each point, and fill polygon
figure(1); clf
plot(x,y,'k','LineWidth',2); hold on
axis tight
axis equal


for j=1:length(x)-1
    j
    %add point j to beginning and end of fetch poly to close it
    fetch_poly_x = zeros(size(theta+2));
    fetch_poly_y = zeros(size(theta+2));
    fetch_poly_x(1) = x(j);
    fetch_poly_y(1) = y(j);
    fetch_poly_x(end) = x(j);
    fetch_poly_y(end) = y(j);
    
    %     for i=2:length(theta)
    % compute a ray
    %         dx = D*cos(theta(i));
    %         dy = D*sin(theta(i));
    dx = D*cos(theta);
    dy = D*sin(theta);
    
    
    % check if a small ray is in or on the polygon
    dx_eps = eps * cos(theta);
    dy_eps = eps * sin(theta);
    
   in_island = zeros(size(theta)); on_island = zeros(size(theta));
   
    if length(shoreline_fetch) > 1
        for ii = find(~bounding_idx) 
            [in_add, on_add] = inpolygon(x(j)+dx_eps,y(j)+dy_eps,shoreline_fetch{ii,1}(:,1),shoreline_fetch{ii,1}(:,2));
            in_island = in_island | in_add; on_island = on_island | on_add;
        end
    end
    
    
    
    %inside of the bounding
    [in_bounding, on_bounding] = inpolygon(x(j)+dx_eps,y(j)+dy_eps,x_bounding,y_bounding);
    
    in = ~in_island & in_bounding; on = ~on_island & on_bounding; % add them both up!
    
    
    %if it's not in or on, make it point j
    %         if ~any(in | on)
    %             fetch_poly_x(i) = x(j);
    %             fetch_poly_y(i) = y(j);
    %             continue
    %         end
    
    % for the ones that aren't in or on, make it point j
    fetch_poly_x(~(in|on)) = x(j);
    fetch_poly_y(~(in|on)) = y(j);
    ind_inon = find(in|on);
    
    
    % find the first intersection between the ray and the polygon,
    % using min distance
    
    for i=2:length(ind_inon)
        x1 = [x(j), x(j) + dx(ind_inon(i))];
        y1 = [y(j), y(j) + dy(ind_inon(i))];
        
        xi = NaN; yi = NaN; %get something to start with
        for ii = 1:1:length(shoreline_fetch)
            [xi_add,yi_add] = polyxpoly(x1,y1,shoreline_fetch{ii,1}(:,1),shoreline_fetch{ii,1}(:,2)); % find intersection between ray and each shoreline polygon
            xi = [xi;xi_add]; %make an array of the x and y of all intersections
            yi = [yi;yi_add];
        end
        xi = xi(2:end); yi = yi(2:end); %remove the NaN
        
        
        [m,k] = min((x(j)-xi).^2 + (y(j)-yi).^2);
        while m < eps^2
            xi(k) = [];
            yi(k) = [];
            [m,k] = min((x(j)-xi).^2 + (y(j)-yi).^2);
        end
        %if there is no intersection (lines are parallel and on top of each
        %other), then use point j
        if isempty(xi)
            fetch_poly_x(ind_inon(i)) = x(j);
            fetch_poly_y(ind_inon(i)) = y(j);
            continue
        end
        
        % point int = points that make up fetch area polygon
        x_int = xi(k);
        y_int = yi(k);
        fetch_poly_x(ind_inon(i)) = x_int;
        fetch_poly_y(ind_inon(i)) = y_int;
%                 plot([x(j) x_int],[y(j) y_int]); drawnow   % plot all intersecting rays
    end
%     clf; plot(x_bounding,y_bounding); hold on

    % get rid of zeros (for debugging)
    idx = (fetch_poly_x ~= 0 | fetch_poly_y ~= 0);
    fetch_poly_x = fetch_poly_x(idx);
    fetch_poly_y = fetch_poly_y(idx);
    
    %plot filled fetch polygons
%     fill(fetch_poly_x, fetch_poly_y,'b'); drawnow
%     hold on
%     scatter(x(j),y(j),'o')
    
    
    % compute fetch area
    FetchArea(j)=polyarea(fetch_poly_x,fetch_poly_y);
    %     break
end

% normalize fetch area by total lake area
normFetchArea = FetchArea./LakeArea;


% % plot norm fetch area on point from which was calculated
figure(2)
axis tight
axis equal
scatter3(x,y,normFetchArea,[],normFetchArea,'.')
colormap(parula)
colorbar
view(2)


% %save to csv
% save=[x y normFetchArea];
% % csvwrite(savename,save)

% figure(3)
% scatter3(x,y,FetchArea,[],FetchArea,'.')
% colormap(parula)
% colorbar
% view(2)
% axis tight
% axis equal




% fetchper=normFetchArea;
% movav=tsmovavg(fetchper,'t',2,1); %window=1/10 of lake points
% figure(4)
% scatter3(x,y,movav,[],movav,'.')
% colormap(parula)
% colorbar
% view(2)
% axis tight
% axis equal

