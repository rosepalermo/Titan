function [FetchArea,normFetchArea] = fetch_xy_original(x,y);

% %get rid of duplicate points (when it goes exactly around a pixel)
% % indices to unique values in column 3
% [~, ind] = unique(M(:, 4:5), 'rows');
% % duplicate indices
% duplicate_ind = setdiff(1:size(M, 1), ind);
% % duplicate values
% duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
% M(duplicate_ind,:)=[];
% x = M(:,4);
% y = M(:,5);
% x=x(1:5:end);
% y=y(1:5:end);

%add the first point to the end to make periodic--make sure to take off
%last point in dowave if you're leaving this.
% x = x'; y = y';
% x = [x;x(1)];
% y = [y;y(1)];

x_ = x - mean(x);
y_ = y - mean(y);
D = 2*max(sqrt(x_.^2 + y_.^2)) + 10;

%find 1/2 the smallest dist from one point to the next
eps = 0.5*sqrt(min((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2));

% angles around from point
theta = linspace(0,2*pi,1000);%1000 steps looked good

%zeros for fetch area
FetchArea = zeros(size(x));
LakeArea = polyarea(x,y);

% % plot lake area, rays that extend from each point, and fill polygon
% figure(1); clf
% plot(x,y,'k','LineWidth',2); hold on
% axis tight
% axis equal


parfor j=1:length(x)-1
    j;
    %add point j to beginning and end of fetch poly to close it
    fetch_poly_x = zeros(size(theta+2));
    fetch_poly_y = zeros(size(theta+2));
    fetch_poly_x(1) = x(j);
    fetch_poly_y(1) = y(j);
    fetch_poly_x(end) = x(j);
    fetch_poly_y(end) = y(j);
    
    for i=2:length(theta)
        % compute a ray
        dx = D*cos(theta(i));
        dy = D*sin(theta(i));

        % check if a small ray is in or on the polygon
        dx_eps = eps*dx/norm([dx;dy]);
        dy_eps = eps*dy/norm([dx;dy]);
        [in, on] = inpolygon(x(j)+dx_eps,y(j)+dy_eps,x,y);
        %if it's not in or on, make it point j
        if ~any(in | on)
            fetch_poly_x(i) = x(j);
            fetch_poly_y(i) = y(j);
            continue
        end
        
        % find the first intersection between the ray and the polygon,
        % using min distance
        x1 = [x(j), x(j) + dx];
        y1 = [y(j), y(j) + dy];
        
        [xi,yi] = polyxpoly(x1,y1,x,y);
        [m,k] = min((x(j)-xi).^2 + (y(j)-yi).^2);
        while m < eps^2
            xi(k) = [];
            yi(k) = [];
            [m,k] = min((x(j)-xi).^2 + (y(j)-yi).^2);
        end
        %if there is no intersection (lines are parallel and on top of each
        %other), then use point j
        if isempty(xi)
            fetch_poly_x(i) = x(j);
            fetch_poly_y(i) = y(j);
            continue
        end
        
        % point int = points that make up fetch area polygon
        x_int = xi(k);
        y_int = yi(k);
        fetch_poly_x(i) = x_int;
        fetch_poly_y(i) = y_int;
%         plot([x(j) x_int],[y(j) y_int])   % plot all intersecting rays
    end
    
    % get rid of zeros (for debugging)
    idx = (fetch_poly_x ~= 0 | fetch_poly_y ~= 0);
    fetch_poly_x = fetch_poly_x(idx);
    fetch_poly_y = fetch_poly_y(idx);
    
    %plot filled fetch polygons
%     fill(fetch_poly_x, fetch_poly_y,'b');% drawnow

    % compute fetch area
    FetchArea(j)=polyarea(fetch_poly_x,fetch_poly_y);
%     break
end

% normalize fetch area by total lake area
normFetchArea = FetchArea./LakeArea;


% % plot norm fetch area on point from which was calculated
% figure(2)
% axis tight
% axis equal
% scatter3(x,y,normFetchArea,eps,normFetchArea,'.')
% colormap(parula)
% colorbar
% view(2)


% %save to csv
% save=[x y normFetchArea];
% % csvwrite(savename,save)

% figure(3)
% scatter3(x,y,normFetchArea,[],normFetchArea,'.')
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

