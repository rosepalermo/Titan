% function [FetchArea] = fetch_mw_cell(shoreline_fetch)

% rng(2)
%
% An example polygon with lots of concavity
th = linspace(0,2*pi,60);
r = 2 + rand(size(th))-0.5;
x_ex = r.*cos(th);
y_ex = r.*sin(th);

shoreline_fetch = cell(2,1);
shoreline_fetch{1,1}(:,1) = x_ex;
shoreline_fetch{1,1}(:,2) = y_ex;
shoreline_fetch{2,1}(:,1) = [0; 0.5; 1; 0.5; 0];
shoreline_fetch{2,1}(:,2) = [0; 0; 0; 1; 0];


% N=100; %1000 data points
% th = linspace(0,2*pi,N)'; %theta from 0 to 2pi
% amp = 2*ones(N,1);
% for i=1:500 %higher numbers here make more higher frequency fluctuations
%     a = rand()-0.5; %random number from -0.5 to 0.5
%     b = rand()-0.5;
%     amp = amp + a/i*cos(i*th) + b/i*sin(i*th); %change circle by a fourier function over it
% end
% x=(amp.*cos(th))'; y = (amp.*sin(th))';

% x_tot = NaN;
% y_tot = NaN;
% for ii = 1:length(shoreline_fetch) %add up all of the shoreline pieces (islands too)
%     x_tot = [x_tot; shoreline_fetch{ii,1}(:,1)];
%     y_tot = [y_tot; shoreline_fetch{ii,1}(:,2)];
% end
% x_tot = x_tot(2:end); y_tot = y_tot(2:end); % get rid of the NaN
% x = x_tot + abs(min(x_tot)); % add min to make all +
% y = y_tot + abs(min(y_tot));


% if the length is greater than 1, there are islands.
length_cell=cellfun(@length,shoreline_fetch);
bounding_idx = length_cell == max(length_cell);
slbounding = shoreline_fetch(length_cell == max(length_cell));
x_bounding = slbounding{1,1}(:,1)'; % apostrophe to make rows
y_bounding = slbounding{1,1}(:,2)';


for obj = 1:length(shoreline_fetch)
    
    x = shoreline_fetch{obj,1}(:,1);
    y = shoreline_fetch{obj,1}(:,2);
    x = x'; y = y';
    x = [x x(1)]; % close it
    y = [y y(1)];
    
    
    nvert = length(x) - 1; % # vertices in polygon, assuming closed
    % Set up rays
    nray = 1000;
    theta = linspace(0,2*pi,nray+1);
    theta = theta(1:end-1); % eliminate duplicate at 2pi
    D = 2*max(sqrt(x.^2 + y.^2));
    eps = D*1e-8;
    % Set up output: light-of-sight intersection point for each vertex/angle
    % pair
    xlos = nan(nray,nvert);
    ylos = nan(nray,nvert);
    % Loop over polygon vertices...
    for iv = 1:nvert
        xend = x(iv) + D.*cos(theta);
        yend = y(iv) + D.*sin(theta);
        dx_eps = (xend - x(iv))*1e-8;
        dy_eps = (yend - y(iv))*1e-8;
        xseg = [ones(size(theta))*x(iv)+dx_eps; xend; nan(size(theta))];
        yseg = [ones(size(theta))*y(iv)+dy_eps; yend; nan(size(theta))];
        [xi,yi] = polyxpoly(xseg(:), yseg(:), x_bounding, y_bounding);
        if length(shoreline_fetch) > 1
            for ii = 2:length(shoreline_fetch)
            [xi_add,yi_add] = polyxpoly(xseg(:), yseg(:), [shoreline_fetch{ii,1}(:,1); shoreline_fetch{ii,1}(1,1)],[shoreline_fetch{ii,1}(:,2);shoreline_fetch{ii,1}(1,2)]);
            xi = [xi;xi_add]; yi = [yi;yi_add];
            end
        end
            
        % Unique intersection points
        xyi = setdiff(unique([xi yi], 'rows'), [x(iv) y(iv)], 'rows');
        xi = xyi(:,1);
        yi = xyi(:,2);
        
        % Figure out which points are the first intersection points along each
        % ray by grouping by angle
        ang = atan2(yi - y(iv), xi - x(iv));
        ang = interp1(wrapToPi(theta), theta, ang, 'nearest'); % remove rounding error
        d = sqrt((yi - y(iv)).^2 + (xi - x(iv)).^2);
        [G, unqang] = findgroups(ang);
        dmin = splitapply(@min, d, G);
        [~,loc] = ismember([unqang dmin], [ang d], 'rows');
        ximin = xi(loc);
        yimin = yi(loc);
        
        %     % Now figure out whether each ray was inside or outside the polygon
        %     % before intersecting.  Can check simply by seeing if bisector is
        %     % inside or outside.
        %     dx = ximin - x(iv);
        %     dy = yimin - y(iv);
        %     isin = inpolygon(x(iv) + dx/2, y(iv) + dy/2, x, y);
        %     ximin = ximin(isin);
        %     yimin = yimin(isin);
        %     unqang = unqang(isin);
        
        % Now figure out whether each ray was inside or outside the polygon
        % before intersecting.  Check if some small dist away from the point is
        % in
        dx = ximin - x(iv);
        dy = yimin - y(iv);
        [in_bounding,on_bounding] = inpolygon(x(iv) + dx/2, y(iv) + dy/2, x_bounding, y_bounding); % is it in the lake?
        
        % if islands, is it in the island?
        in_island = zeros(size(in_bounding)); on_island = zeros(size(on_bounding));
        if length(shoreline_fetch) > 1
            for ii = 2: length(shoreline_fetch)
                [in_add, on_add] = inpolygon(x(iv) + dx*1e-8, y(iv) + dy*1e-8,[shoreline_fetch{ii,1}(:,1); shoreline_fetch{ii,1}(1,1)],[shoreline_fetch{ii,1}(:,2);shoreline_fetch{ii,1}(1,2)]);
                in_island = in_island | in_add; on_island = on_island | on_add;
            end
        end
        
        
        %in lake but not on island
        isin = ~in_island & in_bounding; on = ~on_island & on_bounding;
        
        ximin = ximin(isin);
        yimin = yimin(isin);
        unqang = unqang(isin);
        
        % Replace in proper slot
        [~,loc] = ismember(unqang, theta);
        xlos(loc,iv) = ximin;
        ylos(loc,iv) = yimin;
    end
    
    % For plotting, define each resulting line segment from vertex to
    % line-of-sight point
    xplt = permute(cat(3, x(1:end-1).*ones(nray,1), xlos), [3 1 2]);
    yplt = permute(cat(3, y(1:end-1).*ones(nray,1), ylos), [3 1 2]);
    % To check, let's plot just the rays associated with angle idx=50
    plot(x,y,'k');
    hold on;
    % for ivv = 1:nray
    for ivv = 1:1
        for iv = 1:nvert
            plot(xplt(:,ivv,iv), yplt(:,ivv,iv));
            hold on
        end
    end
    
    figure()
    plot(x,y,'k')
    hold on
    
    
    %replace NaNs with coordinate for each vertex to close the polygon
    x_big = repmat(x(1:end-1),length(theta),1);
    xlos(isnan(xlos)) = x_big(isnan(xlos));
    y_big = repmat(y(1:end-1),length(theta),1);
    ylos(isnan(ylos)) = y_big(isnan(ylos));
    
    
    
    %     for iv = 1:nvert
    %         figure()
    %         plot(x,y)
    %         hold on
    %         fill(xlos(:,iv), ylos(:,iv),'b');
    %         hold on
    %     end
    FetchArea = cell(size(shoreline_fetch));
    for iv = 1:nvert
        FetchArea{obj,1}(:,iv)=polyarea(xlos(:,iv),ylos(:,iv));
    end
end
% title(sprintf('theta = %.2f from each vertex', theta(50)));

% end