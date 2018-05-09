x = x'; y = y'; %make the columns of x and y rows of x and y
x = [x x(1)]; % close it
y = [y y(1)];
nvert = length(x) - 1; % # vertices in polygon, assuming closed
% Set up rays
nray = 1000;
theta = linspace(0,2*pi,nray+1);
theta = theta(1:end-1); % eliminate duplicate at 2pi
D = 2*max(sqrt(x.^2 + y.^2));
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
      [xi,yi] = polyxpoly(xseg(:), yseg(:), x, y);
      
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
      
      % Now figure out whether each ray was inside or outside the polygon
      % before intersecting.  Can check simply by seeing if bisector is
      % inside or outside.
      dx = ximin - x(iv);
      dy = yimin - y(iv);
      isin = inpolygon(x(iv) + dx/2, y(iv) + dy/2, x, y);
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
% plot(x,y,'k');
% hold on;
% for ivv = 1:nray
%     for iv = 1:nvert
%         plot(xplt(:,ivv,iv), yplt(:,ivv,iv));
%         hold on
%     end
% end
% 
% figure()
% plot(x,y,'k')
% hold on


%replace NaNs with coordinate for each vertex to close the polygon
x_big = repmat(x(1:end-1),length(theta),1);
xlos(isnan(xlos)) = x_big(isnan(xlos));
y_big = repmat(y(1:end-1),length(theta),1);
ylos(isnan(ylos)) = y_big(isnan(ylos));



    for iv = 1:nvert
        figure()
        plot(x,y)
        hold on
        fill(xlos(:,iv), ylos(:,iv),'b');
        hold on
    end

for iv = 1:nvert
    FetchArea(iv)=polyarea(xlos(:,iv),ylos(:,iv));
end
% title(sprintf('theta = %.2f from each vertex', theta(50)));