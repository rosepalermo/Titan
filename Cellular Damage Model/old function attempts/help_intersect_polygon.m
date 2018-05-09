% polygon. in this example, a circle
theta = 0 : 0.1 : 2*pi;
radius = 2;
x = radius * cos(theta);
y = radius * sin(theta);
% add the first point to the end to close the polygon
x = [x x(1)];
y = [y y(1)];
plot(x,y,'k','LineWidth',2)
hold on

% find 1/2 the smallest distance from one point to the next
eps = 0.5*sqrt(min((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2));

% find 2 * the maximum distance between 2 points
D = 2*max(sqrt(x.^2 + y.^2));

% angles around from point
theta = linspace(0,2*pi,1000);

for j=1:length(x)
    dx = D*cos(theta);
    dy = D*sin(theta);
    % check if a small ray is in or on the polygon
    dx_eps = eps * cos(theta);
    dy_eps = eps * sin(theta);
    
     %inside of the bounding polygon ( in this case a circle)
    [in, on] = inpolygon(x(j)+dx_eps,y(j)+dy_eps,x,y);
    ind_inon = find(in|on);
    
    % for the points that are in or on the polygon, find the intersection
    % of a line segment in that direction with the polygon
     for i=2:length(ind_inon)
         
         % line segment to find intersections for each angle from point j
        x1 = [x(j), x(j) + dx(ind_inon(i))];
        y1 = [y(j), y(j) + dy(ind_inon(i))];
        
        [xi,yi] = polyxpoly(x1,y1,x,y); % find intersection between ray and each polygon
        
        % find shortest intersection that isn't 0
                [m,k] = min((x(j)-xi).^2 + (y(j)-yi).^2);
        while m < eps^2
            xi(k) = [];
            yi(k) = [];
            [m,k] = min((x(j)-xi).^2 + (y(j)-yi).^2);
        end
         %if there is no intersection (lines are parallel and on top of each
        %other), then use point j because I want a filler in there
        if isempty(xi)
            x_int_save(ind_inon(i)) = x(j);
            y_int_save(ind_inon(i)) = y(j);
            continue
        end
        % save intersecting points
                % point int = points that make up fetch area polygon
        x_int = xi(k);
        y_int = yi(k);
        x_int_save(ind_inon(i)) = x_int;
        y_int_save(ind_inon(i)) = y_int;
        plot([x(j) x_int],[y(j) y_int]); drawnow   % plot all intersecting rays
     end
end

    