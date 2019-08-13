function [WaveArea,FetchArea] = fetch_vis(fetch_sl_cells)

% variables for the vis code
epsilon = 1e-6;
snap_distance = 0;

% this is bc of ccw shoreline..
for i = 1:length(fetch_sl_cells)
    fetch_sl_cells{i} = flipud(fetch_sl_cells{i});
    environment{i} = [fetch_sl_cells{i}(:,:);fetch_sl_cells{i}(1,:)];
end

% % "environment" is the shoreline
% environment = fetch_sl_cells;

for ii = 1:length(fetch_sl_cells)
    %shoreline
    x = fetch_sl_cells{ii}(:,1);
    y = fetch_sl_cells{ii}(:,2);
    
    
    % find a point eps inside the shoreline, because it won't let me use a
    % point on the boundary of the environment to calculate.
    x_wrapped = [x; x(1,:)];
    y_wrapped = [y; y(1,:)];
    dif = [x_wrapped(2:end), y_wrapped(2:end)] - [x, y];
    perp = [dif(:,2), -dif(:,1)];
    normal = sqrt(dif(:,1).^2 + dif(:,2).^2);
    eps = .01;
    perp = eps*perp./normal;
    x_inside = x + perp(:,1);
    y_inside = y + perp(:,2);
    
    % find fetch polygon for each point
    for l = 1:length(x_inside)
        V{ii,l} = visibility_polygon( [x_inside(l) y_inside(l)] , environment , epsilon , snap_distance );
        
        % fetch & wave area & distance & cos(theta-phi)
        FetchArea{ii}(l,1) = polyarea(V{l}(:,1),V{l}(:,2));
        Fetch_dist{ii}(l,1) = sqrt(sum(([x(l),y(l)] - V{l}).^2,2));
        % Wave weighting = (F)*cosang
        weighted_fetch_dist{ii}(l,1) = ([x(l),y(l)] - V{l})*perp(l,:)'/eps;
        % cos(theta - phi) = dot product of slvec and losvec
        cosang{ii}(l,1) = weighted_fetch_dist{l}./Fetch_dist{l};
        cosang{ii}(isnan(cosang{l})) = 0;
        Wavepts{ii}(l,1) = [x(l),y(l)]+(V{l}-[x(l),y(l)]).*cosang{l};
        WaveArea{ii}(l,1) = polyarea(Wavepts{l}(:,1),Wavepts{l}(:,2));
    end
end
end