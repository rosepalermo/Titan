% rayfetch_test.m

% The code now does the following:
% Outputs ray azimuths and coast normal azimuths.
% Skips rays that would correspond to waves that "come from behind".
% Shows how to compute weighted fetch area using the cosine.
% Allows for islands, both as originating points and shoreline
% intersections.
%
% To do:
% Consider finding shoreline intersection points instead of just stopping before intersecting shoreline

%% Main calculation

% load a test lake and the ordered shoreline
load('fetch_input8con_lake','F_lake'); % lake is 1 for water, 0 for land
load('fetch_input8con','fetch_sl_cells'); % x and y are ordered clockwise if i increases downward, first point != last point

shorelines = fetch_sl_cells;
% shorelines{1}(:,1) = flipud(shorelines{1}(:,1));
% shorelines{1}(:,2) = flipud(shorelines{1}(:,2));
lake = F_lake;

cellsize = 62.5; % in units of x and y
nrays = 180; % 360/5number of rays around the circle. Recommend nrays >= 36
delta = 0.05; % 0.5 distance increment along each ray, in cells. Recommend delta <= 0.5

ns = length(shorelines);

tic

% % Convert to matrix units 
% for n = 1:ns
% 
%     shorelines{n}(:,1) = shorelines{n}(:,1)/cellsize;
%     shorelines{n}(:,2) = shorelines{n}(:,2)/cellsize;                
%     
% end
% 
% 
% % % Detect islands.
% % % We reverse the order of island points so that islands are ordered CCW if i increases up)
% island = zeros(ns,1);
% for n = 1:ns
% 
%     sj = shorelines{n}(:,1);
%     si = shorelines{n}(:,2);
%     
%     % Loop through other shorelines in the grid. If any point of shoreline
%     % n is inside any other shoreline polygon, shoreline n is an island.
%     for nn = 1:ns
%         if nn ~= n % skip shoreline n
%             pj = shorelines{nn}(:,1);
%             pi = shorelines{nn}(:,2);
%             
%             if inpolygon(sj(1),si(1),pj,pi)
%                 island(n) = 1;
%             end
%         end
%     end
% 
%     % Convert from x,y units to matrix indices
%     if island(n)
%         shorelines{n}(:,1) = flipud(sj);
%         shorelines{n}(:,2) = flipud(si);        
%     end
%     
% end

[warea,~,wpiws,wpjws,rayaz,coastaz,fpis,fpjs] = GetFetchArea(shorelines,lake,nrays,delta,cellsize);

toc

% for n = 1:ns
% 
%     shorelines{n}(:,1) = shorelines{n}(:,1)*cellsize;
%     shorelines{n}(:,2) = shorelines{n}(:,2)*cellsize;                
%     
% end
%% plot weighted fetch area

figure
ns = length(shorelines);
for n=1:ns
    si = shorelines{n}(:,2); % row indices
    sj = shorelines{n}(:,1); % column indices
    plot(sj,si,'-k')
    hold on
    scatter(sj,si,16,log10(warea{n}),'filled'); axis equal
end
set(gca,'ydir','reverse')
colorbar
title('log10(weighted fetch area in cells)')


%% plot rays for selected points

figure
% imagesc(~lake); axis equal tight; colormap gray
hold on

ns = length(shorelines);

% plot the shorelines, including islands
for s=1:ns
    si = shorelines{s}(:,2); % row indices
    sj = shorelines{s}(:,1); % column indices
    plot(sj,si,'-k')
end

% now plot fetch rays, one point at a time
for s=1:ns
    si = shorelines{s}(:,2); % row indices
    sj = shorelines{s}(:,1); % column indices

    for n = 1:length(si)
        h = zeros(1,nrays);
        hw = zeros(1,nrays);

        
        f = plot(sj(n),si(n),'ok','markersize',16,'markerfacecolor','r');
        hc = plot([sj(n) sj(n)+cos(coastaz{s}(n))],[si(n) si(n)+sin(coastaz{s}(n))],'b');
        for k = 1:nrays
            h(k) = plot([sj(n) fpjs{s}(n,k)],[si(n) fpis{s}(n,k)],'-r'); % unweighted rays
            hw(k) = plot([sj(n) wpjws{s}(n,k)],[si(n) wpiws{s}(n,k)],'-g'); % weighted rays
        end
        pause
        delete(f)
        delete(h)
        delete(hw)
        delete(hc)
    end
end
