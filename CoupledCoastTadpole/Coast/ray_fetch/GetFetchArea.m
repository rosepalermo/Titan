function [warea,farea,wpiws,wpjws,rayazimuths,coastazimuths,fpis,fpjs] = GetFetchArea(shorelines,lake,nrays,delta,cellsize)

ns = length(shorelines);

% Convert to matrix units 
for n = 1:ns
    shorelines{n}(:,1) = shorelines{n}(:,1)/cellsize;
    shorelines{n}(:,2) = shorelines{n}(:,2)/cellsize;                
    
end

% Detect islands.
% We reverse the order of island points so that islands are ordered CCW if i increases up)
island = zeros(ns,1);
for n = 1:ns

    sj = shorelines{n}(:,1);
    si = shorelines{n}(:,2);
    
    % Loop through other shorelines in the grid. If any point of shoreline
    % n is inside any other shoreline polygon, shoreline n is an island.
    for nn = 1:ns
        if nn ~= n % skip shoreline n
            pj = shorelines{nn}(:,1);
            pi = shorelines{nn}(:,2);
            
            if inpolygon(sj(1),si(1),pj,pi)
                island(n) = 1;
            end
        end
    end

%     if island(n)
%         shorelines{n}(:,1) = flipud(sj);
%         shorelines{n}(:,2) = flipud(si);        
%     end
    
end

% lake is 1 for water, 0 for land
% shorelines is a cell array, one element for each closed shoreline or island.
% each cell contains x and y coordinates (column and row indices) ordered clockwise if i increases downward, first point != last point
% That ordering (CW if i increases down, CCW if i increases up) is the same
% for lakes and islands.

ns = length(shorelines);
warea = cell(ns,1);
farea = cell(ns,1);
wpiws = cell(ns,1);
wpjws = cell(ns,1);
fpis = cell(ns,1);
fpjs = cell(ns,1);
rayazimuths = cell(ns,1);
coastazimuths = cell(ns,1);

% What we need to do when we encounter land is (1) check if it's on a
% shoreline (not necessarily the one we're currently working on), and (2)
% find the previous and next points on THAT shoreline. 

% Create one matrix that contains the linear indices into vectors that
% contain the previous and next shoreline points. 

sz = size(lake);

Sidx = zeros(sz); % this will contain the 1-based indices of each shoreline point in Sprev and Snext
siprev = []; % vector of 1-based row indices of previous points in the lake matrix
sinext = []; % vector of 1-based row indices of next points in the lake matrix
sjprev = []; % vector of 1-based column indices of previous points in the lake matrix
sjnext = []; % vector of 1-based column indices of next points in the lake matrix

c = 1; % This is the index where we start adding the next shoreline in Sprev and Snext

% first, loop through the shorelines to build these arrays to pass to
% rayfetch.c
for n = 1:ns

%     % Convert from x,y values to matrix indices
%     si = shorelines{n}(:,2)/cellsize; % row (y) indices
%     sj = shorelines{n}(:,1)/cellsize; % column (x) indices

    si = shorelines{n}(:,2); % row (y) indices
    sj = shorelines{n}(:,1); % column (x) indices
       
    sind = int32(sub2ind(sz,si,sj)); % linear indices
    
    npts = length(sind); % Number of points in this shoreline segment

    Sidx(sind) = c:c+npts-1;
    
    siprev(c:c+npts-1) = [si(end); si(1:end-1)]; % row indices of previous shoreline points
    sinext(c:c+npts-1) = [si(2:end); si(1)]; % col indices of next shoreline points
    sjprev(c:c+npts-1) = [sj(end); sj(1:end-1)]; % row indices of previous shoreline points
    sjnext(c:c+npts-1) = [sj(2:end); sj(1)]; % col indices of next shoreline points
        
    c = c + npts;

end

siprev = siprev(:);
sinext = sinext(:);
sjprev = sjprev(:);
sjnext = sjnext(:);

% Now loop through the shorelines to compute fetch area 
for n = 1:ns

%     si = shorelines{n}(:,2)/cellsize; % row (y) indices
%     sj = shorelines{n}(:,1)/cellsize; % column (x) indices

    si = shorelines{n}(:,2); % row (y) indices
    sj = shorelines{n}(:,1); % column (x) indices

    [fpi,fpj,rayaz,coastaz] = rayfetch(lake,si(:),sj(:),nrays,delta,Sidx,siprev,sinext,sjprev,sjnext);
    
    
    % convert from matrix units back to model units
    si = shorelines{n}(:,2)*cellsize; % row (y) indices
    sj = shorelines{n}(:,1)*cellsize; % column (x) indices
    fpi = fpi*cellsize;
    fpj = fpj*cellsize;
    
    
    % compute weighted fetch polygon areas
    npts = length(si);
    weights = cos( abs(repmat(coastaz,[1 nrays]) - repmat(rayaz,[npts 1])) ); % cosine of the angle between the ray direction and the coast normal
    Li = fpi - si(:); % i distances from shoreline points to ends of fetch rays
    Lj = fpj - sj(:); % j distances
    wpiw = si + weights.*Li; % to get endpoints of weighted fetch rays, add weighted distances
    wpjw = sj + weights.*Lj;
    warea{n} = polyarea(wpjw,wpiw,2);
    fpiw = si + Li; % get endpoints of fetch rays
    fpjw = sj + Lj;
    farea{n} = polyarea(fpjw,fpiw,2);
    
    wpiws{n} = wpiw;
    wpjws{n} = wpjw;
    fpiws{n} = fpiw;
    fpjws{n} = fpjw;
    fpis{n} = fpi;
    fpjs{n} = fpj;    
    rayazimuths{n} = rayaz;
    coastazimuths{n} = coastaz;

end



% return island and associated output order back to original order
% We reverse the order of island points so that islands are ordered CCW if i increases up)
for n = 1:ns

    sj = shorelines{n}(:,1);
    si = shorelines{n}(:,2);

    if island(n)
        shorelines{n}(:,1) = flipud(sj);
        shorelines{n}(:,2) = flipud(si);
        warea{n} = flipud(warea{n});
        farea{n} = flipud(farea{n});
        wpiws{n} = flipud(wpiws{n});
        wpjws{n} = flipud(wpjws{n});
        rayazimuths{n} = flipud(rayazimuths{n});
        coastazimuths{n} = flipud(coastazimuths{n});
        fpis{n} = flipud(fpis{n});
        fpjs{n} = flipud(fpjs{n});
    end
    
end