% taylor_order_shoreline.m

load('wavet2v2_lake.mat', 'lake')

land = double(~lake);
tic
[B,L,N,A] = bwboundaries(land,4);

% As long as there isn't liquid on the boundary of the grid, the first
% object should be the main landmass. The holes of that object will be the
% lakes/seas, and the object children of the main object will be islands --
% and those will have already had their coasts traced!

k = 1; % the index of the object that is the main landmass; we assume it is 1

lakes = find(A(:,k))'; % the indices of the holes that are lakes/seas

% We know that each of the points on the boundary of each lake has at least
% one neighbor that is land. So all we have to do is find one of those land
% points, and then we can use bwtraceboundary to trace the 4-connected land
% coastline.

nextrow = [1,0,0; 1,0,-1; 0,0,-1];
nextcol = [0,-1,-1; 0,1,0; 1,1,0];

coasts = cell(size(B)); % this will hold the ordered (CW) coastlines of the lakes/seas

for l = lakes 
    bdy = B{l};
    r1 = bdy(1,1);
    c1 = bdy(1,2);
    r2 = bdy(2,1);
    c2 = bdy(2,2);
    
    % We know that point 2 is clockwise from point 1 on the 8-connected
    % edge of the "water". So we go CCW one increment (45 degrees) in the
    % 8-connected neighborhood of point 1, and that's our land point.
    dr = r2-r1; % the difference in row index
    dc = c2-c1; % the difference in column index (note that columns increase down!)
    ridx = 2+dr;
    cidx = 2+dc;
    rowCCW = r1 + nextrow(ridx,cidx);
    colCCW = c1 + nextcol(ridx,cidx);
    
    coasts{l} = bwtraceboundary(land,[rowCCW,colCCW],'S',4,Inf,'counterclockwise');
end

% Now we add in the boundaries of the islands, which we already have from
% bwboundaries
islands = 2:N; % the indices of the objects that are islands
coasts(islands) = B(islands); % The non-empty elements of the cell array coasts 
                              % should now have CW-ordered coasts for the land, not the liquid.
toc
% let's see what we got
figure
imagesc(L); axis image
hold on
for c = [islands,lakes]
    coast = coasts{c};
    plot(coast(:,2),coast(:,1),'w','LineWidth',2);
end
