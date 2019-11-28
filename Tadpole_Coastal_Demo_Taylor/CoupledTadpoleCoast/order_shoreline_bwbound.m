function [sl_cell,keepme,cells2trash] = order_shoreline_bwbound(lake);

% taylor_order_shoreline.m
% Rose updated to couple with coastal damage models
%


% This code finds the ordered shoreline points (8-connected) for all lakes
% and islands using bwboundaries and bwtraceboundary. It also identifies
% points that need to be deleted because they have less than 3 connected
% shoreline points in a row, so small islands and lakes, because it messes
% up the fetch code.


% INPUT:        lake            logical array. liquid is 1. land is 0.
% OUTPUT:       sl_cell         cell array of ALL ordered boundaries
%               keepme          indices of sl_cell to keep (larger than 3
%                               in length)
%               cells2trash     shoreline points that are less than 3 in
%                               length. sl_cells(~keep me).


[~,total_lakes,total_islands] = find_first_order_lakes(lake);
land = double(~lake);
[B_land,L_land,~,~] = bwboundaries(land,8);
if length(B_land) ==1 % this occurs when shoreline hits the boundary
    sl_cells = [];
    return
end

[B,L,N,A] = bwboundaries(lake,8);

% As long as there isn't liquid on the boundary of the grid, the first
% object should be the main landmass. The holes of that object will be the
% lakes/seas, and the object children of the main object will be islands --
% and those will have already had their coasts traced!

% Assume the boundary is the first element of B and get rid of it.
% BB = B(2:end);

k = 1; % the index of the object that is the main landmass; we assume it is 1

% lakes = find(A(:,k))'; % the indices of the holes that are lakes/seas
lakes = zeros(N,1);
for i=1:N
    if ~any(A(i,:))
        lakes(i) = i;
    end
end
lakestemp = zeros(N,1); % or N?
lakestemp(lakes) = lakes;
lakes = lakestemp;
if sum(total_islands)>0
    islands = zeros(length(lakes),1);
    islands_temp = zeros(length(lakes),1);
    for i=find(lakes)'
        islands_temp = zeros(size(lakes)) | A(:,i);
        islands(islands_temp) = i;
        clearvars islands_temp
        %     first_order_islands = first_order_islands | A(:,i);
        %     first_order_islands(:,2) = first_order_islands(:,1)*i;
    end
end

% We know that each of the points on the boundary of each lake has at least
% one neighbor that is land. So all we have to do is find one of those land
% points, and then we can use bwtraceboundary to trace the 4-connected land
% coastline.

nextrow = [1,0,0; 1,0,-1; 0,0,-1];
nextcol = [0,-1,-1; 0,1,0; 1,1,0];

% size(B,1)-1 because it includes the boundary of the domain as a boundary
coasts = cell(1,size(B,1)-1); % this will hold the ordered (CW) coastlines of the lakes/seas
for l = find(lakes)
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
    
    coasts{1} = (bwtraceboundary(land,[rowCCW,colCCW],'S',4,Inf,'counterclockwise'));
    coasts{1}(end,:) = [];
end

% Now we add in the boundaries of the islands, which we already have from
% bwboundaries
% islands = 3:N; % the indices of the objects that are islands
if sum(total_islands)>0% if there are islands
    if sum(islands)>0 % if there are islands
    coasts(1,find(islands)) = flipud(B(find(islands))); % The non-empty elements of the cell array coasts
    for island_idx=find(islands)
        coasts{island_idx}(end, :) = []; % remove the last point of each island, duplicate of first point
%         if any(ismember(coasts{island_idx},coasts{1},'rows'))
%             coasts{island_idx} = []; % remove the island that overlaps with a promontory
%         end
    end
    % should now have CW-ordered coasts for the land, not the liquid.
    end
end

% Rose to check if units 3 or less still mess up in the fetch calc.
% get rid of units 3 or less because they mess up in the fetch calc.
length_cells = cellfun(@length,coasts,'uni',false)';
length_cells = cell2mat(length_cells);
keepme = (length_cells) > 3;
cells2trash = cell2mat(coasts(~keepme)');
sl_cell = coasts(keepme);

% test = 1;

% % PLOT
% figure
% imagesc(L); axis image
% hold on
% for c = [lakes]
%     coast = coasts{c};
%     scatter(coast(:,2),coast(:,1),'w');
% end
% 
% %plot the cells we're going to keep in cyan
% for c=1:length(sl_cell)
%     slplot = sl_cell{c};
%     plot(slplot(:,2),slplot(:,1),'c','LineWidth',2);
% end
% %plot the cells we're going to trash in red
% scatter(cells2trash(:,2),cells2trash(:,1),'r','LineWidth',2);

end
