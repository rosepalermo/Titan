clear all
load('test_polygon')

% find the first point for the trace
first_ind = find(edge);
[i,j] = ind2sub(size(polygon),first_ind(1));

edge = 2*double(edge); % times 2 so that you can tell the difference in the figure
[edgeXind,edgeYind] = find(edge);
%add the edge of the polygon that I'm trying to find in order to the polygon
polygon_temp = (polygon+edge);

%trace the whole thing -- problem:lose the details where the polygon has tight angles
edge_bwxy = bwtraceboundary(polygon_temp,[i,j],'NE');

figure()
imagesc(polygon')
title('original polygon')

figure()
imagesc(polygon_temp')
hold on
plot(edgeXind,edgeYind,'rx')
plot(edge_bwxy(:,1),edge_bwxy(:,2),'ro')
title('polygon + edge')

figure()
imagesc(edge')
title('edge, note missing vertices from figure 2')


