function [binfine,bincoarse] = SimplifyShoreline(x,y,tolfine,tolcoarse)

% This version makes two approximations of the shoreline using the
% Douglas-Peucker algorithm: a fine approximation using a tolerance tolfine
% (expressed as a fraction of the max range of x and y) and a coarse
% approximation using tolcoarse. Then, when we're looping through the
% vertices, we will choose which points to draw from the fine approximation
% (within a certain distance of the vertex) and which points to draw from
% the coarse approximation (beyond that distance). The approximations are
% given binary arrays indicating the indices of the points in x and y to retain.

np = length(x);

tolfine = tolfine*max([range(x) range(y)]);
tolcoarse = tolcoarse*max([range(x) range(y)]);

[~,ifine] = simplifyPolygon([x,y], tolfine);
[~,icoarse] = simplifyPolygon([x,y], tolcoarse);

binfine = zeros(1,np);
bincoarse = zeros(1,np);
binfine(ifine) = 1;
bincoarse(icoarse) = 1;
