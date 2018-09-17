clear all
addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Cellular Damage Model/old function attempts')
load('test_polygon')

land = ~polygon;
[edge] = addidshoreline_cardonly(polygon,~polygon);
[x,y] = find(edge); % for test, x and y are indices of shoreline

% Rose you need to make one long list of all shoreline points on the lake
% and islands for the future

vx = ones(size(x))*x(1); vy = ones(size(y))*y(1);
% pcolor(edge)
[X,Y] = meshgrid(min(x):max(x), min(y):max(y));
plot(X,Y,'b.')
hold on

shoreline = [x,y];
plot(x,y,'ro')
plot_rays(shoreline(1,:),shoreline)

bsxfun(@plus,shoreline,-shoreline(1,:))
