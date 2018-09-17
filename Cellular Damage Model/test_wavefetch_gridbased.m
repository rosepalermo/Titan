% This is a test of doing a cellular version of the fetch calc code
% starting with just one vertex
clear all
addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Cellular Damage Model/old function attempts')
load('test_polygon')

[edge] = addidshoreline_cardonly(polygon,~polygon);
[x,y] = find(edge); % for test, x and y are indices of shoreline

% nvert = length(x) - 1; % # vertices in polygon, assuming closed
nvert=1;
% Set up rays
nray = 1000;
theta = linspace(0,2*pi,nray+1);
theta = theta(1:end-1); % eliminate duplicate at 2pi
D = 2*max(sqrt(x.^2 + y.^2));
eps = D*1e-8;
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
%     xseg = [ones(size(theta))*x(iv)+dx_eps; xend; nan(size(theta))];
%     yseg = [ones(size(theta))*y(iv)+dy_eps; yend; nan(size(theta))];
    xseg = [ones(size(theta))*x(iv)+dx_eps; xend; zeros(size(theta))]';
    yseg = [ones(size(theta))*y(iv)+dy_eps; yend; zeros(size(theta))]';
    v1 = [xseg(:,1) yseg(:,1) xseg(:,3)];
    v2 = [xseg(:,2) yseg(:,2) yseg(:,3)];
    pt = [x y zeros(size(x))];
    % now have the line segments, the originating point, and all shoreline
    % points...
    d = zeros(length(pt),length(theta));
    for ii = 1:length(theta)
        for iii = 1:length(pt)
            a = v1(ii,:) - v2(ii,:);
            b = pt(iii) - v2(ii,:);
            d(iii,ii) = norm(cross(a,b)) / norm(a); % distance from points on shoreline to line segments
        end
    end
    
end

