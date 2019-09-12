% fetch.m

clear

addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Tadpole 2')

p.Nx = 400;                 %     p.Nx             Number of grid points in x direction
p.Ny = 400;                 %     p.Ny             Number of grid points in y direction
p.dx = 125/2;                 %     p.dx             Grid spacing in the x direction (m)
p.dy = 125/2;                 %     p.dy             Grid spacing in the y direction (m)

p.beta = 1.6;               %     p.beta           Negative slope of the power spectrum. 0 = white noise, more positive values are "redder" (more variance at longer wavelengths)
p.variance = 10000;         %     p.variance       Variance of elevation (m^2)
p.periodic = 1;             %     p.periodic       Elevations will be periodic at the boundaries (1) or not (0, default)


% land = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);
% 
% % make lowest 10% of elevations 1, otherwise 0
% SL = prctile(land(:),30);
% land = land - SL;
% land(land >= 0) = 0; 
% land(land < 0) = 1; 
% 
% % Make a figure to make sure the liquid doesn't span a boundary.
% figure; imagesc(land);

% volcanic island shield
% xc = (p.Nx-1)/2+1; % center of island
% yc = (p.Ny-1)/2+1;
% r = 0.8*(p.Nx-1)/2; % radius of island
% S = 0.05; % slope of shield
% [X Y] = meshgrid(1:p.Nx,1:p.Ny);
% D = sqrt((X-xc).^2 + (Y-yc).^2);
% shield = (r - D)*p.dx*S; % the initial shield profile
initgaus = get_gaussian_boundary([800 800], 0.1, sqrt(p.variance));

% add a random component
noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);
noise = noise - mean(noise(:));

figure()
imagesc(noise)
title('noise')

% sea = (initgaus + noise);
[H Wss] = Hann2D(noise);
sea = H;
figure()
imagesc(sea)
title('H + noise')

SL = prctile(sea(:),10);
sealine = sea;
sealine(sealine >= SL) = 0; 
sealine(sealine < SL) = 1; 

figure()
imagesc(sealine)
title('sea')


% Find points on the coast (1) or not (0)
coast = bwmorph(sealine,'remove');

close all

figure()
imagesc(coast)

% get a list of locations of all the points on the coast. We will do
% everything in matrix indices (y = row, x = col)
[yc,xc] = find(coast);

% For each point on the coast

% nc = length(yc);
% 
% for i=1:nc
%     % the segments joining point i with each other point on the coast are
%     % the candidate rays. Calculate the slope (as an angle CCW from E) from
%     % point i to each of the other coast points. Now we have the starting
%     % point, end point, and azimuth of each segment.
%     
%     % Next we find the indices of the cells that each segment passes
%     % through (excluding the cells corresponding to the 2 endpoints).
%     
%     % Find the x coordinates of the cols containing the segment (except the
%     % starting point and endpoint). 
%     
%     % Shift those x coords half a cell towards the starting point of the
%     % segment (use signs of unit vector components or direction cosines).
%     
%     % Calculate the y coords on the segment corresponding to the shifted 
%     % coordinates using the equation for the segment (starting point and
%     % slope). Round to the nearest integer to get the row coords. Shift the
%     % x coords back a half cell (away from the starting point) to their 
%     % original values. (Could do the rounding and shifting in same step, so
%     % no need to modify col coords really.)
% end
