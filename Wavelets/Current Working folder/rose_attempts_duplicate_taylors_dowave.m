clear

% load river centerline coordinates
filename = 'lg_all_pts.xls';
M = xlsread(filename);
savename = 'trash.csv'; % savename is commented out in dowave, so not saving anything


%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

x0=M(:,4);
y0=M(:,5);

% CONSIDER NOT INTERPOLATING HERE. WE DON'T NEED EVENLY SPACED POINTS TO
% COMPUTE AZIMUTHS, AND IT CAN CREATE STEPS IN THE AZIMUTH VS. DISTANCE
% SERIES IF INTERPOLATING CREATES MULTIPLE SUCCESSIVE POINTS WITH THE SAME
% AZIMUTH.

% % interpolate
% dist=sqrt((x0(2:end)-x0(1:end-1)).^2+(y0(2:end)-y0(1:end-1)).^2); % distance increment
% dist_cum = [0;cumsum(dist)];
% dist_total = sum(dist);
% dist_interp = linspace(0,dist_total,length(x0));
% x_int = interp1(dist_cum,x0,dist_interp);
% y_int = interp1(dist_cum,y0,dist_interp);
% figure()
% scatter(x0,y0,'r')
hold on
% scatter(x_int,y_int,'k')
% 
% x0 = x_int';
% y0 = y_int';


% rearrange ligeia for debugging
x0 = [x0(23186:end);x0(1:23185)];
y0 = [y0(23186:end);y0(1:23185)];



% add first AND SECOND points to the end for meander
x = [x0;x0(1:2)];
y = [y0;y0(1:2)];



% transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% streamwise distance
[theta, dtheta, deltad] = meander_titan(x,y);
x = x(1:end-2);
y = y(1:end-2);


% do the wavelet transforms, comparing the spectra to those of an 
% autoregressive process of order n (a.k.a. an AR(n) process)
n=2;
% dowave_duplicate(theta,deltad,n,x,y,savename);
dowave(theta,deltad,n,x,y,savename);

% %plot where the first point is
% figure
% % plot(x,y,'k','LineWidth',2)
% hold on
% scatter(x,y,'.','k')
% scatter(x(1),y(1),'*','r')
% scatter(x(30),y(30),'*','b')
% xlabel('X')
% ylabel('Y')
% axis tight
% axis equal
% 
