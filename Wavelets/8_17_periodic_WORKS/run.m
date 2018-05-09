% leif.m

% load river centerline coordinates
filename = 'lg_all_pts.xls';
% imagename = 'jingpo.tif';
% [A,R]=gceotiffread(imagename);
savename = 'trash.csv';

M = xlsread(filename);
%M = csvread(filename);


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

% % interpolate
% dist=sqrt((x0(2:end)-x0(1:end-1)).^2+(y0(2:end)-y0(1:end-1)).^2); % distance increment
% dist_cum = [0;cumsum(dist)];
% dist_total = sum(dist);
% dist_interp = linspace(0,dist_total,length(x0));
% xtest = interp1(dist_cum,x0,dist_interp);
% ytest = interp1(dist_cum,y0,dist_interp);
% figure()
% scatter(x0,y0,'r')
% hold on
% scatter(xtest,ytest,'k')
% 
% x0 = xtest';
% y0 = ytest';

% make rearrange ligeia for debugging
% x0 = [x0(23186:end);x0(1:23185)];
% y0 = [y0(23186:end);y0(1:23185)];


% sebago island data
% xx = [x0(3319+25:3459);x0(2223:2787);x0(3460:4155);x0(5442:7123);x0(87:1664);x0(1941:2222);x0(4156:5044);x0(8329:8703);x0(7124:8328);x0(5381:5440);x0(5045:5379)];
% yy = [y0(3319+25:3459);y0(2223:2787);y0(3460:4155);y0(5442:7123);y0(87:1664);y0(1941:2222);y0(4156:5044);y0(8329:8703);y0(7124:8328);y0(5381:5440);y0(5045:5379)];
% x_island = [x0(1665:1940);x0(2788:3318)];
% y_island = [y0(1665:1940);y0(2788:3318)];
xx = x0;
yy = y0;
% x=xx;y=yy;

%reduce to every other element
% xx=xx(1:2:end);
% yy=yy(1:2:end);
% plot(xx,yy)
% hold on


%add the first point to the end to close--make sure to take off
%last point in dowave if you're leaving this.


% when all of the lake
x = [xx;xx(1)];
y = [yy;yy(1)];

%when only part of the lake
% x = xx; y = yy;





% transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% streamwise distance
[theta, dtheta, deltad] = meander_titan(x,y);
x = x(1:end-1);
y = y(1:end-1);


% do the wavelet transforms, comparing the spectra to those of an 
% autoregressive process of order n (a.k.a. an AR(n) process)
n=2;
dowave_duplicate(theta,deltad,n,x,y,savename);
%dowave(dtheta,deltad,n,xx,yy);

%plot where the first point is
figure
% plot(x,y,'k','LineWidth',2)
hold on
scatter(x,y,'.','k')
scatter(x(1),y(1),'*','r')
scatter(x(30),y(30),'*','b')
xlabel('X')
ylabel('Y')
axis tight
axis equal

% %plot in thirds
% x1=x(1:length(x)/3);
% x2=x(length(x1)+1:length(x)./3*2);
% x3=x(2*length(x2)+1:length(x));
% 
% y1=y(1:length(y)/3);
% y2=y(length(y1)+1:length(y)./3*2);
% y3=y(2*length(y2)+1:length(y));
% 
% 
% figure
% plot(x,y,'k','LineWidth',2)
% hold on
% scatter(x(1),y(1),'*','r')
% plot(x1,y1,'r','LineWidth',2)
% plot(x2,y2,'g','LineWidth',2)
% plot(x3,y3,'b','LineWidth',2)
% scatter(x(1),y(1),'x','r','LineWidth',2)
% %scatter(x(100),y(100),'o')
% xlabel('X')
% ylabel('Y')
% axis tight
% axis equal

