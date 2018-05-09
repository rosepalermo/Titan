% leif.m

% load river centerline coordinates
filename = 'lg_all_pts.xls';
% imagename = 'jingpo.tif';
% [A,R]=gceotiffread(imagename);
savename = 'trash.csv';

M = xlsread(filename);
% M = csvread(filename);


%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

xx = M(:,4);
yy = M(:,5);
% xx=x;yy=y;
%reduce to every 5th element
% xx=x(1:5:end);
% yy=y(1:5:end);
%add the first point to the end to make periodic
% x=ones(length(xx)+1,1);
% y=ones(length(yy)+1,1);
% x(1:end-1) = xx(:,1);
% y(1:end-1) = yy(:,1);
% x(end)=x(1);
% y(end)=y(1);

x = xx;
y = yy;

% total_length = arclength(x,y,'linear')
% [xy]=interparc(0:(200/total_length):1,x,y,'spline');
% x=xy(:,1);y=xy(:,2);
% xx=ones(length(x)+1,1);
% yy=ones(length(y)+1,1);
% xx(1:end-1) = y(1:end);
% yy(1:end-1) = x(1:end);
% xx(end)=xx(1);
% yy(end)=yy(1);



% make the first point (0,0)
% x = x-x(1);
% y = y-y(1);

% transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% streamwise distance
[theta, dtheta, deltad] = meander(x,y);
xx=x;
yy=y;

% do the wavelet transforms, comparing the spectra to those of an 
% autoregressive process of order n (a.k.a. an AR(n) process)
n=1;
dowave(theta,deltad,n,xx,yy,savename);
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
axis square

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

