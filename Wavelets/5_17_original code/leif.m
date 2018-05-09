% leif.m

% load river centerline coordinates
filename = 'lake_at_5_76_automated_xy.csv';
% imagename = 'jingpo.tif';
% [A,R]=gceotiffread(imagename);
savename = 'lake_at_5_76_automated_xy_specslope.csv';

%M = xlsread(filename);
M = csvread(filename);


%get rid of duplicate points (when it goes exactly around a pixel)
% indices to unique values in column 3
[~, ind] = unique(M(:, 4:5), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 4) M(duplicate_ind, 5)];
M(duplicate_ind,:)=[];

%add the first point to the end to make periodic
x=ones(length(M(:,4))+1,1);
y=ones(length(M(:,5))+1,1);
x(1:end-1) = M(:,4);
y(1:end-1) = M(:,5);

%reduce to every 5th element
x=x(1:5:end,:);
y=y(1:5:end,:);
x(end)=x(1);
y(end)=y(1);

% x = M(:,4);
% y = M(:,5);

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
plot(x,y,'k','LineWidth',2)
hold on
scatter(x(1),y(1),'*','r')
scatter(x(1),y(1),'x','r','LineWidth',2)
xlabel('X')
ylabel('Y')
axis tight
axis equal

%plot in thirds
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
% % plot(x1,y1,'r','LineWidth',2)
% % plot(x2,y2,'g','LineWidth',2)
% % plot(x3,y3,'b','LineWidth',2)
% scatter(x(1),y(1),'x','r','LineWidth',2)
% %scatter(x(100),y(100),'o')
% xlabel('X')
% ylabel('Y')
% axis tight
% axis equal
% 
