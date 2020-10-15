% calculate wavelet power spectrum
x0 = x; y0 = y;

% get rid of duplicate points

M = [x0 y0];
% indices to unique values in column 3
[~, ind] = unique(M, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(M, 1), ind);
% duplicate values
duplicate_val = [M(duplicate_ind, 1) M(duplicate_ind, 2)];
M(duplicate_ind,:)=[];
% fetch(duplicate_ind) = [];


x0=M(:,1);
y0=M(:,2);

% add first AND SECOND points to the end for meander
x = [x0;x0(1:2)];
y = [y0;y0(1:2)];



% transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% streamwise distance
[theta, dtheta, deltad] = meander_titan(x,y);
x = x(1:end-2);
y = y(1:end-2);
theta = theta(1:end-1);


n=2;

[period{1},eq14holder] = dowave_greece(theta,deltad,n,x,y,'test',0,x,2);
holdersize= size(eq14holder');
eq14save{1} = eq14holder';

thisone = size(eq14save);