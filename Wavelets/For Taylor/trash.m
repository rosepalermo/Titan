filename = 'Scotland_gp.csv';
M = csvread(filename);
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
lat = 111360;
lon = 60772;
x0 = (x0-min(x0)) * lon;%changed from lon/lat*lon to just lon
y0 = (y0-min(y0)) * lat;
dist = sqrt((x0(2:end) - x0(1:end-1)).^2 +  (y0(2:end) - y0(1:end-1)).^2);
dist_off = find(dist>1000);
