% prep contour for Titan Modeling
% 0.5 contour because 1m contour did not close a loop

close all;clear all


load('example_800x800_beta1_6.mat')
output = solution;
% 1 m contour at t = 1
figure()
tile = repmat(output,3);
imagesc(tile(:,:,1))
hold on

maxc = 50*floor(max(max(tile(:,:,1)))/50);
[C,h] = contour(tile(:,:,1),[0.5,1]);

idx = find(C(1,:) == 1);
Llen = C(2,idx);
for k1 = 1:length(idx)
    conturc{k1,1} = C(:,idx(k1)+1 : idx(k1)+1+Llen(k1)-1);
    conturc{k1,2} = C(:,idx(k1)+1 : idx(k1)+1+Llen(k1)-1);
    if ~any(conturc{k1,1}<1000) | ~any(conturc{k1,1}>2000)
        contourxy{k1,1} = conturc{k1,1};
    end
    if ~any(conturc{k1,2}<1000) | ~any(conturc{k1,2}>2000)
        contourxy{k1,2} = conturc{k1,2};
    end
end

length_cells = cellfun(@length,contourxy,'uni',false);
length_cells = cell2mat(length_cells);
[length_cells_sort,sortind] = sort(length_cells);
length_cells_sort = flipud(length_cells_sort);
sortind = flipud(sortind);
contourxy_sorted = contourxy(sortind,:);


figure()
[C,h] = contour(tile(:,:,1),[0.5,1]);
hold on
plot(contourxy_sorted{3,1}(1,:), contourxy_sorted{3,1}(2,:),'r')
set(gca,'Ydir','reverse')
grid
x_0_5m_t1v1 = contourxy_sorted{1,1}(1,:);
y_0_5m_t1v1 = contourxy_sorted{1,1}(2,:);

% clearvars -except x_0_5m_t1v1 x_0_5m_t1v1 tile

%% 1 m contour at t = 3

figure()

imagesc(tile(:,:,3))
hold on

maxc = 50*floor(max(max(tile(:,:,3)))/50);
[C,h] = contour(tile(:,:,3),[0.5,1]);

idx = find(C(1,:) == 1);
Llen = C(2,idx);
for k1 = 1:length(idx)
    conturc{k1,1} = C(:,idx(k1)+1 : idx(k1)+1+Llen(k1)-1);
    conturc{k1,2} = C(:,idx(k1)+1 : idx(k1)+1+Llen(k1)-1);
    if ~any(conturc{k1,1}<400) | ~any(conturc{k1,1}>1300)
        contourxy{k1,1} = conturc{k1,1};
    end
    if ~any(conturc{k1,2}<200) | ~any(conturc{k1,2}>900)
        contourxy{k1,2} = conturc{k1,2};
    end
end
% get rid of empty cells
contourxy = contourxy(find(~cellfun(@isempty,contourxy)));

% sort cells
length_cells = cellfun(@length,contourxy,'uni',false);
length_cells = cell2mat(length_cells);
[length_cells_sort,sortind] = sort(length_cells);
length_cells_sort = flipud(length_cells_sort);
sortind = flipud(sortind);
contourxy_sorted = contourxy(sortind,:);


figure()
[C,h] = contour(tile(:,:,3),[0.5,1]);
hold on
plot(contourxy_sorted{3,1}(1,:), contourxy_sorted{3,1}(2,:),'r')
set(gca,'Ydir','reverse')
grid
x_0_5m_t3v1 = contourxy_sorted{1,1}(1,:);
y_0_5m_t3v1 = contourxy_sorted{1,1}(2,:);

%% 1 m contour at t = 2

figure()

imagesc(tile(:,:,2))
hold on

maxc = 50*floor(max(max(tile(:,:,2)))/50);
[C,h] = contour(tile(:,:,2),[0.5,1]);

idx = find(C(1,:) == 1);
Llen = C(2,idx);
for k1 = 1:length(idx)
    conturc{k1,1} = C(:,idx(k1)+1 : idx(k1)+1+Llen(k1)-1);
    conturc{k1,2} = C(:,idx(k1)+1 : idx(k1)+1+Llen(k1)-1);
end
% get rid of empty cells
contourxy = conturc;
contourxy = contourxy(find(~cellfun(@isempty,contourxy)));

% find longest cell
val=cellfun(@(x) numel(x),contourxy);
out=contourxy(val==max(val));
figure
[C,h] = contour(tile(:,:,2),[0.5,1]);
hold on
plot(out{1,1}(1,:),out{1,1}(2,:),'r')
set(gca,'Ydir','reverse')


x_0_5m_t2v1 = out{1,1}(1,:);
y_0_5m_t2v1 = out{1,1}(2,:);


save('xy_beta1_6_v1.mat','x_0_5m_t1v1','y_0_5m_t1v1','x_0_5m_t2v1','y_0_5m_t2v1','x_0_5m_t3v1','y_0_5m_t3v1')
