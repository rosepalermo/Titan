% prep contour for Titan Modeling
close all;clear all

load('example_Rose_800x800.mat')

% 1 m contour at t = 1
figure()
tile = repmat(output,2);
imagesc(tile(:,:,1))
hold on

maxc = 50*floor(max(max(tile(:,:,1)))/50);
[C,h] = contour(tile(:,:,1),[1,1]);

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
figure()
[C,h] = contour(tile(:,:,1),[1,1]);
hold on
plot(contourxy{60,1}(1,:), contourxy{60,1}(2,:))
set(gca,'Ydir','reverse')
grid
x_1m_t1 = contourxy{60,1}(1,:);
y_1m_t1 = contourxy{60,1}(2,:);

clearvars -except x_1m_t1 y_1m_t1 tile

%% 1 m contour at t = 3

figure()

imagesc(tile(:,:,3))
hold on

maxc = 50*floor(max(max(tile(:,:,1)))/50);
[C,h] = contour(tile(:,:,3),[1,1]);

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

% find longest cell
val=cellfun(@(x) numel(x),contourxy);
out=contourxy(val==max(val));
figure
[C,h] = contour(tile(:,:,3),[1,1]);
hold on
plot(out{1,1}(1,:),out{1,1}(2,:),'r')
set(gca,'Ydir','reverse')


x_1m_t3 = out{1,1}(1,:);
y_1m_t3 = out{1,1}(2,:);

save('xycontours.mat','x_1m_t1','y_1m_t1','x_1m_t3','y_1m_t3')
