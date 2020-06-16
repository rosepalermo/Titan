% test mean neighbor subtimesteps
load('fetch_for_interp_test.mat')
x = 1:length(lake_save{1});
y = 1:length(lake_save{1});
[X,Y] = meshgrid(x,y);

fetchind = 16;

% plot fetch
figure()
imagesc(fetch_save{fetchind})
hold on
[yind,xind] = ind2sub(size(lake_save{fetchind}),shorelinediff{fetchind});
% plot the points in the next time step that need to be calculated from
% previous fetch
scatter(xind,yind,'r')
indshoreline = sub2ind(size(lake_save{fetchind}),X(~isnan(fetch_save{fetchind})),Y(~isnan(fetch_save{fetchind})));
F_meanneighbor = mean_neighboring_fetch(lake_save{fetchind},shorelinediff{fetchind},indshoreline,fetch_save{fetchind}(~isnan(fetch_save{fetchind})));
% scatter3(xind,yind,vq_meanneighbor,[],vq_meanneighbor)

actualfetch = fetch_save{fetchind+1}(shorelinediff{fetchind});
% scatter3(xind,yind,actualfetch,[],actualfetch)

figure()
plot(actualfetch)
hold on
plot(F_meanneighbor)
legend('actual','meanneighbor')

figure()
hold on
plot(abs(F_meanneighbor-actualfetch))
plot(zeros(length(F_meanneighbor)))
sum(abs(F_meanneighbor-actualfetch))

plot
