% test interp methods for subtimesteps
load('fetch_for_interp_test.mat')
x = 1:length(lake_save{1});
y = 1:length(lake_save{1});
[X,Y] = meshgrid(x,y);


% test to make sure this interp function calculates the same as below
% Yup! yay!
[fetch_interp] = interp_fetch_for_ind(lake_save{18},shorelinediff{18},fetch_save{18});



%%
fetchind = 18;

% plot fetch
% figure()
% imagesc(fetch_save{fetchind})
% hold on
[yind,xind] = ind2sub(size(lake_save{1}),shorelinediff{fetchind});
% plot the points in the next time step that need to be calculated from
% previous fetch
% scatter(xind,yind,'r')

F_nearest = scatteredInterpolant(X(~isnan(fetch_save{fetchind})),Y(~isnan(fetch_save{fetchind})),fetch_save{fetchind}(~isnan(fetch_save{fetchind})),'nearest');
vq_nearest = F_nearest(xind,yind);
% scatter3(xind,yind,vq_nearest,[],vq_nearest)

F_linear = scatteredInterpolant(X(~isnan(fetch_save{fetchind})),Y(~isnan(fetch_save{fetchind})),fetch_save{fetchind}(~isnan(fetch_save{fetchind})),'linear');
vq_linear = F_linear(xind,yind);
% scatter3(xind,yind,vq_linear,[],vq_linear)

F_natural = scatteredInterpolant(X(~isnan(fetch_save{fetchind})),Y(~isnan(fetch_save{fetchind})),fetch_save{fetchind}(~isnan(fetch_save{fetchind})),'natural');
vq_natural = F_natural(xind,yind);
% scatter3(xind,yind,vq_natural,[],vq_natural)

indshoreline = sub2ind(size(lake_save{fetchind}),X(~isnan(fetch_save{fetchind})),Y(~isnan(fetch_save{fetchind})));
F_meanneighbor = mean_neighboring_fetch(lake_save{fetchind},shorelinediff{fetchind},indshoreline,fetch_save{fetchind}(~isnan(fetch_save{fetchind})));

actualfetch = fetch_save{fetchind+1}(shorelinediff{fetchind});
% scatter3(xind,yind,actualfetch,[],actualfetch)

% figure()
% plot(actualfetch,'k','LineWidth',1.5)
hold on
% plot(vq_nearest)
% plot(vq_linear)
plot(vq_natural)
% plot(F_meanneighbor)
% legend('actual','nearest','linear','natural','meanneighbor')

% figure()
% hold on
% plot(abs(vq_nearest-actualfetch))
% plot(abs(vq_linear-actualfetch))
% plot(abs(vq_natural-actualfetch))
% plot(abs(F_meanneighbor-actualfetch))
% legend(sprintf('nearest %f',sum(abs(vq_nearest-actualfetch))),sprintf('linear %f',sum(abs(vq_linear-actualfetch))),sprintf('natural %f',sum(abs(vq_natural-actualfetch))),sprintf('meanneighbor %f',sum(abs(F_meanneighbor-actualfetch))))
% % sum(abs(vq_nearest-actualfetch))
% sum(abs(vq_linear-actualfetch))
% sum(abs(vq_natural-actualfetch))
% sum(abs(F_meanneighbor-actualfetch))


