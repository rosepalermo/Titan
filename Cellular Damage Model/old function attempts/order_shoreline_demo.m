%ordering the shoreline demo

% small test polygon
load test_polygon.mat
lake_small_test = polygon;

% actual model example
load wavet2v2_lake.mat


%% first the test with the small lake
% 
% [shoreline_small_test] = addidshoreline_cardonly(lake_small_test,~lake_small_test);
% [indshoreline_ocw_small,cells2trash] = order_cw_lastpoint(lake_small_test,shoreline_small_test); % ccw ordered ind = indshoreline
% 
% figure()
% imagesc(lake_small_test');
% hold on
% plot(indshoreline_ocw_small{1,1}(:,1),indshoreline_ocw_small{1,1}(:,2),'c','LineWidth',2)

%% test with the model lake

[L,fol] = find_first_order_lakes(lake);
F_lake = (L == fol(5)); % This is the longest one.

shoreline = addidshoreline_cardonly(F_lake,~F_lake); %rewrite shoreline to be card only for fetch.. will damage corners separately later
tic
[indshoreline_ocw,cells2trash] = order_cw_lastpoint(F_lake,shoreline); % ccw ordered ind = indshoreline
toc

figure()
imagesc(F_lake');
hold on
plot(indshoreline_ocw{1,1}(:,1),indshoreline_ocw{1,1}(:,2),'c','LineWidth',2)