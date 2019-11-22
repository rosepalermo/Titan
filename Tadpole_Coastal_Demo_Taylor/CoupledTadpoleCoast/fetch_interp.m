close all; clear all

load('initialtopo.mat')
x = 1:400; y = x;
[X,Y] = meshgrid(x,y);

figure()
ind_all = [NaN];
fetch_all = [NaN];
fetchcorn_all = [NaN];
corners_all = [NaN];

% iterate over 0.1 m and calculate damage for uniform erosion
for SL = 10
    % define the shoreline
    lake = init<SL;
    
    [F_lake_all,~,~,~] = find_first_order_lakes(lake);
    for ff = 1:length(F_lake_all)
        F_lake = F_lake_all{ff};
        if length(find(F_lake))<2
            continue
        end
        [fetch,fetch_corn,corners,indshoreline] = find_fetch(F_lake,X,Y,lake);
        fetch_all = [fetch_all;fetch];
        fetchcorn_all = [fetchcorn_all;fetch_corn];
        ind_all = [ind_all;indshoreline];
        corners_all = [corners_all;corners];
    end
end
ind_all = ind_all(2:end);
fetch_all = fetch_all(2:end);
fetchcorn_all = fetchcorn_all(2:end);
corners_all = corners_all(2:end);

% plot
[y_ind,x_ind] = ind2sub(size(lake),ind_all);
scatter3(x_ind,y_ind,fetch_all,[],fetch_all)

view(2)

caxis([0 10])
axis equal
x_all = x_all(2:end); y_all = y_all(2:end);


%%
SL1 = 1; SL2 = 10;

lake1 = init<SL1;
lake2 = init<SL2;


ind_all1 = [NaN];
fetch_all1 = [NaN];
fetchcorn_all1 = [NaN];
corners_all1 = [NaN];

    [F_lake_all,~,~,~] = find_first_order_lakes(lake1);
    for ff = 1:length(F_lake_all)
        F_lake = F_lake_all{ff};
        if length(find(F_lake))<2
            continue
        end
        [fetch,fetch_corn,corners,indshoreline] = find_fetch(F_lake,X,Y,lake1);
        fetch_all1 = [fetch_all1;fetch];
        fetchcorn_all1 = [fetchcorn_all1;fetch_corn];
        ind_all1 = [ind_all1;indshoreline];
        corners_all1 = [corners_all1;corners];
    end
    
ind_all1 = ind_all1(2:end);
fetch_all1 = fetch_all1(2:end);
fetchcorn_all1 = fetchcorn_all1(2:end);
corners_all1 = corners_all1(2:end);
ind_sl1 = [ind_all1;corners_all1];
fetch_sl1 = [fetch_all1;fetchcorn_all1];


ind_all2 = [NaN];
fetch_all2 = [NaN];
fetchcorn_all2 = [NaN];
corners_all2 = [NaN];

        [F_lake_all,~,~,~] = find_first_order_lakes(lake2);
    for ff = 1:length(F_lake_all)
        F_lake = F_lake_all{ff};
        if length(find(F_lake))<2
            continue
        end
        [fetch,fetch_corn,corners,indshoreline] = find_fetch(F_lake,X,Y,lake2);
        fetch_all2 = [fetch_all2;fetch];
        fetchcorn_all2 = [fetchcorn_all2;fetch_corn];
        ind_all2 = [ind_all1;indshoreline];
        corners_all2 = [corners_all1;corners];
    end

ind_all2 = ind_all1(2:end);
fetch_all2 = fetch_all1(2:end);
fetchcorn_all2 = fetchcorn_all1(2:end);
corners_all2 = corners_all1(2:end);
ind_sl2 = [ind_all2;corners_all2];
fetch_sl2 = [fetch_all2;fetchcorn_all2];
    
[y_ind1,x_ind1] = ind2sub(size(lake1),ind_sl2);
[y_ind2,x_ind2] = ind2sub(size(lake2),ind_sl1);

% figure
% scatter3(x_ind1,y_ind1,dam1,[],dam1)
% view(2)
% hold on
% scatter3(x_ind2,y_ind2,dam2,[],dam2)

[sl2] = addidshoreline(lake2,~lake2); % corners and edges
sl2(sl2>0) = 1;
intermediate = lake2 - lake1 + sl2;
% figure()
% imagesc(intermediate)
[y_indint,x_indint] = ind2sub(size(intermediate),find(intermediate));

fetch_all_sl12 = [fetch_sl1; fetch_sl2];
[y_fetch,x_fetch] = ind2sub(size(lake),[ind_sl1;ind_sl2]);

% daminterp(daminterp=0)=NaN;
Vq = griddata(x_fetch,y_fetch,fetch_all_sl12,x_indint,y_indint);

figure
scatter3(x_indint,y_indint,Vq,[],Vq)
view(2)
caxis([0 10])
axis equal

figure(); imagesc(init); caxis([SL1 SL2])
axis equal
set(gca,'Ydir','normal')