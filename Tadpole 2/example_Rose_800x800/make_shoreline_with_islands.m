% lake with islands example

load('example_800x800_beta1_6.mat');

input = solution(:,:,2);
tile = repmat(input,3);
tile_lake = tile;
tile_lake = tile_lake(484:1297, 678:1563);
tile_lake = tile_lake > 76;
% tile_lake(tile>1) = 0;
% tile_lake(tile<=1) = 1;
% tile_lake = tile_lake(ceil(size(tile_lake,2)/6):4*floor(size(tile_lake,1)/6),ceil(size(tile_lake,2)/6):4*floor(size(tile_lake,1)/6));
% return
[L,fol,foi] = find_first_order_lakes(tile_lake);

F_lake = L==28;
[indshoreline_ocw,keepme,cells2trash] = order_shoreline_bwbound(F_lake);
cells_with_shorelines = find(keepme);
[X,Y] = meshgrid(1:length(tile_lake),1:length(tile_lake));
for l = 1: length(cells_with_shorelines)
    indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{cells_with_shorelines(l),1}(:,1),indshoreline_ocw{cells_with_shorelines(l),1}(:,2));
    fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
    fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
end
fetch_sl_cells{244,2} = 'fol';
for ll = 1:243
fetch_sl_cells{ll,2} = 'foi';
end
