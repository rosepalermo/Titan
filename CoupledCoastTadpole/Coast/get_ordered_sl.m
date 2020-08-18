function [indshoreline_ordered] = get_ordered_sl(lake,p)

% initialize
x = 1:size(lake,2);
y = 1:size(lake,1);
[X,Y] = meshgrid(x,y);
X = X*p.dx; Y = Y*p.dy;

% find number of first order lakes
[F_lake_all,~,~,~] = find_first_order_lakes(lake);

for ff = 1:length(F_lake_all)
    F_lake = F_lake_all{ff};
    if length(find(F_lake))<2
        continue
    end    
    %order the shoreline and islands
    % new order the shoreline code
    [indshoreline_ocw,~,cells2trash_ff,p] = order_shoreline_bwbound(F_lake,p);
    if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
        break
    end
    for l = 1: length(indshoreline_ocw)
        indshoreline_ordered{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
        fetch_sl_cells{l,1}(:,1) = X(indshoreline_ordered{l,1});
        fetch_sl_cells{l,1}(:,2) = Y(indshoreline_ordered{l,1});
    end

    indshoreline_ordered = cell2mat(indshoreline_ordered);

end