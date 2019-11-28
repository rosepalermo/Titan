function [dam,damcorn,corners,indshoreline] = find_fetch(F_lake,X,Y,lake);

disp('fetch')
clearvars fetch_sl_cells indshoreline WaveArea_cell
disp('ordering')
% order the shoreline
[indshoreline_ocw,~,cells2trash] = order_shoreline_bwbound(F_lake);
if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
    return
end
% my original order the shoreline code
%             [indshoreline_ocw,cells2trash] = order_cw_lastpoint(F_lake,shoreline); % ccw ordered ind = indshoreline
%             keepme = 1:length(indshoreline_ocw);
disp('ordered')
for l = 1: length(indshoreline_ocw)
    indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
    fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
    fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
end
disp('calculating wave')
% calculate wave weighted (sqrt(F)*cos(theta-phi))
[~,WaveArea_cell] = fetch_vis_approx(fetch_sl_cells);
% Damage the shoreline
indshoreline = cell2mat(indshoreline);
%         [shoreline] = addidshoreline_cardonly(lake,land); % edges only
[shoreline] = addidshoreline(F_lake,~F_lake); % corners and edges
dam = cell2mat(WaveArea_cell);
corners = setdiff(find(shoreline),indshoreline);
if ~isempty(cells2trash)
    [c2t]=sub2ind(size(X),cells2trash(:,1),cells2trash(:,2)); % find the cells that arent the trash cells
    if exist('c2t') & ~isempty(corners)
        corners = corners(~ismember(corners,c2t));
    end
end

if ~isempty(corners)
    corners = corners(shoreline(corners)<1.5); % if less than 1.5, only a corner. not also a side
    [damcorn] = damagecorners(lake,corners,indshoreline,dam);
else
    corners =[];
    damcorn = [];
end

dam = dam.*shoreline(indshoreline);

end
