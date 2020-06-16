function [adt,dam_matrix,wave_matrix,p] = getCADT(p,g,lake);

%This function finds the damage associated with the maximum fetch in the
%landscape. Based on max damage, it calculates an adaptive timestep
if isfield(p,'boundary')
    adt = nan;
    dam_matrix = nan;
    return
end

x = 1:size(lake,2);
y = 1:size(lake,1);
[X,Y] = meshgrid(x,y);
X = X*p.dx; Y = Y*p.dy;
land = ~lake;

[F_lake_all,~,~,~] = find_first_order_lakes(lake);

erodedind_save = [NaN];
dam_matrix = zeros(size(lake));
wave_matrix = zeros(size(lake));

for ff = 1:length(F_lake_all)
    F_lake = F_lake_all{ff};
    if length(find(F_lake))<2
        continue
    end
    [indshoreline_ocw,~,cells2trash] = order_shoreline_bwbound(F_lake);
    if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
        break
    end
    for l = 1: length(indshoreline_ocw)
        indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
        fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
        fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
    end
    
    [WaveArea_cell,~] = fetch_vis_approx(fetch_sl_cells);
    wave_weighting = cell2mat(WaveArea_cell);
    [shoreline] = addidshoreline(lake,land);
    indshoreline = cell2mat(indshoreline);
    dam = p.dt*p.Kcoast*shoreline(indshoreline).*wave_weighting;
    dam_matrix(indshoreline) = dam;
    wave_matrix(indshoreline) = wave_weighting;
    
    maxdam(ff)= max(dam);
    clearvars indshoreline dam fetch_sl_cells
end

% because damage is normalized to 1, we need as many time steps as if we
% round up the damage. When <1, only one time step, so not really adaptive.
% When >1, we divide by ceil(damage) so that it is never greater than 1.
adt = ceil(max(maxdam));

dam_matrix = dam_matrix./adt;

end

