function [dam_matrix,wave_weight_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p)

% initialize
x = 1:size(lake,2);
y = 1:size(lake,1);
[X,Y] = meshgrid(x,y);
X = X*p.dx; Y = Y*p.dy;
wave_weight_matrix = nan(size(lake));
dam_matrix = zeros(size(lake));
cells2trash = [];

% find the shoreline
[shoreline] = addidshoreline(lake,~lake);

% find number of first order lakes
[F_lake_all,~,~,~] = find_first_order_lakes(lake);

for ff = 1:length(F_lake_all)
    F_lake = F_lake_all{ff};
    if length(find(F_lake))<2
        continue
    end
    clearvars fetch_sl_cells indshoreline_ordered WaveArea_cell
    
    %order the shoreline and islands
    % new order the shoreline code
    [indshoreline_ocw,~,cells2trash_ff,p] = order_shoreline_bwbound(F_lake,p);
    if isfield(p,'boundary')
        dam_matrix= [];
        wave_weight_matrix= [];
        indshoreline_ordered= [];
        cells2trash= [];
        return
    end
    cells2trash = [cells2trash;cells2trash_ff];
    
    if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
        break
    end
    
    for l = 1: length(indshoreline_ocw)
        indshoreline_ordered{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
        fetch_sl_cells{l,1}(:,1) = X(indshoreline_ordered{l,1});
        fetch_sl_cells{l,1}(:,2) = Y(indshoreline_ordered{l,1});
    end
    % calculate wave weighted (sqrt(F)*cos(theta-phi))
    [WaveArea_cell,~] = fetch_vis_approx(fetch_sl_cells);% first is wave, second is fetch!!
    clearvars erodedind
    %         end
    
    % Damage the shoreline
    indshoreline_ordered = cell2mat(indshoreline_ordered);
    %         [shoreline] = addidshoreline_cardonly(lake,land); % edges only
    wave_weighting = cell2mat(WaveArea_cell);
    wave_weighting = cell2mat(WaveArea_cell)./p.Ao;
    wave_weight_matrix(indshoreline_ordered) = wave_weighting;
    dam = p.dt*p.Kwave*shoreline(indshoreline_ordered).*wave_weighting*p.So./p.dxo;
    dam_matrix(indshoreline_ordered) = dam_matrix(indshoreline_ordered)+dam; % have to do this because if a shoreline cell is an island on the big one and has a second order lake in it, we don't want the damage from the second order lake to be the only damage it receives. need to add up damage from both sides.
    
    
    
end