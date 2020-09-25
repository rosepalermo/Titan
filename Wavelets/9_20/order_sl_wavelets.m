% find longest shoreline in each of the time steps and order it

% load data
load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/091820/rand_hwave_luniform.mat')
addpath(genpath('/Users/rosepalermo/Documents/GitHub/Titan/CoupledCoastTadpole'))
% identify shoreline
ind = size(g.output,3);
if isfield(p,'boundary')
    ind = size(g.output,3)-1;
end

for i = 1:ind
    
    lake = g.output(:,:,i)<p.sealevel_init;
% find number of first order lakes
[F_lake_all,~,~,~] = find_first_order_lakes(lake);

for ff = 1:length(F_lake_all)
    length_sl(ff) = length(find(F_lake_all{ff}));
end
    ind_lake = find(length_sl==max(length_sl));

    F_lake = F_lake_all{ind_lake};
    if length(find(F_lake))<2
        continue
    end
    clearvars fetch_sl_cells indshoreline_ordered WaveArea_cell
    
    %order the shoreline and islands
    % new order the shoreline code
    [indshoreline_ocw,~,~,p] = order_shoreline_bwbound(F_lake,p);
    if isfield(p,'boundary')
        dam_matrix= [];
        wave_weight_matrix= [];
        indshoreline_ordered= [];
        cells2trash= [];
        return
    end
    
    if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
        break
    end
    
    sl_ind = indshoreline_ocw{1};
    x = 1:p.dx:p.Nx*p.dx;
    y = 1:p.dy:p.Ny*p.dy;
    [X,Y] = meshgrid(x,y);
    indshoreline_ordered = sub2ind(size(X),sl_ind(:,1),sl_ind(:,2));
    fetch_sl_cells(:,1) = X(indshoreline_ordered);
    fetch_sl_cells(:,2) = Y(indshoreline_ordered);
    
    sl_save{i} = fetch_sl_cells;

end


