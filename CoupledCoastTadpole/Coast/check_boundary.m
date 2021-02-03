function [p] = check_boundary(lake,p)

% check if shoreline is on a boundary
% find number of first order lakes
[F_lake_all,~,~,~] = find_first_order_lakes(lake,p);
for ff = 1:length(F_lake_all)
    F_lake = F_lake_all{ff};
    if length(find(F_lake))<2
        continue
    end    
    % use order_shoerline_bwbound to find if boundary
    [~,~,~,p] = order_shoreline_bwbound(F_lake,p);
    if isfield(p,'boundary') && ~isfield(p,'sl_analysis')
        dam_matrix= [];
        indshoreline= [];
        return
    end
    
end
end