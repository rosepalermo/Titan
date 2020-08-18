function [adt,dam_matrix,p] = getCADT_uniform(p,lake);

%This function finds the damage associated with the maximum perimeter of 
% coastline. Based on max damage, it calculates an adaptive timestep
if isfield(p,'boundary')
    adt = nan;
    dam_matrix = nan;
    return
end

[dam_matrix] = get_dam_uniform(lake,p);

if max(dam_matrix,[],'all')>0.1
adt = ceil(max(dam_matrix,[],'all'))*10;

dam_matrix = dam_matrix./adt;
else
    adt = 1;
end

end

