function [sum_dam_excess,ind_excess,time_excess,dam_matrix] = calc_excess_dam(strength,dam_matrix,p)

% calculate the excess damage
ind_excess = find(strength<0);
dam_excess = -strength(ind_excess);
sum_dam_excess = sum(dam_excess);
if sum_dam_excess>0
    time_excess_temp = dam_excess./(dam_matrix(ind_excess)).*p.dt;
%     dam_matrix(ind_excess) - dam_excess
    dam_matrix(ind_excess) = dam_matrix(ind_excess) - dam_excess;
    time_excess = zeros(size(dam_matrix));
    time_excess(ind_excess) = time_excess_temp;
else
    time_excess = [];
end
if any(any(isinf(time_excess)))
end
end
