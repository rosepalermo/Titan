function        [time_dam_excess,ind_rec_excess,p] = pass_excess_time(time_excess,ind_new,ind_sl_new,ind_excess,lake,p)
% this function takes the excess indices, the excess time, and the new
% indices to find the max excess time of the eight connected neighbors of
% the new indices



[indxsldif,indysldif] = ind2sub(size(lake),ind_new);
for xxx = 1:length(ind_new)
    if indxsldif(xxx) == 1 | indxsldif(xxx) == length(lake) | indysldif(xxx) == 1 | indysldif(xxx) == length(lake)
        time_dam_excess = [];
        ind_rec_excess = [];
        p.boundary = 1; % hit the boundary, end the simulation
        return
    end
end
[indrslexcess,indcslexcess] = ind2sub(size(lake),ind_excess);
time_dam_excess = cell(length(ind_excess),1);

for xxx = 1: length(ind_excess)
    
    D = sub2ind(size(lake),indrslexcess(xxx),indcslexcess(xxx)-1);
    DL = sub2ind(size(lake),indrslexcess(xxx)-1,indcslexcess(xxx)-1);
    L = sub2ind(size(lake),indrslexcess(xxx)-1,indcslexcess(xxx));
    UL = sub2ind(size(lake),indrslexcess(xxx)-1,indcslexcess(xxx)+1);
    U = sub2ind(size(lake),indrslexcess(xxx),indcslexcess(xxx)+1);
    UR = sub2ind(size(lake),indrslexcess(xxx)+1,indcslexcess(xxx)+1);
    R = sub2ind(size(lake),indrslexcess(xxx)+1,indcslexcess(xxx));
    DR = sub2ind(size(lake),indrslexcess(xxx)+1,indcslexcess(xxx)-1);
    

    directions = [U, UL, L, DL, D, DR, R, UR]; % all 8 directions.
    
    ind_rec_excess_temp{xxx} = directions(ismember(directions,ind_sl_new));
    if ~isempty(ind_rec_excess_temp{xxx})
        time_dam_excess_temp{xxx} = time_excess(ind_excess(xxx))./sum(find(ind_rec_excess_temp{xxx}))*ones(size(ind_rec_excess_temp{xxx}));
    else
        time_dam_excess_temp{xxx} = [];
    end
%     time_dam_excess(xxx) = max(time_excess([U, UL, L, DL, D, DR, R, UR]));
end
ind_rec_excess = cell2mat(ind_rec_excess_temp);
time_dam_excess = cell2mat(time_dam_excess_temp);
end
