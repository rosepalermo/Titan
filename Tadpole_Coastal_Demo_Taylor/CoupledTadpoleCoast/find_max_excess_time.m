function        [time_dam_excess,p] = find_max_excess_time(time_excess,ind_new,lake,p)
% this function takes the excess indices, the excess time, and the new
% indices to find the max excess time of the eight connected neighbors of
% the new indices



[indxsldif,indysldif] = ind2sub(size(lake),ind_new);
time_dam_excess = zeros(length(ind_new),1);

for xxx = 1: length(ind_new)
    if indxsldif(xxx) == 1 | indxsldif(xxx) == length(lake) | indysldif(xxx) == 1 | indysldif(xxx) == length(lake)
        time_dam_excess = zeros(length(ind_new),1);
        p.t = p.tf; % hit the boundary, end the simulation
        return
    else
        U = sub2ind(size(lake),indxsldif(xxx),indysldif(xxx)-1);
        UL = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx)-1);
        L = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx));
        DL = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx)+1);
        D = sub2ind(size(lake),indxsldif(xxx),indysldif(xxx)+1);
        DR = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx)+1);
        R = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx));
        UR = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx)-1);
        
        time_dam_excess(xxx) = max(time_excess([U, UL, L, DL, D, DR, R, UR]));
    end
end
end
