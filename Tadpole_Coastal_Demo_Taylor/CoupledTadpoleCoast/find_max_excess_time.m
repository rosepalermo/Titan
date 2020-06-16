function        time_dam_excess = find_max_excess_time(ind_excess,time_excess,ind_new,lake);
% this function takes the excess indices, the excess time, and the new
% indices to find the max excess time of the eight connected neighbors of
% the new indices



[indysldif,indxsldif] = ind2sub(size(lake),ind_new);
time_dam_excess = zeros(length(ind_new),1);


for xxx = 1: length(ind_new)
    
    U = sub2ind(size(lake),indxsldif(xxx),indysldif(xxx)-1);
    UL = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx)-1);
    L = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx));
    DL = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx)+1);
    D = sub2ind(size(lake),indxsldif(xxx),indysldif(xxx)+1);
    DR = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx)+1);
    R = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx));
    UR = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx)-1);
    
    Um = ismember(U,ind_excess);
    ULm = ismember(UL,ind_excess);
    Lm = ismember(L,ind_excess);
    DLm = ismember(DL,ind_excess);
    Dm = ismember(D,ind_excess);
    DRm = ismember(DR,ind_excess);
    Rm = ismember(R,ind_excess);
    URm = ismember(UR,ind_excess);
    
    temporary = 0;
    num = Um + ULm + Lm + DLm + Dm + DRm + Rm + URm;
    if num > 0
        if Um
            temporary = [temporary nansum(time_excess(find(ind_excess == U)))];
        end
        if ULm
            temporary = [temporary nansum(time_excess(find(ind_excess == UL)))];
        end
        if Lm
            temporary = [temporary nansum(time_excess(find(ind_excess == L)))];
        end
        if DLm
            temporary = [temporary nansum(time_excess(find(ind_excess == DL)))];
        end
        if Dm
            temporary = [temporary nansum(time_excess(find(ind_excess == D)))];
        end
        if DRm
            temporary = [temporary nansum(time_excess(find(ind_excess == DR)))];
        end
        if Rm
            temporary = [temporary nansum(time_excess(find(ind_excess == R)))];
        end
        if URm
            temporary = [temporary nansum(time_excess(find(ind_excess == UR)))];
        end
        
        time_dam_excess(xxx) = max(temporary);
    else
        time_dam_excess(xxx) = 0;
    end
end
end
