function [mnf] = mean_neighboring_fetch(lake,shorelinediff,indshoreline,fetch)


% this function finds the mean of the Wave weighting of all of the edge
% points adjacent to each corner point



[indysldif,indxsldif] = ind2sub(size(lake),shorelinediff);
if any(ismember([1,size(lake,2)],indysldif)) | any(ismember([1,size(lake,2)],indxsldif))
    mnf = zeros(length(shorelinediff),1);
    return
end
mnf = zeros(length(shorelinediff),1);


for xxx = 1: length(shorelinediff)
    
    U = sub2ind(size(lake),indxsldif(xxx),indysldif(xxx)-1);
    UL = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx)-1);
    L = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx));
    DL = sub2ind(size(lake),indxsldif(xxx)-1,indysldif(xxx)+1);
    D = sub2ind(size(lake),indxsldif(xxx),indysldif(xxx)+1);
    DR = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx)+1);
    R = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx));
    UR = sub2ind(size(lake),indxsldif(xxx)+1,indysldif(xxx)-1);
    
    Um = ismember(U,indshoreline);
    ULm = ismember(UL,indshoreline);
    Lm = ismember(L,indshoreline);
    DLm = ismember(DL,indshoreline);
    Dm = ismember(D,indshoreline);
    DRm = ismember(DR,indshoreline);
    Rm = ismember(R,indshoreline);
    URm = ismember(UR,indshoreline);
    
    temporary = 0;
    num = Um + ULm + Lm + DLm + Dm + DRm + Rm + URm;
    if num > 0
        if Um
            temporary = temporary + nansum(fetch(indshoreline == U));
        end
        if ULm
            temporary = temporary + nansum(fetch(indshoreline == UL));
        end
        if Lm
            temporary = temporary + nansum(fetch(indshoreline == L));
        end
        if DLm
            temporary = temporary + nansum(fetch(indshoreline == DL));
        end
        if Dm
            temporary = temporary + nansum(fetch(indshoreline == D));
        end
        if DRm
            temporary = temporary + nansum(fetch(indshoreline == DR));
        end
        if Rm
            temporary = temporary + nansum(fetch(indshoreline == R));
        end
        if URm
            temporary = temporary + nansum(fetch(indshoreline == UR));
        end
        
        mnf(xxx) = temporary/num;
    else
        mnf(xxx) = 0;
    end
end
end