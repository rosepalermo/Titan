function [damcorn] = damagecorners(lake,corners,indshoreline,dam)


% this function finds the mean of the Wave weighting of all of the edge
% points adjacent to each corner point



[indxcorners,indycorners] = ind2sub(size(lake),corners);
if any(ismember([1,size(lake,2)],indycorners)) | any(ismember([1,size(lake,2)],indxcorners))
    damcorn = zeros(length(corners),1);
    return
end
damcorn = zeros(length(corners),1);


for xxx = 1: length(corners)
    
    U = sub2ind(size(lake),indxcorners(xxx),indycorners(xxx)-1);
    UL = sub2ind(size(lake),indxcorners(xxx)-1,indycorners(xxx)-1);
    L = sub2ind(size(lake),indxcorners(xxx)-1,indycorners(xxx));
    DL = sub2ind(size(lake),indxcorners(xxx)-1,indycorners(xxx)+1);
    D = sub2ind(size(lake),indxcorners(xxx),indycorners(xxx)+1);
    DR = sub2ind(size(lake),indxcorners(xxx)+1,indycorners(xxx)+1);
    R = sub2ind(size(lake),indxcorners(xxx)+1,indycorners(xxx));
    UR = sub2ind(size(lake),indxcorners(xxx)+1,indycorners(xxx)-1);
    
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
            temporary = temporary + sum(dam(indshoreline == U));
        end
        if ULm
            temporary = temporary + sum(dam(indshoreline == UL));
        end
        if Lm
            temporary = temporary + sum(dam(indshoreline == L));
        end
        if DLm
            temporary = temporary + sum(dam(indshoreline == DL));
        end
        if Dm
            temporary = temporary + sum(dam(indshoreline == D));
        end
        if DRm
            temporary = temporary + sum(dam(indshoreline == DR));
        end
        if Rm
            temporary = temporary + sum(dam(indshoreline == R));
        end
        if URm
            temporary = temporary + sum(dam(indshoreline == UR));
        end
        
        damcorn(xxx) = temporary/num;
    else
        damcorn(xxx) = 0;
    end
end
end