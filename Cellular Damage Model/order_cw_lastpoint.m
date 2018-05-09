function [sl_ord_ccw] = order_cw_lastpoint(lake,shoreline);

% This code will start at a point on the shoreline, look counterclockwise
% at the 8-connected neighbors to find the next shoreline point, and
% continue to order the shoreline points in a counter clockwise LastPoint


%inputs
% lake - logical matrix that defines lake as true and not lake as false
% shoreline - output of addidshoreline_cardonly

[slXind,slYind] = find(shoreline);
slX_orderedind = nan(length(slXind),1);
slY_orderedind = nan(length(slYind),1);
land = ~lake;

% For Plotting
% edge = 2*shoreline;
% plot_sl = edge+land;
% imagesc(plot_sl');hold on;scatter(slXind(1),slYind(1))

% look around the first point to find the second
U = land(slXind(1),slYind(1)-1);
UL = land(slXind(1)-1,slYind(1)-1);
L = land(slXind(1)-1,slYind(1));
DL = land(slXind(1)-1,slYind(1)+1);
D = land(slXind(1),slYind(1)+1);
DR = land(slXind(1)+1,slYind(1)+1);
R = land(slXind(1)+1,slYind(1));
UR = land(slXind(1)+1,slYind(1)-1);

slX_orderedind(1) = slXind(1);
slY_orderedind(1) = slYind(1);

if (U - UL) < 0
    slX_orderedind(2) = slXind(1)-1;
    slY_orderedind(2) = slYind(1)-1;
    
elseif (UL - L) < 0
    slX_orderedind(2) = slXind(1)-1;
    slY_orderedind(2) = slYind(1);
    
elseif (L - DL) < 0
    slX_orderedind(2) = slXind(1)-1;
    slY_orderedind(2) = slYind(1)+1;
    
elseif (DL - D) < 0
    slX_orderedind(2) = slXind(1);
    slY_orderedind(2) = slYind(1)+1;
    
elseif (D - DR) < 0
    slX_orderedind(2) = slXind(1)+1;
    slY_orderedind(2) = slYind(1)+1;
    
elseif (DR - R) < 0
    slX_orderedind(2) = slXind(1)+1;
    slY_orderedind(2) = slYind(1);
    
elseif (R - UR) < 0
    slX_orderedind(2) = slXind(1)+1;
    slY_orderedind(2) = slYind(1)-1;
    
elseif (UR - U) < 0
    slX_orderedind(2) = slXind(1);
    slY_orderedind(2) = slYind(1)-1;
end


for i = 3:length(slX_orderedind)+1
    U = land(slX_orderedind(i-1),slY_orderedind(i-1)-1);
    UL = land(slX_orderedind(i-1)-1,slY_orderedind(i-1)-1);
    L = land(slX_orderedind(i-1)-1,slY_orderedind(i-1));
    DL = land(slX_orderedind(i-1)-1,slY_orderedind(i-1)+1);
    D = land(slX_orderedind(i-1),slY_orderedind(i-1)+1);
    DR = land(slX_orderedind(i-1)+1,slY_orderedind(i-1)+1);
    R = land(slX_orderedind(i-1)+1,slY_orderedind(i-1));
    UR = land(slX_orderedind(i-1)+1,slY_orderedind(i-1)-1);
    
    if [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)-1]
        LastPoint = 2;
    elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)]
        LastPoint = 3;
    elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)+1]
        LastPoint = 4;
    elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1),slY_orderedind(i-1)+1]
        LastPoint = 5;
    elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)+1]
        LastPoint = 6;
    elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)]
        LastPoint = 7;
    elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)-1]
        LastPoint = 8;
    elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1),slY_orderedind(i-1)-1]
        LastPoint = 1;
    end
    
    if LastPoint == 1
        if (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
        end
    end
    
    if LastPoint == 2
        if (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
        end
    end
    
    if LastPoint == 3
        if (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;        
            
        elseif (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
        end
    end
    
    if LastPoint == 4
        if (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
                    
        elseif (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
        end
    end
    
    if LastPoint == 5
        if (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
        end
    end
    
    if LastPoint == 6
        if (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        end
    end
    
    if LastPoint == 7
        if (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
        end
    end
    
    if LastPoint == 8
        if (UR - U) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (U - UL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            
        elseif (UL - L) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (L - DL) < 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DL - D) < 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (D - DR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            
        elseif (DR - R) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            
        elseif (R - UR) < 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
        end
    end
end

sl_ord_cw = [slX_orderedind,slY_orderedind];
sl_ord_ccw = flipud(sl_ord_cw);
