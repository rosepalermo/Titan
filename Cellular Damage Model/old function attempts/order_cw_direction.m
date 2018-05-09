% clear all
% load('test_polygon')
polygon = lake;
[edge] = addidshoreline_cardonly(polygon,~polygon);

% This code will start at a point on the shoreline, look counterclockwise
% at the 8-connected neighbors to find the next shoreline point, and
% continue to order the shoreline points in a counter clockwise direction
slind = find(edge);
[slXind,slYind] = find(edge);
slX_orderedind = nan(length(slXind),1);
slY_orderedind = nan(length(slYind),1);
lake = polygon;
land = ~lake;
shoreline = 2*edge;
plot_sl = shoreline+land;
imagesc(plot_sl');hold on;scatter(slXind(1),slYind(1))

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
direction = 0;
if (U - UL) > 0
    slX_orderedind(2) = slXind(1)-1;
    slY_orderedind(2) = slYind(1)-1;
    direction = 1;
elseif (UL - L) > 0
    slX_orderedind(2) = slXind(1)-1;
    slY_orderedind(2) = slYind(1);
    direction = 2;
elseif (L - DL) > 0
    slX_orderedind(2) = slXind(1)-1;
    slY_orderedind(2) = slYind(1)+1;
    direction = 3;
    
elseif (DL - D) > 0
    slX_orderedind(2) = slXind(1);
    slY_orderedind(2) = slYind(1)+1;
    direction = 4;
    
    
elseif (D - DR) > 0
    slX_orderedind(2) = slXind(1)+1;
    slY_orderedind(2) = slYind(1)+1;
    direction =5;
    
elseif (DR - R) > 0
    slX_orderedind(2) = slXind(1)+1;
    slY_orderedind(2) = slYind(1);
    direction = 6;
    
    
elseif (R - UR) > 0
    slX_orderedind(2) = slXind(1)+1;
    slY_orderedind(2) = slYind(1)-1;
    direction = 7;
    
elseif (UR - U) > 0
    slX_orderedind(2) = slXind(1);
    slY_orderedind(2) = slYind(1)-1;
    direction = 8;
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
    if direction == 1
        if (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        elseif (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
            
            
        elseif (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
            
        elseif (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
            
            
        elseif (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
            
        elseif (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
            
            
        elseif (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
            
        elseif (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        end
    end
    if direction == 2
        
        if (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
            
            
        elseif (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
            
        elseif (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
            
            
        elseif (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
            
        elseif (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
            
            
        elseif (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
            
        elseif (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        elseif (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        end
    end
    if direction == 3
        if (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
            
        elseif (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
            
            
        elseif (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
            
        elseif (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
            
            
        elseif (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
            
        elseif (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        elseif (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        elseif (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
        end
    end
    if direction == 4
        if (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
            
            
        elseif (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
            
        elseif (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
            
            
        elseif (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
            
        elseif (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        elseif (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        elseif (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
        elseif (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
        end
    end
    if direction == 5
        if (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
            
        elseif (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
            
            
        elseif (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
            
        elseif (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        elseif (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        elseif (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
        elseif (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
        elseif (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
        end
    end
    
    if direction == 6
        if (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
            
            
        elseif (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
            
        elseif (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        elseif (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        elseif (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
        elseif (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
        elseif (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
        elseif (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
        end
    end
    
    if direction == 7
        if (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
            
        elseif (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        elseif (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        elseif (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
        elseif (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
        elseif (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
        elseif (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
        elseif (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
        end
    end
    
    if direction == 8
        if (UR - U) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 8;
        elseif (U - UL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 1;
        elseif (UL - L) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 2;
        elseif (L - DL) > 0
            slX_orderedind(i) = slX_orderedind(i-1)-1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 3;
        elseif (DL - D) > 0
            slX_orderedind(i) = slX_orderedind(i-1);
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction = 4;
        elseif (D - DR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)+1;
            direction =5;
        elseif (DR - R) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1);
            direction = 6;
        elseif (R - UR) > 0
            slX_orderedind(i) = slX_orderedind(i-1)+1;
            slY_orderedind(i) = slY_orderedind(i-1)-1;
            direction = 7;
        end
    end
    
    
    
end
