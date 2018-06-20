function [slX_orderedind, slY_orderedind, LastPoint,LP_state] = find2ndpt(land,slXind,slYind)

slX_orderedind = nan(2,1);
slY_orderedind = nan(2,1);


% find 2nd shoreline point
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
else
    slX_orderedind(2) = slXind(1);
    slY_orderedind(2) = slYind(1);
end

% find LastPoint(2)
if [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1)-1,slY_orderedind(1)-1]
    LastPoint(2,1) = 2;
    LP_state(2,:) = [slX_orderedind(1)-1,slY_orderedind(1)-1];
elseif [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1)-1,slY_orderedind(1)]
    LastPoint(2,1) = 3;
    LP_state(2,:) = [slX_orderedind(1)-1,slY_orderedind(1)];
elseif [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1)-1,slY_orderedind(1)+1]
    LastPoint(2,1) = 4;
    LP_state(2,:) = [slX_orderedind(1)-1,slY_orderedind(1)+1];
elseif [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1),slY_orderedind(1)+1]
    LastPoint(2,1) = 5;
    LP_state(2,:) = [slX_orderedind(1),slY_orderedind(1)+1];
elseif [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1)+1,slY_orderedind(1)+1]
    LastPoint(2,1) = 6;
    LP_state(2,:) = [slX_orderedind(1)+1,slY_orderedind(1)+1];
elseif [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1)+1,slY_orderedind(1)]
    LastPoint(2,1) = 7;
    LP_state(2,:) = [slX_orderedind(1)+1,slY_orderedind(1)];
elseif [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1)+1,slY_orderedind(1)-1]
    LastPoint(2,1) = 8;
    LP_state(2,:) = [slX_orderedind(1)+1,slY_orderedind(1)-1];
elseif [slX_orderedind(2),slY_orderedind(2)] == [slX_orderedind(1),slY_orderedind(1)-1]
    LastPoint(2,1) = 1;
    LP_state(2,:) = [slX_orderedind(1),slY_orderedind(1)-1];
else
    LastPoint = zeros(2,1);
    LP_state(2,:) = zeros(2,1);
end



end