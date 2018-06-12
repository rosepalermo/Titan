function [sl_cell] = order_cw_lastpoint(lake,shoreline)

% This code will start at a point on the shoreline, look counterclockwise
% at the 8-connected neighbors to find the next shoreline point, and
% continue to order the shoreline points in a counter clockwise LastPoint


%inputs
% lake - logical matrix that defines lake as true and not lake as false
% shoreline - output of addidshoreline_cardonly

[slXind,slYind] = find(shoreline);
sl_cell = cell(1,1);
land = ~lake;

% For Plotting
% edge = 2*shoreline;
% plot_sl = edge+land;
% imagesc(plot_sl');hold on;scatter(slXind(1),slYind(1))

% look around the first point to find the second and find the "Last Point"
% for point number 2
[slX_orderedind, slY_orderedind, LastPoint] = find2ndpt(land,slXind(1),slYind(1));


state = [slX_orderedind, slY_orderedind, LastPoint];
unq = unique(state,'rows');

ind = [slXind slYind];
ordered_cw = [slX_orderedind slY_orderedind];
ordered_cw_all = [slX_orderedind slY_orderedind];
i = 3;
obj = 1;

while sum(~ismember(ind,ordered_cw_all,'rows')) > 0
    if length(unq(:,1)) == length(state(:,1))
        U = land(slX_orderedind(i-1),slY_orderedind(i-1)-1);
        UL = land(slX_orderedind(i-1)-1,slY_orderedind(i-1)-1);
        L = land(slX_orderedind(i-1)-1,slY_orderedind(i-1));
        DL = land(slX_orderedind(i-1)-1,slY_orderedind(i-1)+1);
        D = land(slX_orderedind(i-1),slY_orderedind(i-1)+1);
        DR = land(slX_orderedind(i-1)+1,slY_orderedind(i-1)+1);
        R = land(slX_orderedind(i-1)+1,slY_orderedind(i-1));
        UR = land(slX_orderedind(i-1)+1,slY_orderedind(i-1)-1);
        
        if [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)-1]
            LastPoint(i) = 2;
        elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)]
            LastPoint(i) = 3;
        elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)+1]
            LastPoint(i) = 4;
        elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1),slY_orderedind(i-1)+1]
            LastPoint(i) = 5;
        elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)+1]
            LastPoint(i) = 6;
        elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)]
            LastPoint(i) = 7;
        elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)-1]
            LastPoint(i) = 8;
        elseif [slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1),slY_orderedind(i-1)-1]
            LastPoint(i) = 1;
        end
        
        if LastPoint(i) == 1
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
            
            
        elseif LastPoint(i) == 2
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
        
        elseif LastPoint(i) == 3
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
        
        elseif LastPoint(i) == 4
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
        
        elseif LastPoint(i) == 5
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
        
        elseif LastPoint(i) == 6
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
        
        elseif LastPoint(i) == 7
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
        
        elseif LastPoint(i) == 8
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
            else
                slX_orderedind(i) = slX_orderedind(i-1);
                slY_orderedind(i) = slY_orderedind(i-1);
            end
        end
        i = i+1;
        state = [slX_orderedind, slY_orderedind, LastPoint];
        unq = unique(state,'rows');
    else
        ordered_cw = [slX_orderedind slY_orderedind];
        ordered_cw_all = [ordered_cw_all;ordered_cw];
        ordered_ccw = flipud(ordered_cw);
        sl_cell{obj,1} = ordered_ccw;
        obj = obj+1;
        % start at point number 1 for the next object
        i = 2;
        ptsleft = ind(find(~ismember(ind,ordered_cw_all,'rows')),:);
        if ~isempty(ptsleft)
            clearvars slX_orderedind slY_orderedind LastPoint
            slX_orderedind(1) = ptsleft(1,1);
            slY_orderedind(1) = ptsleft(1,2);
            [slX_orderedind, slY_orderedind, LastPoint] = find2ndpt(land,slX_orderedind(1),slY_orderedind(1));
        end
        
        
    end
end

%sort the elements of the shoreline by the length in descending order
length_cells = cellfun(@length,sl_cell,'uni',false);
length_cells = cell2mat(length_cells);
[length_cells_sort,sortind] = sort(length_cells);
length_cells_sort = flipud(length_cells_sort);
sortind = flipud(sortind);
sl_cell = sl_cell(sortind,:);

end




