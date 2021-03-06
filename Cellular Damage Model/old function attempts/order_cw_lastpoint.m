function [sl_cell,cells2trash] = order_cw_lastpoint(lake,shoreline)

% This code will start at a point on the shoreline, look counterclockwise
% at the 8-connected neighbors to find the next shoreline point, and
% continue to order the shoreline points in counter clockwise direction.
% Flips so clockwise at the end. Last point is the starting point for look
% direction.
% THIS CODE DOES NOT FIND CORNERS. ONLY CARDINAL DIRECTION.

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
[slX_orderedind, slY_orderedind, LastPoint,LP_state] = find2ndpt(land,slXind(1),slYind(1));


state = [slX_orderedind, slY_orderedind, LP_state];
state_all = [state];
unq = unique(state,'rows');
unq_all = unq;

ind = [slXind slYind];
ordered_cw = [slX_orderedind slY_orderedind];
ordered_cw_all = [slX_orderedind slY_orderedind];
i = 3;
obj = 1;

% imagesc(lake');hold on;

while sum(~ismember(ind,ordered_cw_all,'rows')) > 0

    if length(unq(:,1)) == length(state(:,1)) && ~all([slX_orderedind(1), slY_orderedind(1)] == [slX_orderedind(2), slY_orderedind(2)])
%         scatter(slX_orderedind,slY_orderedind);
        U = land(slX_orderedind(i-1),slY_orderedind(i-1)-1);
        UL = land(slX_orderedind(i-1)-1,slY_orderedind(i-1)-1);
        L = land(slX_orderedind(i-1)-1,slY_orderedind(i-1));
        DL = land(slX_orderedind(i-1)-1,slY_orderedind(i-1)+1);
        D = land(slX_orderedind(i-1),slY_orderedind(i-1)+1);
        DR = land(slX_orderedind(i-1)+1,slY_orderedind(i-1)+1);
        R = land(slX_orderedind(i-1)+1,slY_orderedind(i-1));
        UR = land(slX_orderedind(i-1)+1,slY_orderedind(i-1)-1);
        
        if all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)-1])
            LastPoint(i) = 2;
            LP_state(i,:) = [slX_orderedind(i-1)-1,slY_orderedind(i-1)-1];
        elseif all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)])
            LastPoint(i) = 3;
            LP_state(i,:) = [slX_orderedind(i-1)-1,slY_orderedind(i-1)];
        elseif all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)-1,slY_orderedind(i-1)+1])
            LastPoint(i) = 4;
            LP_state(i,:) = [slX_orderedind(i-1)-1,slY_orderedind(i-1)+1];
        elseif all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1),slY_orderedind(i-1)+1])
            LastPoint(i) = 5;
            LP_state(i,:) = [slX_orderedind(i-1),slY_orderedind(i-1)+1];
        elseif all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)+1])
            LastPoint(i) = 6;
            LP_state(i,:) = [slX_orderedind(i-1)+1,slY_orderedind(i-1)+1];
        elseif all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)])
            LastPoint(i) = 7;
            LP_state(i,:) = [slX_orderedind(i-1)+1,slY_orderedind(i-1)];
        elseif all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1)+1,slY_orderedind(i-1)-1])
            LastPoint(i) = 8;
            LP_state(i,:) = [slX_orderedind(i-1)+1,slY_orderedind(i-1)-1];
        elseif all([slX_orderedind(i-2),slY_orderedind(i-2)] == [slX_orderedind(i-1),slY_orderedind(i-1)-1])
            LastPoint(i) = 1;
            LP_state(i,:) = [slX_orderedind(i-1),slY_orderedind(i-1)-1];
        else
            break
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
        state = [slX_orderedind, slY_orderedind, LP_state];
        unq = unique(state,'rows');
    else
%                 scatter(slX_orderedind,slY_orderedind);

        ordered_cw = [slX_orderedind, slY_orderedind];
        if obj == 1
            ordered_cw_all = ordered_cw;
            state_all = state;
            unq_all = unq;
        else
            ordered_cw_all = [ordered_cw_all;ordered_cw];
            state_all = [state_all;state];
            unq_all = [unq_all;unq];

        end
        ordered_ccw = flipud(ordered_cw);
        sl_cell{obj,1} = ordered_ccw(1:end-3,:);
        obj = obj+1;
        i = 3;
        % start at point number 1 for the next object
        ptsleft = ind(find(~ismember(ind,ordered_cw_all,'rows')),:);
        if ~isempty(ptsleft)
            clearvars slX_orderedind slY_orderedind LastPoint state unq
            slX_orderedind(1) = ptsleft(1,1);
            slY_orderedind(1) = ptsleft(1,2);
            [slX_orderedind, slY_orderedind, LastPoint, LP_state] = find2ndpt(land,slX_orderedind(1),slY_orderedind(1));
            state = [slX_orderedind, slY_orderedind, LP_state];
            unq = unique(state,'rows');
        end
        
        
    end
end

% get rid of units 2 or less because they mess up in the fetch calc.
length_cells = cellfun(@length,sl_cell,'uni',false);
length_cells = cell2mat(length_cells);
keepme = (length_cells) > 3;
cells2trash = cell2mat(sl_cell(~keepme));
sl_cell = sl_cell(keepme);

%sort the elements of the shoreline by the length in descending order
length_cells = cellfun(@length,sl_cell,'uni',false);
length_cells = cell2mat(length_cells);
[length_cells_sort,sortind] = sort(length_cells);
length_cells_sort = flipud(length_cells_sort);
sortind = flipud(sortind);
sl_cell = sl_cell(sortind,:);

end




