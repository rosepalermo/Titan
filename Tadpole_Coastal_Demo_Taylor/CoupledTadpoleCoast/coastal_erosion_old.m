function [lake,strength,erodedind_save,p] = coastal_erosion_old(lake,fetch_on,strength,p,dam_matrix,wave_matrix)

% Titan analogue damage model for coastal erosion of a lake
% Rose Palermo 6-2019
% input is a matrix called "lake" -- this is a logical array of a landscape
% where 1 is lake and 0 is land. takes "lake" and either preforms uniform
% erosion on the shoreline or finds first order lakes


% initialize lake, land
x = 1:size(lake,2);
y = 1:size(lake,1);
[X,Y] = meshgrid(x,y);
X = X*p.dx; Y = Y*p.dy;
land = ~lake;
if length(dam_matrix)==0
    dam_matrix = zeros(size(lake));
end

%% if 0, uniform erosion
if ~fetch_on % if no fetch, the order doesn't matter and we can calc damage.
    [shoreline] = addidshoreline(lake,land); % corners are part of the shoreline!
    indshoreline = find(shoreline);
    uniform_weight = double(shoreline);
    
    
    dam = p.dt*p.Kcoast*uniform_weight(indshoreline);
    dam_matrix(indshoreline) = dam;
    
    % damage the shoreline
    ind_sl_old = find(shoreline);
    strength = strength - dam_matrix;
    lake(strength<=0) = 1;
    land = ~lake;
    
    %update the lake & shoreline & find new sl indices
    [shoreline_new] = addidshoreline(lake,land); % corners and edges
    ind_sl_new = find(shoreline_new);
    uniform_weight_new = double(shoreline_new);
    ind_new = ind_sl_new(~ismember(ind_sl_new,ind_sl_old));
    
    % find the points that failed during this time step and erode the
    % neighbors for the amount of time since the coastline points failed
    [sum_dam_excess,ind_excess,time_excess] = calc_excess_dam(strength,dam_matrix,p);
    %     lake(ind_excess) = 1;
    while sum_dam_excess > 0
        % make strength 0 for excess that I'm about to distribute
        strength(find(strength<0)) = 0;
        % damage the neighbors for the amount of time after the cell failed

        % find max time excess of 8 connected neighbors
        time_dam_excess = find_max_excess_time(ind_excess,time_excess,ind_new,lake);
        % damage shoreline for max excess time of 8 con neighbor
        strength(ind_new) = strength(ind_new) - time_dam_excess.*p.Kcoast.*shoreline_new(ind_new);
        
        % UPDATE LAKE
        lake(strength<=0) = 1;
        land = ~lake;
        
        % UPDATE SHORELINE
        % this will be where I take damage the neighbors of cells with
        % excess damage
        [shoreline_new] = addidshoreline(lake,land); % corners and edges
        ind_sl_new = find(shoreline_new);
        ind_new = ind_sl_new(~ismember(ind_sl_new,ind_sl_old));
        
        % recalculate sum of excess damage (if the neighboring cells erode,
        % then the time will need to pass on to their neighbors in the next loop)
        % end
        [sum_dam_excess,ind_excess,time_excess] = calc_excess_dam(strength);
    end
        
    % THIS IS PROBABLY WRONG-- MISSING POINTS ERODED DURING REDISTRIBUTION
    % OF EXCESS TIME
    erodedind = indshoreline(strength(indshoreline)==0);
    
    if p.doStreamPower % if doing streampower, return strength to IC to represent strength of underlying rock.
        strength(erodedind) = p.strength;
    end
    
    % change land to lake at eroded pts
    lake(erodedind) = true;
    erodedind_save = erodedind;
end

%% if 1, wave erosion

if fetch_on
    erodedind_save = [NaN];
    if sum(dam_matrix,'all')==0
        % calculate dam_matrix unless it already exists
        % find number of first order lakes
        [F_lake_all,~,~,~] = find_first_order_lakes(lake);
        wave_weight_matrix = zeros(size(lake));
        for ff = 1:length(F_lake_all)
            F_lake = F_lake_all{ff};
            if length(find(F_lake))<2
                continue
            end
            clearvars fetch_sl_cells indshoreline WaveArea_cell
            
            %order the shoreline and islands
            disp('ordering')
            % new order the shoreline code
            [indshoreline_ocw,~,cells2trash] = order_shoreline_bwbound(F_lake);
            if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
                break
            end
            % my original order the shoreline code
            disp('ordered')
            for l = 1: length(indshoreline_ocw)
                indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
                fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
                fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
            end
            disp('calculating wave')
            % calculate wave weighted (sqrt(F)*cos(theta-phi))
            [WaveArea_cell,~] = fetch_vis_approx(fetch_sl_cells);% first is wave, second is fetch!!
            disp('wave calculated')
            clearvars erodedind
            %         end
            
            % Damage the shoreline
            indshoreline = cell2mat(indshoreline);
            %         [shoreline] = addidshoreline_cardonly(lake,land); % edges only
            [shoreline] = addidshoreline(lake,land); % corners and edges
            wave_weighting = cell2mat(WaveArea_cell);
            wave_weight_matrix(indshoreline) = wave_weighting;
            dam = p.dt*p.Kcoast*shoreline(indshoreline).*wave_weighting;
            dam_matrix(indshoreline) = dam;
            
            if max(dam)>1  %if the maximum damage is greater than 1, we need the adaptive coastal time step. grid bias will be introduced otherwise.
                p.doAdaptiveCoastalTimeStep = 1;
                max(dam)
            end
        end
    end
    if isfield(p,'boundary')
        erodedind_save = [];
        return
    end
    
    % damage shoreline
    [shoreline] = addidshoreline(lake,land);
    strength = strength - dam_matrix;
    ind_sl_old = find(shoreline);
    lake(strength<=0) = 1;
    land = ~lake;
    
    % this will be where I take damage the neighbors of cells with
    % excess damage
    [shoreline_new] = addidshoreline(lake,land); % corners and edges
    ind_sl_new = find(shoreline_new);
    ind_new = ind_sl_new(~ismember(ind_sl_new,ind_sl_old));

    
    % calculate sum of excess damage
    [sum_dam_excess,ind_excess,time_excess] = calc_excess_dam(strength,dam_matrix,p);
    while sum_dam_excess > 0
        % make strength 0 for excess that I'm about to distribute
        strength(find(strength<0)) = 0;
        % damage the neighbors for the amount of time after the cell failed
        % interpolate wave_weighting for the neighboring cells
        wave_weighting_interp = interp_fetch_for_ind(lake,ind_new,wave_matrix);
        % find max time excess of 8 connected neighbors
        time_dam_excess = find_max_excess_time(ind_excess,time_excess,ind_new,lake);
        % damage shoreline for max excess time of 8 con neighbor
        strength(ind_new) = strength(ind_new) - time_dam_excess.*p.Kcoast.*shoreline_new(ind_new).*wave_weighting_interp;
        
        % UPDATE LAKE
        lake(strength<=0) = 1;
        land = ~lake;
        
        % UPDATE SHORELINE
        % this will be where I take damage the neighbors of cells with
        % excess damage
        [shoreline_new] = addidshoreline(lake,land); % corners and edges
        ind_sl_new = find(shoreline_new);
        ind_new = ind_sl_new(~ismember(ind_sl_new,ind_sl_old));
        
        % recalculate sum of excess damage (if the neighboring cells erode,
        % then the time will need to pass on to their neighbors in the next loop)
        % end
        [sum_dam_excess,ind_excess,time_excess] = calc_excess_dam(strength);
    end
%     strength(strength<0) = 0; % if strength is negative, make it 0 for convenience
    % find eroded points
    erodedind = find(strength<=0);
    % erode points that were only 2-1 cells
    % connected because it messed up the fetch calculations..
    if exist('cells2trash') && ~isempty(cells2trash)
        cells2trash = sub2ind(size(lake),cells2trash(:,1),cells2trash(:,2));
        erodedind =[erodedind;cells2trash];
        strength(cells2trash) = 0;
    end
    
    if p.doStreamPower
        strength(erodedind) = p.strength;
    end
    
    % change land to lake at eroded pts
    lake(erodedind) = true;
    erodedind_save = [erodedind_save;erodedind];
end


end