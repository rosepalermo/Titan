function [lake,strength,erodedind_save] = coastal_erosion(lake,fetch_on,strength,p)

% Titan analogue damage model for coastal erosion of a lake
% Rose Palermo 6-2019
% input is a matrix called "lake" -- this is a logical array of a landscape
% where 1 is lake and 0 is land. takes "lake" and either preforms uniform
% erosion on the shoreline or finds first order lakes


% initialize lake, land
x = 1:size(lake,2);
y = 1:size(lake,1);
[X,Y] = meshgrid(x,y);
land = ~lake;
% if fetch_on
% %     strength = 10000000*double(land);
%     strength = 50000*double(land);
% else
%     strength = 10*double(land);
% end

% initialize save variable
% shoreline_save = cell(1,1);


        
%% if 0, uniform erosion
    if ~fetch_on % if no fetch, the order doesn't matter and we can calc damage.
        [shoreline] = addidshoreline(lake,land); % corners are part of the shoreline!
        indshoreline = find(shoreline);
        dam = double(shoreline);
        %damage the shoreline
        strength(indshoreline) = strength(indshoreline)-dam(indshoreline);
        strength(strength<0) = 0;
        erodedind = indshoreline(strength(indshoreline)<=0);
        strength(erodedind) = 0;
        erodedX = X(erodedind);
        erodedY = Y(erodedind);
%         erodedi = cat(2,erodedX,erodedY);

if p.doStreamPower
    strength(erodedind) = p.strength;
end
        
        % change land to lake at eroded pts
        lake(erodedind) = true;
        land = ~lake;
        erodedind_save = erodedind;
%         [shoreline] = addidshoreline(lake,land);
%         shoreline_save{ff} = find(shoreline);
%         dam_save{ff} = dam;
%         lake_save{ff} = lake;
%         eroded = eroded(2:end,:);
    end
    
    %% if 1, wave erosion
    
    if fetch_on
        % if fetch_on 
       
        
        % loop over # of objects
        % find number of first order lakes
%         [L,fol] = find_first_order_lakes(lake_tile);
        [F_lake_all,~,~,~] = find_first_order_lakes(lake);
        % remove first order lakes that touch the boundaries
%         boundaries = [L(1,:)';L(end,:)';L(:,1);L(:,end)];
%         boundaries = boundaries(boundaries>0);
%         boundaries = unique(boundaries);
        erodedind_save = [NaN];
        for ff = 1:length(F_lake_all)
            F_lake = F_lake_all{ff};
            if length(find(F_lake))<2
                continue
            end
%             F_lake = (L == fol(ff));
% if (exist('erodedind','var')) | (ff==1)
            disp('fetch')
            clearvars fetch_sl_cells indshoreline WaveArea_cell
            
            %order the shoreline and islands
%             shoreline = addidshoreline_cardonly(lake,land); %rewrite shoreline to be card only for fetch.. will damage corners separately later
%             shoreline = addidshoreline_cardonly(F_lake,~F_lake); %rewrite shoreline to be card only for fetch.. will damage corners separately later            
            disp('ordering')
            % new order the shoreline code
            [indshoreline_ocw,~,cells2trash] = order_shoreline_bwbound(F_lake);
            if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
                break
            end
            % my original order the shoreline code
%             [indshoreline_ocw,cells2trash] = order_cw_lastpoint(F_lake,shoreline); % ccw ordered ind = indshoreline
%             keepme = 1:length(indshoreline_ocw);
            disp('ordered')
            for l = 1: length(indshoreline_ocw)
                indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
                fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
                fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
            end
            disp('calculating wave')
            % calculate wave weighted (sqrt(F)*cos(theta-phi))
            [~,WaveArea_cell] = fetch_vis_approx(fetch_sl_cells);% first is wave, second is fetch!!
            %             [WaveArea_cell,~] = fetch_wavefield_cell(fetch_sl_cells);
            %         [WaveArea_cell] = {ones(size(fetch_sl_cells{1,1},1),1)}; % ones to test debugging with
            disp('wave calculated')
            clearvars erodedind
%         end
        
        % go back to normal lake, not lake tile for damage?
        
        % Damage the shoreline
        indshoreline = cell2mat(indshoreline);
        %         [shoreline] = addidshoreline_cardonly(lake,land); % edges only
        [shoreline] = addidshoreline(F_lake,~F_lake); % corners and edges
        dam = cell2mat(WaveArea_cell);
%         SHOULDNT NEED TO DO THIS ANYMORE BECAUSE NEW ORDER THE SHORELINE
%         INCLUDES CORNERS
%         % find corners and damage if they exist alone (not in the cells to
%         % trash)
        corners = setdiff(find(shoreline),indshoreline);
        if ~isempty(cells2trash)
            [c2t]=sub2ind(size(X),cells2trash(:,1),cells2trash(:,2)); % find the cells that arent the trash cells
            if exist('c2t') & ~isempty(corners)
                corners = corners(~ismember(corners,c2t));
            end
        end
%         
%         % find the mean of the damage for points next to the corners. make that
%         % the damage for that corner
        if ~isempty(corners)
            corners = corners(shoreline(corners)<1.5); % if less than 1.5, only a corner. not also a side
            [damcorn] = damagecorners(lake,corners,indshoreline,dam);
            %         strength(indshoreline) = strength(indshoreline) - ones(length(indshoreline),1).*dam;
            strength(corners) = strength(corners) - shoreline(corners).*damcorn;
        end
        
        strength(strength<0) = 0; % if strength is negative, make it 0 for convenience

%   Find the corners and change the damage to sum corners* sqrt2/2 * wave
%   weighting
        [sl_nocorners] = addidshoreline_cardonly(F_lake,~F_lake);
        corners = setdiff(find(shoreline),find(sl_nocorners));
        cornind = ismember(indshoreline,corners);
%         dam(cornind) = dam(cornind).*shoreline(indshoreline(cornind));
        dam = dam.*shoreline(indshoreline);

        % DAMAGE THE COASTLINE
        strength(indshoreline) = strength(indshoreline) - shoreline(indshoreline).*dam;
        
        % find eroded points
        erodedind = indshoreline(strength(indshoreline)<=0);
        % erode points that weren't a corner and were only 2-1 cells
        % connected because it messed up the fetch calculations..
        if ~isempty(cells2trash)
            cells2trash = sub2ind(size(lake),cells2trash(:,1),cells2trash(:,2));
%             erodedind_12 = cells2trash(find(~ismember(cells2trash,corners)));
            erodedind =[erodedind;cells2trash];
        end
        
        
        strength(erodedind) = 0;
        if p.doStreamPower
            strength(erodedind) = p.strength;
        end
        erodedX = X(erodedind);
        erodedY = Y(erodedind);
        erodedi = cat(2,erodedX,erodedY);
%         eroded{ff} = cat(1,eroded,erodedi);
        
        % change land to lake at eroded pts
        lake(erodedind) = true;
        erodedind_save = [erodedind_save;erodedind];
%         land = ~lake;
        % update shoreline
%         [shoreline] = addidshoreline(lake,land);
%         ordered_sl_save{ff} = fetch_sl_cells;
%         corners_save{ff} = corners;
%         damcorners_save{ff} = damcorn;
%         eroded{ff} = eroded{ff}(2:end,:);
        end
    end
    
% lake_tile = lake_tile(ceil(size(lake_tile,2)/6):4*floor(size(lake_tile,1)/6),ceil(size(lake_tile,2)/6):4*floor(size(lake_tile,1)/6));
end