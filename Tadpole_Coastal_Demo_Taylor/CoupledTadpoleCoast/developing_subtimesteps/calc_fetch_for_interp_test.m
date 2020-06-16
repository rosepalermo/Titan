% test interp vs nearest neighbor subtimesteps
% this script will calculate the wave area for two time steps that will be
% used in the interp vs nearest neighbor comparison

% load data
runname = 'uniform_river_tf1e5_Kf_2e-05';
addpath('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/ModelingAGU19');
addpath(genpath('/Users/rosepalermo/Documents/GitHub/Titan/Tadpole_Coastal_Demo_Taylor'));
load(runname)

p.doStreamPower = 0;

x = p.dx*(0:p.Nx-1);
y = p.dy*(p.Ny-1:-1:0);

p.Kcoast = 1e-10;

% find the lake
elev = output(:,:,5);
sl = sealevel(5)/5;
time = t(5);
xx = 1:length(elev);
yy = xx;
[X,Y] = meshgrid(xx,yy);

lake = elev<sl;
land = ~lake;
strength = double(land);
lake_save{1} = lake;
strength_save{1} = strength;


erodedind_save = {};
erodi = 0;
for i = 2:20
    
    fetch_matrix = nan(size(lake));
    
    
    % calculate fetch at ts1, erode, then calc fetch at ts2    
        [F_lake_all,~,~,~] = find_first_order_lakes(lake);
        erodedind_save = [NaN];
        for ff = 1:length(F_lake_all)
            F_lake = F_lake_all{ff};
            if length(find(F_lake))<2
                continue
            end
            clearvars fetch_sl_cells indshoreline WaveArea_cell
            
            %order the shoreline and islands
            % new order the shoreline code
            [indshoreline_ocw,~,cells2trash] = order_shoreline_bwbound(F_lake);
            if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
                break
            end
            
            for l = 1: length(indshoreline_ocw)
                indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
                fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
                fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
            end
            
            % calculate wave weighted (sqrt(F)*cos(theta-phi))
            [WaveArea_cell,~] = fetch_vis_approx(fetch_sl_cells);% first is wave, second is fetch!!
            clearvars erodedind
                
        % Damage the shoreline
        indshoreline = cell2mat(indshoreline);
        %         [shoreline] = addidshoreline_cardonly(lake,land); % edges only
        [shoreline] = addidshoreline(lake,land); % corners and edges
        dam = cell2mat(WaveArea_cell);
        strength(indshoreline) = strength(indshoreline) - p.dt*p.Kcoast*shoreline(indshoreline).*dam; % Taylor's modified line that depends on a rate constant
        adt = p.dt*p.Kcoast*shoreline(indshoreline).*dam;
        fetch_matrix(indshoreline) = adt;

        if max(adt)>1  %if the maximum damage is greater than 1, we need the adaptive coastal time step. grid bias will be introduced otherwise.
            p.doAdaptiveCoastalTimeStep = 1;
            max(adt)
        end
        
%         strength(strength<0) = 0; % if strength is negative, make it 0 for convenience
%         
        % find eroded points
        erodedind = indshoreline(strength(indshoreline)<=0);
        % erode points that weren't a corner and were only 2-1 cells
        % connected because it messed up the fetch calculations..
        if ~isempty(cells2trash)
            cells2trash = sub2ind(size(lake),cells2trash(:,1),cells2trash(:,2));
%             erodedind_12 = cells2trash(find(~ismember(cells2trash,corners)));
            erodedind =[erodedind;cells2trash];
            strength(cells2trash) = 0;
        end
        
        
%         strength(erodedind) = 0;
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
        
        
        end
    lake_save{i} = lake;
    strength_save{i} = strength;
    fetch_save{i-1} = fetch_matrix;
    clearvars fetch_matrix
end

clearvars erodedind
for ii = 2:length(lake_save);
    erodedind{ii-1} = setdiff(find(lake_save{ii}),find(lake_save{ii-1}));
end

for ii = 2:length(fetch_save);
    shorelinediff{ii-1} = setdiff(find(~isnan(fetch_save{ii})),find(~isnan(fetch_save{ii-1})));
end

save('fetch_for_interp_test.mat','lake_save','fetch_save','erodedind','shorelinediff','strength_save','p')