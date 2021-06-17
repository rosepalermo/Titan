function [p,g] = waveerosion(p,g)

if isfield(p,'boundary')
    g.nLakeCells = length(find(g.U<g.sealevel));
    g.dam_wave = zeros(size(g.U));
    g.wave_matrix = zeros(size(g.U));
    g.adt = 0;
    return
end

% g.U is the elevation
% g.wave_input = g.U<g.sealevel(p.n); % a meter above 0 I'm calling the coastline. 
g.wave_input = g.U<g.sealevel; % a meter above 0 I'm calling the coastline. 
p.dt_save = p.dt;
g.nLakeCells = length(find(g.wave_input));

if p.doAdaptiveCoastalTimeStep
    i = 0;
    [adt,dam_matrix,wave_matrix,cells2trash,p] = getCADT_wave(p,g.wave_input);
    p.dt = p.dt_save./adt; %adaptive time step
    
%     if p.t> 1100
%         disp('testing')
%     end
    
    while i<adt
        i = i+1;    
        if i == 1
            [g.wave_output,g.Strength,p,dam_matrix,wave_matrix] = coastal_erosion(g.wave_input,1,g.Strength,p,dam_matrix,wave_matrix,cells2trash);
        elseif p.doAdaptiveCoastalTimeStep_interpwave
            [g.wave_output,g.Strength,p,dam_matrix,wave_matrix] = coastal_erosion(g.wave_input,1,g.Strength,p,dam_matrix,wave_matrix,[]);
        else
            [g.wave_output,g.Strength,p,dam_matrix,wave_matrix] = coastal_erosion(g.wave_input,1,g.Strength,p,[],[],[]);
        end
        erodedind = find(g.wave_output - g.wave_input);
        
        % update interpolated wave matrix
        [shoreline_old] = addidshoreline(g.wave_input,~g.wave_input,p);
        ind_sl_old = find(shoreline_old);
        [shoreline_new] = addidshoreline(g.wave_input,~g.wave_input,p);
        ind_sl_new = find(shoreline_new);
        ind_interp  = ind_sl_new(~ismember(ind_sl_new,ind_sl_old));
        wave_weighting_interp = interp_fetch_for_ind(g.wave_output,ind_interp,wave_matrix); % this already includes
        wave_matrix(ind_interp) = wave_weighting_interp;
        dam_matrix(ind_interp) = p.dt*p.Kwave*shoreline_new(ind_interp).*wave_matrix(ind_interp)*p.So./p.dxo./p.Ao;

        g.wave_input = g.wave_output; % update the input for the loop!!
        % g.U(erodedind) = g.sealevel(p.n)-0.5;
        g.U(erodedind) = g.sealevel-0.5;
        
    end
    p.dt = p.dt_save; %return the time step to what it was.
else
    
    % call wave erosion function, output updates U and strength
    [g.wave_output,g.Strength,p,dam_matrix,wave_matrix] = coastal_erosion(g.wave_input,1,g.Strength,p,[],[],[]);
    erodedind = find(g.wave_output - g.wave_input);
    % make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
    % should be.
    % g.U(erodedind) = g.sealevel(p.n)-p.deptheroded;
    g.U(erodedind) = g.sealevel-p.deptheroded;
end

% imagesc(g.wave_output)
g.nLakeCells = length(find(g.wave_output));
g.dam_wave = dam_matrix;
g.wave_matrix = wave_matrix;
g.adt = adt;
end