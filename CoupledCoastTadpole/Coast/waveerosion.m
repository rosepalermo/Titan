function [p,g] = waveerosion(p,g)

if isfield(p,'boundary')
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
            [g.wave_output,g.Strength,p,dam_matrix] = coastal_erosion(g.wave_input,1,g.Strength,p,dam_matrix,wave_matrix,cells2trash);
        else
            [g.wave_output,g.Strength,p,dam_matrix] = coastal_erosion(g.wave_input,1,g.Strength,p,[],[],[]);
        end
        erodedind = find(g.wave_output - g.wave_input);
        g.wave_input = g.wave_output; % update the input for the loop!!
        % make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
        % should be.
        % g.U(erodedind) = g.sealevel(p.n)-0.5;
        g.U(erodedind) = g.sealevel-0.5;
    end
    p.dt = p.dt_save; %return the time step to what it was.
else
    
    % call wave erosion function, output updates U and strength
    [g.wave_output,g.Strength,p,dam_matrix] = coastal_erosion(g.wave_input,1,g.Strength,p,[],[],[]);
    erodedind = find(g.wave_output - g.wave_input);
    % make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
    % should be.
    % g.U(erodedind) = g.sealevel(p.n)-0.5;
    g.U(erodedind) = g.sealevel-0.5;
end

% imagesc(g.wave_output)
g.nLakeCells = length(find(g.wave_output));
g.dam_wave = dam_matrix;
end