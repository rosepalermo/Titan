function [p,g] = waveerosion(p,g)

% g.U is the elevation
% g.wave_input = g.U<g.sealevel(p.n); % a meter above 0 I'm calling the coastline. after looking at a bunch of contour maps, it's close enough. Could even go to 2 m probably.
g.wave_input = g.U<g.sealevel; % a meter above 0 I'm calling the coastline. after looking at a bunch of contour maps, it's close enough. Could even go to 2 m probably.
p.dt_save = p.dt;

if p.doAdaptiveCoastalTimeStep
    dt_temp = p.dt;
    i = 0;
    [adt,dam_matrix,wave_matrix,p] = getCADT(p,g,g.wave_input);
    p.dt = dt_temp./adt; %adaptive time step
    
    while i<=adt
        i = i+1;
        if i == 1
            [g.wave_output,g.Strength,erodedind] = coastal_erosion(g.wave_input,1,g.Strength,p,dam_matrix,wave_matrix);
        else
            [g.wave_output,g.Strength,erodedind] = coastal_erosion(g.wave_input,1,g.Strength,p,[],[]);
        end
        erodedind(1) = [];
        % make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
        % should be.
        % g.U(erodedind) = g.sealevel(p.n)-0.5;
        g.U(erodedind) = g.sealevel-0.5;
    end
    p.dt = dt_temp; %return the time step to what it was.
else
    
    % call wave erosion function, output updates U and strength
    [g.wave_output,g.Strength,erodedind,p] = coastal_erosion(g.wave_input,1,g.Strength,p,[]);
    
    erodedind(1) = [];
    % make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
    % should be.
    % g.U(erodedind) = g.sealevel(p.n)-0.5;
    g.U(erodedind) = g.sealevel-0.5;
end

if isfield(p,'boundary')
    return
end


end