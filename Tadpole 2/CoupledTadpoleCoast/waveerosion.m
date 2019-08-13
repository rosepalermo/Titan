function [p g] = waveerosion(p,g)

% g.U is the elevation
g.wave_input = g.U<g.sealevel; % a meter above 0 I'm calling the coastline. after looking at a bunch of contour maps, it's close enough. Could even go to 2 m probably.


% call wave erosion function, output updates U and strength
[g.wave_output,g.strength] = coastal_erosion(g.wave_input,p.doWaveErosion,g.Strength);


% make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
% should be.
g.U(g.wave_output<=0) = g.sealevel-0.5;


% make all subaqueous points fixed points (not just new in case of sl
% fall)
p.F(g.U<g.sealevel) = 1;


end