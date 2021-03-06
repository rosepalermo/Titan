function [p g] = waveerosion(p,g)

% g.U is the elevation
g.wave_input = g.U<g.sealevel(p.n); % a meter above 0 I'm calling the coastline. after looking at a bunch of contour maps, it's close enough. Could even go to 2 m probably.


% call wave erosion function, output updates U and strength
[g.wave_output,g.strength,erodedind] = coastal_erosion(g.wave_input,1,g.Strength,p);

erodedind(1) = [];
% make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
% should be.
g.U(erodedind) = g.sealevel(p.n)-0.5;


% make all subaqueous points fixed points (not just new in case of sl
% fall)
p.F(g.U<g.sealevel(p.n)) = 1;


end