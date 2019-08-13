function [p g] = uniformerosion(p,g)

if n>1
% g.U is the elevation
g.Coast_input = g.U<g.sealevel; % a meter above 0 I'm calling the coastline. after looking at a bunch of contour maps, it's close enough. Could even go to 2 m probably.
end

% call wave erosion function, output updates U and strength
[g.uniform_output,g.Strength] = coastal_erosion(g.Coast_input,p.doUniformErosion,g.Strength);


% make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
% should be.
g.U(g.uniform_outpu<=0) = g.sealevel-0.5;


% make all subaqueous points fixed points (not just new in case of sl
% fall)
p.F(g.U<g.sealevel) = 1;


end