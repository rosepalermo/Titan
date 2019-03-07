function [p g] = input_coastalerosion(p,g)

% g.U is the elevation
g.Coast_input = g.U<g.sealevel; % a meter above 0 I'm calling the coastline. after looking at a bunch of contour maps, it's close enough. Could even go to 2 m probably.


% need to do something fancy here to be able to do coastal erosion on
% multiple lakes...

% call coastal erosion function here
% output from coastal erosion should be g.Coast_output
% g.Coast_output should be lake in my coastal erosion code.
% take the eroded points as an output also.
% also output the new strength matrix as g.Strength


% make new subaqueous points elevation 0.5? Talk to Andrew about what this depth
% should be.
g.U(coastal_erodedpoints) = 0.5;

% make all subaqueous points fixed points (not just new in case of sl
% fall)
p.F(g.U<g.sealevel) = 1;


end