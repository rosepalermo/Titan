function [p,g] = SeaLevel(p,g)

% if p.n>1
%     g.sealevel(p.n) = g.sealevel(p.n-1) + p.dt*p.SLR; % will need to change this for sinusoidal
% else
%     g.sealevel(1) = p.sealevel_init;
% end

% g.sealevel(p.n) = p.sealevel_all(p.n);

% Case 1: Constant rate of SL change
% g.sealevel = g.sealevel + p.dt*p.SLR;

% Case 2: time-dependent function

SLperiod = 2e4; % sea level period (yr)
SLamp = 100; % sea level amplitude (m)
g.sealevel = p.sealevel_init + SLamp * sin(2*pi*p.t/SLperiod);
