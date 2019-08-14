function g = Source(p,g)

% Source terms or perturbations (Uplift/subsidence, etc.)

% surface uplift relative to boundaries
g = Uplift(p,g);

if p.n>1
    g.sealevel(p.n) = g.sealevel(p.n-1) + p.dt.*p.SLR; % will need to change this for sinusoidal
else
    g.sealevel(1) = p.sealevel_init;
end