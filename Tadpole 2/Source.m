function g = Source(p,g)

% Source terms or perturbations (Uplift/subsidence, etc.)

% surface uplift relative to boundaries
g = Uplift(p,g);

if n>1
    g.sealevel(n) = g.sealevel(n-1) + p.dt.*p.SLR;
else
    g.sealevel(1) = g.sealevel_init;
end