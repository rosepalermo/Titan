function [init] = get_IC(p,rfactor,idx)

%inputs initial matrix p and outputs initial conditions for coupled tadpole
% and coastal erosion simulation;
if exist('idx')
    rng(idx)
end
% create a random component
noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);

relief = max(noise(:)) - min(noise(:));

% create a depression with a desired depth (expressed as a multiple of the noise relief)
depth = rfactor*relief;
depression = -1 * depth * Hann2D(ones(size(noise))); % this has a max of zero and a min of -depth

% add the depression to the noise to create the initial condition.
init = depression + noise;

% adjust the elevations so pctwet % of the domain is below initial SL
pctwet = 10;
Zshift = prctile(init(:),pctwet);
init = init - Zshift + p.sealevel_init;

end
