function [init] = get_IC(p);

%inputs initial matrix p and outputs initial conditions for coupled tadpole
% and coastal erosion simulation



% add a random component
noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);
noise = noise - mean(noise(:));

% sea = (initgaus + noise);
[H Wss] = Hann2D(noise);
init = H + noise;
