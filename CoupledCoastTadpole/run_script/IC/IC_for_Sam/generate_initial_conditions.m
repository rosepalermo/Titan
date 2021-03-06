% generate initial condition

p.Ny = 200; %Number of x grid points
p.Nx = 200; %Number of y grid points
p.beta = 1.6; %Negative slope of the power spectrum. 0 = white noise, more positive values are "redder" (more variance at longer wavelengths)
p.variance = 10000; %Variance of elevation (m^2)
p.periodic = 1; %Elevations will be periodic at the boundaries (1) or not (0, default)
p.sealevel_init = 1; % Initial sea level

idx = 1; % random seed - this is an optional input
rfactor = 0.25;
[init] = get_IC(p,rfactor,idx);

