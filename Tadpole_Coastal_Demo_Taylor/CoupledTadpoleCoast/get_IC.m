function [init] = get_IC(p,rfactor)

%inputs initial matrix p and outputs initial conditions for coupled tadpole
% and coastal erosion simulation

% Rose code:

% rng(0)
% 
% % add a random component
% noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);
% % noise = noise - mean(noise(:));
% % noise = noise - min(noise(:));
% % noise = noise - prctile(noise(:),10);
% 
% % sea = (initgaus + noise);
% [H Wss] = Hann2D(noise);
% init = H;
% 
% % A = diag(2*ones(100,1)) + diag(-1*ones(99,1),-1) + diag(-1*ones(99,1),1);
% % plot(.1*mvnrnd(zeros(10,100), inv(A))')

% Taylor code:

rng(0)

% create a random component
noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);

relief = max(noise(:)) - min(noise(:));

% create a depression with a desired depth (expressed as a multiple of the noise relief) 
depth = rfactor*relief;
depression = -1 * depth * Hann2D(ones(size(noise))); % this has a max of zero and a min of -depth

% add the depression to the noise to create the initial condition.
init = depression + noise;

% A = diag(2*ones(100,1)) + diag(-1*ones(99,1),-1) + diag(-1*ones(99,1),1);
% plot(.1*mvnrnd(zeros(10,100), inv(A))')