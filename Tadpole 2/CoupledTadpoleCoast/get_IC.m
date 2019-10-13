function [init] = get_IC(p);

%inputs initial matrix p and outputs initial conditions for coupled tadpole
% and coastal erosion simulation

rng(0)

% add a random component
noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);
% noise = noise - mean(noise(:));
% noise = noise - min(noise(:));
% noise = noise - prctile(noise(:),10);

% sea = (initgaus + noise);
[H Wss] = Hann2D(noise);
init = H;

% A = diag(2*ones(100,1)) + diag(-1*ones(99,1),-1) + diag(-1*ones(99,1),1);
% plot(.1*mvnrnd(zeros(10,100), inv(A))')