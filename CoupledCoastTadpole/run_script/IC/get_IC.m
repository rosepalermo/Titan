function [init,p] = get_IC(p,rfactor,idx)

%inputs initial matrix p and outputs initial conditions for coupled tadpole
% and coastal erosion simulation;
if isfield(p,'rand_gen')
    rng(idx)
end

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
% if p.Nx ==400
% rng(0)
% elseif p.Nx ==200
%     rng(1)
% end

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
depression = depression -Zshift + p.sealevel_init;
noise = noise-Zshift+p.sealevel_init;
    if p.Nx == 400
        p.Ao = 8.9298e+07;
        p.Ao_cells = 30368;
    else
        diff_ = min(min(depression))-min(min(noise)); % diff_o is 173 -- rng(0)
        Ao = depression<ceil(diff_); % diff_o is 173 -- rng(0)
        [shoreline] = addidshoreline(Ao,~Ao,p);
        indsl = find(shoreline);
        x = p.dx*(1:size(Ao,2));
        y = p.dy*(1:size(Ao,1));
        [X,Y] = meshgrid(x,y);
        radius = 0.5*(max(X(indsl))-min(X(indsl)));
        if isempty(radius)
            radius = 1;
        end
        p.Ao = 3/4*pi*radius^2;
        p.Ao_cells = length(find(Ao));
    end

end
