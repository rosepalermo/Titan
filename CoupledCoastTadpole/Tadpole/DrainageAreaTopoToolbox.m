function g = DrainageAreaTopoToolbox(p,g)

% pad to create periodic boundaries for TopoToolbox. We assume Nx and Ny
% are even throughout!
elevtile = repmat(g.U,[3,3]);
elevtile = elevtile(p.Ny/2+1:3*p.Ny-p.Ny/2,p.Nx/2+1:3*p.Nx-p.Nx/2);

snkstile = repmat(~g.C,[3,3]); % logical matrix with sinks = 1
snkstile = snkstile(p.Ny/2+1:3*p.Ny-p.Ny/2,p.Nx/2+1:3*p.Nx-p.Nx/2);


[nr,nc] = size(elevtile);

x = p.dx*(0:nc-1);
y = p.dy*(nr-1:-1:0);

elevtile = GRIDobj(x,y,elevtile);
snkstile = GRIDobj(x,y,snkstile); % make it a GRIDobj

FD = FLOWobj(elevtile,'preprocess','carve','sinks',snkstile,'mex',true); % flow direction object with steepest descent and carving through closed depressions that aren't sinks
A  = flowacc(FD); % flow accumulation in # of cells

Acells = A.Z(p.Ny/2+1:2*p.Ny-p.Ny/2,p.Nx/2+1:2*p.Nx-p.Nx/2); % Extract the original size grid

Acells = Acells.*g.C; % Assign drainage area of zero to any sinks

% multiply by cell area
g.A = p.dx*p.dy*Acells;

% We also need a grid indicating which points were flooded!
elevtilefs = fillsinks(elevtile,snkstile); % sink-filled GRIDobj

elevfs = elevtilefs.Z(p.Ny/2+1:2*p.Ny-p.Ny/2,p.Nx/2+1:2*p.Nx-p.Nx/2); % Extract the original size sink-filled grid for use below 

elevtilefs = elevtilefs - elevtile; % GRIDobj of sink depth
elevtilefs.Z(elevtilefs.Z > 0) = 1; % flooded areas are 1, non-flooded areas are 0

g.fl = elevtilefs.Z(p.Ny/2+1:2*p.Ny-p.Ny/2,p.Nx/2+1:2*p.Nx-p.Nx/2); % Extract the original size grid

% finally, we need the weights array for the upwind differencing. This
% could probably be constructed from TopoToolbox output, but this is fast
% enough for now.
doflood = 0;
[~, ~, ~, g.W] = mexD8(elevfs,p.dy/p.dx,1*(~g.C),p.bvec,doflood);
