function [dam_matrix, indshoreline, p] = get_dam_uniform(lake,p)
% initialize dam matrix
dam_matrix = zeros(size(lake));

% find the shoreline and shoreline indices
[shoreline] = addidshoreline(lake,~lake,p);
indshoreline = find(shoreline);

% uniform weighting
uniform_weight = double(shoreline);

%calc damage matrix for uniform erosion
dam = p.dt*p.Kuniform*uniform_weight(indshoreline)*p.So./p.dxo*p.dx;

dam_matrix(indshoreline) = dam;

%% check if shoreline is on a boundary
[p] = check_boundary(lake,p);

end