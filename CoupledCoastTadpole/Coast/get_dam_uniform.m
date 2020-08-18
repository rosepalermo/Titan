function [dam_matrix, indshoreline] = get_dam_uniform(lake,p)
% initialize dam matrix
dam_matrix = zeros(size(lake));

% find the shoreline and shoreline indices
[shoreline] = addidshoreline(lake,~lake);
indshoreline = find(shoreline);

% uniform weighting
uniform_weight = double(shoreline);

%calc damage matrix for uniform erosion
dam = p.dt*p.Kuniform*uniform_weight(indshoreline)*p.So./p.dxo;
dam_matrix(indshoreline) = dam;

end