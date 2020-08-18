% fetch interpolation for specified points
function [wave_interp] = interp_fetch_for_ind(lake,ind_calc,wave_matrix)

% make X and Y grid for scattered interpolation
x = 1:length(lake);
y = 1:length(lake);
[X,Y] = meshgrid(x,y);

[xind,yind] = ind2sub(size(lake),ind_calc);

F_natural = scatteredInterpolant(X(~isnan(wave_matrix)),Y(~isnan(wave_matrix)),wave_matrix(~isnan(wave_matrix)),'natural');
wave_interp = F_natural(xind,yind);

end
