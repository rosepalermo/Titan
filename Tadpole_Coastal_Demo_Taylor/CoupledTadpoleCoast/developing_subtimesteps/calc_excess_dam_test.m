% calculate the excess damage

load('fetch_for_interp_test.mat')
x = 1:length(lake_save{1});
y = 1:length(lake_save{1});
[X,Y] = meshgrid(x,y);

for i = 1:length(strength_save)-1
    ind_excess_total{i} = find(strength_save{i+1}<0);
    dam_excess_total{i} = strength_save{i+1}(ind_excess_total{i});
    if i == 1
        ind_excess{i} = find(strength_save{i+1}<0);
    else
        ind_excess{i} = ind_excess_total{i}(~ismember(ind_excess_total{i},ind_excess_total{i-1}));
    end
    dam_excess{i} = strength_save{i+1}(ind_excess{i});
    sum_dam_excess(i) = sum(dam_excess{i});
    time_excess{i} = -dam_excess{i}./(fetch_save{i}(ind_excess{i})).*p.dt;
end

