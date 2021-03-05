% prepare for analysis of many model runs
pmin = 512;
pmax = 1024;

load('idx_list_v1')
idx_list(33) = NaN; idx_list(46) = NaN; idx_list(isnan(idx_list)) = [];
%initial conditions
for runs = 1:length(idx_list)
file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_uniform_K_0_1.mat'];
load(file_name)
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,1,true);

[~,eq_14]= wavelets(x*p.dx,y*p.dx,wave_weight,pmin,pmax,['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_1'],1);

save_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_init_results.mat'];
save(save_name,'x','y','wave_weight','eq_14')
end

% model results
for runs = 1:length(idx_list)
file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_uniform_K_0_1.mat'];
load(file_name)
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,'end',true);

[~,eq_14]= wavelets(x*p.dx,y*p.dx,wave_weight,pmin,pmax,['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_end'],1);

save_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_uniform_K_0_1_end_results.mat'];
save(save_name,'x','y','wave_weight','eq_14')
end

for runs = 1:length(idx_list)
file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_wave_K_0_1.mat'];
load(file_name)
[x,y,wave_weight,xyz_name] = calc_xyw(file_name,'end',true);

[~,eq_14]= wavelets(x*p.dx,y*p.dx,wave_weight,pmin,pmax,['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_end'],1);

save_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_wave_K_0_1_end_results.mat'];
save(save_name,'x','y','wave_weight','eq_14')
end