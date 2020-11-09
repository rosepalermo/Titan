% calculate x,y and fetch

function [x,y,wave_weight] = calc_xyw(file_name,i,save_on)

load(file_name)
sn = '_xyw';
iter = num2str(i);
save_name = strcat(file_name(1:end-4),'_',iter,sn,file_name(end-3:end));

[x,y] = plotisl(p,g,i);
if i == 'end'
    lake = g.output(:,:,end)<= p.sealevel_init;
else
    lake = g.output(:,:,i)<=p.sealevel_init;
end

if ~isfield(p,'Kwave')
    p.Kwave = p.Kuniform;
end
[dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p);

ind = sub2ind(size(lake),y,x);
wave_weight = wave_weight_matrix(ind);

if save_on
    save(save_name,'x','y','wave_weight','wave_weight_matrix','fetch_matrix','save_name')
end

end