% calculate x,y and fetch

function [x,y,wave_weight,save_name] = calc_xyw(file_name,i,save_on,lake,p)
temp = exist('lake');
if ~temp
load(file_name)
end
p.sl_analysis = 1;
sn = '_xyw';
iter = num2str(i);
save_name = strcat(file_name(1:end-4),'_',iter,sn,file_name(end-3:end));
if ~temp
        [x,y] = plotisl(p,g,i,lake);
    if i == 'end'
        lake = g.output(:,:,end)<= p.sealevel_init;
    else
        lake = g.output(:,:,i)<=p.sealevel_init;
    end
else
    [x,y] = plotisl(p,[],i,lake);
end
if ~isfield(p,'Kwave')
    p.Kwave = p.Kuniform;
end
[~,wave_weight_matrix,fetch_matrix,~,~,~] = get_dam_wave(lake,p);

ind = sub2ind(size(lake),y,x);
wave_weight = wave_weight_matrix(ind);

if save_on
    save(save_name,'x','y','wave_weight','wave_weight_matrix','fetch_matrix','save_name')
end

end
