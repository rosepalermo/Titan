

rfactor = 0.25; % 0.25; % depth of the depression as a function of relief of the noise surface
init_temp = cell(0);
ii = 1;
for idx = 1:10000
load('p_init_v1.mat')
p.dtmax = 1e6;
p.Kf = 5e-6;
p.rand_gen = 1;
p.size_final = 2.4;
p.saveint = 40;
% p.sealevel_init = 20;
[init_temp,p] = get_IC(p,rfactor,idx);

solution = Tadpole(init_temp,p);

init{idx} = solution(:,:,end);

lake = init{idx}<40;

% does the SL hit a boundary?
[p] = check_boundary(lake,p);

% if it doesn't hit the boundary, then we're good
if ~isfield(p,'boundary')
   idx_save(ii) = idx;
   ii = ii+1;
end
clearvars p
end

save('ic_generated_v1','init','idx_save')
