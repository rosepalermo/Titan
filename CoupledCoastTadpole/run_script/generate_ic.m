cd ..

folder = fileparts(which('getpath_CCT.m'));
addpath(genpath(folder));

rfactor = 0.25; % 0.25; % depth of the depression as a function of relief of the noise surface
init_temp = cell(0);
% ii = 1;

for idx = 1000:2000
    load('p_init_v1.mat')
    p.relief_final = 0.92;
    p.size_final = 1;
    p.dtmax = 1e6;
    p.Kf = 5e-6;
    p.rand_gen = 1;
    % p.size_final = 2.4;
    p.saveint = 100;
    % p.sealevel_init = 20;
    [init_temp,p] = get_IC(p,rfactor,idx);
    p.relief_init = mean(init_temp,'all')-min(init_temp,[],'all');
    solution = Tadpole(init_temp,p);
    
    init = solution(:,:,end);
    
    lake = init<40;
    
    % does the SL hit a boundary?
    [p] = check_boundary(lake,p);
    
    % if it doesn't hit the boundary, then we're good
    if ~isfield(p,'boundary')
        %     idx_save(ii) = idx;
        %     init_save{ii} = init;
        %     ii = ii+1;
        folder = '/home/rpalermo/TitanModelOutput/generate_init/';
        savename = ['mrfn92_idx_',num2str(idx)];
        f_savename = [folder,savename];
        save(f_savename,'init','idx','p')
    end
    
    clearvars p
end
