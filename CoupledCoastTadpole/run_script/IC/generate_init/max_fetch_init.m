% measure max and mean fetch for all IC

% list of file names
foldername = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/generate_init/'];
filepattern = fullfile(foldername,'mrfn94*.mat');
files = dir(filepattern);

for i = 1:length(files)
    load([files(i).folder,'/',files(i).name])
    lake = init<50;
    p.Kwave = 1e-1;
    p.nrays = 180;
    p.delta = 0.05;
    p.dt = 1;
    idx_list(i) = idx;
%     [~,~,fetch_matrix,~,~,~] = get_dam_wave(lake,p);
%     max_fetch(i) = max(fetch_matrix,[],'all');
%     fetch_90(i) = quantile(fetch_matrix,0.9,'all');
end

% Ao = nanmean(max_fetch);
% Ao_cells = nanmean(max_fetch)./p.dx./p.dy;