% input k and folder

function [Kc_,folder] = inputs_k_folder(cluster);

% cluster -- 1 running on cluster, 2 running locally

if cluster
    folder = '/home/rpalermo/TitanModelOutput/093020/';
    
elseif ~cluster
    folder = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/093020/';
    
end

% K
% Kf_ = [5e-10 5e-8 5e-6]; % rivers
Kc_ = [1e-2 1.5e-2 1.5e-1]; % uniform/wave
% Kc_ = [1e-4; 1e-3; 1e-2];% 5e-13 2e-13];
% Kc_ = 1.5e-8*1000;%for circle, wave p.Kcoast = 1.5e-8
% Kc_ = p.Kf;
end