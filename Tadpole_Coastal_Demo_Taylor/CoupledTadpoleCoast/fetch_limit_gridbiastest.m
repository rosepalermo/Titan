% wave grid bias testing
% at what fetch length does everything just get eroded and lose weighting
% to grid bias?
load('circle_init.mat')
% for the circle case, dt=1 & p.Kcoast = 1
dt=1;
Kcoast=1;

% UNIFORM ONLY
% 8-connected -- polygon approximation
lake_8con = ~init;
strength_8con = init*10;

for i = 1:2000
    [sl_8con] = addidshoreline(lake_8con,~lake_8con); % gives me shoreline with weighting
    strength_8con= strength_8con-sl_8con;
    lake_8con(strength_8con<0) = 1;
end
figure()
imagesc(lake_8con)


% Wave with no fetch limitation
lake_nolim = ~init;
strength_8con = init*10;

for i = 1:2000
    [sl_nolim] = addidshoreline(lake_nolim,~lake_nolim); % gives me shoreline with weighting
    [F_lake_all_nolim,~,~,~] = find_first_order_lakes(lake_nolim);
    for ff = 1:length(F_lake_all_nolim)
        F_lake = F_lake_all_nolim{ff};
        if length(find(F_lake))<2
            continue
        end
        [fetch,fetch_corn,corners,indshoreline] = find_fetch(F_lake,X,Y,lake);
        strength(indshoreline) = strength(indshoreline) - p.dt*p.Kcoast*shoreline(indshoreline).*dam; % Taylor's modified line that depends on a rate constant

    end
    strength_nolim= strength_nolim-sl_nolim;
    lake_nolim(strength_nolim<0) = 1;
end
figure()
imagesc(lake_nolim)

test = dt*Kcoast*shoreline(indshoreline).*dam;
figure()
plot(test)
testing = ones(size(shoreline));
testing(indshoreline) = testing(indshoreline)-test;
imagesc(testing)

% Wave with limit of 1 unit
fetch_limit = 1;


% Wave with limit of 10 units
fetch_limit = 10;


% Wave with limit of 100 units
fetch_limit = 100;

% Wave with limit of 1000 units
fetch_limit = 1000;
