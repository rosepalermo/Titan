% testing grid bias

% first 4-con vs 8-con, no uniform weighting

load('circle_init.mat')

%% NO UNIFORM WEIGHTING

% 4-connected -- DIAMOND
lake_4con = ~init;
for i = 1:200
    [sl_4con] = addidshoreline_cardonly(lake_4con,~lake_4con); % gives me shoreline with weighting
    sl_4con(sl_4con<0) = 1; % remove weighting
    lake_4con(find(sl_4con))= 1;
end
figure()
imagesc(lake_4con)

% 8-connected -- SQUARE
lake_8con = ~init;
for i = 1:200
    [sl_8con] = addidshoreline(lake_8con,~lake_8con); % gives me shoreline with weighting
    sl_8con(sl_8con<0) = 1; % remove weighting
    lake_8con(find(sl_8con))= 1;
end
figure()
imagesc(lake_8con)

%% UNIFORM WEIGHTING

% 4-connected -- square with facets
lake_4con = ~init;
strength_4con = init*10;
for i = 1:2000
    [sl_4con] = addidshoreline_cardonly(lake_4con,~lake_4con); % gives me shoreline with weighting
    strength_4con= strength_4con-sl_4con;
    lake_4con(strength_4con<0) = 1;
end
figure()
imagesc(lake_4con)

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