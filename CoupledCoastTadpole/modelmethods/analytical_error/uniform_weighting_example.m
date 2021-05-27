% example uniform weighting on a shoreline

load('lake_uniformweighting_example.mat')
p.con8 = 1;
[shoreline] = addidshoreline(lake,~lake,p);
imagesc(shoreline)