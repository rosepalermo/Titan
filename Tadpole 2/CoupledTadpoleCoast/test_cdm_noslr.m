% test_cdm_noslr
clear
% Small test lake for cdm_no_slr
addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Cellular Damage Model')
% loadtestlake
load('testlake.mat')
fetch_on = 1;
savename = 'test_cdm_noslr';
cdm_no_sealevelrise(lake,fetch_on,savename)
