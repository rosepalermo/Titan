% addpath('/home/rpalermo/Titan2/Tadpole 2/example_Rose_800x800')
% addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Tadpole 2/example_Rose_800x800')
% addpath('/Users/rosepalermo/Dropbox (MIT)/Titan/taylor_fetch/Rosemakingapproxfn')
% addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Cellular Damage Model')
% contours
load('xy_beta1_6_v1.mat')
lakex_mr{1,1}(1,:) = x_0_5m_t1v1;
lakey_mr{1,1}(1,:) = y_0_5m_t1v1;
eps_mr(1) = 100;
dx_mr(1) = 0.25; dy_mr(1) = 0.25;
nametemp{1} = 't1v1';

lakex_mr{2,1}(1,:) = x_0_5m_t2v1;
lakey_mr{2,1}(1,:) = y_0_5m_t2v1;
eps_mr(2) = 100;
dx_mr(2) = 0.25; dy_mr(2) = 0.25;
nametemp{2} = 't2v1';

load('xy_beta1_6_v2.mat')
lakex_mr{3,1}(1,:) = x_0_5m_t1v2;
lakey_mr{3,1}(1,:) = y_0_5m_t1v2;
eps_mr(3) = 100;
dx_mr(3) = 0.25; dy_mr(3) = 0.25;
nametemp{3} = 't1v2';

lakex_mr{4,1}(1,:) = x_0_5m_t2v2;
lakey_mr{4,1}(1,:) = y_0_5m_t2v2;
eps_mr(4) = 100;
dx_mr(4) = 0.25; dy_mr(4) = 0.25;
nametemp{4} = 't2v2';

lakex_mr{5,1}(1,:) = x_0_5m_t3v2;
lakey_mr{5,1}(1,:) = y_0_5m_t3v2;
eps_mr(5) = 100;
dx_mr(5) = 0.25; dy_mr(5) = 0.25;
nametemp{5} = 't3v2';

load('xy_beta1_6_v1.mat')
lakex_mr{6,1}(1,:) = x_0_5m_t1v1;
lakey_mr{6,1}(1,:) = y_0_5m_t1v1;
eps_mr(6) = 200;
dx_mr(6) = 10; dy_mr(6) = 10;
nametemp{6} = 'test_t1v1';


for fo = 2
    fo_ = [1 0];
    fetch_on = fo_(fo);
    savetemp = {'wave' 'uniform'};
    for mr = 2
        savename = [savetemp{1,fo},nametemp{1,mr}]
        [lake,~,~] = gridlake(lakex_mr{mr,1}(1,:),lakey_mr{mr,1}(1,:),dx_mr(mr),dy_mr(mr),eps_mr(mr));
        cdm_no_sealevelrise(lake,fetch_on,savename)
        close all
    end
end