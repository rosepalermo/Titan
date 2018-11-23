% contours
load('xy_beta1_6_v1.mat')
lakex_mr{1,1}(1,:) = x_0_5m_t1v1;
lakey_mr{1,1}(1,:) = y_0_5m_t1v1;
eps_mr(1) = 200;
dx_mr(1) = 1; dy_mr(1) = 1;
nametemp{1} = 't1v1';

lakex_mr{2,1}(1,:) = x_0_5m_t3v1;
lakey_mr{2,1}(1,:) = y_0_5m_t3v1;
eps_mr(2) = 200;
dx_mr(2) = 1; dy_mr(2) = 1;
nametemp{2} = 't3v1';

load('xy_beta1_6_v2.mat')
lakex_mr{3,1}(1,:) = x_0_5m_t1v2;
lakey_mr{3,1}(1,:) = y_0_5m_t1v2;
eps_mr(3) = 200;
dx_mr(3) = 1; dy_mr(3) = 1;
nametemp{3} = 't1v2';

lakex_mr{4,1}(1,:) = x_0_5m_t3v2;
lakey_mr{4,1}(1,:) = y_0_5m_t3v2;
eps_mr(4) = 200;
dx_mr(4) = 1; dy_mr(4) = 1;
nametemp{4} = 't3v2';

for fo = 1:2
    fo_ = [1 0];
    fetch_on = fo_(fo);
    savetemp = {'wave' 'uniform'};
    for modelrun = 1:4
        lakex = lakex_mr{modelrun,1}(1,:);
        lakey = lakey_mr{modelrun,1}(1,:);
        eps = eps_mr(modelrun);
        dx = dx_mr(modelrun);
        dy = dy_mr(modelrun);
        savename = [savetemp{1,fo},nametemp{1,modelrun}];
        cdm_Titan(lakex,lakey,eps,dx,dy,modelrun,fetch_on,savename)
    end
end