% contours
load('xy_beta1_6_v1.mat')
lakex_mr{1,1}(1,:) = x_0_5m_t1v1;
lakey_mr{1,1}(1,:) = y_0_5m_t1v1;
eps_mr(1) = 100;
dx_mr(1) = 1; dy_mr(1) = 1;
nametemp{1} = 't1v1';

lakex_mr{2,1}(1,:) = x_0_5m_t2v1;
lakey_mr{2,1}(1,:) = y_0_5m_t2v1;
eps_mr(2) = 100;
dx_mr(2) = 1; dy_mr(2) = 1;
nametemp{2} = 't2v1';

load('xy_beta1_6_v2.mat')
lakex_mr{3,1}(1,:) = x_0_5m_t1v2;
lakey_mr{3,1}(1,:) = y_0_5m_t1v2;
eps_mr(3) = 100;
dx_mr(3) = 1; dy_mr(3) = 1;
nametemp{3} = 't1v2';

lakex_mr{4,1}(1,:) = x_0_5m_t2v2;
lakey_mr{4,1}(1,:) = y_0_5m_t2v2;
eps_mr(4) = 100;
dx_mr(4) = 1; dy_mr(4) = 1;
nametemp{4} = 't2v2';

lakex_mr{5,1}(1,:) = x_0_5m_t3v2;
lakey_mr{5,1}(1,:) = y_0_5m_t3v2;
eps_mr(5) = 100;
dx_mr(5) = 1; dy_mr(5) = 1;
nametemp{5} = 't3v2';

load('xy_beta1_6_v1.mat')
lakex_mr{6,1}(1,:) = x_0_5m_t1v1;
lakey_mr{6,1}(1,:) = y_0_5m_t1v1;
eps_mr(6) = 200;
dx_mr(6) = 5; dy_mr(6) = 5;
nametemp{6} = 'test_t1v1';

for fo = 1:2
    fo_ = [1 0];
    fetch_on = fo_(fo);
    savetemp = {'wave' 'uniform'};
    for mr = 6
        mr_ = [1 2 3 4 5 6];
        modelrun = mr_(mr);
        lakex = lakex_mr{modelrun,1}(1,:);
        lakey = lakey_mr{modelrun,1}(1,:);
%         figure()
%         plot(lakex,lakey)
        eps = eps_mr(modelrun);
        dx = dx_mr(modelrun);
        dy = dy_mr(modelrun);
        savename = [savetemp{1,fo},nametemp{1,modelrun}]
        cdm_Titan(lakex,lakey,eps,dx,dy,modelrun,fetch_on,savename)
        close all
    end
end