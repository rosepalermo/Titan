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

rng(2)
N=1000; %1000 data points
th = linspace(0,2*pi,N)'; %theta from 0 to 2pi
amp = 2*ones(N,1);
for i=1:500 %higher numbers here make more higher frequency fluctuations
    a = rand()-0.5; %random number from -0.5 to 0.5
    b = rand()-0.5;
    amp = amp + a/i*cos(i*th) + b/i*sin(i*th); %change circle by a fourier function over it
end
lakex=amp.*cos(th); lakey = amp.*sin(th);
lakex_mr{7,1}(1,:) = lakex;
lakey_mr{7,1}(1,:) = lakey;
eps_mr(7) = 1;
dx_mr(7) = .01; dy_mr(7) = .01;
nametemp{7} = 'test_random';


for fo = 1:2
    fo_ = [1 0];
    fetch_on = fo_(fo);
    savetemp = {'wave' 'uniform'};
    for mr = 2:5
        mr_ = [1 2 3 4 5 6 7];
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