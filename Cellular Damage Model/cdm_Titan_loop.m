% contours
load('xycontours.mat')
lakex_mr{1,1}(1,:) = x_1m_t1;
lakey_mr{1,1}(1,:) = y_1m_t1;
eps_mr(1) = 200;
dx_mr(1) = 1; dy_mr(1) = 1;

load('xycontours.mat')
lakex_mr{2,1}(1,:) = x_1m_t2;
lakey_mr{2,1}(1,:) = y_1m_t2;
eps_mr(2) = 200;
dx_mr(2) = 1; dy_mr(2) = 1;


load('xycontours.mat')
lakex_mr{3,1}(1,:) = x_1m_t3;
lakey_mr{3,1}(1,:) = y_1m_t3;
eps_mr(3) = 200;
dx_mr(3) = 1; dy_mr(3) = 1;

% th = linspace(0,2*pi,60);
% r = 2 + rand(size(th))-0.5 ;
% lakex_mr{4,1}(1,:) = r.*cos(th );
% lakey_mr{4,1}(1,:) = r.*sin(th );
% eps_mr(4) = 10;
% dx_mr(4) = 0.05; dy_mr(4) = 0.05;
% load('testlake_rand.mat')

% rng(2)
% N=100; %1000 data points
% th = linspace(0,2*pi,N)'; %theta from 0 to 2pi
% amp = 2*ones(N,1);
% for i=1:500 %higher numbers here make more higher frequency fluctuations
%     a = rand()-0.5; %random number from -0.5 to 0.5
%     b = rand()-0.5;
%     amp = amp + a/i*cos(i*th) + b/i*sin(i*th); %change circle by a fourier function over it
% end
% lakex_mr{3,1}(1,:)=amp.*cos(th); lakey_mr{3,1}(1,:) = amp.*sin(th);
% eps_mr(3) = 5;
% dx_mr(3) = 0.01; dy_mr(3) = 0.01;

for modelrun = 1:3
    lakex = lakex_mr{modelrun,1}(1,:);
    lakey = lakey_mr{modelrun,1}(1,:);
    eps = eps_mr(modelrun);
    dx = dx_mr(modelrun);
    dy = dy_mr(modelrun);
    cdm_Titan(lakex,lakey,eps,dx,dy,modelrun)
end