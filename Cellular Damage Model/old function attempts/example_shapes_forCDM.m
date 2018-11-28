% % % Generate Data
% rng(2)
% N=1000; %1000 data points
% th = linspace(0,2*pi,N)'; %theta from 0 to 2pi
% amp = 2*ones(N,1);
% for i=1:500 %higher numbers here make more higher frequency fluctuations
%     a = rand()-0.5; %random number from -0.5 to 0.5
%     b = rand()-0.5;
%     amp = amp + a/i*cos(i*th) + b/i*sin(i*th); %change circle by a fourier function over it
% end
% lakex=amp.*cos(th); lakey = amp.*sin(th);


% circle
% theta = 0 : 0.1 : 2*pi;
% radius = 2;
% lakex = radius * cos(theta);
% lakey = radius * sin(theta);

% triangle
% lakey = cat(2,0:0.01:1,1-.01:-0.01:0,zeros(1,length(1:-.01:-1)));
% lakex = cat(2,-1:0.01:1,1:-.01:-1);

% square
% lakex = .5*cat(2, -1*ones(1,length(-1:0.01:1)),-1:0.01:1,ones(1,length(-1:0.01:1)),1:-0.01:-1);
% lakey = .5*cat(2, -1:0.01:1, ones(1,length(-1:0.01:1)), 1:-0.01:-1, -1*ones(1,length(-1:0.01:1)));
