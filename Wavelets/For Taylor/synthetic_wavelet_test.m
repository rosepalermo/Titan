% test synthetic

save_on = false;
% HERE IS A SYNTHETIC SIGNAL WITH THE SAME LENGTH AS theta


th = linspace(0,2*pi,60);
r = 2 + rand(size(th))-0.5;
x = r.*cos(th);
y = r.*sin(th);
N=100; %1000 data points
th = linspace(0,2*pi,N)'; %theta from 0 to 2pi
amp = 2*ones(N,1);
for i=1:5 %higher numbers here make more higher frequency fluctuations
    a = rand()-0.5; %random number from -0.5 to 0.5
    b = rand()-0.5;
    amp = amp + a/i*cos(i*th) + b/i*sin(i*th); %change circle by a fourier function over it
end
x0=(amp.*cos(th)); y0 = (amp.*sin(th));
% add first AND SECOND points to the end for meander
x = [x0;x0(1:2)];
y = [y0;y0(1:2)];

% transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% streamwise distance
[theta, dtheta, deltad] = meander_titan(x,y);
x = x(1:end-2);
y = y(1:end-2);

i=2;
n=2;
dowave(theta,deltad,n,x,y,[],save_on,[],i);

%SECOND SYNTHETIC SIGNAL

t=deltad*(0:length(theta)-1);
t = t(:);
T = t(end)+deltad;
hwin = 0.5*(1-cos(2*pi*t/T));
synth = 4*sin(2*pi*t/(T/4)) + hwin.*( 1*sin(2*pi*t/(T/64)) );

n=2; i = 2;
dowave(synth,deltad,n,x,y,[],save_on,[],i);