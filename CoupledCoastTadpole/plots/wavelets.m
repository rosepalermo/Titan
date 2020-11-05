function [period, eq14] = wavelets(x,y,fetch,pmin,pmax,save_name)

% calculate wavelet power spectrum
% input is x and y coordinates

% transform into azimuth and d(azimuth)/d(distance), evenly spaced in
% streamwise distance
% [theta, dtheta, deltad] = meander_titan(x,y);
[theta,curvi,deltad] = CoastalAzimuth(x,y);


n=2;

[period,eq14] = dowave_greece(theta,deltad,n,x,y,save_name,1,fetch,pmin,pmax);

end