% titan_wavelet_test.m

clear 

load meander_test_wave x y

% The first and last points in x and y should not be the same. Periodic 
% endpoints will be taken care of inside the function.

[azi,curvi,deltad] = CoastalAzimuth(x,y);