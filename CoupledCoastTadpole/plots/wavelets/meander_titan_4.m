function [azi curvi deltad] = meander_titan4(x,y)
% meander_titan.m
% changed to subtract a line with slope 2*pi

% the workspace should include two vectors, x and y, which contain the x
% and y coordinates along the river in the downstream direction.

% determine the azimuth from each point to the next point
dx=diff(x);
dy=diff(y);
az=atan2(dy,dx);

% the azimuth applies to the midpoints of the segments defined by x and y
xaz=0.5*(x(1:end-1) + x(2:end));
yaz=0.5*(y(1:end-1) + y(2:end));

% calculate cumulative streamwise distance
dxaz=diff(xaz);
dyaz=diff(yaz);
dd=sqrt(dxaz.^2+dyaz.^2); % distance increment
d=[0; cumsum(dd)];
% d=d-min(d); % cumulative streamwise distance
deltad=mean(dd); % the mean distance spacing


% determine whether the river makes a left or right turn at each point by
% computing the direction of the cross product of the vectors from 
% (x1,y1)-->(x2,y2) and (x2,y2)-->(x3,y3). Note that this excludes the 
% first and last points.
% 
% If xprod > 0  left turn
%    xprod = 0  straight
%    xprod < 0  right turn
%
% Note that we assume the first link is straight %ROSE: CHECK THAT THIS IS
% NOT THE SOURCE OF SHORT-WAVELENGTH POWER AT ENDPOINTS

x1=x(1:end-2);
x2=x(2:end-1);
x3=x(3:end);

y1=y(1:end-2);
y2=y(2:end-1);
y3=y(3:end);

xprod = [0; (x2 - x1).*(y3 - y1) - (y2 - y1).*(x3 - x1)];


% to get around the ambiguity of the arctan and the periodicity of the
% azimuth, we want to make sure that if the river turns left (CCW), the 
% azimuth increases, and if it turns right (CW), the azimuth decreases.

for i=2:length(az)
    if xprod(i)>0 % left turn
        if az(i)<az(i-1)
            az(i:end)=az(i:end)+2*pi;
        end
    elseif xprod(i)<0 % right turn
        if az(i)>az(i-1)
            az(i:end)=az(i:end)-2*pi;
        end
    % otherwise there's no turn, it's straight        
    end
end

% subtract the mean, which sets the mean flow direction to zero
az = az - mean(az);

% NEW CODE FOR TITAN LAKES/SEAS
% By definition, the first and last azimuths should differ by 2*pi.
% We want to make the data series periodic.

% Find whether the data series has a positive or negative background slope
% slopesign = sign(az(end)-az(1));
bgslope = (az(end)-az(1))/(d(end) - d(1));
% bgslope = slopesign * 2*pi / (d(end) - d(1));
bgx = d - mean(d);

% Detrend the data series by subtracting a line with a range of 2*pi from the data.
az = az - bgslope*bgx;
az = az - mean(az);
% The azimuths should now be periodic with a mean of zero.

% END NEW CODE FOR TITAN LAKES/SEAS

% calculate rate of change of azimuth (d(azimuth)/d(distance)). This will exclude the first and last
% points of the original river path
dazdd = diff(az)./dd;

% interpolate azimuth and rate of change of azimuth at points evenly 
% spaced in the downstream direction
di=(0:(length(d)-1))*deltad;
di=di(:);
azi = interp1(d,az,di);
curvi = interp1(d(1:end-1),dazdd,di(1:end-1));
curvi = curvi(~isnan(curvi));
