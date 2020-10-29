function [azi,curvi,deltamean] = CoastalAzimuth(x,y)
% The two vectors x and y contain the x and y coordinates along the coast,
% ordered CCW (if y is positive up). The first and last elements in x and y
% should not be identical (i.e., the starting point should not be 
% repeated at the end).

% Now we make x and y periodic by repeating the first 2 points at the end.
% This will also make our azimuths periodic (first and last points in the 
% azimuth vectors will be identical), which is what we want.
x = [x; x(1:2)];
y = [y; y(1:2)];

% determine the azimuth from each point to the next point
dx=diff(x);
dy=diff(y);
az=atan2(dy,dx);

% the azimuth applies to the midpoints of the segments defined by x and y
xaz=0.5*(x(1:end-1) + x(2:end));
yaz=0.5*(y(1:end-1) + y(2:end));

% calculate cumulative distance
deltad = sqrt(dx.^2 + dy.^2); % distances between shoreline points
cumdist = [0; cumsum(deltad)]; % cumulative distance of shoreline points
daz = 0.5*(cumdist(1:end-1)+cumdist(2:end)); % cumulative distance of azimuth points
daz = daz - daz(1); % make the starting distance for azimuth points zero
deltadaz = diff(daz); % distances between azimuth points
deltamean = daz(end)/(length(daz)-1); % average distance between azimuth points

% There's a problem here: if the shoreline doubles back on itself, the
% two points are the same, so the distance will be calculated as zero. Need
% to handle this special case.

% determine whether the coast makes a left or right turn at each point by
% computing the direction of the cross product of the vectors from 
% (x1,y1)-->(x2,y2) and (x2,y2)-->(x3,y3). 
% 
% If xprod > 0  left turn
%    xprod = 0  straight
%    xprod < 0  right turn
%

x1=x(1:end-1);
x2=x(2:end);
x3=x([3:end 1]);

y1=y(1:end-1);
y2=y(2:end);
y3=y([3:end 1]);

xprod = (x2 - x1).*(y3 - y1) - (y2 - y1).*(x3 - x1);
% The first element tells us about the turn from azimuth point 1 to azimuth
% point 2, etc. The last element tells us about the turn from the last
% azimuth point back to the first one.

% to get around the ambiguity of the arctan and the periodicity of the
% azimuth, we want to make sure that if the river turns left (CCW), the 
% azimuth increases, and if it turns right (CW), the azimuth decreases.

for i=2:length(az)
    if xprod(i-1)>0 % left turn from i-1 to i
        if az(i)<az(i-1)
            az(i:end)=az(i:end)+2*pi;
        end
    elseif xprod(i-1)<0 % right turn from i-1 to i
        if az(i)>az(i-1)
            az(i:end)=az(i:end)-2*pi;
        end
    elseif (xaz(i)==xaz(i-1) && yaz(i)==yaz(i-1))
        % If this az point and the previous one are identical, the
        % shoreline doubles back on itself. The shoreline is ordered CCW, so if
        % this happens, it must be a right (CW) turn, because we only have
        % doubled-back promontories, not doubled-back embayments.
        if az(i)>az(i-1)
            az(i:end)=az(i:end)-2*pi;
        end
    % otherwise xprod(i-1)==0 and there's no turn; it's straight        
    end
end

% subtract the mean, which sets the mean azimuth to zero
az = az - mean(az);

% By definition, the first and last azimuths should differ by 2*pi, BECAUSE
% THE FIRST AND LAST AZIMUTHS REPRESENT THE SEGMENT FROM (x1,y1) to (x2,y2).
% We want to make the data series periodic with a mean of zero.

% Find whether the data series has a positive or negative background slope
slopesign = sign(az(end)-az(1));
bgslope = slopesign * 2*pi / (daz(end) - daz(1));
bgx = daz - mean(daz);

% Detrend the data series by subtracting a line with a range of 2*pi from the data.
az = az - bgslope*bgx;
az = az - mean(az); % THIS ISN'T REALLY NECESSARY, BUT REMOVES A BIT OF RESIDUAL IMPRECISION
% The azimuths should now be periodic with a mean of zero, AND THE FIRST
% AND LAST POINTS SHOULD HAVE THE SAME AZIMUTH. WE WILL TRIM THE LAST POINT
% LATER.

% calculate rate of change of azimuth (d(azimuth)/d(distance)). 
% dazdd = diff([az; az(1)])./[deltadaz; deltadaz(1)]; % THIS WILL MAKE IT PERIODIC TOO
dazdd = diff(az)./deltadaz;

% interpolate azimuth and rate of change of azimuth at points evenly 
% spaced in the downstream direction
di=(0:(length(daz)-1))*deltamean;
di=di(:);
azi = interp1(daz,az,di);
% curvi = interp1(daz,dazdd,di);
curvi = interp1(daz(1:end-1),dazdd,di(1:end-1));

% FINALLY, TRIM AZIMUTHS AND CURVATURES 

% TRIM AZIMUTH BY REMOVING THE LAST POINT (WHICH IS IDENTICAL TO THE FIRST
% POINT, EVEN AFTER INTERPOLATION)
azi = azi(1:end-1);
% curvi = curvi(1:end-1); % NOT SURE THIS IS CORRECT. FIRST AND LAST POINT DON'T MATCH BEFORE TRIMMING.
