function [shoreline] = addidshoreline_cardonly(lake,land)

%identify the shoreline by checking neighbors (boundaries btw land and lake) in logical
%
% input: logical matrix that describes the location of a lake (true = lake; false = land) 
% input: land. land = ~lake; -- if only one lake.
% output: the location of the land around the shoreline.
%


% separate lake matrix into 2 to check against each other for above and
% below neighbors
lakea = lake(:,1:end-1);
lakeb = lake(:,2:end);
edgeb = xor(lakea,lakeb); % true when boundary with neighbor below is land/lake

% separate lake matrix into 2 to check against each other for left and
% right neighbors
lakel = lake(1:end-1,:);
laker = lake(2:end,:);
edger = xor(lakel,laker); % true when boundary with neighbor to right is land/lake

% bottom left corner
lakea = lake(1:end-1,1:end-1);
lakeb = lake(2:end,2:end);
edgebl = xor(lakea,lakeb); % true when boundary with neighbor below is land/lake

% bottom right corner
lakea = lake(1:end-1,2:end);
lakeb = lake(2:end,1:end-1);
edgebr = xor(lakea,lakeb); % true when boundary with neighbor below is land/lake


%% add 1234
% 
% shorelinecard = zeros(size(lake));
% shorelinecorn = zeros(size(lake));
% 
% shoreline = zeros(size(lake));
% shoreline(:,1:end-1) = land(:,1:end-1) & edgeb; %is this cell land and boundary below?
% shoreline(:,2:end) = shoreline(:,2:end) + (land(:,2:end) & edgeb);% is this cell land and boundary above?
% shoreline(1:end-1,:) = shoreline(1:end-1,:) + (land(1:end-1,:) & edger); % is this cell land and boundary right?
% shoreline(2:end,:) = shoreline(2:end,:) + (land(2:end,:) & edger); % is this cell land and boundary left?
% shoreline(1:end-1,1:end-1) = shoreline(1:end-1,1:end-1) + (land(1:end-1,1:end-1) & edgebl); % is land and bottom left boundary of lake?
% shoreline(2:end,2:end) = shoreline(2:end,2:end) + (land(2:end,2:end) & edgebl); %is land and top right boundary of lake?
% shoreline(1:end-1,2:end) = shoreline(1:end-1,2:end) + (land(1:end-1,2:end) & edgebr); %is land and bottom right boundary of lake?
% shoreline(2:end,1:end-1) = shoreline(2:end,1:end-1) + (land(2:end,1:end-1) & edgebr); % is land and top left boundary of lake?



%% #edges + #corners/sqrt2

shorelinecard = zeros(size(lake));
shorelinecard(:,1:end-1) = land(:,1:end-1) & edgeb; %is this cell land and boundary below?
shorelinecard(:,2:end) = shorelinecard(:,2:end) + (land(:,2:end) & edgeb);% is this cell land and boundary above?
shorelinecard(1:end-1,:) = shorelinecard(1:end-1,:) + (land(1:end-1,:) & edger); % is this cell land and boundary right?
shorelinecard(2:end,:) = shorelinecard(2:end,:) + (land(2:end,:) & edger); % is this cell land and boundary left?
% shorelinecorn = zeros(size(lake));
% shorelinecorn(1:end-1,1:end-1) = land(1:end-1,1:end-1) & edgebl; % is land and bottom left boundary of lake?
% shorelinecorn(2:end,2:end) = shorelinecorn(2:end,2:end) + (land(2:end,2:end) & edgebl); %is land and top right boundary of lake?
% shorelinecorn(1:end-1,2:end) = shorelinecorn(1:end-1,2:end) + (land(1:end-1,2:end) & edgebr); %is land and bottom right boundary of lake?
% shorelinecorn(2:end,1:end-1) = shorelinecorn(2:end,1:end-1) + (land(2:end,1:end-1) & edgebr); % is land and top left boundary of lake?

shoreline = shorelinecard;% + shorelinecorn./sqrt(2);


% figure()
% imagesc(X,Y,double(shoreline))
% shading flat