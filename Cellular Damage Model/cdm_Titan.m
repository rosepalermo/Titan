% Titan analogue damage model for coastal erosion of a lake
% Rose Palermo 2-2018
%
%
% Input or generate a lake. Makes a grid larger than this lake.
% decides water and land based on if the grid cell is in the polygon or
% not. then incurs damage based on geometry and based on fetch.
%
%

close all;clear all;

tic

% fetch weighted?
fetch_on = false;

% run time
tmax = 100;

% when creating a gif
plot_now = true;
gif_on = true;
filename = 'fetch_5_2018_example.gif';

load('xycontours.mat')
lakex = x_1m_t1;
lakey = y_1m_t1;

% % Test lake
% th = linspace(0,2*pi,60);
% r = 2 + rand(size(th))-0.5 ;
% lakex = r.*cos(th );
% lakey = r.*sin(th );


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

figure()
plot(lakex,lakey)
axis square
LakeArea = polyarea(lakex,lakey);

%make a grid larger than lake by eps
eps = 200;
dx = 0.5; dy = 0.5;
x = (min(lakex)-eps):dx:(max(lakex)+eps);
y = (min(lakey)-eps):dy:(max(lakey)+eps);
[X,Y] = meshgrid(x,y);
Xinon = reshape(X,[],1);
Yinon = reshape(Y,[],1);


%points in and on the polygon are the lake
[in, on] = inpolygon(Xinon,Yinon,lakex,lakey);
lake = in + on;
lake = reshape(lake,length(y),length(x));
% figure()
% imagesc(x,y,lake)
% shading interp

%give land some sort of strength that will be damaged and destroyed
land = ~lake;
strength = 10*double(land);

% identify the shoreline
% [shoreline,shorelinecard,shorelinecorn] = idshoreline(lake,land);
[shoreline] = addidshoreline_cardonly(lake,land);
% indshoreline_last = find(shoreline);


% % plot the shoreline
% figure()
% imagesc(x,y,double(shoreline))
% shading flat

%% loop here for chemical weathering- style model
% ii = 1;
eroded = nan(1,2);
h = figure;
for i = 1:tmax
    i
    %     if ~fetch_on % if no fetch, the order is fine like this.
    % create array of uniform damage and subtract
    indshoreline = find(shoreline);
    [slr,slc] = find(shoreline);
    %     end
    
    if fetch_on
        if (exist('erodedind','var')) | (i == 1)
            disp('fetch')
            %order the shoreline
            [indshoreline] = order_cw_lastpoint(lake,shoreline); % ccw ordered ind = indshoreline

            % calculate fetch
            [FetchArea] = fetch_mw(X(indshoreline(:,1)),Y(indshoreline(:,2)));
            normfetch = FetchArea./LakeArea; % divide fetch area by original lake area
            
            clearvars erodedind
        end
    end
    
    
    
    % to damage shoreline by shoreline matrix
    if fetch_on
        dam = ones(1,length(indshoreline)) .* normfetch;
        strength(indshoreline) = strength(indshoreline) - dam';
    end
    if ~fetch_on
        dam = double(shoreline);
        strength(indshoreline) = strength(indshoreline)-dam(indshoreline);
    end
    
    
    % find eroded points
    erodedind = indshoreline(strength(indshoreline)<=0);
    strength(erodedind) = 0;
    erodedX = X(erodedind);
    erodedY = Y(erodedind);
    erodedi = cat(2,erodedX,erodedY);
    eroded = cat(1,eroded,erodedi);
    
    % change land to lake at eroded pts
    lake(erodedind) = true;
    land = ~lake;
    % update shoreline
    %     [shoreline,shorelinecard,shorelinecorn] = idshoreline(lake,land);
    [shoreline] = addidshoreline_cardonly(lake,land);
    
    
    
    
    %plot
    if plot_now
        drawnow
        imagesc(x,y,double(shoreline))
        shading flat
        axis square
        %         axis([-2 -0.5 0.5 2])
        str = sprintf('Time step = %d',i);
        title(str)
        hold on
        %              scatter(eroded(:,1),eroded(:,2),'c')
    end
    if gif_on
        if i == tmax
            hold on
            plot(lakex,lakey,'w')
        end
        
        if (length(eroded(:,1))>1 && rem(i,1)==0) || i == 1
            %Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            %Write to the GIF File
            if i == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
        end
    end
    
    save(['test',num2str(i),'.mat'])
    
end


%% plot
eroded = eroded(2:end,:);
% plot initial shoreline and final shoreline
figure()
imagesc('XData',x,'YData',y,'CData',strength)
shading flat
colormap((parula))
hold on
plot(lakex,lakey,'w')
axis square
% scatter(eroded(:,1),eroded(:,2),'c')

toc