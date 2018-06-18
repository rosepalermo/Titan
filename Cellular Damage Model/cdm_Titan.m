function cdm_Titan(lakex,lakey,eps,dx,dy,modelrun)

% Titan analogue damage model for coastal erosion of a lake
% Rose Palermo 2-2018
%
%
% Input or generate a lake. Makes a grid larger than this lake.
% decides water and land based on if the grid cell is in the polygon or
% not. then incurs damage based on geometry and based on fetch.
%
%

% close all;clear all;

% tic

% fetch weighted?
fetch_on = true;

% run time
tmax = 1000;

% when creating a gif
plot_now = true;
gif_on = false;
save_on = false;
shoreline_save = cell(1,1);
filename = [num2str(modelrun),'fetch_5_2018_example.gif'];

% load('xycontours.mat')
% lakex = x_1m_t1;
% lakey = y_1m_t1;

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

% figure()
% plot(lakex,lakey)
% axis square
LakeArea = polyarea(lakex,lakey); % original lake area-- does NOT change throughout simulation

%make a grid larger than lake by eps
% eps = 200;
% dx = 0.5; dy = 0.5;
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
% [shoreline] = addidshoreline_cardonly(lake,land);
[shoreline] = addidshoreline(lake,land);

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
    if ~fetch_on % if no fetch, the order is fine like this.
        indshoreline = find(shoreline);
        dam = double(shoreline);
        %damage the shoreline
        strength(indshoreline) = strength(indshoreline)-dam(indshoreline);
        strength(strength<0) = 0;
    end
    
    if fetch_on
        if (exist('erodedind','var')) | (i == 1)
            disp('fetch')
            clearvars fetch_sl_cells indshoreline WaveArea_cell
            %order the shoreline
            shoreline = addidshoreline_cardonly(lake,land);
            [indshoreline,cells2trash] = order_cw_lastpoint(lake,shoreline); % ccw ordered ind = indshoreline
            for l = 1: length(indshoreline)
                indshoreline{l,1} = sub2ind(size(X),indshoreline{l,1}(:,1),indshoreline{l,1}(:,2));
                fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
                fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
            end
            % calculate wave weighted (sqrt(F)*cos(theta-phi))
            [WaveArea_cell] = fetch_wavefield_cell(fetch_sl_cells);
%         [WaveArea_cell] = {ones(size(fetch_sl_cells{1,1},1),1)}; % ones to test debugging with
            
            clearvars erodedind
        end
        % Damage the shoreline
        indshoreline = cell2mat(indshoreline);
        %         [shoreline] = addidshoreline_cardonly(lake,land); % edges only
        [shoreline] = addidshoreline(lake,land); % corners and edges
        corners = setdiff(find(shoreline),indshoreline);
        corners = corners(shoreline(corners)<1.5);
        dam = cell2mat(WaveArea_cell);
        [damcorn] = damagecorners(lake,corners,indshoreline,dam);
        %         strength(indshoreline) = strength(indshoreline) - ones(length(indshoreline),1).*dam;
        strength(indshoreline) = strength(indshoreline) - shoreline(indshoreline).*dam;
        strength(corners) = strength(corners) - shoreline(corners).*damcorn;
        strength(strength<0) = 0;
    end
    
    


    
    
    % find eroded points 
    erodedind = indshoreline(strength(indshoreline)<=0);
    if fetch_on
        % erode points that weren't a corner and were only 2-1 cells
        % connected because it messed up the fetch calculations..
        findsl = find(shoreline);
        if ~isempty(cells2trash)
            cells2trash = sub2ind(size(lake),cells2trash(:,1),cells2trash(:,2));
            erodedind_12 = cells2trash(find(~ismember(cells2trash,corners)));
            erodedind =[erodedind;erodedind_12];
        end
    end

         
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
%     [shoreline] = addidshoreline_cardonly(lake,land); % update to card only because only card for fetch part
    [shoreline] = addidshoreline(lake,land);
    if ~fetch_on
        [shoreline] = addidshoreline(lake,land);
    end
    
    
    
    %plot
    if plot_now
        drawnow
        %         imagesc(x,y,double(shoreline))
        imagesc(x,y,(strength))
        colormap('gray')
        shading flat
        axis square
        %         axis([-2 -0.5 0.5 2])
%         if modelrun == 1
%             axis([570 650 740 780])
%         end
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
    
    
    
    %     if save_on
    %         if ~fetch_on
    %             save([num2str(modelrun),'uniform',num2str(i),'.mat'])
    %         elseif fetch_on
    %             save([num2str(modelrun),'wave',num2str(i),'.mat'])
    %         end
    %     end
    
    shoreline_save{i,1} = find(shoreline);
    
    if save_on
        saveas(gcf,['C:\Users\Rose Palermo\Documents\Titan\Modeling\6_17_pregeneralsfigs\',num2str(modelrun),'wave',num2str(i),'.fig'])
    end
    
end

if save_on
    save(['C:\Users\Rose Palermo\Documents\Titan\Modeling\6_17_pregeneralsfigs\',num2str(modelrun),'wave','.mat'],'shoreline_save')
end
%% plot
eroded = eroded(2:end,:);
% plot initial shoreline and final shoreline
figure()
imagesc('XData',x,'YData',y,'CData',strength)
shading flat
colormap((gray))
hold on
plot(lakex,lakey,'w')
axis square
% scatter(eroded(:,1),eroded(:,2),'c')

% toc
end