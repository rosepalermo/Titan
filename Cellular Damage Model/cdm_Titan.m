
function cdm_Titan(lake,X,Y,modelrun,fetch_on,savename)

% Titan analogue damage model for coastal erosion of a lake
% Rose Palermo 2-2018
% Last update to this code was 6-2019. Used this code to make the
% waveerosion and uniformerosion functions for the coupled tadpole model.
%
% Input or generate a lake. Makes a grid larger than this lake.
% decides water and land based on if the grid cell is in the polygon or
% not. then incurs damage based on geometry and based on fetch.
%


% close all;clear all;



% tic

% fetch weighted?
% fetch_on = true;

% run time

tmax = 25;

% when creating a gif
% savefolder = 'D:\Titan\Modeling\river_and_wave_1_2019\';
% savefolder = '/home/rpalermo/titan_models';
savefolder = '/Users/rosepalermo/Documents/Research/Titan/River_and_wave_7_19';
plot_now = false;
gif_on = false;
save_on = true;
shoreline_save = cell(1,1);
% filename = [num2str(modelrun),'riverandwave.gif'];


x = X(1,:); y = Y(:,1);
% figure()
% plot(lakex,lakey)
% axis square
% LakeArea = polyarea(lakex,lakey); % original lake area-- does NOT change throughout simulation

% %make a grid larger than lake by eps
% % eps = 5;
% % dx = 0.05; dy = 0.05;
% x = (min(lakex)-eps):dx:(max(lakex)+eps);
% y = (min(lakey)-eps):dy:(max(lakey)+eps);
% [X,Y] = meshgrid(x,y);
% Xinon = reshape(X,[],1);
% Yinon = reshape(Y,[],1);
% 
% 
% %points in and on the polygon are the lake
% [in, on] = inpoly([Xinon,Yinon]',[lakex;lakey]);
% lake = in + on;
% lake = reshape(lake,length(y),length(x));
% % figure()
% % imagesc(x,y,lake)
% % shading interp

%give land some sort of strength that will be damaged and destroyed
land = ~lake;
if fetch_on
%     strength = 10000000*double(land);
    strength = 50000*double(land);
else
    strength = 10*double(land);
end

% identify the shoreline
% [shoreline,shorelinecard,shorelinecorn] = idshoreline(lake,land);
% [shoreline] = addidshoreline_cardonly(lake,land);
[shoreline] = addidshoreline(lake,land); % corners are part of the shoreline!

% indshoreline_last = find(shoreline);


% % plot the shoreline
% figure()
% imagesc(x,y,double(shoreline))
% shading flat

%% loop here for chemical weathering- style model
% ii = 1;
eroded = nan(1,2);
% h = figure;
% h.Position = [0,0,1000,500]

for i = 1:tmax
    i
% figure();imagesc(strength)
        
        
    if ~fetch_on % if no fetch, the order doesn't matter and we can calc damage.
        indshoreline = find(shoreline);
        dam = double(shoreline);
        %damage the shoreline
        strength(indshoreline) = strength(indshoreline)-dam(indshoreline);
        strength(strength<0) = 0;
        erodedind = indshoreline(strength(indshoreline)<=0);
        strength(erodedind) = 0;
        erodedX = X(erodedind);
        erodedY = Y(erodedind);
        erodedi = cat(2,erodedX,erodedY);
        eroded = cat(1,eroded,erodedi);
        
        % change land to lake at eroded pts
        lake(erodedind) = true;
        land = ~lake;
        [shoreline] = addidshoreline(lake,land);
    end
    
    if fetch_on
        
        % if fetch_on % loop over # of objects
        % find number of first order lakes
        [F_lake_all] = find_first_order_lakes(lake);
        for ff = 1:length(F_lake_all)
            
            F_lake = F_lake_all{ff};
        if (exist('erodedind','var')) | (i == 1) | (ff>1)
            disp('fetch')
            clearvars fetch_sl_cells indshoreline WaveArea_cell
            
            %order the shoreline and islands
%             shoreline = addidshoreline_cardonly(lake,land); %rewrite shoreline to be card only for fetch.. will damage corners separately later
            shoreline = addidshoreline_cardonly(F_lake,~F_lake); %rewrite shoreline to be card only for fetch.. will damage corners separately later            
            disp('ordering')
            [indshoreline_ocw,cells2trash] = order_cw_lastpoint(F_lake,shoreline); % ccw ordered ind = indshoreline
            disp('ordered')
            for l = 1: length(indshoreline_ocw)
                indshoreline{l,1} = sub2ind(size(X),indshoreline_ocw{l,1}(:,1),indshoreline_ocw{l,1}(:,2));
                fetch_sl_cells{l,1}(:,1) = X(indshoreline{l,1});
                fetch_sl_cells{l,1}(:,2) = Y(indshoreline{l,1});
            end
            disp('calculating wave')
            % calculate wave weighted (sqrt(F)*cos(theta-phi))
%             [WaveArea_cell,~] = fetch_wavefield_cell(fetch_sl_cells);
            [WaveArea_cell] = fetch_vis_approx(fetch_sl_cells);
 
%         [WaveArea_cell] = {ones(size(fetch_sl_cells{1,1},1),1)}; % ones to test debugging with
            disp('wave calculated')
            clearvars erodedind
        end
        
        % Damage the shoreline
        indshoreline = cell2mat(indshoreline);
        %         [shoreline] = addidshoreline_cardonly(lake,land); % edges only
        [shoreline] = addidshoreline(lake,land); % corners and edges
        dam = cell2mat(WaveArea_cell);
%         dam = WaveArea_cell';
        strength(indshoreline) = strength(indshoreline) - shoreline(indshoreline).*dam;
        
        % find corners and damage if they exist alone (not in the cells to
        % trash)
        corners = setdiff(find(shoreline),indshoreline);
        if ~isempty(cells2trash)
            [c2t]=sub2ind(size(X),cells2trash(:,1),cells2trash(:,2)); % find the cells that arent the trash cells
            if exist('c2t') & ~isempty(corners)
                corners = corners(~ismember(corners,c2t));
            end
        end
        
        % find the mean of the damage for points next to the corners. make that
        % the damage for that corner
        if ~isempty(corners)
            corners = corners(shoreline(corners)<1.5); % if less than 1.5, only a corner. not also a side
            [damcorn] = damagecorners(lake,corners,indshoreline,dam);
            %         strength(indshoreline) = strength(indshoreline) - ones(length(indshoreline),1).*dam;
            strength(corners) = strength(corners) - shoreline(corners).*damcorn;
        end
        
        strength(strength<0) = 0; % if strength is negative, make it 0 for convenience
        
        
        % find eroded points
        erodedind = indshoreline(strength(indshoreline)<=0);
        % erode points that weren't a corner and were only 2-1 cells
        % connected because it messed up the fetch calculations..
        if ~isempty(cells2trash)
            cells2trash = sub2ind(size(lake),cells2trash(:,1),cells2trash(:,2));
            erodedind_12 = cells2trash(find(~ismember(cells2trash,corners)));
            erodedind =[erodedind;erodedind_12];
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
        ordered_sl_save{i,ff} = fetch_sl_cells;
        corners_save{i,ff} = corners;
        damcorners_save{i,ff} = damcorn;
        end
    end
    
    
    
    
    %plot
    if plot_now
        drawnow
        %         imagesc(x,y,double(shoreline))
        p1 = subplot(1,2,1)
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
    dam_save{i,1} = dam;
    lake_save{i,1} = lake;
%     if fetch_on
%         ordered_sl_save{i,1} = fetch_sl_cells;
%         corners_save{i,1} = corners;
%         damcorners_save{i,1} = damcorn;
    
%     slplot = cell2mat(fetch_sl_cells);
%     p2 = subplot(1,2,2)
%     cla(p2)
%     scatter3(slplot(:,1),slplot(:,2),dam,[],dam)
%     hold on
%     scatter3(X(corners),Y(corners),damcorn,[],damcorn)
%     axis equal
%     axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y))])
%     colorbar
%     colormap(jet)
%     view(2)
%     end
    
    
    if save_on && plot_now
        %         saveas(gcf,['C:\Users\Rose Palermo\Documents\Titan\Modeling\6_17_pregeneralsfigs\',num2str(modelrun),'wave',num2str(i),'.fig'])
        %         saveas(gcf,['C:\Users\Rose Palermo\Documents\Titan\Modeling\6_17_pregeneralsfigs\','wave_rednoise',num2str(i),'.fig'])
        
%         saveas(gcf,['D:\Titan\Modeling\river_and_wave_1_2019\',savename,num2str(i),'.fig'])
        saveas(gcf,[savefolder, savename,num2str(i),'.fig'])

    end
    
    if fetch_on
        if save_on
            %     save(['C:\Users\Rose Palermo\Documents\Titan\Modeling\6_17_pregeneralsfigs\',num2str(modelrun),'wave','.mat'],'shoreline_save')
%             save(['D:\Titan\Modeling\river_and_wave_1_2019\',savename,'.mat'],'shoreline_save','ordered_sl_save','dam_save','corners_save','damcorners_save','X','Y')
            save([savefolder, savename,'.mat'],'shoreline_save','ordered_sl_save','dam_save','corners_save','damcorners_save','X','Y','lake_save')
            
        end
    else
%         save(['D:\Titan\Modeling\river_and_wave_1_2019\',savename,'.mat'],'shoreline_save','dam_save','X','Y','lake_save')
        save([savefolder, savename,'.mat'],'shoreline_save','dam_save','X','Y','lake_save')

    end
end

%% plot
eroded = eroded(2:end,:);
% plot initial shoreline and final shoreline
% figure()
% imagesc('XData',x,'YData',y,'CData',strength)
% shading flat
% colormap((gray))
% hold on
% plot(lakex,lakey,'w')
% axis square
% scatter(eroded(:,1),eroded(:,2),'c')

% toc
end