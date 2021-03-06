
function cdm_no_sealevelrise(lake,fetch_on,savename)

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

tmax = 20;
p.doStreamPower = 0;

% when creating a gif
% savefolder = 'D:\Titan\Modeling\river_and_wave_1_2019\';
% savefolder = '/home/rpalermo/titan_models/';
savefolder = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/River_and_wave_9_19/';
plot_now = false;
gif_on = false;
save_on = false;
lake_save = cell(tmax,1);
lake_save{1} = lake;
strength = cell(tmax,1);
% filename = [num2str(modelrun),'riverandwave.gif'];


%give land some sort of strength that will be damaged and destroyed
land = ~lake;
if fetch_on
        strength{1} = 500000000*double(land);
%     strength{1} = 50000*double(land);
else
    strength{1} = 10*double(land);
end


for i = 2:tmax
    i
    % figure();imagesc(strength)
    
    
    [lake_save{i},strength{i}] = coastal_erosion(lake_save{i-1},fetch_on,strength{i-1},p);
    
    %plot
    if plot_now
        drawnow
        %         imagesc(x,y,double(shoreline))
        p1 = subplot(1,2,1)
        imagesc((strength{i}))
        colormap('gray')
        shading flat
        axis square
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
    
    
    
    
%     shoreline_save{i,1} = find(shoreline);
%     lake_save{i,1} = lake;
    
    % save a plot if there is one
    if save_on && plot_now
        %         saveas(gcf,['C:\Users\Rose Palermo\Documents\Titan\Modeling\6_17_pregeneralsfigs\',num2str(modelrun),'wave',num2str(i),'.fig'])
        %         saveas(gcf,['C:\Users\Rose Palermo\Documents\Titan\Modeling\6_17_pregeneralsfigs\','wave_rednoise',num2str(i),'.fig'])
        
        %         saveas(gcf,['D:\Titan\Modeling\river_and_wave_1_2019\',savename,num2str(i),'.fig'])
%         saveas(gcf,[savefolder, savename,num2str(i),'.fig'])
        
    end
    
    %save the data
    if save_on
        save([savefolder, savename,'_',num2str(i),'.mat'],'lake_save','strength')
    end
end

end