% close all
% clear
addpath('/Users/rosepalermo/Documents/GitHub/Titan2/Wavelets/Current Working folder')
save_on = false;
%     load('river_coastal_sl4wavelets.mat')
load('original_river_incised.mat')
for i = 3
    
%     clearvars -except period global_Save i save_on eq14save sl_cell_wave sl_cell_uniform
    
    % 2 = REDNOISE
    savename{1} = 'test';
    % 3 = WAVE t1v1_30
    savename{2} = 'test';
    
    fetch = [];
    if save_on
        savename = savename{i};
    end
    
    if i == 1
        sl = sl_cell_wave;
    elseif i == 2
        sl = sl_cell_uniform;
    elseif i == 3
        sl = sl_river;
    end
    
    for ii = 1:1
        %     xx = A(5).cord{1,1}(:,2); yy = A(5).cord{1,1}(:,1);
        %     A(5).cord{1,1}(:,1) = yy;
        %     A(5).cord{1,1}(:,2) = xx;
        x0 = sl{ii}(:,1);
        y0 = sl{ii}(:,2);
        
        % fetch = WaveArea_save{1, i}{1, 1};
        M = [x0 y0];
        % indices to unique values in column 3
        [~, ind] = unique(M, 'rows');
        % duplicate indices
        duplicate_ind = setdiff(1:size(M, 1), ind);
        % duplicate values
        duplicate_val = [M(duplicate_ind, 1) M(duplicate_ind, 2)];
        M(duplicate_ind,:)=[];
        % fetch(duplicate_ind) = [];
        
        
        x0=M(:,1);
        y0=M(:,2);
        
        
        
        % add first AND SECOND points to the end for meander
        x = [x0;x0(1:2)];
        y = [y0;y0(1:2)];
        
        
        % transform into azimuth and d(azimuth)/d(distance), evenly spaced in
        % streamwise distance
        [theta, dtheta, deltad] = meander_titan(x,y);
        x = x(1:end-2);
        y = y(1:end-2);
        deltad = deltad*10;
        
        % do the wavelet transforms, comparing the spectra to those of an
        % autoregressive process of order n (a.k.a. an AR(n) process)
        n=2;
        % dowave_duplicate(theta,deltad,n,x,y,savename);
        
        % theta = [theta(length(theta)*2/3:end);theta;theta(1:length(theta)*1/3)];
        % x = [x(length(x)*2/3:end);x;x(1:length(x)*1/3)];
        % y = [y(length(y)*2/3:end);y;y(1:length(y)*1/3)];
        
        clear eq14holder
        
        [period{i},eq14holder] = dowave_greece(theta,deltad,n,x,y,savename,save_on,fetch,i);
        
        
        holdersize= size(eq14holder');
        eq14save{i} = eq14holder';
        
        thisone = size(eq14save)
    end
    
end
