% Visibility algorithm/find fetch demo
% Rose Palermo
% This demo shows a calculation of the fetch around a lake using the
% visilibity1 algorithm

%% before running the demo, need to mex things in the src folder
% mex -setup
% mex -v in_environment.cpp visilibity.o
% mex -v shortest_path.cpp visilibity.o
% mex -v visibility_graph.cpp visilibity.o
% mex -v visibility_polygon.cpp visilibity.o
%%
%add your path to the visibility folder here:
% addpath('/Users/rosepalermo/Documents/GitHub/VisiLibity1/src')

% clear
%
% % load a test lake (wavet2v2) and the ordered shoreline:
% load test_vis_alg.mat fetch_sl_cells
function [WaveArea,FetchArea] = fetch_vis_approx(fetch_sl_cells);

slccw = fetch_sl_cells;
% for i = 1:length(slccw)
%     slccw{i} = flipud[fetch_sl_cells{i}(:,:);fetch_sl_cells{i}(1,:)];
% end
% slccw{1} = flipud(fetch_sl_cells{1});
np = 0;
x_all = [];
y_all = [];
% parameters for shoreline simplification
tolfine = 0.0001; % 0.0001 works well
tolcoarse = 0.005; % 0.005 works well
nhood = 2; % This is a tradeoff between choosing a larger value that is more likely to work everywhere vs. iterating more times at a few points. It seems to have only a secondary effect on execution time.
dist = 20;%0.1*max([range(x) range(y)]); % Typical value is ~10-100 in units of x

% parameters for visibility functions
eps = 0.1;
epsilon = 1e-4;
snap_distance = 0.05; % Looks like we get some invalid (empty array) visibility polygons if snap_distance >= eps.
% disp('simplify')
parfor ii = 1:length(fetch_sl_cells)
    % shoreline. x and y are ordered clockwise, first point != last point
    x = slccw{ii}(:,1);
    y = slccw{ii}(:,2);
%     np = np + length(x);
%     x_all = [x_all;x];
%     y_all = [y_all;y];
        
    % Approximate the shoreline using the Douglas-Peucker algorithm
    [binfine{ii},bincoarse{ii}] = SimplifyShoreline(x,y,tolfine,tolcoarse); % binary arrays where 1 indicates that a point lies on a fine or coarse approximations of the shoreline
    
end

% find fetch polygon for each point
WaveArea = cell(length(slccw),1);
tic
for k = 1:length(slccw)
    WaveArea{k,1} = zeros(length(slccw{k}),1);
    for l = 1:length(slccw{k})
%             l
        % For each point, subsample the shoreline at a desired interval.
%         disp('subsample')
        [Pobs,env,~,~,nbi] = SubsampleShorelineislands(slccw,k,l,binfine,bincoarse,dist,nhood,eps,epsilon);
%         disp('vis')
        
        V = visibility_polygon(Pobs, env, epsilon, snap_distance);
        

%         disp('vis calcl')
        % fetch & wave area & distance & cos(theta-phi)
        FetchArea{k,1}(l,1) = polyarea(V(:,1),V(:,2));
        Fetch_dist = sqrt(sum(([slccw{k}(l,1),slccw{k}(l,2)] - V).^2,2));
        minFetch_dist = min(Fetch_dist,200); % this is how we can limit fetch eventually
        % Wave weighting = (F)*cosang
        weighted_fetch_dist = ([slccw{k}(l,1),slccw{k}(l,2)] - V)*nbi'; % magnitude of fetch * magnitude of normal vector * cosang
        % cos(theta - phi) = dot product of slvec and losvec
        cosang = -weighted_fetch_dist./Fetch_dist; %[mag_fetch*mag_norm(which is 1)*cosang]/mag_fetch = cosang
        cosang(isnan(cosang)) = 0;
        Wavepts = [slccw{k}(l,1),slccw{k}(l,2)]+(V-[slccw{k}(l,1),slccw{k}(l,2)]).*cosang;
        WaveArea{k,1}(l,1) = polyarea(Wavepts(:,1),Wavepts(:,2));
        
    end
end

runtime = toc;
% disp(['dist=' num2str(dist) ' took ' num2str(runtime) ' sec'])

%% Calculate fetch polygon areas and compare with exact calculation
% NEED TO UPDATED EXACT CALCULATION TO INCLUDE ISLANDS

% Areas_approx = zeros(np,1);
%
% tic
% for l = 1:np
%
%     xs = V{k,l}(:,1);
%     ys = V{k,l}(:,2);
%     Areas_approx(l) = polyarea([xs;xs(1)],[ys;ys(1)]);
%
% end
% toc
%
% % figure
% % hist(log10(AreasAppx))
%
% load('Fetch_exact.mat', 'Areas')
% Areas_exact = Areas;
% Amin=min(Areas_exact);
% Amax=max(Areas_exact);
%
% figure
% loglog(Areas_exact,Areas_approx,'.k')
% hold on
% set(gca,'tickdir','out')
% xlabel('Exact Fetch Area' )
% ylabel('Approximated Fetch Area' )
% title(['tfine=' num2str(tolfine) ', tcoarse=' num2str(tolcoarse) ', nhood=' num2str(nhood) ', runtime=' num2str(runtime) 's'])
% set(gca,'tickdir','out')
% axis equal square
%
% plot([Amin Amax],2*[Amin Amax],'r')
% plot([Amin Amax],[Amin Amax]/2,'r')
% plot([Amin Amax],1.5*[Amin Amax],'b')
% plot([Amin Amax],[Amin Amax]/1.5,'b')
% plot([Amin Amax],[Amin Amax]/1.1,'m')
% plot([Amin Amax],1.1*[Amin Amax],'m')

% %% plot an example
%
% pt = 100; % choose a point
%
% % Plot the shoreline
% figure
% plot(x,y,'-k')
% axis equal
% alpha 0.5
% hold on
%
%
% % Compute and plot the actual visibility polygon
% Pobs = SubsampleShoreline(x,y,pt,binfine,bincoarse,dist,nhood,eps,epsilon);
% environment = {flipud([fetch_sl_cells{1,1};fetch_sl_cells{1,1}(1,:)])};
%
% Vtest{1} = visibility_polygon(Pobs, environment, epsilon, snap_distance);
% patch( Vtest{1}(:,1) , Vtest{1}(:,2) , [0.4 0.4 0.4] , 'FaceAlpha',0.5 );
%
% %Plot observer
% plot( Pobs(1) , Pobs(2) , 'o' , 'Markersize' , 9 , 'MarkerEdgeColor' , 'k' , 'MarkerFaceColor' , [0.4 0.4 0.4] );
%
%
% % Compute and plot the subsampled visibility polygon
% [Pobs,environment,xsub,ysub] = SubsampleShoreline(x,y,pt,binfine,bincoarse,dist,nhood,eps,epsilon);
% plot([xsub;xsub(1)],[ysub;ysub(1)],'-r')
%
% Vtest{2} = visibility_polygon(Pobs, environment, epsilon, snap_distance);
% patch( Vtest{2}(:,1) , Vtest{2}(:,2) , 'r' ,'FaceAlpha',0.5 );
%
% %Plot observer
% plot( Pobs(1) , Pobs(2) , 's' , 'Markersize' , 9 , 'MarkerEdgeColor' , 'k' , 'MarkerFaceColor' , 'r' );
end