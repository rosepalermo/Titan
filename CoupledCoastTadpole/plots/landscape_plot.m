% CoupledModelPlot.m

% file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_15.mat');
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_wave_Kc0_1.mat');
% file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111820/river_mwave_muniform.mat');
load(file_name)
% close(p.fighandle)
v = VideoWriter(file_name(1:end-4));
x = p.dx*(0:p.Nx-1);
y = p.dy*(p.Ny-1:-1:0);
open(v)
i2plot = 248; % saved frames to plot
figure()
for i = 1:i2plot
    for ii = 1
    elev = g.output(:,:,i);
    sl = g.sealevelsave(i);
    time = g.t(i);
    
    % set sea level to elevation 0
    elev = elev - sl;
    
    DEM = GRIDobj(x,y,elev);
%     cmapsea  = [0  0 1;  0.3 0.75 1; 1 1 1];
    cmapsea  = [0  0 1;  0.3 0.75 1];
%     cmapland = flipud([0.93 0.69 0.13; .8 .6 0; 1  1 .8 ]);
%     cmapland = [0 0.4 0.2; 0.87 0.98 0.76; 0.88 0.98 0.76;  0.4  0.17 0 ];
    cmapland = [0.93 0.69 0.13; 1  1 .8 ];
%     [cmap,climits] = demcmap(DEM.Z);
    [cmap,climits] = demcmap(DEM.Z,256,cmapsea,cmapland);

%     [cmap,zlimits] = ttcmap(DEM,'cmap','gmtrelief');    
%     imageschs(DEM,[],'caxis',zlimits,'colormap',cmap)   

    
%     f = figure('WindowStyle','docked');
%     imageschs(DEM,[],'caxis',climits,'colormap',cmap)
    imageschs(DEM,[],'caxis',[climits(1)+1 climits(2)],'colormap',cmap)
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
    colorbar off
%     imageschs(DEM,[],'colorbar',false); % default colormap
%     imageschs(DEM,[],'colormap',[.9 .9 .9],'colorbar',false); % just hillshade
    
    hold on
    contour(x,y,DEM.Z,0*[1,1],'k','linewidth',2)
    hold off

    % stream networks!
    elevtile = repmat(elev,[3,3]);
    elevtile = elevtile(p.Ny/2+1:3*p.Ny-p.Ny/2,p.Nx/2+1:3*p.Nx-p.Nx/2);
    snks = elev < 0;
    snkstile = repmat(snks,[3,3]); % logical matrix with sinks = 1
    snkstile = snkstile(p.Ny/2+1:3*p.Ny-p.Ny/2,p.Nx/2+1:3*p.Nx-p.Nx/2);
    [nr,nc] = size(elevtile);
    xtile = p.dx*(0:nc-1);
    ytile = p.dy*(nr-1:-1:0);
    elevtile = GRIDobj(xtile,ytile,elevtile);
    snkstile = GRIDobj(xtile,ytile,snkstile); % make it a GRIDobj
    FD = FLOWobj(elevtile,'preprocess','carve','sinks',snkstile,'mex',true);
    Amincells = 500; % minimum drainage area in cells
    S = STREAMobj(FD,'minarea',Amincells);
    S = modify(S,'upstreamto',elevtile>=0); % only want rivers above SL. Really this should be >=0, but the way TopoToolbox trims the network makes the outlets stop short of the coast.
    % now we need to crop or at least shift the stream network
    S.x = S.x - p.dx*p.Nx/2;
    S.y = S.y - p.dy*p.Ny/2;
    
    hold on
%     plot(S,'color',[.5 .5 .5],'linewidth',2)
    plot(S,'k','linewidth',2)
    hold off
    
    set(gca,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)])
    
%     savefig([runname '_t' num2str(time) '_rivers.fig'])
%     print([runname '_t' num2str(time) '_rivers'],'-dpng','-r300');
        frame = getframe(gcf);
    writeVideo(v,frame);
%     close(f)
    end
end