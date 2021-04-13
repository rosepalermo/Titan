% Ligeia Mare calculate wave weighting using visilibity vs fetch

% calculate fetch using 1000 intersecting rays for ligeia Mare
load('Ligeia.mat')
 histogram(fetch)
%%
% rayfetch
[lake,X,Y] = gridlake(x_test',y_test',62.5,62.5,62.5*2);
load('p_init_v1.mat')
p.nrays = 180;
p.delta = 0.05;
p.Ao = 1; p.dt = 1; p.Kwave = 1;
[~,wave_weight_matrix,fetch_matrix,~,~,~] = get_dam_wave_pts(lake,p,X,Y);
%%
load('Ligeia_test.mat')
[indshoreline_ocw,~,~,~] = order_shoreline_bwbound(sebago_lake,p);
load('Ligeia_indshoreline.mat')
load('Ligeia_vis.mat')
ind = sub2ind(size(lake),indshoreline(:,1),indshoreline(:,2));
fetch_vislib = fetch_matrix_vis(ind);
        figure()
    scatter3(indshoreline(:,2),indshoreline(:,1),fetch_vislib,[],fetch_vislib,'filled')
    view(2)
    grid off
    axis tight
    axis equal
    h = colorbar;
    ylabel(h, 'fetch area')
    xlabel('meters');ylabel('meters');
    set(gca,'FontSize',16)


% figure()
% histogram(fetch_matrix)

% visilibity
% [~,wave_weight_matrix_vis,fetch_matrix_vis,~,~,~] = get_dam_wave_pts_vislib(lake,p,X,Y);

