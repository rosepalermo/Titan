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
addpath('/Users/rosepalermo/Documents/GitHub/VisiLibity1/src')

% variables for the vis code
epsilon = 1e-6;%1e-5;
snap_distance = 0;%0.05;

% load a test lake (wavet2v2) and the ordered shoreline: 
load test_vis_alg.mat

%shoreline
x = fetch_sl_cells{1,1}(:,1);
y = fetch_sl_cells{1,1}(:,2);
% "environment" is the shoreline
environment = {flipud([fetch_sl_cells{1,1};fetch_sl_cells{1,1}(1,:)])};

% find a point eps inside the shoreline, because it won't let me use a
% point on the boundary of the environment to calculate.
x_wrapped = [x; x(1,:)];
y_wrapped = [y; y(1,:)];
diff = [x_wrapped(2:end), y_wrapped(2:end)] - [x, y];
perp = [diff(:,2), -diff(:,1)];
norm = sqrt(diff(:,1).^2 + diff(:,2).^2);
eps = .001;
perp = eps*perp./norm;
x_inside = x + perp(:,1);
y_inside = y + perp(:,2);

% find fetch polygon for each point
for l = 1:10%length(x_inside)
    tic
V{l} = visibility_polygon( [x_inside(l) y_inside(l)] , environment , epsilon , snap_distance );
toc

figure()
%Plot Environment
patch( environment{1}(:,1) , environment{1}(:,2) , 0.1*ones(1,length(environment{1}(:,1)) ) , ...
       'w' , 'linewidth' , 1.5 );
% for i = 2 : size(environment,2)
%     patch( environment{i}(:,1) , environment{i}(:,2) , 0.1*ones(1,length(environment{i}(:,1)) ) , ...
%            'k' , 'EdgeColor' , [0 0 0] , 'FaceColor' , [0.8 0.8 0.8] , 'linewidth' , 1.5 );
% end
hold on
%Plot observer
    plot3( x_inside(l) , y_inside(l) , 0.3 , ...
           'o' , 'Markersize' , 9 , 'MarkerEdgeColor' , 'y' , 'MarkerFaceColor' , 'k' );

%Compute and plot visibility polygon
V{l} = visibility_polygon( [x_inside(l) y_inside(l)] , environment , epsilon , snap_distance );
    patch( V{l}(:,1) , V{l}(:,2) , 0.1*ones( size(V{l},1) , 1 ) , ...
           'r' , 'linewidth' , 1.5 );
    plot3( V{l}(:,1) , V{l}(:,2) , 0.1*ones( size(V{l},1) , 1 ) , ...
           'b*' , 'Markersize' , 5 );
end

%% plot an example (point 100)

figure()
%Plot Environment
patch( environment{1}(:,1) , environment{1}(:,2) , 0.1*ones(1,length(environment{1}(:,1)) ) , ...
       'w' , 'linewidth' , 1.5 );
for i = 2 : size(environment,2)
    patch( environment{i}(:,1) , environment{i}(:,2) , 0.1*ones(1,length(environment{i}(:,1)) ) , ...
           'k' , 'EdgeColor' , [0 0 0] , 'FaceColor' , [0.8 0.8 0.8] , 'linewidth' , 1.5 );
end
hold on
%Plot observer
    plot3( x_inside(100) , y_inside(100) , 0.3 , ...
           'o' , 'Markersize' , 9 , 'MarkerEdgeColor' , 'y' , 'MarkerFaceColor' , 'k' );

%Compute and plot visibility polygon
V{1} = visibility_polygon( [x_inside(100) y_inside(100)] , environment , epsilon , snap_distance );
    patch( V{1}(:,1) , V{1}(:,2) , 0.1*ones( size(V{1},1) , 1 ) , ...
           'r' , 'linewidth' , 1.5 );
    plot3( V{1}(:,1) , V{1}(:,2) , 0.1*ones( size(V{1},1) , 1 ) , ...
           'b*' , 'Markersize' , 5 );