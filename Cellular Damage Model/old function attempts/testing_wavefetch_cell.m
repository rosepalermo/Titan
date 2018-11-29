% An example polygon with lots of concavity
th = linspace(0,2*pi,60);
r = 2 + rand(size(th))-0.5;
x_ex = r.*cos(th);
y_ex = r.*sin(th);

shoreline_fetch = cell(1,1);
shoreline_fetch{1,1}(:,1) = x_ex;
shoreline_fetch{1,1}(:,2) = y_ex;
shoreline_fetch{2,1}(:,1) = [0; 0.5; 1; 0.5; 0];
shoreline_fetch{2,1}(:,2) = [0; 0; 0; 1; 0];
shoreline_fetch{3,1}(:,1) = [0; -0.5; -1; -0.5; 0];
shoreline_fetch{3,1}(:,2) = [0; 0; 0; -1; 0];

% [testwave,~] = fetch_wavefield_cell(shoreline_fetch);
% 
% figure()
% plot(x_ex,y_ex,'k')
% hold on
% scatter3(x_ex,y_ex,testwave{1,1},[],testwave{1,1})
% view(2)
% scatter3(shoreline_fetch{2,1}(:,1),shoreline_fetch{2,1}(:,2),testwave{2,1},[],testwave{2,1})
% plot(shoreline_fetch{2,1}(:,1),shoreline_fetch{2,1}(:,2),'k')
% scatter3(shoreline_fetch{3,1}(:,1),shoreline_fetch{3,1}(:,2),testwave{3,1},[],testwave{3,1})
% plot(shoreline_fetch{3,1}(:,1),shoreline_fetch{3,1}(:,2),'k')

lakex = [x_ex'; shoreline_fetch{2,1}(:,1); shoreline_fetch{3,1}(:,1)];
lakey = [y_ex'; shoreline_fetch{2,1}(:,2); shoreline_fetch{3,1}(:,2)];
eps = 2;
dx = 0.01; dy = 0.01
modelrun = 1; fetch_on = 1; savename = 'test';
cdm_Titan(lakex,lakey,eps,dx,dy,modelrun,fetch_on,savename)