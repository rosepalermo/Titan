% example_Titan.m
%
% sample script for running Tadpole


clear

folder = fileparts(which('example_Titan_coupled.m'));
addpath(genpath(folder));

%% SET PARAMETERS %%

% --------------- space and time resolution ------------------------------- 

p.Nx = 400;                 %     p.Nx             Number of grid points in x direction
p.Ny = 400;                 %     p.Ny             Number of grid points in y direction
p.dx = 125/2;                 %     p.dx             Grid spacing in the x direction (m)
p.dy = 125/2;                 %     p.dy             Grid spacing in the y direction (m)

p.doAdaptiveTimeStep = 1;   % p.doAdaptiveTimeStep Turn adaptive time step based on Courant number on (1) or off (0). If set to off, time step is p.dtmax
p.doAdaptiveCoastalTimeStep = 1;
p.dtmax = 300;%1e4;              %     p.dtmax          maximum time step (yr)
p.Courant = 0.9;            %     p.Courant        maximum Courant number

p.tf = 3e10;                 %     p.tf             Total time of the simulation (yr)


% ----- boundary conditions, source terms, and flow routing ---------------

p.bdy.left  = 'periodic';   %     p.bdy            a struct with fields 'left' 'right' 'lower' 'upper'
p.bdy.right = 'periodic';   %                      specifying boundary condition:
p.bdy.upper = 'periodic';   % 
p.bdy.lower = 'periodic';   %                      'fixed'    --> constant elevation (derivatives set equal to zero)
                            %                      'mirror'   --> "mirror" boundary (centered differencing using boundary
                            %                                     point and one interior point)
                            %                      'periodic' --> periodic boundary (centered differencing using boundary
                            %                                     point, one interior point, and one
                            %                                     point from opposite boundary)

p.E = 1e-10;                %     p.E              Rate of surface uplift or base level lowering (m/yr)
p.routing = 'Dms';          %     p.routing        Choose which flow routing method to use: 'D8' (steepest descent), 'Dinf' (Tarboton's D-infinity), or 'Dms' (multi-slope)
p.flood = 1;                %     p.flood          1=route flow through local minima, 0=don't

p.F = zeros(p.Ny,p.Nx);     %     p.F              Optional matrix of fixed points, in addition to boundary conditions above. Points with p.F == 1 will have constant elevation
                            
% ------------------ plotting and output ----------------------------------

p.doDrawPlot = 0;           %     p.doDrawPlot     Display solution as the run progresses
p.plotint = 1;%100;            %     p.plotint        Plot will be redrawn every plotint iterations
p.plottype = 'elevation';             %     p.plottype       1=perspective view, 2=drainage area map, 3=curvature map, 4=elevation map, 5=contour map, 6=shaded relief, 7=colored shaded relief
                            %
p.doSaveOutput = 0;         %     p.SaveOutput     Save model output to a .mat file
p.saveint = 5; %1000;              %     p.saveint        Elevation grid will be saved every saveint iterations
p.runname = 'trash';        %     p.runname:       Character string naming the run. If specified 
                            %                      (and if p.saveint~=0), the parameters and elevations at each 
                            %                      save interal will be saved in a binary .MAT file called <runname>.mat
                           
 
                            
% ------------------ hillslope processes ----------------------------------                            
                                                    
p.doDiffusion = 0;          %     p.doDiffusion    Turn hillslope diffusion on (1) or off (0)
p.D = 0.005;                %     p.D              Hillslope diffusivity (m^2/yr)
                            %
p.doLandslides = 0;         %     p.doLandslides   Turn landslides on (1) or off (0)
p.Sc = 0.6;                 %     p.Sc             Critical slope (m/m)


% ---------------- bedrock channel incision -------------------------------                           

p.doStreamPower = 0;        %     p.doStreamPower  Turn bedrock channel incision on (1) or off (0)
p.doChannelDiffusion = 0;   %     p.doChannelDiffusion Turn diffusion in channels on (1) or off (0)
p.Kf = 1e-5; % 5e-6;                %     p.Kf             Coefficient in stream power incision law (kg m^(1+2m) yr^-2)
p.m = 0.5;                  %     p.m              Drainage area exponent in stream power law
p.w = 1.0;                  %     p.w              Slope exponent in stream power law
p.Kw = 0.5*p.dx;                 %     p.Kw             Coefficient relating channel width to drainage area
p.wexp = 0;                 %     p.wexp           Exponent relating channel width to drainage area
p.thetac = 0;               %     p.thetac         Threshold for fluvial incision

% ---------------- coastal erosion -------------------------------                           

p.doWaveErosion = 0;        %     p.doWaveErosion  Turn fetch based coastal erosion on (1) or off (0)
p.doUniformErosion = 1;     %     p.doUniformErosion  Turn uniform coastal erosion on (1) or off (0)
% p.SLR = 50/p.tf;                  %     p.SLR            Rate of sea level rise (m/yr)
p.sealevel_init = 1;        %     p.sealevel_init  Initial sea level
if p.doUniformErosion
    p.strength = 1;         %     p.strength       Initial strength of the bedrock
    p.Kcoast = 1e-4;        %     p.Kcoast         Coastal erosion rate constant (damage * strength^-1 * yr^-1)
elseif p.doWaveErosion

    p.strength = 1;         %     p.strength       Initial strength of the bedrock
    p.Kcoast = 1.3e-12;     %1/p.dtmax/((p.Nx*p.dx).^2)/4; % maximum damage that could occur on an island that sees the whole domain
    %               5.5e-8;       %     p.Kcoast         Coastal erosion rate constant (damage * strength^-1 * yr^-1)
        % the range between 4e-8 and 5.5e-8 was hard to choose from. I
        % chose 4e-8 because it produced a reasonable range (showing difference 
        % between embayments and open coast) at 100m sea level (our amplitude). 
        % Less erosion will occur below that because less of the coastline
        % will be eroded total, but still a difference between embayment
        % and open coast

else
    p.strength = 0;
end

% no sea level change 1, sinusoidal sea level change 0
p.noSLR = 1;

% ------------------ initial conditions -----------------------------------                           

p.beta = 1.6; % 1.6;             %     p.beta           Negative slope of the power spectrum. 0 = white noise, more positive values are "redder" (more variance at longer wavelengths)
p.variance = 10000;         %     p.variance       Variance of elevation (m^2)
p.periodic = 1;             %     p.periodic       Elevations will be periodic at the boundaries (1) or not (0, default)


%% CREATE INITIAL SURFACE %%

% noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);

% initgaus = get_gaussian_boundary([800 800], 0.3, 10);
% init = (initgaus + noise);
rfactor = 0.25; % 0.25; % depth of the depression as a function of relief of the noise surface
init = get_IC(p,rfactor);

% adjust the elevations so pctwet % of the domain is below initial SL
pctwet = 10;
Zshift = prctile(init(:),pctwet);
init = init - Zshift + p.sealevel_init;


% set fixed points according to initial sea level. Note that p.F will need
% to be updated each time step according to changes in elevations,
% coastal positions, and sea level. --> No, that is what g.C is for (that's
% why it's a "grid" in g and p.F is a parameter in F). 

% p.F(init < p.sealevel_init) = 1; % I forget if you decided that points with elevations equal to SL would be considered land or submerged. Here I assumed they are land; if submerged, this line should be <= instead of <

%test circle
[init] = 10*test_circle(p);
p.F = zeros(size(init));
% p.F(init < p.sealevel_init) = 1; % I forget if you decided that points with elevations equal to SL would be considered land or submerged. Here I assumed they are land; if submerged, this line should be <= instead of <

% make lowest 10% of elevations fixed points
% SL = prctile(init(:),10);
% init = init - SL;
% % % init(init < 0) = 0; 
% p.F(init <= 0) = 1;

%% RUN THE MODEL %%

% run the model, storing the final elevation grid in solution
Kf_ = [5e-10 5e-8 5e-6];% 5e-5 5e-7];%[5e-6 4.5e-6 4e-6 5.5e-6 6e-6 6.5e-6];
folder = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/River_and_wave_3_20/';
run = 'wave_river_tf1e5_Kf_';
for i = 1:1%:length(Kf_)
    p.Kf = Kf_(i);
    kf__ = num2str(Kf_(i));
    p.runname = strcat(folder,run,kf__);
    solution = Tadpole(init,p);
end
% save('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/ModelingAGU19/uniform1.mat',solution,p,init)
