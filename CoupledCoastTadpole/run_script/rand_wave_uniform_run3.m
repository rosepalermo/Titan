% example_Titan.m
%
% sample script for running Tadpole


clear

cd ..

folder = fileparts(which('getpath_CCT.m'));
addpath(genpath(folder));

init_circle =0;
init_square = 0;
rand_IC = 1;
river_IC =0;
%% SET PARAMETERS %%

% --------------- space and time resolution ------------------------------- 

p.Nx = 400;                 %     p.Nx             Number of grid points in x direction
p.Ny = 400;                 %     p.Ny             Number of grid points in y direction
p.dx = 125/2;                 %     p.dx             Grid spacing in the x direction (m)
p.dy = 125/2;                 %     p.dy             Grid spacing in the y direction (m)

p.doAdaptiveTimeStep = 1;   % p.doAdaptiveTimeStep Turn adaptive time step based on Courant number on (1) or off (0). If set to off, time step is p.dtmax
p.doAdaptiveCoastalTimeStep = 1;
p.dtmax = 100;%1e4;              %     p.dtmax          maximum time step (yr)
p.Courant = 0.9;            %     p.Courant        maximum Courant number

% p.tf = 1e5;                 %     p.tf             Total time of the simulation (yr)
p.size_final = 1.2;

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
p.doSaveOutput = 1;         %     p.SaveOutput     Save model output to a .mat file
p.saveint = 1;              %     p.saveint        Elevation grid will be saved every saveint iterations
% p.runname = 'trash';        %     p.runname:       Character string naming the run. If specified 
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

p.doWaveErosion = 1;        %     p.doWaveErosion  Turn fetch based coastal erosion on (1) or off (0)
p.doUniformErosion = 1;     %     p.doUniformErosion  Turn uniform coastal erosion on (1) or off (0)
p.So = 1;
p.dxo = 100;
% p.SLR = 50/p.tf;                  %     p.SLR            Rate of sea level rise (m/yr)
p.sealevel_init = 1;        %     p.sealevel_init  Initial sea level
if p.doUniformErosion
    p.strength = p.So*p.dx/p.dxo;         %     p.strength       Initial strength of the bedrock
    p.Kuniform = 1e-4;        %     p.Kcoast         Coastal erosion rate constant (damage * strength^-1 * yr^-1)
elseif p.doWaveErosion

    p.strength = p.So*p.dx/p.dxo;         %     p.strength       Initial strength of the bedrock
    p.Kwave = 1.3e-12;     %1/p.dtmax/((p.Nx*p.dx).^2)/4; % maximum damage that could occur on an island that sees the whole domain
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
if rand_IC
    % initgaus = get_gaussian_boundary([800 800], 0.3, 10);
    % init = (initgaus + noise);
    rfactor = 0.25; % 0.25; % depth of the depression as a function of relief of the noise surface
    [init,p] = get_IC(p,rfactor);
elseif river_IC
    load('riverIC.mat')
    init = riverIC;
    p.Ao = 8.9298e+07;  
    p.Ao_cells = 30368;
%test circle
elseif init_circle
    [init,p] = test_circle(p);
elseif init_square
    [init,p] = test_square(p);
end

% % set fixed points
p.F = zeros(size(init));
p.F(init < p.sealevel_init) = 1; % I forget if you decided that points with elevations equal to SL would be considered land or submerged. Here I assumed they are land; if submerged, this line should be <= instead of <

%% RUN THE MODEL %%

% run the model, storing the final elevation grid in solution
% Kf_ = [5e-10 5e-8 5e-6]; % rivers
Kc_ = [1e-4 1.5e-4 1.5e-3]; % uniform/wave
% Kc_ = [1e-4; 1e-3; 1e-2];% 5e-13 2e-13];
% Kc_ = 1.5e-8*1000;%for circle, wave p.Kcoast = 1.5e-8
% Kc_ = p.Kf;
% p.folder = '/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/riverIC/';
p.folder = '/home/rpalermo/TitanModelOutput/08_2020/results1/save_more/';
p.run = 'rand_mwave_muniform';
time = 'time';
for i = 1
    p.tf = 1e5;
    tstart = tic;
    if p.doUniformErosion
        p.Kuniform = Kc_(2);
    end
    if p.doWaveErosion
        p.Kwave = Kc_(2);
    end
    p.runname = strcat(p.folder,p.run);
    solution = Tadpole(init,p);
    time_end(i) = toc(tstart);
    timename = strcat(p.folder,p.run,p.run2,time);
    save(timename,'time_end','p','i');
end
