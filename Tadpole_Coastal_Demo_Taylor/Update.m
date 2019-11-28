function [p,g] = Update(p,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE elevations using operator splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANNEL INCISION
if p.doStreamPower
    [p,g] = FT(p,g); % Forward-time explicit
    % Note that FT performs two additional steps:
    % 1. adjusts time step according to the Courant number
    % 2. records locations of channels
    
%     [p,g] = BoundaryMat(p,g); % We don't really need to do this here
%     because points shouldn't erode below SL. Small errors could make that
%     happen (channel points erode below SL), but we can ignore it for now
%     and let SeaLevel() pick them up at the next iteration.
end


% LANDSLIDING
if p.doLandslides
    g = Landslide(p,g);
%     [p,g] = BoundaryMat(p,g); % Landslides shouldn't lower any cells
%     below SL
end


% DIFFUSION
if p.doDiffusion
    if (p.doStreamPower && ~p.doChannelDiffusion) || (p.doStreamPower && p.doAdaptiveTimeStep)
        g = SetUpADI(p,g); % update ADI matrices to exclude channel points
        % in the future, we could do this more efficiently by constructing
        % archetypal ADI matrices at the beginning, and then zeroing out
        % rows or multiplying by dtnew/dtold as necessary, depending on 
        % where the channels are and what the new time step is
    end
    g = ADI(p,g); % alternating-direction implicit (Crank-Nicolson in 2D)
    
    [p,g] = BoundaryMat(p,g);
end

% SOURCE TERMS / PERTURBATIONS 
g = Source(p,g);
% [p,g] = BoundaryMat(p,g); % Subsidence or other source terms could put
% some cells below SL, but as long as we call SeaLevel() next and BoundaryMat() after that, don't need BoundaryMat() here. 

% SEA LEVEL CHANGE
[p,g] = SeaLevel(p,g);
[p,g] = BoundaryMat(p,g); % Update boundary matrices, including g.C, which tracks what is submerged

% COASTAL EROSION
if p.doWaveErosion
    [p,g] = waveerosion(p,g);
    [p,g] = BoundaryMat(p,g);
end

if p.doUniformErosion
    [p,g] = uniformerosion(p,g);
    [p,g] = BoundaryMat(p,g);
end

