function [p,g] = TadpoleRun(p,g)

% TadpoleRun.m
%
% Performs main iteration loop of Tadpole

n=0;
p.lastsave = 0;
while p.t < p.tf || g.nLakeCells < p.Ao_cells*p.size_final % if either the final time is reached or the max lake size-- keeping final time to keep my kill switches
    
    n = n + 1;
    p.n = n;
    %%%%%%%%%%%%%%%%%%%%%% UPDATE ELEVATIONS %%%%%%%%%%%%%%%%%%%%%%%
    
    [p,g] = Update(p,g);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%% INCREMENT TIME %%%%%%%%%%%%%%%%%%%%%%%%%
    
    p.t = p.t + p.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if it is time to redraw the plot, do so.
    if p.doDrawPlot
        if ~rem(n,p.plotint)
            DrawPlot(n,p,g)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if p.doSaveOutput
        %         if ~rem(n,p.saveint)
        if ~p.noSLR
            if p.n == 1 || ((sign(g.SL_slope(p.n)) == -1) && (sign(g.SL_slope(p.n-1)) == 1)) % if first ts iteration or if highstand
                %             p.lastsave = n/p.saveint + 1;
                p.lastsave = p.lastsave+1;
                g.output(:,:,p.lastsave) = g.U;
                g.t(p.lastsave) = p.t;
                g.sealevelsave(p.lastsave) = g.sealevel;
                save(p.runname, '-v7.3', 'p', 'g');
            end
        else
            if p.n == 1 || ~rem(n,p.saveint)
                %             p.lastsave = n/p.saveint + 1;
                p.lastsave = p.lastsave+1;
                g.output(:,:,p.lastsave) = g.U;
                g.t(p.lastsave) = p.t;
                g.sealevelsave(p.lastsave) = g.sealevel;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

p.iterations = n;