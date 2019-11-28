function [p,g] = TadpoleFinalize(p,g)

% TadpoleFinalize.m
%
% Performs finalization steps for Tadpole

if p.doSaveOutput % if we're saving results

    g.output(:,:,p.lastsave+1)=g.U;

    g.t(p.lastsave+1)=p.t;
    
    g.sealevelsave(p.lastsave+1)=g.sealevel;

    output = g.output;
    
    t = g.t;
    
    sealevel = g.sealevelsave;
    
%     p = rmfield (p,'fighandle');
    
    save(p.runname, '-v7.3', 'p', 't', 'output', 'sealevel');

    % take a snapshot of the final topography and save it as an image
    Snapshot(p,g);
    
else
    
    g.output = g.U;

end

if p.doDrawPlot % if we're plotting the solution
    DrawPlot(p.iterations,p,g)
end
