function [Pobs,env,xsub,ysub,nbi] = SubsampleShoreline(slccw,k,l,binfine,bincoarse,dist,nhood,eps,epsilon)

x_l = slccw{k}(l,1);
y_l = slccw{k}(l,2);
np = length(slccw{k}(:,1));

for i = 1:length(slccw)
    
    x = slccw{i}(:,1);
    y = slccw{i}(:,2);
    
%     if length(x)>20
        
        np_x = length(x);
        % compute squared euclidean distances to every other shoreline point
        dx = x_l-x;
        dy = y_l-y;
        dsq = dx.*dx + dy.*dy;
        
        % include points that are within a distance dist to point l.
        binclose = zeros(1,np_x);
        binclose(dsq <= dist*dist) = 1;
        
        % now we just need to know which elements to draw from fine, and which from
        % coarse
        % Only subsample if the coastline is > 20 (i.e. don't subsample small islands or lakes)
        if length(x) > 20
            binsimp = binfine{i}.*binclose + bincoarse{i}.*(~binclose);
        else
            binsimp = ones(size(binfine{i}));
        end
        observer_is_valid = 0;
        
        while ~observer_is_valid
            
            % Add in point l and the specified neighborhood
            neighb = (l-nhood):(l+nhood);
           
            
            % periodic bc's
            neighb(neighb<1) = neighb(neighb<1)+np;
            neighb(neighb>np) = neighb(neighb>np)-np;
            
            if length(x) > 20
                binsimp(neighb) = 1;
            end
            
            % subsample
            ind = find(binsimp);
            lind = find(ind==l); % what is the index of point l in binsimp?
            xsub{i} = x(ind);
            ysub{i} = y(ind);
            
            
            % Find the interior observer point, which lies a distance eps inside
            % the shoreline along the mid-vector defined by the unit normals to the
            % two sides adjacent to point l (the mid-vector bisects the interior
            % angle made by the two adjacent sides)
            % ROSE: WILL ISLANDS REQUIRE SPECIAL HANDLING TO MAKE SURE THE OBSERVER
            % POINT IS IN THE LIQUID?
            
            nsub = length(xsub{i});
            
            if lind==1
                lprev = nsub;
            else
                lprev = lind-1;
            end
            
            if lind==nsub
                lnext = 1;
            else
                lnext = lind+1;
            end
            
            xp = xsub{i}(lprev); % previous point
            yp = ysub{i}(lprev);
            x0 = xsub{i}(lind); % point l
            y0 = ysub{i}(lind);
            xn = xsub{i}(lnext); % next point
            yn = ysub{i}(lnext);
            
            if i == k % only making 1 observer. other closed loops are islands
                if k == 1 % this is the lake
                    [Pobs,nbi] = Pint([xn,yn],[x0,y0],[xp,yp],eps);
                else % this is an island
                    [Pobs,nbi] = Pint([xp,yp],[x0,y0],[xn,yn],eps);
                end
            end
            
            % Visilibity wants the "environment" to be a CCW-ordered polygon
            
            % I found that there are some errors if the first point == last point.
            % Better to not repeat this point.
            % But either way, for some reason the polygon isn't quite right unless it's sorted to put point l first. Can't figure out why.
            %         env{i} = [xsub{i}([lind:nsub, 1:lind-1]),ysub{i}([lind:nsub, 1:lind-1])];
            
            % if the polygon we're working on is the one that has l, make l
            % the first point, if not, just use the ordered subsampled xy
            if i == k
                env{i} = flipud([xsub{i}([lind:nsub, 1:lind-1]),ysub{i}([lind:nsub, 1:lind-1])]);
            else
                env{i} = [xsub{i},ysub{i}];
            end
            
            if i == k % if we're on the polygon that the observer point is on
                if in_environment(Pobs, env, epsilon) % If the selected point is in the environment, we're done.
                    observer_is_valid = 1;
                else % If not, we probably have a polygon that intersects itself. Increase nhood and try again until the observer is in the environment.
                    nhood = nhood+1;
                    %         disp(['Iterating on point' num2str(l) '. nhood=' num2str(nhood)])
                end
            else
                observer_is_valid = 1;
            end
        end
%     else
%         % Do the same thing as above but don't actually do any sub
%         xsub{i} = x;
%         ysub{i} = y;
%         nsub = length(xsub{i});
% 
%         if lind==1
%             lprev = nsub;
%         else
%             lprev = lind-1;
%         end
% 
%         if lind==nsub
%             lnext = 1;
%         else
%             lnext = lind+1;
%         end
% 
%         xp = xsub{i}(lprev); % previous point
%         yp = ysub{i}(lprev);
%         x0 = xsub{i}(lind); % point l
%         y0 = ysub{i}(lind);
%         xn = xsub{i}(lnext); % next point
%         yn = ysub{i}(lnext);
% 
%         if i == k % only making 1 observer. other closed loops are islands
%             if k == 1 % this is the lake
%                 Pobs = Pint([xp,yp],[x0,y0],[xn,yn],eps);
%             else % this is an island
%                 Pobs = Pint([xn,yn],[x0,y0],[xp,yp],eps);
%             end
%         end
% 
%         % Visilibity wants the "environment" to be a CCW-ordered polygon
% 
%         % I found that there are some errors if the first point == last point.
%         % Better to not repeat this point.
%         % But either way, for some reason the polygon isn't quite right unless it's sorted to put point l first. Can't figure out why.
%         %         env{i} = [xsub{i}([lind:nsub, 1:lind-1]),ysub{i}([lind:nsub, 1:lind-1])];
% 
%         % if the polygon we're working on is the one that has l, make l
%         % the first point, if not, just use the ordered subsampled xy
%         if i == k
%             env{i} = flipud([xsub{i}([lind:nsub, 1:lind-1]),ysub{i}([lind:nsub, 1:lind-1])]);
%         else
%             env{i} = [xsub{i},ysub{i}];
%         end
        
%     end
    
end