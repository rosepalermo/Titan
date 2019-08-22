function [Pobs,env,xsub,ysub,nbi] = SubsampleShorelineislands(slccw,k,l,binfine,bincoarse,dist,nhood,eps,epsilon)

x_l = slccw{k}(l,1);
y_l = slccw{k}(l,2);

% The size with which to increase the distance each iteration, this
% parameter might need tuning.
delta_dist = 1.0;

binclose = cell(length(slccw), 1);
binsimp = cell(length(slccw), 1);
ind = cell(length(slccw), 1);
env = cell(length(slccw), 1);

% Precompute the squared distance
dsq = cell(length(slccw), 1);
for i=1:length(slccw)
    dsq{i} = (slccw{i}(:,1) - x_l).^2 + (slccw{i}(:,2) - y_l).^2;
end

% Increase the distance until all of the points are contained within it.
observer_is_valid = 0;
while ~observer_is_valid
    % Update binclose and binsimp for all of the bodies.
    for i=1:length(slccw)
        % Populate binclose
        binclose{i} = dsq{i}' < dist^2;
        
        binsimp{i} = binclose{i} + bincoarse{i}.*(~binclose{i});
        ind{i} = find(binsimp{i});
    end
    
    % Get the next and previous points.
    l_ind = find(ind{k} == l);
    if l_ind > 1
        l_ind_prev = l_ind-1;
    else
        l_ind_prev = length(ind{k});
    end
    
    if l_ind < length(ind{k})
        l_ind_next = l_ind+1;
    else
        l_ind_next = 1;
    end
    
    xp = slccw{k}(ind{k}(l_ind_prev), 1);
    yp = slccw{k}(ind{k}(l_ind_prev), 2);
    xn = slccw{k}(ind{k}(l_ind_next), 1);
    yn = slccw{k}(ind{k}(l_ind_next), 2);

    % Get the observer position.
    if k == 1 % this is the lake
        [Pobs,nbi] = Pint([xn,yn],[x_l,y_l],[xp,yp],eps);
    else % this is an island
        [Pobs,nbi] = Pint([xp,yp],[x_l,y_l],[xn,yn],eps);
    end
    
    % Populate the environment
    for i=1:length(slccw)
        xsub = slccw{i}(ind{i},1);
        ysub = slccw{i}(ind{i},2);
        nsub = length(ind{i});
        
        % if the polygon we're working on is the one that has l, make l
        % the first point, if not, just use the ordered subsampled xy
        if i == k
            env{i} = flipud([xsub([l_ind:nsub, 1:l_ind-1]),ysub([l_ind:nsub, 1:l_ind-1])]);
        else
            env{i} = [xsub,ysub];
        end
    end
    
    % Check if the point is inside the environment.
    if in_environment(Pobs, env, epsilon)
        % If the selected point is in the environment, we're done.
        observer_is_valid = 1;
    else        
        % Increase the distance for the next iteration.
        dist = dist + delta_dist;
    end
end
