function [P,Pobs,nbi] = island_area(P0,Pobs,l_ind,eps)
P = zeros(size(P0));

% make the island epsilon bigger
for l = 1:length(P0)
    if l > 1
        l_ind_prev = l-1;
    else
        l_ind_prev = length(P0);
    end
    
    if l < length(P0)
        l_ind_next = l+1;
    else
        l_ind_next = 1;
    end
    
    xp = P0((l_ind_prev), 1);
    yp = P0((l_ind_prev), 2);
    xn = P0((l_ind_next), 1);
    yn = P0((l_ind_next), 2);
    [P(l,:),~] = Pint([xn,yn],[P0(l,1),P0(l,2)],[xp,yp],eps);
   
end

% find the new observer point for the slightly bigger island if i == k
if size(Pobs)>0
l = l_ind;
    if l > 1
        l_ind_prev = l-1;
    else
        l_ind_prev = length(P);
    end
    
    if l < length(P)
        l_ind_next = l+1;
    else
        l_ind_next = 1;
    end
    
    xp = P((l_ind_prev), 1);
    yp = P((l_ind_prev), 2);
    xn = P((l_ind_next), 1);
    yn = P((l_ind_next), 2);
[Pobs,nbi] = Pint([xn,yn],[P(l,1),P(l,2)],[xp,yp],eps/2);
end
end