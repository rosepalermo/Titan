% old corners script. No longer need to do this because new definition of
% shoreline includes the corners

% %         SHOULDNT NEED TO DO THIS ANYMORE BECAUSE NEW ORDER THE SHORELINE
% %         INCLUDES CORNERS
% %         % find corners and damage if they exist alone (not in the cells to
% %         % trash)
%         corners = setdiff(find(shoreline),indshoreline);
%         if ~isempty(cells2trash)
%             [c2t]=sub2ind(size(X),cells2trash(:,1),cells2trash(:,2)); % find the cells that arent the trash cells
%             if exist('c2t') & ~isempty(corners)
%                 corners = corners(~ismember(corners,c2t));
%             end
%         end
% %         
% %         % find the mean of the damage for points next to the corners. make that
% %         % the damage for that corner
%         if ~isempty(corners)
%             corners = corners(shoreline(corners)<1.5); % if less than 1.5, only a corner. not also a side
%             [damcorn] = damagecorners(lake,corners,indshoreline,dam);
%             %         strength(indshoreline) = strength(indshoreline) - ones(length(indshoreline),1).*dam;
%             strength(corners) = strength(corners) - p.dt*p.Kcoast*shoreline(corners).*damcorn;
%         end

% %   Find the corners and change the damage to sum corners* sqrt2/2 * wave
% %   weighting
%         [sl_nocorners] = addidshoreline_cardonly(F_lake,~F_lake);
%         corners = setdiff(find(shoreline),find(sl_nocorners));
%         cornind = ismember(indshoreline,corners);
% %         dam(cornind) = dam(cornind).*shoreline(indshoreline(cornind));
%         dam = dam.*shoreline(indshoreline);
%  

%             erodedind_12 = cells2trash(find(~ismember(cells2trash,corners)));
