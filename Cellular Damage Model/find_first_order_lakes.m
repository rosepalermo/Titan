function [L,fol] = find_first_order_lakes(lake)


% load('test_polygon_island.mat')
% lake = repmat(polygon,2);
% load('wavet2v2_lake_025.mat')

% figure()
[bound, L, n, A]  = bwboundaries(lake);
% imagesc(lake');
% hold on;
for i = 1:length(bound)
    scatter(bound{i,1}(:,1),bound{i,1}(:,2));
end

n = length(bound);
first_order_lakes = zeros(n,1);
for i=1:n
    if ~any(A(i,:))
        first_order_lakes(i) = 1;
    end
end

first_order_islands = zeros(n,1);
for i=find(first_order_lakes)'
    first_order_islands = first_order_islands | A(:,i);
end

second_order_lakes = zeros(n,1);
for i=find(first_order_islands)'
    second_order_lakes = second_order_lakes | A(:,i);
end


fol = find(first_order_lakes);
% imagesc(L == fol(1))

% 
% enclosing_boundary  = find(A(6,:));
% enclosed_boundaries = find(A(:,6));

end
