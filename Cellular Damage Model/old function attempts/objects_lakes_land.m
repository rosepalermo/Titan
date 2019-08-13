% load('test_polygon_island.mat')
% lake = repmat(polygon,2);
% load('wavet2v2_lake_025.mat')
load('test_lake_objects.mat')
% lake = L;
figure()
[bound, L, ~, A]  = bwboundaries(lake,8);
n = length(bound);
imagesc(lake');
hold on;
for i = 1:n
    plot(bound{i,1}(:,1),bound{i,1}(:,2));
end

first_order_lakes = zeros(n,1);
for i=1:n
    if ~any(A(i,:))
        first_order_lakes(i) = i;
    end
end

first_order_islands = zeros(n,1);
islands_temp = zeros(n,1);
for i=find(first_order_lakes)'
    islands_temp = zeros(size(first_order_lakes)) | A(:,i);
    first_order_islands(islands_temp) = i;
    clearvars islands_temp
%     first_order_islands = first_order_islands | A(:,i);
%     first_order_islands(:,2) = first_order_islands(:,1)*i;
end

second_order_lakes = zeros(n,1);
for i=find(first_order_islands)'
    sol_temp = zeros(size(second_order_lakes)) | A(:,i);
    second_order_lakes_i(sol_temp) = find(sol_temp);
    second_order_lakes(sol_temp) = i;
    clearvars sol_temp
end


second_order_islands = zeros(n,1);
for i=find(second_order_islands)'
    second_order_islands = second_order_islands | A(:,i);
end

fol = first_order_lakes(find(first_order_lakes>0));
foi = first_order_islands(find(first_order_islands>0));
imagesc(L == fol(1))

total_lakes = first_order_lakes+second_order_lakes_i;
tl = total_lakes(find(total_lakes>0));

total_islands = first_order_islands+second_order_islands;
ti = total_islands(find(total_islands>0));
%%
close all

L_lake=cell(length(tl),1);
for i = 1:length(fol)
    L_lake{i} = zeros(size(L));
    F_lake{i} = zeros(size(L));
    L_lake{i}(find(L==fol(i))) = L(find(L == fol(i)));
    F_lake{i}(find(L==fol(i))) = ones(size(find(L == fol(i))));
    foi_temp = find(first_order_islands==fol(i));
    if ~isempty(foi_temp) % are there any islands in this lake?
        for ii = 1:length(foi_temp) % if so, find them
            L_lake{i}(find(L==foi_temp(ii))) = L(find(L==foi_temp(ii)));
            F_lake{i}(find(L==foi_temp(ii))) = zeros(size(find(L==foi_temp(ii))));
            %             imagesc(F_lake{i})
            sol_temp = find(second_order_lakes==foi_temp(ii));
            if ~isempty(sol_temp) % are there any lakes in these islands?
                for iii = 1:length(sol_temp)
                    L_lake{i}(find(L==sol_temp(iii))) = L(find(L==sol_temp(iii)));
                    F_lake{i}(find(L==sol_temp(iii))) = ones(size(find(L==sol_temp(iii))));
                    %             imagesc(F_lake{i})
                end
            end
            
        end
    end

    
%     figure()
%     imagesc(F_lake{i})
end

% for i = 1:10%length(fol)
%     F_lake{i} = zeros(size(L));
%     F_lake{i}(find(L==first_order_lakes(tl(i)))) = L(find(L == first_order_lakes(tl(i)))); 
%     foi_temp = find(first_order_islands==tl(i));
%     if ~isempty(foi_temp)
%     for ii = 1:length(foi_temp)
%         F_lake{i}(find(L==foi_temp(ii))) = L(find(L==foi_temp(ii)));
%     end
%     end
%     sol_temp = find(
%     figure()
%     imagesc(F_lake{i})
% end




[B_test, L_test, ~, A_test]  = bwboundaries(F_lake{1},8);
figure()
imagesc(L_test)


% [B,L,n,A] = bwboundaries(tile);
enclosing_boundary  = find(A(6,:));
enclosed_boundaries = find(A(:,6));
% enclosing_boundary  = find(A(2,:));
% enclosed_boundaries = find(A(:,2));
% enclosing_boundary  = find(A(3,:));
% enclosed_boundaries = find(A(:,3));
% enclosing_boundary  = find(A(4,:));
% enclosed_boundaries = find(A(:,4));
% enclosing_boundary  = find(A(5,:));
% enclosed_boundaries = find(A(:,5));
