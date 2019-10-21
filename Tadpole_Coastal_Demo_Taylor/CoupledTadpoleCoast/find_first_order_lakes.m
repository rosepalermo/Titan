function [F_lake,total_lakes,total_islands,L] = find_first_order_lakes(lake)

[bound, L, ~, A]  = bwboundaries(lake,8);
n = length(bound);

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
second_order_lakes_i = zeros(n,1);
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
% close all

L_lake=cell(length(tl),1); %first order lake and everything inside with numbers from L
F_lake=cell(length(tl),1); %first order lake and everything inside with ones and zeros (1 = lake, 0 = land)
for i = 1:length(tl)
    L_lake{i} = zeros(size(L));
    F_lake{i} = zeros(size(L));
    L_lake{i}(find(L==tl(i))) = L(find(L == tl(i)));
    F_lake{i}(find(L==tl(i))) = ones(size(find(L == tl(i))));
    foi_temp = find(first_order_islands==tl(i));
    if ~isempty(foi_temp) % are there any first order islands in this lake?
        for ii = 1:length(foi_temp) % if so, find them
            L_lake{i}(find(L==foi_temp(ii))) = L(find(L==foi_temp(ii)));
            F_lake{i}(find(L==foi_temp(ii))) = zeros(size(find(L==foi_temp(ii))));
        end
    end
    soi_temp = find(second_order_islands==tl(i));
    if ~isempty(soi_temp) % are there any second order islands in this lake?
        for ii = 1:length(soi_temp) % if so, find them
            L_lake{i}(find(L==soi_temp(ii))) = L(find(L==soi_temp(ii)));
            F_lake{i}(find(L==soi_temp(ii))) = zeros(size(find(L==soi_temp(ii))));
        end
    end

end



end
