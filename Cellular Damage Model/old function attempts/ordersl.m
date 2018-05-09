function [xsl_ord,ysl_ord,ind_save] = ordersl(xsl,ysl)
% this didn't work for some places where there is a cell big embayment. need to check and manaully fix after. Less than ideal

sl_trash = [xsl ysl];

ind = 1;
ind_save = zeros(length(xsl),1);
ind_save(1) = ind;
sl(1,:) = sl_trash(1,:);
for i = 2:length(xsl)
    dist = sl_twsrash(ind,:) - sl_trash;
    sl_trash(ind,:) = NaN;
    dist(ind,:) = NaN;
    [~,ind] = nanmin(sum(dist.^2,2));
    ind_save(i) = ind;
    sl(i,:) = sl_trash(ind,:);
end
xsl_ord = sl(:,1);
ysl_ord = sl(:,2);