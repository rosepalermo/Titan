function plot_rays(start,final)
    num_points = size(final,1);
    px = nan(3*num_points,1);
    py = nan(3*num_points,1);
    px(1:3:end) = start(1);
    py(1:3:end) = start(2);
    px(2:3:end) = final(:,1);
    py(2:3:end) = final(:,2);
    plot(px,py)
end