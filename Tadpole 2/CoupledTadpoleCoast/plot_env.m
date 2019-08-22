function plot_env(env)
hold on
for i=1:length(env)
    plot(env{i}(:,1), env{i}(:,2),'k','LineWidth',2)
end