function plot_env(env,Pobs,V)
hold on
for i=1:length(env)
    plot(env{i}(:,1), env{i}(:,2),'k','LineWidth',2)
end

if size(V) >0
patch( V(:,1) , V(:,2) , 0.1*ones( size(V,1) , 1 ) , ...
           'b' , 'linewidth' , 1.5 );
       alpha(.3)
    plot3( V(:,1) , V(:,2) , 0.1*ones( size(V,1) , 1 ) , ...
           'b.' , 'Markersize' , 5 );
end
       scatter(Pobs(1),Pobs(2),40,'r*')
