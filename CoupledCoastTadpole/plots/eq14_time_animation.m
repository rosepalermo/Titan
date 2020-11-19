% animation of eq14 plotted on shoreline through time
% uniform
file_name = ('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/111120/river_uniform_Kc0_1.mat');
load(file_name)
v = VideoWriter(file_name(1:end-4));
open(v);
ind = 1:14;
for i = 1:length(ind)
    for kk = 1:10
%     [x,y,wave_weight] = calc_xyw(file_name,7,true);
drawnow
%     set(gca,'Ydir','reverse')

[x,y] = plotisl(p,g,i);
wavelets(x,y,[],8,16,'trash',0)
xlim([0 p.Nx])
ylim([0 p.Ny])
    time_n = strrep(num2str(g.t(ind(i))),'.','_');
    title(string(strcat('time',' ',string(time_n))))
    frame = getframe(gcf);
    writeVideo(v,frame);
    end
end
close(v)

% for i = 1:14
%     drawnow
%     set(gca,'Ydir','reverse')
% xlim([0 p.Nx])
% ylim([0 p.Ny])
% [x,y] = plotisl(p,g,i);
% plot(x,y)
% 
% hold on
% end