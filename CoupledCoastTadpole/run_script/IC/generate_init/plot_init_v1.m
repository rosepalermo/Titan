% plot AGU IC
load('river_wideria_200x200.mat')
figure()
subplot(1,2,1); imagesc(init); axis equal
subplot(1,2,2); imagesc(init<40); axis equal

load('/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/results1/cluster/generate_init/ic_generated_1.mat')
init_test1 = cell(84,1);
for i =1:84
init_test1{i} = init{idx_save(i)};
end

drawnow
figure()
for i = 1:84
    subplot(1,2,1)
    imagesc(init_test1{i})
    axis equal
    subplot(1,2,2)
    imagesc(init_test1{i}<40)
    axis equal
    i
    pause
end