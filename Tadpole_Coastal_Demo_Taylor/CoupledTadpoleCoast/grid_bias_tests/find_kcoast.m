%load('test_pkcoast_sl_1.mat')
load('test_pkcoast_sl_100.mat')
test = p.dt*p.Kcoast*shoreline(indshoreline).*dam;
figure()
plot(test)
testing = ones(size(shoreline));
testing(indshoreline) = testing(indshoreline)-test;
imagesc(testing)

% decided upon 4e-8, but could have maybe gone up to 5e-8, depending on if
% we care about the bottom or rhte top. I did the bottom