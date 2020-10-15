% test different meander_titan scripts to find differences
load('meander_test.mat')
[azi1,curvi1,deltad1] = meander_titan(x,y);
[azi2,curvi2,deltad2] = meander_titan_2(x,y);
[azi3,curvi3,deltad3] = meander_titan_3(x,y);
azi1 = azi1(1:end-1);
azi2 = azi2(1:end-1);
azi3 = azi3(1:end-1);

figure()
plot(azi1)
hold on
plot(azi2)
plot(azi3)
legend('1','2','3')