% measuring the sensitivity to resolution and time step

% take the total strength loss and divide it by the number of shoreline
% cells to get a noramlized strength loss

% input shoreline that is not gridded

load('fetch_input8con','fetch_sl_cells'); % x and y are ordered clockwise if i increases downward, first point != last point
x = fetch_sl_cells{1}(:,1);
y = fetch_sl_cells{1}(:,2);

p.Kuniform = 0.001;
p.Kwave = 0.001;
p.con8 = 1;
p.So = 1;
p.dxo = 100;
p.doStreamPower = 0;
p.doLandslides = 0;
p.doUniformErosion = 1;
p.doWaveErosion = 0;

% grid shoreline and calculate strength loss over 100 years with time steps
% varying from 5 to 100 years
sts = 1:1:10;
tend = 1000;
p.dt = 100;
p.dx = 62.5;
totalSloss_frac = zeros(length(dx),length(dt));
for i = 1:length(sts)
        [lake,~,~] = gridlake(flipud(x)',flipud(y)',dx(i),dx(i),400);
%     [init,~] = test_circle(p,dx(i),200);
%     [init,~] = test_square(p,dx(i),200);
%     lake = ~init;
        [~,totalSloss_frac(i,t)] = calcstrengthloss_sts(lake,tend,p,sts,0);
end
%%
figure()
imagesc(dt,dx,totalSloss_frac*100)
colormap gray
colorbar
xlabel('dt')
ylabel('dx')
set(gca,'FontSize',12)

figure()
plot(dx,totalSloss_frac(:,1)*100)
xlabel('dx')
ylabel('total Strength loss/ domain initial strength')