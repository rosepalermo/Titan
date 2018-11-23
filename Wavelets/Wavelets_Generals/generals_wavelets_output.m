% wavelet characteristics
model = { 'Initial Conditions'	'Fluvial'	'Uniform'	'Wave'	'Lake Powell'	'Isles of North and South Uist'	'Isles west coast'	'Isles east coast'	'Ligeia Mare'};

% variance
v = [116.9	12.2	111.9	88.6	66.8	78.1	82	47.9	84.4];
% mean
me = [20	20	20	20	20	20	13.8	23	20];
% median
md = [17.6	19.6	17.3	18.6	18.8	20	11.8	22.9	19.5];
% skewness
sk = [1.2	0.4	1.9	1.4	0.9	0.1	1.1	0.3	0.3];

figure()
% subplot(1,2,1)
scatter(md(1),sk(1),1000,'k','.')
hold on
scatter(md(2),sk(2),1000,'b','.')
scatter(md(3),sk(3),1000,'g','.')
scatter(md(4),sk(4),1000,'r','.')
xlabel('Median');ylabel('Skewness')
set(gca,'FontSize',14)
legend(model(1:4),'location','southwest')
figure()
subplot(1,2,2)
scatter(md(1),sk(1),'k')
hold on
scatter(md(2),sk(2),'b')
scatter(md(3),sk(3),'g')
scatter(md(4),sk(4),'r')
scatter(md(5),sk(5),'x','MarkerEdgeColor',[0.7 0.7 0.7]) % Lake Powell
scatter(md(6),sk(6),'d','MarkerEdgeColor',[0.7 0.7 0.7]) % Isles
scatter(md(7),sk(7),'.','MarkerEdgeColor',[0.7 0.7 0.7]) % Isles west
scatter(md(8),sk(8),'s','MarkerEdgeColor',[0.7 0.7 0.7]) % Isles east
scatter(md(9),sk(9),'*','MarkerEdgeColor',[0.7 0.7 0.7]) % LGM
xlabel('Median');ylabel('Variance')
set(gca,'FontSize',14)
legend(model,'location','southwest')
