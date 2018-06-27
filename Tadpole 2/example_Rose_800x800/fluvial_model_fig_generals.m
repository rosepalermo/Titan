load('example_Rose_800x800.mat')
tile = repmat(output,2);
load('xycontours.mat')

figure()
subplot(2,2,1)
imagesc(tile(:,:,1))
axis equal
set(gca,'Xlim',[400 1200]); set(gca,'Ylim',[150 950]);set(gca,'clim',[0 200]);
set(gca,'Ydir','normal')

subplot(2,2,2)
imagesc(tile(:,:,2))
axis equal
set(gca,'Xlim',[400 1200]); set(gca,'Ylim',[150 950]);set(gca,'clim',[0 200]);
set(gca,'Ydir','normal')

subplot(2,2,3)
imagesc(tile(:,:,3))
axis equal
set(gca,'Xlim',[400 1200]); set(gca,'Ylim',[150 950]);set(gca,'clim',[0 200]);
set(gca,'Ydir','normal')

subplot(2,2,4)
plot(x_1m_t1/1000,y_1m_t1/1000,'k')
hold on
plot(x_1m_t2/1000,y_1m_t2/1000,'Color',[0.3 0.3 0.3])
plot(x_1m_t3/1000,y_1m_t3/1000,'Color',[0.7 0.7 0.7])
axis equal tight


figure()
subaxis(1,2,1,'Spacing',0.05,'Margin',0.05)
plot(x_1m_t1/1000,y_1m_t1/1000,'k')
hold on
plot(x_1m_t2/1000,y_1m_t2/1000,'Color',[0.3 0.3 0.3])
plot(x_1m_t3/1000,y_1m_t3/1000,'Color',[0.7 0.7 0.7])
axis equal tight
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])

subaxis(1,2,2,'Spacing',0.05,'Margin',0.05)
plot(x_1m_t1/1000,y_1m_t1/1000,'k')
hold on
plot(x_1m_t2/1000,y_1m_t2/1000,'Color',[0.3 0.3 0.3])
plot(x_1m_t3/1000,y_1m_t3/1000,'Color',[0.7 0.7 0.7])
axis equal tight
set(gca,'XLim',([0.45 0.75])); set(gca,'YLim',([0.6 0.9]))
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])