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
plot(x_1m_t1,y_1m_t1,'k')
hold on
plot(x_1m_t2,y_1m_t2,'Color',[0.3 0.3 0.3])
plot(x_1m_t3,y_1m_t3,'Color',[0.7 0.7 0.7])
axis equal tight