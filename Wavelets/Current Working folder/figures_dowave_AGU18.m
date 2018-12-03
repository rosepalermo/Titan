% dowave figures

figure() 
plot(t/100,eq14)
xlabel('distance along coast (km)')
title('sum of eq14')

figure()
h = histogram(eq14,10,'Normalization','probability')
% h = findobj(gca,'Type','patch');
h.FaceColor = 'k';
% h.FaceColor = [0.6 0.6 0.6];
h.EdgeColor = 'w';
title('sum of eq14')

if i == 12
    figure()
    scatter3(t,y,eq14,[],eq14,'.')
    title('sum of eq14')
    view(2)
    colorbar
    set(gca,'Clim',[0 0.00025])
    axis equal tight
elseif i == 2|i == 3|i == 4|i == 5
    figure()
    scatter3(xx/1e3,yy/1e3,eq14,[],eq14,'.')
    view(2)
    axis equal tight
    set(gca,'XLim',([0.5 0.8])); set(gca,'YLim',([0.5 0.8])); %set(gca,'Clim',[0 mean(rness)+2*std(rness)])
    set(gca,'Clim',[0 0.00025])
    set(gca,'FontSize',14)
    set(gca,'xtick',[],'ytick',[])
    set(gca,'xticklabel',[],'yticklabel',[])
    text(xx(1:500:end-500)/1e3,yy(1:500:end-500)/1e3,plotidx,'FontSize',14);
    
    figure()
    scatter3(xx,yy,eq14,[],eq14,'.')
    title('sum of eq14')
    view(2)
    colorbar
    set(gca,'Clim',[0 0.00025])
    axis equal tight
    text(xx(1:500:end-500),yy(1:500:end-500),plotidx,'FontSize',14);
    
    else
    figure()
    scatter3(xx,yy,eq14,[],eq14,'.')
    title('sum of eq14')
    view(2)
    colorbar
    set(gca,'Clim',[0 0.00025])
    axis equal tight
    text(xx(1:500:end-500),yy(1:500:end-500),plotidx,'FontSize',14);
end

figure()
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
imagesc(t,log2(period(pband1)),log2(power(pband1,:)));  %*** uncomment for 'image' plot
colormap jet
set(gca,'FontSize',12)
xlabel('alongshore distance (m)')
ylabel('Period')
title('Wavelet Power Spectrum')
set(gca,'XLim',xlim(:))
set(gca,'CLim',[-12 8])



