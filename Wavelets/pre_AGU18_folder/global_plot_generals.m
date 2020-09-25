load('global_generals.mat')

figure()
semilogx(global_Save{2},log2(period{2}),'k') %red noise
hold on
semilogx(global_Save{3},log2(period{3}),'r') % wave
semilogx(global_Save{4},log2(period{4}),'b') % rivers
semilogx(global_Save{5},log2(period{5}),'g') % uniform

% Yticks = 2.^(fix(log2(min(period{5}))):fix(log2(max(period{5}))));

xlabel('Power (amplitude^2)')
ylabel('log2(Period)')
set(gca,'YDir','reverse')
% set(gca,'XLim',[0,1.25*max(global_Save{5})])
set(gca,'FontSize',12)
legend('Red Noise','Wave Modified','River Incised','Uniform Erosion')

figure()
semilogy(log2(period{2}),global_Save{2},'k') %red noise
hold on
% semilogy(log2(period{3}),global_Save{3},'r') % wave
% semilogy(log2(period{4}),global_Save{4},'b') % rivers
% semilogy(log2(period{5}),global_Save{5},'g') % uniform

% Yticks = 2.^(fix(log2(min(period{5}))):fix(log2(max(period{5}))));

ylabel('Power (radians^2)')
xlabel('Wavelength (log2(meters))')
% set(gca,'XLim',[0,1.25*max(global_Save{5})])
set(gca,'FontSize',12)
% legend('Red Noise','Wave Modified','River Incised','Uniform Erosion','location','northwest')

figure()
semilogy(log2(period{2}),global_Save{2},'k') %red noise
hold on
semilogy(log2(period{3}),global_Save{3},'r') % wave
semilogy(log2(period{4}),global_Save{4},'b') % rivers
semilogy(log2(period{5}),global_Save{5},'g') % uniform
ylabel('Power (radians^2)')
xlabel('Wavelength (log2(meters))')
set(gca,'XLim',[2 7])
set(gca,'FontSize',12)
legend('Red Noise','Wave Modified','River Incised','Uniform Erosion','location','northwest')