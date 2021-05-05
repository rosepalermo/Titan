%% Plot fetch on shoreline
    hh = figure();
%     subplot('position',[457 113 857 692]);
    subplot(2,2,2)
    scatter3(x,y,eq_14',[],eq_14','filled')
    view(2)
    grid off
    axis tight
    axis equal
%     xlim([0 200])
%     ylim([20 180])
%         set(gca,'Ydir','reverse')
    h = colorbar;
    ylabel(h, 'azimuthal variance')
%     set(gca,'Clim',[0 3e-4])
    xlabel('meters');ylabel('meters');
%     load('clim_eq14.mat')
%     set(gca,'Clim',clim_eq14)
    set(gca,'FontSize',16)
%     set(gca,'ColorScale','log')
%     set(gca,'YTickLabel',[])
%     set(gca,'XTickLabel',[])
%     set(gca,'fontsize',18)
    grid off
%     if save_on
%         fig = '.eps'; fig_suf ='_eq14_sl'; 
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
subplot(2,2,1)
scatter3(x,y,wave_weight,[],wave_weight,'filled')
view(2)
grid off
axis tight
axis equal
%     xlim([0 200])
%     ylim([20 180])
%         set(gca,'Ydir','reverse')
h = colorbar;
ylabel(h, 'fetch area')
%     set(gca,'Clim',[1e9 1e10])
xlabel('meters');ylabel('meters');
%     load('clim_eq_14.mat')
%     set(gca,'Clim',clim_eq_14)
set(gca,'FontSize',16)
%     set(gca,'ColorScale','log')
%     set(gca,'YTickLabel',[])
%     set(gca,'XTickLabel',[])
%     set(gca,'fontsize',18)
%     if save_on
%         fig = '.eps'; fig_suf ='_ww_sl';
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end

%% Plot fetch vs roughness
% figure()
% % scatter(log(wave_weight),eq_14')
% X = [log(wave_weight),eq_14'];
% X(X==inf) = NaN;
% X(X==-inf) = NaN;
% hist3(X,'CdataMode','auto'); view(2)
% xlabel('log fetch')
% ylabel('azimuthal variance (radians^2)')
% set(gca,'FontSize',16)
%     if save_on
%         fig = '.eps'; fig_suf ='fvr';
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end

%% Plot fetch vs roughness
%% linear
% figure()
% % scatter(log(wave_weight),eq_14')
% [B,~,idx] = histcounts(wave_weight);
% % plot(B);
% % xlabel('bin')
% % ylabel('N')
%     if save_on
%         fig = '.eps'; fig_suf ='hist_lin';
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
% figure()
% idx = idx+1; % this is because it starts at 0
% meaneq_14 = accumarray(idx(:),eq_14,[],@mean);
% meaneq_14(meaneq_14==0)=NaN;
% meanfetch = accumarray(idx(:),wave_weight,[],@mean);
% meanfetch(meanfetch==0)=NaN;
% medianeq_14 = accumarray(idx(:),eq_14,[],@median);
% medianeq_14(medianeq_14==0)=NaN;
% medianfetch = accumarray(idx(:),wave_weight,[],@median);
% medianfetch(medianfetch==0)=NaN;
% stdeq_14 = accumarray(idx(:),eq_14,[],@std);
% stdeq_14(stdeq_14==0)=NaN;
% B = [1 B]';
% SEM = stdeq_14./sqrt(B);                         % Standard Error Of The Mean
% CI95 = SEM .* tinv(0.975, B-1);              % 95% Confidence Intervals
% stdfetch = accumarray(idx(:),wave_weight,[],@std);
% stdfetch(stdfetch==0)=NaN;
%
% plot(meanfetch,meaneq_14,'k','LineWidth',2)
% hold on
% plot(wave_weight,eq_14,'.','Color',[0.8 0.8 0.8])
% % scatter(medianfetch,medianeq_14,'k*')
% errorbar(meanfetch,meaneq_14,CI95,'k')
% % legend('mean','median')
% xlabel('weighted fetch area')
% ylabel('azimuthal variance (radians^2)')
% set(gca,'FontSize',16)
%     if save_on
%         fig = '.eps'; fig_suf ='fvr_mean';
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
%% log
%     figure()
subplot(2,2,3)
scatter(log(wave_weight),eq_14')
[B,~,idx] = loghistcounts(wave_weight);
% plot(B);
% xlabel('bin')
% ylabel('N')
%     if save_on
%         fig = '.eps'; fig_suf ='hist_log';
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
% figure()
idx = idx+1; % this is because it starts at 0
meaneq_14 = accumarray(idx(:),eq_14,[],@mean);
meaneq_14(meaneq_14==0)=NaN;
meanfetch = accumarray(idx(:),wave_weight,[],@mean);
meanfetch(meanfetch==0)=NaN;
medianeq_14 = accumarray(idx(:),eq_14,[],@median);
medianeq_14(medianeq_14==0)=NaN;
medianfetch = accumarray(idx(:),wave_weight,[],@median);
medianfetch(medianfetch==0)=NaN;
stdeq_14 = accumarray(idx(:),eq_14,[],@std);
stdeq_14(stdeq_14==0)=NaN;
B = [1 B]';
SEM = stdeq_14./sqrt(B);                         % Standard Error Of The Mean
CI95 = SEM .* tinv(0.975, B-1);              % 95% Confidence Intervals
stdfetch = accumarray(idx(:),wave_weight,[],@std);
stdfetch(stdfetch==0)=NaN;

p_2 = semilogx(meanfetch(B>1),meaneq_14(B>1),'k','LineWidth',2);
hold on
p2_2 = plot(wave_weight,eq_14,'.','Color',[0.8 0.8 0.8]);
% plot(wave_weight,eq_14,'.','Color','g')
% scatter(medianfetch,medianeq_14,'k*')
p3_2 = errorbar(meanfetch(B>1),meaneq_14(B>1),CI95(B>1),'k');
% legend('mean','median')
% ylim([0 5e-4])
xlabel('weighted fetch area')
ylabel('azimuthal variance (radians^2)')
% legend('mean','data','95% CI')
set(gca,'FontSize',16)
%     if save_on
%         fig = '.eps'; fig_suf ='fvr_mean_log';
%         range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%         figname = strcat(savename,range,fig_suf,'.jpg');
%         saveas(gcf,figname)
%     end
pp = polyfit(log(meanfetch(B>1)),meaneq_14(B>1),1);
slope = pp(1)

% Plot fetch vs roughness
%  figure()
subplot(2,2,4)
scatter(log(wave_weight),eq_14')
% scatter(log(wave_weight),eq_14')
X = [log(wave_weight),eq_14'];
X(X==inf) = NaN;
X(X==-inf) = NaN;
hist3(X,'CdataMode','auto','Nbins',[20 20]); view(2)
xlabel('log fetch')
ylabel('wavelet variance')
set(gca,'FontSize',14)
colorbar
%      if save_on
%          fig = '.eps'; fig_suf ='fvr'; figname = strcat(savename,fig_suf);
%          fig = '.eps'; fig_suf ='fvr';
%          range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
%          figname = strcat(savename,range,fig_suf,'.jpg');
%          saveas(gcf,figname)
%      end
hh.Position = [457 113 857 692];
savename = strcat(file_name(1:end-4),'_');

fig = '.eps'; fig_suf ='mapv_fvr'; figname = strcat(savename,fig_suf);
range = strcat('_min',num2str(pmin),'_max',num2str(pmax));
figname = strcat(savename,range,fig_suf,'.jpg');
saveas(gcf,figname)

