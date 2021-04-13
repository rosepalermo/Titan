% automate plots
load('idx_list_mrfn94'); idx_list(33) = NaN; idx_list(46) = NaN; idx_list(isnan(idx_list)) = [];
% load('idx_list_mrfn92'); idx_list(1) = NaN; idx_list(5:6) = NaN; idx_list(21) = NaN; idx_list(isnan(idx_list)) = [];
figure()
hold on
wave_edges = linspace(5e-4,1e0,100);
for runs = 1:length(idx_list)
    subplot(3,1,1)

file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_init_results.mat'];
load(file_name)
[B,~,idx] = loghistcounts(wave_weight);

idx = idx+1; % this is because it starts at 0
meaneq_14 = accumarray(idx(:),eq_14,[],@mean);
meaneq_14(meaneq_14==0)=NaN;
meanwave_weight = accumarray(idx(:),wave_weight,[],@mean);
meanwave_weight(meanwave_weight==0)=NaN;
medianeq_14 = accumarray(idx(:),eq_14,[],@median);
medianeq_14(medianeq_14==0)=NaN;
medianwave_weight = accumarray(idx(:),wave_weight,[],@median);
medianwave_weight(medianwave_weight==0)=NaN;
stdeq_14 = accumarray(idx(:),eq_14,[],@std);
stdeq_14(stdeq_14==0)=NaN;
B = [1 B]';
SEM = stdeq_14./sqrt(B);                         % Standard Error Of The Mean
CI95 = SEM .* tinv(0.975, B-1);              % 95% Confidence Intervals
stdwave_weight = accumarray(idx(:),wave_weight,[],@std);
stdwave_weight(stdwave_weight==0)=NaN;

pp = polyfit(log(meanwave_weight(B>1)),meaneq_14(B>1),1);
slope_init(runs) = pp(1);

semilogx(meanwave_weight(B>1),meaneq_14(B>1),'k','LineWidth',2);
hold on
% plot(wave_weight,eq_14,'.','Color',[0.8 0.8 0.8]);

errorbar(meanwave_weight(B>1),meaneq_14(B>1),CI95(B>1),'k');
% legend('mean','median')
% ylim([0 8e-4])
% xlim([1e-4 1e0])
xlabel('weighted weighted fetch area')
ylabel('azimuthal variance (radians^2)')
% legend('mean','data','95% CI')
set(gca,'FontSize',16)
title('Initial Conditions')

subplot(3,1,2)
file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_uniform_K_0_1_end_results.mat'];
load(file_name)
[B,~,idx] = loghistcounts(wave_weight);

idx = idx+1; % this is because it starts at 0
meaneq_14 = accumarray(idx(:),eq_14,[],@mean);
meaneq_14(meaneq_14==0)=NaN;
meanwave_weight = accumarray(idx(:),wave_weight,[],@mean);
meanwave_weight(meanwave_weight==0)=NaN;
medianeq_14 = accumarray(idx(:),eq_14,[],@median);
medianeq_14(medianeq_14==0)=NaN;
medianwave_weight = accumarray(idx(:),wave_weight,[],@median);
medianwave_weight(medianwave_weight==0)=NaN;
stdeq_14 = accumarray(idx(:),eq_14,[],@std);
stdeq_14(stdeq_14==0)=NaN;
B = [1 B]';
SEM = stdeq_14./sqrt(B);                         % Standard Error Of The Mean
CI95 = SEM .* tinv(0.975, B-1);              % 95% Confidence Intervals
stdwave_weight = accumarray(idx(:),wave_weight,[],@std);
stdwave_weight(stdwave_weight==0)=NaN;

pp = polyfit(log(meanwave_weight(B>1)),meaneq_14(B>1),1);
slope_uni(runs) = pp(1);

semilogx(meanwave_weight(B>1),meaneq_14(B>1),'k','LineWidth',2);
hold on
% plot(wave_weight,eq_14,'.','Color',[0.8 0.8 0.8]);

errorbar(meanwave_weight(B>1),meaneq_14(B>1),CI95(B>1),'k');
% legend('mean','median')
ylim([0 8e-4])
xlim([1e-4 1e0])
xlabel('weighted weighted fetch area')
ylabel('azimuthal variance (radians^2)')
% legend('mean','data','95% CI')
set(gca,'FontSize',16)
title('Uniform')


subplot(3,1,3)
file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_wave_K_0_1_end_results.mat'];
load(file_name)
[B,~,idx] = loghistcounts(wave_weight);

idx = idx+1; % this is because it starts at 0
meaneq_14 = accumarray(idx(:),eq_14,[],@mean);
meaneq_14(meaneq_14==0)=NaN;
meanwave_weight = accumarray(idx(:),wave_weight,[],@mean);
meanwave_weight(meanwave_weight==0)=NaN;
medianeq_14 = accumarray(idx(:),eq_14,[],@median);
medianeq_14(medianeq_14==0)=NaN;
medianwave_weight = accumarray(idx(:),wave_weight,[],@median);
medianwave_weight(medianwave_weight==0)=NaN;
stdeq_14 = accumarray(idx(:),eq_14,[],@std);
stdeq_14(stdeq_14==0)=NaN;
B = [1 B]';
SEM = stdeq_14./sqrt(B);                         % Standard Error Of The Mean
CI95 = SEM .* tinv(0.975, B-1);              % 95% Confidence Intervals
stdwave_weight = accumarray(idx(:),wave_weight,[],@std);
stdwave_weight(stdwave_weight==0)=NaN;

pp = polyfit(log(meanwave_weight(B>1)),meaneq_14(B>1),1);
slope_wave(runs) = pp(1);

semilogx(meanwave_weight(B>1),meaneq_14(B>1),'k','LineWidth',2);
hold on
% plot(wave_weight,eq_14,'.','Color',[0.8 0.8 0.8]);

errorbar(meanwave_weight(B>1),meaneq_14(B>1),CI95(B>1),'k');
% legend('mean','median')
ylim([0 8e-4])
xlim([1e-4 1e0])
xlabel('weighted weighted fetch area')
ylabel('azimuthal variance (radians^2)')
% legend('mean','data','95% CI')
set(gca,'FontSize',16)
title('Wave')
whos x
end
%%
figure()
slope_all = [slope_init slope_uni slope_wave];
slope_edges = linspace(min(slope_all),max(slope_all),50);
histogram(slope_init,slope_edges); 
hold on; 
histogram(slope_uni,slope_edges); 
histogram(slope_wave,slope_edges); 
legend('init','uniform','wave')

dslopeu = slope_uni-slope_init;
dslopew = slope_wave-slope_init;
dslope_edges = linspace(min([dslopeu dslopew]),max([dslopeu dslopew]),50);

figure()
histogram(dslopeu,dslope_edges)
hold on
histogram(dslopew,dslope_edges)
xlabel('change in slope')
legend('uniform - init','wave - init')

%%
% edges_eq14 = linspace(1e-5,9e-4)
% for runs = 1:length(idx_list)
%     figure()
% file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_init_results.mat'];
% load(file_name)
% histogram(eq_14)
% hold on
% file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_uniform_K_0_1_end_results.mat'];
% load(file_name)
% histogram(eq_14)
% file_name = ['/Users/rosepalermo/Documents/Research/Titan/ModelOutput/paper1/022021/r_idx_',num2str(idx_list(runs)),'_wave_K_0_1_end_results.mat'];
% load(file_name) 
% histogram(eq_14)
% 
% pause
% end
