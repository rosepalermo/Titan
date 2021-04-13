% plot fetch vs azimuth

figure()
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
% title('Initial Conditions')