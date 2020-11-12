
load('test_histcounts.mat')

figure()
% scatter(log(fetch),eq14')
[B,~,idx] = loghistcounts(fetch);
idx = idx+1; % this is because it starts at 0
meaneq14 = accumarray(idx(:),eq14,[],@mean);
meaneq14(meaneq14==0)=NaN;
meanfetch = accumarray(idx(:),fetch,[],@mean);
meanfetch(meanfetch==0)=NaN;
medianeq14 = accumarray(idx(:),eq14,[],@median);
medianeq14(medianeq14==0)=NaN;
medianfetch = accumarray(idx(:),fetch,[],@median);
medianfetch(medianfetch==0)=NaN;
stdeq14 = accumarray(idx(:),eq14,[],@std);
stdeq14(stdeq14==0)=NaN;
stdfetch = accumarray(idx(:),fetch,[],@std);
stdfetch(stdfetch==0)=NaN;

semilogx(meanfetch,meaneq14,'k','LineWidth',2)
hold on
scatter(medianfetch,medianeq14,'k*')
errorbar(meanfetch,meaneq14,stdeq14)
legend('mean','median')
xlabel('Wave weighting')
ylabel('Wavelet variance')
set(gca,'FontSize',14)