%% Plot HMM summary statistics
% summary plots:
fig = setup_figure([],1,2.25);
subplot(3,1,1)
distributionPlot(hmm_1stlevel.FO,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
set(gca,'YLim',[0 1.1*max(FracOcc(:))]);
xlabel('RSN-state'), ylabel('Proportion'), title('Fractional Occupancy')
grid on;

subplot(3,1,2)
distributionPlot(hmm_1stlevel.LT_mu ./ config.sample_rate * 1000,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Life Times'); xlabel('RSN-state'), ylabel('Time (ms)'),
grid on;
YL = 1.1*max(hmm_1stlevel.LT_mu(:))./ config.sample_rate * 1000;
set(gca,'YLim',[0 YL]);

subplot(3,1,3);
% NOTE: if we want to use a logscale in the axis, enable this
%{
distributionPlot(log10(ITmerged ./ config.sample_rate * 1000),'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Interval Times');xlabel('RSN-state'), ylabel('Time (ms)');
grid on
YL(2) =1.5* max(mean(log10(ITmerged ./ config.sample_rate * 1000)));
YL(1) = min(squash(log10(ITmerged ./ config.sample_rate * 1000)));
set(gca,'YLim',YL)
set(gca,'YTick',(log10(1000*[0.05,0.1,0.5,1,5,10])))
y_labels = get(gca,'YTickLabel');
for i=1:length(y_labels)
  y_labels{i}=num2str(10.^(str2num(y_labels{i})),1);
end
set(gca,'YTickLabels',y_labels);
%}
% print([config.figdir '0_temporalstats_IT_logscale'],'-depsc')

distributionPlot(hmm_1stlevel.IT_mu ./ config.sample_rate*1000,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Interval Times');xlabel('RSN-state'), ylabel('Time (ms)');grid on
YL(2) =1.5* max(mean((hmm_1stlevel.IT_mu ./ config.sample_rate * 1000)));
YL(1) = 0;
set(gca,'YLim',YL)

save_figure(fig, [config.figdir, 'figure_supp_hmm_stats/', '0_HMM_summary_statistics'])