%% Figure 2 supplement:  analyse by quintiles
percentiles = 0:20:100;
fig=setup_figure([],2,2);
if whichstudy==4
  alpha_thresh = 0.0000001;
else
  alpha_thresh = 0.05;
end
clear ax
for ip=1:length(percentiles)-1
  [FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));
  mean_direction_p = squeeze(mean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
  FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
  rotational_momentum_p = squeeze(imag(nansum(nansum(angleplot.*FO_assym_p))));
  
  ax(ip,1) = axes('Position', [-.05, 1.02-0.2*ip, 0.45, 0.13]);
  cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(alpha_thresh/bonf_ncomparisons));
  title({sprintf('Rotational momentum = %0.3f', mean(rotational_momentum_p, 'omitnan')/hmm_1stlevel.max_theoretical_rotational_momentum),''})

  
  ax(ip,2) = axes('Position', [0.45, 1.03-0.196*ip, 0.5, 0.15]);
  ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
  distributionPlot(sqrt(ITmerged),'showMM',2,'color',{color_scheme{1:size(FO,2)}});
  set(ax(ip,2),'YLim',[0 1.1*max(sqrt(ITmerged(:)))])
  if ip==5
    xlabel('RSN-State')
  end
  grid on; ylabel('Interval Times (s)')
  for ii=1:length(ax(ip,2).YTickLabel)
    ax(ip,2).YTickLabel{ii} = num2str(str2num(ax(ip,2).YTickLabel{ii})^2);
  end
end
set_font(10, {'title', 'label'})

save_figure([config.figdir,'2supp_Cyclicalpatterns_percentiled']);
close
hmm_1stlevel.FO_quartile = FO_p;

% Make an extra plot with median interval time vs rotational momentum?
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));

percentiles = 5:5:95;
for ip=1:length(percentiles)
  [FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,[percentiles(ip)-4 percentiles(ip)+5]);
  FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
  ITmean(ip) = mean(mean(cellfun(@mean,t_intervals_p)));ITmean(ip) = ITmean(ip)./config.sample_rate;
  % ITmedian(ip) = cellfun(@median,t_intervals_p);ITmedian(ip) = ITmedian(ip)./config.sample_rate;
  rotational_momentum_p(ip,:) = squeeze(imag(nansum(nansum(angleplot.*FO_assym_p))));
end

fig = setup_figure([],2,0.6);
subplot(1,2,1),
shadedErrorBar(percentiles,mean(rotational_momentum_p,2), std(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
xlabel('Percentile'), ylabel('Rotational Momentum')
yticks([0])
subplot(1,2,2)
shadedErrorBar((ITmean),mean(rotational_momentum_p,2), std(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
yticks([0])
xlabel('Mean IT (s)'), ylabel('Rotational Momentum')
save_figure([config.figdir,'2supp_Cyclicalpatterns_rotationalmomentum_percentiled']);
close