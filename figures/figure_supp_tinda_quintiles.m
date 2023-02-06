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
  [FO_p,FO_pvals_p,t_intervals_p, FO_stat_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));
  hmm_1stlevel.control.quintile.assym_permtest{ip} = FO_stat_p;
  hmm_1stlevel.control.quintile.assym_permtest{ip}.pvals = FO_pvals_p;
  a=[];
  for i=1:K
    for j=1:K
      [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(FO_p(i,j,1,:)), squeeze(FO_p(i,j,2,:)));
    end
  end
  hmm_1stlevel.control.quintile.assym_ttest{ip} = a;
  t_int(ip,:) = round(mean(mean(cellfun(@mean, t_intervals_p),2)));
  mean_direction_p = squeeze(mean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
  FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
  rotational_momentum_p = compute_rotational_momentum(angleplot, FO_assym_p); squeeze(imag(nansum(nansum(angleplot.*FO_assym_p))));
  tmp(:,ip) = rotational_momentum_p;
  
  ax(ip,1) = axes('Position', [-.05, 1.02-0.2*ip, 0.45, 0.13]);
  cyclicalstatesubplot(bestseq,mean_direction_p,a.pvals<hmm_1stlevel.assym_ttest.alpha_thresh);
  title({sprintf('Rotational momentum = %0.3f', mean(rotational_momentum_p, 'omitnan')/hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum),''})

  
  ax(ip,2) = axes('Position', [0.45, 1.03-0.196*ip, 0.5, 0.15]);
  ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
  tmp2(ip) = median(mean(ITmerged,2));
  distributionPlot(sqrt(ITmerged),'showMM',2,'color',{color_scheme{1:size(FO,2)}});
  set(ax(ip,2),'YLim',[0 1.1*max(sqrt(ITmerged(:)))])
  if ip==5
    xlabel('RSN-State')
  end
  grid on; ylabel('Interval Times (s)')
  for ii=1:length(ax(ip,2).YTickLabel)
    ax(ip,2).YTickLabel{ii} = num2str(str2num(ax(ip,2).YTickLabel{ii})^2);
  end
  hmm_1stlevel.control.quintile.FO_intervals(:,:,:,ip) = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO_p,3));
end
hmm_1stlevel.control.quintile.rotational_momentum = tmp;
set_font(10, {'title', 'label'})

save_figure([config.figdir, 'figure_supp_tinda_quintiles/','2supp_Cyclicalpatterns_percentiled']);
close

fig = setup_figure([],1.5,0.5);
boxplot_with_scatter(tmp./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum)
xlabel('Mean interval time (ms)'), ylabel('M')
title('Cycles instantiate over longer temporal intervals')
xticklabels(num2str(t_int));
box off;
save_figure([config.figdir, 'figure_supp_tinda_quintiles/','2supp_rotational_momentum_percentiled']);

% Make an extra plot with median interval time vs rotational momentum?
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));

percentiles = 5:5:95;
clear rotational_momentum_p
for ip=1:length(percentiles)
  [FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,[percentiles(ip)-4 percentiles(ip)+5]);
  FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
  ITmean(ip) = mean(mean(cellfun(@mean,t_intervals_p)));ITmean(ip) = ITmean(ip)./config.sample_rate;
  % ITmedian(ip) = cellfun(@median,t_intervals_p);ITmedian(ip) = ITmedian(ip)./config.sample_rate;
  rotational_momentum_p(ip,:) = compute_rotational_momentum(angleplot, FO_assym_p);
  hmm_1stlevel.control.ntile.FO_intervals(:,:,:,ip) = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO_p,3));
end
hmm_1stlevel.control.ntile.rotational_momentum = rotational_momentum_p;

fig = setup_figure([],2,0.6);
subplot(1,2,1),
shadedErrorBar(percentiles,nanmean(rotational_momentum_p,2), nanstd(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
xlabel('Percentile'), ylabel('Rotational Momentum')
yticks([0])
subplot(1,2,2)
shadedErrorBar((ITmean),nanmean(rotational_momentum_p,2), nanstd(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
yticks([0])
xlabel('Mean IT (s)'), ylabel('Rotational momentum')
save_figure([config.figdir, 'figure_supp_tinda_quintiles/', '2supp_Cyclicalpatterns_rotationalmomentum_percentiled']);
close