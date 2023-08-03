%% Figure Supplement 2: analyse by intervals >2 heartbeats long

% assume HR are greater than 50BPM:
maxHR = (60/50*2); % 2.4 Hz is twice the lowest HB
percentile = [maxHR*config.sample_rate,NaN]; % shortest permissible intervals:

[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentile);
mean_direction_p = squeeze(nanmean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
rotational_momentum_p = compute_rotational_momentum(angleplot, FO_assym_p);
ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
a=[];
for i=1:K
  for j=1:K
    [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(FO_p(i,j,1,:)), squeeze(FO_p(i,j,2,:)));
  end
end


fig = setup_figure([], 2,0.4); clear ax
ax(1) = axes('Position', [-.03, 0.1, 0.45, 0.7]);
cyclicalstatesubplot(bestseq,mean_direction_p,a.pvals<hmm_1stlevel.assym_ttest.alpha_thresh);
title({sprintf('S = %0.3f', mean(rotational_momentum_p, 'omitnan')/hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum),''})

ax(2) = axes('Position', [0.5, 0.15, 0.45, 0.7]);
distributionPlot(ITmerged,'showMM',2,'color',{color_scheme{1:size(FO,2)}});

set(gca,'YLim',[0 1.1*max(ITmerged(:))])
title('Interval times > 2 sec');xlabel('State'), ylabel('Interval time (sec)');grid on;
set_font(10, {'title', 'label'})

save_figure([config.figdir, 'figure_supp_tinda_heartbeat/','2supp_Cyclicalpatterns_greaterthan2sec'])
close

% save metrics:
hmm_1stlevel.control.cardiac.FO_intervals = FO_p;
hmm_1stlevel.control.cardiac.rotational_momentum = rotational_momentum_p;
hmm_1stlevel.control.cardiac.FO_assym = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO_p,3));
