%% Figure 2 supplement: analyse intervals <half a respiratory cycle

% assume resp cycle of 16BPM;
minRR = 60/16*0.5;
percentile = [NaN,minRR*config.sample_rate]; % shortest permissible intervals:
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentile);
mean_direction_p = squeeze(nanmean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
rotational_momentum_p = squeeze(imag(sum(sum(angleplot.*FO_assym_p))));
ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;

fig = setup_figure([], 2,0.4); clear ax
ax(1) = axes('Position', [-.03, 0.1, 0.45, 0.7]);
cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(0.05/bonf_ncomparisons));
title({'Rotational momentum', sprintf('= %0.3f', mean(rotational_momentum_p, 'omitnan')/hmm_1stlevel.max_theoretical_rotational_momentum),''})

ax(2) = axes('Position', [0.5, 0.15, 0.45, 0.7]);
distributionPlot(ITmerged,'showMM',2,'color',{color_scheme{1:size(FO,2)}});

set(gca,'YLim',[0 1.1*max(ITmerged(:))])
title('Interval times < 2.5 sec');xlabel('RSN-State'), ylabel('Interval time (sec)');grid on;
set_font(10, {'title', 'label'})
save_figure([config.figdir,'2supp_Cyclicalpatterns_lessthan2.5sec.png'])

hmm_1stlevel.FO_respcontrol = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO,3));
