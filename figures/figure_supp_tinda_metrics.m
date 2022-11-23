%% make a boxplot seperately for each state (vs simulation)
clear stat
fig = setup_figure([],2,.75);
for whichstate=1:K+1
  if whichstate==K+1
    fig = setup_figure([],1,.6); hold on;
    tmp = [simulation_average.TIDA, hmm_1stlevel.cycle_metrics.TIDA];
  else
    subplot(4,3,whichstate), hold on
    tmp = [simulation_average.TIDA_perstate(:, whichstate), hmm_1stlevel.cycle_metrics.TIDA_perstate(:,whichstate)];
  end
  [stat{whichstate}.H, stat{whichstate}.P, stat{whichstate}.CI, stat{whichstate}.stats] = ttest(tmp(:,2),tmp(:,1), 'tail', 'right');
  clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
  for k=1:2
    scatter(ones(size(tmp(:,k))).*(k+(rand(size(tmp(:,k)))-0.5)/2),tmp(:,k),'filled', 'MarkerFaceAlpha',0.6)
  end
  boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000) % put this on top
  if stat{whichstate}.H
    sigstar({[1,2]}, stat{whichstate}.P)
  end
  box off
  % ylim([-.36, 0.36])
  set(gca, 'XTickLabels', {'simulated', 'observed'})
  ylabel({'Mean FO asym'})
  if whichstate==K
    fname=[config.figdir, 'figure_supp_tinda_metrics/','2supp_TIDA_perstate'];
    save_figure(fname);
  elseif whichstate==K+1
    fname=[config.figdir, 'figure_supp_tinda_metrics/','2supp_TIDA'];
    save_figure(fname);
  end
end

%% TODO: make figures for all perstate and per subject measures. Also for simulations