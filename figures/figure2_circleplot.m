%% Figure 2: plot circular diagram

if whichstudy~=4
  sigpoints = FO_pvals<(0.05/bonf_ncomparisons);
else
  sigpoints = FO_pvals<0.0000001*(0.05/bonf_ncomparisons);
end



for ext = {'rotationalmomentum'}
  %%
  fig = setup_figure([],2,0.5);
  clear ax
  ax(1) = axes('Position', [0.025, 0, .3, 1]);
  cyclicalstateplot(bestseq,mean_direction, sigpoints,color_scheme,false);
  title({'Observed rotational', sprintf('momentum = %0.3f', mean(hmm_1stlevel.rotational_momentum)./hmm_1stlevel.max_theoretical_rotational_momentum), ''})

  ax(2) = axes('Position', [0.375, 0, .3, 1]);
  if whichstudy==4
    cyclicalstateplot(simulation_average.bestsequencemetrics{2},simulation_average.mean_direction, simulation_average.FO_pvals<0.0000001*(0.05/bonf_ncomparisons),color_scheme,false);
  else
    cyclicalstateplot(simulation_average.bestsequencemetrics{2},simulation_average.mean_direction, simulation_average.FO_pvals<(0.05/bonf_ncomparisons),color_scheme,false);
  end
title({'Simulated rotational', sprintf('momentum = %0.3f', mean(simulation_average.rotational_momentum)./hmm_1stlevel.max_theoretical_rotational_momentum), ''})
  
  if strcmp(ext, 'FOasym') || strcmp(ext{1}, 'FOasym')
    tmp = [simulation_average.TIDA, hmm_1stlevel.TIDA];
    ylb = 'Mean FO asymmetry';
    [stat.H, stat.P, stat.CI, stat.stats] = ttest(tmp(:,2),tmp(:,1), 'tail', 'right');
  elseif strcmp(ext, 'rotationalmomentum') || strcmp(ext{1}, 'rotationalmomentum')
    tmp = [simulation_average.rotational_momentum, hmm_1stlevel.rotational_momentum]./hmm_1stlevel.max_theoretical_rotational_momentum;
    ylb = 'Rotational momentum';
    [stat.H, stat.P, stat.CI, stat.stats] = ttest(tmp(:,2),tmp(:,1));
  end
  ax(3) = axes('Position', [0.755, 0.15, 0.25, 0.7]); hold on
  vline(0, '--k')
  
  clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
  for k=1:2
    scatter(ones(size(tmp(:,k))).*(k+(rand(size(tmp(:,k)))-0.5)/2),tmp(:,k), 'filled', 'MarkerFaceColor', clr{k}, 'MarkerFaceAlpha',0.6)
  end
  h=boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this back on top
  set(h, 'linew', 2)
  box off
  if stat.H
    sigstar({[1,2]}, stat.P)
  end
%   title({'Difference between', 'observed data and', 'simulations', ''})
  set(ax(3), 'XTickLabels', {'simulated', 'observed'})
  ylabel(ylb)
  set_font(10, {'title', 'label'})
  fname = [config.figdir,'2_cyclicalpattern_', ext{1}];
  save_figure(fname);
  save([fname, 'stat'], 'stat')
end


% also plot the *real* circle plot individually
cyclicalstateplot(bestseq,mean_direction, sigpoints,color_scheme);
save_figure([config.figdir,'2ind_Cyclicalpattern']);



% as well as the simulated
cyclicalstateplot(simulation_average.bestsequencemetrics{2},simulation_average.mean_direction, simulation_average.FO_pvals<(0.05/bonf_ncomparisons), color_scheme);
save_figure([config.figdir,'2ind_Cyclicalpattern_simulated']);

