%% Figure 2: plot circular diagram
if ~exist('simulation', 'var'), simulation = hmm_1stlevel.simulation; end
if ~exist('simulation_average', 'var'), simulation_average = hmm_1stlevel.simulation_average; end
alpha_thresh = 0.05/(K.^2-K);

if whichstudy==4
  alpha_thresh = 0.0000001 * alpha_thresh;
end


clr = [{[0 0.4470 0.7410]}, {[0.8500 0.3250 0.0980]}, {[0.9290 0.6940 0.1250]}];
study = {'MEGUK (N=55)', 'MEGUK (N=55)', 'HCP (N=79)', 'Cam-CAN (N=600)'};
sigpoints = hmm_1stlevel.assym_ttest.sigpoints;

%% Figure 2 (per study)
fig = setup_figure([],2,.45);
clear ax

% asymmetry matrix
ax(1) = axes('Position', [0.05, 0.0, 0.225, .775]); hold on
CL = max(abs(hmm_1stlevel.cycle_metrics.mean_assym(:)))*[-1 1];
cmap = flipud(brewermap(256, 'RdBu'));
tmp = hmm_1stlevel.cycle_metrics.mean_assym(bestseq, bestseq);
tmp(isnan(tmp))=0;
imagesc(tmp, CL), colormap(cmap)
set(gca, 'YDir', 'reverse')
xticks(1:K), yticks(1:K), xticklabels(bestseq), yticklabels(bestseq), colorbar('Location', 'SouthOutside')
ylabel('From state')
xlabel('To state')
title('Mean FO asymmetry')
xlim([0.5 12.5])
ylim([0.5 12.5])

% cycle plot
ax(2) = axes('Position', [1/3, 0, .3, .9]); %ax(1) = axes('Position', [0.025, 0, .2, 1]);
cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_direction, hmm_1stlevel.assym_ttest.sigpoints,color_scheme,false);
title({study{whichstudy}, sprintf('M = %0.3f', mean(hmm_1stlevel.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), ''})

% rotational momentum vs perm
if whichstudy==4
  facealpha=0.2;
else
  facealpha=0.6;
end
ax(3) = axes('Position', [0.725, 0.05, 0.25, 0.8])
tmp1=hmm_1stlevel.cycle_metrics.rotational_momentum./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
boxplot_with_scatter([tmp1*nan, tmp1], clr([2,1]), facealpha)
tmp2 = mean(hmm_1stlevel.perm.rotational_momentum,2)'./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
scatter(ones(size(tmp2)).*((2-1)+(rand(size(tmp2))-0.5)/2),tmp2,25, brewermap(length(tmp2), 'Dark2'), 'filled', 'MarkerFaceAlpha',0.6)
h1 = boxplot([tmp2;nan(1,length(tmp2))]', 'Colors', 'k', 'Width', .75, 'whisker', 1000)
% set(h, {'linew'}, {2})
ylim(1.1*[min([min(tmp1), min(tmp2)]), max([max(tmp1), max(tmp2)])])
sigstar({[1,2]}, hmm_1stlevel.perm.obs_vs_perm)
xlim([-1 3.5])
box off
xticks([1 2])
ylabel('M')
xticklabels({'perm', 'obs'})
view([90,-90])
c1=brewermap(length(tmp2), 'Dark2');
for k=1:length(c1)
  cmap1{k} = c1(k,:);
end

if whichstudy==1
  y1 = -.1;
  x2 = 0.785;
elseif whichstudy==3
  y1 = -0.09;
  x2 = 0.78;
elseif whichstudy==4
  y1 = -0.075;
  x2 = 0.778;
end
text(3.2, y1+.02, 'individual subjects')
plot(3.2, y1-.01, '.', 'Color', clr{1}, 'MarkerSize', 10)
text(2.8, y1+.02, 'individual perm')
plot(2.8, y1-.02, '.', 'Color', cmap1{5}, 'MarkerSize', 10)
plot(2.8, y1-.01, '.', 'Color', cmap1{2}, 'MarkerSize', 10)
plot(2.8, y1, '.', 'Color', cmap1{6}, 'MarkerSize', 10)
title('Rotational Momentum')

% inset
ax(4) = axes('Position', [x2, .16, 0.15, 0.19]), hold on
facealpha = 0.6;
scatter(ones(size(tmp2)).*((2-1)+(rand(size(tmp2))-0.5)/2),tmp2,25, c1, 'filled', 'MarkerFaceAlpha',0.6)
h = boxplot([tmp2]', 'Colors', 'k', 'whisker', 1000, 'Width', .75)
set(h, {'linew'}, {2})
xlim([0.5 1.5])
xticklabels({'zoom'})
ax(4).YAxis.Exponent = 0;
view([90,-90])

set_font(10, {'title', 'label'})
fname = [config.figdir, 'figure2_circleplot/','2_cyclicalpattern'];
save_figure(fname);


%% Supplements
for ext = {'rotationalmomentum'}
 
  fig = setup_figure([],2,0.7);
  clear ax
  ax(1) = axes('Position', [0.025, .45, .3, .4]); %ax(1) = axes('Position', [0.025, 0, .2, 1]);
  cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_direction, sigpoints,color_scheme,false);
  title({sprintf('M = %0.3f', mean(hmm_1stlevel.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), ''})
  
  ax(2) = axes('Position', [0.35, .45, .3, .4]); %  ax(2) = axes('Position', [0.275, 0, .2, 1]);
  cyclicalstateplot(simulation{1}.bestsequencemetrics{1},simulation{1}.cycle_metrics.mean_direction, simulation{1}.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh,color_scheme,false);
  title({sprintf('M = %0.3f', mean(simulation{1}.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), ''})
  
  ax(4) = axes('Position', [0.675, .45, .3, .4]);%ax(4) = axes('Position', [0.525, 0, .2, 1]);
  cyclicalstateplot(simulation_average.bestsequencemetrics{1},simulation_average.cycle_metrics.mean_direction, simulation_average.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh,color_scheme,false);
  title({sprintf('M = %0.3f', mean(simulation_average.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), ''})
  
  if strcmp(ext, 'FOasym') || strcmp(ext{1}, 'FOasym')
    tmp = [hmm_1stlevel.cycle_metrics.TIDA, simulation{1}.TIDA, simulation_average.TIDA];
    ylb = 'Mean FO asymmetry';
  elseif strcmp(ext, 'rotationalmomentum') || strcmp(ext{1}, 'rotationalmomentum')
    tmp = [hmm_1stlevel.cycle_metrics.rotational_momentum, simulation{1}.cycle_metrics.rotational_momentum, simulation_average.cycle_metrics.rotational_momentum]./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
    ylb = {'                     clockwise \bf \leftarrow M \rightarrow \rm counterclockwise'};
  end
  %   ax(3) = axes('Position', [0.8, 0.15, 0.15, 0.7]); hold on
  ax(3) = axes('Position', [0.1, 0.05, 0.45, .3]); hold on
  vline(0, '--k')
  
  clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
  for k=1:3
    scatter(ones(size(tmp(:,k))).*(k+(rand(size(tmp(:,k)))-0.5)/2),tmp(:,k), 'filled', 'MarkerFaceColor', clr{k}, 'MarkerFaceAlpha',0.6)
  end
  h=boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this back on top
  set(h, 'linew', 2)
  box off
  if strcmp(ext, 'FOasym') || strcmp(ext{1}, 'FOasym')
    sigstar({[1,2], [1,3]}, [hmm_1stlevel.metric_vs_sim.TIDA.prob, hmm_1stlevel.metric_vs_sim_avg.TIDA.prob])
  elseif strcmp(ext, 'rotationalmomentum') || strcmp(ext{1}, 'rotationalmomentum')
    sigstar({[1,2], [1,3]}, [hmm_1stlevel.metric_vs_sim.rotational_momentum.prob, hmm_1stlevel.metric_vs_sim_avg.rotational_momentum.prob])
  end
  %   title({'Difference between', 'observed data and', 'simulations', ''})
  set(ax(3), 'XTickLabels', {'obs', '1 sim', '100 sim'})
  %   xtickangle(45)
  ylabel(ylb)
  view([-270 90])
  title('Rotational momentum vs. simulations')
  
  
  ax(5) = axes('Position', [0.65, 0.075, 0.3, .265]); hold on
  CL = max(abs(hmm_1stlevel.cycle_metrics.mean_assym(:)))*[-1 1];
  cmap = flipud(brewermap(256, 'RdBu'));
  tmp = hmm_1stlevel.cycle_metrics.mean_assym(bestseq, bestseq);
  tmp(isnan(tmp))=0;
  imagesc(tmp, CL), colormap(cmap)
  set(gca, 'YDir', 'reverse')
  xticks(1:K), yticks(1:K), xticklabels(bestseq), yticklabels(bestseq), colorbar
  ylabel('From state')
  xlabel('To state')
  title('Mean FO asymmetry')
  xlim([0.5 12.5])
  ylim([0.5 12.5])
  
  set_font(10, {'title', 'label'})
  fname = [config.figdir, 'figure2_circleplot/','2_cyclicalpattern_', ext{1}];
  save_figure(fname);
  
end
%%

% also plot the *real* circle plot individually
fig = setup_figure([],1,1);
cyclicalstateplot(bestseq,hmm_1stlevel.cycle_metrics.mean_direction, sigpoints,color_scheme, false);
title({sprintf('M = %0.3f', mean(hmm_1stlevel.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), ''})
save_figure([config.figdir, 'figure2_circleplot/','2ind_Cyclicalpattern']);

% as well as the simulated
fig = setup_figure([],1,1);
cyclicalstateplot(simulation{1}.bestsequencemetrics{1},simulation{1}.cycle_metrics.mean_direction, simulation{1}.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);
title({sprintf('M = %0.3f', mean(simulation{1}.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), ''})
save_figure([config.figdir, 'figure2_circleplot/','2ind_Cyclicalpattern_simulation1']);

% and the simulation average
fig = setup_figure([],1,1);
cyclicalstateplot(simulation_average.bestsequencemetrics{1},simulation_average.cycle_metrics.mean_direction, simulation_average.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);
title({sprintf('M = %0.3f', mean(simulation_average.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), ''})
save_figure([config.figdir, 'figure2_circleplot/','2ind_Cyclicalpattern_simulationAvg']);


% Asym matrix
fig = setup_figure([],1,1);
CL = max(abs(hmm_1stlevel.cycle_metrics.mean_assym(:)))*[-1 1];
cmap = flipud(brewermap(256, 'RdBu'));
tmp = hmm_1stlevel.cycle_metrics.mean_assym(bestseq, bestseq);
tmp(isnan(tmp))=0;
imagesc(tmp, CL), colormap(cmap)
set(gca, 'YDir', 'reverse')
xticks(1:K), yticks(1:K), xticklabels(bestseq), yticklabels(bestseq), colorbar
ylabel('From state')
xlabel('To state')
title('Mean FO asymmetry')
xlim([0.5 12.5])
ylim([0.5 12.5])

set_font(10, {'title', 'label'})
fname = [config.figdir, 'figure2_circleplot/','2supp_FOassym_matrix'];
save_figure(fname);


%% Rotational momentum vs simulations
fig = setup_figure([],1,1); hold on
ext{1} = 'rotationalmomentum';
if strcmp(ext, 'FOasym') || strcmp(ext{1}, 'FOasym')
    tmp = [hmm_1stlevel.cycle_metrics.TIDA, simulation{1}.TIDA, simulation_average.TIDA];
    ylb = 'Mean FO asymmetry';
  elseif strcmp(ext, 'rotationalmomentum') || strcmp(ext{1}, 'rotationalmomentum')
    tmp = [hmm_1stlevel.cycle_metrics.rotational_momentum, simulation{1}.cycle_metrics.rotational_momentum, simulation_average.cycle_metrics.rotational_momentum]./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
    ylb = {'                     clockwise \bf \leftarrow M \rightarrow \rm counterclockwise'};
  end
  %   ax(3) = axes('Position', [0.8, 0.15, 0.15, 0.7]); hold on
%   ax(3) = axes('Position', [0.1, 0.05, 0.45, .3]); hold on
  vline(0, '--k')
  
  clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
  for k=1:3
    scatter(ones(size(tmp(:,k))).*(k+(rand(size(tmp(:,k)))-0.5)/2),tmp(:,k), 'filled', 'MarkerFaceColor', clr{k}, 'MarkerFaceAlpha',0.6)
  end
  h=boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this back on top
  set(h, 'linew', 2)
  box off
  if strcmp(ext, 'FOasym') || strcmp(ext{1}, 'FOasym')
    sigstar({[1,2], [1,3]}, [hmm_1stlevel.metric_vs_sim.TIDA.prob, hmm_1stlevel.metric_vs_sim_avg.TIDA.prob])
  elseif strcmp(ext, 'rotationalmomentum') || strcmp(ext{1}, 'rotationalmomentum')
    sigstar({[1,2], [1,3]}, [hmm_1stlevel.metric_vs_sim.rotational_momentum.prob, hmm_1stlevel.metric_vs_sim_avg.rotational_momentum.prob])
  end
  %   title({'Difference between', 'observed data and', 'simulations', ''})
  xticklabels( {'obs', '1 sim', '100 sim'})
  %   xtickangle(45)
  ylabel(ylb)
  view([-270 90])
  title('Rotational Momentum vs. simulations')
  set_font(10, {'title', 'label'})
fname = [config.figdir, 'figure2_circleplot/','2supp_RotationalMomentum_obs_vs_sim'];
save_figure(fname);
