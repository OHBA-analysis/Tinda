%% Figure 2: plot circular diagram
if ~exist('simulation', 'var'), simulation = hmm_1stlevel.simulation; end
alpha_thresh = 0.05/(K.^2-K);

if whichstudy==4
  alpha_thresh = 0.0000001 * alpha_thresh;
end


clr = [{[0 0.4470 0.7410]}, {[0.8500 0.3250 0.0980]}, {[0.9290 0.6940 0.1250]}];
study = {'MEG UK (N=55)', 'MEG UK (N=55)', 'HCP (N=79)', 'Cam-CAN (N=600)', '', 'Cam-CAN (N=612)', 'WakeHen (N=19)', 'WakeHen (N=19)'};
sigpoints = hmm_1stlevel.assym_ttest.sigpoints;

%% Circle plot comparing with and without ordering
fig = setup_figure([],2,.4);

ax(1) = axes('Position', [0.05, 0, .275, .9]); %ax(1) = axes('Position', [0.025, 0, .2, 1]);
cyclicalstateplot([1 12:-1:2], hmm_1stlevel.cycle_metrics_no_ordering.mean_direction, hmm_1stlevel.assym_ttest.sigpoints,color_scheme,false);
title({'Original ordering', ''})

ax(2) = axes('Position', [0.375, 0, .275, .9]); %ax(1) = axes('Position', [0.025, 0, .2, 1]);
cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_direction, hmm_1stlevel.assym_ttest.sigpoints,color_scheme,false);
title({'Optimal ordering', '' })

ax(3) = axes('Position', [0.745, 0.15, .25, .7]); %ax(1) = axes('Position', [0.025, 0, .2, 1]);
boxplot_with_scatter([hmm_1stlevel.cycle_metrics_no_ordering.rotational_momentum, hmm_1stlevel.cycle_metrics.rotational_momentum]./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum)
xticks(1:2)
xticklabels({'original', 'optimal'})
ylabel('Cycle strength (a.u.)')
xlabel('Ordering')
box off

sigstar([1 2], hmm_1stlevel.metric_vs_no_ordering.rotational_momentum.pvals)


fname = [config.figdir, 'figure2_circleplot/','2_cyclicalpattern_bestseq_vs_no_ordering'];
save_figure(fname, false);

%% Figure 2 (per study)
fig = setup_figure([],2,.45);
clear ax

% asymmetry matrix
ax(1) = axes('Position', [0.06, 0.0, 0.235, .79]); hold on
CL = max(abs(hmm_1stlevel.cycle_metrics.mean_assym(:)))*[-1 1];
cmap = flipud(brewermap(256, 'RdBu'));
bestseq_cw = 1:12;%[bestseq(1), bestseq(end:-1:2)];
tmp = hmm_1stlevel.cycle_metrics.mean_assym(bestseq_cw, bestseq_cw);
tmp(isnan(tmp))=0;
imagesc(tmp, CL), colormap(cmap)
set(gca, 'YDir', 'reverse')
xticks(1:K), yticks(1:K), xticklabels(bestseq_cw), yticklabels(bestseq_cw), 
cb=colorbar('Location', 'SouthOutside')
cb.Ticks = [ceil(100*cb.Limits(1))/100, 0,floor(100*cb.Limits(2))/100];
% cb.Position(3) = cb.Position(3)/2;
% cb.Position(1) = cb.Position(1)+cb.Position(3)/2;

xtickangle(0)
ylabel('State interval')
xlabel({'State FO asymmetry'})
title({'Mean FO asymmetry',''})
xlim([0.5 12.5])
ylim([0.5 12.5])

% add asteriks in matrix
mask = hmm_1stlevel.assym_ttest.sigpoints(bestseq_cw,bestseq_cw);
mask2 = (hmm_1stlevel.assym_ttest.pvals(bestseq_cw, bestseq_cw) < 0.05/132) & ~mask;

[R, C] = ndgrid(1:12, 1:12);
R = R(:)+.1; C = C(:)-.15;
%rows are Y values, columns are X values !
vals = mask(:);
vals2 = mask2(:); % The points that are significant with less stringent threshold;
for k=1:length(vals)
    if vals(k)
        tx{k} = '*';
    else
        tx{k}='';
    end
    if vals2(k)
        tx2{k} = '*';
    else
        tx2{k} = '';
    end
end
text(C, R, tx, 'color', 'k', 'FontSize', 6)
% text(C, R, tx2, 'color', [.7 .7 .7], 'FontSize', 6)


% cycle plot
ax(2) = axes('Position', [.34, 0, .3, .9]); %ax(1) = axes('Position', [0.025, 0, .2, 1]);
cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_direction, hmm_1stlevel.assym_ttest.sigpoints,color_scheme,false);
title({study{whichstudy}, ''})

% rotational momentum vs perm
if whichstudy==4
  facealpha=0.2;
else
  facealpha=0.6;
end
ax(3) = axes('Position', [0.73, 0.05, 0.25, 0.8])
tmp1=hmm_1stlevel.cycle_metrics.rotational_momentum./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
boxplot_with_scatter([tmp1*nan, tmp1], clr([2,1]), facealpha)
tmp2 = mean(hmm_1stlevel.perm.rotational_momentum,2)'./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
scatter(ones(size(tmp2)).*((2-1)+(rand(size(tmp2))-0.5)/2),tmp2,25, brewermap(length(tmp2), 'Dark2'), 'filled', 'MarkerFaceAlpha',0.6)
h1 = boxplot([tmp2;nan(1,length(tmp2))]', 'Colors', 'k', 'Width', .75, 'whisker', 1000)
% set(h, {'linew'}, {2})
ylim(1.2*[min([min(tmp1), min(tmp2)]), max([max(tmp1), max(tmp2)])])
sigstar({[1,2]}, hmm_1stlevel.perm.obs_vs_perm)
xlim([-1.5 3.5])
box off
xticks([1 2])
ylabel('S (a.u.)')
xticklabels({'perm', 'obs'})
view([90,-90])
c1=brewermap(length(tmp2), 'Dark2');
for k=1:length(c1)
  cmap1{k} = c1(k,:);
end

if whichstudy==1
  y1 = 0
  x2 = 0.75;
  plot([0.25 0.75], [-0.06, mean(tmp2)], '--k')
  plot([0.25 0.75], [.11, mean(tmp2)], '--k')
elseif whichstudy==3
  y1 = 0.01;
  x2 = 0.73;
    plot([0.4 0.75], [-0.053, mean(tmp2)], '--k')
  plot([0.4 0.75], [.08, mean(tmp2)], '--k')
elseif whichstudy==4
  y1 = 0.075;
  x2 = 0.778;
elseif whichstudy==6
    y1 = 0.015;
  x2 = 0.726;
  plot([0.36 0.72], [-0.1, mean(tmp2)], '--k')
  plot([0.36 0.72], [.113, mean(tmp2)], '--k')
elseif whichstudy==7
      y1 = 0.015;
  x2 = 0.78;
  plot([0.25 0.72], [-0.02, mean(tmp2)], '--k')
  plot([0.25 0.72], [.145, mean(tmp2)], '--k')
end
text(3.2, y1-.02, 'individual subjects')
plot(3.2, y1-.03, '.', 'Color', clr{1}, 'MarkerSize', 10)
text(2.8, y1-.02, 'individual perm')
plot(2.8, y1-.03, '.', 'Color', cmap1{5}, 'MarkerSize', 10)
plot(2.8, y1-.04, '.', 'Color', cmap1{2}, 'MarkerSize', 10)
plot(2.8, y1-.05, '.', 'Color', cmap1{6}, 'MarkerSize', 10)
title('Cycle strength')
if whichstudy==6
    ylim([-.1 .25])
elseif whichstudy ==7
  ylim([-.05 .2])
else
    ylim(1.1*ylim)
end

ax(4) = axes('Position', [x2, .18, 0.15, 0.19]), hold on
facealpha = 0.6;
scatter(ones(size(tmp2)).*((2-1)+(rand(size(tmp2))-0.5)/2),tmp2,25, c1, 'filled', 'MarkerFaceAlpha',0.6)
h = boxplot([tmp2]', 'Colors', 'k', 'whisker', 1000, 'Width', .75)
set(h, {'linew'}, {2})
xlim([0.5 1.5])
ax(4).YAxis.Exponent = 0;
view([90,-90])
if whichstudy==3
    yticks([0.0075 0.01])
end
xticks([])
set_font(10, {'title', 'label'})
fname = [config.figdir, 'figure2_circleplot/','2_cyclicalpattern'];
ax(3).Children(9).Position(2) = ax(3).Children(9).Position(2)+0.02

save_figure(fname, false);


%% Supplements
%
for ext = {'rotationalmomentum'}
 
  fig = setup_figure([],2,0.7);
  clear ax
  ax(1) = axes('Position', [0.025, .45, .3, .4]); %ax(1) = axes('Position', [0.025, 0, .2, 1]);
  cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_direction, hmm_1stlevel.assym_ttest.sigpoints,color_scheme,false);
  title({'Observed', ''})
  
  ax(2) = axes('Position', [0.35, .45, .3, .4]); %  ax(2) = axes('Position', [0.275, 0, .2, 1]);
  cyclicalstateplot(bestseq,hmm_1stlevel.simulation{1}.cycle_metrics.mean_direction, hmm_1stlevel.simulation{1}.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh,color_scheme,false);
  title({'Simulated from transprob', ''})
  
  ax(4) = axes('Position', [0.675, .45, .3, .4]);%ax(4) = axes('Position', [0.525, 0, .2, 1]);
  cyclicalstateplot(bestseq,hmm_1stlevel.simulation_average.cycle_metrics.mean_direction, hmm_1stlevel.simulation_average.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh,color_scheme,false);
  title({'Simulated 100x from transprob', ''})
  
  if strcmp(ext, 'FOasym') || strcmp(ext{1}, 'FOasym')
    tmp = [hmm_1stlevel.cycle_metrics.TIDA, simulation{1}.TIDA, simulation_average.TIDA];
    ylb = 'Mean FO asymmetry';
  elseif strcmp(ext, 'rotationalmomentum') || strcmp(ext{1}, 'rotationalmomentum')
    tmp = [hmm_1stlevel.cycle_metrics.rotational_momentum, hmm_1stlevel.simulation{1}.cycle_metrics.rotational_momentum, hmm_1stlevel.simulation_average.cycle_metrics.rotational_momentum]./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
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
%}
% also plot the *real* circle plot individually
fig = setup_figure([],1,1);
cyclicalstateplot(bestseq,hmm_1stlevel.cycle_metrics.mean_direction, hmm_1stlevel.assym_ttest.sigpoints,color_scheme, false);
title({'Observed', ''})
save_figure([config.figdir, 'figure2_circleplot/','2ind_Cyclicalpattern']);

% as well as the simulated
fig = setup_figure([],1,1);
cyclicalstateplot(bestseq,hmm_1stlevel.simulation{1}.cycle_metrics.mean_direction, hmm_1stlevel.simulation{1}.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);
title({'Simulated from transprob', ''})
save_figure([config.figdir, 'figure2_circleplot/','2ind_Cyclicalpattern_simulation1']);
%
% and the simulation average
fig = setup_figure([],1,1);
cyclicalstateplot(bestseq,hmm_1stlevel.simulation_average.cycle_metrics.mean_direction, hmm_1stlevel.simulation_average.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);
title({'Simulated 100x from transprob', ''})
save_figure([config.figdir, 'figure2_circleplot/','2ind_Cyclicalpattern_simulationAvg']);
%}

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
    tmp = [hmm_1stlevel.cycle_metrics.rotational_momentum, hmm_1stlevel.simulation{1}.cycle_metrics.rotational_momentum, hmm_1stlevel.simulation_average.cycle_metrics.rotational_momentum]./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
    ylb = {'S'};
  end
  % vline(h, '--k')
  
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
    a=gca;
    a.Children(1).Position(2) = a.Children(1).Position(2)*1.1;
    a.Children(4).YData = a.Children(4).YData*.9;
    a.Children(3).Position(2) = a.Children(3).Position(2)*1;
  xticklabels( {'obs', '1 sim', '100 sim'})
    ylim([-.1 .25])
  ylabel(ylb)
  view([-270 90])
  title('Cycle strength vs. simulations')
  set_font(10, {'title', 'label'})
fname = [config.figdir, 'figure2_circleplot/','2supp_RotationalMomentum_obs_vs_sim'];
save_figure(fname, false);
