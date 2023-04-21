% Run TINDA on the vpath simulated from the transprob matrix.
%% Subject level FO assym
n_sim_perm = 100;
simulation=cell(n_sim_perm,1);
for iperm=1:n_sim_perm
  iperm
  for k=1:length(vpath)
    simulation_vpath{k} = simulateVpath(vpath{k},hmmT{k},K);
  end
  [simulation{iperm}.FO_intervals,FO_pvals,~, FO_stat] = computeLongTermAsymmetry(simulation_vpath,hmmT,K);
  simulation{iperm}.assym_permtest = [];
  simulation{iperm}.assym_permtest.stat = FO_stat;
  simulation{iperm}.assym_permtest.pvals = FO_pvals;
  simulation{iperm}.assym_ttest=[];
  K=12;
  a=[];
  for i=1:K
    for j=1:K
      [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(simulation{iperm}.FO_intervals(i,j,1,:)), squeeze(simulation{iperm}.FO_intervals(i,j,2,:)));
    end
  end
  simulation{iperm}.assym_ttest = a;
  
  
  simulation{iperm}.bestsequencemetrics = optimiseSequentialPattern(simulation{iperm}.FO_intervals);
  angleplot_sim = circle_angles(simulation{iperm}.bestsequencemetrics{1});
  simulation{iperm}.cycle_metrics = compute_tinda_metrics(config, simulation{iperm}.bestsequencemetrics{1}, angleplot_sim, simulation{iperm}.FO_intervals, simulation{iperm}.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);
  
end
hmm_1stlevel.simulation = simulation;


clear q;
for iperm=1:n_sim_perm
q(iperm,:) = hmm_1stlevel.simulation{iperm}.cycle_metrics.rotational_momentum./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
end

fig = setup_figure([],1,1);
shadedErrorBar(1:100, mean(q,2), std(q,[],2)./sqrt(config.nSj))
hline(mean(hmm_1stlevel.cycle_metrics.rotational_momentum./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), '--k')
xlabel('Simulation')
ylabel('M')
yl = ylim;
ylim([yl(1)*1.1 0])

text(0.4, .1, 'Observed', 'Units', 'normalized')
title({'Mean (+SEM) rotational momentum', 'from Markovian model'})
box off

save_figure([config.figdir, 'figure_supp_tinda_metrics/','2supp_rotational_momentum_markovian']);



%%
% Is there something in the HMM transprob matrix when we look at the
% simulation aggregate?
simulation_average=[];
for k=1:n_sim_perm
  simulation_average.FO_intervals(:,:,:,:,k) = simulation{k}.FO_intervals;
end
simulation_average.FO_intervals = mean(simulation_average.FO_intervals,5);
[FO_pvals, FO_stat] = FO_permutation_test(simulation_average.FO_intervals, K, config.nSj);
simulation_average.assym_permtest = [];
simulation_average.assym_permtest.stat = FO_stat;
simulation_average.assym_permtest.pvals = FO_pvals;
simulation_average.assym_ttest=[];
K=12;
a=[];
for i=1:K
for j=1:K
[a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(simulation_average.FO_intervals(i,j,1,:)), squeeze(simulation_average.FO_intervals(i,j,2,:)));
end
end
simulation_average.assym_ttest = a;
simulation_average.bestsequencemetrics = optimiseSequentialPattern(simulation_average.FO_intervals);
simulation_average.cycle_metrics = compute_tinda_metrics(config, simulation_average.bestsequencemetrics{1}, angleplot, simulation_average.FO_intervals, simulation_average.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);

% compare the observed circularity with the simulated one
% dat1=[];
% dat1.dimord = 'rpt_chan_time';
% dat1.label{1} = 'circularity';
% dat1.time=1;
% dat2=dat1;
% dat1.trial =  hmm_1stlevel.circularity_subject;
% dat2.trial = hmm_1stlevel.FO_stats_simulation_average.circularity_subject;
%
% cfg=[];
% cfg.method = 'montecarlo';
% cfg.statistic = 'depsamplesT';
% cfg.design = [ones(1,config.nSj), 2*ones(1,config.nSj); 1:config.nSj, 1:config.nSj];
% cfg.ivar = 1;
% cfg.uvar = 2;
% cfg.numrandomization = 100000;
% cfg.tail = 1;
% stat_c = ft_timelockstatistics(cfg, dat1, dat2);
% simulation_average.cycle_metrics.circularity_stat_obs_vs_perm = stat_c;


hmm_1stlevel.simulation_average = simulation_average;

%% Group level FO assym

% Alternatively, do the FO assym on the group level

% we can either simulate the vpath based on the group level transprob, or
% the individual vpath
% let's do the group:
nsim=100;
for sim=1:nsim
  sim
  clear hmmT_sim vpath_sim
  if 1
    if size(hmmT{1},2)>size(hmmT{1},1)
      hmmT_sim = cat(2,hmmT{:})';
    else
      hmmT_sim = cat(1,hmmT{:});
    end
    vpath_sim = simulateVpath(cat(1,vpath{:}),hmmT_sim,K);
  else
    hmmT_sim = hmmT;
    for k=1:length(vpath)
      vpath_sim{k} = simulateVpath(vpath{k},hmmT_sim{k},K);
    end
    vpath_sim = cat(1,vpath_sim{:});
    hmmT_sim = cat(1,hmmT_sim{:});
  end
  [FO_group_sim{sim}, ~, ~] = computeLongTermAsymmetry({vpath_sim},{hmmT_sim},K);
end
FO_group_sim = cat(4,FO_group_sim{:});
hmm_1stlevel.FO_simulation_group = FO_group_sim;

%% Compare the observed metrics with the simulated ones
% per subject measures
cfg=[];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.design = [ones(1,config.nSj), 2*ones(1,config.nSj); 1:config.nSj, 1:config.nSj];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 100000;
cfg.correcttail = 'prob';

dat1=[];
dat1.dimord = 'rpt_chan_time';
dat1.label{1} = 'metric';

measures = {'FO_assym_rv_coef', 'FO_assym_subject_fit', 'TIDA', 'rotational_momentum', 'circularity_subject', 'TIDA_perstate', 'rotational_momentum_perstate'};
for im = measures
  m = im{1};
  if strcmp(m ,'FO_assym_rv_coef') || ...
      strcmp(m, 'FO_assym_subject_fit') || strcmp(m, 'TIDA') ||... 
      strcmp(m, 'TIDA_perstate') || strcmp(m, 'circularity') % these are all positive numbers
    cfg.tail = 1;
  else
    cfg.tail = -1; % rotational momentum should have a tail of -1
  end
  
  dat1.time=1:size(hmm_1stlevel.cycle_metrics.(m),2);
  dat1.trial = [];
  dat2=dat1;
  
  dat1.trial(:,1,:) = hmm_1stlevel.cycle_metrics.(m);
  dat2.trial(:,1,:) = hmm_1stlevel.simulation{1}.cycle_metrics.(m);
  
  hmm_1stlevel.metric_vs_sim.(m) = ft_timelockstatistics(cfg, dat1, dat2);
  
  dat2.trial(:,1,:) = hmm_1stlevel.simulation_average.cycle_metrics.(m);
  hmm_1stlevel.metric_vs_sim_avg.(m) = ft_timelockstatistics(cfg, dat1, dat2);
end
