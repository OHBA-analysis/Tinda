% Run TINDA on the vpath simulated from the transprob matrix.
%% Subject level FO assym
n_sim_perm = 100;
simulation=cell(n_sim_perm,1);
for iperm=1:n_sim_perm
  iperm
  for k=1:length(vpath)
    simulation_vpath{k} = simulateVpath(vpath{k},hmmT{k},K);
  end
  [simulation{iperm}.FO_intervals,simulation{iperm}.FO_pvals,simulation{iperm}.t_intervals, simulation{iperm}.FO_stat] = computeLongTermAsymmetry(simulation_vpath,hmmT,K);
  simulation{iperm}.bestsequencemetrics = optimiseSequentialPattern(simulation{iperm}.FO_intervals);
  simulation{iperm}.cycle_metrics = compute_tinda_metrics(config, simulation{iperm}.bestsequencemetrics{2}, angleplot, simulation{iperm}.FO_intervals, simulation{iperm}.FO_pvals<alpha_thresh, color_scheme, false);
end
hmm_1stlevel.simulation = simulation;


% Is there something in the HMM transprob matrix when we look at the
% simulation aggregate? 
simulation_average=[];
for k=1:n_sim_perm
  simulation_average.FO_intervals(:,:,:,:,k) = simulation{k}.FO_intervals;
end
simulation_average.FO_intervals = mean(simulation_average.FO_intervals,5);
[simulation_average.FO_pvals, simulation_average.FO_stat] = FO_permutation_test(simulation_average.FO_intervals, K, config.nSj);
simulation_average.bestsequencemetrics = optimiseSequentialPattern(simulation_average.FO_intervals);
simulation_average.cycle_metrics = compute_tinda_metrics(config, simulation_average.bestsequencemetrics{2}, angleplot, simulation_average.FO_intervals, simulation_average.FO_pvals<alpha_thresh, color_scheme, false);

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
nsim=1000;
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
FO_group_sim = cat(4,FO_all_sim{:});
hmm_1stlevel.FO_simulation_group = FO_group_sim;
