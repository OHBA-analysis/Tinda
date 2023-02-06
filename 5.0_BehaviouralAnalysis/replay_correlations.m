%{
for whichstudy = [1, 2]
  is=[];
  if whichstudy==1
    MFStartup;
    is.dirname = 'Replaydata4Cam/';
      is.sessionsOfInterest = [2,3,7,9]; % lci, lci, rst, rst
  else
    MFStartup_studyII;
    is.dirname='StrLearn_MEGexp/';
      is.sessionsOfInterest = [1,2,3,4,8]; % rst, lci, lci, lci, rst
  end
  load([is.rootBehav, 'Exptrial.mat'])
  is.usePrecomputed = true;
  is.iRun = [1,2];
  is.topPercentile=1;
  is.lambda = 5; % hardcoded best lambda for all subjects
  is.whenInMS = is.tss*is.msPerSample;
  is.TOI = 200; % 200ms post stimulus.
%   is.nStim = 8;
%   is.ntrlAll=120; % per session (includes upside down images)
%   is.ntrlValid=96; % per session
  is.nSes = 3;
  is.nChan = 273;
  is.nTime = 401;
  is.nPerm = 1000;
  is.K = 12;
  is.plot=false;
  fsample=250;
  is.t_window=fsample/2;
  is.useMask=0;
  config.resultsdir = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study1/';
  color_scheme = colorscheme(1);
  config.figdir = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study1/figures/replay/';
  
  
  %% Get the Study specific HMM results
  fname = fullfile(is.studydir, 'Neuron2020Analysis/', sprintf('Study%d',whichstudy), "hmm_1to45hz", sprintf("hmm%dusingtemplate_parc_giles_symmetric__pcdim80_voxelwise_embed13_K12_big1_dyn_modelhmm.mat", whichstudy+1));
  replayflag=1;
  load(fname, 'hmm')
  if replayflag% use the replay ordering
    load(fname, 'new_state_ordering')
  else
    load([config.resultsdir, 'coherence_state_ordering.mat'])
  end
  hmm = hmm_permutestates(hmm,new_state_ordering);
  Gamma = hmm.gamma;
  K=hmm.K;
  
  % and the timings
  fname=fullfile(is.studydir,'Neuron2020Analysis/', sprintf('Study%d',whichstudy), 'hmm_1to45hz/hmm_parc_giles_symmetric__pcdim80_voxelwise_embed13.mat');
  load(fname,'hmmT','subj_inds')
  
  % load in the subject masks
  datadir = fullfile(is.studydir,'Neuron2020Analysis/', sprintf('Study%d',whichstudy), 'bfnew_1to45hz/');
  [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datadir);
  gs{whichstudy}=goodsamples;
  tp{whichstudy}=triggerpoints;
  
  nSes = length(goodsamples);
  K=size(Gamma,2);
  
  if length(triggerpoints{1})>40000
    Fs = 250;
  else
    Fs = 100;
  end
  
  t=[-is.t_window:is.t_window]./Fs;
  method=1;
  
  
  scan_T = cell2mat(hmmT);
  scan_T_studies{whichstudy}=scan_T;
  R = [[1;1+cumsum(scan_T(1:end-1))'],cumsum(scan_T(1:end))'];
  if ~all(R(2:end,1)==R(1:end-1,2)+1)
    % correct any uncorrected offsets in R:
    R(2:end,1) = R(1:end-1,2)+1;
  end

  for iSes=1:nSes
    iSes
    % convert Gamma to padded timecourse
    Gam_long = NaN(length(goodsamples{iSes}),K);
    Gam_long(goodsamples{iSes},:) = Gamma(R(iSes,1):R(iSes,2),:);
    Gam_long = Gam_long(triggerpoints{iSes},:);
    
    G{whichstudy}{iSes} = Gam_long - nanmean(Gam_long);
    
    temp1 = zeros(length(goodsamples{iSes}),1);
    temp1(goodsamples{iSes}) = hmm.statepath(R(iSes,1):R(iSes,2));
    temp1 = temp1(triggerpoints{iSes},:);
    vpath_studies{whichstudy}{iSes} = temp1;
    T_studies{whichstudy}{iSes} = [length(vpath_studies{whichstudy}{iSes})];
  end
  tmp = load([is.resultsdir, 'replayidx/replay_idx.mat'])
  if exist('topidx', 'var') && ~isempty(topidx)
    topidx = cat(2, topidx, tmp.topidx);
  else
    topidx = tmp.topidx;
  end
  hmm_study{whichstudy}=hmm;
  wd = '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/';
      if whichstudy==1
        load([wd,'GenericReplayData/STUDYI_ReplayOnset/StrReplayOnset'],'ToRall'); % the replay scores for first resting state
        replayScores{whichstudy}(1:2:(2*21),:) = ToRall;
        load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
        replayScores{whichstudy}(2:2:(2*21),:) = ToRall;
    else
        load([wd,'GenericReplayData/STUDYII_ReplayOnset/STUDYII_PreplayOnset'],'ToRall'); % the replay scores for first resting state
        replayScores{whichstudy}(1:2:(2*22),:) = ToRall;
        load([wd,'GenericReplayData/STUDYII_ReplayOnset/STUDYII_ReplayOnset'],'ToRall');
        replayScores{whichstudy}(2:2:(2*22),:) = ToRall;
      end
end
%}
get_replayData


%% Combine the two studies.
% Load the replay timing info
rep_scores = cat(1,rep_scores{:});
vpath = cat(2,vpath_studies{:});
T = cat(2, T_studies{:});
G = cat(2, G{:});

nSes = length(vpath);
nSj = nSes/2;
config.nSj=nSj;


%% Get the sequence from the canonical RS
replayflag=0;
if replayflag
  optimalseqfile = [is.studydir, 'Neuron2020Analysis/CanonicalRS/250Hz/hmm_1to45hz/', 'bestseq',int2str(whichstudy),'.mat'];
else
  optimalseqfile = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study1/bestseq1_coherence.mat';
end
if ~isfile(optimalseqfile)
  bestsequencemetrics = optimiseSequentialPattern(FO);
  save(optimalseqfile,'bestsequencemetrics');
else
  load(optimalseqfile);
end
bestseq = bestsequencemetrics{1};
% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
dp_2d = zeros(12,2);
for i3=1:12
  disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
  dp_2d(bestseq(i3),:) = [real(disttoplot_manual(bestseq(i3))), imag(disttoplot_manual(bestseq(i3)))];
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));


%% create a circle path projection
disttoplot_manual = zeros(12,2);
for i=1:12
  temp = exp(sqrt(-1)*(i+2)/12*2*pi);
  disttoplot_manual(bestseq(i),:) = [real(temp),imag(temp)];
end

for iSj=1:nSj
  task_evoked_circle(:,:,iSj) = G_evoked(:,:,iSj)*disttoplot_manual;
end
mu_2d = 3*mean(task_evoked_circle,3);
%plot(mu_2d,'k');
%
figure();
cyclicalstatesubplot(bestseq,zeros(12,12),zeros(12,12))
hold on;
t_ofinterest = -0.5:1/250:0.5;
for t=1:size(mu_2d,1)
  plot(mu_2d(t,1),mu_2d(t,2),'.','Color',[1 1 1]*t./(length(t_ofinterest)));
  pause(0.001)
end

config.resultsdir = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study1/';
config.figdir = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study1/figures/replay/';

%% Compute TINDA on the resting state
% run tinda on the continuous data - also run when accounting for replay
% evoked gammas (using the glm)
[FO_replay,pvals_replay,~,stat_replay] = computeLongTermAsymmetry(vpath,T,K);

% Run tinda statistics without correcting for evoked response
FO_replay = (FO_replay(:,:,:,1:2:end)+FO_replay(:,:,:,2:2:end))/2;
[pvals_replay, stat_replay] = FO_permutation_test(FO_replay, K, nSj);
a=[];
for i=1:K
  for j=1:K
    [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(FO_replay(i,j,1,:)), squeeze(FO_replay(i,j,2,:)));
  end
end
replay.FO_intervals = FO_replay;
replay.assym_ttest = a;
replay.assym_permtest.stat = stat_replay;
replay.assym_permtest.pvals = pvals_replay;


% find the best ordering
replay.bestseq.modelhmm = bestseq;
replay.angleplot.modelhmm = angleplot;
bestsequencemetrics = optimiseSequentialPattern(FO_replay);
replay.bestseq.replay = bestsequencemetrics{1};
disttoplot_manual = zeros(12,1);
for i3=1:12
  disttoplot_manual(replay.bestseq.replay(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
replay.angleplot.replay = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));

% get some metrics
seq = {'modelhmm', 'replay'};
for k=1:2
  replay.cycle_metrics.(seq{k}) = compute_tinda_metrics(config, replay.bestseq.(seq{k}), replay.angleplot.(seq{k}), replay.FO_intervals, replay.assym_ttest.pvals<0.05, color_scheme);
end

%{
% and with correctiong for evoked response
FO_replay_residual = (FO_replay_residual(:,:,:,1:2:end)+FO_replay_residual(:,:,:,2:2:end))/2;
[pvals_replay_residual, stat_replay_residual] = FO_permutation_test(FO_replay_residual(:,:,:,:), K, nSj);
hmm_1stlevel.replay_res.FO_intervals = FO_replay_residual;
hmm_1stlevel.replay_res.FO_stat = stat_replay_residual;
hmm_1stlevel.replay_res.FO_pvals = pvals_replay_residual;
bestsequencemetrics = optimiseSequentialPattern(FO_replay_residual);
bestseq_replay_residual = bestsequencemetrics{1};
disttoplot_manual = zeros(12,1);
for i3=1:12
  disttoplot_manual(bestseq_replay_residual(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot_replay_residual = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));
hmm_1stlevel.replay_res.cycle_metrics = compute_tinda_metrics(config, bestseq_replay_residual, angleplot, FO_replay_residual, pvals_replay_residual<0.05, color_scheme);
%}
%% Make cycle plots
for k=1:2
  setup_figure([], 2,0.4)
  ax(1) = axes('Position', [0.05 0.05 0.3 .8]), box off
  cyclicalstateplot(replay.bestseq.(seq{k}),replay.cycle_metrics.(seq{k}).mean_direction,replay.assym_ttest.pvals<0.05, color_scheme, false);
  title({'p<0.05 (uncorrected)', ''})
  ax(1) = axes('Position', [0.425 0.1 0.2 .8]), box off, hold on
  clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
  tmp = replay.cycle_metrics.(seq{k}).rotational_momentum./replay.cycle_metrics.(seq{k}).max_theoretical_rotational_momentum;
  scatter(ones(size(tmp(:,1))).*(1+(rand(size(tmp(:,1)))-0.5)/2),tmp(:,1), 'filled', 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha',0.6)
  h=boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this back on top
  set(h, 'linew', 2)
  box off
  ylabel('M')
  xticks([])
  title({'', sprintf('M = %0.3f', mean(replay.cycle_metrics.(seq{k}).rotational_momentum)./replay.cycle_metrics.(seq{k}).max_theoretical_rotational_momentum), ''})
  ax(3) = axes('Position', [0.65 0.05 0.3 .8]), box off
  cyclicalstateplot(replay.bestseq.(seq{k}),replay.cycle_metrics.(seq{k}).mean_direction,replay.assym_ttest.pvals<0.01, color_scheme, false);
  title({'p<0.01 (uncorrected)', ''})
  save_figure([config.figdir, sprintf('2_Cyclicalpattern_order_%s', seq{k})],[],false);
end




%% Plot the distribution of replay events on the cycle plot
for perc=[1,5]
  topidx = replay_scores2index(vpath, rep_scores, [], [], perc);
  
  
  G_evoked = [];
  for iSes=1:nSes
      tmp = zeros(251, 12, length(topidx{iSes}));
      for k=1:length(topidx{iSes})
          ix = topidx{iSes}(k);
          tmp(:,:,k) = G{iSes}(ix-125:ix+125,:);
      end
      G_evoked(:,:, iSes) = nanmean(tmp,3);
  end
  G_evoked = squeeze(nanmean(reshape(G_evoked, 251, 12, 2, []),3));

  
  % figure;
  q{1} = squeeze(mean(betas_norm{1}))'; %squeeze(betas_norm{1}(126,:,:))';
  q{2} = squeeze(mean(betas_norm{2}))'; %squeeze(betas_norm{2}(126,:,:))';
  q{3} = [q{1}; q{2}];
  ext = {'study1', 'study2', 'bothstudies'};
  for ik=1:2
    setup_figure([], 2,0.35);
    for k=1:3
      ax(k) = axes('Position', [0.05+(k-1)*0.3 0.05 0.25 0.8]); hold on
      CL = [min(mean(q{k})), max(mean(q{k}))];
      % for two sequence orders
      
      cyclicalstatesubplot(replay.bestseq.(seq{ik}), zeros(12), zeros(12), color_scheme)
      cyclicalstate_distributionplot(replay.bestseq.(seq{ik}), mean(q{k})', CL, [], false,[],[],1)
      title(sprintf('%s', ext{k}))
      plot(0, 0, '.k', 'MarkerSize', 20)
    end
    suptitle('Replay distribution')
    save_figure([config.figdir, sprintf('2_spiderplot_order_%s_perc%d',(seq{ik}), perc)],[],false);
  end
  
  
  %% Redo TINDA where Replay is state 13
  color_scheme_K13 = [color_scheme, {[0 0.4470 0.7410]}];
  vpath_K13 = vpath;
  for k=1:length(vpath)
    vpath_K13{k}(topidx{k}) = K+1;
  end
  [FO_replay_orig,~,t_intervals,~] = computeLongTermAsymmetry(vpath_K13,T,K+1);
  
  FO_replay = (FO_replay_orig(:,:,:,1:2:end)+FO_replay_orig(:,:,:,2:2:end))/2;
  [pvals_replay, stat_replay] = FO_permutation_test(FO_replay, K+1, nSj);
  replay.K13.FO_intervals=FO_replay;
  a=[];
  for i=1:K+1
    for j=1:K+1
      [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(FO_replay(i,j,1,:)), squeeze(FO_replay(i,j,2,:)));
    end
  end
  replay.K13.assym_ttest = a;
  replay.K13.assym_permtest.stat = stat_replay;
  replay.K13.assym_permtest.pvals = pvals_replay;
  
  % Find best ordering by finding the best gap for state 13
  assym = squeeze((replay.K13.FO_intervals(:,:,1,:)-replay.K13.FO_intervals(:,:,2,:))./mean(replay.K13.FO_intervals,3));
  for ik=1:2 % two sequence orders
    bestseq_K13 = [];
    minR = 0;
    for k=1:K
      bs = replay.bestseq.(seq{ik});
      bs = [bs(1:k), 13, bs(k+1:end)];
      disttoplot_manual = zeros(13,1);
      for i3=1:13
        disttoplot_manual(bs(i3)) = exp(sqrt(-1)*i3/13*2*pi);
      end
      ap = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));
      tmp_rotmom = compute_rotational_momentum(ap, assym);
      if mean(tmp_rotmom)<minR
        minR = mean(tmp_rotmom);
        replay.K13.bestseq.(seq{ik}) = bs;
        replay.K13.angleplot.(seq{ik}) = ap;
      end
    end
    replay.K13.cycle_metrics.(seq{ik}) = compute_tinda_metrics(config, [], replay.K13.angleplot.(seq{ik}), replay.K13.FO_intervals, replay.K13.assym_ttest.pvals<0.05, color_scheme_K13);
  end
  
  %% K13 cycle plot
  for k=1:2
    setup_figure([], 2,0.4)
    ax(1) = axes('Position', [0.05 0.05 0.3 .8]), box off
    cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,replay.K13.assym_ttest.pvals<0.05, color_scheme_K13, false);
    title({'p<0.05 (uncorrected)', ''})
    ax(1) = axes('Position', [0.425 0.1 0.2 .8]), box off, hold on
    clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
    tmp = replay.K13.cycle_metrics.(seq{k}).rotational_momentum./replay.K13.cycle_metrics.(seq{k}).max_theoretical_rotational_momentum;
    scatter(ones(size(tmp(:,1))).*(1+(rand(size(tmp(:,1)))-0.5)/2),tmp(:,1), 'filled', 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha',0.6)
    h=boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this back on top
    set(h, 'linew', 2)
    box off
    ylabel('M')
    xticks([])
    title({'', sprintf('M = %0.3f', mean(replay.K13.cycle_metrics.(seq{k}).rotational_momentum)./replay.K13.cycle_metrics.(seq{k}).max_theoretical_rotational_momentum), ''})
    ax(3) = axes('Position', [0.65 0.05 0.3 .8]), box off
    cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,replay.K13.assym_ttest.pvals<0.01, color_scheme_K13, false);
    title({'p<0.01 (uncorrected)', ''})
    save_figure([config.figdir, sprintf('2_Cyclicalpattern_K13_order_%s_perc%d', seq{k}, perc)], [], false);
    
    
    % also plot just the arrows for state 13
    fig = setup_figure([],1,1);
    sigpoints = [zeros(12), ones(12,1); ones(1,13)].*(replay.K13.assym_ttest.pvals<0.05/12);
    cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,sigpoints, color_scheme_K13, false);
    save_figure([config.figdir, sprintf('2_Cyclicalpattern_K13_replay_order_%s_perc%d', seq{k}, perc)], [], false);
    
    % and with state 13 highlighted
        fig = setup_figure([],1,1);
    sigpoints1 = (replay.K13.assym_ttest.pvals<0.01);
    sigpoints2 = (replay.K13.assym_ttest.pvals<0.05/12);
    cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,[ones(12), zeros(12,1); zeros(1,13)].*sigpoints1, color_scheme_K13, false,[],[],[.8 .8 .8]);
    cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,[zeros(12), ones(12,1); ones(1,13)].*sigpoints2, color_scheme_K13, false);
    save_figure([config.figdir, sprintf('2_Cyclicalpattern_K13_replay_highlighted_order_%s_perc%d', seq{k}, perc)], [], false);
    
  end
  
  %% Replay figure
    fig = setup_figure([],1.5,.5);
    k=2;
    ax(1) = axes('Position',[0.05, 0.05, 0.4, 0.9])
    cyclicalstateplot(replay.bestseq.(seq{k}),replay.cycle_metrics.(seq{k}).mean_direction,replay.assym_ttest.pvals<0.01, color_scheme, false);
    
    sigpoints1 = (replay.K13.assym_ttest.pvals<0.01);
    sigpoints2 = (replay.K13.assym_ttest.pvals<0.05/12);
    ax(2) = axes('Position',[0.55, 0.05, 0.4, 0.9]), hold on

    cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,[ones(12), zeros(12,1); zeros(1,13)].*sigpoints1, color_scheme_K13, false,[],[],[.8 .8 .8]);
    cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,[zeros(12), ones(12,1); ones(1,13)].*sigpoints2, color_scheme_K13, false);
    save_figure([config.figdir, sprintf('4_Cyclicalpattern_replay_%s_perc%d', seq{k}, perc)], [], false);
    
    for ii=1:12, dum{ii}=ii; end
    figure; subplot(2,1,1), boxplot_with_scatter(squeeze(replay.K13.cycle_metrics.replay.FO_assym(13,:,:))'), title('tinda on replay intervals')
    sigstar(dum, replay.K13.assym_ttest.pvals(13,1:12)*24)
    subplot(2,1,2), boxplot_with_scatter(squeeze(replay.K13.cycle_metrics.replay.FO_assym(:,13,:))'), title('tinda on state intervals')
    sigstar(dum, replay.K13.assym_ttest.pvals(1:12,13)*24)
    suptitle('Bonferroni corrected p<(0.05/24)')
    save_figure([config.figdir, sprintf('4_replay_FO_assym_%s_perc%d', seq{k}, perc)], [], false);

  
  
  %% Replay spider plot
  ttl = {{'TINDA on' 'replay intervals'}, {'TINDA on' 'state intervals'}};
  CM = (brewermap(128,'RdBu'));
  ax=[];
  q=[];
  for k=1:12
    q{1}(k) = replay.K13.assym_ttest.stat(k,13).tstat;% tinda on state intervals
    q{2}(k) = -replay.K13.assym_ttest.stat(13,k).tstat; % tinda on replay intervals (positive is replay leads into state)
  end
  CL = max(abs(cat(2,q{:})))'*[-1 1];
  
  for ik = 1:2
    setup_figure([], 1.5,0.5)
    for j=1:2
      CL = [0, max(abs(q{j}))];
      ax(j) = axes('Position', [0.1+(j-1)*0.5 0.0, 0.33 0.9]);
      hold on
      clr = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
      q1 = q{j}; q1(q1<0)=0;
      q2 = q{j}; q2(q2>0)=0;
      lim=130;
      cyclicalstatesubplot(replay.bestseq.(seq{ik}), zeros(12), zeros(12), color_scheme)
      cyclicalstate_distributionplot(replay.bestseq.(seq{ik}), q1', CL, CM, false, clr{1},lim)
      cyclicalstate_distributionplot(replay.bestseq.(seq{ik}), -q2', CL, CM, false, clr{2}, lim)
      plot(0, 0, '.k', 'MarkerSize', 20)
      title(ttl{j})
    end
    ax(3) = axes('Position', [0.525 0.45, 0.1 0.5]);
    plot(1:10, rand(3,10)), xlim([-1 0]);
    legend({['state leads',newline, 'into replay'], ['replay leads',newline, 'into state']}, 'Location', 'north')
    legend('boxoff')
    axis off
    save_figure([config.figdir, sprintf('2_spiderplot_asymmetry_%s_perc%d', seq{ik}, perc)],[],false);
  end
  
  
  
  %% Histogram
  percentiles = 5:5:95;
  for i=1:length(percentiles)
      [FO_replay_p(:,:,:,:,i),~,t_intervals_p{i},~] = computeLongTermAsymmetry(vpath_K13,T,K+1, [percentiles(i)-4, percentiles(i)+5]);
  end
  
  
  
  
  save([config.resultsdir, 'tinda_replay_perc', num2str(perc)], 'replay')
end




%% How long before/after a replay event is a state active?
for ik=1:12
  for iSes=1:nSes
    ix{iSes} = [];
    for k=1:length(topidx{iSes})
      try
        ix{iSes}(k) = find(vpath{iSes}(topidx{iSes}(k):end)==ik,1)./250;
        ix2{iSes}(k) = (topidx{iSes}(k)-max(find(vpath{iSes}(1:topidx{iSes}(k))==ik)))./250;
      catch
        ix{iSes}(k) = nan;
        ix2{iSes}(k) = nan;
      end
    end
  end
  tim(ik,:) = cellfun(@nanmean, ix);
  tim2(ik,:) = cellfun(@nanmean, ix2);
end
tim=(tim(:,1:2:end) + tim(:,2:2:end))./2;
tim2=(tim2(:,1:2:end) + tim2(:,2:2:end))./2;

%}