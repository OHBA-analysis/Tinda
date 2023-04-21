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
    
    [FO_replay_orig,~,t_intervals,~] = computeLongTermAsymmetry(vpath,T,K,[],[],topidx,false);
    
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
        sigpoints = [zeros(12), ones(12,1); ones(1,13)].*(replay.K13.assym_ttest.pvals<0.05/24);
        cyclicalstateplot(replay.K13.bestseq.(seq{k}),replay.K13.cycle_metrics.(seq{k}).mean_direction,sigpoints, color_scheme_K13, false);
        save_figure([config.figdir, sprintf('2_Cyclicalpattern_K13_replay_order_%s_perc%d', seq{k}, perc)], [], false);
        
        % and with state 13 highlighted
        fig = setup_figure([],1,1);
        sigpoints1 = (replay.K13.assym_ttest.pvals<0.01);
        sigpoints2 = (replay.K13.assym_ttest.pvals<0.05/24);
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
    
    
    
    %% TINDA run on quartiles
    k=2;
    percentiles = 10:10:90;
    percentiles
    for i=1:length(percentiles)
        [FO_p(:,:,:,:,i),~,t_intervals_p{i},~] = computeLongTermAsymmetry(vpath,T,K, [percentiles(i)-9, percentiles(i)+10],[],topidx,false,6);% do it in percentiles
    end
    % if we want to look at absolute sum of occurrence (not normalized by
    % interval length) - do the following
    if 0
        for k=1:9
            for k2=1:86
                for k3=1:13
                    FO_p(k3,:,:,k2,k) = FO_p(k3,:,:,k2,k)*sum(intervals{k}{k2,k3});
                end
            end
        end
    end
    
    
    FO_replay_p = (FO_p(:,:,:,1:2:end,:) + FO_p(:,:,:,2:2:end,:))/2;
    
    
    for i=1:length(percentiles)
        t_int(:,i) = mean(cellfun(@nanmean, t_intervals_p{i}));
    end
    replay.K13.bin_perc.FO = FO_p;
    replay.K13.bin_perc.percentiles=percentiles;
    replay.K13.bin_perc.t_intervals=t_intervals_p;
    replay.K13.bin_perc.t_intervals_mean=t_int;
    
    % plot group mean
    q = squeeze(mean(FO_replay_p,4));
    fig=setup_figure([],2,1);  suptitle({'Binned and percentiled FO',''})
    cmap = inferno;
    bins = 1:size(FO_replay_p,3);
    subplot(3,2,1); imagesc(bins, 1:9, squeeze(q(1,2,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title('PAN FO in DMN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(1,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','4','5','6','\bfDMN\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,2); imagesc(bins, 1:9, squeeze(q(13,2,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title('PAN FO in replay intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(13,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,3); imagesc(bins, 1:9, squeeze(q(2,1,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title('DMN FO in PAN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(2,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfPAN\rm', '1','2','3','4','5','6','\bfPAN\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,4); imagesc(bins, 1:9, squeeze(q(13,1,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title('DMN FO in replay intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(13,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,5); imagesc(bins, 1:9, squeeze(q(1,13,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title('Replay FO in DMN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(1,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','4','5','6','\bfDMN\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,6); imagesc(bins, 1:9, squeeze(q(2,13,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title('Replay FO in PAN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(2,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfPAN\rm', '1','2','3','4','5','6','\bfPAN\rm'}), xtickangle(45), vline(3.5, '-k')
    save_figure([config.figdir, sprintf('2_FO_binned_percentiled_perc%d', perc)],[],false);
    
    % normalize group mean by row
    qq = squeeze(nanmean(FO_replay_p./mean(FO_replay_p,3),4));
    cmap = flipud(brewermap(256, 'RdBu'));
    fig=setup_figure([],2,1);  suptitle({'Binned and percentiled FO (normalized by percentile FO)', ''})
    subplot(3,2,1); imagesc(bins, 1:9, squeeze(qq(1,2,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title('PAN FO in DMN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(1,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','4','5','6','\bfDMN\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,2); imagesc(bins, 1:9, squeeze(qq(13,2,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title('PAN FO in replay intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(13,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,3); imagesc(bins, 1:9, squeeze(qq(2,1,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title('DMN FO in PAN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(2,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfPAN\rm', '1','2','3','4','5','6','\bfPAN\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,4); imagesc(bins, 1:9, squeeze(qq(13,1,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title('DMN FO in replay intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(13,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,5); imagesc(bins, 1:9, squeeze(qq(1,13,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title('Replay FO in DMN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(1,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','4','5','6','\bfDMN\rm'}), xtickangle(45), vline(3.5, '-k')
    subplot(3,2,6); imagesc(bins, 1:9, squeeze(qq(2,13,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title('Replay FO in PAN intervals'), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(2,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfPAN\rm', '1','2','3','4','5','6','\bfPAN\rm'}), xtickangle(45), vline(3.5, '-k')
    save_figure([config.figdir, sprintf('2_FO_binned_percentiled_normalized_perc%d', perc)],[],false);
    
    figure; cnt=1;
    for k=1:13
        for ik=1:13
            if ik~=k
                subplot(13,13,cnt)
                imagesc(bins, 1:9, squeeze(qq(k,ik,:,:))');  c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, colormap(cmap)
                yticks(1:9), yticklabels(round(t_int(k,:)))
            end
            cnt=cnt+1;
        end
    end
    save_figure([config.figdir, sprintf('2_FO_binned_percentiled_normalized_allstates_perc%d', perc)],[],false);
    
    % Zoom to the longest interval times
    figure; cnt=1;
    for k=1:13
        for ik=1:13
            if ik~=k
                subplot(13,13,cnt)
                imagesc(bins, 1:4, squeeze(qq(k,ik,:,6:9))');  c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, colormap(cmap)
                yticks(1:4), yticklabels(round(t_int(k,6:9))), xticks(1:6)
            end
            cnt=cnt+1;
        end
    end
    save_figure([config.figdir, sprintf('2_FO_binned_percentiled_normalized_allstates_zoomed_perc%d', perc)],[],false);
    
    % Just the replay intervals
    fig=setup_figure([],2,1); cnt=1;
    k=13;
    for ik=1:12
        if ik~=k
            subplot(4,3,cnt)
            imagesc(bins, 1:4, squeeze(qq(k,ik,:,6:9))');  c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, colormap(cmap)
            yticks(1:4), yticklabels(round(t_int(k,6:9))),
            xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45)
        end
        title(['State ', num2str(cnt)]), ylabel('mean IT (ms)')
        cnt=cnt+1;
    end
    suptitle('state FO in replay interval')
    save_figure([config.figdir, sprintf('2_FO_binned_percentiled_normalized_replayIT_zoomed_perc%d', perc)],[],false);
    
    
    % Just the state intervals
    fig=setup_figure([],2,1); cnt=1;
    for k=1:12;
        ik=13;
        if ik~=k
            subplot(4,3,cnt)
            imagesc(bins, 1:4, squeeze(qq(k,ik,:,6:9))');  c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, colormap(cmap)
            yticks(1:4), yticklabels(round(t_int(k,6:9))),
            xticks([0.5 1:6, 6.5]),xticklabels({'\bfState\rm', '1','2','3','4','5','6','\bfState\rm'}), xtickangle(45)
        end
        title(['State ', num2str(cnt)]), ylabel('mean IT (ms)')
        cnt=cnt+1;
    end
    suptitle('Replay FO in state interval')
    save_figure([config.figdir, sprintf('2_FO_binned_percentiled_normalized_stateIT_zoomed_perc%d', perc)],[],false);
    
    
    
    %% Tinda run on specific interval durations
    k=2;
    lags =  [0,4,16,50,100,200,500,1000,5000]/1000; % ms
    lags_samples = lags*config.sample_rate;

    [FO_p,~,t_intervals_p,~] = computeLongTermAsymmetry(vpath,T,K, lags_samples,[],topidx,false,[],'abs');% do it in absolute terms

    % if we want to look at absolute sum of occurrence (not normalized by
    % interval length) - do the following
    if 0
        for k=1:9
            for k2=1:86
                for k3=1:13
                    FO_p(k3,:,:,k2,k) = FO_p(k3,:,:,k2,k)*sum(intervals{k}{k2,k3});
                end
            end
        end
    end
    
    
    FO_replay_p = (FO_p(:,:,:,1:2:end,:) + FO_p(:,:,:,2:2:end,:))/2;
    a=[];
    for i=1:K+1
        for j=1:K+1
            for ip=1:length(lags)-1
                [a.h(i,j,ip), a.pvals(i,j,ip), a.ci(i,j,ip,:), a.stat(i,j,ip)] = ttest(squeeze(FO_replay_p(i,j,1,:,ip)), squeeze(FO_replay_p(i,j,2,:,ip)));
                tstat(i,j,ip) = a.stat(i,j,ip).tstat;
            end
        end
    end
    replay.K13.perc.assym_ttest = a;
    
    replay.K13.perc.FO = FO_p;
    replay.K13.perc.lags=lags;
    replay.K13.perc.t_intervals=t_intervals_p;

    
    % get some extra variables for creating a bubblechart in a newer matlab
    % version
    bubble.mean_direction = squeeze(mean(FO_replay_p(:,:,1,:,:) - FO_replay_p(:,:,2,:,:),4));
    bestseq_nott = load([config.basedir, 'Study1/bestseq1_coherence.mat']);
    bestseq_nott = bestseq_nott.bestsequencemetrics{1};
    bubble.bestseq_nott = bestseq_nott;
    bubble.angleplot_nott = circle_angles(bestseq_nott); 
    bubble.tstat=tstat;
    bubble.cmap = flipud(brewermap(256,'RdBu'));
    replay.K13.perc.bubbleplot = bubble;
   
    
    save([config.resultsdir, 'tinda_replay_perc', num2str(perc)], 'replay')
end


%% Look at intervals between replay/DMN/PAN
states=[1,2,13];
sum_occ = cell(length(states), length(states), nSes);
sum_occ_bin=sum_occ; intervals = sum_occ;
FO_interplay_ses = nan(length(states), length(states), nSes);
nbins=3;
for k1=states
    k1
    for k2=setdiff(states,k1)
        k2
        k3=setdiff(states, [k1,k2]);
        for iSes = 1:nSes
            v = vpath{iSes};
            rep = zeros(size(v));
            rep(topidx{iSes})=1;
            t = T{iSes};
            tmp_sum_occ=[];tmp_sum_occ_bin=[];
            intv=[];
            ix=0;
            w=0;
            while w==0
                % times are in samples (Fs=1/250)
                if k1>K
                    k1_ontime = ix + find(rep(ix+1:end)==1,1);
                    k1_offtime = k1_ontime + find(diff(rep(k1_ontime:end))~=0,1);
                else
                    k1_ontime = ix + find(v(ix+1:end)==k1,1);
                    k1_offtime = k1_ontime + find(diff(v(k1_ontime:end))~=0,1);
                end
                
                if k2>K
                    k2_ontime = k1_offtime  + find(rep(k1_offtime:end)==1,1) - 1;
                else
                    k2_ontime = k1_offtime  + find(v(k1_offtime:end)==k2,1) - 1;
                end
                
                if k2_ontime == k1_offtime %k2 comes right after k1
                    intv=[intv; nan];
                    tmp_sum_occ = [tmp_sum_occ; 0];
                    tmp_sum_occ_bin = [tmp_sum_occ_bin; zeros(1,nbins)];
                elseif (~any(v(k1_offtime:k2_ontime)==k1) || ~any(rep(k1_offtime:k2_ontime)==1)) && ...% if k1 appears again before k2 does - abort
                        ~isempty(k2_ontime)  % or if k2 doesn't appear anymore - also abort
                    intv=[intv; k2_ontime-k1_offtime];
                    bin_size = floor((k2_ontime-k1_offtime)./nbins);
                    residual_bin = floor(((k2_ontime-k1_offtime) - floor((k2_ontime-k1_offtime)./nbins)*nbins)./(nbins-1)); % %binsize*nbins might nog equal the interval. Leave residual in between bins
                    
                    if k3>K
                        tmp_sum_occ = [tmp_sum_occ; sum(rep(k1_offtime:k2_ontime)==1)];
                        tmp_sum_occ_bin = [tmp_sum_occ_bin; sum(rep(k1_offtime:k1_offtime+bin_size-1)==1), sum(rep(k1_offtime+bin_size+residual_bin:k1_offtime+2*bin_size+residual_bin-1)==1), sum(rep(k2_ontime-bin_size:k2_ontime-1)==1)];
                    else
                        tmp_sum_occ = [tmp_sum_occ; sum(v(k1_offtime:k2_ontime)==k3)];
                        tmp_sum_occ_bin = [tmp_sum_occ_bin; sum(v(k1_offtime:k1_offtime+bin_size-1)==setdiff(states, [k1,k2])), sum(v(k1_offtime+bin_size+residual_bin:k1_offtime+2*bin_size+residual_bin-1)==setdiff(states, [k1,k2])), sum(v(k2_ontime-bin_size:k2_ontime-1)==setdiff(states, [k1,k2]))];
                    end
                end
                ix = k1_offtime;
                
                if isempty(k1_ontime) || isempty(k1_offtime) || isempty(k2_ontime)
                    w=1;
                end
            end
            
            sum_occ{k1,k2,iSes} = tmp_sum_occ;
            sum_occ_bin{k1,k2,iSes} = tmp_sum_occ_bin; %
            intervals{k1,k2,iSes} = intv;
            FO_interplay_ses(k1,k2,iSes) = mean(tmp_sum_occ(~isnan(intv))./intv(~isnan(intv)));%./mean(v==k1);
        end
    end
end
sum_occ=sum_occ(states,states,:);
sum_occ_bin=sum_occ_bin(states,states,:);
intervals=intervals(states,states,:);
FO_interplay=(FO_interplay_ses(states,states,1:2:end)+FO_interplay_ses(states,states,2:2:end))./2;

% Now seperate into percentiles
percentiles = 10:10:90;
nperc=length(percentiles);
for k1=1:3
    for k2=setdiff(1:3,k1)
        for iSes=1:nSes
            tmp1 = intervals{k1,k2,iSes};
            tmp2 = [intervals{k1,k2,iSes}; intervals{k2,k1,iSes}];tmp2(isnan(tmp2))=[];% take percentile over both i-j and j-i
            for k3=1:nperc
                ix = tmp1>=percentile(tmp2, percentiles(k3)-9) & tmp1<percentile(tmp2, percentiles(k3)+10);
                sum_occ_bin_perc(k1,k2,k3,:,iSes) = nanmean(sum_occ_bin{k1,k2,iSes}(ix,:));%nanmean(sum_occ_bin{k1,k2,iSes}(ix,:)./intervals{k1,k2,iSes}(ix));
                intervals_perc(k1,k2,k3,iSes) = percentile(tmp2, percentiles(k3))./250;
            end
        end
    end
end
sum_occ_bin_perc = nanmean(sum_occ_bin_perc,5); % mean over subjects

% and normalize by interval percentile (rows):
for k1=1:3
    for k2=setdiff(1:3,k1)
        perc_avg(k1,k2,:,:,:) = repmat(mean((sum_occ_bin_perc(k1,k2,:,:,:) + sum_occ_bin_perc(k2,k1,:,:,:))./2,4), [1,1,1,nbins,1]);
    end
end
sum_occ_bin_perc_norm = nanmean(sum_occ_bin_perc./perc_avg,5); % normalized mean over subjects
intervals_perc = nanmean(intervals_perc,4)*1000;



cmap=inferno;
fig=setup_figure([],1.5,2);
subplot(3,1,1),hold on, imagesc(1:3,1:nperc,squeeze(sum_occ_bin_perc(1,2,:,:))), axis xy, imagesc(4:6, 1:nperc, squeeze(sum_occ_bin_perc(2,1,:,:))), xlim([0.5, 6.5]), ylim([0.5 nperc+.5])
c=caxis; caxis([0,1].*max(abs(c))); ylabel('mean IT'), xlabel('Interval time bin')
colorbar, yticks(1:nperc), yticklabels(round(squeeze(intervals_perc(1,2,:)))), colormap(cmap), xticks([0.5 1:3 3.5 4:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','\bfPAN\rm', '4','5','6','\bfDMN\rm'}),
xtickangle(45), vline(3.5, '-k'), title('Replay Occurrence')

subplot(3,1,2),hold on, imagesc(1:3,1:nperc,squeeze(sum_occ_bin_perc(1,3,:,:))), axis xy, imagesc(4:6, 1:nperc, squeeze(sum_occ_bin_perc(3,1,:,:))), xlim([0.5, 6.5]), ylim([0.5 nperc+.5])
c=caxis; caxis([0,1].*max(abs(c))); ylabel('mean IT'), xlabel('Interval time bin')
colorbar, yticks(1:nperc), yticklabels(round(squeeze(intervals_perc(1,3,:)))), colormap(cmap), xticks([0.5 1:3 3.5 4:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','\bfReplay\rm', '4','5','6','\bfDMN\rm'}),
xtickangle(45), vline(3.5, '-k'), title('PAN Occurrence')

subplot(3,1,3),hold on, imagesc(1:3,1:nperc,squeeze(sum_occ_bin_perc(2,3,:,:))), axis xy, imagesc(4:6, 1:nperc, squeeze(sum_occ_bin_perc(3,2,:,:))), xlim([0.5, 6.5]), ylim([0.5 nperc+.5])
c=caxis; caxis([0,1].*max(abs(c))); ylabel('mean IT'), xlabel('Interval time bin')
colorbar, yticks(1:nperc), yticklabels(round(squeeze(intervals_perc(2,3,:)))), colormap(cmap), xticks([0.5 1:3 3.5 4:6, 6.5]),xticklabels({'\bfPAN\rm', '1','2','3','\bfReplay\rm', '4','5','6','\bfPAN\rm'}),
xtickangle(45), vline(3.5, '-k'), title('DMN Occurrence')
suptitle({'Binned and percentiled FO', ''})
save_figure([config.figdir, sprintf('2_FO_DMN_PAN_Replay_binned_percentiled_perc%d', perc)],[],false);


cmap = flipud(brewermap(256, 'RdBu'));
fig=setup_figure([],1.5,2);
subplot(3,1,1),hold on, imagesc(1:3,1:nperc,squeeze(sum_occ_bin_perc_norm(1,2,:,:))), axis xy, imagesc(4:6, 1:nperc, squeeze(sum_occ_bin_perc_norm(2,1,:,:))), xlim([0.5, 6.5]), ylim([0.5 nperc+.5])
c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); ylabel('mean IT'), xlabel('Interval time bin')
colorbar, yticks(1:nperc), yticklabels(round(squeeze(intervals_perc(1,2,:)))), colormap(cmap), xticks([0.5 1:3 3.5 4:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','\bfPAN\rm', '4','5','6','\bfDMN\rm'}),
xtickangle(45), vline(3.5, '-k'), title('Replay Occurrence')

subplot(3,1,2),hold on, imagesc(1:3,1:nperc,squeeze(sum_occ_bin_perc_norm(1,3,:,:))), axis xy, imagesc(4:6, 1:nperc, squeeze(sum_occ_bin_perc_norm(3,1,:,:))), xlim([0.5, 6.5]), ylim([0.5 nperc+.5])
c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); ylabel('mean IT'), xlabel('Interval time bin')
colorbar, yticks(1:nperc), yticklabels(round(squeeze(intervals_perc(1,3,:)))), colormap(cmap), xticks([0.5 1:3 3.5 4:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','\bfReplay\rm', '4','5','6','\bfDMN\rm'}),
xtickangle(45), vline(3.5, '-k'), title('PAN Occurrence')

subplot(3,1,3),hold on, imagesc(1:3,1:nperc,squeeze(sum_occ_bin_perc_norm(2,3,:,:))), axis xy, imagesc(4:6, 1:nperc, squeeze(sum_occ_bin_perc_norm(3,2,:,:))), xlim([0.5, 6.5]), ylim([0.5 nperc+.5])
c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); ylabel('mean IT'), xlabel('Interval time bin')
colorbar, yticks(1:nperc), yticklabels(round(squeeze(intervals_perc(2,3,:)))), colormap(cmap), xticks([0.5 1:3 3.5 4:6, 6.5]),xticklabels({'\bfPAN\rm', '1','2','3','\bfReplay\rm', '4','5','6','\bfPAN\rm'}),
xtickangle(45), vline(3.5, '-k'), title('DMN Occurrence')
suptitle({'Binned and percentiled FO (normalized by percentile FO)', ''})
save_figure([config.figdir, sprintf('2_FO_DMN_PAN_Replay_binned_percentiled_normalized_perc%d', perc)],[],false);

%%



% make bar plot for various time intervals
t = [0,1,50, 200, 2000, Inf]./1000; % in s
w_avg = zeros(length(states), length(states), length(t)-1, nSes);
avg = zeros(length(states), length(states), length(t)-1, nSes);
total = zeros(length(states), length(states), length(t)-1, nSes);
for k1=1:length(states)
    for k2=setdiff(1:length(states),k1)
        for it = 1:length(t)-1
            for iSes=1:nSes
                ix = intervals{k1,k2,iSes}./250>=t(it) & intervals{k1,k2,iSes}./250 < t(it+1);
                w_avg(k1,k2,it,iSes) = mean(sum_occ{k1,k2,iSes}(ix)./intervals{k1,k2,iSes}(ix)./250);
                avg(k1,k2,it,iSes) = mean(sum_occ{k1,k2,iSes}(ix));
                total(k1,k2,it,iSes) = sum(sum_occ{k1,k2,iSes}(ix));
            end
        end
    end
end
w_avg = (w_avg(:,:,:,1:2:end) + w_avg(:,:,:,2:2:end))./2;
W=nanmean(w_avg,4);
avg = (avg(:,:,:,1:2:end) + avg(:,:,:,2:2:end))./2;
Av=nanmean(avg,4);
total = (total(:,:,:,1:2:end) + total(:,:,:,2:2:end));
TO = nanmean(total,4);

% Boxplot?
% histogram (re time)?
% where in the interval - interval bins?


fig=setup_figure([],1.5,2);
subplot(3,1,1), hold on
boxplot_with_scatter([squeeze(w_avg(1,2,:,:));-ones(3,43);squeeze(w_avg(2,1,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)]),
ylim([0,1.05*max([squash(w_avg(1,2,:,:)); squash(w_avg(2,1,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfDMN\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfPAN\rm','0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfDMN\rm'})
xtickangle(45), ylabel('FO'), title('Replay occurrence'),
subplot(3,1,2)
boxplot_with_scatter([squeeze(w_avg(1,3,:,:));-ones(3,43);squeeze(w_avg(3,1,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)]),
ylim([0,1.05*max([squash(w_avg(1,3,:,:)); squash(w_avg(3,1,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfDMN\rm',  '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfReplay\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfDMN\rm'})
xtickangle(45), ylabel('FO'), title('PAN occurrence')
subplot(3,1,3)
boxplot_with_scatter([squeeze(w_avg(2,3,:,:));-ones(3,43);squeeze(w_avg(3,2,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)])
ylim([0,1.05*max([squash(w_avg(2,3,:,:)); squash(w_avg(3,2,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfPAN\rm',  '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfReplay\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfPAN\rm'})
xtickangle(45), ylabel('FO'), title('DMN occurrence')
suptitle('sum of occurrence divided by interval length')
save_figure([config.figdir, sprintf('2_FO_DMN_PAN_Replay_1_perc%d', perc)],[],false);

fig=setup_figure([],1.5,2);
subplot(3,1,1), hold on
boxplot_with_scatter([squeeze(avg(1,2,:,:));-ones(3,43);squeeze(avg(2,1,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)]),
ylim([0,1.05*max([squash(avg(1,2,:,:)); squash(avg(2,1,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfDMN\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfPAN\rm','0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfDMN\rm'})
xtickangle(45), ylabel('FO'), title('Replay occurrence'),
subplot(3,1,2)
boxplot_with_scatter([squeeze(avg(1,3,:,:));-ones(3,43);squeeze(avg(3,1,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)]),
ylim([0,1.05*max([squash(avg(1,3,:,:)); squash(avg(3,1,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfDMN\rm',  '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfReplay\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfDMN\rm'})
xtickangle(45), ylabel('FO'), title('PAN occurrence')
subplot(3,1,3)
boxplot_with_scatter([squeeze(avg(2,3,:,:));-ones(3,43);squeeze(avg(3,2,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)])
ylim([0,1.05*max([squash(avg(2,3,:,:)); squash(avg(3,2,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfPAN\rm',  '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfReplay\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfPAN\rm'})
xtickangle(45), ylabel('FO'), title('DMN occurrence')
suptitle('sum of occurrence divided by number of intervals')
save_figure([config.figdir, sprintf('2_FO_DMN_PAN_Replay_2_perc%d', perc)],[],false);

fig=setup_figure([],1.5,2);
subplot(3,1,1), hold on
boxplot_with_scatter([squeeze(total(1,2,:,:));-ones(3,43);squeeze(total(2,1,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)]),
ylim([0,1.05*max([squash(total(1,2,:,:)); squash(total(2,1,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfDMN\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfPAN\rm','0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfDMN\rm'})
xtickangle(45), ylabel('FO'), title('Replay occurrence'),
subplot(3,1,2)
boxplot_with_scatter([squeeze(total(1,3,:,:));-ones(3,43);squeeze(total(3,1,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)]),
ylim([0,1.05*max([squash(total(1,3,:,:)); squash(total(3,1,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfDMN\rm',  '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfReplay\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfDMN\rm'})
xtickangle(45), ylabel('FO'), title('PAN occurrence')
subplot(3,1,3)
boxplot_with_scatter([squeeze(total(2,3,:,:));-ones(3,43);squeeze(total(3,2,:,:))]', [repmat({[0 0.4470 0.7410]}, 8,1); repmat({[0.8500 0.3250 0.0980]}, 5,1)])
ylim([0,1.05*max([squash(total(2,3,:,:)); squash(total(3,2,:,:))])]),xlim([0,14]), h=gca;h.TickLabelInterpreter = 'tex';
xticks([0:5, 7, 9:14]), xticklabels({'\bfPAN\rm',  '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfReplay\rm', '0ms', '1-50ms', '50-200ms', '.2-2s', '2+s', '\bfPAN\rm'})
xtickangle(45), ylabel('FO'), title('DMN occurrence')
suptitle('sum of occurrence')
save_figure([config.figdir, sprintf('2_FO_DMN_PAN_Replay_3_perc%d', perc)],[],false);

%% Second level HMM
% create a sliding window DMN FO
w = 2*ceil(0.2*250/2); % in samples
for iSes=1:nSes
    X_poiss_DMN{iSes} = zeros(length(vpath{iSes}),1);
    for i = w/2+1:length(vpath{iSes})-w/2
        X_poiss_DMN{iSes}(i,1) = sum(vpath{iSes}(i-w/2:i+w/2-1)==1);
    end
    %     X_poiss_DMN{iSes}=zscore(X_poiss_DMN{1});
end

% Run 2nd level HMM
n_runs=5;
for i_run = 1:n_runs
    options = [];
    options.K = 2;
    options.distribution = 'poisson';
    %options.Pstructure = eye(options.K) + diag(ones(1,options.K-1),1);
    %ptions.Pstructure(options.K,1) = 1;
    options.initrep = 4; % this the number of parallel cores
    options.useParallel = false;
    options.DirichletDiag=1;
    
    
    Sjs_init = randperm(config.nSj,2);
    X_init = []; T_init = [];
    for iSj = Sjs_init
        X_init = [X_init; X_poiss_DMN{iSj}];
        T_init = [T_init; T{iSj}];
    end
    [hmminit,gammainit] = hmmmar(X_init,T_init,options);
    
    options = rmfield(options,'initrep');
    options.hmm = hmminit;
    
    options.decodeGamma = false;
    options.standardise = false;
    [hmmtemp,Gammatemp,~,~,~,~,fehist] = hmmmar(cat(1,X_poiss_DMN{:}),cat(1,T{:}),options);
    
    if i_run==1 || fehist(end)<lowestfe
        hmmPoiss = hmmtemp;
        GammaPoiss = Gammatemp;
        lowestfe = fehist(end);
    end
    % record a few stats to get a sense of how deviant these really are:
    feall(i_run,1) = fehist(1);
    feall(i_run,2) = fehist(end);
    gamsum(i_run,:) = mean(Gammatemp);
end

save([config.resultsdir, sprintf('tinda_replay_perc%d_dmn_burst', perc)], 'options', 'feall', 'w', 'X_poiss_DMN')

% state 1 is DMN off
% state 2 is DMN on
options.dropstates=false;
for iSes=1:nSes
    LT = getStateLifeTimes(vpath_DMN{iSes},T{iSes},options);
    LTmedian(iSes,:) = cellfun(@median,LT);
    LTmu(iSes,:) = cellfun(@mean,LT);
    FracOcc(iSes,:) = getFractionalOccupancy(vpath_DMN{iSes},sum(T{iSes}),options);
    IT = getStateIntervalTimes(vpath_DMN{iSes},T{iSes},options);
    ITmedian(iSes,:) = cellfun(@median,IT);
    ITmu(iSes,:) = cellfun(@mean,IT);
    
end
LTmedian = (LTmedian(1:2:end,:) + LTmedian(2:2:end,:))./2;
ITmedian = (ITmedian(1:2:end,:) + ITmedian(2:2:end,:))./2;
LTmu = (LTmu(1:2:end,:) + LTmu(2:2:end,:))./2;
ITmu = (ITmu(1:2:end,:) + ITmu(2:2:end,:))./2;
FracOcc = (FracOcc(1:2:end,:) + FracOcc(2:2:end,:))./2;

dmn_burst = [];
dmn_burst.FO = FracOcc;
dmn_burst.LT = LT;
dmn_burst.LTmedian = LTmedian;
dmn_burst.LTmu = LTmu;
dmn_burst.IT = IT;
dmn_burst.ITmedian = ITmedian;
dmn_burst.ITmu = ITmu;


figure;
subplot(1,3,1), distributionPlot(FracOcc), title('Fractional Occupancy')
subplot(1,3,2), distributionPlot(ITmedian), title('Interval Times')
subplot(1,3,3), distributionPlot(LTmedian), title('Life Times')
save_figure([config.figdir, '2_DMN_burst_stats'],[],false);

% TINDAsum_occ_bin_perc
k=2;
percentiles = [0:20:100];
for i=1:length(percentiles)-1
    [FO_p_masked(:,:,:,:,i),~,t_intervals_p{i},~] = computeLongTermAsymmetry(vpath_DMN,T,2, [percentiles(i), percentiles(i+1)],[],topidx,false,6);
end
FO_p_masked = (FO_p_masked(:,:,:,1:2:end,:) + FO_p_masked(:,:,:,2:2:end,:))/2;
sum_occ_bin_perc = squeeze(mean(FO_p_masked,4));
t_int = squeeze(mean(cellfun(@mean, cat(3,t_intervals_p{:})),1));
fig=setup_figure([],2,1);
cmap = inferno;
bins = 1:size(FO_p_masked,3);
subplot(2,2,1); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc(1,3,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title({'Replay when DMN is ON', '(DMN OFF interval)'}), xlabel('DMN bursting'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(1,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfDMN silence \rm', '1','2','3','4','5','6','\bfDMN silence\rm'}), xtickangle(45), vline(3.5, '-k')
subplot(2,2,2); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc(3,2,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title({'DMN burst in', 'replay intervals'}), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(3,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')
subplot(2,2,3); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc(2,3,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title({'Replay when DMN is OFF', '(DMN ON interval)'}), xlabel('DMN not bursting'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(2,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfDMN burst\rm', '1','2','3','4','5','6','\bfDMN burst\rm'}), xtickangle(45), vline(3.5, '-k')
subplot(2,2,4); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc(3,1,:,:))'); c=caxis; caxis([0,1].*max(abs(c))); axis xy, title({'DMN silence in', 'replay intervals'}), xlabel('Interval time bin'), ylabel('mean IT'), colorbar, yticks(1:9), yticklabels(round(t_int(3,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')
suptitle({'Binned and percentiled FO',''})
save_figure([config.figdir, sprintf('2_DMN_burst_FO_binned_percentiled_perc%d', perc)],[],false);

% Normalize FO by row
sum_occ_bin_perc_norm = squeeze(nanmean(FO_p_masked./mean(FO_p_masked,3),4));
cmap = flipud(brewermap(256, 'RdBu'));
fig=setup_figure([],2,1);
subplot(2,2,1); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc_norm(1,3,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title({'Replay when DMN is ON', '(DMN OFF interval)'}), xlabel('DMN bursting'), ylabel('mean IT (ms)'), colorbar, yticks(1:9), yticklabels(round(t_int(1,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfnoDMN\rm', '1','2','3','4','5','6','\bfnoDMN\rm'}), xtickangle(45), vline(3.5, '-k')
subplot(2,2,2); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc_norm(3,2,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title({'DMN burst in', 'replay intervals'}), xlabel('Interval time bin'), ylabel('mean IT (ms)'), colorbar, yticks(1:9), yticklabels(round(t_int(3,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')
subplot(2,2,3); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc_norm(2,3,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title({'Replay when DMN is OFF', '(DMN ON interval)'}), xlabel('DMN not bursting'), ylabel('mean IT (ms)'), colorbar, yticks(1:9), yticklabels(round(t_int(2,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfDMN\rm', '1','2','3','4','5','6','\bfDMN\rm'}), xtickangle(45), vline(3.5, '-k')
subplot(2,2,4); imagesc(bins, 1:5, squeeze(sum_occ_bin_perc_norm(3,1,:,:))'); c=caxis; c=max(abs(c-1)); caxis([1-c 1+c]); axis xy, title({'DMN silence in', 'replay intervals'}), xlabel('Interval time bin'), ylabel('mean IT (ms)'), colorbar, yticks(1:9), yticklabels(round(t_int(3,:))), colormap(cmap), xticks([0.5 1:6, 6.5]),xticklabels({'\bfReplay\rm', '1','2','3','4','5','6','\bfReplay\rm'}), xtickangle(45), vline(3.5, '-k')

pause(1); suptitle({'Binned and percentiled FO (normalized by percentile FO)'})
save_figure([config.figdir, sprintf('2_DMN_burst_FO_binned_percentiled_normalized_perc%d', perc)],[],false);

dmn_burst.FO_intervals = FO_p_masked;
dmn_burst.t_intervals = t_intervals_p;
dmn_burst.sum_occ_bin_perc = sum_occ_bin_perc;
dmn_burst.sum_occ_bin_perc_norm = sum_occ_bin_perc_norm;
save([config.resultsdir, sprintf('tinda_replay_perc%d_dmn_burst', perc)], 'dmn_burst', '-append')



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