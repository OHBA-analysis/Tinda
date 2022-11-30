whichstudy = 2;
if whichstudy==1
  MFStartup;
else
  MFStartup_studyII;
  load([is.rootBehav, 'Exptrial.mat'])
  is.dirname='StrLearn_MEGexp/';
end
if ~isfolder(is.AnalysisPath),   mkdir(is.AnalysisPath), end
is.usePrecomputed = true;
is.iRun = 2;
is.topPercentile=1;
is.lambda = 5; % hardcoded best lambda for all subjects
is.whenInMS = is.tss*is.msPerSample;
is.TOI = 200; % 200ms post stimulus.
is.sessionsOfInterest = [1,2,3,4,8]; % rst, lci, lci, lci, rst
is.nStim = 8;
is.ntrlAll=120; % per session (includes upside down images)
is.ntrlValid=96; % per session
is.nSes = 3;
is.nChan = 273;
is.nTime = 401;
is.iRun = 2;
is.nPerm = 1000;
is.K = 12;
is.plot=false;
fsample=250;
is.t_window=fsample/2;


%% Get the Study specific HMM results
fname = fullfile(is.studydir, 'Neuron2020Analysis/', sprintf('Study%d',whichstudy), "hmm_1to45hz", "hmm5usingtemplate_parc_giles_symmetric__pcdim80_voxelwise_embed13_K12_big1_dyn_modelhmm.mat");
replayflag=1;
if replayflag% use the replay ordering
  load(fname, 'hmm', 'new_state_ordering')
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

  G(iSes,:,:) = Gam_long;

  temp1 = zeros(length(goodsamples{iSes}),1);
  temp1(goodsamples{iSes}) = hmm.statepath(R(iSes,1):R(iSes,2));
  temp1 = temp1(triggerpoints{iSes},:);
  vpath{iSes} = temp1;
  T{iSes} = [length(vpath{iSes})];
end

% vpath=vpath(is.iRun:2:end);
% hmmT = hmmT(is.iRun:2:end);
%% Get the sequence from the canonical RS
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
bestseq = bestsequencemetrics{2};
% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
    disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));


%% Run TINDA on replay events
% Load the replay timing info
load([is.resultsdir, 'replayidx/replay_idx.mat'])
for k=1:K
  for i=1:44
    vpath_replay{i} = zeros(size(vpath{i}));
    vpath_replay{i}(vpath{i}==k) = 1;
    vpath_replay{i}(topidx{i}) = 2;
  end
  [FO_intervals{k},FO_pvals{k},t_intervals{k},FO_stat{k}] = computeLongTermAsymmetry(vpath_replay,T,2);
  %redo the stats because we first want to average over sessions within a
  %subject
  FO_intervals{k} = (FO_intervals{k}(:,:,:,1:2:end) + FO_intervals{k}(:,:,:,2:2:end))/2;
  [FO_pvals{k}, FO_stat{k}] = FO_permutation_test(FO_intervals{k}(:,:,:,:), 2, 22);
  
  % only select the asymmetry of replay events (1st entry)
  FO_intervals{k} = FO_intervals{k}(1,2,:,:);
  FO_pvals{k} = FO_pvals{k}(1,2);
  FO_stat{k}.prob = FO_stat{k}.prob(1);
  FO_stat{k}.cirange = FO_stat{k}.cirange(1);
  FO_stat{k}.mask = FO_stat{k}.mask(1);
  FO_stat{k}.stat = FO_stat{k}.stat(1);
  FO_stat{k}.ref = FO_stat{k}.ref(1);
  FO_stat{k}.time = 1;
  FO_assym{k} = squeeze(FO_intervals{k}(:,:,1,:)) - squeeze(FO_intervals{k}(:,:,2,:));
end


% plot distribution of absolute asymetry
q  = abs(mean(cat(2,FO_assym{:})));
CL = [0 max((q))];
cyclicalstate_distributionplot(bestseq, (q)', CL)



%% Run TINDA on the replay events
[CL(1), CL(2)] = bounds(squash(q(f>13 & f<40,:)));tmp = max(abs(CL-1));CL = [1-tmp, 1+tmp]; cyclicalstate_distributionplot(bestseq, mean(q(f<40 & f>13, :),1), CL, hotcold(6:5:end-5,:))








%% First check for HMM cycles in the resting state data (irrespective of replay)
[FO,pvals,t_intervals] = computeLongTermAsymmetry(vpath,T,K);
FO = (FO(:,:,:,1:2:end) + FO(:,:,:,2:2:end))/2;
[pvals, stat] = FO_permutation_test(FO(:,:,:,:), K, 22);

bonf_ncomparisons = 2*(K.^2-K);
sigpoints = pvals<(0.05);

mean_direction = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
mean_assym = squeeze(mean((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3),4));

hmm_1stlevel.FO_intervals = FO;
hmm_1stlevel.FO_assym = squeeze((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3));

% This is one way to define the strength of the circularity.
rotational_momentum = imag(sum(sum(angleplot.*hmm_1stlevel.FO_assym)));
hmm_1stlevel.rotational_momentum = squeeze(rotational_momentum);


% plot as circular diagram:
cyclicalstateplot(bestseq,mean_direction,sigpoints);
gcf;
print([is.resultsdir, 'Figures/', 'CyclicalPattern_uncorrected'],'-dpng');


[circularity, circle_pval, ~, ~, fig] = geometric_circularity(mean_direction(bestseq, bestseq), sigpoints(bestseq, bestseq));
hmm_1stlevel.circularity = circularity;
hmm_1stlevel.pval = circle_pval;
tmp = hmm_1stlevel.FO_assym; tmp(isnan(tmp))=0;
for k=1:size(tmp,3)
  [hmm_1stlevel.circularity_subject(k,1), hmm_1stlevel.circularity_pval_subject(k,1), ~, ~, ~] = geometric_circularity(tmp(bestseq, bestseq,k), sigpoints(bestseq, bestseq),[],[],0);
end
gcf;
print([is.resultsdir, 'Figures/', 'Circularity_uncorrected'],'-dpng');


%% Find the Replay locked v-path
% shape the vpath in a way in which we get trials concatenated. Also make a
% T matrix.
fname = [is.resultsdir, 'replayidx/replay_idx'];
if is.useMask
  fname = [fname, sprintf('_mask%d')];
end
tmp=load([fname, '.mat']);
topidx=tmp.topidx;
mult=4; % 2 seconds on each side
window = is.t_window*mult;
time = (-window:window)/fsample;
for iSj=1:length(vpath)
  v = vpath{iSj};
  ix = topidx{iSj};
  v_tl = zeros(2*window+1, length(ix));
  rem=[];
  for i=1:length(ix)
    try
      v_tl(:,i) = v(ix(i)-window:ix(i)+window);
    catch
      rem = [rem, i];
    end
  end
  v_tl(:, rem) = [];
  vpath_replay{iSj} = reshape(v_tl, [], 1);
  
  T_replay{iSj} = (2*window+1)*ones(size(v_tl,2),1);
end

hmm_1stlevel.Replay.window = 2*window+1;
hmm_1stlevel.Replay.time = time;

%% Compute the replay locked assymetry (TINDA) 
% CAN WE USE THE NEXT LINE?
[FO_replay,pvals_replay,~,stat_replay, FO_replay_residual, pvals_replay_residual, stat_residual] = computeLongTermAsymmetry(vpath,T,hmm.K,[],1, topidx{iSj});

[FO_replay,pvals_replay,~,stat_replay, FO_replay_residual, pvals_replay_residual, stat_residual] = computeLongTermAsymmetry(vpath_replay,T_replay,hmm.K,[],1);
FO_replay = (FO_replay(:,:,:,1:2:end)+FO_replay(:,:,:,2:2:end))/2;
FO_replay_residual = (FO_replay_residual(:,:,:,1:2:end)+FO_replay_residual(:,:,:,2:2:end))/2;

[pvals_replay, stat_replay] = FO_permutation_test(FO_replay(:,:,:,:), K, 22);

[pvals_replay_residual, stat_residual] = FO_permutation_test(FO_replay_residual(:,:,:,:), K, 22);


% first without correction for the evoked response
sigpoints = pvals_replay<(0.05);

mean_direction = squeeze(nanmean(FO_replay(:,:,1,:)-FO_replay(:,:,2,:),4));
mean_assym = squeeze(mean((FO_replay(:,:,1,:)-FO_replay(:,:,2,:))./mean(FO_replay,3),4));

hmm_1stlevel.Replay.raw.FO_intervals = FO_replay;
hmm_1stlevel.Replay.raw.FO_assym = squeeze((FO_replay(:,:,1,:)-FO_replay(:,:,2,:))./mean(FO_replay,3));
hmm_1stlevel.Replay.raw.FO_assym_pval = pvals_replay;

rotational_momentum = imag(nansum(nansum(angleplot.*hmm_1stlevel.Replay.raw.FO_assym)));
hmm_1stlevel.Replay.raw.rotational_momentum = squeeze(rotational_momentum);

cyclicalstateplot(bestseq,mean_direction,sigpoints);
gcf;
print([is.resultsdir, 'Figures/', 'CyclicalPattern_ReplayRaw_uncorrected'],'-dpng');
[circularity_raw, circle_pval_raw, ~, ~, fig] = geometric_circularity(mean_direction(bestseq, bestseq), sigpoints(bestseq, bestseq));
hmm_1stlevel.Replay.raw.circularity = circularity_raw;
hmm_1stlevel.Replay.raw.circularity_pval = circle_pval_raw;
gcf;
print([is.resultsdir, 'Figures/', 'Circularity_ReplayRaw_uncorrected'],'-dpng');

tmp = hmm_1stlevel.Replay.raw.FO_assym; tmp(isnan(tmp))=0;
for k=1:size(tmp,3)
  [hmm_1stlevel.Replay.raw.circularity_subject(k,1), hmm_1stlevel.Replay.raw.circularity_pval_subject(k,1), ~, ~, ~] = geometric_circularity(tmp(bestseq, bestseq,k), sigpoints(bestseq, bestseq),[],[],0);
end

% Now correct for the evoked response
sigpoints_residual = pvals_replay_residual<(0.05);

mean_direction_residual = squeeze(nanmean(FO_replay_residual(:,:,1,:)-FO_replay_residual(:,:,2,:),4));
mean_assym_residual = squeeze(mean((FO_replay_residual(:,:,1,:)-FO_replay_residual(:,:,2,:))./mean(FO_replay_residual,3),4));

hmm_1stlevel.Replay.FO_intervals = FO_replay_residual;
hmm_1stlevel.Replay.FO_assym = squeeze((FO_replay_residual(:,:,1,:)-FO_replay_residual(:,:,2,:))./mean(FO_replay_residual,3));

rotational_momentum_residual = imag(nansum(nansum(angleplot.*hmm_1stlevel.Replay.FO_assym)));
hmm_1stlevel.Replay.rotational_momentum = squeeze(rotational_momentum_residual);

cyclicalstateplot(bestseq,mean_direction_residual,sigpoints_residual);
gcf;
print([is.resultsdir, 'Figures/', 'CyclicalPattern_ReplayResidual_uncorrected'],'-dpng');
[circularity_residual, circle_pval_residual, ~, ~, fig] = geometric_circularity(bestseq, mean_direction_residual, sigpoints_residual);
hmm_1stlevel.Replay.raw.circularity = circularity_residual;
hmm_1stlevel.Replay.raw.circularity_pval = circle_pval_residual;
gcf;
print([is.resultsdir, 'Figures/', 'Circularity_ReplayResidual_uncorrected'],'-dpng');
tmp = hmm_1stlevel.Replay.FO_assym; tmp(isnan(tmp))=0;
for k=1:size(tmp,3)
  [hmm_1stlevel.Replay.circularity_subject(k,1), hmm_1stlevel.Replay.circularity_pval_subject(k,1), ~, ~, ~] = geometric_circularity(tmp(bestseq, bestseq,k), sigpoints(bestseq, bestseq),[],[],0);
end

save([is.resultsdir, 'tinda/HMMsummarymetrics.mat'], 'hmm_1stlevel', '-append')





%% Vpath Circle analysis
triggerpoints=triggerpoints(2:2:end);
goodsamples=goodsamples(2:2:end);
for i=1:22
  temp1 = getCircleVpath(vpath{i}, bestseq);
  vpathcircle{i} = double(temp1);
  trialonset = zeros(75000,1);
  trialonset(topidx{i})=1;
  [tmp, vpathcircle_freq, vpathcircle_time] = fourieranalysis_circleVpath(vpathcircle{i}, trialonset);
  circlepow_evoked(i,:,:) = squeeze(nanmean(abs(tmp).^2));
  circlespctrm_evoked(i,:,:) = squeeze(nanmean(tmp));
end

for i=1:22
for k=1:length(topidx{i})
Q{i}(k,:) = vpathcircle{i}(topidx{i}(k)-250:topidx{i}(k)+250);
end
end
for k=1:22
Q2(k,:) = nanmean(Q{k});
end


%% Within a window, sum over the positions in the unit circle and do spectral analysis on this
% get replay scores
% observed data

is.doPermutation=false;
if whichstudy==1
  tmp=load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
  replayScores = tmp.ToRall;
else
  StimAll = load([is.AnalysisPath, 'classifiers/Sequence_by_Training_4Cell/StimAll.mat']);
  sfAll = StimAll.sfAll;
  sbAll = StimAll.sbAll;
  [~, ~, Mreplay, ~] = find_bestLambda(is, sfAll, sbAll);
  ToRall = compute_replayTimecourse(is, [], Mreplay, is.iRun);
  replayScores = ToRall;
end

cnt=1;
for iSes=is.iRun:2:nSes 
    temp1 = zeros(length(goodsamples{iSes}),1);
    temp1(goodsamples{iSes}) = hmm.statepath(R(iSes,1):R(iSes,2));
    temp1 = temp1(triggerpoints{iSes},:);
    vpath{cnt} = temp1;
    cnt=cnt+1;
end


% get top percentile probabilities
topperc = replayScores > prctile(replayScores(:),100-is.topPercentile);
topperc = [zeros(size(topperc,1),1),diff(topperc,[],2)]==1;% eliminate adjacent points
topperc(:,1:is.t_window)=0;topperc(:,end-is.t_window:end)=0;% eliminate border points:



fs=250;
W = 125; %set arbitrary half second window for averaging
Noverlap = 1; % number of points to overlap over successive windows
t_last = 0;

% get the circle plot positions
optimalseqfile = [is.matsdir,'tinda/', 'bestseq.mat'];
load(optimalseqfile);
bestseq = bestsequencemetrics{2};

% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
  disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end

for iSj=1:22
 if 0
   GoodChannel = find_goodChannel(iSj, is); % get good channels
   gnf = train_classifier(iSj,is,Exptrial, GoodChannel); % get decoding models
   Rreds = apply_toRest(iSj, is, gnf, GoodChannel, is.iRun); % get all resting state decoding results
   RredsAll = Rreds{1,37,1,is.lambda}; % select the decoding results of the correct model
   replayScores =  1./(1+exp(RredsAll'));
   topperc = replayScores > prctile(replayScores(:),100-is.topPercentile);
   topperc = [zeros(size(topperc,1),1),diff(topperc,[],2)]==1;% eliminate adjacent points
   topperc(:,1:is.t_window)=0;topperc(:,end-is.t_window:end)=0;% eliminate border points:
   topidx = upsample_replayTimes(topperc, zeros(1,1,75000), is.t_window);
 else

  topidx = upsample_replayTimes(topperc(iSj,:), zeros(1,1,75000), is.t_window);
 end
  vpcircle = vpath{iSj};
  for k=1:12
    vpcircle(vpcircle==k) = disttoplot_manual(k);
  end
  
  % smooth vpcircle with half second window
  vpcircle=smoothdata(vpcircle,'movmean',W);

  % get replay locked phase information
  vpcircle_timelocked{iSj} = zeros(length(topidx), Fs+1);
  for i=1:length(topidx)
    vpcircle_timelocked{iSj}(i,:) = vpcircle(topidx(i)-W:topidx(i)+W);
  end
  vpcircle_timelocked_group(iSj,:) = mean(vpcircle_timelocked{iSj});
end


figure; for k=1:22
subplot(3,8,k), plot(t,angle(mean(vpcircle_timelocked{k}))), ylim([-pi, pi]), title(sprintf('subj %d', k)), xlabel('time (s)')
end
suptitle('cyclical state phase timelocked to replay')
figure; plot(t, angle(mean(vpcircle_timelocked_group)))
title('cyclical state phase timelocked to replay - GROUP AVERAGE'), ylim([-pi,pi]), xlabel('time (s)')


