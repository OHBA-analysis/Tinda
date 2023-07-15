% script called to load viterbi paths inferred and hmm objects and run
% post-hoc sequence analysis:
if ~exist('whichstudy','var')
    whichstudy = 1; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = colorscheme(whichstudy);
hotcold = cmap_hotcold;
clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};

% include option for simulations here - results to be saved in adjacent
% folder:
simtests = false;
if simtests
    config.figdir = strrep(config.figdir,['Study',int2str(whichstudy)],['Study',int2str(whichstudy),'_simtest']);
    mkdir(config.figdir);
end
use_WB_nnmf=true; % whether or not to use the wide band NNMF Diego's Nature Comms to select power and coherence (alternative is selecting 1-30 Hz)
useMT = true; % use the MT results instead of Cam's wavelet approach

%% Load HMM results
% first load data and plot basic temporal statistics:
clear new_state_ordering
if isempty(config.reordering_states)
    new_state_ordering=1:hmm.K;
else
    fname_stateorder = [config.reordering_states, '_state_ordering'];
    if useMT
        if use_WB_nnmf
            fname_stateorder = [fname_stateorder, '_MT_nnmf'];
        else
            fname_stateorder = [fname_stateorder, '_MT'];
        end
        load([config.resultsdir, fname_stateorder])
    end
end

if whichstudy>=6
    hmm.K=12;
else
    temp = load(fullfile(config.hmmfolder,config.hmmfilename));

    hmm = temp.hmm;
    if ~isfield(hmm,'gamma') && whichstudy<4
        hmm.gamma = temp.Gamma;
        hmm.statepath = temp.vpath;
    end
    if whichstudy<3
        %load(config.prepdatafile,'hmmT','subj_inds');
        if strcmp(config.reordering_states, 'replay')
            new_state_ordering = temp.new_state_ordering;
        end

        hmm = hmm_permutestates(hmm, new_state_ordering);
        for i=1:config.nSj
            hmmT{i} = sum(hmm.subj_inds==i);
        end
    elseif whichstudy==3
        hmm = hmm_permutestates(hmm, new_state_ordering);
        hmmT = temp.T_all;
        hmm.subj_inds = zeros(size(hmm.statepath));
        t_offset = 0;
        for i=1:length(hmmT)
            t_length = sum(hmmT{i}) - length(hmmT{i})*(length(hmm.train.embeddedlags)-1);
            hmm.subj_inds(t_offset + [1:t_length]) = i;
            t_offset = t_offset + t_length;
        end
        hmm.subj_inds = ceil(hmm.subj_inds/3); % account for multiple runs per subj
        hmmTold = reshape(hmmT,3,config.nSj);
        hmmTsubj = cell(config.nSj,1);
        for i=1:config.nSj
            hmmTsubj{i,:} = [hmmTold{1,i},hmmTold{2,i},hmmTold{3,i}]-(length(hmm.train.embeddedlags)-1);
        end
        hmmT = hmmTsubj;
        clear hmmTsubj hmmTold;
    elseif whichstudy==4
        hmm = hmm_permutestates(hmm, new_state_ordering);
        hmmT = temp.T_all;
        % correct for embedded lags:
        for i=1:length(hmmT)
            hmmT{i} = hmmT{i} - (length(hmm.train.embeddedlags)-1);
        end
        load(config.matfilelist);
        if ~isfolder(mat_files_orth{1})
            for k=1:length(mat_files_orth)
                tmp1 = fileparts(config.matfilelist);
                [~, tmp2] = fileparts(mat_files_orth{k});
                mat_files_orth{k} = fullfile(tmp1, tmp2);
            end
        end
    end
    clear temp vpath;
end
%% First level HMM statistics
if whichstudy<4
    Gamma = hmm.gamma;
end
K = hmm.K;
FO = zeros(K,K,2,config.nSj);
opts = [];
opts.K = 12;
opts.Fs = config.sample_rate;
opts.dropstates=0;
if whichstudy>=6
    q=load([config.resultsdir, 'vpath']);
    vpath=cellfun(@(x) x+1, cellfun(@double, cellfun(@transpose, q.vpath, 'UniformOutput', false), 'UniformOutput',false), 'UniformOutput', false);
    for k=1:numel(vpath)
        k
        tmp = vpath{k}+100;
        for k2=1:K
            tmp(tmp==new_state_ordering(k2)+100) = k2;
        end
        vpath{k} = tmp;
        if ~isfield(q, 'T') || length(q.T)~=length(q.vpath)
            q.T(k) = length(vpath{k});
        end
    end
    hmmT = num2cell(q.T);
    clear q
    if whichstudy==7 || whichstudy==8
        hmmT_ses = hmmT;
        vpath_ses = vpath;
        tmp = reshape(vpath, [config.nSes, config.nSj]);
        clear vpath hmmT
        for k=1:size(tmp,2)
            vpath{k} = cat(1,tmp{:,k});
            hmmT{k} = length(vpath{k});
        end
    end
end
for subnum=1:config.nSj
    if whichstudy<6
        fprintf(['\nProcessing subj ',int2str(subnum)]);
        if whichstudy~=4
            vpath{subnum} = hmm.statepath(hmm.subj_inds==subnum);
        else
            temp = load(mat_files_orth{subnum},'vpath');
            vpath_old = temp.vpath;
            vpath{subnum} = zeros(size(vpath_old));
            for k=1:K
                vpath{subnum}(vpath_old==new_state_ordering(k)) = k;
            end
        end
        
    end
    if simtests
        try % note unusual syntax here is just to catch very rare precision errors in numeric simulation
            vpath{subnum} = simulateVpath(vpath{subnum},hmmT{subnum},hmm.K);
        catch
            vpath{subnum} = simulateVpath(vpath{subnum},hmmT{subnum},hmm.K);
        end
    end
    LT = getStateLifeTimes(vpath{subnum},hmmT{subnum},opts);
    LTmerged(subnum,:) = cellfun(@mean,LT);
    LTmedian(subnum,:) = cellfun(@median,LT);
    LTstd(subnum,:) = cellfun(@std,LT);
    FracOcc(subnum,:) = getFractionalOccupancy(vpath{subnum},sum(hmmT{subnum}),opts);
        IT = getStateIntervalTimes(vpath{subnum},hmmT{subnum},opts);
    ITmerged(subnum,:) = cellfun(@mean,IT);
    ITmedian(subnum,:) = cellfun(@median,IT);
    ITstd(subnum,:) = cellfun(@std,IT);
    if whichstudy==7 || whichstudy==8
      for iSes=1:config.nSes
        FracOcc_ses(subnum,iSes,:) = getFractionalOccupancy(vpath_ses{(subnum-1)*config.nSes+iSes},sum(hmmT_ses{(subnum-1)*config.nSes+iSes}),opts);
        LT_ses = getStateLifeTimes(vpath_ses{(subnum-1)*config.nSes+iSes},hmmT_ses{(subnum-1)*config.nSes+iSes},opts);
        LTmerged_ses(subnum,iSes,:) = cellfun(@mean,LT_ses);
        LTmedian_ses(subnum,iSes,:) = cellfun(@median,LT_ses);
        LTstd_ses(subnum,iSes,:) = cellfun(@std,LT_ses);
        IT_ses = getStateIntervalTimes(vpath_ses{(subnum-1)*config.nSes+iSes},hmmT_ses{(subnum-1)*config.nSes+iSes},opts);
        ITmerged_ses(subnum,iSes,:) = cellfun(@mean,IT_ses);
        ITmedian_ses(subnum,iSes,:) = cellfun(@median,IT_ses);
        ITstd_ses(subnum,iSes,:) = cellfun(@std,IT_ses);
      end
    end

end

% save key metrics to look at covariates later:
hmm_1stlevel = [];
hmm_1stlevel.FO = FracOcc;
hmm_1stlevel.LT_mu = LTmerged;
hmm_1stlevel.LT_med = LTmedian;
hmm_1stlevel.LT_std = LTstd;
hmm_1stlevel.IT_mu = ITmerged;
hmm_1stlevel.IT_med = ITmedian;
hmm_1stlevel.IT_std = ITstd;
if whichstudy==7 ||whichstudy==8
  hmm_1stlevel.per_ses.FO = FracOcc_ses;
  hmm_1stlevel.per_ses.LT_mu = LTmerged_ses;
  hmm_1stlevel.per_ses.LT_med = LTmedian_ses;
  hmm_1stlevel.per_ses.LT_std = LTstd_ses;
  hmm_1stlevel.per_ses.IT_mu = ITmerged_ses;
  hmm_1stlevel.per_ses.IT_med = ITmedian_ses;
  hmm_1stlevel.per_ses.IT_std = ITstd_ses;
end
%% load MT PSDs for each state:
diagselect = find(eye(config.parc.n_parcels));
offdiagselect = find(~eye(config.parc.n_parcels));
if whichstudy==1
    fname = [config.resultsdir, 'coherence_state_ordering'];
else
    fname = [config.resultsdir, 'study1matched_state_ordering';];
end

loadHMMspectra_MT

% prepare sepctra and topos
tmp1 = nanmean(psd,4);
C = coh(:,:,:,offdiagselect);
C = nanmean(C,4);
for whichstate = 1:K
    pow_state_freq{whichstate} = squeeze(tmp1(:,whichstate,:));
    coh_state_freq{whichstate} = squeeze(C(:,whichstate,:));
    if use_WB_nnmf
        pow_state_topo{whichstate} = squeeze(nanmean((psd_wb(:,whichstate,:)),1));
        coh_state_topo{whichstate} = squeeze(mean(coh_wb(:,whichstate,:,:),1));
    else
        pow_state_topo{whichstate} = squeeze(nanmean(nanmean((psd(:,whichstate,:,:)),3),1));
        coh_state_topo{whichstate} = squeeze(mean(mean(coh(:,whichstate,:,:,:),3),1));
    end
end
clear C tmp1

group_avg_coh_wb = squeeze(mean(coh_wb));
group_avg_psd_wb = squeeze(mean(psd_wb));
save([config.resultsdir, 'group_avg_spectral_maps.mat'], 'group_avg_coh_wb', 'group_avg_psd_wb')

%% Compute long term assymetry:
[FO_intervals,FO_pvals,t_intervals,FO_stat] = computeLongTermAsymmetry(vpath,hmmT,K);

hmm_1stlevel.FO_intervals = FO_intervals;
hmm_1stlevel.assym_permtest = FO_stat;
hmm_1stlevel.assym_permtest.pvals = FO_pvals;
a=[];
for i=1:K
    for j=1:K
        [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(FO_intervals(i,j,1,:)), squeeze(FO_intervals(i,j,2,:)));
    end
end
hmm_1stlevel.assym_ttest = a;

% correcting for the number of tests. Correcting for the two tails in the
% FO_assym is done inside the permutation test (FO_permutation_test)
bonf_ncomparisons = (K.^2-K);
alpha_thresh = (0.05/bonf_ncomparisons);
hmm_1stlevel.assym_permtest.alpha_thresh = alpha_thresh;
hmm_1stlevel.assym_permtest.sigpoints = hmm_1stlevel.assym_permtest.pvals<alpha_thresh;
if whichstudy==3
    alpha_thresh = 0.1*alpha_thresh;
elseif whichstudy==4
    alpha_thresh = 0.000000001 * alpha_thresh;
elseif whichstudy==6
    alpha_thresh=10^-20 * alpha_thresh;
end
hmm_1stlevel.assym_ttest.alpha_thresh = alpha_thresh;
hmm_1stlevel.assym_ttest.sigpoints = hmm_1stlevel.assym_ttest.pvals<alpha_thresh;


% Run TINDA on the group level.
[FO_group,~,~] = computeLongTermAsymmetry({cat(1,vpath{:})},{squash(cat(2,hmmT{:}))},K);
hmm_1stlevel.FO_intervals_group = FO_group;


%% Find the optimial ordering
% this script determines the optimal state ordering for a circular plot; it
% then determines whether such a sequentially organised network could arise
% by chance by random shuffles of the rows of the transmat
if strcmp(config.reordering_states, 'coherence')
    optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_coherence' ,'.mat'];
elseif strcmp(config.reordering_states, 'study1matched')
    optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_study1matched' ,'.mat'];
else
    optimalseqfile = [config.hmmfolder,'bestseq',int2str(whichstudy),'.mat'];
end
if ~isfile(optimalseqfile)
    bestsequencemetrics = optimiseSequentialPattern(FO_intervals);
    save(optimalseqfile,'bestsequencemetrics');
else
    load(optimalseqfile);
end
bestseq = bestsequencemetrics{1};

if strcmp(config.reordering_states, 'replay') && whichstudy<3
    % put the lowest coherence state at 12 o'clock
    bestseq=circshift(bestseq, 12-find(bestseq==12)+1);
end
angleplot = circle_angles(bestseq);

%% Compute TINDA metrics

hmm_1stlevel.cycle_metrics = compute_tinda_metrics(config, bestseq, angleplot, FO_intervals, hmm_1stlevel.assym_ttest.sigpoints, color_scheme);

hmm_1stlevel.cycle_metrics_no_ordering = compute_tinda_metrics(config, [1 12:-1:2], circle_angles([1 12:-1:2]), FO_intervals, hmm_1stlevel.assym_ttest.sigpoints, color_scheme);

clear a
hmm_1stlevel.metric_vs_no_ordering=[];
[a.h, a.pvals, a.ci, a.stat] = ttest(hmm_1stlevel.cycle_metrics.rotational_momentum, hmm_1stlevel.cycle_metrics_no_ordering.rotational_momentum, 'Tail', 1);
hmm_1stlevel.metric_vs_no_ordering.rotational_momentum = a;

clear a
% compute cycle_metrics per session
if whichstudy==3 
    for iSes=1:config.nSes
        for iSj=1:config.nSj
            tmp = cumsum(hmmT{iSj});
            vpath_ses{iSj, iSes} = vpath{iSj}(1+tmp(iSes)-hmmT{iSj}(iSes) : tmp(iSes));
            hmmT_ses{iSj, iSes} = hmmT{iSj}(iSes);
        end
        a=[];
        [a.FO_intervals,a.assym_permtest.FO_pvals,a.t_intervals,a.assym_permtest.FO_stat] = computeLongTermAsymmetry(vpath_ses(:, iSes),hmmT_ses(:,iSes),K);
        
        hmm_1stlevel.tinda_per_ses{iSes} = a;
        b=[];
        for i=1:K
            for j=1:K
                [b.h(i,j), b.pvals(i,j), b.ci(i,j,:), b.stat(i,j)] = ttest(squeeze(a.FO_intervals(i,j,1,:)), squeeze(a.FO_intervals(i,j,2,:)));
            end
        end
        hmm_1stlevel.tinda_per_ses{iSes}.assym_ttest = b;
        hmm_1stlevel.tinda_per_ses{iSes}.cycle_metrics = compute_tinda_metrics(config, [], angleplot, a.FO_intervals, b.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme);
    end
elseif whichstudy==7
  for iSes=1:config.nSes
    [FO_intervals_ses,~,t_intervals_ses,~] = computeLongTermAsymmetry(vpath_ses(iSes:config.nSes:end),hmmT_ses(iSes:config.nSes:end),K,[],[],[],false);
    hmm_1stlevel.per_ses.FO_intervals(:,:,:,iSes,:) = FO_intervals_ses;
    hmm_1stlevel.per_ses.t_intervals(iSes,:,:) = t_intervals_ses;
    hmm_1stlevel.per_ses.cycle_metrics(iSes) = compute_tinda_metrics(config, bestseq, angleplot, FO_intervals_ses, hmm_1stlevel.assym_ttest.sigpoints, color_scheme);
  end
end



%% Run TINDA on permuted state labels
rng(42);
nperm = 1000;
for iperm=1:nperm
    iperm
    for k=1:length(vpath)
        vpath_perm{k} = zeros(size(vpath{k}));
        perm_ordering = randperm(K);
        for i=1:K
            vpath_perm{k}(vpath{k}==i)=perm_ordering(i);
        end
    end
    [FO_intervals_perm{iperm},FO_pvals_perm{iperm},~,FO_stat_perm{iperm}] = computeLongTermAsymmetry(vpath_perm,hmmT,K);
    bestsequencemetrics_perm{iperm} = optimiseSequentialPattern(FO_intervals_perm{iperm});
    angleplot_perm = circle_angles(bestsequencemetrics_perm{iperm}{1});
    FO_assym_perm = squeeze((FO_intervals_perm{iperm}(:,:,1,:)-FO_intervals_perm{iperm}(:,:,2,:))./mean(FO_intervals_perm{iperm},3));
    rotational_momentum_perm(iperm,:) = compute_rotational_momentum(angleplot_perm, FO_assym_perm);
end

hmm_1stlevel.perm = [];
hmm_1stlevel.perm.FO_intervals = FO_intervals_perm;
hmm_1stlevel.perm.assym_permtest.pvals = FO_pvals_perm;
hmm_1stlevel.perm.assym_permtest.stat = FO_stat_perm;
for k=1:nperm
    a=[];
    for i=1:K
        for j=1:K
            [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(hmm_1stlevel.perm.FO_intervals{k}(i,j,1,:)), squeeze(hmm_1stlevel.perm.FO_intervals{k}(i,j,2,:)));
        end
    end
    hmm_1stlevel.perm.assym_ttest{k} = a;
end
hmm_1stlevel.perm.bestsequencemetrics = bestsequencemetrics_perm;
hmm_1stlevel.perm.rotational_momentum = rotational_momentum_perm;
hmm_1stlevel.perm.obs_vs_perm = max([1./(nperm+1), (sum(mean(hmm_1stlevel.perm.rotational_momentum ,2)>mean(hmm_1stlevel.cycle_metrics.rotational_momentum)))/nperm]);

hmm_1stlevel.perm.obs_no_ordering_vs_perm = max([1./(nperm+1), (sum(mean(hmm_1stlevel.perm.rotational_momentum ,2)>mean(hmm_1stlevel.cycle_metrics_no_ordering.rotational_momentum)))/nperm]);


%% Run TINDA on the vpath simulated from the transprob matrix.
SequenceAnalysis_transprob_simulations
% saved in hmm_1stlevel.FO_state_simulation

%% Run TINDA with regressing out trial effects
if whichstudy==7 || whichstudy==8
  load([config.resultsdir, 'wakehen_events_badseg_adapted.mat'])

  % trigger info:
  % 5-7: Famous face
  % 13-15: Unfamiliar face
  % 17-19: Scrambled face
  
  tmp = cat(1,events{:});
  event_codes = unique(tmp(:,3));
  vis_events = [5:7, 13:15, 17:19];
  mot_events = setdiff(event_codes, vis_events)';

  % get reaction times and state that was active at moment of visual
  % stimulus and button press
  wakeman_henson_correlation;

  
  % average reaction time 932 ms
  % we want minimal overlap between the two windows
  event_twindow = {[0, 125], [-125, 75]};
  baseline_twindow = {[-26, -1], [-151 -126]}; 
  event_grouping = {vis_events, mot_events}; % first one is visual, second one button press
  for j=1:length(event_grouping)
      event_id = event_grouping{j};
      for k=1:length(events)
          % only take events with these IDs 5:7, 13:15, 17:19
          idx = find(any(events{k}(:,3)==event_id,2));
          grouped_events{j,k} = events{k}(idx, [1 3]);
            
          % remove events at the border
          id = grouped_events{j,k}(:,1)<(-event_twindow{j}(1)) | grouped_events{j,k}(:,1)>=length(vpath_ses{k})-event_twindow{j}(2);
          grouped_events{j,k} = grouped_events{j,k}(~id,:);
      end
  end
  trialdef = [];
  trialdef.events = grouped_events;
  trialdef.twindow = event_twindow;
  trialdef.baseline = baseline_twindow;
    
  % correct for both visual and motor erf
  [FO_intervals,~,~,~, FO_residual] = computeLongTermAsymmetry(vpath_ses,hmmT_ses,K,[],true, trialdef,false);
  
  hmm_1stlevel.per_ses.trial_glm.FO_intervals = reshape(FO_residual, [K,K,2,config.nSes, config.nSj]);
  hmm_1stlevel.trial_glm.FO_intervals = squeeze(mean(hmm_1stlevel.per_ses.trial_glm.FO_intervals,4));
  a=[];
  for i=1:K
    for j=1:K
      [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(hmm_1stlevel.trial_glm.FO_intervals(i,j,1,:)), squeeze(hmm_1stlevel.trial_glm.FO_intervals(i,j,2,:)));
    end
  end
  hmm_1stlevel.trial_glm.assym_ttest = a;
  hmm_1stlevel.trial_glm.assym_ttest.alpha_thresh = alpha_thresh;
  hmm_1stlevel.trial_glm.assym_ttest.sigpoints = hmm_1stlevel.trial_glm.assym_ttest.pvals<alpha_thresh;
  hmm_1stlevel.trial_glm.bestsequencemetrics = optimiseSequentialPattern(hmm_1stlevel.trial_glm.FO_intervals);
  hmm_1stlevel.trial_glm.cycle_metrics = compute_tinda_metrics(config, hmm_1stlevel.trial_glm.bestsequencemetrics{1}, circle_angles(hmm_1stlevel.trial_glm.bestsequencemetrics{1}), hmm_1stlevel.trial_glm.FO_intervals, hmm_1stlevel.trial_glm.assym_ttest.sigpoints, color_scheme);
end

%% plot HMM summary statistics
figure_supp_hmm_stats

%% plot HMM state spectra vs freq (averaged over all parcels):
figure_supp_hmm_spectra


%% create the TINDA example figure in a seperate script
figure1_tinda_example


%% Figure 1 Supplement: Plot each figure separately with power and coherence maps
% also creates a figure with all states' circle plots separately
figure_supp_tinda_states


%% Figure 2: plot circular diagram
figure2_circleplot


%% Plot circle metrics
figure_supp_tinda_metrics


%% Figure 2 supplement:  analyse by quintiles
figure_supp_tinda_quintiles

if whichstudy==1
%% Figure Supplement 2: analyse by intervals >2 heartbeats long
figure_supp_tinda_heartbeat


%% Figure 2 supplement: analyse intervals <half a respiratory cycle
figure_supp_tinda_respiration
end

%% save metrics
if exist([config.resultsdir, 'HMMsummarymetrics'])
    save([config.resultsdir, 'HMMsummarymetrics'],'hmm_1stlevel','-append');
else
    save([config.resultsdir, 'HMMsummarymetrics'],'hmm_1stlevel')
end


%% Figure 3: Spectral information in the circle plot

figure3_spectral_circle

if whichstudy==7 || whichstudy==8
  dotrialglm=true;
  figure3_spectral_circle
end


%% Figure 3 supplement: TINDA Movie
%{
outname = [config.resultsdir, 'TINDA.avi'];
tinda_movie(bestseq, mean_direction, sigpoints, f, psd, coh, outname)
%}

%% Compare studies
figure_supp_state_pathways

% states are matched to study 1. How well do the circle plots reproduce?

