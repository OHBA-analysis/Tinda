% script called to load viterbi paths inferred and hmm objects and run
% post-hoc sequence analysis:
if ~exist('whichstudy','var')
  whichstudy = 1; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = colorscheme(whichstudy);
hotcold = cmap_hotcold;

% include option for simulations here - results to be saved in adjacent
% folder:
simtests = false;
if simtests
  config.figdir = strrep(config.figdir,['Study',int2str(whichstudy)],['Study',int2str(whichstudy),'_simtest']);
  mkdir(config.figdir);
end
use_WB_nnmf=false; % whether or not to use the wide band NNMF as in Higgins 2020 to select power and coherence (alternative is selecting 1-30 Hz)

%% Load HMM results
% first load data and plot basic temporal statistics:
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
  elseif strcmp(config.reordering_states, 'coherence')
    load([config.figdir, 'coherence_state_ordering.mat'])
  else
    new_state_ordering=1:hmm.K;
  end
  hmm = hmm_permutestates(hmm, new_state_ordering);
  for i=1:config.nSj
    hmmT{i} = sum(hmm.subj_inds==i);
  end
elseif whichstudy==3
  if strcmp(config.reordering_states, 'coherence')
    load([config.figdir, 'coherence_state_ordering.mat'])
  else
    new_state_ordering=1:hmm.K;
  end
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
  hmmTsubj = cell(config.nSj);
  for i=1:config.nSj
    hmmTsubj{i} = [hmmTold{1,i},hmmTold{2,i},hmmTold{3,i}]-(length(hmm.train.embeddedlags)-1);
  end
  hmmT = hmmTsubj;
  clear hmmTsubj hmmTold;
elseif whichstudy==4
  if strcmp(config.reordering_states, 'coherence')
    load([config.figdir, 'coherence_state_ordering.mat'])
  else
    new_state_ordering=1:hmm.K;
  end
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
for subnum=1:config.nSj
  fprintf(['\nProcessing subj ',int2str(subnum)]);
  if whichstudy~=4
    vpath{subnum} = hmm.statepath(hmm.subj_inds==subnum);
  else
    temp = load(mat_files_orth{subnum},'vpath');
    vpath{subnum} = temp.vpath;
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

%% load wavelet PSDs for each state:
diagselect = find(eye(config.parc.n_parcels));
offdiagselect = find(~eye(config.parc.n_parcels));
fname = [config.figdir, 'coherence_state_ordering.mat'];

if whichstudy==3
  % for HCP need to recompute run indices (each subject has multiple runs)
  run_inds = zeros(size(hmm.statepath));
  t_offset = 0;
  for i=1:length(hmmT)
    t_length = sum(hmmT{i}) - length(hmmT{1})*(length(hmm.train.embeddedlags)-1);
    run_inds(t_offset + [1:t_length]) = i;
    t_offset = t_offset + t_length;
  end
  if strcmp(config.reordering_states, 'coherence')
    if ~isfile(fname)
      [pow_no_ordering, coh_no_ordering, f] = loadHMMspectra(config,whichstudy,hmm,run_inds,[],false, false);
      [~, new_state_ordering] = sort(nanmean(nanmean(nanmean(coh_no_ordering(:,:,1:nearest(f,30),offdiagselect),4),3),1), 'descend');
      save(fname, 'new_state_ordering')
      P = pow_no_ordering(:, new_state_ordering,:,:,:);
      coh = coh_no_ordering(:, new_state_ordering,:,:,:);
    else
      [P, coh, f] = loadHMMspectra(config,whichstudy,hmm,run_inds,[], false);
    end
  else
    [P, coh, f] = loadHMMspectra(config,whichstudy,hmm,run_inds,[], false);
  end
  [P,coh,f] = loadHMMspectra(config,whichstudy,hmm,run_inds,[],false);
elseif whichstudy==4
  % compute FO per subj:
  for i=1:length(vpath)
    for k=1:K
      FO_subj(i,k) = mean(vpath{i}==k);
    end
  end
  if strcmp(config.reordering_states, 'coherence')
    if ~isfile(fname)
      [pow_no_ordering, coh_no_ordering, f] = loadHMMspectra(config,whichstudy,hmm,[],FO_subj,false, false);
      [~, new_state_ordering] = sort(nanmean(nanmean(nanmean(coh_no_ordering(:,:,1:nearest(f,30),offdiagselect),4),3),1), 'descend');
      save(fname, 'new_state_ordering')
      P = pow_no_ordering(:, new_state_ordering,:,:,:);
      coh = coh_no_ordering(:, new_state_ordering,:,:,:);
    else
      [P,coh,f] = loadHMMspectra(config,whichstudy,hmm,[],FO_subj,false);
    end
  else
    [P,coh,f] = loadHMMspectra(config,whichstudy,hmm,[],FO_subj,false);
  end
else
  % Find the coherence state ordering (low to high coherence)
  if strcmp(config.reordering_states, 'coherence')
    if ~isfile(fname)
      [pow_no_ordering, coh_no_ordering, f] = loadHMMspectra(config,whichstudy,hmm,hmm.subj_inds,[],false, false);
      [~, new_state_ordering] = sort(nanmean(nanmean(nanmean(coh_no_ordering(:,:,1:nearest(f,30),offdiagselect),4),3),1), 'descend');
      save(fname, 'new_state_ordering')
    end
  end
  [P,coh,f] = loadHMMspectra(config,whichstudy,hmm,hmm.subj_inds,[],false);
end

psd = abs(P(:,:,1:nearest(f,30), diagselect));
coh = coh(:,:,1:nearest(f,30), :,:);
coh(:,:,:,diagselect)=0;
f = f(1:nearest(f, 30));
sqrtf=sqrt(f);
%
% get the static power and coherence, i.e. weighted by FO
sz=size(psd);
static_pow = repmat(mean(hmm_1stlevel.FO), [sz(1),1, sz([3 4])]) .* psd;
sz=size(coh);
static_coh = repmat(mean(hmm_1stlevel.FO), [sz(1),1, sz([3 4,5])]) .* coh;
% now get the parcel/frequency averages for plotting
powAvg_freq = nanmean(squeeze(sum(nanmean(static_pow,4),2)));
powAvg_topo = squeeze(nanmean(sum(nanmean(static_pow,3),2),1));
cohAvg_freq = nanmean(squeeze(sum(nanmean(static_coh(:,:,:,offdiagselect),4),2)));
cohAvg_topo = squeeze(nanmean(sum(nanmean(static_coh,3),2),1));

%% Compute long term assymetry:

[FO_intervals,FO_pvals,t_intervals,FO_stat] = computeLongTermAsymmetry(vpath,hmmT,K);
hmm_1stlevel.FO_intervals = FO_intervals;
hmm_1stlevel.FO_stat = FO_stat;
hmm_1stlevel.FO_pvals = FO_pvals;

% correcting for the number of tests. Correcting for the two tails in the 
% FO_assym is done inside the permutation test (FO_permutation_test)
bonf_ncomparisons = (K.^2-K); 
alpha_thresh = 0.05;
if whichstudy==4
  alpha_thresh = 0.0000001 * alpha_thresh;
end
alpha_thresh = (alpha_thresh/bonf_ncomparisons);
sigpoints = FO_pvals<alpha_thresh;


% Run TINDA on the group level.
[FO_group,~,~] = computeLongTermAsymmetry({cat(1,vpath{:})},{squash(cat(1,hmmT{:}))},K);
hmm_1stlevel.group_FO = FO_group;

%% Find the optimial ordering
% this script determines the optimal state ordering for a circular plot; it
% then determines whether such a sequentially organised network could arise
% by chance by random shuffles of the rows of the transmat
if strcmp(config.reordering_states, 'coherence')
  optimalseqfile = [config.figdir,'bestseq',int2str(whichstudy),'_coherence' ,'.mat'];
else
  optimalseqfile = [config.hmmfolder,'bestseq',int2str(whichstudy),'.mat'];
end
if ~isfile(optimalseqfile)
  bestsequencemetrics = optimiseSequentialPattern(FO_intervals);
  save(optimalseqfile,'bestsequencemetrics');
else
  load(optimalseqfile);
end
bestseq = bestsequencemetrics{2};

if strcmp(config.reordering_states, 'replay') && whichstudy<3
  % put the lowest coherence state at 12 o'clock
  bestseq=circshift(bestseq, 12-find(bestseq==12)+1);
end
angleplot = circle_angles(bestseq);

%% Compute TINDA metrics

hmm_1stlevel.cycle_metrics = compute_tinda_metrics(config, bestseq, angleplot, FO_intervals, sigpoints, color_scheme);


%% Run TINDA on the vpath simulated from the transprob matrix.
SequenceAnalysis_transprob_simulations
% saved in hmm_1stlevel.FO_state_simulation


% and compare the observed metrics with the simulated ones
% TODO

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


% also make a boxplot seperately for each state (vs simulation)
figure_supp_tinda_metrics



%% Figure 2 supplement:  analyse by quintiles
figure_supp_tinda_quintiles


%% Figure Supplement 2: analyse by intervals >2 heartbeats long
figure_supp_tinda_heartbeat


%% Figure 2 supplement: analyse intervals <half a respiratory cycle
figure_supp_tinda_respiration


%% save metrics
if exist([config.figdir, 'HMMsummarymetrics'])
  save([config.figdir, 'HMMsummarymetrics'],'hmm_1stlevel','-append');
else
  save([config.figdir, 'HMMsummarymetrics'],'hmm_1stlevel')
end


%% Figure 3: Spectral information in the circle plot
figure3_spectral_circle



%% Figure 3 supplement: TINDA Movie

outname = [config.figdir, 'TINDA.avi'];
tinda_movie(bestseq, mean_direction, sigpoints, f, psd, coh, outname)
