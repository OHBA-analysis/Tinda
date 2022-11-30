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
use_WB_nnmf=true; % whether or not to use the wide band NNMF Diego's Nature Comms to select power and coherence (alternative is selecting 1-30 Hz)
useMT = true; % use the MT results instead of Cam's wavelet approach

%% Load HMM results
% first load data and plot basic temporal statistics:
temp = load(fullfile(config.hmmfolder,config.hmmfilename));

fname_stateorder = 'coherence_state_ordering';
if useMT
  if use_WB_nnmf
    fname_stateorder = [fname_stateorder, '_MT_nnmf'];
  else
    fname_stateorder = [fname_stateorder, '_MT'];
  end
end

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
    load([config.resultsdir, fname_stateorder])
  else
    new_state_ordering=1:hmm.K;
  end
  hmm = hmm_permutestates(hmm, new_state_ordering);
  for i=1:config.nSj
    hmmT{i} = sum(hmm.subj_inds==i);
  end
elseif whichstudy==3
  if strcmp(config.reordering_states, 'coherence')
    load([config.resultsdir, fname_stateorder])
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
    load([config.resultsdir, fname_stateorder])
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

%% load MT PSDs for each state:
diagselect = find(eye(config.parc.n_parcels));
offdiagselect = find(~eye(config.parc.n_parcels));
fname = [config.resultsdir, 'coherence_state_ordering'];

loadHMMspectra_MT

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
hmm_1stlevel.FO_intervals_group = FO_group;

%% Find the optimial ordering
% this script determines the optimal state ordering for a circular plot; it
% then determines whether such a sequentially organised network could arise
% by chance by random shuffles of the rows of the transmat
if strcmp(config.reordering_states, 'coherence')
  optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_coherence' ,'.mat'];
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

measures = {'FO_assym_subject_fit', 'TIDA', 'rotational_momentum', 'circularity_subject', 'TIDA_perstate', 'rotational_momentum_perstate'};
for im = measures
  m = im{1};
  if strcmp(m, 'FO_assym_subject_fit') || strcmp(m, 'TIDA') ||... 
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


%% Figure Supplement 2: analyse by intervals >2 heartbeats long
figure_supp_tinda_heartbeat


%% Figure 2 supplement: analyse intervals <half a respiratory cycle
figure_supp_tinda_respiration


%% save metrics
if exist([config.resultsdir, 'HMMsummarymetrics'])
  save([config.resultsdir, 'HMMsummarymetrics'],'hmm_1stlevel','-append');
else
  save([config.resultsdir, 'HMMsummarymetrics'],'hmm_1stlevel')
end


%% Figure 3: Spectral information in the circle plot
figure3_spectral_circle



%% Figure 3 supplement: TINDA Movie
%{
outname = [config.resultsdir, 'TINDA.avi'];
tinda_movie(bestseq, mean_direction, sigpoints, f, psd, coh, outname)
%}