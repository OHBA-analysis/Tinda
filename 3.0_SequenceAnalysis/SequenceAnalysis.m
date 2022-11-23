% script called to load viterbi paths inferred and hmm objects and run
% post-hoc sequence analysis:
if ~exist('whichstudy','var')
  whichstudy = 1; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
if whichstudy==1 || whichstudy==2
  color_scheme = set1_cols();
elseif whichstudy==3 || whichstudy==5
  tmp = circshift([0.400000000000000,0.760784313725490,0.647058823529412;0.988235294117647,0.552941176470588,0.384313725490196;0.552941176470588,0.627450980392157,0.796078431372549;0.905882352941177,0.541176470588235,0.764705882352941;0.650980392156863,0.847058823529412,0.329411764705882;1,0.850980392156863,0.184313725490196;0.898039215686275,0.768627450980392,0.580392156862745;0.701960784313725,0.701960784313725,0.701960784313725;0.400000000000000,0.260784313725490,0.647058823529412;0.738235294117647,0.352941176470588,0.384313725490196;0.552941176470588,0.827450980392157,0.946078431372549;0.605882352941177,0.541176470588235,0.164705882352941],-2,1);
  for k=1:12
    color_scheme{k} = tmp(k,:);
  end
elseif whichstudy==4
  tmp = brewermap(12, 'Set3');
  for k=1:12
    color_scheme{k} = tmp(k,:);
  end
end
cm = colormap(hot);
cm2 = cm;
cm2(:,1) = cm(:,3);
cm2(:,3) = cm(:,1);
hotcold=[flipud(cm2);cm];

% include option for simulations here - results to be saved in adjacent
% folder:
simtests = false;
if simtests
  config.figdir = strrep(config.figdir,['Study',int2str(whichstudy)],['Study',int2str(whichstudy),'_simtest']);
  mkdir(config.figdir);
end
use_WB_nnmf=false; % whether or not to use the wide band NNMF as in Higgins 2020 to select power and coherence (alternative is selecting 1-30 Hz)

%% part 1: First level HMM statistics
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

%% plot HMM summary statistics
figure_supp_hmm_stats

%% Compute long term assymetry:

[FO_intervals,FO_pvals,t_intervals,FO_stat] = computeLongTermAsymmetry(vpath,hmmT,K);

bonf_ncomparisons = 2*(K.^2-K);

mean_direction = squeeze(mean(FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:),4));
mean_assym = squeeze(mean((FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:))./mean(FO_intervals,3),4));

hmm_1stlevel.FO_intervals = FO_intervals;
hmm_1stlevel.FO_assym = squeeze((FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:))./mean(FO_intervals,3));
hmm_1stlevel.FO_stat = FO_stat;
hmm_1stlevel.FO_pvals = FO_pvals;

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
% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
  disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));

% This is one way to define the strength of the circularity.
rotational_momentum = imag(sum(sum(angleplot.*hmm_1stlevel.FO_assym)));
hmm_1stlevel.rotational_momentum = squeeze(rotational_momentum);
hmm_1stlevel.max_theoretical_rotational_momentum = imag(sum(sum(angleplot.*sign(imag(angleplot)))));

hmm_1stlevel.TIDA = nanmean(reshape(abs(hmm_1stlevel.FO_assym), K*K,[]))';
% also compute this per state
for k=1:K
  hmm_1stlevel.TIDA_perstate(:,k) = nanmean(abs([squeeze(hmm_1stlevel.FO_assym(:,k,:));squeeze(hmm_1stlevel.FO_assym(k,:,:))]));
end

% Run TINDA on the vpath simulated from the transprob matrix.
n_sim_perm = 100;
simulation=cell(n_sim_perm,1);
for iperm=1:n_sim_perm
  iperm
  for k=1:length(vpath)
    simulation_vpath{k} = simulateVpath(vpath{k},hmmT{k},K);
  end
  [simulation{iperm}.FO_intervals,simulation{iperm}.pvals,simulation{iperm}.t_intervals] = computeLongTermAsymmetry(simulation_vpath,hmmT,K);
  simulation{iperm}.mean_direction = squeeze(mean(simulation{iperm}.FO_intervals(:,:,1,:)-simulation{iperm}.FO_intervals(:,:,2,:),4));
  simulation{iperm}.mean_assym = squeeze(nanmean((simulation{iperm}.FO_intervals(:,:,1,:)-simulation{iperm}.FO_intervals(:,:,2,:))./mean(simulation{iperm}.FO_intervals,3),4));
  simulation{iperm}.FO_assym = squeeze((simulation{iperm}.FO_intervals(:,:,1,:)-simulation{iperm}.FO_intervals(:,:,2,:))./mean(simulation{iperm}.FO_intervals,3));
  simulation{iperm}.rotational_momentum = squeeze(imag(sum(sum(angleplot.*simulation{iperm}.FO_assym))));
  simulation{iperm}.TIDA = nanmean(abs(reshape((simulation{iperm}.FO_assym), K*K,[])))';
  for k=1:K
    simulation{iperm}.TIDA_perstate(:,k) = nanmean(abs([squeeze(simulation{iperm}.FO_assym(:,k,:));squeeze(simulation{iperm}.FO_assym(k,:,:))]));
  end
  simulation{iperm}.bestsequencemetrics = optimiseSequentialPattern(simulation{iperm}.FO_intervals);
end
hmm_1stlevel.FO_stats_simulation = simulation;

% potentially average over simulations before calculating all these
% measures
simulation_average=[];
for k=1:n_sim_perm
  simulation_average.FO_intervals(:,:,:,:,k) = simulation{k}.FO_intervals;
end
simulation_average.FO_intervals = mean(simulation_average.FO_intervals,5);
[simulation_average.FO_pvals, simulation_average.FO_stat] = FO_permutation_test(simulation_average.FO_intervals, K, config.nSj);
simulation_average.mean_direction = squeeze(mean(simulation_average.FO_intervals(:,:,1,:)-simulation_average.FO_intervals(:,:,2,:),4));
simulation_average.mean_assym = squeeze(nanmean((simulation_average.FO_intervals(:,:,1,:)-simulation_average.FO_intervals(:,:,2,:))./mean(simulation_average.FO_intervals,3),4));
simulation_average.FO_assym = squeeze((simulation_average.FO_intervals(:,:,1,:)-simulation_average.FO_intervals(:,:,2,:))./mean(simulation_average.FO_intervals,3));
simulation_average.rotational_momentum = squeeze(imag(sum(sum(angleplot.*simulation_average.FO_assym))));
simulation_average.TIDA = nanmean(abs(reshape((simulation_average.FO_assym), K*K,[])))';
for k=1:K
  simulation_average.TIDA_perstate(:,k) = nanmean(abs([squeeze(simulation_average.FO_assym(:,k,:));squeeze(simulation_average.FO_assym(k,:,:))]));
end
simulation_average.bestsequencemetrics = optimiseSequentialPattern(simulation_average.FO_intervals);
hmm_1stlevel.FO_stats_simulation_average = simulation_average;

% Alternatively, do the FO assym on the group level
[FO_group,~,~] = computeLongTermAsymmetry({cat(1,vpath{:})},{squash(cat(1,hmmT{:}))},K);
hmm_1stlevel.group_FO = FO_group;
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
hmm_1stlevel.group_FO_simulation = FO_group_sim;

hmm_1stlevel.FO_group_stats = [];
hmm_1stlevel.FO_group_stats.FO_group_asym=FO_group(:,:,1)-FO_group(:,:,2);
hmm_1stlevel.FO_group_stats.FO_group_asym_perm = squeeze(FO_group_sim(:,:,1,:)-FO_group_sim(:,:,2,:));
% find whether there is a difference with the simulation
pval_group = [];
for k=1:12
  for l=1:12
    hmm_1stlevel.FO_group_stats.pval_group_vs_sim(k,l) = 1-sum(FO_all_asym(k,l)>squeeze(FO_all_asym_perm(k,l,:)))/nsim;
  end
end




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


%% plot average state spectra vs freq (averaged over all parcels):
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



[circularity, circle_pval, ~, ~, fig] = geometric_circularity(bestseq, mean_direction, sigpoints,[],[],[],color_scheme);
set_font(10, {'title', 'label'})
save_figure([config.figdir,'2supp_CyclicalpatternVsPermutations']);

hmm_1stlevel.circularity = circularity;
hmm_1stlevel.circularity_pval = circle_pval;
tmp = hmm_1stlevel.FO_assym; tmp(isnan(tmp))=0;
for k=1:config.nSj
  [hmm_1stlevel.circularity_subject(k,1), hmm_1stlevel.circularity_pval_subject(k,1), ~, ~, ~] = geometric_circularity(bestseq, tmp(:, :,k), sigpoints,[],[],0, color_scheme);
end

% also for the simulation average
[circularity_sim, circle_pval_sim, ~, ~, fig] = geometric_circularity(simulation_average.bestsequencemetrics{2}, simulation_average.mean_direction, simulation_average.FO_pvals<(0.05/bonf_ncomparisons));
hmm_1stlevel.FO_stats_simulation_average.circularity = circularity_sim;
hmm_1stlevel.FO_stats_simulation_average.circularity_pval = circle_pval_sim;
set_font(10, {'title', 'label'})
save_figure([config.figdir,'2supp_CyclicalpatternVsPermutations_simulation']);
tmp = simulation_average.FO_assym; tmp(isnan(tmp))=0;
for k=1:config.nSj
  [hmm_1stlevel.FO_stats_simulation_average.circularity_subject(k,1), hmm_1stlevel.FO_stats_simulation_average.circularity_pval_subject(k,1), ~, ~, ~] = geometric_circularity(simulation_average.bestsequencemetrics{2}, tmp(:,:,k), simulation_average.FO_pvals<(0.05/bonf_ncomparisons),[],[],0, color_scheme);
end

% compare the observed circularity with the simulated one
dat1=[];
dat1.dimord = 'rpt_chan_time';
dat1.label{1} = 'circularity';
dat1.time=1;
dat2=dat1;
dat1.trial =  hmm_1stlevel.circularity_subject;
dat2.trial = hmm_1stlevel.FO_stats_simulation_average.circularity_subject;

cfg=[];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.design = [ones(1,config.nSj), 2*ones(1,config.nSj); 1:config.nSj, 1:config.nSj];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 100000;
cfg.tail = 1;
stat_c = ft_timelockstatistics(cfg, dat1, dat2);
hmm_1stlevel.FO_stats_simulation_average.circularity_stat_obs_vs_perm = stat_c;

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
