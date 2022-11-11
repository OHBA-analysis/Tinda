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
%%
% summary plots:
fig = setup_figure([],1,2.25);
subplot(3,1,1)
distributionPlot(FracOcc,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
set(gca,'YLim',[0 1.1*max(FracOcc(:))]);
xlabel('RSN-state'), ylabel('Proportion'), title('Fractional Occupancy')
grid on;

subplot(3,1,2)
distributionPlot(LTmerged ./ config.sample_rate * 1000,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Life Times'); xlabel('RSN-state'), ylabel('Time (ms)'),
grid on;
YL = 1.1*max(LTmerged(:))./ config.sample_rate * 1000;
set(gca,'YLim',[0 YL]);

subplot(3,1,3);
% NOTE: if we want to use a logscale in the axis, enable this
%{
distributionPlot(log10(ITmerged ./ config.sample_rate * 1000),'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Interval Times');xlabel('RSN-state'), ylabel('Time (ms)');
grid on
YL(2) =1.5* max(mean(log10(ITmerged ./ config.sample_rate * 1000)));
YL(1) = min(squash(log10(ITmerged ./ config.sample_rate * 1000)));
set(gca,'YLim',YL)
set(gca,'YTick',(log10(1000*[0.05,0.1,0.5,1,5,10])))
y_labels = get(gca,'YTickLabel');
for i=1:length(y_labels)
  y_labels{i}=num2str(10.^(str2num(y_labels{i})),1);
end
set(gca,'YTickLabels',y_labels);
%}
% print([config.figdir '0_temporalstats_IT_logscale'],'-depsc')

distributionPlot(ITmerged ./ config.sample_rate*1000,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Interval Times');xlabel('RSN-state'), ylabel('Time (ms)');grid on
YL(2) =1.5* max(mean((ITmerged ./ config.sample_rate * 1000)));
YL(1) = 0;
set(gca,'YLim',YL)

save_figure(fig, [config.figdir '0_HMM_summary_statistics'])

%% and actually compute long term assymetry:

[FO,FO_pvals,t_intervals,FO_stat] = computeLongTermAsymmetry(vpath,hmmT,K);

bonf_ncomparisons = 2*(K.^2-K);

mean_direction = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
mean_assym = squeeze(mean((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3),4));

hmm_1stlevel.FO_intervals = FO;
hmm_1stlevel.FO_assym = squeeze((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3));
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
if 1%~isfile(optimalseqfile)
  bestsequencemetrics = optimiseSequentialPattern(FO);
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
  [simulation{iperm}.FO,simulation{iperm}.pvals,simulation{iperm}.t_intervals] = computeLongTermAsymmetry(simulation_vpath,hmmT,K);
  simulation{iperm}.mean_direction = squeeze(mean(simulation{iperm}.FO(:,:,1,:)-simulation{iperm}.FO(:,:,2,:),4));
  simulation{iperm}.mean_assym = squeeze(nanmean((simulation{iperm}.FO(:,:,1,:)-simulation{iperm}.FO(:,:,2,:))./mean(simulation{iperm}.FO,3),4));
  simulation{iperm}.FO_assym = squeeze((simulation{iperm}.FO(:,:,1,:)-simulation{iperm}.FO(:,:,2,:))./mean(simulation{iperm}.FO,3));
  simulation{iperm}.rotational_momentum = squeeze(imag(sum(sum(angleplot.*simulation{iperm}.FO_assym))));
  simulation{iperm}.TIDA = nanmean(abs(reshape((simulation{iperm}.FO_assym), K*K,[])))';
  for k=1:K
    simulation{iperm}.TIDA_perstate(:,k) = nanmean(abs([squeeze(simulation{iperm}.FO_assym(:,k,:));squeeze(simulation{iperm}.FO_assym(k,:,:))]));
  end
  simulation{iperm}.bestsequencemetrics = optimiseSequentialPattern(simulation{iperm}.FO);
end
hmm_1stlevel.FO_stats_simulation = simulation;

% potentially average over simulations before calculating all these
% measures
simulation_average=[];
for k=1:n_sim_perm
  simulation_average.FO(:,:,:,:,k) = simulation{k}.FO;
end
simulation_average.FO = mean(simulation_average.FO,5); 
[simulation_average.FO_pvals, simulation_average.FO_stat] = FO_permutation_test(simulation_average.FO, K, config.nSj);
simulation_average.mean_direction = squeeze(mean(simulation_average.FO(:,:,1,:)-simulation_average.FO(:,:,2,:),4));
simulation_average.mean_assym = squeeze(nanmean((simulation_average.FO(:,:,1,:)-simulation_average.FO(:,:,2,:))./mean(simulation_average.FO,3),4));
simulation_average.FO_assym = squeeze((simulation_average.FO(:,:,1,:)-simulation_average.FO(:,:,2,:))./mean(simulation_average.FO,3));
simulation_average.rotational_momentum = squeeze(imag(sum(sum(angleplot.*simulation_average.FO_assym))));
simulation_average.TIDA = nanmean(abs(reshape((simulation_average.FO_assym), K*K,[])))';
for k=1:K
  simulation_average.TIDA_perstate(:,k) = nanmean(abs([squeeze(simulation_average.FO_assym(:,k,:));squeeze(simulation_average.FO_assym(k,:,:))]));
end
simulation_average.bestsequencemetrics = optimiseSequentialPattern(simulation_average.FO);
hmm_1stlevel.FO_stats_simulation_average = simulation_average;

% Alternatively, do the FO assym on the group level
[FO_group,~,~] = computeLongTermAsymmetry({cat(1,vpath{:})},{squash(cat(1,hmmT{:}))},K);
hmm_1stlevel.group_FO = FO_group;
% we can either simulate the vpath based on the group level transprob, or
% the individual vpath
% let do the group:
nsim=1000;
for sim=1:nsim
  sim
  clear hmmT_sim vpath_sim
if 1
  hmmT_sim = cat(1,hmmT{:});
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

%% plot average state spectra vs freq (averaged over all parcels):

Pmean = squeeze(mean(mean(psd,1),4));

% relorabs = 'rel'; % 'rel' or 'abs'
for relorabs = {'abs', 'rel'}
  fig = setup_figure([],1,1);
if strcmp(relorabs, 'abs')
  sup = '';
elseif strcmp(relorabs, 'rel')
  sup = '_relative';
end
ls = {'-','--'};
for i=1:hmm.K
  if strcmp(relorabs, 'rel')
    plot(sqrtf,Pmean(i,:),'Color',color_scheme{i},'LineStyle',ls{1+(i>6)},'LineWidth',2);
    set_sqrt_ax(f)
    xlim(sqrt([f(1), 30]))
  else
    plot(f,Pmean(i,:),'Color',color_scheme{i},'LineStyle',ls{1+(i>6)},'LineWidth',2);
    xlim([f(1), 30])
  end
  hold on;
  leglabels{i} = ['State ',int2str(i)];
end

xlabel('Freq (Hz)'), ylabel('PSD');
yticks([])
legend(leglabels, 'Location', 'NorthOutside', 'NumColumns', K/4)
set_font(10, {'title', 'label'})
save_figure([config.figdir,'1supp_PSDperstate', sup]);
end

%% Make a power vs coherence plot

fig = setup_figure([],1,.75); hold on
for k=1:K
    scatter((squeeze(log10(nanmean(nanmean((psd(:,k,:,:)),4),3)))), log10(squeeze(nanmean(nanmean((coh(:,k,:,offdiagselect)),4),3))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
  l{k} = sprintf('State %d', k);
end
% axis off
box off
yticks([])
ylabel('Coherence')
xticks([])
xlabel('PSD')
legend(l, 'Location', 'EastOutside', 'NumColumns', 1)
box off
axis square
set_font(10, {'label', 'title'})
xlim(log10([min(min((squeeze(nanmean(nanmean((psd),4),3)))))*0.95, max(max((squeeze(nanmean(nanmean((psd),4),3)))))*1.05]))
ylim(log10([min(min((squeeze(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),3)))))*1.05, max(max((squeeze(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),3)))))*0.95]))

save_figure([config.figdir,'1supp_PSDvsCoh'] )


%% fit NNMF over space and find modes:

if use_WB_nnmf
P2 = reshape(squeeze(mean(psd,1)),[12*size(Pmean,2),config.parc.n_parcels]);

maxP = 6;
if strcmp(config.reordering_states, 'coherence')
  nnmffile = [config.hmmfolder,'nnmf_spatialfit_coherence.mat'];
else
  nnmffile = [config.hmmfolder,'nnmf_spatialfit.mat'];
end
if isfile(nnmffile)
  load(nnmffile)
else
  for iter=1:10
    fprintf(['Run ',int2str(iter),'\n'])
    [atemp btemp]=nnmf_mww(P2,maxP,'replicates',500,'algorithm','als');
    res = mean(mean((P2 - (atemp*btemp)).^2));
    if iter==1 || res<lastres
      a=atemp;
      b=btemp;
      lastres = res;
    end
  end
  clear atemp btemp
  a = reshape(a,[K,size(P2,1)/K,maxP]);
  save(nnmffile,'a','b','maxP', 'f');
end

end

%% Figure 1: Plot TINDA example
if whichstudy==1
  iSj=5;
  if strcmp(config.reordering_states, 'coherence')
    whichstate=1;
  else
    whichstate=2;
  end
  fprintf(['\n Subject: ',int2str(iSj), ', State: ' int2str(whichstate)]);
  
  % load raw data
  D = spm_eeg_load(replace(hmm.data_files{iSj}, '/Users/chiggins/data/YunzheData/Replaydata4Cam/WooliePipeline/spm/', '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/'));
  
  % get vpath
  Gamma_subj = Gamma(hmm.subj_inds==iSj,:);
  [~,vpath_subj] = max(Gamma_subj,[],2);
  
  
  % parcellation
  parc=config.parc;
  mni_coords = config.parc.roi_centers;
  nparcels = length(parc.labels);
  
  % Figure setup
  fig = setup_figure([],2,0.6); clear ax
  cmap = colormap('inferno');
  local_clim = true; % for the power maps, create a clim based on that state's range
  
  %%%%%%%%%%%%%%%%%%%
  % DATA TRACE PLOT %
  %%%%%%%%%%%%%%%%%%%
  ax(1) = axes('Position', [.1 .55 .25 .4]); hold on
  
  tlength=720;
  if (whichstate==1 && strcmp(config.reordering_states, 'replay')) || (whichstate==2 && strcmp(config.reordering_states, 'coherence'))
    tstart=34300;
    % this is the interval we'll plot TINDA for
    t1 = 106;
    t2 = 668;
  elseif (whichstate==2 && strcmp(config.reordering_states, 'replay')) || (whichstate==1 && strcmp(config.reordering_states, 'coherence'))
    tstart=30300;
    t1 = 76;
    t2 = 697;
  end
  t_segment = tstart + [1:tlength];
  thalf = mean([t1, t2]);
  t_segment1 = tstart + [t1:tlength/2];
  t_segment2 = tstart + [tlength/2+1:t2];
  
  q = vpath_subj(t_segment)==whichstate;
  s1 = [find(diff(q)==1)+1 find(diff(q)==-1)];
  mn = min(min(D(1:8,t_segment)));
  mx = max(max(D(1:8,t_segment)));
  
%   plot_Gamma(Gamma_subj(t_segment,:), t_segment, true,false,[]), colororder(cat(1,color_scheme{:}))  
%   a = gca;
%   for k=1:length(a.Children)
%     a.Children(k).FaceAlpha=0.2;
%     a.Children(k).LineWidth=0.001;
%   end
%   box off
%   ylim([.02 .99])

  % plot raw trace
  plot((1/250):(1/250):(length(t_segment)/250), D(1:8,t_segment(1:tlength))')
  % annotate visits to whichstate
  for jj=1:size(s1,1)
    fill([s1(jj,1) s1(jj,1) s1(jj,2) s1(jj,2)]./250, 1.3*[mn mx mx mn], color_scheme{whichstate}, 'FaceAlpha', .5, 'EdgeColor', 'none');
  end
  set(ax(1),'xcolor','none')
  ax(1).XAxis.Label.Color=[0 0 0];
  ax(1).XAxis.Label.Visible='on';
  ax(1).YAxis.Label.Color=[0 0 0];
  set(ax(1),'YTick',[]);
  axis off
  text(ax(1).Position(1)-0.3, ax(1).Position(2)+0.05, {'Resting', 'state', 'data'}, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10)

  %%%%%%%%%%%%%%
  % VPATH PLOT %
  %%%%%%%%%%%%%%
  ax(2) = axes('Position', [0.1 0.1 0.25 0.4]);
  hold on
  
  % annotate intervals
  cb = [256,193,1; 201, 204, 231]/256; % interval annotation colors
  fill([t1-1 t1-1 thalf thalf]./250, [0 13 13 0], cb(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none') % interval T1
  fill([thalf thalf t2 t2]./250, [0 13 13 0], cb(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none')% interval T2
  
  yl = [0 13];
  % plot vpath
  for k=1:12
    myline = NaN(length(t_segment),1);
    myline(vpath_subj(t_segment)==k)=(13-k);
    plot((1/250):(1/250):(length(t_segment)/250),myline,'LineWidth',20,'Color',color_scheme{k});
    yaxislabels{13-k} = ['State ',int2str(k)];
  end
  set(ax(2),'YTick',[1:12]);
  yaxislabels{13-whichstate} = ['\color{red} ' yaxislabels{13-whichstate}];
  set(ax(2),'YTickLabels',yaxislabels);
  grid off
  set(ax(2),'xcolor','none')
  ax(2).XAxis.Label.Color=[0 0 0];
  ax(2).XAxis.Label.Visible='on';
  xlabel('Time');
  ax(2).YAxis.Label.Color=[0 0 0];
  ylim(yl)
  title('HMM state timecourse', 'HorizontalAlignment', 'center')
%   title('FO', 'HorizontalAlignment', 'center')
%   title('\boldmath$FO^T$', 'Interpreter', 'Latex', 'FontSize', 10)

  text(0.28, -.06, sprintf('%s', sym("T_1")), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold')
  text(0.72, -.06, sprintf('%s', sym("T_2")), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold')

  
  %%%%%%%%%%%%%%%%
  % T1 BAR GRAPH %
  %%%%%%%%%%%%%%%%
  ax(3) = axes('Position', [0.375 0.1 0.075 0.4]);
  hold on
  
  
  for ii=1:12
    T1(ii) = sum(vpath_subj(t_segment1)==ii);
    h=barh(ii, T1(ii));
    set(h,'FaceColor', color_scheme{ii})
  end
  set(gca, 'Ydir', 'reverse')
  axis off
  ylim(yl)
%   title({'$$\ \ \ \ \ \ \ \forall_{k\neq{1}}\sum_{T_{1}|T_{2}}Z_{t}==k$$'},'Interpreter', 'Latex', 'FontSize', 10, 'HorizontalAlignment', 'center')
  title('\boldmath$FO^{T_1}$', 'Interpreter', 'Latex', 'FontSize', 10, 'FontWeight', 'bold')

  %%%%%%%%%%%%%%%%
  % T2 BAR GRAPH %
  %%%%%%%%%%%%%%%%
  ax(4) = axes('Position', [0.475 0.1 0.075 0.4]);
  hold on
  for ii=1:12
    T2(ii) = sum(vpath_subj(t_segment2)==ii);
    h=barh(ii, T2(ii));
    set(h,'FaceColor', color_scheme{ii})
  end
  set(gca, 'Ydir', 'reverse')
  axis off
  ylim(yl)
  title('\boldmath$FO^{T_2}$', 'Interpreter', 'Latex', 'FontSize', 10, 'FontWeight', 'bold')

  %%%%%%%%%%%%%%%%%%%%%
  % DISTRIBUTION PLOT %
  %%%%%%%%%%%%%%%%%%%%%
  sigpoints = pvals<(0.05/bonf_ncomparisons);
  for ii=1:12
    ax(4+ii) = axes('Position', [0.575 0.12+(12-ii)*0.032 0.1 0.025]);
    hold on
  end
  
  clear d
  cb = [256,193,1; 201, 204, 231]/256;
  for ii=1:12
    axes(ax(4+ii))
    if ii==1
%         title({'$$\forall_{k\neq{1}} \sum_{i=1}^{N}Q_{i,T_{1}}-Q_{i,T_{2}}$$'},'Interpreter', 'Latex', 'FontSize', 10, 'HorizontalAlignment', 'center')
      title('FO asym', 'FontSize', 10)
    end
    if any(setdiff(1:12, whichstate)==ii)
      d{1} = squeeze(hmm_1stlevel.FO_intervals(whichstate,ii,1,:));
      d{2} = squeeze(hmm_1stlevel.FO_intervals(whichstate,ii,2,:));
      % exagerate significant differences
      if sigpoints(ii,whichstate)
        if mean_direction(ii,whichstate)<0
          d{1} = d{1}+0.5*mean(d{2});
        elseif mean_direction(ii,whichstate)>0
          d{2} = d{2}+0.5*mean(d{2});
        end
      end
      h2 = raincloud_plot(d{2}, 'box_on', 0, 'color', cb(2,:), 'alpha', 0.5);
      h1 = raincloud_plot(d{1}, 'box_on', 0, 'color', cb(1,:), 'alpha', 0.5);
      y2 = get(ax(4+ii), 'YLim');
      ylim([0 y2(2)*1.3])
    end
    axis off
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
  end

  %%%%%%%%%%%%%%%%%%%%%
  % TINDA CIRCLE PLOT %
  %%%%%%%%%%%%%%%%%%%%%
  ax(25) = axes('Position', [0.70 0.3500 0.3 0.145]); axis off
  title('TINDA', 'FontSize', 10)
  ax(17) = axes('Position', [0.70 0.0500 0.3 0.4]);
  cyclicalstateplot_perstate(bestseq,mean_direction,sigpoints,find(bestseq==whichstate),false);
% seq=[12:-1:1];
%   cyclicalstateplot_perstate(seq,mean_direction,sigpoints,find(seq==whichstate),false);


  %%%%%%%%%%%
  % PSD MAP %
  %%%%%%%%%%%


  ax(18) = axes('Position',[0.34        0.8 0.175 0.2] ); % top left
  ax(19) = axes('Position',[0.34+0.185  0.8  0.175 0.2] ); % top right
  ax(20) = axes('Position',[0.34+0.025  0.6  0.175 0.2] ); % bottom left
  ax(21) = axes('Position',[0.34+0.16   0.6  0.175 0.2] ); % bottom right
  ax(26) = axes('Position', [0.37 0.80 0.3 0.15]); axis off
  title('PSD', 'FontSize', 10)

  
  net_mean = zeros(nparcels,size(Gamma,2));
  thresh_mean = zeros(size(Gamma,2),1);
  for k = 1:K
    if use_WB_nnmf
      net_mean(:,k) = squeeze(nnmfWB_res.nnmf_psd_maps(k,1,:))';
    else
      net_mean(:,k) = squeeze(nanmean(nanmean((psd(:,k,:,:)),3),1))';
    end
  end
  net_mean=log10(net_mean);
  toplot = net_mean(:,whichstate);%-mean(net_mean,2);
  psdthresh = min(net_mean(:)); % lowest value
  if local_clim
    CL = [min(squash(toplot(:,:))) max(squash(toplot(:,:)))];
  else
    CL = [min(squash(net_mean(:,:))) max(squash(net_mean(:,:)))];
  end
  f2 = plot_surface_4way(parc,toplot,0,false,'trilinear',[],psdthresh,CL,ax(18:21));
   
  %%%%%%%%%%%%%%%%%
  % COHERENCE MAP %
  %%%%%%%%%%%%%%%%%
  ax(22) = axes('Position',[0.63+0.02 0.765 0.24 0.24] ); % left
  ax(23) = axes('Position',[0.63+0.18  0.765 0.24 0.24] ); % right
  ax(24) = axes('Position',[.74 0.6 0.22 0.22]);% bottom
 
  
  if use_WB_nnmf
    graph = abs(squeeze(nnmfWB_res.nnmf_coh_maps(whichstate,1,:,:)));
  else
    graph = squeeze(nanmean(nanmean((coh(:,whichstate,:,:,:)),3),1));
    graph(diagselect)=0;
  end
  graph_mean = squeeze(nanmean(nanmean(nanmean((coh),3),1),2));
  graph_mean(diagselect)=0;
  [~, ax(22:24), ~] = plot_coh_topo(ax(22:24), mni_coords, graph, graph_mean);
   ax(27) = axes('Position', [0.7 0.80 0.3 0.15]); axis off
  title('Coh', 'FontSize', 10, 'HorizontalAlignment', 'center')
  
  
  % change colormap for power
  for ii=18:21
    colormap(ax(ii), 'inferno')
  end

  %%%%%%%%%%%%
  % SAVE FIG %
  %%%%%%%%%%%%%
  set_font(10, {'label', 'title'})
  save_figure(fig, [config.figdir '/1_tinda_example_state',int2str(whichstate)]);
%   close
end
%% Figure 1 Supplement: Plot each figure separately with power and coherence maps
parc=config.parc;
nparcels=config.parc.n_parcels;
local_clim=1;
cmap = colormap('inferno');
mni_coords = config.parc.roi_centers;
if whichstudy==4
alpha = 0.0000001;
else
  alpha = 0.05;
end

net_mean = zeros(nparcels,K);
thresh_mean = zeros(K,1);
for k = 1:K
  if use_WB_nnmf
    net_mean(:,k) = squeeze(nnmfWB_res.nnmf_psd_maps(k,1,:))';
  else
    net_mean(:,k) = log10(squeeze(nanmean(nanmean((psd(:,k,:,:)),3),1))');
  end
end
if whichstudy==1
  statelabels = {'DMN', 'parietal alpha', 'fronto-temporal delta/theta', 'visual alpha', 'temporal theta', 'fronto-parietal', 'sensorimotor beta', 'parietal', 'anterior', 'parieto-temporal', 'sensorimotor', 'visual'};
end

graph_mean = squeeze(nanmean(nanmean(nanmean((coh),1),2),3));
graph_mean_triu=squash(triu(graph_mean));
for whichstate =1:K

  fig = setup_figure([],2,0.75); axis off
  suptitle(sprintf('RSN state %d - %s', whichstate, statelabels{whichstate}))
% TINDA
  ax(11) = axes('Position', [0.35 0.74 0.3 0.145]); axis off
  title('TINDA', 'FontSize', 10)
  ax(1) = axes('Position', [0.325 0.5, 0.35, 0.35]);
  cyclicalstateplot_perstate(bestseq,mean_direction,pvals<(alpha/bonf_ncomparisons),find(bestseq==whichstate),false);

  ax(9) = axes('Position', [0.05 0.53, 0.25, 0.4]);cla; hold on
  if use_WB_nnmf
    pow = nnmfWB_res.nnmf_psd_specs(1,1:nearest(f,30));
    f=f(1:nearest(f,30));
  else
    powAvg = log10(squeeze(nanmean(nanmean(psd(:,:,:,:),4),2)));
    pow = log10(squeeze(nanmean((psd(:,whichstate,:,:)),4)));
  end
  plot(f,mean(powAvg), '--', 'Color', [.8 .8 .8])
  shadedErrorBar(f,mean(pow,1), std(pow,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  set(gca, 'YTick', [])
  set(gca, 'Xtick', nearest(f,30)/4:nearest(f,30)/4:nearest(f,30))
  xlabel('Frequency (Hz)')
  ylabel('PSD')
  [yl(1) yl(2)] = bounds(log10(squash(squeeze(nanmean(nanmean((psd(:,:,nearest(f,2):end,:)),1),4)))));
  ylim([0.95 1.05].*yl)
  box off

  ax(2) = axes('Position',[0    0.25  0.21 0.21]); % top left
  ax(3) = axes('Position',[0.23  0.25  0.21 0.21]); % top right
  ax(4) = axes('Position',[0.03  0.05 0.21 0.21]);% bottom left
  ax(5) = axes('Position',[0.2   0.05 0.21 0.21]); % bottom right
  
  
  toplot = net_mean(:,whichstate);%-mean(net_mean,2);
  psdthresh = min(net_mean(:))*0.9; % lowest value
  if local_clim
    CL = [min(squash(toplot(:,:))) max(squash(toplot(:,:)))];
  else
    CL = [min(squash(net_mean(:,:))) max(squash(net_mean(:,:)))];
  end
  f2 = plot_surface_4way(parc,toplot,0,false,'trilinear',[],psdthresh,CL,ax(2:5));
  
  
  % coherence topo
  if use_WB_nnmf
    graph = abs(squeeze(nnmfWB_res.nnmf_coh_maps(whichstate,1,:,:)));
  else
    graph = squeeze(nanmean(nanmean(coh(:,whichstate,:,:,:),3),1));
  end
  ax(6) = axes('Position',[0.5 0.225  0.24 0.24]);cla
  ax(7) = axes('Position',[0.75  0.225  0.24 0.24]);cla
  ax(8) = axes('Position',[0.625 0.05 0.24 0.24]);cla;
  [~, ax(6:8), ~] = plot_coh_topo(ax(6:8), mni_coords, graph, graph_mean);

  
  ax(10) = axes('Position', [0.73 0.53, 0.25, 0.4]);cla; hold on
  if use_WB_nnmf
    C = nanmean(nnmfWB_res.nnmf_coh_maps(whichstate,1,offdiagselect),3) * nnmfWB_res.nnmf_coh_specs(1,1:nearest(f,30));
  else
    C = coh(:,:,:,offdiagselect);
    Cavg = squeeze(nanmean(nanmean(C(:,:,:,:),4),2));
    C = squeeze(nanmean(C(:,whichstate,:,:),4));
  end
  plot(f, mean(Cavg), '--', 'Color', [.8 .8 .8])
  
  shadedErrorBar(f,mean(C,1), std(C,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  set(gca, 'YTick', [])
  set(gca, 'Xtick', nearest(f,30)/4:nearest(f,30)/4:nearest(f,30))
  xlabel('Frequency (Hz)')
  ylabel('Coherence')
  [yl(1) yl(2)] = bounds(squash(squeeze(nanmean(nanmean(coh(:,:,nearest(f,2):end,offdiagselect),1),4))));
  ylim([0.95, 1.05].*yl)
  box off
  
  
  % colormap for power
  for ii=2:5
    colormap(ax(ii), 'inferno')
  end
  set_font(10, {'title', 'label'})
  save_figure(fig, [config.figdir '/1supp_tinda_state',int2str(whichstate)],false);
  
  % also save the one with relative x axis
  axes(ax(9))
  cla
  plot(sqrtf,mean(powAvg), '--', 'Color', [.8 .8 .8])
  shadedErrorBar(sqrtf,mean(pow,1), std(pow,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  set_sqrt_ax(f)
  xlim(sqrtf([1 end]))
  
  axes(ax(10))
  cla
  plot(sqrtf, mean(Cavg), '--', 'Color', [.8 .8 .8])
  shadedErrorBar(sqrtf,mean(C,1), std(C,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  set_sqrt_ax(f)
  xlim(sqrtf([1 end]))
  save_figure(fig, [config.figdir '/1supp_tinda_state',int2str(whichstate), '_relative'],false);

end


%% Fig 1 supplement: plot as multiple individual state plots:
fig=setup_figure([],2,1);
if whichstudy<4
  cyclicalstateplot_perstate(bestseq,mean_direction,pvals<(0.05/bonf_ncomparisons),[],false,color_scheme);
else
  cyclicalstateplot_perstate(bestseq,mean_direction,pvals<0.0000001*(0.05/bonf_ncomparisons),[],false, color_scheme);
end
save_figure([config.figdir,'1supp_StatePathways']);

% also print for legend:
figure('Position', [440 579 114 219]);
quiver(0,0,1,0,'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.8);hold on;
quiver(0,1,1,0,'Color',[0 0 0]+0.8,'LineWidth',2,'MaxHeadSize',0.8);hold on;
axis off;
print([config.figdir,'1supp_StatePathways_legend'], '-dpng')





%% Figure 2: plot circular diagram
if whichstudy~=4
  sigpoints = FO_pvals<(0.05/bonf_ncomparisons);
else
  sigpoints = FO_pvals<0.0000001*(0.05/bonf_ncomparisons);
end



for ext = {'rotationalmomentum', 'FOasym'}  
fig = setup_figure([],2,0.5);
clear ax
ax(1) = axes('Position', [0.025, 0, .3, 1]);
cyclicalstateplot(bestseq,mean_direction, sigpoints,color_scheme,false);
text(-.3, 1.3, 'observed')

ax(2) = axes('Position', [0.375, 0, .3, 1]);
if whichstudy==4
  cyclicalstateplot(simulation_average.bestsequencemetrics{2},simulation_average.mean_direction, simulation_average.FO_pvals<0.0000001*(0.05/bonf_ncomparisons),color_scheme,false);
else
  cyclicalstateplot(simulation_average.bestsequencemetrics{2},simulation_average.mean_direction, simulation_average.FO_pvals<(0.05/bonf_ncomparisons),color_scheme,false);
end
text(-.3, 1.3, 'simulated')

  if strcmp(ext, 'FOasym')
    tmp = [simulation_average.TIDA, hmm_1stlevel.TIDA];
    ylb = 'Mean FO asymmetry';
    [stat.H, stat.P, stat.CI, stat.stats] = ttest(tmp(:,2),tmp(:,1), 'tail', 'right');
  elseif strcmp(ext, 'rotationalmomentum')
    tmp = [simulation_average.rotational_momentum, hmm_1stlevel.rotational_momentum];
    ylb = 'Rotational momentum';
    [stat.H, stat.P, stat.CI, stat.stats] = ttest(tmp(:,2),tmp(:,1));
  end
ax(3) = axes('Position', [0.755, 0.15, 0.25, 0.7]); hold on
vline(0, '--k')

clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
for k=1:2
  scatter(ones(size(tmp(:,k))).*(k+(rand(size(tmp(:,k)))-0.5)/2),tmp(:,k), 'filled', 'MarkerFaceColor', clr{k}, 'MarkerFaceAlpha',0.6)
end
h=boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this back on top
set(h, 'linew', 2)
box off
if stat.H
  sigstar({[1,2]}, stat.P)
end
title({'Difference between', 'observed data and', 'simulations'})
set(ax(3), 'XTickLabels', {'simulated', 'observed'})
ylabel(ylb)
set_font(10, {'title', 'label'})
fname = [config.figdir,'2_yclicalpattern_', ext{1}];
save_figure(fname);
save([fname, 'stat'], 'stat')
end


% also make a boxplot seperately for each state (vs simulation)
clear stat
fig = setup_figure([],2,.75);
for whichstate=1:K+1
  if whichstate==K+1
    fig = setup_figure([],1,.6); hold on;
    tmp = [simulation_average.TIDA, hmm_1stlevel.TIDA];
  else
    subplot(4,3,whichstate), hold on
    tmp = [simulation_average.TIDA_perstate(:, whichstate), hmm_1stlevel.TIDA_perstate(:,whichstate)];
  end
  [stat{whichstate}.H, stat{whichstate}.P, stat{whichstate}.CI, stat{whichstate}.stats] = ttest(tmp(:,2),tmp(:,1), 'tail', 'right');
  clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
for k=1:2
  scatter(ones(size(tmp(:,k))).*(k+(rand(size(tmp(:,k)))-0.5)/2),tmp(:,k),'filled', 'MarkerFaceAlpha',0.6)
end
boxplot(tmp, 'Colors', 'k', 'Width', .75, 'whisker', 1000) % put this on top
if stat{whichstate}.H
  sigstar({[1,2]}, stat{whichstate}.P)
end
box off
% ylim([-.36, 0.36])
set(gca, 'XTickLabels', {'simulated', 'observed'})
ylabel({'Mean FO asym'})
if whichstate==K
  fname=[config.figdir,'2supp_TIDA_perstate'];
  save_figure(fname);
  save([fname, 'stat'], 'stat')
elseif whichstate==K+1
  fname=[config.figdir,'2supp_TIDA'];
  save_figure(fname);
  stat=stat{K+1};
  save([fname, 'stat'], 'stat')
end
end


% also plot the *real* circle plot individually
cyclicalstateplot(bestseq,mean_direction, sigpoints,color_scheme);
save_figure([config.figdir,'2ind_Cyclicalpattern']);

% as well as the simulated
cyclicalstateplot(simulation_average.bestsequencemetrics{2},simulation_average.mean_direction, simulation_average.FO_pvals<(0.05/bonf_ncomparisons), color_scheme);
save_figure([config.figdir,'2ind_Cyclicalpattern_simulated']);


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

%% Figure 2 supplement:  analyse by quartiles
percentiles = 0:20:100;
fig=setup_figure([],2,2);
if whichstudy==4
  alpha_thresh = 0.0000001;
else
  alpha_thresh = 0.05;
end
clear ax
for ip=1:length(percentiles)-1
  [FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));
  mean_direction_p = squeeze(mean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
  ax(ip,1) = axes('Position', [-.05, 1.02-0.2*ip, 0.45, 0.15]);
  cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(alpha_thresh/bonf_ncomparisons));
  ax(ip,2) = axes('Position', [0.45, 1.03-0.196*ip, 0.5, 0.15]);
  ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
  distributionPlot(sqrt(ITmerged),'showMM',2,'color',{color_scheme{1:size(FO,2)}});
  set(ax(ip,2),'YLim',[0 1.1*max(sqrt(ITmerged(:)))])
  if ip==5
    xlabel('RSN-State')
  end
  grid on; ylabel('Interval Times (s)')
  for ii=1:length(ax(ip,2).YTickLabel)
    ax(ip,2).YTickLabel{ii} = num2str(str2num(ax(ip,2).YTickLabel{ii})^2);
  end
end
set_font(10, {'title', 'label'})

save_figure([config.figdir,'2supp_Cyclicalpatterns_percentiled']);
close
hmm_1stlevel.FO_quartile = FO_p;

% Make an extra plot with median interval time vs rotational momentum?
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));

percentiles = 5:5:95;
for ip=1:length(percentiles)
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,[percentiles(ip)-4 percentiles(ip)+5]);
FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
ITmean(ip) = mean(mean(cellfun(@mean,t_intervals_p)));ITmean(ip) = ITmean(ip)./config.sample_rate;
% ITmedian(ip) = cellfun(@median,t_intervals_p);ITmedian(ip) = ITmedian(ip)./config.sample_rate;
rotational_momentum_p(ip,:) = squeeze(imag(nansum(nansum(angleplot.*FO_assym_p))));
end

fig = setup_figure([],2,0.6); 
subplot(1,2,1), 
shadedErrorBar(percentiles,mean(rotational_momentum_p,2), std(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
xlabel('Percentile'), ylabel('Rotational Momentum')
yticks([0])
subplot(1,2,2)
shadedErrorBar((ITmean),mean(rotational_momentum_p,2), std(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
yticks([0])
xlabel('Mean IT (s)'), ylabel('Rotational Momentum')
save_figure([config.figdir,'2supp_Cyclicalpatterns_rotationalmomentum_percentiled']);
close

%% Figure Supplement 2: analyse by intervals >2 heartbeats long

% assume HR are greater than 50BPM:
maxHR = (60/50*2);
percentile = [maxHR*config.sample_rate,NaN]; % shortest permissible intervals:

[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentile);
mean_direction_p = squeeze(nanmean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;

fig = setup_figure([], 2,0.4); clear ax
ax(1) = axes('Position', [-.03, 0.1, 0.45, 0.8]);
cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(0.05/bonf_ncomparisons));

ax(2) = axes('Position', [0.5, 0.15, 0.45, 0.7]);
distributionPlot(ITmerged,'showMM',2,'color',{color_scheme{1:size(FO,2)}});

set(gca,'YLim',[0 1.1*max(ITmerged(:))])
title('Interval times');xlabel('RSN-State'), ylabel('Interval time (sec)');grid on;
set_font(10, {'title', 'label'})

save_figure([config.figdir,'2supp_Cyclicalpatterns_greaterthan2sec'])
close

% save metrics:
hmm_1stlevel.FO_cardiaccontrol = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO,3));

%% Figure 2 supplement: analyse intervals <half a respiratory cycle

% assume resp cycle of 16BPM;
minRR = 60/16*0.5;
percentile = [NaN,minRR*config.sample_rate]; % shortest permissible intervals:
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentile);
mean_direction_p = squeeze(nanmean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;

fig = setup_figure([], 2,0.4); clear ax
ax(1) = axes('Position', [-.03, 0.1, 0.45, 0.8]);
cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(0.05/bonf_ncomparisons));

ax(2) = axes('Position', [0.5, 0.15, 0.45, 0.7]);
distributionPlot(ITmerged,'showMM',2,'color',{color_scheme{1:size(FO,2)}});

set(gca,'YLim',[0 1.1*max(ITmerged(:))])
title('Interval times');xlabel('RSN-State'), ylabel('Interval time (sec)');grid on;
set_font(10, {'title', 'label'})
save_figure([config.figdir,'2supp_Cyclicalpatterns_lessthan2.5sec'])

hmm_1stlevel.FO_respcontrol = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO,3));


% save metrics for later analysis:
if exist([config.figdir, 'HMMsummarymetrics'])
  save([config.figdir, 'HMMsummarymetrics'],'hmm_1stlevel','-append');
else
  save([config.figdir, 'HMMsummarymetrics'],'hmm_1stlevel')
end


%% Figure 3: Spectral information in the circle plot

diffmode = {'rel'}; % whether to plot psd/coh relative to average over states
use_sqrt_f = [true, false] ; % whether to plot PSD/coh with sqrt of f axis
statecolor = [true, false]; % whether to plot PSD/coh in the state color
for sc = statecolor 
  if sc
    sup2 = '_statecolor';
  else
    sup2 = '_samecolor';
  end
for ii = use_sqrt_f
  if ii
  sup1 = '_sqrtax';
  else
  sup1 = '_linax';
  end

pos = zeros(12,2);
for i=1:12
    temp = exp(sqrt(-1)*(i+2)/12*2*pi);
    pos(i,:) = [real(temp),imag(temp)];
end
% manually shift some positions
pos = (pos/4 + 0.25)*2;
pos(:,1) = pos(:,1)*0.88;
pos(pos(:,2)>0.55,2) = pos(pos(:,2)>0.55,2) - 0.05;
pos(pos(:,2)<0.4,2) = pos(pos(:,2)<0.4,2) + 0.05;
pos([2 6],1) = pos([2 6],1) - [0.075;0.075];
pos([8 12],1) = pos([8 12],1) + [0.075;0.075];

pos(bestseq,:) = pos;
fig = setup_figure([],2,1);
ax(13,1) = axes('Position', [0.295 0.19 0.38 .6]); hold on
axes(ax(13,1))
cyclicalstateplot(bestseq,mean_direction, zeros(12), color_scheme, false,[],true)
coh_avg = squeeze(nanmedian(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),1),2));
% pow_avg = squeeze(nanmedian(nanmean(nanmean(log10(psd),4),1),2));
pow_avg = squeeze(nanmedian(nanmean(nanmean((psd),4),1),2));

for k=1:12
  ax(k,2) = axes('Position', [0.125 0.025 0 0]+[0.85 0.85 1 1].*[pos(k,1), pos(k,2), 0.1 0.1]); hold on
  yyaxis left
  yyaxis right
  if sc
    cl{1} = color_scheme{k};
    cl{2} = cl{1};
    linestyle = ':';
    ax(k,2).YAxis(1).Color = 'k';
    ax(k,2).YAxis(2).Color = 'k';
  else
    cl{1} = [0 0.4470 0.7410];
    cl{2} = [0.8500 0.3250 0.0980];
    linestyle = '-';
  end
coh_state = squeeze(nanmean(nanmean((coh(:,k,:,offdiagselect)),4),1));
pow_state = squeeze(nanmean(nanmean((psd(:,k,:,:)),4),1));

yyaxis left
tmpP = pow_state;
if strcmp(diffmode, 'rel')
  tmpP = log10(tmpP./pow_avg);
  [ylt(1) ylt(2)] = bounds(squash(log10(squeeze(nanmean(nanmean((psd(:,:,:,:)),4),1))./pow_avg')));
  [ylt(3) ylt(4)] = bounds(squash(log10(squeeze(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),1))./coh_avg')));
  [yl(1) yl(2)] = bounds(ylt);
else
  [yl(1) yl(2)] = bounds(squash(squeeze(nanmean(nanmean(log10(psd(:,:,nearest(f,3):end,:)),1),4))));
end

if ii
  plot(sqrtf,tmpP, 'LineWidth',2, 'Color', cl{1});
  set_sqrt_ax(f)
else
  plot(f,tmpP, 'LineWidth',2, 'Color', cl{1});
end
  hline(0)

% ylim([0.95,1.05].*yl)
ylim([1.05,1.05].*yl)

yticks([])
if sc
  ylabel(['PSD ' char(8212)]);
else
  ylabel('PSD');
end
yyaxis right
tmpC = coh_state;
if strcmp(diffmode, 'rel')
  tmpC = log10(tmpC./coh_avg);
else
  [yl(1) yl(2)] = bounds(squash(squeeze(nanmean(nanmean((coh(:,:,nearest(f,3):end,offdiagselect)),1),4))));
end

if ii
  plot(sqrtf,tmpC, 'LineWidth',2, 'Color', cl{2}, 'LineStyle', linestyle);
  xlim(sqrt([f(1) 30]))
  set_sqrt_ax(f)
else
 plot(f,tmpC, 'LineWidth',2, 'Color', cl{2}, 'LineStyle', linestyle);
 xlim([f(1) 30])
end
ylim([1.05 1.05].*yl)
hline(0)

if sc
  ylabel(['Coh ' '---']);
else
  ylabel('Coh');
end
yticks([])
box off

toplot = log10(squeeze(nanmean(nanmean(psd(:,k,:,:),3),1)));
coords_left = mni_coords(:,1)<0;
thresh = prctile(toplot,75);
% project the right side to the left side (all non thresholded values are shown on the left projection)
toplot(1:2:end) = max([toplot(1:2:end), toplot(2:2:end)],[],2, 'omitnan');
toplot(toplot<thresh)=NaN;

[CL(1), CL(2)] = bounds(squash(log10(squeeze((nanmean(nanmean(psd,3),1))))));
ax(k,1) = axes('Position', [0 0.025 0 0]+[0.85 0.85 1 1].*[pos(k,1), pos(k,2), 0.1 0.1]);
plot_surface_4way(parc,toplot,1,true,'trilinear', [],CL(1)*0.9,CL, ax(k,1))
cmap = colormap(inferno);
colormap(cmap(50:end,:))

end

ax(11,1) = axes('Position', [0.375, 0.38, 0.25, 0.25]); hold on
clear l
for k=1:K
    scatter(log10(squeeze(nanmean(nanmean((psd(:,k,:,:)),4),3))), log10(squeeze(nanmean(nanmean(coh(:,k,:,offdiagselect),4),3))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
  l{k} = sprintf('State %d', k);
end
% axis off
box off
yticks([])
ylabel('Coherence')
xticks([])
xlabel('PSD')
box off
axis square
xlim(log10([min(min((squeeze(nanmean(nanmean((psd),4),3)))))*0.95, max(max((squeeze(nanmean(nanmean((psd),4),3)))))*1.05]))
ylim(log10([min(min((squeeze(nanmean(nanmean(coh(:,:,:,offdiagselect),4),3)))))*1.05, max(max((squeeze(nanmean(nanmean(coh(:,:,:,offdiagselect),4),3)))))*0.95]))

set_font(10, {'label', 'title'})
save_figure([config.figdir,'3_Spectral_circle', sup1, sup2], false)
end
end




%% Figure 3 supplement: TINDA Movie
if whichstudy==1
  outname = [config.figdir, 'TINDA.avi'];
  if use_WB_nnmf
    nnmfWB_res.coh_freq(1,:,:) = squeeze(mean(mean(mean(coh,5),4),1));
    tinda_movie(bestseq, mean_direction, sigpoints, f, nnmfWB_res, nnmfWB_res, outname)
  else
    tinda_movie(bestseq, mean_direction, sigpoints, f, psd, coh, outname)
  end
end



%% refit earlier spatial modes to camcan data:

if whichstudy==4
  MEGUKfolder = '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/CanonicalRS/250Hz/hmm_1to45hz/'
  nnmffile_MEGUK = [MEGUKfolder,'nnmf_spatialfit.mat'];
  load(nnmffile_MEGUK)
  
  opt = statset('maxiter',1);
  [b_new,a_new]=nnmf_mww(P2',maxP,'replicates',500,'algorithm','als','w0',b','options',opt);
  b_new = b_new';
  a_new = a_new';
  a_new = reshape(a_new,[12,size(P2,1)/12,maxP]);
  
  figure();
  bestorder = [1:6];
  plotCyclicalTimeFreqPattern(seq,config.parc,a_new(:,:,bestorder),b_new(bestorder,:),f, config)
  print([config.figdir,'2CPSD_modes_samefit'],'-dpng');
  
end

%% load age data from camcan and check for relationship with FO assymetry:
if whichstudy==4
  info = camcan_getparticipantinfo(config);
  
  % check how age predicts each FO assymetry measure:
  assym = (FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3);
  for k1=1:12
    for k2=1:12
      Y = permute(assym(k1,k2,1,:),[4,1,2,3]);
      X = info(:,1);
      [R,P] = corrcoef(X,Y);
      rho(k1,k2) = R(1,2);
      pval(k1,k2) = P(1,2);
    end
  end
  mask = pval<(0.05/bonf_ncomparisons);
  cyclicalstateplot(bestseq,rho,pval<(0.05/bonf_ncomparisons));
  
  print([config.figdir,'1D_Cyclicalpattern_AgeEffect'],'-dpng');
  
  % does overall assymetry increase/ decrease with age?
  temp = (permute(assym,[4,1,2,3]));
  metric = sum(sum(abs(temp(:,find(~eye(12)))),3),2);
  [R,P] = corrcoef(metric,X);
  figure('Position',[440 516 960 282]);
  subplot(1,2,1);
  scatter(X,metric,'filled');
  plot4paper('Age','Sequential assymetry')
  title(['rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
  % sequential patterns become slightly stronger with age
  % check if this due to outliers:
  outlier = abs(demean(metric))>2*std(metric);
  
  subplot(1,2,2);
  [R,P] = corrcoef(metric(~outlier),X(~outlier));
  scatter(X(~outlier),metric(~outlier),'filled');
  plot4paper('Age','Sequential assymetry')
  title(['Outliers removed: rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
  print([config.figdir,'1E_SequenceAssym_AgeEffect'],'-dpng');
  
  % or gender?
  p = anova1(metric(~outlier),info(~outlier,3));
  % significant effect for gender also
  
  % look at age within a gender group:
  figure('Position',[440 516 960 282]);
  subplot(1,2,1);
  outlier = abs(demean(metric))>2*std(metric) | info(:,3)==1;
  [R,P] = corrcoef(metric(~outlier),X(~outlier));
  scatter(X(~outlier),metric(~outlier),'filled');
  plot4paper('Age','Sequential assymetry')
  title(['Male group: rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
  subplot(1,2,2);
  outlier = abs(demean(metric))>2*std(metric) | info(:,3)==2;
  [R,P] = corrcoef(metric(~outlier),X(~outlier));
  scatter(X(~outlier),metric(~outlier),'filled');
  plot4paper('Age','Sequential assymetry')
  title(['Female group: rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
end


%% Circular Vpath analysis
%{
W=125;
W_half = (W-1)/2;
for i=1:length(vpath)
  vpathcircle{i}= getCircleVpath(vpath{i}, bestseq);
  tmp = zeros((W-1)/2,1);
  vpathcircle_window{i} = tmp;
  for k=W_half+1:length(vpathcircle{i})-W_half
    vpathcircle_window{i} = [vpathcircle_window{i}; mean(vpathcircle{i}(k-W_half:k+W_half))];
  end
  vpathcircle_window{i} = [vpathcircle_window{i}; tmp];
  
  
  % do Fourier analysis on the circle vpath
  [tmp, circlefreq, circletime] = fourieranalysis_circleVpath(vpathcircle{i}, []);
  circlepow(i,:,:) = squeeze(nanmean(abs(tmp).^2));
  circlespctrm(i,:,:) = squeeze(nanmean(tmp));
  
  [tmp, circlefreq, circletime] = fourieranalysis_circleVpath(vpathcircle_window{i}, []);
  circlepow_window(i,:,:) = squeeze(nanmean(abs(tmp).^2));
  circlespctrm_window(i,:,:) = squeeze(nanmean(tmp));
end

figure; plot(circlefreq, nanmean(circlepow,3)), hold on, plot(circlefreq, nanmean(nanmean(circlepow,3)), 'k', 'LineWidth', 2)
title('Circle vpath PSD per subject, CamCan Rest');
plot4paper('Frequenzy (Hz)', 'PSD');
print([config.figdir,'1F1_CircleFreq'],'-dpng');

figure; plot(circlefreq, nanmean(circlepow_window,3)), hold on, plot(circlefreq, nanmean(nanmean(circlepow_window,3)), 'k', 'LineWidth', 2)
title('Windowed circle vpath PSD per subject, CamCan Rest');
plot4paper('Frequenzy (Hz)', 'PSD');
print([config.figdir,'1F2_CircleFreq_windowed'],'-dpng');
%}
