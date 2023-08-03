% Load in spectra
if whichstudy==4 || whichstudy==6
  coh = h5read([config.resultsdir, 'spectra/spectra.h5'], '/coh');
  coh = permute(coh, [5,4,3,1,2]);
  f = h5read([config.resultsdir, 'spectra/spectra.h5'], '/f');
  psd = h5read([config.resultsdir, 'spectra/spectra.h5'], '/psd');
  if whichstudy==4
    psd = double(permute(psd.r, [4,3,2,1]));
  elseif whichstudy==6
    psd = double(permute(psd, [4,3,2,1]));
  end
  subj_weight = h5read([config.resultsdir, 'spectra/spectra.h5'], '/w')';
  wb_comp = h5read([config.resultsdir, 'spectra/spectra.h5'], '/wb_comp');
  nnmf_psd = h5read([config.resultsdir, 'spectra/spectra.h5'], '/nnmf_psd');
  nnmf_psd = permute(nnmf_psd, [4,3,2,1]);
  nnmf_coh = h5read([config.resultsdir, 'spectra/spectra.h5'], '/nnmf_coh');
  nnmf_coh = permute(nnmf_coh, [5,4,3,1,2]);
else
  if exist([config.resultsdir, 'spectra/spectra.mat'])
    load([config.resultsdir, 'spectra/spectra.mat']);
    subj_weight = w; % fraction of total training data per subject (for weighting FO)
  else
    % put the stuff for loadHMMspectra here
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
      elseif strcmp(config.reordering_states, 'study1matched')
        % load template coherence
        config1 = getStudyDetails(1);
        template = load([config1.resultsdir, 'group_avg_spectral_maps.mat']);
        match_states(template.coh, group_avg_coh_wb)
        % manually match with study1 coherence

        save(fname, 'new_state_ordering')
        P = pow_no_ordering(:, new_state_ordering,:,:,:);
        coh = coh_no_ordering(:, new_state_ordering,:,:,:);
      else
        [P, coh, f] = loadHMMspectra(config,whichstudy,hmm,run_inds,[], false);
      end
      [P,coh,f] = loadHMMspectra(config,whichstudy,hmm,run_inds,[],false);
    elseif whichstudy==4 || whichstudy==6
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

    powAvg_topo = squeeze(mean(mean(psd_wb,2)));
    cohAvg_topo = squeeze(mean(mean(coh_wb,2)));

    clear static_pow static_coh
  end
end

%% permute state numbers
if use_WB_nnmf
  fname = [fname, '_MT_nnmf'];
else
  fname = [fname, '_MT'];
end

if exist([fname, '.mat'], 'file')
  load(fname)
else
  if strcmp(config.reordering_states, 'coherence')
    if use_WB_nnmf
      [~, new_state_ordering] = sort(mean(reshape(permute(nnmf_coh(:,:,1,:,:), [2,1,3,4,5]), 12,[]),2),'descend');
    else
      [~, new_state_ordering] = sort(mean(reshape(permute(coh, [2,1,3,4,5]), 12,[]),2),'descend');
    end
  elseif strcmp(config.reordering_states, 'study1matched') % matched to MT_WB_nnmf
    if whichstudy==3 % manual ordering because a different parcellation was used
      new_state_ordering = [12, 11, 10, 9, 8, 3, 1, 6, 4, 7, 5, 2];
    elseif whichstudy==7
      new_state_ordering = [3     8     2    10     7     1     4    12     6     9    11     5];
    else
      if whichstudy==8
        config1 = getStudyDetails(6);
      else
        config1 = getStudyDetails(1);
      end
      template = load([config1.resultsdir, 'group_avg_spectral_maps.mat']);
      new_state_ordering = match_states(template.group_avg_coh_wb, squeeze(mean(nnmf_coh(:,:,1,:,:),1)));
    end
  end
  save(fname, 'new_state_ordering')
end

psd = psd(:,new_state_ordering,:,:);
coh = coh(:,new_state_ordering,:,:,:);
nnmf_psd = nnmf_psd(:,new_state_ordering,:,:);
nnmf_coh = nnmf_coh(:,new_state_ordering,:,:,:);

f_orig = f;
coh(:,:,:,diagselect)=0;
sqrtf=sqrt(f);
psd_wb = squeeze(nnmf_psd(:,:,1,:));
coh_wb = squeeze(nnmf_coh(:,:,1,:,:));

% get the static power and coherence, i.e. weighted by FO
sz=size(psd);
if whichstudy==7 || whichstudy==8
  [s0,s1,s2,s3] = size(psd);
  group_FO = subj_weight*transpose(reshape(permute(hmm_1stlevel.per_ses.FO, [3,2,1]), s1, []));
else
  group_FO = subj_weight*hmm_1stlevel.FO;
end
static_pow = repmat(group_FO, [sz(1),1, sz([3 4])]) .* psd;
static_pow_wb = repmat(group_FO, [sz(1),1, sz(4)]) .* psd_wb;

sz=size(coh);
static_coh = repmat(group_FO, [sz(1),1, sz([3 4,5])]) .* coh;
static_coh_wb = repmat(group_FO, [sz(1),1, sz(4), sz(4)]) .* coh_wb;

% now get the parcel/frequency averages for plotting
powAvg_freq = nanmean(squeeze(sum(nanmean(static_pow,4),2)));
cohAvg_freq = nanmean(squeeze(sum(static_coh,2)));
cohAvg_freq = nanmean(cohAvg_freq(:,:,offdiagselect),3); % breaking up in two lines helps with memory

if use_WB_nnmf
  powAvg_topo = squeeze(nanmean(sum(static_pow_wb,2),1));
  cohAvg_topo = squeeze(nanmean(sum(static_coh_wb,2),1));
else
  powAvg_topo = squeeze(nanmean(sum(nanmean(static_pow,3),2),1));
  cohAvg_topo = squeeze(nanmean(sum(nanmean(static_coh,3),2),1));
end
if whichstudy==7
  % average over sessions
  coh = permute(squeeze(mean(reshape(permute(coh, [2:5,1]), s1,s2,s3,s3,config.nSes, config.nSj),5)), [5 1:4]);
  psd = permute(squeeze(mean(reshape(permute(psd, [2:4,1]), s1,s2,s3,config.nSes, config.nSj),4)), [4 1:3]);
  coh_wb = permute(squeeze(mean(reshape(permute(coh_wb, [2:4,1]), s1,s3,s3,config.nSes, config.nSj),4)), [4 1:3]);
  psd_wb = permute(squeeze(mean(reshape(permute(psd_wb, [2:3,1]), s1,s3,config.nSes, config.nSj),3)), [3 1:2]);
end

clear static_pow static_coh static_pow_wb static_coh_wb nnmf_coh nnmf_psd

