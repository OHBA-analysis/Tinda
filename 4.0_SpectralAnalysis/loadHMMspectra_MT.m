if useMT
  if whichstudy==4
    f = h5read([config.resultsdir, 'spectra/spectra.h5'], '/f');
    psd = h5read([config.resultsdir, 'spectra/spectra.h5'], '/psd');
    coh = h5read([config.resultsdir, 'spectra/spectra.h5'], '/coh');
    subj_weight = h5read([config.resultsdir, 'spectra/spectra.h5'], '/w');
    wb_comp = h5read([config.resultsdir, 'spectra/spectra.h5'], '/wb_comp');
    nnmf_psd = h5read([config.resultsdir, 'spectra/spectra.h5'], '/nnmf_psd');
    nnmf_coh = h5read([config.resultsdir, 'spectra/spectra.h5'], '/nnmf_coh');
  else
    load([config.resultsdir, 'spectra/spectra.mat']);
    subj_weight = w; % fraction of total training data per subject (for weighting FO)
  end
  % permute state numbers
  if use_WB_nnmf
    fname = [fname, '_MT_nnmf'];
    if exist([fname, '.mat'], 'file')
      load(fname)
    else
      [~, new_state_ordering] = sort(mean(reshape(permute(nnmf_coh(:,:,1,:,:), [2,1,3,4,5]), 12,[]),2),'descend');
      save(fname, 'new_state_ordering')
    end
  else
    fname = [fname, '_MT'];
    if exist([fname, '.mat'], 'file')
      load(fname)
    else
      [~, new_state_ordering] = sort(mean(reshape(permute(coh, [2,1,3,4,5]), 12,[]),2),'descend');
      save(fname, 'new_state_ordering')
    end
  end
  
  psd = psd(:,new_state_ordering,:,:);
  coh = coh(:,new_state_ordering,:,:,:);
  nnmf_psd = nnmf_psd(:,new_state_ordering,:,:);
  nnmf_coh = nnmf_coh(:,new_state_ordering,:,:,:);
  
  f_orig = f;
%   psd = abs(psd(:,:,1:nearest(f,30), :));
%   coh = coh(:,:,1:nearest(f,30), :,:);
%   f = f(1:nearest(f, 30));
  coh(:,:,:,diagselect)=0;
  sqrtf=sqrt(f);
  psd_wb = squeeze(nnmf_psd(:,:,1,:));
  coh_wb = squeeze(nnmf_coh(:,:,1,:,:));
  
  % get the static power and coherence, i.e. weighted by FO
  sz=size(psd);
  group_FO = subj_weight*hmm_1stlevel.FO;
  static_pow = repmat(group_FO, [sz(1),1, sz([3 4])]) .* psd;
  static_pow_wb = repmat(group_FO, [sz(1),1, sz(4)]) .* psd_wb;
  
  sz=size(coh);
  static_coh = repmat(group_FO, [sz(1),1, sz([3 4,5])]) .* coh;
  static_coh_wb = repmat(group_FO, [sz(1),1, sz(4), sz(4)]) .* coh_wb;
  
  % now get the parcel/frequency averages for plotting
  powAvg_freq = nanmean(squeeze(sum(nanmean(static_pow,4),2)));
  cohAvg_freq = nanmean(squeeze(sum(nanmean(static_coh(:,:,:,offdiagselect),4),2)));
  
  if use_WB_nnmf
    powAvg_topo = squeeze(nanmean(sum(static_pow_wb,2),1));
    cohAvg_topo = squeeze(nanmean(sum(static_coh_wb,2),1));
  else
    powAvg_topo = squeeze(nanmean(sum(nanmean(static_pow,3),2),1));
    cohAvg_topo = squeeze(nanmean(sum(nanmean(static_coh,3),2),1));
  end
  
else
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
  
  powAvg_topo = squeeze(mean(mean(psd_wb,2)));
  cohAvg_topo = squeeze(mean(mean(coh_wb,2)));
end