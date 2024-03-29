% this script fits a second level sequential HMM to an already inferred
% viterbi path sequence
if ~exist('whichstudy','var')
  whichstudy = 1; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = set1_cols();

%% Load HMM results
use_WB_nnmf=true; % whether or not to use the wide band NNMF Diego's Nature Comms to select power and coherence (alternative is selecting 1-30 Hz)
useMT = true; % use the MT results instead of Cam's wavelet approach


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
  temp = load(fullfile(config.resultsdir,config.hmmfilename));

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

end

%% Set Poisson window length to average state lifetime:
load(config.metricfile, 'hmm_1stlevel')
try
  Poisswindow = ceil(mean(hmm_1stlevel.LT_mu(:)));
catch
  opts = [];
  opts.K = 12;
  opts.Fs = config.sample_rate;
  opts.dropstates=0;
  for k=1:length(vpath)
    lt = getStateLifeTimes(vpath{k}, hmmT{k}, opts);
    lt_mu(k,:) = cellfun(@mean, lt);
  end
    Poisswindow = ceil(mean(lt_mu(:)));
end

if whichstudy==7
  vpath = vpath_ses;
  hmmT = hmmT_ses;
  config.nSj = 114;
end

W=Poisswindow;
Noverlap = 1; % number of points to overlap over successive windows
T_poiss = [];
X_poiss = [];
Poiss_subj_inds = [];
t_last = 0;

if Noverlap~=0
  overlapstring='_overlappingWindows';
end
if whichstudy<4 || ~exist([config.resultsdir,'Poissdata_',int2str(W),overlapstring,'/filelist.mat'])
  for iSj=1:config.nSj
    vpTemp = vpath{iSj};
    t_last = t_last + sum(hmmT{iSj});
    fprintf(['\nProcessing subject ',int2str(iSj)]);
    if whichstudy==4 || whichstudy==6 %reset for each subject and save, rather than concatenate
      X_poiss = [];
      T_poiss = [];
    end
    for i=1:length(hmmT{iSj})
      % Event count over window of length W:
      if Noverlap==0
        T_poiss = [T_poiss;floor(hmmT{iSj}(i)/W)];
        on_t = 1 + sum(hmmT{iSj}(1:(i-1)));
      else
        % set for overlapping windows
        T_poiss = [T_poiss;hmmT{iSj}(i)-W+1];
        on_t = sum(hmmT{iSj}(1:(i-1)));
      end


      if T_poiss(end)>1
        if Noverlap==0
          temp = reshape(vpTemp(on_t:(W*T_poiss(end)+on_t-1),:),[W,T_poiss(end)]);
        else
          temp = zeros(W,T_poiss(end));
          for i2=1:(T_poiss(end))
            temp(:,i2) = vpTemp(on_t + [i2:(i2+W-1)]);
          end
        end
        temp2 = zeros(T_poiss(end),hmm.K);
        for k=1:hmm.K
          temp2(:,k)=sum(temp==k,1);
        end
        X_poiss = [X_poiss;temp2];
        Poiss_subj_inds = [Poiss_subj_inds;repmat(iSj,T_poiss(end),1)];
      else
        T_poiss(end)=[];
      end
    end
    T_poiss(T_poiss==0)=[];
    if whichstudy==4 || whichstudy==6 || whichstudy==7
      if ~isfolder([config.resultsdir,'Poissdata',int2str(W),overlapstring])
        mkdir([config.resultsdir,'Poissdata',int2str(W),overlapstring]);
      end
      mat_files_poiss{iSj} = [config.resultsdir,'Poissdata',int2str(W),overlapstring,'/PoissDataSj',int2str(iSj)];
      X = X_poiss;
      T = T_poiss;
      if exist(mat_files_poiss{iSj})
        error('careful - about to overwrite files containing STC info!')
      end
      save(mat_files_poiss{iSj},'X','T');
    elseif whichstudy==5

      mat_files_poiss{iSj} = strrep(mat_files_orth{iSj},'_orth',['_Poissdata',int2str(W)]);
      %mat_files_poiss{iSj} = [config.resultsdir,'Poissdata_',int2str(W),overlapstring,'/PoissDataSj',int2str(iSj)];
      X = X_poiss;
      T = T_poiss;

      if exist(mat_files_poiss{iSj})
        %error('careful - about to overwrite files containing STC info!')
      end
      save(mat_files_poiss{iSj},'X','T');

    end
  end


  if whichstudy==4 || whichstudy==6
    save([config.resultsdir,'Poissdata_',int2str(W),overlapstring,'/filelist.mat'],'mat_files_poiss');
  elseif whichstudy==5
    mkdir([config.resultsdir,'TASK_Poissdata_',int2str(W),overlapstring])
    save([config.resultsdir,'TASK_Poissdata_',int2str(W),overlapstring,'/filelist.mat'],'mat_files_poiss');
  end

else
  if whichstudy==4 || whichstudy==6
    load([config.resultsdir,'Poissdata_',int2str(W),overlapstring,'/filelist.mat'],'mat_files_poiss');
  elseif whichstudy==5
    load([config.resultsdir,'TASK_Poissdata_',int2str(W),overlapstring,'/filelist.mat'],'mat_files_poiss');
  end

end


%% We want to construct a state timecourse in quadrants of the cycle:

% normalise weights by FO:
origmetrics = load(config.metricfile);
FO_mean = mean(origmetrics.hmm_1stlevel.FO,1);
K2=4; % 4 state meta HMM
if whichstudy==1
  n_runs = 50; % do permutations for invalid circles
else
  n_runs=4;
end

if whichstudy~=5 %do not fit a new hmm for HCP task data, use the resting state one
  for i_run = 1:n_runs
    options = [];
    options.K = K2; % Note this didn't work with 4 states; one is knocked out and results in wierd behaviour
    options.distribution = 'poisson';
    options.Pstructure = eye(options.K) + diag(ones(1,options.K-1),1);
    options.Pstructure(options.K,1) = 1;
    %options.initrep = 4; % this the number of parallel cores
    options.useParallel = false;

    % load optimal sequence file:
    if strcmp(config.reordering_states, 'coherence')
      optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_coherence' ,'.mat'];
    elseif strcmp(config.reordering_states, 'study1matched')
      optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_study1matched' ,'.mat'];
    else
      optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'.mat'];
    end
    load(optimalseqfile);
    bestseq = bestsequencemetrics{1};

    % this is where we can shuffle entries for null test:
    if i_run<5
      bestseq = circshift(bestseq,i_run-1);
    else
      inds = randperm(12);
      bestseq = bestseq(inds);
    end

    disttoplot_manual = zeros(12,2);
    for i=1:12
      temp = exp(sqrt(-1)*(i+2)/12*2*pi);
      disttoplot_manual(bestseq(i),:) = [real(temp),imag(temp)];
    end

    circleposition = disttoplot_manual(:,1) + sqrt(-1)*disttoplot_manual(:,2);
    metastateposition = (2^-0.5)*exp(sqrt(-1)*(pi/2-[0:(K2-1)]*2*pi/K2));

    for i1=1:12
      for i2=1:K2
        FOweighting(i2,i1) = real(circleposition(i1))*real(metastateposition(i2)) + ...
          imag(circleposition(i1))*imag(metastateposition(i2));
      end
    end


    FOweighting = 1 + FOweighting;
    FO_metastate = repmat(FO_mean,K2,1).*FOweighting;
    FO_metastate = FO_metastate ./ repmat(sum(FO_metastate,2),1,12);

    % load initialised hmm:
    load([config.resultsdir, 'hmm_flat.mat'])
    if K2==4
      hmm_flat.state(K2) = hmm_flat.state(3);
      hmm_flat.K = K2;
      tmp = 0.01*ones(K2); tmp(find(eye(K2)))=1;
      hmm_flat.P = options.Pstructure.*(tmp);
      hmm_flat.Pi = ones(1,K2)./K2;
      hmm_flat.Dir_alpha = [hmm_flat.Dir_alpha mean(hmm_flat.Dir_alpha)];
      hmm_flat.Dir2d_alpha = hmm_flat.P;
      hmm_flat.prior.Dir2d_alpha = hmm_flat.Dir2d_alpha;
      hmm_flat.prior.Dir_alpha = ones(1,K2);
      hmm_flat.train.active = ones(1,K2);
      hmm_flat.train.K = K2;
      hmm_flat.train.Pstructure = options.Pstructure;
    end
    for k=1:K2
      hmm_flat.state(k).W.W_shape = W*hmm_flat.state(k).W.W_rate*FO_metastate(k,:);
      hmm_flat.state(k).W.W_mean = hmm_flat.state(k).W.W_shape./hmm_flat.state(k).W.W_rate;
    end

    options.hmm = hmm_flat;

    options.decodeGamma = false;
    options.standardise = false;
    options.cyc = 1; % we only allow this to run a single iteration - not to converge to global free energy minima
    if whichstudy~=4 && whichstudy~=6
      [hmmtemp,Gammatemp,~,~,~,~,fehist] = hmmmar(X_poiss,T_poiss,options);

      if i_run==1 || fehist(end)<lowestfe
        hmmPoiss = hmmtemp;
        GammaPoiss = Gammatemp;
        lowestfe = fehist(end);
      end
      % record a few stats to get a sense of how deviant these really are:
      feall(i_run,1) = fehist(1);
      feall(i_run,2) = fehist(end);
      gamsum(i_run,:) = mean(Gammatemp);
    else

      options.cyc=1; % limit the number of inference cycles to run on each noisy update
      % Hack to implemnet stochastic inference:
      n_batch = 100;
      for i=1
        thisbatch = sort(randperm(config.nSj,n_batch));
        X_poiss = []; T_poiss = [];
        % loading batch files:
        fprintf(['\nLoading batch files, batch',int2str(i)]);
        for j=1:n_batch
          temp = load(mat_files_poiss{thisbatch(j)});
          X_poiss = [X_poiss;temp.X];
          T_poiss = [T_poiss;temp.T];
        end
        fprintf(['\nFitting hmm to new batch:']);
        [hmmtemp,Gammatemp,~,~,~,~,fehist] = hmmmar(X_poiss,T_poiss,options);
  allhmm.hmm{i_run} = hmmtemp;
  allhmm.gamma{i_run} = Gammatemp;
        if i_run==1 || fehist(end)<lowestfe
          hmmPoiss = hmmtemp;
          GammaPoiss = Gammatemp;
          lowestfe = fehist(end);
        end
        feall(i_run,1) = fehist(1);
        feall(i_run,2) = fehist(end);
        gamsum(i_run,:) = mean(Gammatemp);
      end
    end
  end
end

%% plot free energy of valid vs invalid, ie null dist of models:

if whichstudy==1
fig=setup_figure([],1,1);
labels = [repmat({'valid'},4,1);repmat({'invalid'},46,1)];
[p,anovatab,stats] = anova1(feall(:,2),labels);
plot4paper('Model','Free Energy')
save([config.resultsdir, '2ndlevel_HMM_vs_permutations'], 'p', 'anovatab', 'stats','feall', 'gamsum','labels')
end

%% plot metastate profile:
ttl{1} = 'MEG UK';
ttl{3} = 'HCP';
ttl{6} = 'Cam-CAN';
ttl{7} = 'WakeHen';

FO_2ndlevel = mean(GammaPoiss);
%FO_2ndlevel = mean(Gammatemp);
statedist_all = zeros(1,12);
for i=1:K2
  statedist_all = statedist_all + FO_2ndlevel(i)*hmmPoiss.state(i).W.W_mean;
end
statedist_all = statedist_all ./ 125; % ratio rather than integer
for i=1:K2
  colorweights(:,i) = hmmPoiss.state(i).W.W_mean ./ 125 - statedist_all;
end
colorweights = (colorweights-min(colorweights(:))+0.01);

colorweights = log(colorweights) - min(log(colorweights(:))) + 0.01;
colorweights = colorweights./(max(colorweights(:)+0.1));


fig=setup_figure([], 1.5, 0.35);
CM = colormap(fig,hot(100));
load(optimalseqfile);
bestseq = bestsequencemetrics{1};

for i=1:K2
  if K2==4
    ax(i) = axes('Position', [0.035+(i-1)*0.225 0.1 0.175, 0.8]);
  else
    ax(i) = axes('Position', [0.035+(i-1)*0.3 0.1 0.225, 0.8]);
  end
  for i2=1:12
    CW{i2} = CM(ceil(length(CM)*colorweights(i2,i)),:);
  end
  cyclicalstatesubplot(bestseq,zeros(12),zeros(12),CW);
  if i==2
    title({ttl{whichstudy}, ''})
  end
end
ax(5) = axes('Position', [.9 0.175 0.1, .6]); axis off
h=colorbar;
% save_figure([config.figdir,'figure4_correlations/4supp_metastate_profile'],false);


%% Infer state path for each subject
if whichstudy<4
  hmmPoiss.gamma = GammaPoiss;
  disttoplot = plotMDS_states(hmmPoiss);

  save([config.resultsdir,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(options.K),overlapstring,'.mat'],'hmmPoiss','feall','GammaPoiss','T_poiss','Poiss_subj_inds');

  % plot convergence details;
  figure();
  subplot(3,1,1);
  plot(feall,'LineWidth',2);
  plot4paper('Init run','Free energy');
  legend('On init','On convergence')
  subplot(3,1,2);
  plot(std(gamsum,[],2),'*');
  plot4paper('Init run','Gamma std');
  subplot(3,1,3);
  plot(min(gamsum,[],2),'*');
  plot4paper('Init run','Min Gamma');
  save_figure([config.resultsdir,'figure/hmm_2ndlevel_ConvergenceRecord_window',int2str(W)],false);
else
  if whichstudy==4 || whichstudy==6 || whichstudy==7
    save([config.resultsdir,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(options.K),overlapstring,'.mat'],'hmmPoiss');
  else
    options = [];
    options.K = K2; % Note this didn't work with 4 states; one is knocked out and results in wierd behaviour
    options.distribution = 'poisson';
    options.Pstructure = eye(options.K) + diag(ones(1,options.K-1),1);
    options.Pstructure(options.K,1) = 1;
    options.useParallel = false;
    options.standardise = false;
    load([config.resultsdir,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(options.K),overlapstring,'.mat'],'hmmPoiss');

  end

  % and infer each subject's associated state timecourse:
  options.updateObs = 0;
  options.decodeGamma = 1;
  options.cyc = 1;
  options.hmm = hmmPoiss;

  mkdir([config.resultsdir,'Poissdata_',int2str(W),overlapstring])
  for i=1:length(mat_files_poiss)
    fprintf(['\n inferring STC for sj ',int2str(i)]);
    temp = load(mat_files_poiss{i});
    [~,Gamma] = hmmmar(temp.X,temp.T,options);
    mat_files_poiss{i} = [config.resultsdir,'Poissdata',int2str(W),overlapstring, '/PoissDataSj',int2str(i)];
    save(mat_files_poiss{i},'Gamma', '-append');
  end
  try
    save([config.resultsdir,'Poissdata',int2str(W),overlapstring, '/filelist.mat'],'mat_files_poiss', '-append')
  catch
    save([config.resultsdir,'Poissdata',int2str(W),overlapstring, '/filelist.mat'],'mat_files_poiss')
  end
end

%% finally, save summary stats to file for later analysis:
K=K2;
overlapstring='_overlappingWindows';

% infer cycle times:
if whichstudy<4
  load([config.resultsdir,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(K),overlapstring,'.mat'],'hmmPoiss','feall','GammaPoiss','T_poiss','Poiss_subj_inds');

  samp_minute = (config.sample_rate*60); % split into minute by minute chunks
  for subnum = 1:config.nSj

    Gamtemp = hmmPoiss.gamma(Poiss_subj_inds==subnum,:);
    if whichstudy==3
      cycletimes = getStateCycleTimes(Gamtemp,T_poiss([1:3]+(subnum-1)*3))./config.sample_rate;
      lifetimes_meta{subnum} = getStateLifeTimes(Gamtemp,T_poiss([1:3]+(subnum-1)*3));
    else
      cycletimes = getStateCycleTimes(Gamtemp,T_poiss(subnum))./config.sample_rate;
      lifetimes_meta{subnum} = getStateLifeTimes(Gamtemp,T_poiss(subnum));
    end
    cycletime_mu(subnum,:) = mean(cycletimes);
    cycletime_std(subnum,:) = std(cycletimes);
    cycletime_med(subnum,:) = median(cycletimes);
    FO_meta(subnum,:) = getFractionalOccupancy(Gamtemp,length(Gamtemp));
    cyctimes{subnum} = cycletimes;

    % how many microstates are typically traversed in each cycle?
    vpTemp = vpath{subnum};
    if length(T_poiss)~=config.nSj
      error('Need to realign state timecourses for sub segments')
    else
      vpTemp = vpTemp(ceil(W/2):end-floor(W/2));
      if length(vpTemp)~=length(Gamtemp)
        error('Sizes misaligned')
      end
    end
    if whichstudy==3
      num_transitions = getMicroStateCycleVisits(Gamtemp,T_poiss([1:3]+(subnum-1)*3),vpTemp);
    else
      [num_transitions,microsequences{subnum}] = getMicroStateCycleVisits(Gamtemp,T_poiss(subnum),vpTemp);
    end
    microtransitions_all{subnum} = num_transitions;
    microtransitions_mu(subnum,:) = mean(num_transitions);
    microtransitions_median(subnum,:) = median(num_transitions);

    if whichstudy==3 % get session by session detail:
      for subnum = 1:config.nSj
        for isession=1:3
          t0 = 1+sum(T_poiss(1:((subnum-1)*3)+(isession-1)));tend = sum(T_poiss(1:((subnum-1)*3)+(isession)));
          Gamtemp = hmmPoiss.gamma(t0:tend,:);
          cycletimes = getStateCycleTimes(Gamtemp,tend-t0+1)./config.sample_rate;
          cycletime_mu_sess(subnum,isession) = mean(cycletimes);
          cycletime_std_sess(subnum,isession) = std(cycletimes);
          cycletime_med_sess(subnum,isession) = median(cycletimes);
          FO_meta_sess(subnum,isession,:) = getFractionalOccupancy(Gamtemp,length(Gamtemp));
        end
      end
    end
  end
else
  load([config.resultsdir,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(options.K),overlapstring,'.mat'],'hmmPoiss');
  if isfield(config, 'Poiss_dir') && ~contains(config.Poiss_dir,'overlappingWindows')
    clear FO
    for subnum=1:config.nSj
      Gamtemp = hmmPoiss.gamma(Poiss_subj_inds==subnum,:);
      cycletimes = getStateIntervalTimes(Gamtemp,length(Gamtemp))./config.sample_rate;
      cycletime_mu(subnum,:) = cellfun(@mean,cycletimes);
      cycletime_std(subnum,:) = cellfun(@std,cycletimes);
      cycletime_med(subnum,:) = cellfun(@median,cycletimes);
      FO(subnum,:) = getFractionalOccupancy(Gamtemp,length(Gamtemp));
      cyctimes{subnum} = cycletimes;
      lifetimes_meta{subnum} = getStateLifeTimes(Gamtemp,length(Gamtemp));
      for ik=1:length(lifetimes_meta{subnum})
        lifetimes_meta{subnum}{ik} = lifetimes_meta{subnum}{ik}./config.sample_rate;
      end
      

      % how many microstates are typically traversed in each cycle?
      vpTemp = load(mat_files_orth{subnum},'vpath','T');
      vpTemp = vpTemp.vpath;
      if length(T_poiss)~=config.nSj
        error('Need to realign state timecourses for sub segments')
      else
        vpTemp = vpTemp(ceil(W/2):end-floor(W/2));
        if length(vpTemp)~=length(Gamtemp)
          error('Sizes misaligned')
        end
      end
      if whichstudy==3
        num_transitions = getMicroStateCycleVisits(Gamtemp,T_poiss([1:3]+(subnum-1)*3),vpTemp);
      else
        [num_transitions,microsequences{subnum}] = getMicroStateCycleVisits(Gamtemp,T_poiss(subnum),vpTemp);
      end
      microtransitions_all{subnum} = num_transitions;
      microtransitions_mu(subnum,:) = mean(num_transitions);
      microtransitions_median(subnum,:) = median(num_transitions);
    end  

  else
    figdir = [config.figdir,'4_covariates_W',int2str(W),'_overlappingWindows/'];
    mkdir(figdir);
    clear cycletimes cycletime_mu cycletime_std cycletime_med FO cyctimes lifetimes cycletime_mu_min cycletime_med_min cycletime_std_min
    load([config.resultsdir,'Poissdata',int2str(W),overlapstring, '/filelist.mat'])
    samp_2minute = config.sample_rate*2*60;
    for subnum=1:config.nSj
      fprintf(['\nSubj: ',int2str(subnum)]);
      load(mat_files_poiss{subnum},'Gamma');

      T = length(Gamma);

      %cycletimes = getStateIntervalTimes(Gamma,T);
      cycletimes{1} = getStateCycleTimes(Gamma,T)./config.sample_rate;
      cycletime_mu(subnum,:) = cellfun(@mean,cycletimes);
      cycletime_std(subnum,:) = cellfun(@std,cycletimes);
      cycletime_med(subnum,:) = cellfun(@median,cycletimes);
      FO_meta(subnum,:) = getFractionalOccupancy(Gamma,length(Gamma));
      cyctimes{subnum} = cycletimes;
      lifetimes_meta{subnum} = getStateLifeTimes(Gamma,T);
      for ik=1:length(lifetimes_meta{subnum})
        lifetimes_meta{subnum}{ik} = lifetimes_meta{subnum}{ik}./config.sample_rate;
      end

      % how many microstates are typically traversed in each cycle?
      vpTemp = load(mat_files_orth{subnum},'vpath','T');

      % correct Gamma for points that were erronesouly enteres:
      temp2 = load(mat_files_poiss{subnum});
      to_remove = sum(temp2.X')==0;
      Gamtemp2 = Gamma(~to_remove,:);
      T_amended = vpTemp.T-14; %-W+1
      if length(vpTemp.vpath)~=sum(T_amended)
        error('Not aligned!!!')
      end
      vpTemp = vpTemp.vpath;
      t_start = ceil(W/2);
      t_end = -ceil(W/2)+1;
      vpTemp2 = [];
      for i=1:length(T_amended)
        t_end = t_end + T_amended(i);
        vpTemp2 = [vpTemp2; vpTemp(t_start:t_end)];
        t_start = t_start + T_amended(i);
      end
      T_amended = T_amended-W+1;
      if length(vpTemp2)~=sum(T_amended)
        error('Need to realign state timecourses for sub segments')
      end
      [num_transitions,microsequences{subnum}] = getMicroStateCycleVisits(Gamtemp2,T_amended,vpTemp2);

      microtransitions_all{subnum} = num_transitions;
      microtransitions_mu(subnum,:) = mean(num_transitions);
      microtransitions_median(subnum,:) = median(num_transitions);

      % also get minute by minute detail:
      for imin=1:5
        startseg = find(cumsum(T)>(imin-1)*samp_2minute,1)-1;
        endseg = find(cumsum(T)>(imin)*samp_2minute,1)-1;
        if ~isempty(startseg) && ~isempty(endseg)
          if startseg==endseg
            T_sub = samp_2minute;
          else
            T_sub = sum(T(1:startseg+1)) - (imin-1)*samp_2minute;
            T_sub = [T_sub;T(startseg+2:endseg)];
            T_sub = [T_sub;samp_2minute - sum(T_sub)];
          end
          try
            Gamtemp = Gamma((imin-1)*samp_2minute + [1:samp_2minute],:);
          catch
            Gamtemp = Gamma((imin-1)*samp_2minute :end,:);
            T_sub = T_sub(cumsum(T_sub)<length(Gamtemp));
            T_sub = [T_sub;length(Gamtemp)-sum(T_sub)];
          end
          temp = getStateCycleTimes(Gamtemp,T_sub)./config.sample_rate;
          cyctime{subnum,imin} = temp;
          cycletime_mu_min(subnum,imin) = mean(temp);
          cycletime_med_min(subnum,imin) = median(temp);
          cycletime_std_min(subnum,imin) = std(temp);
        end
      end
    end
  end
end

hmm_2ndlevel = [];
hmm_2ndlevel.cyctime_mu = cycletime_mu;
hmm_2ndlevel.cyctime_med = cycletime_med;
hmm_2ndlevel.cyctime_std = cycletime_std;
hmm_2ndlevel.FO = FO_meta;
hmm_2ndlevel.cyctimes_all = cyctimes;
hmm_2ndlevel.LT_all = lifetimes_meta;

hmm_2ndlevel.microtransitions_all = microtransitions_all;
hmm_2ndlevel.microtransitions_mu = microtransitions_mu;
hmm_2ndlevel.microtransitions_med = microtransitions_median;
hmm_2ndlevel.microsequences = microsequences;

if whichstudy==3
  hmm_2ndlevel.cycletime_mu_sess = cycletime_mu_sess;
  hmm_2ndlevel.cycletime_std_sess = cycletime_std_sess;
  hmm_2ndlevel.cycletime_med_sess = cycletime_med_sess;
  hmm_2ndlevel.FO_meta_sess = FO_meta_sess;
elseif whichstudy==4 || whichstudy==6
  hmm_2ndlevel.cyctime_min = cyctime;
  hmm_2ndlevel.cycletime_mu_min=cycletime_mu_min;
  hmm_2ndlevel.cycletime_med_min=cycletime_med_min;
  hmm_2ndlevel.cycletime_std_min=cycletime_std_min;
end

save(config.metricfile,'hmm_2ndlevel','-append');

%% run some second level metastate stats:

N = length(hmm_2ndlevel.cyctime_mu);
%plot(hmm_2ndlevel.cyctimes_all{randi(N)})
outliers = abs(1./hmm_2ndlevel.cyctime_mu(:,1) - mean(1./hmm_2ndlevel.cyctime_mu(:,1))) > 2*std(1./hmm_2ndlevel.cyctime_mu(:,1));

figure('Position',[588 326 274 471]);
distributionPlot(hmm_2ndlevel.microtransitions_mu(~outliers,:),'showMM',2,'color',{color_scheme{1:size(hmm_1stlevel.FO,2)}})
%title('Expected state visits per cycle');
plot4paper('','Expected Number of state visits per cycle');grid on;
set(gca,'XTickLabel',{'Total state visits','Unique state visits'})
set(gca,'XTickLabelRotation',45)
%set(gca,'YLim',[0 YL],'FontSize',fontsize,'FontSize',fontsize);
print([config.figdir,'Fig4A_StateVisitsPerCycle_window',int2str(W),'_K',int2str(K)],'-dpng');

figure('Position',[588 63 412 735]);
color_scheme{13} = [0 0 0]
%subplot(211);
distributionPlot([hmm_1stlevel.LT_mu(~outliers,:),hmm_2ndlevel.cyctime_mu(~outliers,:)] ./ config.sample_rate * 1000,'showMM',2,'color',{color_scheme{:}})
title('Life Times');plot4paper('','Time (ms)');grid on;
YL = 1.1*prctile(hmm_1stlevel.LT_mu(:),97.5)./ config.sample_rate * 1000;
for i=1:12;labels{i} = ['State',int2str(i)];end
labels{13} = 'CycleTime';
set(gca,'XTickLabel',labels)
print([config.figdir,'Fig4B_CycleTimesVsLifeTimes',int2str(W),'_K',int2str(K)],'-dpng');

%% how many microstates are typically visited in each cycle?

boxplot(microtransitions_median)

%% do consecutive traversals tend to visit the same states?
clear stats s

for subnum = 1:config.nSj
    fprintf(['\nRunning subj', int2str(subnum),'\n'])
    for T = 1:20
        for i_perm=1:100
            temp = hmm_2ndlevel.microsequences{subnum};
            N = length(temp);
            if i_perm==1
                ordering = 1:N;
            else
                ordering = randperm(N);
            end
            L_common = zeros(N-T,1);
            L_unique = zeros(N-T,1);
            for n = 1:(N-T)
                L_common(n) = length(intersect(temp{ordering(n)},temp{ordering(n+T)}));
                L_unique(n) = length(union(temp{ordering(n)},temp{ordering(n+T)}));
            end
            stats(subnum,i_perm,T) = mean(L_common ./ L_unique);
        end
    end
end

%% 
figure('Position',[588 63 412 735]);
temp = {stats(~outliers,1,1),squash(stats(~outliers,2:end,1))};
p = anova1([temp{1};temp{2}(:)],[ones(length(temp{1}),1);2*ones(length(temp{2}(:)),1)],'off')
distributionPlot(temp,'showMM',2,'color',{color_scheme{[1,3]}})
title(['p=',num2str(p,3)]);plot4paper('','Proportion common states visited');grid on;
set(gca,'XTickLabel',{'Consecutive cycles','Nonconsecutive cycles'})
print([config.figdir,'Fig4D_ConsecutivePattern1',int2str(W),'_K',int2str(K)],'-dpng');

%%
figure;
shadedErrorBar(1:20,squeeze(nanmean(stats(~outliers,1,:),1)),squeeze(nanstd(stats(~outliers,1,:),[],1)./sqrt(subnum)));
hold on;
for i_perm = 2:100
    s(i_perm-1,:) = squeeze(nanmean(stats(~outliers,i_perm,:),1));
end
plot(1:20,repmat(max(s(:)),1,20),'r--');
h(1) = plot(nan,nan,'k-');
h(2) = plot(1:20,repmat(min(s(:)),1,20),'r--');
plot4paper('Cycles N distance apart','Proportion common states visited')
legend(h,{'Observed statistic','p<5e-4'})
print([config.figdir,'Fig4D_ConsecutivePattern2',int2str(W),'_K',int2str(K)],'-dpng');

%%

if whichstudy==1
    [info,varnames] = MEGUK_getparticipantinfo();
    subj_age = info(:,11);
    subj_gender = info(:,12); % 1 denotes male, 2 denotes female
    subj_RTs = info(:,3);
elseif whichstudy==3
    [info,varnames] = HCP_getparticipantinfo(config);
    subj_age = info(:,4);
    subj_gender = info(:,3); % 1 denotes female
    subj_RTs = info(:,end);
elseif whichstudy==4
    info = camcan_getparticipantinfo(config);
    subj_age = info(:,1);
    subj_gender = info(:,3); % ??
    subj_RTs = info(:,5);
    
end
info = [subj_age,subj_gender,subj_RTs];

figure();
i=3
teststat = stats(~outliers,1,1) - max(stats(~outliers,2:end,1),[],2)
scatter(info(~outliers,i),teststat)
[C,P] = corr(info(~outliers,i),teststat)

%% histogram of state transitions in phase space:

bestseq = bestsequencemetrics{1}; % MATS: flagging you to check this bit please. Can't remember which of the optimals we chose in the end
    
disttoplot_manual = zeros(12,2);
for i=1:12
    temp = exp(sqrt(-1)*(i+2)/12*2*pi);
    disttoplot_manual(bestseq(i),:) = [real(temp),imag(temp)];
end
    %circleposition = exp(j*(pi/2-[0:11]*2*pi/12));
    %circleposition = circleposition(bestseq);
    circleposition = disttoplot_manual(:,1) + sqrt(-1)*disttoplot_manual(:,2);
phaseshiftcounter = zeros(config.nSj,12);
phases = linspace(-pi,5*pi/6,12);
jointcounter = zeros(12,12,config.nSj);
for subnum=1:config.nSj

    if whichstudy==4
        clear vpath;
        vpTemp = load(mat_files_orth{subnum},'vpath','T');
        vpath{subnum} = vpTemp.vpath;
    end
    transitions = find(diff(vpath{subnum})~=0);
    phaseshift = zeros(length(transitions),1);
    for i=1:length(transitions)
        phase_t = angle(circleposition(vpath{subnum}(transitions(i))));
        phase_tplus1 = angle(circleposition(vpath{subnum}(transitions(i)+1)));
        phaseshift(i) = phase_tplus1 - phase_t;
    end
    % unwrap phases onto same scale:
    phaseshift(phaseshift > pi-0.01) = phaseshift(phaseshift > pi-0.01)-2*pi;
    phaseshift(phaseshift < -pi-0.01) = phaseshift(phaseshift < -pi-0.01)+2*pi;
    phaseshiftcounter(subnum,:) = hist(phaseshift,phases);

    % track phase shifts for consecutive transitions:
    for i=1:12
        pointsin = phaseshift(1:end-1) >= phases(i) - 0.1 & phaseshift(1:end-1) <= phases(i) + 0.1;
        pointsin = find(pointsin);
        for i2 = 1:12
            jointcounter(i,i2,subnum) = sum(phaseshift(pointsin + 1) >= phases(i2) - 0.1 & phaseshift(pointsin + 1) <= phases(i2) + 0.1);
        end
    end
end
phaseshiftcounter = phaseshiftcounter ./ repmat(sum(phaseshiftcounter,2),1,12);
figure();
bar(phases,mean(phaseshiftcounter))
plot4paper('Phase Shift','Proportion of microstate transitions')
labs = {'-\pi','-5\pi/6','-2\pi/3','-\pi/2','-\pi/3','-\pi/6','0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6'}
set(gca,'XTick',phases,'XTickLabel',labs)
title('Histogram of Microstate Transition Phase Shifts')
print([config.figdir,'Fig4C_PhaseShifts',int2str(W),'_K',int2str(K)],'-dpng');

%% show effect of age as well:
%B = pinv([ones(config.nSj,1),subj_age])*phaseshiftcounter
for i=1:12;
    phaselabels{i} = num2str(phases(i),2);
end
%tbl = table([phaseshiftcounter,subj_age],'VariableNames',{phaselabels{:},'Age'})
clear B
for i=1:12
    tbl = table(phaseshiftcounter(:,i),subj_age,'VariableNames',{'Phase','Age'})
    mdl = fitlm(tbl,'Phase ~ Age')
    pvals(i) = mdl.Coefficients.pValue(2);
    B(i) = mdl.Coefficients.Estimate(2);     
end
figure('Position',[200 324 815 376]);
subplot(121);
for i=1:12
    if pvals(i)<0.05/12
        col = 'green';
    else
        col = 'blue';
    end
    bar(phases(i),B(i),0.8*pi/6,'green');
    hold on
end
h(1) = bar(nan,nan,'green');
h(2) = bar(nan,nan,'blue');
legend(h,{'p<5e-3','N.S.'})
plot4paper('Phase Shift','Age Cofficient')
labs = {'-\pi','-5\pi/6','-2\pi/3','-\pi/2','-\pi/3','-\pi/6','0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6'}
set(gca,'XTick',phases,'XTickLabel',labs)
%title('Histogram of Microstate Transition Phase Shifts')
% fit multinomial model - can we predict age:
tbl = table(phaseshiftcounter(:,1),subj_age,'VariableNames',{'Phase1','Age'})
for i=2:12
    name = ['Phase',num2str(i)];
    tbl.(name) = phaseshiftcounter(:,i)
end
mdl = fitlm(tbl,'Age ~ Phase1 + Phase2 + Phase3 + Phase4 + Phase5 + Phase6 + Phase7 + Phase8 + Phase9 + Phase10 + Phase11 + Phase12');

subplot(122)
scatter(subj_age,mdl.predict,'filled')
XL = xlim();hold on;
plot(XL,XL,'Color',[1 1 1]*0.8,'LineWidth',0.8)
plot4paper('Actual Age','Predicted Age')
title(['p=6.2e-15'])
var_expl = [var(subj_age) - var(mdl.predict-subj_age)] ./ var(subj_age)
print([config.figdir,'Fig4C_PhaseShifts_Age_',int2str(W),'_K',int2str(K)],'-dpng');

%% Gender:

clear B
for i=1:12
    tbl = table(phaseshiftcounter(:,i),normalize(subj_gender),'VariableNames',{'Phase','Gender'})
    mdl = fitlm(tbl,'Phase ~ Gender')
    pvals(i) = mdl.Coefficients.pValue(2);
    B(i) = mdl.Coefficients.Estimate(2);     
end
figure('Position',[200 324 815 376]);
subplot(121);
bar(phases,B)
plot4paper('Phase Shift','Gender effect')

set(gca,'XTick',phases,'XTickLabel',labs)
%title('Histogram of Microstate Transition Phase Shifts')
% fit multinomial model - can we predict age:
tbl = table(phaseshiftcounter(:,1),normalize(subj_gender),'VariableNames',{'Phase1','Gender'})
for i=2:12
    name = ['Phase',num2str(i)];
    tbl.(name) = phaseshiftcounter(:,i)
end
mdl = fitlm(tbl,'Gender ~ Phase1 + Phase2 + Phase3 + Phase4 + Phase5 + Phase6 + Phase7 + Phase8 + Phase9 + Phase10 + Phase11 + Phase12');

subplot(122)
scatter(normalize(subj_gender),mdl.predict,'filled')
plot4paper('Gender','Predicted Age')

%% Histogram phase shifts t and t+1 timesteps apart:
figure();
imagesc(sum(jointcounter,3))

%%
[A,B] = pca(phaseshiftcounter)

figure()
for i=1:4
    subplot(1,4,i)
    bar(A(:,i))
end

[C,P] = corr(B(:,2),subj_age)

