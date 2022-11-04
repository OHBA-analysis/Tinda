function [FO,pvals,tintervalsout, stat, FO_residual,pvals_residual, stat_res] = computeLongTermAsymmetry(vpath,T,K,intervalpercentiles,dotaskglm)
% computes a single subject's interval FO assymettry matrix. Note will
% return NaN if there are no intervals
if ~iscolumn(vpath)
  vpath = vpath';
  if ~iscolumn(vpath)
    error('Viterbi path must be passed in as column vector');
  end
end
runonquartiles = false; % implement this later to account for sub-elements of intervals (ie not just first half vs second half)
if nargin<4 || isempty(intervalpercentiles)
  intervalpercentiles = [0,100];
end
if nargin<5
  dotaskglm = false;
end
if ~iscell(vpath) || ~iscell(T)
  error('vpath and T should be cells, with one entry per subject');
end
nSj = length(vpath);
FO = zeros(K,K,2,nSj,length(intervalpercentiles)-1);
FO_all = [];
if dotaskglm, FO_residual = zeros(K,K,2,nSj,length(intervalpercentiles)-1);
else, FO_residual=[]; end

for iSj=1:nSj
  if sum(T{iSj})~=length(vpath{iSj})
    error('Dimension mismatch between T and vpath');
  end
  T_breaks = cumsum(T{iSj}); % indices of HMM chain breaks to be excluded
  if dotaskglm
    % get the demeaned vpath data to use in the GLM later
    v=zeros(length(vpath{iSj}), K);
    for ik=1:K
      v(:,ik) = vpath{iSj}==ik;
    end
    % reshape v to get trials
    v=reshape(v,[], length(T{iSj}), K);
    v_demean = demean(v, 2);
    v_demean = reshape(v_demean, [], K);
  end
  for ik=1:K
    intervalDesign = zeros(length(vpath{iSj}),1);
    tempaway = [];
    tempto = [];
    ontimes = find(diff([0;vpath{iSj}==ik;0])==1)-1;
    offtimes = find(diff([0;vpath{iSj}==ik;0])==-1);
    %         if isempty(percentiles)
    %             for t=1:length(offtimes)-1
    %                 % assert interval doesn't cross segment boundaries:
    %                 if find(ontimes(t+1)<T_breaks,1) == find(offtimes(t)<T_breaks,1)
    %                     t_interval = ontimes(t+1)-offtimes(t);
    %                     tempaway = [tempaway;vpath(offtimes(t):offtimes(t)-1+ceil(t_interval/2),:)];
    %                     tempto = [tempto;vpath(ontimes(t+1)+1-ceil(t_interval/2):ontimes(t+1),:)];
    %                     if nargout>1
    %                         tintervals{ik}(t) = t_interval+1;
    %                     end
    %                 end
    %             end
    %             for ik2=1:K
    %                 FO(ik,ik2,1) = mean(tempaway==ik2);
    %                 FO(ik,ik2,2) = mean(tempto==ik2);
    %             end
    %         else
    tintervals{iSj,ik} = ontimes(2:end)-offtimes(1:end-1);
    for ip=1:length(intervalpercentiles)-1
      if ~any(isnan(intervalpercentiles(ip:ip+1)))
        p_low = prctile(tintervals{iSj,ik},intervalpercentiles(ip));
        p_high = prctile(tintervals{iSj,ik},intervalpercentiles(ip+1));
      elseif isnan(intervalpercentiles(ip))
        p_low = 0;
        p_high = intervalpercentiles(ip+1);
      elseif isnan(intervalpercentiles(ip+1))
        p_low = intervalpercentiles(1);
        p_high = max(tintervals{iSj,ik})+1;
      end
      tempaway = [];
      tempto = [];
      for t=1:length(offtimes)-1
        t_interval = ontimes(t+1)-offtimes(t);
        if ~any(ismember(T_breaks,[offtimes(t):ontimes(t+1)]))
          if t_interval>=p_low && t_interval<=p_high
            if ~runonquartiles
              tempaway = [tempaway;vpath{iSj}(offtimes(t):offtimes(t)+floor(t_interval/2),:)];
              tempto = [tempto;vpath{iSj}(ontimes(t+1)-floor(t_interval/2):ontimes(t+1),:)];
              intervalDesign(offtimes(t):offtimes(t)+floor(t_interval/2), 1) = 1;
              intervalDesign(ontimes(t+1)-floor(t_interval/2):ontimes(t+1), 1) = 2;
            else
              tempaway = [tempaway;vpath{iSj}(offtimes(t):offtimes(t)+floor(t_interval/4),:)];
              tempto = [tempto;vpath{iSj}(ontimes(t+1)-floor(t_interval/4):ontimes(t+1),:)];
              intervalDesign(offtimes(t):offtimes(t)+floor(t_interval/4), 1) = 1;
              intervalDesign(ontimes(t+1)-floor(t_interval/4):ontimes(t+1), 1) = 2;
            end
            %intervaltracker{ip,ik}(length(intervaltracker{ip,ik})+1) = t_interval;
          end
        end
      end
      for ik2=1:K
        if ~isempty(tempaway)
          FO(ik,ik2,1,iSj,ip) = mean(tempaway==ik2);
          FO(ik,ik2,2,iSj,ip) = mean(tempto==ik2);
        end
      end
      tintervalsout{iSj,ik,ip} = tintervals{iSj,ik}(tintervals{iSj,ik}>=p_low & tintervals{iSj,ik}<=p_high);
      
      if dotaskglm
        des = zeros(length(v_demean), 2);
        des(:,1) = intervalDesign==1;
        des(:,2) = intervalDesign==2;
        for ik2=setdiff(1:K,ik)
          b = glmfit(des, v_demean(:,ik2));
          FO_residual(ik,ik2,:,iSj,ip) = b(2:3);
        end
      end
    end
    %         end
  end
end

% Do Permutation test on FO (instead of t-test)

dat1=[];
dat1.dimord = 'rpt_chan_time';
dat1.label{1} = 'FO asym';
dat1.time=1:(K^2-K);
dat2=dat1;
cfg=[];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.design = [ones(1,nSj), 2*ones(1,nSj); 1:nSj, 1:nSj];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 100000;

pvals = zeros(K,K,length(intervalpercentiles)-1);
if dotaskglm
  pvals_residual = zeros(K,K,length(intervalpercentiles)-1);
else
  pvals_residual=[]; stat_res=[]; 
end

for ip=1:length(intervalpercentiles)-1
  % do permutation test of FO
  
  tmp = permute(squeeze(FO(:,:,1,:,ip)), [3,1,2]);
  dat1.trial(:,1,:) = tmp(:, ~eye(K));
  
  tmp = permute(squeeze(FO(:,:,2,:,ip)), [3,1,2]);
  dat2.trial(:,1,:) = tmp(:, ~eye(K));
  
  stat{ip} = ft_timelockstatistics(cfg, dat1, dat2);
  tmp = ones(K);
  tmp(~eye(K)) = stat{ip}.prob;
  pvals(:,:,ip) = tmp;
  
  if dotaskglm
    tmp = permute(squeeze(FO_residual(:,:,1,:,ip)), [3,1,2]);
    dat1.trial(:,1,:) = tmp(:, ~eye(K));
    
    tmp = permute(squeeze(FO_residual(:,:,2,:,ip)), [3,1,2]);
    dat2.trial(:,1,:) = tmp(:, ~eye(K));
    
    stat_res{ip} = ft_timelockstatistics(cfg, dat1, dat2);
    tmp = ones(K);
    tmp(~eye(K)) = stat_res{ip}.prob;
    pvals_residual(:,:,ip) = tmp;
  end
end
if length(stat)==1
  stat = stat{1};
end
if length(stat_res)==1
  stat_res = stat_res{1};
end