function [FO,pvals,tintervalsout, stat, FO_residual,pvals_residual, stat_res] = computeLongTermAsymmetry(vpath,T,K,intervalpercentiles,dotaskglm, replayidx)
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
    if exist('replayidx') && ~isempty(replayidx)
      v=zeros(length(vpath{iSj}), K);
      for ik=1:K
        v(:,ik) = vpath{iSj}==ik;
      end
      v_mean = zeros(251, length(replayidx),K);
      rem=[];
      for k=1:length(replayidx)
        ix = replayidx(k);
        try
          v_mean(:, k,:) = v(ix-125:ix+125,:);
        catch
          rem = [rem, k];
        end
      end
      v_mean(:, rem,:) = [];
      v_mean = squeeze(mean(v_mean,2));
      replayidx(rem)=[];
      v_demean = v;
      for k=1:length(replayidx)
        ix = replayidx(k);
        v_demean(ix-125:ix+125,:) = v(ix-125:ix+125,:)-v_mean;
      end
      
      
    else
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

pvals = zeros(K,K,length(intervalpercentiles)-1);
stat=[];
if dotaskglm
  pvals_residual = zeros(K,K,length(intervalpercentiles)-1);
else
  pvals_residual=[]; stat_res=[];
end
if nSj>1
  for ip=1:length(intervalpercentiles)-1
    % do permutation test of FO
    
    [pvals(:,:,ip), stat{ip}] = FO_permutation_test(FO(:,:,:,:,ip), K, nSj);
    
    if dotaskglm
      [pvals_residual(:,:,ip), stat_res{ip}] = FO_permutation_test(FO_residual(:,:,:,:,ip), K, nSj);
    end
  end
  if length(stat)==1
    stat = stat{1};
  end
  if length(stat_res)==1
    stat_res = stat_res{1};
  end
else
  warning('The FO asymmetry cannot not be statistically tested because there is only one subject')
end