function [FO,pvals,tintervalsout] = computeLongTermAsymmetry(vpath,T,K,intervalpercentiles)
% computes a single subject's interval FO assymettry matrix. Note will
% return NaN if there are no intervals
if ~iscolumn(vpath)
  vpath = vpath';
  if ~iscolumn(vpath)
    error('Viterbi path must be passed in as column vector');
  end
end
runonquartiles = false; % implement this later to account for sub-elements of intervals (ie not just first half vs second half)
if nargin<4
  intervalpercentiles = [0,100];
end
if ~iscell(vpath) || ~iscell(T)
  error('vpath and T should be cells, with one entry per subject');
end
nSj = length(vpath);
FO = zeros(K,K,2,nSj,length(intervalpercentiles)-1);
for iSj=1:nSj
  if sum(T{iSj})~=length(vpath{iSj})
    error('Dimension mismatch between T and vpath');
  end
  T_breaks = cumsum(T{iSj}); % indices of HMM chain breaks to be excluded
  for ik=1:K
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
            else
              tempaway = [tempaway;vpath{iSj}(offtimes(t):offtimes(t)+floor(t_interval/4),:)];
              tempto = [tempto;vpath{iSj}(ontimes(t+1)-floor(t_interval/4):ontimes(t+1),:)];
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
    end
    %         end
  end
end

% paired t-tests to evaluate significance:
pvals = zeros(K,K,length(intervalpercentiles)-1);
for ip=1:length(intervalpercentiles)-1
  for ik1=1:12
    for ik2=1:12
      [~,pvals(ik1,ik2,ip)] = ttest(squeeze(FO(ik1,ik2,1,:,ip)-FO(ik1,ik2,2,:,ip)));
    end
  end
end
end