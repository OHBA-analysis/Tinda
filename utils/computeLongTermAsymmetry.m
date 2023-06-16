function [FO,pvals,tintervalsout, stat, FO_residual,pvals_residual, stat_res] = computeLongTermAsymmetry(vpath,T,K,intervalbins, dotaskglm, replayidx, dostat, ntiles, intervalbin_mode)
% computes a single subject's interval FO assymettry matrix. Note will
% return NaN if there are no intervals
if ~iscolumn(vpath)
    vpath = vpath';
    if ~iscolumn(vpath)
        error('Viterbi path must be passed in as column vector');
    end
end
if ~exist('ntiles', 'var') || isempty(ntiles)
    ntiles = 2; % implement this later to account for sub-elements of intervals (ie not just first half vs second half)
end
if nargin<4 || isempty(intervalbins)
    intervalbins = [0,100];
end
if ~exist('intervalbin_mode', 'var') || isempty(intervalbin_mode)
    intervalbin_mode = 'percentile'; % can be 'abs', 'percentile'
end
if ~exist('dotaskglm', 'var') || isempty(dotaskglm)
    dotaskglm = false;
end
if ~iscell(vpath) || ~iscell(T)
    error('vpath and T should be cells, with one entry per subject');
end
if ~exist('numperm', 'var') || isempty(numperm)
    numperm = 10^5;
end

if ~exist('dostat', 'var') || isempty(dostat)
    dostat=1;
end

if exist('replayidx', 'var') && ~isempty(replayidx)
    replay_state=1;
else
    replay_state=0;
end    


nSj = length(vpath);

if ntiles>2
    FO = zeros(K+replay_state,K+replay_state,2,nSj,length(intervalbins)-1, ceil(ntiles/2));
    dostat=false;
else
    FO = zeros(K+replay_state,K+replay_state,2,nSj,length(intervalbins)-1);
end
FO_all = [];
if dotaskglm, FO_residual = zeros(K+replay_state,K+replay_state,2,nSj,length(intervalbins)-1);
else, FO_residual=[]; end

for iSj=1:nSj
    if replay_state
        replay = zeros(size(vpath{iSj}));
        replay(replayidx{iSj})=1;
    end
    if sum(T{iSj})~=length(vpath{iSj})
        error('Dimension mismatch between T and vpath');
    end
    T_breaks = cumsum(T{iSj}); % indices of HMM chain breaks to be excluded
    if dotaskglm
        if exist('replayidx', 'var') && ~isempty(replayidx)
            v=zeros(length(vpath{iSj}), K);
            for ik=1:K
                v(:,ik) = vpath{iSj}==ik;
            end
            v_mean = zeros(251, length(replayidx{iSj}),K);
            rem=[];
            for k=1:length(replayidx{iSj})
                ix = replayidx{iSj}(k);
                try
                    v_mean(:, k,:) = v(ix-125:ix+125,:);
                catch
                    rem = [rem, k];
                end
            end
            v_mean(:, rem,:) = [];
            v_mean = squeeze(mean(v_mean,2));
            replayidx{iSj}(rem)=[];
            v_demean = v;
            for k=1:length(replayidx{iSj})
                ix = replayidx{iSj}(k);
                v_demean(ix-125:ix+125,:) = v(ix-125:ix+125,:)-v_mean;
            end
            
            
        else % this will be used when we have a trial x time matrix.
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
    for ik=1:K+replay_state
        intervalDesign = zeros(length(vpath{iSj}),1);
        tempaway = [];
        tempto = [];
        
        if ik>K % this means there is a replay state
            ontimes = find(diff([0;replay;0])==1)-1;
            offtimes = find(diff([0;replay;0])==-1);
        else
            ontimes = find(diff([0;vpath{iSj}==ik;0])==1)-1;
            offtimes = find(diff([0;vpath{iSj}==ik;0])==-1);
        end
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
        for ip=1:length(intervalbins)-1
            if ~any(isnan(intervalbins(ip:ip+1)))
                if strcmp(intervalbin_mode, 'abs')
                    p_low = intervalbins(ip);
                    p_high = intervalbins(ip+1)-1; % this forces it to be smaller than this value
                else
                    
                if 0 % looks at percentages without zero interval lengths
                    tmp = tintervals{iSj,ik};
                    tmp(tmp==0)=[];
                    p_low = prctile(tmp,intervalbins(ip));
                    p_high = prctile(tmp,intervalbins(ip+1));
                else
                    p_low = prctile(tintervals{iSj,ik},intervalbins(ip));
                    p_high = prctile(tintervals{iSj,ik},intervalbins(ip+1));
                end
                end
            elseif isnan(intervalbins(ip))
                p_low = 0;
                p_high = intervalbins(ip+1);
            elseif isnan(intervalbins(ip+1))
                p_low = intervalbins(1);
                p_high = max(tintervals{iSj,ik})+1;
            end
            tempaway = [];
            tempto = [];
            for ntile=1:ntiles/2 % this means we can now seperate each interval into however many chunks we want. 
                % if ntiles>2, FO(ik1, ik2, 1, iSj, ip, 1:ntile/2) will contain
                % 0:ntile:ntiles/2 and FO(ik1, ik2, 2, iSj, ip, 1:ntile/2)
                % will contain -ntiles/2:ntile:0.
                tempaway = [];
                tempto = [];
                if replay_state
                    tempaway_replay = [];
                    tempto_replay = [];
                end
                for t=1:length(offtimes)-1
                    t_interval = ontimes(t+1)-offtimes(t); 
                    if ~any(ismember(T_breaks,[offtimes(t):ontimes(t+1)]))
                        if t_interval>=p_low && t_interval<=p_high
                            if ntiles==2 % keep this because it's the historic way.
                                tempaway = [tempaway;vpath{iSj}(offtimes(t):offtimes(t)+floor(t_interval/2),:)];
                                tempto = [tempto;vpath{iSj}(ontimes(t+1)-floor(t_interval/2):ontimes(t+1),:)];
                                if replay_state
                                    tempaway_replay = [tempaway_replay;replay(offtimes(t):offtimes(t)+floor(t_interval/2),:)];
                                    tempto_replay = [tempto_replay; replay(ontimes(t+1)-floor(t_interval/2):ontimes(t+1),:)];
                                end
                            else
                                t_interval = t_interval+1;
                                tempaway = [tempaway;vpath{iSj}(offtimes(t)+(ntile-1)*floor(t_interval/ntiles):offtimes(t)+ntile*floor(t_interval/ntiles)-1,:)];
                                tempto = [tempto;vpath{iSj}(ontimes(t+1)-((ntiles/2)-(ntile-1))*floor(t_interval/ntiles)+1:ontimes(t+1)-((ntiles/2)-(ntile))*floor(t_interval/ntiles),:)];
                                if replay_state
                                    tempaway_replay = [tempaway_replay;replay(offtimes(t)+(ntile-1)*floor(t_interval/ntiles):offtimes(t)+ntile*floor(t_interval/ntiles)-1,:)];
                                    tempto_replay = [tempto_replay;replay(ontimes(t+1)-((ntiles/2)-(ntile-1))*floor(t_interval/ntiles)+1:ontimes(t+1)-((ntiles/2)-(ntile))*floor(t_interval/ntiles),:)];
                                end
                                
                            end
                            intervalDesign(offtimes(t):offtimes(t)+floor(t_interval/ntiles), 1) = 1;
                            intervalDesign(ontimes(t+1)-floor(t_interval/ntiles):ontimes(t+1), 1) = 2;
                            %intervaltracker{ip,ik}(length(intervaltracker{ip,ik})+1) = t_interval;
                        end
                    end
                end
                for ik2=1:K+replay_state
                    if ~isempty(tempaway)
                        if ik2>K % replay state
                            FO(ik,ik2,1,iSj,ip,ntile) = mean(tempaway_replay==1);
                            FO(ik,ik2,2,iSj,ip,ntile) = mean(tempto_replay==1);
                        else
                            FO(ik,ik2,1,iSj,ip,ntile) = mean(tempaway==ik2);
                            FO(ik,ik2,2,iSj,ip,ntile) = mean(tempto==ik2);
                        end
                    end
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

if ntiles>2
    % reshape such that the third dimension is ntiles in the order of one
    % state activation to the next
    [s1,s2,s3,s4,s5,s6] = size(FO);
    FO = permute(reshape(permute(FO, [1,2,4,5,6,3]), [s1,s2,s4,s5,s6*s3]), [1,2,5,3,4]);
end

if dostat
    % Do Permutation test on FO (instead of t-test)
    pvals = zeros(K,K,length(intervalbins)-1);
    stat=[];
    if dotaskglm
        pvals_residual = zeros(K,K,length(intervalbins)-1);
    else
        pvals_residual=[]; stat_res=[];
    end
    if nSj>1
        for ip=1:length(intervalbins)-1
            % do permutation test of FO
            
            [pvals(:,:,ip), stat{ip}] = FO_permutation_test(FO(:,:,:,:,ip), K, nSj, numperm);
            
            if dotaskglm
                [pvals_residual(:,:,ip), stat_res{ip}] = FO_permutation_test(FO_residual(:,:,:,:,ip), K, nSj, numperm);
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
else
    stat=[];
    pvals=[];
end