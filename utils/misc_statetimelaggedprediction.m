%% try multinomial fit:
lags = [10:10:250, 500, 750, 1000, 1250];
acc_mnr = zeros(numel(lags),config.nSj);

for subnum=1:config.nSj
    subnum
    vpath_onehot = zeros(size(vpath{subnum},1),K);
    for istate = 1:K
        vpath_onehot(:,istate) = vpath{subnum} == istate;
    end
    cnt=0;
    for t=lags
        cnt=cnt+1;
        preds = zeros(length(vpath_onehot)-t,K);
        for istate = 1:K
            predictors = vpath_onehot(1:end-t,:);
            var = vpath_onehot(1+t:end,istate);
            % do 2 fold cross validation:
            T = floor(length(predictors)/2);
            for iCV=1:2
                if iCV==1
                    mdl = fitglm(predictors(1:T,:),var(1:T),"Distribution","binomial",'link','logit');
                    preds(T+1:end,istate) = mdl.predict(predictors(T+1:end,:));
                    %preds = preds > 0.5;
                    %acc(t,istate,subnum,iCV) = mean(~xor(preds,var(T+1:end)));
                else
                    mdl = fitglm(predictors(T+1:end,:),var(T+1:end),"Distribution","binomial",'link','logit');
                    preds(1:T,istate) = mdl.predict(predictors(1:T,:));
                    %preds = preds > 0.5;
                    %acc(t,istate,subnum,iCV) = mean(~xor(preds,var(1:T)));
                end
            end
        end
        % find max likelihood state:
        vpath_preds = zeros(length(predictors),1);
        for i=1:length(predictors)
            [~,vpath_preds(i)] = max(preds(i,:));
        end
        acc_mnr(cnt,subnum) = mean(vpath_preds==vpath{subnum}(1+t:end));
        
    end
end
figure();
shadedErrorBar(lags, mean(acc_mnr,2), std(acc_mnr,[],2))
xticks(1:length(lags)), xticklabels(lags)



prob=zeros(12,12,250, 55);
for subnum=1:config.nSj
    subnum
    vpath_onehot = zeros(size(vpath{subnum},1),K);
    for istate = 1:K
        vpath_onehot(:,istate) = vpath{subnum} == istate;
    end
    V = vpath_onehot;
    V(V==0)=nan;
    for t=[500,1250]%1:250
        for k1=1:12
            for k2=1:12
                
                prob(k1,k2,t, subnum) = mean(V(1:end-t,k1)==V(1+t:end, k2))./mean(vpath_onehot(:,k2));
            end
        end
    end
end
hmm_1stlevel.control.timelagged.transprob = prob;
%{
for tlag = 10:10:250
    for iSj=1:config.nSj
        P_sim = prob(:,:,tlag,iSj)';
        P_sim = cumsum(P_sim,2);
        P_sim(:, end) = 1;
        t0=0;
        T=hmmT{iSj};
        vpath_new = zeros(size(vpath{iSj}));
        for ik=1:K
            FO(ik) = mean(vpath{iSj}==ik);
        end
        
        for iseg=1:length(T)
            z = rand(T(iseg),1);
            vpath_new(t0+1) = find(z(1)<cumsum(FO,2),1);
            for iT = 2:T(iseg)
                try % wierd syntax here just to catch very rare precision errors
                    vpath_new(t0+iT) = find(z(iT)<P_sim(vpath_new(iT-1),:),1);
                catch
                    vpath_new(t0+iT) = find(z(iT)<P_sim(vpath_new(iT-1),:),1);
                end
            end
            t0 = t0+T(iseg);
        end
        vpath_sim{iSj} = vpath_new;
    end
    [simulation{tlag/10}.FO_intervals,~,~, ~] = computeLongTermAsymmetry(vpath_sim,hmmT,K,[],[],[],false);
    a=[];
    for i=1:K
        for j=1:K
            [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(simulation{tlag/10}.FO_intervals(i,j,1,:)), squeeze(simulation{tlag/10}.FO_intervals(i,j,2,:)));
        end
    end
    simulation{tlag/10}.assym_ttest = a;
    
    simulation{tlag/10}.bestsequencemetrics_sim = optimiseSequentialPattern(simulation{tlag/10}.FO_intervals);
    simulation{tlag/10}.bestseq = simulation{tlag/10}.bestsequencemetrics_sim{1};
    angleplot_sim = circle_angles(simulation{tlag/10}.bestseq);
    
    simulation{tlag/10}.cycle_metrics = compute_tinda_metrics(config, simulation{tlag/10}.bestseq, simulation{tlag/10}.FO_intervals, simulation{tlag/10}.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);
    simulation{tlag/10}.lag = tlag/1000;

end
hmm_1stlevel.control.timelagged.time = (1:250)./1000;
hmm_1stlevel.control.timelagged.simulation = simulation;
%}

%%


cnt=0;
lags = [2,5,12,25,50,125,250,500,1250];
init_length = 2*max(lags);
for tlag = lags
    cnt=cnt+1;
    for iSj=1:config.nSj
        P_sim = prob(:,:,tlag,iSj)';
        P_sim = cumsum(P_sim,2);
        P_sim(:, end) = 1;
        for ik=1:K
            FO(ik) = mean(vpath{iSj}==ik);
        end
        T = length(vpath{iSj})+init_length;
        vpath_new = zeros(T,1);
        % run an initial 5000 samples based on a random starting point in
        % of vpath
        ix_init = randperm(length(vpath{iSj})-tlag-1,1);
        vpath_new(1:tlag) = vpath{iSj}(ix_init:ix_init+tlag-1);
        
        z = rand(T,1);
        for iT =tlag+1:T
            try % wierd syntax here just to catch very rare precision errors
                vpath_new(iT) = find(z(iT)<P_sim(vpath_new(iT-tlag),:),1);
            catch
                vpath_new(iT) = find(z(iT)<P_sim(vpath_new(iT-tlag),:),1);
            end
        end
        
        vpath_sim{iSj} = vpath_new(init_length+1:end);
    end
    [simulation{cnt}.FO_intervals,~,~, ~] = computeLongTermAsymmetry(vpath_sim,hmmT,K,[],[],[],false);
    a=[];
    for i=1:K
        for j=1:K
            [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(simulation{cnt}.FO_intervals(i,j,1,:)), squeeze(simulation{cnt}.FO_intervals(i,j,2,:)));
        end
    end
    simulation{cnt}.assym_ttest = a;
    
    simulation{cnt}.bestsequencemetrics_sim = optimiseSequentialPattern(simulation{cnt}.FO_intervals);
    simulation{cnt}.bestseq = simulation{cnt}.bestsequencemetrics_sim{1};
    
    simulation{cnt}.cycle_metrics = compute_tinda_metrics(config, simulation{cnt}.bestseq, simulation{cnt}.FO_intervals, simulation{cnt}.assym_ttest.pvals<hmm_1stlevel.assym_ttest.alpha_thresh, color_scheme, false);
    simulation{cnt}.lag = tlag/1000;

end
hmm_1stlevel.control.timelagged.time = (1:250)./1000;
hmm_1stlevel.control.timelagged.simulation = simulation;




%%
fig = setup_figure([],1,1)
clear q;
for k=1:length(hmm_1stlevel.control.timelagged.simulation)
    q(:,k) = (hmm_1stlevel.control.timelagged.simulation{k}.cycle_metrics.rotational_momentum)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
end
shadedErrorBar(1:length(lags), mean(q), std(q)./sqrt(config.nSj))
hline(mean(hmm_1stlevel.cycle_metrics.rotational_momentum./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum), '--k')
yl=ylim;
ylim([1.1*yl(1) 0])
xlabel('Time lag (s)')
ylabel('M')
text(0.4, .1, 'Observed', 'Units', 'normalized')
title({'Mean (+SEM) rotational momentum', 'from time-lagged models'})
box off
xticks(1:length(lags)), xticklabels(lags)
xlim([1 length(lags)])

save_figure([config.figdir, 'figure_supp_tinda_metrics/','2supp_rotational_momentum_timelagged']);

