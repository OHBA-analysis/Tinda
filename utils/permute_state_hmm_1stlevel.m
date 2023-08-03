hmm_1stlevel.FO = hmm_1stlevel.FO(:,new_state_ordering);
hmm_1stlevel.LT_mu = hmm_1stlevel.LT_mu(:,new_state_ordering);
hmm_1stlevel.LT_med = hmm_1stlevel.LT_med(:,new_state_ordering);
hmm_1stlevel.LT_std = hmm_1stlevel.LT_std(:,new_state_ordering);
hmm_1stlevel.IT_mu = hmm_1stlevel.IT_mu(:,new_state_ordering);
hmm_1stlevel.IT_med = hmm_1stlevel.IT_med(:,new_state_ordering);
hmm_1stlevel.IT_std = hmm_1stlevel.IT_std(:,new_state_ordering);

hmm_1stlevel.FO_intervals = hmm_1stlevel.FO_intervals(new_state_ordering, new_state_ordering,:,:);
  
optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_coherence' ,'.mat'];
bestsequencemetrics = optimiseSequentialPattern(hmm_1stlevel.FO_intervals);
save(optimalseqfile,'bestsequencemetrics');

%% permtest
hmm_1stlevel.assym_permtest.pvals = hmm_1stlevel.assym_permtest.pvals(new_state_ordering, new_state_ordering);
hmm_1stlevel.assym_permtest.sigpoints = hmm_1stlevel.assym_permtest.sigpoints(new_state_ordering,new_state_ordering);
tmp = {'prob', 'cirange', 'mask', 'stat','ref'};
for k=1:length(tmp)
tmp2 = zeros(12);
tmp2(find(~eye(12))) = hmm_1stlevel.assym_permtest.stat.(tmp{k});
tmp2 = tmp2(new_state_ordering, new_state_ordering);
tmp2 = tmp2(find(~eye(12)));
hmm_1stlevel.assym_permtest.stat.(tmp{k}) = tmp2';
end

%% ttest
hmm_1stlevel.assym_ttest.h = hmm_1stlevel.assym_ttest.h(new_state_ordering,new_state_ordering);
hmm_1stlevel.assym_ttest.pvals = hmm_1stlevel.assym_ttest.pvals(new_state_ordering,new_state_ordering);
hmm_1stlevel.assym_ttest.ci = hmm_1stlevel.assym_ttest.ci(new_state_ordering,new_state_ordering,:);
hmm_1stlevel.assym_ttest.stat = hmm_1stlevel.assym_ttest.stat(new_state_ordering,new_state_ordering);
hmm_1stlevel.assym_ttest.sigpoints = hmm_1stlevel.assym_ttest.sigpoints(new_state_ordering,new_state_ordering);

%% cycle metrics
tmp = {'mean_direction', 'mean_assym','P'};
tmp2 = {'rotational_momentum_perstate', 'TIDA_perstate'};
for k=1:length(tmp);
    if k<3
        hmm_1stlevel.cycle_metrics.(tmp2{k}) = hmm_1stlevel.cycle_metrics.(tmp2{k})(:, new_state_ordering);
    end
    hmm_1stlevel.cycle_metrics.(tmp{k}) = hmm_1stlevel.cycle_metrics.(tmp{k})(new_state_ordering,new_state_ordering);
end

hmm_1stlevel.cycle_metrics.FO_assym = hmm_1stlevel.cycle_metrics.FO_assym(new_state_ordering,new_state_ordering,:);

%% simulations
for isim=1:100
hmm_1stlevel.simulation{isim}.FO_intervals=hmm_1stlevel.simulation{isim}.FO_intervals(new_state_ordering,new_state_ordering,:,:);

% permtest
hmm_1stlevel.simulation{isim}.assym_permtest.pvals = hmm_1stlevel.simulation{isim}.assym_permtest.pvals(new_state_ordering, new_state_ordering);
tmp = {'prob', 'cirange', 'mask', 'stat','ref'};
for k=1:length(tmp)
tmp2 = zeros(12);
tmp2(find(~eye(12))) = hmm_1stlevel.simulation{isim}.assym_permtest.stat.(tmp{k});
tmp2 = tmp2(new_state_ordering, new_state_ordering);
tmp2 = tmp2(find(~eye(12)));
hmm_1stlevel.simulation{isim}.assym_permtest.stat.(tmp{k}) = tmp2';
end

% ttest
hmm_1stlevel.simulation{isim}.assym_ttest.h = hmm_1stlevel.simulation{isim}.assym_ttest.h(new_state_ordering,new_state_ordering);
hmm_1stlevel.simulation{isim}.assym_ttest.pvals = hmm_1stlevel.simulation{isim}.assym_ttest.pvals(new_state_ordering,new_state_ordering);
hmm_1stlevel.simulation{isim}.assym_ttest.ci = hmm_1stlevel.simulation{isim}.assym_ttest.ci(new_state_ordering,new_state_ordering,:);
hmm_1stlevel.simulation{isim}.assym_ttest.stat = hmm_1stlevel.simulation{isim}.assym_ttest.stat(new_state_ordering,new_state_ordering);

% cycle metrics
tmp = {'mean_direction', 'mean_assym','P'};
tmp2 = {'rotational_momentum_perstate', 'TIDA_perstate'};
for k=1:length(tmp);
    if k<3
        hmm_1stlevel.simulation{isim}.cycle_metrics.(tmp2{k}) = hmm_1stlevel.simulation{isim}.cycle_metrics.(tmp2{k})(:, new_state_ordering);
    end
    hmm_1stlevel.simulation{isim}.cycle_metrics.(tmp{k}) = hmm_1stlevel.simulation{isim}.cycle_metrics.(tmp{k})(new_state_ordering,new_state_ordering);
end

hmm_1stlevel.simulation{isim}.cycle_metrics.FO_assym = hmm_1stlevel.simulation{isim}.cycle_metrics.FO_assym(new_state_ordering,new_state_ordering,:);
end
hmm_1stlevel.simulation{1}.bestsequencemetrics = optimiseSequentialPattern(hmm_1stlevel.simulation{1}.FO_intervals);


%% simulation average
hmm_1stlevel.simulation_average.FO_intervals=hmm_1stlevel.simulation_average.FO_intervals(new_state_ordering,new_state_ordering,:,:);

% permtest
hmm_1stlevel.simulation_average.assym_permtest.pvals = hmm_1stlevel.simulation_average.assym_permtest.pvals(new_state_ordering, new_state_ordering);
tmp = {'prob', 'cirange', 'mask', 'stat','ref'};
for k=1:length(tmp)
tmp2 = zeros(12);
tmp2(find(~eye(12))) = hmm_1stlevel.simulation_average.assym_permtest.stat.(tmp{k});
tmp2 = tmp2(new_state_ordering, new_state_ordering);
tmp2 = tmp2(find(~eye(12)));
hmm_1stlevel.simulation_average.assym_permtest.stat.(tmp{k}) = tmp2';
end

% ttest
hmm_1stlevel.simulation_average.assym_ttest.h = hmm_1stlevel.simulation_average.assym_ttest.h(new_state_ordering,new_state_ordering);
hmm_1stlevel.simulation_average.assym_ttest.pvals = hmm_1stlevel.simulation_average.assym_ttest.pvals(new_state_ordering,new_state_ordering);
hmm_1stlevel.simulation_average.assym_ttest.ci = hmm_1stlevel.simulation_average.assym_ttest.ci(new_state_ordering,new_state_ordering,:);
hmm_1stlevel.simulation_average.assym_ttest.stat = hmm_1stlevel.simulation_average.assym_ttest.stat(new_state_ordering,new_state_ordering);

% cycle metrics
tmp = {'mean_direction', 'mean_assym','P'};
tmp2 = {'rotational_momentum_perstate', 'TIDA_perstate'};
for k=1:length(tmp)
    if k<3
        hmm_1stlevel.simulation_average.cycle_metrics.(tmp2{k}) = hmm_1stlevel.simulation_average.cycle_metrics.(tmp2{k})(:, new_state_ordering);
    end
    hmm_1stlevel.simulation_average.cycle_metrics.(tmp{k}) = hmm_1stlevel.simulation_average.cycle_metrics.(tmp{k})(new_state_ordering,new_state_ordering);
end

hmm_1stlevel.simulation_average.cycle_metrics.FO_assym = hmm_1stlevel.simulation_average.cycle_metrics.FO_assym(new_state_ordering,new_state_ordering,:);
hmm_1stlevel.simulation_average.bestsequencemetrics = optimiseSequentialPattern(hmm_1stlevel.simulation_average.FO_intervals);


%% rest

hmm_1stlevel.FO_simulation_group = hmm_1stlevel.FO_simulation_group(new_state_ordering,new_state_ordering,:,:);





