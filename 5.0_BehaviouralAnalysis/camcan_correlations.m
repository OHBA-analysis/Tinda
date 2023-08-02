whichstudy = 6; % 1 denotes the hmm model run in Higgins2020_neuron
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = colorscheme(whichstudy);

info = camcan_getparticipantinfo(config);
subj_age = info(:,1);
subj_gender = info(:,3); % 1 denotes male, 2 female
subj_RTs = info(:,5);
info = [subj_age,subj_gender,subj_RTs];

% load the cyclical summary metrics
load(config.metricfile);

cycletime_mu = hmm_2ndlevel.cyctime_mu;
cycletime_med = hmm_2ndlevel.cyctime_med;
cycletime_std = hmm_2ndlevel.cyctime_std;

cyclerate_mu = 1./cycletime_mu; % this is more gaussian distributed

X = [cyclerate_mu, hmm_1stlevel.cycle_metrics.rotational_momentum];
pred_label = {'constant', 'cycle rate', 'rotational momentum'};
for k=1:12
  pred_label{3+k} = sprintf('FO state %d', k);
end

resp_label = {'age', 'sex'};
% find outliers
outliers = any(transpose(abs(X-mean(X,1)) > 3*std(X,[],1)));
for k = 1:2
  if k==1
    y = subj_age;
    distribution = 'normal';
  elseif k==2
    y = subj_gender-1;
    distribution = 'binomial';
  end
  [b{k},dev{k},stats{k}] = glmfit(X(~outliers,:), y(~outliers), distribution);
end

save([config.resultsdir, 'correlations_demographics'], 'b', 'dev', 'stats', 'pred_label', 'resp_label')
%% plot age correlation
clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

tmp1 = [{subj_age(~outliers)}, {subj_gender(~outliers)}];
tmp2 = [{cyclerate_mu(~outliers)}, {hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum}];
ylb = {'Cycle rate (Hz)', 'Cycle strength'};
ylb2 = {'cycrate', 'rotmom'};
xlb = {'Age', 'Sex'};
ttl = {'Cycle rate vs. age', 'Cycle strength vs. age', 'Cycle rate vs. sex', 'Cycle strength vs. sex'};
pos = {[105 2.45], [105 -0.02], [1.5 2.4], [1.5 -0.018]};
xl = {[0 119], [0.5, 2.5]};
cnt=1;

for k1=1
    for k2=1:2
        fig = setup_figure([], 1,1); hold on
        scatter(tmp1{k1}, (tmp2{k2}), 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{1}, 'MarkerEdgeAlpha', 0.2)
        h1 = lsline()
        h1.Color = 'k';
        h1.LineWidth = 2;
        xlabel(xlb{k1}), ylabel(ylb{k2})
        if k1==2
            xticks([1,2])
            xticklabels({'Male', 'Female'})
        end
        rho = corr(tmp1{k1}, tmp2{k2});
        text(pos{cnt}(1), pos{cnt}(2),{sprintf('R=%.02f', rho), sprintf('p=%0.04f', stats{k1}.p(1+k2))}, 'HorizontalAlignment', 'center')
        
        xlim(xl{k1})
        box off
        title(ttl{cnt})
        cnt=cnt+1;
        save_figure([config.figdir sprintf('figure4_correlations/4_correlation_%s_%s', xlb{k1}, ylb2{k2})])
        
    end
end

% find out how much % the rate increases every 10 years
% 10*b{1}(2)/log(10)


%%
group=subj_gender(~outliers);
for k1=2
    for k2=1:2
        fig = setup_figure([], 1,1); hold on
        
        boxplot_with_scatter((tmp2{k2}), [], [], group)
        xlabel(xlb{k1}), ylabel(ylb{k2})
        if k1==2
            xticks([1,2])
            xticklabels({'Male', 'Female'})
        end
        mu = [mean(tmp2{k2}(group==1)), mean(tmp2{k2}(group==2))];
        sig = [std(tmp2{k2}(group==1)), std(tmp2{k2}(group==2))];
        sig_pool = sqrt(((sum(group==1)-1)*sig(1).^2 + (sum(group==2)-1)*sig(2).^2)/(length(group)-2));
        cohend = diff(mu)./sig_pool;
        h=sigstar([1,2], 4*stats{k1}.p(1+k2))
%         text(pos{cnt}(1), pos{cnt}(2),{sprintf("Cohen's d=%.02f", cohend), sprintf('p=%0.04f', 4*stats{k1}.p(1+k2))}, 'HorizontalAlignment', 'center')
        xlim(xl{k1})
        box off
        title(ttl{cnt})
        cnt=cnt+1;

        save_figure([config.figdir sprintf('figure4_correlations/4_correlation_%s_%s', xlb{k1}, ylb2{k2})])
        
    end
end



%% plot in a single figure together with the heritability
clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};
load([config.basedir, 'Study3/HMMsummarymetrics.mat'], 'heritability')

fig = setup_figure([], 2,1),

ax(1) = axes('Position', [0.125 0.575 0.25 0.35]);
k1=1; %age, sex
k2=1; % rate, strength
cnt=1;
scatter(tmp1{k1}, (tmp2{k2}), 'MarkerFaceColor', clr{3}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{3}, 'MarkerEdgeAlpha', 0.2);
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;
xlabel('\bfAge'), 
rho = corr(tmp1{k1}, tmp2{k2});
title({sprintf('R = %.02f', rho), sprintf('p = %0.04f', stats{k1}.p(1+k2))})
xlim([10 95])
xticks(20:20:80)
ylim([1.5 5])
box off
ylabel({'\bf \fontsize{15} Cycle rate (Hz) \rm', ''})


cnt=2; k1=2; k2=1;
ax(cnt) = axes('Position', [0.4 0.575 0.25 0.3675]);
boxplot_with_scatter((tmp2{k2}), [], .3, group)
xlabel('\bfSex'), 
ax(cnt).YAxis.Visible = 'off'; 
xticks([1,2])
xticklabels({'Male', 'Female'})
mu = [mean(tmp2{k2}(group==1)), mean(tmp2{k2}(group==2))];
sig = [std(tmp2{k2}(group==1)), std(tmp2{k2}(group==2))];
sig_pool = sqrt(((sum(group==1)-1)*sig(1).^2 + (sum(group==2)-1)*sig(2).^2)/(length(group)-2));
cohend = diff(mu)./sig_pool;
h=sigstar([1,2], stats{k1}.p(1+k2))
        title({sprintf("Cohen's d = %.02f", cohend), sprintf('p = %0.04f', stats{k1}.p(1+k2))})
xlim(xl{k1})
ylim([1.5 5])
box off

ax(5) = axes('Position', [0.725 0.575 0.25, 0.37]);
im=1;
absdiff = [abs((mean(heritability.data{im}(heritability.pairs{1,1},:),2) - mean(heritability.data{im}(heritability.pairs{1,2},:),2)));...
    abs((mean(heritability.data{im}(heritability.pairs{2,1},:),2) - mean(heritability.data{im}(heritability.pairs{2,2},:),2)));...
    abs((mean(heritability.data{im}(heritability.pairs{3,1},:),2) - mean(heritability.data{im}(heritability.pairs{3,2},:),2)))];
G = []; for k=1:3; G=[G;k*ones(size(heritability.pairs{k,1}))]; end

c={[0.3020    0.6863    0.2902], [0.5961    0.3059    0.6392], [0.5020    0.6941    0.8275]};
boxplot_with_scatter(absdiff, c, 0.6, G)
box off
xticklabels({'MZ', 'DZ', 'unrelated'})
ylabel('\Delta pairs')
title({sprintf('h^2 = %.02f', heritability.cycle_rate.Ests(1,1)), sprintf('p = %0.04f', heritability.cycle_rate.Ps(1,1))})
xlabel('\bfRelatedness')



cnt=3; k1=1;k2=2;
ax(cnt) = axes('Position', [0.125 0.075 0.25 0.35]);
scatter(tmp1{k1}, (tmp2{k2}), 'MarkerFaceColor', clr{3}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{3}, 'MarkerEdgeAlpha', 0.2);
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;
xlabel('\bfAge'), 
ylabel(ylb{k2})
rho = corr(tmp1{k1}, tmp2{k2});
title({sprintf('R = %.02f', rho), sprintf('p = %0.04f', stats{k1}.p(1+k2))})
xlim([10 95])
xticks(20:20:80)
ylim([-.06 .2])
box off
% text(0.8,1.2, '\bfCycle strength\rm', 'Units', 'Normalized')
ylabel({'\bf \fontsize{15} Cycle strength (a.u.) \rm', ''})


cnt=2; k1=2; k2=2;
ax(cnt) = axes('Position', [0.4 0.075 0.25 0.3675]);
boxplot_with_scatter((tmp2{k2}), [], .3, group)
xlabel('\bfSex'), 
ax(cnt).YAxis.Visible = 'off'; 
xticks([1,2])
xticklabels({'Male', 'Female'})
mu = [mean(tmp2{k2}(group==1)), mean(tmp2{k2}(group==2))];
sig = [std(tmp2{k2}(group==1)), std(tmp2{k2}(group==2))];
sig_pool = sqrt(((sum(group==1)-1)*sig(1).^2 + (sum(group==2)-1)*sig(2).^2)/(length(group)-2));
cohend = diff(mu)./sig_pool;
h=sigstar([1,2], 4*stats{k1}.p(1+k2))
        title({sprintf("Cohen's d = %.02f", cohend), sprintf('p = %0.04f (n.s.)', stats{k1}.p(1+k2))})
xlim(xl{k1})
ylim([-.06 .2])
box off



ax(6) = axes('Position', [0.725 0.075 0.25, 0.37]);
im=2;
absdiff = [abs((mean(heritability.data{im}(heritability.pairs{1,1},:),2) - mean(heritability.data{im}(heritability.pairs{1,2},:),2)));...
    abs((mean(heritability.data{im}(heritability.pairs{2,1},:),2) - mean(heritability.data{im}(heritability.pairs{2,2},:),2)));...
    abs((mean(heritability.data{im}(heritability.pairs{3,1},:),2) - mean(heritability.data{im}(heritability.pairs{3,2},:),2)))];
G = []; for k=1:3; G=[G;k*ones(size(heritability.pairs{k,1}))]; end

c={[0.3020    0.6863    0.2902], [0.5961    0.3059    0.6392], [0.5020    0.6941    0.8275]};
boxplot_with_scatter(absdiff, c, 0.6, G)
box off
xticklabels({'MZ', 'DZ', 'unrelated'})
ylabel('\Delta pairs')
title({sprintf('h^2 = %.02f', heritability.rotational_momentum.Ests(1,1)), sprintf('p = %0.04f (n.s.)', heritability.rotational_momentum.Ps(1,1))})
xlabel('\bfRelatedness')


save_figure([config.figdir 'figure4_correlations/4_correlations_all'],[],false)

%% Correlation between cycle rate and strength
age = info(:,1);

correlation_cycle_rate_strength=[];
correlation_cycle_rate_strength.cycle_strength_corrected = regress_out(hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers), age(~outliers));
correlation_cycle_rate_strength.cycle_rate_corrected = regress_out(cyclerate_mu(~outliers), age(~outliers));

[correlation_cycle_rate_strength.R,correlation_cycle_rate_strength.pval]=corr(correlation_cycle_rate_strength.cycle_rate_corrected, correlation_cycle_rate_strength.cycle_strength_corrected);

save([config.resultsdir, 'correlations_demographics'],'correlation_cycle_rate_strength', '-append')


fig=setup_figure([],1,1);
scatter(cyclerate_mu(~outliers), hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum, 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{1}, 'MarkerEdgeAlpha', 0.2)
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;
xlabel('\bfCycle rate (Hz)'), 
ylabel('\bfCycle strength (a.u.)')
rho = corr(cyclerate_mu(~outliers), hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum);
title({sprintf('R = %.02f', correlation_cycle_rate_strength.R), sprintf('p = %0.01f * 10^{-19}', 10^19*correlation_cycle_rate_strength.pval)})
box off
save_figure([config.figdir, 'figure4_correlations/4_correlation_cycle_rate_strength'])


%% Correlation between age and mean absolute asymmetry

correlation_tida_age=[];
correlation_tida_age.tida = squeeze(nanmean(nanmean(abs(hmm_1stlevel.cycle_metrics.FO_assym(:,:,~outliers)),2),1));
correlation_tida_age.age = age(~outliers);

[correlation_tida_age.R,correlation_tida_age.pval]=corr(correlation_tida_age.tida, correlation_tida_age.age);

% see if this holds when correcting for cycle strength
correlation_tida_age.rotational_momentum = hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers);
[correlation_tida_age.R_corrected_for_M,correlation_tida_age.pval_corrected_for_M]=partialcorr(correlation_tida_age.tida, correlation_tida_age.age, correlation_tida_age.rotational_momentum);

% or the other way around
[correlation_tida_age.R_M_corrected_for_tida,correlation_tida_age.pval_M_corrected_for_tida]=partialcorr(correlation_tida_age.age, correlation_tida_age.rotational_momentum, correlation_tida_age.tida);

save([config.resultsdir, 'correlations_demographics'],'correlation_tida_age', '-append')


fig=setup_figure([],1,1);
scatter(correlation_tida_age.age, correlation_tida_age.tida, 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{1}, 'MarkerEdgeAlpha', 0.2)
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;
xlabel('\bfAge'), 
ylabel('\bfMean absolute FO asymmetry')
rho = corr(correlation_tida_age.age, correlation_tida_age.tida);
xlim([10 95])
xticks(20:20:80)
title({sprintf('R = %.02f', correlation_tida_age.R), sprintf('p = %0.04f', correlation_tida_age.pval)})
box off
save_figure([config.figdir, 'figure4_correlations/4_correlation_tida_age'])

%% Correlation between age and switching rate

correlation_switching_rate_age=[];
for iSj=1:config.nSj
    correlation_switching_rate_age.switching_rate(iSj,1) = mean(diff(vpath{iSj})~=0);
end
correlation_switching_rate_age.switching_rate = correlation_switching_rate_age.switching_rate(~outliers)
correlation_switching_rate_age.age = age(~outliers);

[correlation_switching_rate_age.R,correlation_switching_rate_age.pval]=corr(correlation_switching_rate_age.switching_rate, correlation_switching_rate_age.age);

save([config.resultsdir, 'correlations_demographics'],'correlation_switching_rate_age', '-append')


fig=setup_figure([],1,1);
scatter(correlation_switching_rate_age.age, correlation_switching_rate_age.switching_rate, 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{1}, 'MarkerEdgeAlpha', 0.2)
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;
xlabel('\bfAge'), 
ylabel('\bfMean state-switching rate')
rho = corr(correlation_switching_rate_age.age, correlation_switching_rate_age.switching_rate);
xlim([10 95])
xticks(20:20:80)
title({sprintf('R = %.02f', correlation_switching_rate_age.R), sprintf('p = %0.04f', correlation_switching_rate_age.pval)})
box off
save_figure([config.figdir, 'figure4_correlations/4_correlation_switching_rate_age'])


%% difference in bestseq - Cycle metric correlations with age are not caused by differences in bestseq
[~, ix] = sort(age, 'ascend');
ix_young=ix(1:100);
ix_old = ix(end-99:end);


bestsequencemetrics_young = optimiseSequentialPattern(hmm_1stlevel.FO_intervals(:,:,:,ix_young));
bestsequencemetrics_old = optimiseSequentialPattern(hmm_1stlevel.FO_intervals(:,:,:,ix_old));



angleplot = circle_angles(bestsequencemetrics_young{1});
M_young = compute_rotational_momentum(angleplot, hmm_1stlevel.cycle_metrics.FO_assym);
[R_young, pval_young] = corr(age(~outliers), M_young(~outliers));

angleplot = circle_angles(bestsequencemetrics_old{1});
M_old = compute_rotational_momentum(angleplot, hmm_1stlevel.cycle_metrics.FO_assym);
[R_old, pval_old] = corr(age(~outliers), M_old(~outliers));

bestseq_old_young=[];
bestseq_old_young.ix_old = ix_old;
bestseq_old_young.ix_young = ix_young;
bestseq_old_young.bestseq_young.bestsequencemetrics = bestsequencemetrics_young;
bestseq_old_young.bestseq_young.rotational_momentum = M_young;
bestseq_old_young.bestseq_young.R = R_young;
bestseq_old_young.bestseq_young.pval = pval_young;

bestseq_old_young.bestseq_old.bestsequencemetrics = bestsequencemetrics_old;
bestseq_old_young.bestseq_old.rotational_momentum = M_old;
bestseq_old_young.bestseq_old.R = R_old;
bestseq_old_young.bestseq_old.pval = pval_old;

save([config.resultsdir, 'correlations_demographics'],'bestseq_old_young', '-append')


%% Correlation between connectivity in states and cycle metrics/age
correlation_cyle_metrics_coherence=[];
[correlation_cyle_metrics_coherence.R_strength_coh, correlation_cyle_metrics_coherence.p_strength_coh] = corr(nanmean(coh_wb(~outliers,:,:),3), hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers));

[correlation_cyle_metrics_coherence.R_rate_coh, correlation_cyle_metrics_coherence.p_rate_coh] = corr(nanmean(coh_wb(~outliers,:,:),3), cyclerate_mu(~outliers));
[correlation_cyle_metrics_coherence.R_coh_age,correlation_cyle_metrics_coherence.p_coh_age] = corr(nanmean(coh_wb(~outliers,:,:),3), age(~outliers));

save([config.resultsdir, 'correlations_demographics'],'correlation_cyle_metrics_coherence', '-append')

theta = (0:1/12:.99)*2*pi; % define angles
theta = circshift(theta,-3); % set first angle to 12 o'clock
for k=1:K % reshuffle angles in bestseq order
tmp(k) = theta(bestseq==k);
end
theta = tmp;

xpos = cos(theta);
xpos = (xpos+1)/2.6+0.085;
ypos = sin(theta);
ypos = (ypos+1)/2.5+0.05;

% fig=setup_figure([],1.5,1);
fig=setup_figure([], 1.5, 2)
for k=1:12
% ax(k,1) = axes('Position', [xpos(k)-0.05, ypos(k), 0.075, 0.075]);
ax(k,1) = axes('Position', [0.06,.92-0.98*(k/13), 0.25, 1/18]);
scatter(nanmean(coh_wb(~outliers,k,:),3), hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum, 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{1}, 'MarkerEdgeAlpha', 0.2)
xticks([])
yticks([-.1 -.05 0 .05])
title(sprintf('State %d', k))

if k==1
    if correlation_cyle_metrics_coherence.p_strength_coh(k)<0.05/12
        title({'Cycle strength', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k))})
    else
        title({'Cycle strength', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k))})
    end
else
    if correlation_cyle_metrics_coherence.p_strength_coh(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k)))
    end
end
if k==12
   xlabel('Coherence') 
end
hline(0)
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;

% ax(k,2) = axes('Position', [xpos(k)+0.05, ypos(k), 0.075, 0.075]);
ax(k,2) = axes('Position', [0.38,.92-0.98*(k/13), 0.25, 1/18]);
scatter(nanmean(coh_wb(~outliers,k,:),3), cyclerate_mu(~outliers), 'MarkerFaceColor', clr{2}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{2}, 'MarkerEdgeAlpha', 0.2)
xticks([])
yticks([1,2,3])

if k==1
    if correlation_cyle_metrics_coherence.p_rate_coh(k)<0.05/12
        title({'Cycle rate', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k))})
    else
        title({'Cycle rate', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k))})
    end
else
    if correlation_cyle_metrics_coherence.p_rate_coh(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k)))
    end
end
if k==12
   xlabel('Coherence') 
end

h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;

ax(k,3) = axes('Position', [0.71,.92-0.98*(k/13), 0.25, 1/18]);
scatter(nanmean(coh_wb(~outliers,k,:),3), age(~outliers), 'MarkerFaceColor', clr{3}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{3}, 'MarkerEdgeAlpha', 0.2)
xticks([])
% yticks()

if k==1
    if correlation_cyle_metrics_coherence.p_coh_age(k)<0.05/12
        title({'Age', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_coh_age(k))})
    else
        title({'Age', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_coh_age(k))})
    end
else
    if correlation_cyle_metrics_coherence.p_coh_age(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_coh_age(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_coh_age(k)))
    end
end
if k==12
   xlabel('Coherence') 
end

h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;




end
sgtitle({'Correlation between', 'coherence and cycle metrics', ''})
save_figure([config.figdir, 'figure4_correlations/4_correlation_coherence'], false)

%%
correlation_cyle_metrics_pow=[];
[correlation_cyle_metrics_pow.R_strength_pow, correlation_cyle_metrics_pow.p_strength_pow] = corr(nanmean(psd_wb(~outliers,:,:),3), hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers));

[correlation_cyle_metrics_pow.R_rate_pow, correlation_cyle_metrics_pow.p_rate_pow] = corr(nanmean(psd_wb(~outliers,:,:),3), cyclerate_mu(~outliers));
[correlation_cyle_metrics_pow.R_pow_age,correlation_cyle_metrics_pow.p_pow_age] = corr(nanmean(psd_wb(~outliers,:,:),3), age(~outliers));

save([config.resultsdir, 'correlations_demographics'],'correlation_cyle_metrics_pow', '-append')

theta = (0:1/12:.99)*2*pi; % define angles
theta = circshift(theta,-3); % set first angle to 12 o'clock
for k=1:K % reshuffle angles in bestseq order
tmp(k) = theta(bestseq==k);
end
theta = tmp;

xpos = cos(theta);
xpos = (xpos+1)/2.6+0.085;
ypos = sin(theta);
ypos = (ypos+1)/2.5+0.05;

% fig=setup_figure([],1.5,1);
fig=setup_figure([], 1.5, 2)
for k=1:12
% ax(k,1) = axes('Position', [xpos(k)-0.05, ypos(k), 0.075, 0.075]);
ax(k,1) = axes('Position', [0.06,.92-0.98*(k/13), 0.25, 1/18]);
scatter(nanmean(psd_wb(~outliers,k,:),3), hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum, 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{1}, 'MarkerEdgeAlpha', 0.2)
xticks([])
yticks([-.1 -.05 0 .05])
title(sprintf('State %d', k))

if k==1
    if correlation_cyle_metrics_pow.p_strength_pow(k)<0.05/12
        title({'Cycle strength', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_pow.R_strength_pow(k))})
    else
        title({'Cycle strength', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_pow.R_strength_pow(k))})
    end
else
    if correlation_cyle_metrics_pow.p_strength_pow(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_pow.R_strength_pow(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_pow.R_strength_pow(k)))
    end
end
if k==12
   xlabel('Power') 
end
hline(0)
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;

% ax(k,2) = axes('Position', [xpos(k)+0.05, ypos(k), 0.075, 0.075]);
ax(k,2) = axes('Position', [0.38,.92-0.98*(k/13), 0.25, 1/18]);
scatter(nanmean(psd_wb(~outliers,k,:),3), cyclerate_mu(~outliers), 'MarkerFaceColor', clr{2}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{2}, 'MarkerEdgeAlpha', 0.2)
xticks([])
yticks([1,2,3])

if k==1
    if correlation_cyle_metrics_pow.p_rate_pow(k)<0.05/12
        title({'Cycle rate', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_pow.R_rate_pow(k))})
    else
        title({'Cycle rate', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_pow.R_rate_pow(k))})
    end
else
    if correlation_cyle_metrics_pow.p_rate_pow(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_pow.R_rate_pow(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_pow.R_rate_pow(k)))
    end
end
if k==12
   xlabel('Power') 
end

h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;

ax(k,3) = axes('Position', [0.71,.92-0.98*(k/13), 0.25, 1/18]);
scatter(nanmean(psd_wb(~outliers,k,:),3), age(~outliers), 'MarkerFaceColor', clr{3}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{3}, 'MarkerEdgeAlpha', 0.2)
xticks([])
% yticks()

if k==1
    if correlation_cyle_metrics_pow.p_pow_age(k)<0.05/12
        title({'Age', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_pow.R_pow_age(k))})
    else
        title({'Age', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_pow.R_pow_age(k))})
    end
else
    if correlation_cyle_metrics_pow.p_pow_age(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_pow.R_pow_age(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_pow.R_pow_age(k)))
    end
end
if k==12
   xlabel('Power') 
end

h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;




end
sgtitle({'Correlation between', 'power and cycle metrics', ''})
save_figure([config.figdir, 'figure4_correlations/4_correlation_power'], false)

%%

theta = (0:1/12:.99)*2*pi; % define angles
theta = circshift(theta,-3); % set first angle to 12 o'clock
for k=1:K % reshuffle angles in bestseq order
tmp(k) = theta(bestseq==k);
end
theta = tmp;

xpos = cos(theta);
xpos = (xpos+1)/2.6+0.085;
ypos = sin(theta);
ypos = (ypos+1)/2.5+0.05;

% fig=setup_figure([],1.5,1);
fig=setup_figure([], 1,2.8)
for k=1:12
% ax(k,1) = axes('Position', [xpos(k)-0.05, ypos(k), 0.075, 0.075]);
ax(k,1) = axes('Position', [0.095,.92-0.98*(k/13), 0.4, 1/18]);
scatter(nanmean(coh_wb(~outliers,k,:),3), hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers)./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum, 'MarkerFaceColor', clr{1}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{1}, 'MarkerEdgeAlpha', 0.2)
xticks([])
yticks([-.1 -.05 0 .05])
title(sprintf('State %d', k))

if k==1
    if correlation_cyle_metrics_coherence.p_strength_coh(k)<0.05/12
        title({'Cycle strength', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k))})
    else
        title({'Cycle strength', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k))})
    end
else
    if correlation_cyle_metrics_coherence.p_strength_coh(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_strength_coh(k)))
    end
end
if k==12
   xlabel('Coherence') 
end
hline(0)
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;

% ax(k,2) = axes('Position', [xpos(k)+0.05, ypos(k), 0.075, 0.075]);
ax(k,1) = axes('Position', [0.57,.92-0.98*(k/13), 0.4, 1/18]);
scatter(nanmean(coh_wb(~outliers,k,:),3), cyclerate_mu(~outliers), 'MarkerFaceColor', clr{2}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{2}, 'MarkerEdgeAlpha', 0.2)
xticks([])
yticks([1,2,3])

if k==1
    if correlation_cyle_metrics_coherence.p_rate_coh(k)<0.05/12
        title({'Cycle rate', sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k))})
    else
        title({'Cycle rate', sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k))})
    end
else
    if correlation_cyle_metrics_coherence.p_rate_coh(k)<0.05/12
        title(sprintf('State %d* | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k)))
    else
        title(sprintf('State %d | R=%.02f', k, correlation_cyle_metrics_coherence.R_rate_coh(k)))
    end
end
if k==12
   xlabel('Coherence') 
end

h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;
end
sgtitle({'Correlation between', 'coherence and cycle metrics', ''})
save_figure([config.figdir, 'figure4_correlations/4_correlation_cycle_metrics_coherence'], false)

        %% CCA with intelligence
[cogdata,cogdatalabels] = camcan_getCognitiveData(config);

% remove subjID: 
cogdata = cogdata(:,1:15);
cogdatalabels = cogdatalabels(1:15);
% 1) FldIn   : fluid intelligence           f   ex
% 2) FacReg  : face recognition             f   emo
% 3) EmoRec  : emotional recognition        f   emo
% 4) MltTs   : multitask                    f   ex      (negatively correlated with others - indicates reaction times)
% 5) PicName : picture naming               f   lang
% 6) ProV    : proverb comprehension        c   lang
% 7) MRSp    : motor speed                  f   mot     (negatively correlated with others - indicates reaction times) 
% 8) MRCv    : motor speed varianve         f   mot
% 9) SntRec  : sentence comprehension       c/f lang 
% 10) VSTM   : visual short-term memory     f   mem
% 11) StrRec : story recall                 f   mem
% 12) StW    : spot the word                c   lang
% 13) VrbFl  : verbal fluency               f   lang
% col1: M/F; col2: age; cols3-15: cognitive measures; col16: CCID
cogdatalabels = {'Sex', 'Age', 'Fluid intelligence', 'Face recognition', ...
    'Emotional recognition', 'Multitasking', 'Picture naming', ...
    'Proverb comprehension', 'Motor speed', 'Motor speed var', ...
    'Sentence comprehension', 'Visual short-term memory', 'Story recall', ...
    'Spot the word', 'Verbal fluency'};

% if using nonparametric stats
%{
for k=1:100000
  a = hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers);
  a = a(randperm(length(a)));
  [~,~,~,~,~,stats_perm] = canoncorr(a,cogdata(~outliers,:));
  F_perm(k) = stats_perm.F;
end
%}
% regress out age and sex from all data.

regr = zscore(cogdata(:,1:2));
cycrate_mu_corrected =  demean(regress_out(zscore(1./hmm_2ndlevel.cyctime_mu), regr)); 
rotational_momentum_corrected =  demean(regress_out(zscore(hmm_1stlevel.cycle_metrics.rotational_momentum), regr));   
FO_corrected =  demean(regress_out(zscore(hmm_2ndlevel.FO), regr));   
cogdata_corrected =  demean(regress_out(zscore(cogdata(:,3:end)), regr));   
% find outlier in cogdata Motor Speed (7)
outliers = find(cogdata_corrected(:,7)>5);
keep = setdiff(1:length(cycrate_mu_corrected), outliers);

% remove outlier and renormalize
cogdata_corrected = demean(cogdata_corrected(keep, :));
cycrate_mu_corrected = demean(cycrate_mu_corrected(keep,:));
rotational_momentum_corrected = demean(rotational_momentum_corrected(keep,:));
FO_corrected = demean(FO_corrected(keep,:));


%% CCA
% without correction for age, sex:
% [A,B,r,U,V,stats] = canoncorr(hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers),cogdata(~outliers,3:end));
% with correction:
[A,B,r,U,V,stats] = canoncorr([rotational_momentum_corrected, cycrate_mu_corrected],cogdata_corrected);
%% permutations
nperm=10000
for k=1:nperm
  [~,~,r_perm(k,:),~,~,~] = canoncorr([rotational_momentum_corrected, cycrate_mu_corrected],cogdata_corrected(randperm(length(cogdata_corrected)),:));
end
cca.r = r;
cca.A=A;
cca.B=B;
cca.U=U;
cca.V=V;
cca.stats=stats;
cca.r_perm = r_perm;
cca.cogdatalabels = cogdatalabels;
cca.pval = sum(r_perm>r)./nperm;
cca.nperm=nperm;

save([config.resultsdir, 'correlations_demographics'], 'cca', '-append')

tbl=table;
[B_sorted, ix] = sort(B(:,2), 'descend');
tbl.weight = B_sorted;
tbl.size = abs(B_sorted);
tmp = cogdatalabels(3:end)';
tbl.labels=tmp(ix);

for k=1:length(B_sorted)
  if B_sorted(k)<0
    c{k} = clr{1};
  else
    c{k} = clr{2};
  end
end
%%
fig=setup_figure([],1,1);
ax(1) = axes('Position', [0.05 0.05 0.6 0.9])
w=wordcloud(tbl, 'labels', 'size', 'Color',cat(1,c{:}), 'Shape', 'rectangle')

ax(2) = axes('Position', [0.7 0.05 0.25 0.9])
tbl2=table;

tbl2.weight = A(:,2);
tbl2.size = [1;1] * max(A(:,2));

tbl2.labels={'Cycle strength', 'Cycle rate'}';
for k=1:length(A(:,2))
  if A(k,2)<0
    c2{k} = clr{1};
  else
    c2{k} = clr{2};
  end
end
w=wordcloud(tbl2, 'labels', 'size', 'Color',cat(1,c2{:}), 'Shape', 'rectangle')
save_figure([config.figdir 'figure4_correlations/4_CCA_cycle_metrics_wordcloud'])

%%





tmp = [B];
fig=setup_figure([],2,0.5);
ax(1) = axes('Position', [0.075 0.35 0.2 0.5])
bar(A(:,2))
xticklabels({'Cycle rate', 'Cycle strength'})
title('Cycle metric')
ylabel('Coefficient strength')
box off
% subplot(1,4,[2 4])
ax(2) = axes('Position', [0.35 0.35 0.6 0.5])
bar(B(:,2))
box off
xticklabels(cogdatalabels(3:end))
xtickangle(45)

title('Cognitive test')
sgtitle(sprintf('R=%0.2f, F(%d|%d)=%0.1f, p=%0.3f', r(2), stats.df1(2), stats.df2(2), stats.F(2), stats.p(2)))
save_figure([config.figdir 'figure4_correlations/4_CCA_cycle_metrics'])

%%
resp_label = {'cycle rate', 'rotational momentum'};
% find outliers
for k = 1:2
  if k==1
    y = cycrate_mu_corrected;
  elseif k==2
    y = rotational_momentum_corrected;
  end
  [b{k},dev{k},stats{k}] = glmfit(cogdata_corrected, y, 'normal');
end







%% old Cam's stuff
%{
reglabels = {'Age','Gender'};
for ireg=[1:2];
  labs = {'Mean Cycle ','Median Cycle ','Std Cycle '};
  figure('Position',[7 168 1147 630]);
  for itype = 1:2
    x = cycletime_mu(~outliers,1);
    label1 = 'time';
    if itype==2
      x = 1./x;
      label1 = 'rate';
    end
    subplot(2,3,1+(itype-1)*3);
    scatter(info(~outliers,ireg),x,'filled')
    [R,P] = corrcoef(info(~outliers,ireg),x);
    title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
    plot4paper(reglabels{ireg},[labs{1},label1]);
    H=lsline(gca)
    set(H,'LineWidth',2);
    axis square;
    
    x = cycletime_med(~outliers,1);
    if itype==2
      x = 1./x;
    end
    subplot(2,3,2+(itype-1)*3);
    scatter(info(~outliers,ireg),x,'filled');
    [R,P] = corrcoef(info(~outliers,ireg),x);
    title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
    plot4paper(reglabels{ireg},[labs{2},label1]);
    H=lsline(gca)
    set(H,'LineWidth',2);
    axis square;
    
    subplot(2,3,3+(itype-1)*3);
    scatter(info(~outliers,ireg),cycletime_std(~outliers,1),'filled')
    [R,P] = corrcoef(info(~outliers,ireg),cycletime_std(~outliers,1));
    title(['STD: rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
    plot4paper(reglabels{ireg},[labs{3},label1]);
    H=lsline(gca)
    set(H,'LineWidth',2);
    axis square;
    
    print([figdir '4A_',outlier_string,'metastateCycleCorr_',reglabels{ireg}],'-depsc')
  end
end


%}

