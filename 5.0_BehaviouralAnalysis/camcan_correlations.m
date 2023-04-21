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

cycletime_mu = hmm_2ndlevel.cyctime_mu./1000;
cycletime_med = hmm_2ndlevel.cyctime_med./1000;
cycletime_std = hmm_2ndlevel.cyctime_std./1000;

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
  elseif k==2
    y = subj_gender;
  end
  [b{k},dev{k},stats{k}] = glmfit(X(~outliers,:), y(~outliers), 'normal');
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
        text(pos{cnt}(1), pos{cnt}(2),{sprintf('R=%.02f', rho), sprintf('p=%0.04f', 4*stats{k1}.p(1+k2))}, 'HorizontalAlignment', 'center')
        
        xlim(xl{k1})
        box off
        title(ttl{cnt})
        cnt=cnt+1;
%         save_figure([config.figdir sprintf('figure4_correlations/4_correlation_%s_%s', xlb{k1}, ylb2{k2})])
        
    end
end

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



%%
% and also save in a single figure together with the heritability
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
xlabel(xlb{k1}), 
rho = corr(tmp1{k1}, tmp2{k2});
title({sprintf('R = %.02f', rho), sprintf('p = %0.04f', 4*stats{k1}.p(1+k2))})
xlim([10 95])
xticks(20:20:80)
ylim([1 3.5])
box off
ylabel({'\bf \fontsize{15} Cycle rate \rm', ''})


cnt=2; k1=2; k2=1;
ax(cnt) = axes('Position', [0.4 0.575 0.25 0.3675]);
boxplot_with_scatter((tmp2{k2}), [], .3, group)
xlabel(xlb{k1}), 
ax(cnt).YAxis.Visible = 'off'; 
xticks([1,2])
xticklabels({'Male', 'Female'})
mu = [mean(tmp2{k2}(group==1)), mean(tmp2{k2}(group==2))];
sig = [std(tmp2{k2}(group==1)), std(tmp2{k2}(group==2))];
sig_pool = sqrt(((sum(group==1)-1)*sig(1).^2 + (sum(group==2)-1)*sig(2).^2)/(length(group)-2));
cohend = diff(mu)./sig_pool;
h=sigstar([1,2], 4*stats{k1}.p(1+k2))
        title({sprintf("Cohen's d = %.02f", cohend), sprintf('p = %0.04f', 4*stats{k1}.p(1+k2))})
xlim(xl{k1})
ylim([1 3.5])
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
title({sprintf('h^2 = %.02f', heritability.cycle_rate.Ests(1,1)), sprintf('p = %0.04f', 2*heritability.cycle_rate.Ps(1,1))})
xlabel('Relatedness')



cnt=3; k1=1;k2=2;
ax(cnt) = axes('Position', [0.125 0.075 0.25 0.35]);
scatter(tmp1{k1}, (tmp2{k2}), 'MarkerFaceColor', clr{3}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', clr{3}, 'MarkerEdgeAlpha', 0.2);
h1 = lsline()
h1.Color = 'k';
h1.LineWidth = 2;
xlabel(xlb{k1}), ylabel(ylb{k2})
rho = corr(tmp1{k1}, tmp2{k2});
title({sprintf('R = %.02f', rho), sprintf('p = %0.04f', 4*stats{k1}.p(1+k2))})
xlim([10 95])
xticks(20:20:80)
ylim([-.2 0.1])
box off
% text(0.8,1.2, '\bfCycle strength\rm', 'Units', 'Normalized')
ylabel({'\bf \fontsize{15} Cycle strength \rm', ''})


cnt=2; k1=2; k2=2;
ax(cnt) = axes('Position', [0.4 0.075 0.25 0.3675]);
boxplot_with_scatter((tmp2{k2}), [], .3, group)
xlabel(xlb{k1}), 
ax(cnt).YAxis.Visible = 'off'; 
xticks([1,2])
xticklabels({'Male', 'Female'})
mu = [mean(tmp2{k2}(group==1)), mean(tmp2{k2}(group==2))];
sig = [std(tmp2{k2}(group==1)), std(tmp2{k2}(group==2))];
sig_pool = sqrt(((sum(group==1)-1)*sig(1).^2 + (sum(group==2)-1)*sig(2).^2)/(length(group)-2));
cohend = diff(mu)./sig_pool;
h=sigstar([1,2], 4*stats{k1}.p(1+k2))
        title({sprintf("Cohen's d = %.02f", cohend), sprintf('p = %0.04f', 4*stats{k1}.p(1+k2))})
xlim(xl{k1})
ylim([-.2 .1])
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
title({sprintf('h^2 = %.02f', heritability.rotational_momentum.Ests(1,1)), sprintf('p = %0.04f', 2*heritability.rotational_momentum.Ps(1,1))})
xlabel('Relatedness')


        save_figure([config.figdir 'figure4_correlations/4_correlations_all'],[],false)

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
%%%%%%%%%%%%%%%%%%%%%%%
% rotational momentum %
%%%%%%%%%%%%%%%%%%%%%%%
% without correction for age, sex:
% [A,B,r,U,V,stats] = canoncorr(hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers),cogdata(~outliers,3:end));
% with correction:
[A,B,r,U,V,stats] = canoncorr(rotational_momentum_corrected,cogdata_corrected);
tmp = [B];
stats1 = stats;
stats1.r = r;
fig=setup_figure([],2,0.5)
bar(B)
xticklabels(cogdatalabels(3:end))
xtickangle(45)
ylabel('Coefficient strength')
title(sprintf('R=%0.2f, F(%d|%d)=%0.1f, p=%s', r, stats.df1, stats.df2, stats.F, stats.p))
save_figure([config.figdir 'figure4_correlations/4_CCA_rotational_momentum'])

%%%%%%%%%%%%%%
% cycle rate %
%%%%%%%%%%%%%%
% without correction for age, sex:
% [A,B,r,U,V,stats] = canoncorr(1./hmm_2ndlevel.cyctime_mu(~outliers),cogdata(~outliers,:));
% with correction:
[A,B,r,U,V,stats] = canoncorr(cycrate_mu_corrected,cogdata_corrected);
tmp = [tmp,B];
stats2 = stats;
stats2.r = r;
setup_figure([],2,0.5)
bar(B)
xticklabels(cogdatalabels(3:end))
xtickangle(45)
ylabel('Coefficient strength')
title(sprintf('R=%0.2f, F(%d|%d)=%0.1f, p=%s', r, stats.df1, stats.df2, stats.F, stats.p))
save_figure([config.figdir 'figure4_correlations/4_CCA_cycle_rate'])


setup_figure([],2,0.5)
bar(tmp)
xticklabels(cogdatalabels(3:end))
xtickangle(45)
ylabel('Coefficient strength')
legend({'Rotational momentum', 'Cycle rate'}, 'Location', 'southwest')
title({sprintf('Rotational momentum: R=%0.2f, p=%0.3f', stats1.r, stats1.p),...
  sprintf('Cycle rate: R=%s, p=%s', stats2.r, stats2.p)})
save_figure([config.figdir 'figure4_correlations/4_CCA_combined'])

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

