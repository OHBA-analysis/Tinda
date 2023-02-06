whichstudy = 4; % 1 denotes the hmm model run in Higgins2020_neuron
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

for k1=1:2
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
regr = cogdata(:,1:2);
cycrate_mu_corrected =  regress_out(1./hmm_2ndlevel.cyctime_mu, regr); 
rotational_momentum_corrected =  regress_out(hmm_1stlevel.cycle_metrics.rotational_momentum, regr);   
FO_corrected =  regress_out(hmm_2ndlevel.FO, regr);   
cogdata_corrected =  regress_out(cogdata(:,3:end), regr);   

%%%%%%%%%%%%%%%%%%%%%%%
% rotational momentum %
%%%%%%%%%%%%%%%%%%%%%%%
% without correction for age, sex:
% [A,B,r,U,V,stats] = canoncorr(hmm_1stlevel.cycle_metrics.rotational_momentum(~outliers),cogdata(~outliers,3:end));
% with correction:
[A,B,r,U,V,stats] = canoncorr(rotational_momentum_corrected(~outliers,1),cogdata_corrected(~outliers,:));
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
[A,B,r,U,V,stats] = canoncorr(cycrate_mu_corrected(~outliers,1),cogdata_corrected(~outliers,:));
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



% Do the same but with Bayesian Partial Least Squares (Vidaurre 2013) as in
% Vidaurre PNAS 2018






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

