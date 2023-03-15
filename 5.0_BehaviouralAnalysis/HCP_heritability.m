whichstudy = 3;
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = colorscheme(whichstudy);

[~, ~, subj_age, subj_gender subj_RTs] = getparticipantinfo(whichstudy);
info = [subj_age,subj_gender+1,subj_RTs];

load(config.metricfile);

%% Get detailed HCP participants info

subjdata = readtable([config.participantcovariates, 'unrestricted_aquinn501_4_7_2017_9_4_13.csv']);
rstr_data = readtable([config.resultsdir, 'RESTRICTED_matsvanes_1_30_2023_9_38_27.csv']);
temp = readtable([config.participantcovariates, 'MEGfnames.csv']);
subj_ids = [];
for i=1:size(temp,1)
  subj_id = str2num(temp{i,1}{1}(7:12));
  if ~isempty(subj_id) && length(intersect(subj_id,subj_ids))==0
    subj_ids = [subj_ids;subj_id];
  end
end

inds = [];
for i=1:length(subj_ids)
  inds(i) = find(subjdata.Subject == subj_ids(i));
end
subjdata = subjdata(inds,:);
rstr_data = rstr_data(inds,:);

nSj=numel(inds);
InfMx = table;
InfMx.SubjectID(1:nSj,1) = 1:nSj;%rstr_data.Subject;
InfMx.MotherID = rstr_data.Mother_ID;
InfMx.FatherID = rstr_data.Father_ID;
Zygosity = [];
for k=1:nSj
    if isempty(rstr_data.ZygosityGT{k})
        if isempty(rstr_data.ZygositySR{k})
            Zygosity{k,1} = 'NotTwin';
        else
            Zygosity{k,1} = rstr_data.ZygositySR{k};
        end
    elseif strcmp(rstr_data.ZygosityGT{k}, 'DZ')
        Zygosity{k,1} = 'NotMZ';
    else
        Zygosity{k,1} = rstr_data.ZygosityGT{k};
    end
end
InfMx.Zygosity = Zygosity;
fname = [config.resultsdir, 'APACE_InfMx.csv'];
if ~exist(fname, 'file')
    writetable(InfMx, fname)
end

% also load more detailed data and align:

subjdata_detailed = readtable([config.participantcovariates, 'vars.txt']);
%headers = readtable('/Users/chiggins/data/HCPAnalysis/behav/column_headers.txt');
clear headers;
fid = fopen([config.participantcovariates, 'column_headers.txt']);
tline = fgetl(fid);
i=1;
while ischar(tline)
  temp = strrep(tline,' ','');
  temp = strrep(temp,'-','');
  headers{i,1} = temp;
  tline = fgetl(fid);i=i+1;
end
fclose(fid);

inds = [];
for i=1:length(subj_ids)
  inds(i) = find(subjdata_detailed.Var1 == subj_ids(i));
end
subjdata_detailed = subjdata_detailed(inds,:);
subjdata_detailed.Properties.VariableNames = headers;


%% Set up twin data
twindata = readmatrix([config.participantcovariates, 'twins.txt']);
twinmask = false(size(twindata,1),1);
for i=2:size(twindata)
  if ismember(twindata(1,i),subj_ids)
    twinmask(i) = true;
  end
end
twinstructure_MEG = twindata(twinmask,twinmask);
labels_twins = {'same','unrelated','monozygotic','dizygotic'};



%% Use APACE to determine heritability
measure={'cycle_rate', 'rotational_momentum', 'FO'};

for im = 1:3
  ACEfit_Par=[];
  ACEfit_Par.Model = 'ACE';
  if strcmp(measure{im}, 'cycle_rate')
      data = zscore(1000*1./hmm_2ndlevel.cycletime_mu_sess, [],'all');
      regr = [info(:,1:2)];
  elseif strcmp(measure{im}, 'rotational_momentum')
      data = zscore([hmm_1stlevel.tinda_per_ses{1}.cycle_metrics.rotational_momentum, hmm_1stlevel.tinda_per_ses{2}.cycle_metrics.rotational_momentum, hmm_1stlevel.tinda_per_ses{3}.cycle_metrics.rotational_momentum], [], 'all');
      regr = [info(:,1:2)];
  elseif strcmp(measure{im}, 'FO')
      data = zscore(hmm_2ndlevel.FO, [],'all');
      regr = [info(:,1:2), hmm_1stlevel.FO];
  end
  % regress out age and sex
  data_corrected = demean(regress_out(data, regr));
  ACEfit_Par.P_nm =data_corrected';
  ACEfit_Par.InfMx = [config.resultsdir, 'APACE_InfMx.csv'];
%     ACEfit_Par.InfMx = [config.resultsdir, 'heritability/twins.csv'];

  ACEfit_Par.ResDir = [config.resultsdir, 'heritability/', measure{im}, '/'];

  ACEfit_Par.Dsnmtx = [];
  ACEfit_Par.Nlz = 1;
  ACEfit_Par.AggNlz = 0;
  ACEfit_Par.ContSel = [];
  ACEfit_Par.NoImg = 1;
  
  ACEfit_Par = PrepData(ACEfit_Par)
  
  ACEfit_Par.alpha_CFT = [];
  ACEfit_Par = ACEfit(ACEfit_Par)
  
  ACEfit_Par.nPerm = 10000;
  ACEfit_Par.nBoot = 10000;
  ACEfit_Par.nParallel = 1;
  PrepParallel(ACEfit_Par);
  
  load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
  RunID = 1;
  ACEfit_Perm_Parallel(ACEfit_Par,RunID)
  ACEfit_Perm_Parallel_Results(ACEfit_Par)
  
  ACEfit_Results(ACEfit_Par)
    
  load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
  RunID = 1;
  ACEfit_Boot_Parallel(ACEfit_Par,RunID);
  
  load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
  ACEfit_Boot_Parallel_Results(ACEfit_Par);
  
  load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
  Boot_CIs(ACEfit_Par)
  
  if size(ACEfit_Par.P_nm,1)>1
      AgHe_Method(ACEfit_Par)
  end
  a=[];
  [a.Ests, a.CIs, a.Ps, a.Names] = APACEsummary(ACEfit_Par,'ResultSummary')
  
  heritability.(measure{im}) = a;
end

% save(config.metricfile, 'heritability', '-append')


%% plot correlations twins
[mono1,mono2] = find(twinstructure_MEG==1);
[dy1, dy2] = find(twinstructure_MEG==2);
[non1, non2] = find(twinstructure_MEG==0);


trait = 1./hmm_2ndlevel.cycletime_mu_sess'*1000;
figure; 
plot([trait(i); trait(j)], '.-', 'Color', [0 0 0 0.2]), xlim([0.5 2.5])
clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};





%% CCA with cognitive data
diego_labels = {'inhibitory control', 'processing speed', ...
  'cognitive flexibility', 'attention TN', 'reading', 'life satisfaction', ...
  'agreeableness', 'spatial orientation', 'emotional support', 'friendship', ...
  'positive affect', 'vocabulary', 'instrumental support', 'episodic memory', ...
  'mean purpose', 'fluid intelligence accuracy', 'extraversion', 'emotion recognition', ...
  'attention TP', 'fluid intelligence speed', 'verbal episodic memory', 'fear-affect',...
  'fear-somatic', 'anger-aggression', 'anger-hostility', 'perceived stress', ...
  'perceived hostility', 'sadness', 'anger-affect', 'loneliness', 'perceived rejection', ...
  'openness to experience', 'neuroticism', 'working memory', 'self-efficacy', 'conscientiousness'};

hcp_labels = {'Flanker_Unadj', 'ProcSpeed_Unadj', 'CardSort_Unadj', 'SCPT_TN',...
  'ReadEng_Unadj', 'LifeSatisf_Unadj', 'NEOFAC_A', 'VSPLOT_TC', 'EmotSupp_Unadj', ...
  'Friendship_Unadj', 'PosAffect_Unadj', 'PicVocab_Unadj', 'InstruSupp_Unadj', ...
  'PicSeq_Unadj', 'MeanPurp_Unadj', 'PMAT24_A_CR', 'NEOFAC_E', 'ER40_CR', ...
  'SCPT_TP', 'PMAT24_A_RTCR', 'IWRD_TOT', 'FearAffect_Unadj', 'FearSomat_Unadj', ...
  'AngAggr_Unadj', 'AngAggr_Unadj', 'PercStress_Unadj', 'PercHostil_Unadj', ...
  'Sadness_Unadj', 'AngAffect_Unadj', 'Loneliness_Unadj', 'PercReject_Unadj', ...
  'NEOFAC_O', 'NEOFAC_O', 'ListSort_Unadj', 'SelfEff_Unadj', 'NEOFAC_C'};% note: unsure whether VSPLOT is correct (true correct or reaction time?)




hcp_labels = [diego_labels; hcp_labels];

for k=1:size(hcp_labels,2)
  hcp_ix(k) = find(strcmp(headers, hcp_labels{2,k}));
end

% hcp_camcan_labels = hcp_ix([16, nan, 18, nan, 12, 5, nan, nan, nan, 14, 21, nan]);
hcp_cogix = hcp_ix([1:5, 8, 12, 14, 16, 19, 20, 21, 34]);
cogdata = table2array(subjdata_detailed(:, hcp_cogix));

rotational_momentum_corrected =  regress_out(hmm_1stlevel.cycle_metrics.rotational_momentum, regr);
[A,B,r,U,V,stats] = canoncorr(rotational_momentum_corrected, cogdata);
tmp = [B];
stats1 = stats;
stats1.r = r;
setup_figure([],2,0.5)
bar(B)
% xticklabels(hcp_labels(1,:))
% xtickangle(45)
ylabel('Coefficient strength')
title(sprintf('R=%0.2f, F(%d|%d)=%0.1f, p=%s', r, stats.df1, stats.df2, stats.F, stats.p))
save_figure([config.figdir 'figure4_correlations/4_CCA_rotational_momentum'])

%%%%%%%%%%%%%%
% cycle rate %
%%%%%%%%%%%%%%
cycle_rate_corrected =  regress_out(1./hmm_2ndlevel.cyctime_mu, regr);
[A,B,r,U,V,stats] = canoncorr(cycle_rate_corrected, cogdata);
tmp = [tmp,B];
stats2 = stats;
stats2.r = r;
setup_figure([],2,0.5)
bar(B)
% xticklabels(hcp_labels(1,:))
% xtickangle(45)
ylabel('Coefficient strength')
title(sprintf('R=%0.2f, F(%d|%d)=%0.1f, p=%s', r, stats.df1, stats.df2, stats.F, stats.p))
save_figure([config.figdir 'figure4_correlations/4_CCA_cycle_rate'])


setup_figure([],2,0.5)
bar(tmp)
% xticklabels(hcp_labels(1,:))
% xtickangle(45)
ylabel('Coefficient strength')
legend({'Rotational momentum', 'Cycle rate'}, 'Location', 'southwest')
title({sprintf('Rotational momentum: R=%0.2f, p=%0.3f', stats1.r, stats1.p),...
  sprintf('Cycle rate: R=%s, p=%s', stats2.r, stats2.p)})
save_figure([config.figdir 'figure4_correlations/4_CCA_combined'])















