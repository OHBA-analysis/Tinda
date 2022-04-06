% CamCan script to check for hierarchical variation with age:
% Part1: quantify age variation with lower level HMM states:

%% compare to single state FO:
% script called to load viterbi paths inferred and hmm objects and run
% post-hoc sequence analysis:
if ~exist('whichstudy','var')
    whichstudy = 4; % 4 denotes camcan data
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = set1_cols();


%% part 1: Viterbi path assymetry analysis
% first load data and plot basic temporal statistics:
temp = load(fullfile(config.hmmfolder,config.hmmfilename));

hmm = temp.hmm;
if ~isfield(hmm,'gamma') && whichstudy<4
    hmm.gamma = temp.Gamma;
    hmm.statepath = temp.vpath;
end
if whichstudy<3
    %load(config.prepdatafile,'hmmT','subj_inds');
    hmm = hmm_permutestates(hmm,temp.new_state_ordering);
    for i=1:config.nSj
        hmmT{i} = sum(hmm.subj_inds==i);
    end
elseif whichstudy==3 
    hmmT = temp.T_all;
    hmm.subj_inds = zeros(size(hmm.statepath));
    t_offset = 0;
    for i=1:length(hmmT)
        t_length = sum(hmmT{i}) - length(hmmT{i})*(length(hmm.train.embeddedlags)-1);
        hmm.subj_inds(t_offset + [1:t_length]) = i;
        t_offset = t_offset + t_length;
    end
    hmm.subj_inds = ceil(hmm.subj_inds/3); % account for multiple runs per subj
    hmmTold = reshape(hmmT,3,config.nSj);
    hmmTsubj = cell(config.nSj);
    for i=1:config.nSj
        hmmTsubj{i} = [hmmTold{1,i},hmmTold{2,i},hmmTold{3,i}]-(length(hmm.train.embeddedlags)-1);
    end
    hmmT = hmmTsubj;
    clear hmmTsubj hmmTold;
elseif whichstudy==4
    hmmT = temp.T_all;
    % correct for embedded lags:
    for i=1:length(hmmT)    
        hmmT{i} = hmmT{i} - (length(hmm.train.embeddedlags)-1);
    end
    load(config.matfilelist);
end
clear temp vpath;

if whichstudy<4
    Gamma = hmm.gamma;
end
K = hmm.K;
FO = zeros(K,K,2,config.nSj);
opts = [];
opts.K = 12;
opts.Fs = config.sample_rate;
for subnum=1:config.nSj
    if whichstudy~=4
        vpath{subnum} = hmm.statepath(hmm.subj_inds==subnum);
    else
        temp = load(mat_files_orth{subnum},'vpath');
        vpath{subnum} = temp.vpath;
    end
    if simtests
        try % note unusual syntax here is just to catch very rare precision errors in numeric simulation
            vpath{subnum} = simulateVpath(vpath{subnum},hmmT{subnum},hmm.K);
        catch
            vpath{subnum} = simulateVpath(vpath{subnum},hmmT{subnum},hmm.K);
        end
    end
    LT = getStateLifeTimes(vpath{subnum},hmmT{subnum},opts);
    LTmerged(subnum,:) = cellfun(@mean,LT);
    FracOcc(subnum,:) = getFractionalOccupancy(vpath{subnum},sum(hmmT{subnum}),opts);
    IT = getStateIntervalTimes(vpath{subnum},hmmT{subnum},opts);
    ITmerged(subnum,:) = cellfun(@mean,IT);
end
info = camcan_getparticipantinfo(config);


%% now check for age correlation:

% this now moved to 'allstudies_correlations' as it is generic to all the
% studies covered (ie looking at first level hmm correlations with age and
% gender)

%% PART TWO: hierarchical structure

%%% note much of this analysis is replicated in allstudies_correlations,
%%% but keep here for the ticvmm results...


load(config.secondlevelmodelfile)
windowlength = 125;

%%
if ~contains(config.Poiss_dir,'overlappingWindows')
    figdir = [config.figdir,'4_covariates_W',in2str(W),'/'];
    mkdir(figdir);
    clear FO 
    for subnum=1:600
        Gamtemp = hmmPoiss.gamma(Poiss_subj_inds==subnum,:);
        cycletimes = getStateIntervalTimes(Gamtemp,length(Gamtemp));
        cycletime_mu(subnum,:) = cellfun(@mean,cycletimes);
        cycletime_std(subnum,:) = cellfun(@std,cycletimes);
        cycletime_med(subnum,:) = cellfun(@median,cycletimes);
        FO(subnum,:) = getFractionalOccupancy(Gamtemp,length(Gamtemp));
        cyctimes{subnum} = cycletimes;
        lifetimes{subnum} = getStateLifeTimes(Gamtemp,length(Gamtemp));
    end
    
else
    figdir = [config.figdir,'4_covariates_W',int2str(windowlength),'_overlappingWindows/'];
    mkdir(figdir);
    clear cycletimes cycletime_mu cycletime_std cycletime_med FO cyctimes lifetimes cycletime_mu_min cycletime_med_min cycletime_std_min
    load([config.Poiss_dir,'filelist.mat'])
    samp_2minute = config.sample_rate*2*60;
    for subnum=1:600
        fprintf(['\nSubj: ',int2str(subnum)]);
        load(mat_files_poiss{subnum},'Gamma','T');
        %cycletimes = getStateIntervalTimes(Gamma,T);
        cycletimes{1} = getStateCycleTimes(Gamma,T);
        cycletime_mu(subnum,:) = cellfun(@mean,cycletimes)./config.sample_rate;
        cycletime_std(subnum,:) = cellfun(@std,cycletimes)./config.sample_rate;
        cycletime_med(subnum,:) = cellfun(@median,cycletimes)./config.sample_rate;
        FO(subnum,:) = getFractionalOccupancy(Gamma,length(Gamma))./config.sample_rate;
        cyctimes{subnum} = cycletimes;
        lifetimes{subnum} = getStateLifeTimes(Gamma,T);
        % also get minute by minute detail:
        for imin=1:5
            startseg = find(cumsum(T)>(imin-1)*samp_2minute,1)-1;
            endseg = find(cumsum(T)>(imin)*samp_2minute,1)-1;
            if ~isempty(startseg) && ~isempty(endseg)
                if startseg==endseg
                    T_sub = samp_2minute;
                else
                    T_sub = sum(T(1:startseg+1)) - (imin-1)*samp_2minute;
                    T_sub = [T_sub;T(startseg+2:endseg)];
                    T_sub = [T_sub;samp_2minute - sum(T_sub)];
                end
                try
                    Gamtemp = Gamma((imin-1)*samp_2minute + [1:samp_2minute],:);
                catch
                    Gamtemp = Gamma((imin-1)*samp_2minute :end,:);
                    T_sub = T_sub(cumsum(T_sub)<length(Gamtemp));
                    T_sub = [T_sub;length(Gamtemp)-sum(T_sub)];
                end
                temp = getStateCycleTimes(Gamtemp,T_sub);
                cyctime{subnum,imin} = temp;
                cycletime_mu_min(subnum,imin) = mean(temp);
                cycletime_med_min(subnum,imin) = median(temp);
                cycletime_std_min(subnum,imin) = std(temp);
            end
        end 
    end
end


%% check for correlation:
outliers = abs(cycletime_mu(:,1)-mean(cycletime_mu(:,1))) > 2*std(cycletime_mu(:,1));
    
for ireg=[1]%,3,4];
    
figure('Position',[7 409 993 389]);
subplot(1,3,1);
scatter(info(~outliers,ireg),cycletime_mu(~outliers,1),'filled')
[R,P] = corrcoef(info(~outliers,ireg),cycletime_mu(~outliers,1));
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper(reglabels{ireg},'Mean Cycle Time');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

subplot(1,3,2);
scatter(info(~outliers,ireg),cycletime_med(~outliers,1),'filled');
[R,P] = corrcoef(info(~outliers,ireg),cycletime_med(~outliers,1));
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper(reglabels{ireg},'Median Cycle Time');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

subplot(1,3,3);
scatter(info(~outliers,ireg),cycletime_std(~outliers,1),'filled')
[R,P] = corrcoef(info(~outliers,ireg),cycletime_std(~outliers,1));
title(['STD: rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper(reglabels{ireg},'Cycle Time std');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

print([figdir '4A_metastateCycleCorr_',reglabels{ireg}],'-depsc')
end
%% check for metatstate FO correlation:
for ireg = [1,3,4];
figure('Position',[6 435 1289 363]);
for k =1:3
    subplot(1,3,k);
    scatter(info(~outliers,ireg),FO(~outliers,k),'filled')
    [R,P] = corrcoef(info(~outliers,ireg),FO(~outliers,k));
    title(['Metastate ',int2str(k),', rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
    plot4paper(reglabels{ireg},['Metastate ',int2str(k),' FO']);
    H=lsline(gca)
set(H,'LineWidth',2);
axis square;
end
print([figdir '4B_metastateFOCorr_',reglabels{ireg}],'-depsc')
end

%% side check: is this structure significant after regressing out all effect of lower level FO?

B = pinv(FracOcc)*FO;
FOmeta_resid = FO - FracOcc*B;
for ireg = [1,3,4];
figure('Position',[6 435 1289 363]);
for k =1:3
    subplot(1,3,k);
    scatter(info(~outliers,ireg),FOmeta_resid(~outliers,k),'filled')
    [R,P] = corrcoef(info(~outliers,ireg),FOmeta_resid(~outliers,k));
    title(['Metastate ',int2str(k),', rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
    plot4paper(reglabels{ireg},['Metastate ',int2str(k),' FO']);
    axis square;
end
print([figdir '4B_metastateFOCorr_regressout_',reglabels{ireg}],'-depsc')
end
% note these remain significant as function of age even after lower level
% FO is regressed out - suggests the hierarchical structure contains more
% info than just that contained in the lower level HMM...

%% load cognitive data:

[cogdata,cogdatalabels] = camcan_getCognitiveData(config);

% remove subjID: 
cogdata = cogdata(:,1:15);
cogdatalabels = cogdatalabels(1:15);
% run multiple regression - can I predict cycle times from these data?

[B,~,~,~,stats] = regress(1./cycletime_mu(~outliers),[ones(sum(~outliers),1),normalise(cogdata(~outliers,:))]);
%
[B,~,~,~,stats] = regress(cycletime_mu(~outliers),[ones(sum(~outliers),1),normalise(cogdata(~outliers,:)),normalise(info(~outliers,1))]);

% what if we first regress out age effects?
[B,~,~,~,stats] = regress((cycletime_mu(~outliers)),[ones(sum(~outliers),1),normalise(info(~outliers,[1,3,4]))]);
[B,~,~,~,stats] = regress((cycletime_mu(~outliers)) - [ones(sum(~outliers),1),normalise(info(~outliers,[1,3,4]))]*B,[ones(sum(~outliers),1),normalise(cogdata(~outliers,3:end))]);

stats(3)

%% try non parametric instead of f-tests:
[B,~,~,~,stats] = regress(cycletime_mu(~outliers),[ones(sum(~outliers),1),normalise(info(~outliers,1))]);
y = cycletime_mu(~outliers) - [ones(sum(~outliers),1),normalise(info(~outliers,1))]*B;
X = [ones(sum(~outliers),1),normalise(cogdata(~outliers,:))];
beta = regress(y,X);
ssq_true = sum(sum((y-X*beta).^2));
n_perms = 1000;
for iperm=1:n_perms
    perms = randperm(length(y));
    beta = regress(y,X(perms,:));
    ssq_null(iperm) = sum(sum((y-X(perms,:)*beta).^2));
end
figure();
hist(ssq_null,50);
line([ssq_true,ssq_true],ylim(gca));
nonparam_pvalue = 1 - sum(ssq_true<ssq_null)/n_perms
%%
y = cycletime_mu(~outliers);
X = [ones(sum(~outliers),1),info(~outliers,1)];
beta = regress(y,X);
ssq_true = sum(sum((y-X*beta).^2));
n_perms = 1000;clear ssq_null
for iperm=1:n_perms
    perms = randperm(length(y));
    beta = regress(y,X(perms,:));
    ssq_null(iperm) = sum(sum((y-X(perms,:)*beta).^2));
end
figure();
hist(ssq_null,50);
line([ssq_true,ssq_true],ylim(gca));
nonparam_pvalue = 1 - sum(ssq_true<ssq_null)/n_perms

%% try regression going other way:

[B,~,~,~,stats] = regress(cycletime_mu(~outliers),[ones(sum(~outliers),1),normalise(info(~outliers,1))]);
X = cycletime_mu(~outliers) - [ones(sum(~outliers),1),normalise(info(~outliers,1))]*B;
clear pvals_predict
for i=1:15
    y = normalise(cogdata(~outliers,i));
    [beta,~,~,~,stats] = regress(y,[ones(sum(~outliers),1),X]);
    pvals_predict(i) = stats(3);
end

%% fit RoniTibon's CCA model with one mode:

X = [ones(sum(~outliers),1),normalise([cogdata(~outliers,:)])];
y = 1./cycletime_mu(~outliers);

%B_control = regress(y,[ones(sum(~outliers),1),FracOcc(~outliers,1:end-1)]);
%y = y-FracOcc(~outliers,:)*B_control;

[B,~,~,~,stats] = regress(y,X);

stats(3)
figure();bar(B(2:end))
set(gca,'XTick',1:size(X,2)-1);
%set(gca,'XTickLabel',fullmodellabels{1});
xticklabels(cogdatalabels);
set(gca,'XTickLabelRotation',45);
plot4paper('','Regressor strength');

title(['Full model: p=',num2str(stats(3))]);
print([figdir,'4D_fullregressionmodel'],'-dpng');

%% try cardio measures correlation:
fname = [config.participantcovariates, '/CardioMeasures_summary.xlsx'];
headerlines = 0;
DATA = importdata(fname,' ',headerlines);
[info,subjid] = camcan_getparticipantinfo(config);
HRmean = zeros(600,3);
for i=1:600
    thissub = find(contains(DATA.textdata(:,1),subjid{i}(5:end)));
    if isempty(thissub)
        HRmean(i,:) = NaN;
    else
        HRmean(i,:) = DATA.data(thissub-1,17:19);
    end
end
labels = DATA.textdata(1,19:21);
labels = strrep(labels,'_',' ');

figure('Position',[1 408 999 390]);
for k=1:3
    subplot(1,3,k)
    temp = HRmean(:,k);
    outliers = isnan(temp) | abs(temp-nanmean(temp))>2*nanstd(temp) | abs(cycletime_mu-nanmean(cycletime_mu))>1.5*nanstd(cycletime_mu);
    scatter(temp(~outliers),1./cycletime_mu(~outliers,1),'filled')
    [R,P] = corrcoef(temp(~outliers),1./cycletime_mu(~outliers,1));
    title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
    plot4paper(labels{k},'Cycle Time');
    H=lsline(gca)
    set(H,'LineWidth',2);
    axis square;
end
print([figdir,'4E_Cardiac'],'-dpng');
%% load simple reaction time measure:
fname = [config.participantcovariates, 'RTsimple_summary.xlsx'];
headerlines = 0;
DATA = importdata(fname,' ',headerlines);
[info,subjid] = camcan_getparticipantinfo(config);
for i=1:600
    thissub = find(contains(DATA.textdata(:,1),subjid{i}(5:end)));
    if isempty(thissub)
        RTmean(i,1) = NaN;
    else
        RTmean(i,1) = DATA.data(thissub-1,8);
    end
end
RT_outliers = outliers | isnan(RTmean);% | RTmean>0.7 |cycletime_mu>0.5;
figure();
scatter(1./RTmean(~RT_outliers),1./cycletime_mu(~RT_outliers,1),'filled')
[R,P] = corrcoef(1./RTmean(~RT_outliers),1./cycletime_mu(~RT_outliers,1));
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper('Reaction time','Cycle Time');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;
print([figdir,'4E_RT'],'-dpng');
%% Fit one big regression model, then reduce param by param
X = [cogdata,info(:,[2,4]),HRmean,RTmean];

missing_data = any(isnan(X),2);
fullmodel_outliers = outliers | missing_data;

X = X(~fullmodel_outliers,:);
X = [ones(length(X),1),normalise(X)];

[B,~,~,~,stats] = regress(1./cycletime_mu(~fullmodel_outliers),X);
fullmodellabels = {cogdatalabels{:},'handedness','ticvmm',labels{:},'RT'};
figure();bar(B(2:end))
set(gca,'XTick',1:size(X,2)-1);
%set(gca,'XTickLabel',fullmodellabels{1});
xticklabels(fullmodellabels);
set(gca,'XTickLabelRotation',45);
plot4paper('','Regressor strength');

title(['Full model: p=',num2str(stats(3))]);
%print([figdir,'4D_fullregressionmodel'],'-dpng');


%% eliminate individual groups:
figure('Position',[1 323 1315 475]);
groups = {[1:2,16],[17:20],[3:15,21]};
groupnames = {'Age/handedness','Physiolog','Cognitive'}
for i=1:3
    control_var = setdiff(1:size(X,2),groups{i}+1);
    [B] = regress(1./cycletime_mu(~fullmodel_outliers),X(:,control_var));
    y = 1./cycletime_mu(~fullmodel_outliers) - X(:,control_var)*B;
    [B,~,~,~,stats] = regress(y,X(:,[1,groups{i}+1]));
    
    subplot(1,3,i);
    bar(B(2:end));
    set(gca,'XTick',1:size(X,2)-1);
    %set(gca,'XTickLabel',fullmodellabels{1});
    xticklabels(fullmodellabels(groups{i}));
    set(gca,'XTickLabelRotation',45);
    plot4paper('','Regressor strength');

    title([groupnames{i},'model: p=',num2str(stats(3))]);
end
print([figdir,'4D_reducedregressionmodel'],'-dpng');
%% evaluate temporal patterns over course of scan:

% needs more restructive outlier removal:
to_remove = abs(cycletime_mu_min-nanmean(cycletime_mu_min))>2*nanstd(cycletime_mu_min);
to_remove = any(to_remove(:,1:4),2);
to_remove = to_remove | any(cycletime_mu_min(:,1:4)==0,2) | any(isnan(cycletime_mu_min(:,1:4)),2);
x = (cycletime_mu_min(~to_remove,1:4)./config.sample_rate);

clear p;
p(1) = anova1(x);
p(2) = anova1(x.^-1);

figure('Position',[440 445 1001 353]);
subplot(1,2,1);
distributionPlot(x,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}},'globalNorm',1);
set(gca,'FontSize',fontsize)
title(['Cycletime vs minutes in scan: p=',num2str(p(1),2)]);plot4paper('Minutes in scan','Cycle time');grid on;
set(gca,'XTickLabel',[1:5]*2);
ylim([0,5]);

x = x.^-1;


subplot(1,2,2);
distributionPlot(x,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}},'globalNorm',1);
set(gca,'FontSize',fontsize)
title(['Cyclerate vs minutes in scan: p=',num2str(p(2),2)]);plot4paper('Minutes in scan','Cycle rate');grid on;
set(gca,'XTickLabel',[1:5]*2);
ylim([0,1]);
print([figdir,'4F_cyclesperminute'],'-dpng');
%% how consistent are cycletimes from same vs different subjects?
clear D


D(:,1) = (x(:,1)-x(:,4)).^2;
shuffles = randperm(length(x));
D(:,2) = (x(:,1)-x(shuffles,4)).^2;
P = anova1((sqrt(D)))
ax = gca;
ax.LineWidth=2;
plot4paper('','first minute minus last minute cycle rate');
set(gca,'XTickLabels',{'Same subject','Different subjects'})
title(['ANOVA: p=',num2str(P,3)])
print([figdir,'4F_cyclerateconsistency'],'-dpng');

%% now check for structure over consecutive cycles:
LT_conseccorr = zeros(600,3);
cyc_conseccorr = zeros(600,1);
for subnum=1:600
    for k=1:3
        R = corrcoef(lifetimes{subnum}{k}(1:end-1),lifetimes{subnum}{k}(2:end));
        LT_conseccorr(subnum,k) = R(2,1);
        
        
    end
    R = corrcoef(1./cyctimes{subnum}{1}(1:end-1),1./cyctimes{subnum}{1}(2:end));
    cyc_conseccorr(subnum,1) = R(2,1);
end

% fisher zscore for hypothesis tests:
z_LT = atanh(LT_conseccorr);
z_cyc = atanh(cyc_conseccorr);
for k=1:3
    [~,p_LT(k)] = ttest(z_LT(:,k));
    
end
[~,p_cyc] = ttest(z_cyc);

figure('Position',[7 409 993 389]);
subplot(1,2,1);
distributionPlot(LT_conseccorr,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
set(gca,'FontSize',fontsize)
title(['Consecutive Lifetime Correlation: p<',num2str(max(p_LT))]);plot4paper('MetaState','Correlation coefficient');grid on;

subplot(1,2,2);
distributionPlot(cyc_conseccorr,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
set(gca,'FontSize',fontsize)
title(['Consecutive cyclerate Correlation: p=',num2str(p_cyc)]);plot4paper('MetaState','Correlation coefficient');grid on;
%print([config.figdir '0_temporalstats_FO'],'-depsc')
print([figdir '4F_metastate_consec_corr'],'-depsc')

%% test for correlations with age:
for ireg=1
    for k=1:3
        [R,P] = corrcoef(info(:,ireg),LT_conseccorr(:,k));
        p_consec(ireg,k,1) = P(2,1);
        
    end
    [R,P] = corrcoef(info(:,ireg),cyc_conseccorr(:,1));
        p_consec(ireg,1,2) = P(2,1);
end

%% how does cycle length relate to mean state duration, mean interval time?
figure();
subplot(1,2,1);
[R,P] = corrcoef(mean(LTmerged(~outliers,:),2),cycletime_mu(~outliers,:));
scatter(mean(LTmerged(~outliers,:),2),cycletime_mu(~outliers,1),'filled')
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper('Mean State Lifetime','Cycle Time');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

subplot(1,2,2);
[R,P] = corrcoef(1./mean(LTmerged(~outliers,:),2),1./cycletime_mu(~outliers,:));
scatter(1./mean(LTmerged(~outliers,:),2),1./cycletime_mu(~outliers,1),'filled')
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper('Mean State Lifetime','Cycle Time');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

%% load simple reaction time measure:
% fname = ['/Users/chiggins/Documents/CamCan/RTsimple_summary_amended.txt']
% [fid,errmsg] = fopen(fname,'r','l','UTF16-LE');
% 
% [RT] = textscan(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%s\n');
% %RT = textscan(fid)
% fclose(fid);
% RT



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%BELOW IS REDUNDANT%%%%%%%%%%%%%%%
% 
% %% now check for structure over consecutive cycles:
% LT_conseccorr = zeros(600,3);
% cyc_conseccorr = zeros(600,3);
% for subnum=1:600
%     for k=1:3
%         R = corrcoef(lifetimes{subnum}{k}(1:end-1),lifetimes{subnum}{k}(2:end));
%         LT_conseccorr(subnum,k) = R(2,1);
%         
%         R = corrcoef(cyctimes{subnum}{k}(1:end-1),cyctimes{subnum}{k}(2:end));
%         cyc_conseccorr(subnum,k) = R(2,1);
%     end
% end
% 
% % fisher zscore for hypothesis tests:
% z_LT = atanh(LT_conseccorr);
% z_cyc = atanh(cyc_conseccorr);
% for k=1:3
%     [~,p_LT(k)] = ttest(z_LT(:,k));
%     [~,p_cyc(k)] = ttest(z_cyc(:,k));
% end
% 
% figure('Position',[7 409 993 389]);
% subplot(1,2,1);
% distributionPlot(LT_conseccorr,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
% set(gca,'FontSize',fontsize)
% title('Consecutive Lifetime Correlation');plot4paper('MetaState','Correlation coefficient');grid on;
% 
% subplot(1,2,2);
% distributionPlot(cyc_conseccorr,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
% set(gca,'FontSize',fontsize)
% title('Consecutive cycletime Correlation');plot4paper('MetaState','Correlation coefficient');grid on;
% %print([config.figdir '0_temporalstats_FO'],'-depsc')
% print([figdir '4_metastate_consec_corr'],'-depsc')
% 
% % test for correlations with age:
% for ireg=1
%     for k=1:3
%         [R,P] = corrcoef(info(:,ireg),LT_conseccorr(:,k));
%         p_consec(ireg,k,1) = P(2,1);
%         [R,P] = corrcoef(info(:,ireg),cyc_conseccorr(:,k));
%         p_consec(ireg,k,2) = P(2,1);
%     end
% end
% % only very weakly significant
% 
% %% now check for structure in cycle times within scan (by minutes):
% windowlength = 17;
% samp_2minute = 2*round(250*60/windowlength);
% warning off;
% for subnum=1:600
%     Gamtemp = hmmPoiss.gamma(Poiss_subj_inds==subnum,:);
%     T = length(Gamtemp);
%     T_sub = [repmat(samp_2minute,1,floor(T/samp_2minute)),mod(T,samp_2minute)];
%     for imin=1:length(T_sub)
%         temp = getStateIntervalTimes(Gamtemp(sum(T_sub(1:imin-1))+[1:T_sub(imin)],:),T_sub(imin));
%         cyctime{subnum,imin} = temp{1};
%         cycletime_mu_min(subnum,imin) = mean(cellfun(@mean,temp));
%         cycletime_med_min(subnum,imin) = mean(cellfun(@median,temp));
%         cycletime_std_min(subnum,:) = sum(cellfun(@std,temp));
%     end 
% end
% warning on;
% 
% %% plot:
% figure();
% distributionPlot(cycletime_mu_min(~outlierss,1:5),'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
% set(gca,'FontSize',fontsize)
% title('Cycletime vs minutes in scan');plot4paper('Minutes in scan','Cycle time');grid on;
% %print([config.figdir '0_temporalstats_FO'],'-depsc')
% 
% 
% 
% 
% %% quick check for outliers effect:
% ireg=1;
% figure('Position',[7 409 993 389]);
% outliers = abs(cycletime_mu(:,1)-mean(cycletime_mu(:,1)))>1.5*std(cycletime_mu(:,1));
% subplot(1,3,1);
% scatter(info(~outliers,ireg),log(cycletime_mu(~outliers,1)),'filled')
% [R,P] = corrcoef(info(~outliers,ireg),log(cycletime_mu(~outliers,1)));
% title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
% plot4paper(reglabels{ireg},'Mean Cycle Time');
% H=lsline(gca)
% set(H,'LineWidth',2);
% axis square;
% subplot(1,3,2);
% outliers = abs(cycletime_med(:,1)-mean(cycletime_med(:,1)))>2*std(cycletime_med(:,1));
% 
% scatter(info(~outliers,ireg),cycletime_med(~outliers,1),'filled');
% [R,P] = corrcoef(info(~outliers,ireg),cycletime_med(~outliers,1));
% title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
% plot4paper(reglabels{ireg},'Median Cycle Time');
% H=lsline(gca)
% set(H,'LineWidth',2);
% axis square;
% 
% subplot(1,3,3);
% scatter(info(:,ireg),cycletime_std(:,1),'filled')
% [R,P] = corrcoef(info(:,ireg),cycletime_std(:,1));
% title(['STD: rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
% plot4paper(reglabels{ireg},'Cycle Time std');
% H=lsline(gca)
% set(H,'LineWidth',2);
% axis square;
% 
