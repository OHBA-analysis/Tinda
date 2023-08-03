% generic correlations to plot over all studies:

if ~exist('whichstudy','var')
  whichstudy = 1; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = colorscheme(whichstudy);

% include option for simulations here - results to be saved in adjacent
% folder:
simtests = false;
if simtests
  config.figdir = strrep(config.figdir,['Study',int2str(whichstudy)],['Study',int2str(whichstudy),'_simtest']);
  mkdir(config.figdir);
end

if whichstudy==1
  [info,varnames] = MEGUK_getparticipantinfo();
  subj_age = info(:,11);
  subj_gender = info(:,12); % 1 denotes male, 2 denotes female
  subj_RTs = info(:,3);
elseif whichstudy==3
  [info,varnames] = HCP_getparticipantinfo(config);
  subj_age = info(:,4);
  subj_gender = info(:,3); % 1 denotes female
  subj_RTs = info(:,end);
elseif whichstudy==4
  info = camcan_getparticipantinfo(config);
  subj_age = info(:,1);
  subj_gender = info(:,3); % 1 denotes male, 2 female
  subj_RTs = info(:,5);
  
end
info = [subj_age,subj_gender,subj_RTs];
%% start with standard correlation analyses to be replicated for all studies:

load(config.metricfile);
%% Plot first level HMM statsitics correlation with age and gender:
fontsize = 18;
for ireg = 1:2
  if ireg==1
    reglabel = 'Age';
  else
    reglabel = 'Gender';
  end
  
  %% FRACTIONAL OCCUPANCY
  setup_figure([],1,2);%'Position',[588 63 412 735]);
  subplot(211);
  distributionPlot(hmm_1stlevel.FO,'showMM',2,'color',{color_scheme{1:size(hmm_1stlevel.FO,2)}});
  set(gca,'YLim',[0 1.1*prctile(hmm_1stlevel.FO(:),97.5)])
  title('Fractional Occupancy');xlabel('RSN-State'), ylabel('Proportion');grid on;
  subplot(212);
  pvals = [0.05, 0.01, 0.001];
  pvals_line = {':k', '--k', '-.k'};
  for k=1:12
    to_keep = ~isnan(info(:,ireg));
    [R,P] = corrcoef(hmm_1stlevel.FO(to_keep,k),info(to_keep,ireg));
    rho_FO(k) = R(2,1);
    pvals_FO(k) = P(2,1);
    if pvals_FO(k)>(0.05/12)
      b = bar(k,rho_FO(k));hold on;
      set(b,'FaceColor',color_scheme{k});
    else
      b = bar(k,rho_FO(k));
      hold on;
      set(b,'FaceColor',color_scheme{k});
      set(b,'LineWidth',2);
    end
  end
  
  % add significance asterisks:
  yl = ylim;
  for k=1:12
    sigstar({k}, 12*pvals_FO(k));
  end
  set(gca,'XTick',1:12);
  xlabel('State'), ylabel('Correlation');
%   print([config.figdir '0_temporalstats_FO_',reglabel],'-depsc')
  
%% LIFE TIMES
  figure('Position',[588 63 412 735]);
  subplot(211);
  distributionPlot(hmm_1stlevel.LT_mu ./ config.sample_rate * 1000,'showMM',2,'color',{color_scheme{1:size(hmm_1stlevel.FO,2)}})
  title('Life Times');plot4paper('RSN-State','Time (ms)');grid on;
  YL = 1.1*prctile(hmm_1stlevel.LT_mu(:),97.5)./ config.sample_rate * 1000;
  set(gca,'YLim',[0 YL],'FontSize',fontsize,'FontSize',fontsize);
  subplot(212);
  for k=1:12
    to_keep = ~isnan(info(:,ireg));
    [R,P] = corrcoef(hmm_1stlevel.LT_mu(to_keep,k),info(to_keep,ireg));
    rho_LT(k) = R(2,1);
    pvals_LT(k) = P(2,1);
    if pvals_LT(k)>(0.05/12)
      b = bar(k,rho_LT(k));hold on;
      set(b,'FaceColor',color_scheme{k});
    else
      b = bar(k,rho_LT(k));
      hold on;
      set(b,'FaceColor',color_scheme{k});
      set(b,'LineWidth',2);
    end
  end
  % add significance asterisks:
  for k=1:12
    sigstar({k}, 12*pvals_LT(k));
  end
  set(gca,'XTick',1:12);
  xlabel('State'), ylabel('Correlation');
%   print([config.figdir '0_temporalstats_LT_',reglabel],'-depsc')
  
%% INTERVAL TIMES
  figure('Position',[588 63 412 735]);
  subplot(211);
  distributionPlot(hmm_1stlevel.IT_mu ./ config.sample_rate * 1000,'showMM',2,'color',{color_scheme{1:size(hmm_1stlevel.FO,2)}})
  title('Interval Times');plot4paper('RSN-State','Time (ms)');grid on;
  YL = 1.1*prctile(hmm_1stlevel.IT_mu(:),97.5)./ config.sample_rate * 1000;
  set(gca,'YLim',[0 YL],'FontSize',fontsize,'FontSize',fontsize);
  subplot(212);
  for k=1:12
    to_keep = ~isnan(info(:,ireg));
    [R,P] = corrcoef(hmm_1stlevel.IT_mu(to_keep,k),info(to_keep,ireg));
    rho_IT(k) = R(2,1);
    pvals_IT(k) = P(2,1);
    if pvals_IT(k)>(0.05/12)
      b = bar(k,rho_IT(k));hold on;
      set(b,'FaceColor',color_scheme{k});
    else
      b = bar(k,rho_IT(k));
      hold on;
      set(b,'FaceColor',color_scheme{k});
      set(b,'LineWidth',2);
    end
  end
  % add significance asterisks:
  yl = ylim;
  for k=1:12
    sigstar({k}, 12*pvals_IT(k));
  end
  set(gca,'XTick',1:12);
  xlabel('State'), ylabel('Correlation');
%   print([config.figdir '0_temporalstats_IT_',reglabel],'-depsc')
  
  % correlation with all lifetimes:
  figure();
  subplot(2,2,1);
  to_keep = ~isnan(info(:,ireg));
  scatter(info(to_keep,ireg),mean(hmm_1stlevel.LT_mu(to_keep,:),2),'filled')
  [R,P] = corrcoef(info(to_keep,ireg),mean(hmm_1stlevel.LT_mu(to_keep,:),2));
  title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
  plot4paper(reglabel,'Mean State LT');
  H = lsline(gca);
  set(H,'LineWidth',2);
  
  subplot(2,2,2);
  scatter(info(to_keep,ireg),median(hmm_1stlevel.LT_med(to_keep,:),2),'filled')
  [R,P] = corrcoef(info(to_keep,ireg),median(hmm_1stlevel.LT_med(to_keep,:),2));
  title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
  plot4paper(reglabel,'Median State LT');
  H = lsline(gca);
  set(H,'LineWidth',2);
  
  subplot(2,2,3);
  scatter(info(to_keep,ireg),mean(hmm_1stlevel.IT_mu(to_keep,:),2),'filled')
  [R,P] = corrcoef(info(to_keep,ireg),mean(hmm_1stlevel.IT_mu(to_keep,:),2));
  title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
  plot4paper(reglabel,'Mean State IT');
  H = lsline(gca);
  set(H,'LineWidth',2);
  
  subplot(2,2,4);
  scatter(info(to_keep,ireg),median(hmm_1stlevel.IT_med(to_keep,:),2),'filled')
  [R,P] = corrcoef(info(to_keep,ireg),median(hmm_1stlevel.IT_med(to_keep,:),2));
  title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
  plot4paper(reglabel,'Median State IT');
  H = lsline(gca);
  set(H,'LineWidth',2);
%   print([config.figdir '0_temporalstats_aggregate_',reglabel],'-depsc')
  
end

%% plot population clustering structure:

figure('Position',[1 378 999 420]);
subplot(1,3,1);
imagesc(corrcoef(hmm_1stlevel.FO))
plot4paper('State FO','State FO');
title('FO correlation over subjects');
mkdir([config.figdir 'supp/'])
axis square;
colorbar;

subplot(1,3,2);
imagesc(corrcoef(hmm_1stlevel.LT_mu))
plot4paper('State LT','State LT');
title('LT correlation over subjects');
axis square;
colorbar;

subplot(1,3,3);
imagesc(corrcoef(hmm_1stlevel.IT_mu))
plot4paper('State IT','State IT');
title('IT correlation over subjects');
axis square;
colorbar;
% print([config.figdir 'supp/','1_populationFOcorrelation'],'-dpng')

% extract metric to show which of the two clusters a participant falls
% into:
[a,b] = pca(hmm_1stlevel.FO,'NumComponents',1);
if sum(a(5:9))>0
  a = a*-1;
  b = b*-1;
end
highpowersubjs = b>0;

%% load second level HMM and plot metastate profile:

W=125;K=3;
overlapstring='_overlappingWindows';
if whichstudy==4
  stochstring = '_stoch';
else
  stochstring = '';
end
load([config.hmmfolder,'secondLevelHMM',stochstring,'_Poiss_window',num2str(W),'_K',int2str(K),overlapstring,'.mat'],'hmmPoiss','feall','GammaPoiss','T_poiss','Poiss_subj_inds');

figdir = [config.figdir,'4_covariates_W',int2str(W),'_overlappingWindows/'];
mkdir(figdir);
%%
FO_2ndlevel = mean(GammaPoiss);
statedist_all = zeros(1,12);
for i=1:3
  statedist_all = statedist_all + FO_2ndlevel(i)*hmmPoiss.state(i).W.W_mean;
end
statedist_all = statedist_all ./ 125; % ratio rather than integer
for i=1:3
  colorweights(:,i) = hmmPoiss.state(i).W.W_mean ./ 125 - statedist_all;
end
colorweights = (colorweights-min(colorweights(:))+0.01);

colorweights = log(colorweights) - min(log(colorweights(:))) + 0.01;
colorweights = colorweights./(max(colorweights(:)+0.1));

% make plots:
optimalseqfile = [config.hmmfolder,'bestseq',replace(int2str(whichstudy), '5', '3'),'.mat'];
load(optimalseqfile);
bestseq = bestsequencemetrics{1};
%%
fig=setup_figure([], 1, 0.35);
CM = colormap(fig,hot(20));
for i=1:3
  ax(i) = axes('Position', [0.035+(i-1)*0.3 0.1 0.225, 0.8])
  for i2=1:12
    CW{i2} = CM(ceil(length(CM)*10.^colorweights(i2,i)/10),:);
  end
  cyclicalstatesubplot(bestseq,zeros(12),zeros(12),CW);
  set_font(8)
end
ax(4) = axes('Position', [0.035+(i)*0.275 0.15 0.1, .6])
h=colorbar;
h.Ticks = [0 1];
h.TickLabels = {'low', 'high'}
text(1.1, 1.2, 'FO')
axis off
 save_figure([config.figdir,'figure5_correlations/5_metastate_profile']);
%% check for correlation with cycle rates:
W = 125;



cycletime_mu = hmm_2ndlevel.cyctime_mu;
cycletime_med = hmm_2ndlevel.cyctime_med;
cycletime_std = hmm_2ndlevel.cyctime_std;

for i1=1:4 % run through different outlier / removal strategies
  
  outliers = abs(cycletime_mu(:,1)-mean(cycletime_mu(:,1))) > 2*std(cycletime_mu(:,1));
  outlier_string = '';
  if i1==2
    % option to instead remove outliers based on GMM model fit stats:
    outliers = supp.GMMoutliers;
    outlier_string = '_GMMoutliers';
  elseif i1==3
    outliers = supp.FOoutliers;
    outlier_string = '_HPSubjonly';
  elseif i1==4
    outliers = ~supp.FOoutliers;
    outlier_string = '_LPSubjonly';
  end
  
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
  
  %% check for metastate FO correlation:
  for ireg = [1:2];
    figure('Position',[6 435 1289 363]);
    for k =1:3
      subplot(1,3,k);
      scatter(info(~outliers,ireg),hmm_2ndlevel.FO(~outliers,k),'filled')
      [R,P] = corrcoef(info(~outliers,ireg),hmm_2ndlevel.FO(~outliers,k));
      title(['Metastate ',int2str(k),', rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
      plot4paper(reglabels{ireg},['Metastate ',int2str(k),' FO']);
      H=lsline(gca)
      set(H,'LineWidth',2);
      axis square;
    end
    print([figdir '4B_',outlier_string,'metastateFOCorr_',reglabels{ireg}],'-depsc')
  end
  
  %% side check: is this structure significant after regressing out all effect of lower level FO?
  
  B = pinv(hmm_1stlevel.FO)*hmm_2ndlevel.FO;
  FOmeta_resid = hmm_2ndlevel.FO - hmm_1stlevel.FO*B;
  for ireg = [1:2];
    figure('Position',[6 435 1289 363]);
    for k =1:3
      subplot(1,3,k);
      scatter(info(~outliers,ireg),FOmeta_resid(~outliers,k),'filled')
      [R,P] = corrcoef(info(~outliers,ireg),FOmeta_resid(~outliers,k));
      title(['Metastate ',int2str(k),', rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
      plot4paper(reglabels{ireg},['Metastate ',int2str(k),' FO']);
      axis square;
    end
    print([figdir '4B_',outlier_string,'metastateFOCorr_regressout_',reglabels{ireg}],'-depsc')
  end
  % note these remain significant as function of age even after lower level
  % FO is regressed out - suggests the hierarchical structure contains more
  % info than just that contained in the lower level HMM...
end

%% now check for patterns in 1st level metrics:
outliers = abs(cycletime_mu(:,1)-mean(cycletime_mu(:,1))) > 2*std(cycletime_mu(:,1));
outliers = outliers | isnan(info(:,3));
outlier_string = '';
% predict age from FOassym matrix:
x = permute(hmm_1stlevel.FO_assym,[3,1,2]);

x(isnan(x)) = 0;
x = normalise(x(:,find(1-eye(12))));
[a,b] = pca(x);
x = b(:,1:40);

[B,~,~,~,stats] = regress(info(~outliers,3),[ones(sum(~outliers),1),x(~outliers,:)]);
stats(3)
%%
%assym = hmm_1stlevel.FO_assym;
%X = permute(assym,[3,1,2]);
% compare to standard hmm metrics:
x = [normalise([hmm_1stlevel.FO,hmm_1stlevel.IT_mu,hmm_1stlevel.LT_mu])];
[B,~,~,~,stats] = regress(info(~outliers,3),[ones(sum(~outliers),1),x(~outliers,:)]);
stats(3)
%%
%%%%%%%%%%%%%%%%%%%%%% EXPLORATORY ANALYSES %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check quality of subject fits:

figure();

X = [hmm_1stlevel.FO,hmm_1stlevel.rotational_momentum,squeeze(nansum(nansum(abs(hmm_1stlevel.FO_assym)))),...
  hmm_2ndlevel.cyctime_mu,hmm_2ndlevel.cyctime_med,hmm_2ndlevel.cyctime_std,hmm_2ndlevel.FO];
imagesc(corr(X));

%labels = {'rotational_mom','absFOassym','firstlevelFOK1',

%% What is going on with fano factor outliers???

figure('Position',[1 457 1440 341]);
subplot(1,4,1);
x = hmm_1stlevel.rotational_momentum;
y = hmm_2ndlevel.cyctime_mu;
scatter(x,y)
lsline;
[R,P] = corrcoef(x,y);
title(['p=',num2str(P(2,1))]);
plot4paper('RotationalMomentum','cycletime');

% subplot(1,4,2);
% x = supp.fanofactor(:,2000);
% y = hmm_2ndlevel.cyctime_mu;
% scatter(x,y)
% lsline;
% [R,P] = corrcoef(x,y);
% title(['p=',num2str(P(2,1))]);
% plot4paper('FanoFactor','cycletime');
%title('IT correlation over subjects');
%axis square;

subplot(1,4,3);
x = squeeze(nansum(nansum(abs(hmm_1stlevel.FO_assym),2)));
y = hmm_2ndlevel.cyctime_mu;
scatter(x,y)
lsline;
[R,P] = corrcoef(x,y);
title(['p=',num2str(P(2,1))]);
plot4paper('AbsFOAssym','cycletime');

subplot(1,4,4);
x = squeeze(nansum(nansum(abs(hmm_1stlevel.FO_assym),2)));
y = hmm_1stlevel.rotational_momentum;
scatter(x,y)
lsline;
[R,P] = corrcoef(x,y);
title(['p=',num2str(P(2,1))]);
plot4paper('AbsFOAssym','RotationalMomentum');


print([config.figdir 'supp/','MetricCorrelation'],'-dpng')

% results from HCP suggest two groups - those that fit the global sequence,
% and those that really don't (have some other pattern)

%% Fit GMM with two components to distinguish groups:
rng(10);
x = squeeze(nansum(nansum(abs(hmm_1stlevel.FO_assym),2)));
y = hmm_1stlevel.rotational_momentum;
gmmfit = fitgmdist([x,y],2);
prob_fit = posterior(gmmfit,[x,y]);
i = sum(prob_fit);
if i(2)-i(1)>0
  prob_fit = fliplr(prob_fit);
end
cluster_assignment = prob_fit(:,1)>0.5;
figure();subplot(1,2,1);
scatter(x(cluster_assignment),y(cluster_assignment),'filled');
hold on;
scatter(x(~cluster_assignment),y(~cluster_assignment),'filled')
plot4paper('AbsFOAssym','RotationalMomentum');

GMMoutliers = ~cluster_assignment;
subplot(1,2,2);
stem(GMMoutliers);
title('GMM assigned bad subjects');
print([config.figdir 'supp/','PossibleExclusionCluster'])

% note this is NOT very consistent - if this turns out to be useful beyond
% HCP, use better ways of ensuring consistency

%% or try clustering participants by their FO (ie separate high and low variance subjects)

if exist(config.metricfile)
  load(config.metricfile,'supp');
  supp.GMMoutliers = GMMoutliers;
  supp.FOoutliers = highpowersubjs;
  save(config.metricfile,'supp','-append')
else
  error('Cant find config file')
end

%% are these determined by FO in state 1 and 12?
figure();
scatter(prob_fit(:,1),hmm_1stlevel.FO(:,1)+hmm_1stlevel.FO(:,12),'filled')
lsline;
% suggests the outliers are the ones that have more variance (ie swap
% between state 1 and 12 more)...does this generalise?

%%
x = squeeze(nansum(nansum(abs(hmm_1stlevel.FO_assym),2)));
y = hmm_1stlevel.rotational_momentum;
figure();
scatter(x,y,[],hmm_1stlevel.FO(:,1)+hmm_1stlevel.FO(:,12))

lsline;
[R,P] = corrcoef(x,y);
title(['p=',num2str(P(2,1))]);
plot4paper('AbsFOAssym','RotationalMomentum');

%%
%{
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
  fprintf(['\nPorcessing subj ',int2str(subnum)]);
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
  
end

%% do interval analysis on sub-groups:

[FO1,pvals1,t_intervals1] = computeLongTermAsymmetry(vpath(~supp.GMMoutliers),hmmT(~supp.GMMoutliers),K);
[FO2,pvals2,t_intervals2] = computeLongTermAsymmetry(vpath(supp.GMMoutliers),hmmT(supp.GMMoutliers),K);

% option to exclude state based outliers:
[FO2,pvals2,t_intervals2] = computeLongTermAsymmetry(vpath(supp.GMMoutliers & ~supp.FOoutliers),hmmT(supp.GMMoutliers & ~supp.FOoutliers),K);

bonf_ncomparisons = K.^2-K;

mean_direction1 = squeeze(mean(FO1(:,:,1,:)-FO1(:,:,2,:),4));
mean_direction2 = squeeze(mean(FO2(:,:,1,:)-FO2(:,:,2,:),4));

%%
% this script determines the optimal state ordering for a circular plot; it
% then determines whether such a sequentially organised network could arise
% by chance by random shuffles of the rows of the transmat
temp = squeeze(nansum(nansum(abs(squeeze((FO2(:,:,1,:)-FO2(:,:,2,:))./nanmean(FO2,3))),1),2));
to_exclude = temp>20;

bestsequencemetrics = optimiseSequentialPattern(FO1);
bestseq1 = bestsequencemetrics{1};
bestsequencemetrics = optimiseSequentialPattern(FO2(:,:,~to_exclude));
bestseq2 = bestsequencemetrics{1};
cyclicalstateplot(bestseq1,mean_direction1,pvals1<0.0000001*(0.05/bonf_ncomparisons));
cyclicalstateplot(bestseq1,mean_direction2,pvals2<(0.05/bonf_ncomparisons));

%%

assym1 = squeeze((FO1(:,:,1,:)-FO1(:,:,2,:))./nanmean(FO1,3));
assym2 = squeeze((FO2(:,:,1,:)-FO2(:,:,2,:))./nanmean(FO2,3));

disttoplot_manual2 = zeros(12,1);
disttoplot_manual1 = zeros(12,1);
for i3=1:12
  disttoplot_manual2(bestseq2(i3)) = exp(sqrt(-1)*i3/12*2*pi);
  disttoplot_manual1(bestseq1(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot1 = exp(sqrt(-1)*(angle(disttoplot_manual1.')-angle(disttoplot_manual1)));
angleplot2 = exp(sqrt(-1)*(angle(disttoplot_manual2.')-angle(disttoplot_manual2)));

rotational_momentum1 = imag(sum(sum(angleplot1.*assym1)));
rotational_momentum2 = imag(sum(sum(angleplot1.*assym2)));
figure();
subplot(1,2,1);
scatter(squeeze(nansum(nansum(abs(assym1),2),1)),-rotational_momentum1,'filled');

subplot(1,2,2);
scatter(squeeze(nansum(nansum(abs(assym2),2),1)),rotational_momentum2,'filled');


%% interim summary: on camcan, the second group (n=126) has some differences. these basically entirely disappear if FO outliers are also removed -
% appears these differences are driven by a few 'wierd' subjects

load('/Users/chiggins/Documents/CamCan/CogDatAll.mat')
%}
