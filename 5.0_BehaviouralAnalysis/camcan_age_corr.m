w%% load age data from camcan and check for relationship with FO assymetry:
if ~exist('whichstudy','var')
  whichstudy = 4; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);

% load in the TINDA results
load([config.resultsdir, 'HMMsummarymetrics'])

% set some parameters
bonf_ncomparisons = 132;


if whichstudy==4
  info = camcan_getparticipantinfo(config);
  
  % check how age predicts each FO assymetry measure:
  for k1=1:12
    for k2=1:12
      Y = permute(hmm_1stlevel.cycle_metrics.FO_assym(k1,k2,1,:),[4,1,2,3]);
      X = info(:,1);
      [R,P] = corrcoef(X,Y);
      rho(k1,k2) = R(1,2);
      pval(k1,k2) = P(1,2);
    end
  end
  mask = pval<(0.05/bonf_ncomparisons);
  cyclicalstateplot(bestseq,rho,pval<(0.05/bonf_ncomparisons));
  
  print([config.figdir,'1D_Cyclicalpattern_AgeEffect'],'-dpng');
  
  % does overall assymetry increase/ decrease with age?
  temp = (permute(hmm_1stlevel.FO_assym,[4,1,2,3]));
  metric = sum(sum(abs(temp(:,find(~eye(12)))),3),2);
  [R,P] = corrcoef(metric,X);
  figure('Position',[440 516 960 282]);
  subplot(1,2,1);
  scatter(X,metric,'filled');
  plot4paper('Age','Sequential assymetry')
  title(['rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
  % sequential patterns become slightly stronger with age
  % check if this due to outliers:
  outlier = abs(demean(metric))>2*std(metric);
  
  subplot(1,2,2);
  [R,P] = corrcoef(metric(~outlier),X(~outlier));
  scatter(X(~outlier),metric(~outlier),'filled');
  plot4paper('Age','Sequential assymetry')
  title(['Outliers removed: rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
  print([config.figdir,'1E_SequenceAssym_AgeEffect'],'-dpng');
  
  % or gender?
  p = anova1(metric(~outlier),info(~outlier,3));
  % significant effect for gender also
  
  % look at age within a gender group:
  figure('Position',[440 516 960 282]);
  subplot(1,2,1);
  outlier = abs(demean(metric))>2*std(metric) | info(:,3)==1;
  [R,P] = corrcoef(metric(~outlier),X(~outlier));
  scatter(X(~outlier),metric(~outlier),'filled');
  plot4paper('Age','Sequential assymetry')
  title(['Male group: rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
  subplot(1,2,2);
  outlier = abs(demean(metric))>2*std(metric) | info(:,3)==2;
  [R,P] = corrcoef(metric(~outlier),X(~outlier));
  scatter(X(~outlier),metric(~outlier),'filled');
  plot4paper('Age','Sequential assymetry')
  title(['Female group: rho=',num2str(R(2,1),2),',p=',num2str(P(2,1),1)])
end