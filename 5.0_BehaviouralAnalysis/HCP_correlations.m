% HCP correlations:

whichstudy = 3;
config = getStudyDetails(whichstudy);

% other preliminary setup for plotting etc:
color_scheme = set1_cols();
% include option for simulations here - results to be saved in adjacent
% folder:
simtests = false;
if simtests
    config.figdir = strrep(config.figdir,['Study',int2str(whichstudy)],['Study',int2str(whichstudy),'_simtest']);
    mkdir(config.figdir);
end
%% Load HMM summary stats:

load(config.metricfile);

%% Load HCP participant info:

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

InfMx = [];
InfMx.SubjectID = rstr_data.Subject;
InfMx.MotherID = rstr_data.Mother_ID;
InfMx.FatherID = rstr_data.Father_ID;
Zygosity = [];
for k=1:numel(inds)
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
save([config.resultsdir, 'APACE_InfMx'], 'InfMx')
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

%%

twindata = readmatrix([config.participantcovariates, 'twins.txt']);
twinmask = false(size(twindata,1),1);
for i=2:size(twindata)
    if ismember(twindata(1,i),subj_ids)
        twinmask(i) = true;
    end
end
twinstructure_MEG = twindata(twinmask,twinmask);
labels_twins = {'same','unrelated','monozygotic','dizygotic'};


%% correlation with age:

y = hmm_2ndlevel.cyctime_mu ./ config.sample_rate;
outliers = abs(y-mean(y))>2*std(y);
x = subjdata_detailed.('age');
figure('Position',[7 409 993 389]);
subplot(1,3,1);
scatter(x(~outliers),y(~outliers),'filled')
[R,P] = corrcoef(x(~outliers),y(~outliers));
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper('Age','Mean Cycle Time');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

y = hmm_2ndlevel.cyctime_med ./ config.sample_rate;
outliers = abs(y-mean(y))>2*std(y);

subplot(1,3,2);
scatter(x(~outliers),y(~outliers),'filled')
[R,P] = corrcoef(x(~outliers),y(~outliers));
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper('Age','Median Cycle Time');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

y = hmm_2ndlevel.cyctime_std ./ config.sample_rate;
outliers = abs(y-mean(y))>2*std(y);

subplot(1,3,3);
scatter(x(~outliers),y(~outliers),'filled')
[R,P] = corrcoef(x(~outliers),y(~outliers));
title(['rho=',num2str(R(2,1),2),', p=',num2str(P(2,1),2)]);
plot4paper('Age','Cycle Time std');
H=lsline(gca)
set(H,'LineWidth',2);
axis square;

%% twin structure:

cyctimediffmat = (1./hmm_2ndlevel.cyctime_mu - 1./hmm_2ndlevel.cyctime_mu').^2;
cyctimeL1diffmat = sqrt(cyctimediffmat);

offdiagselect = triu(true(79),1);

%anova1(cyctimeL1diffmat(offdiagselect),twinstructure_MEG(offdiagselect));
offdiagselect = triu(true(79));

cyctimediffmat = [];grouplabs = [];
for i1=1:2
    for i2=i1+1:3
        temp = sqrt((hmm_2ndlevel.cycletime_mu_sess(:,i1) - hmm_2ndlevel.cycletime_mu_sess(:,i2)').^2);
        cyctimediffmat = [cyctimediffmat;temp(offdiagselect)];
        grouplabs = [grouplabs;twinstructure_MEG(offdiagselect)];
    end
end
outliers = abs(cyctimediffmat-mean(cyctimediffmat))>2*std(cyctimediffmat);
anova1(cyctimediffmat(~outliers),grouplabs(~outliers))

%%

outliers = abs(cyctimediffmat-mean(cyctimediffmat))>2*std(cyctimediffmat);
selected = grouplabs==0 | grouplabs ==-1;
anova1(cyctimediffmat(~outliers & selected),grouplabs(~outliers& selected))

%% Look at corrrelation with Diego's fMRI metastates:

temp = load([config.fmri_metastates, 'BigHMM_820_uniqueP1_K12_mean.mat']);

temp3 = reshape(temp.BigGamma,4800,820,12);
temp4 = squeeze(mean(temp3,1));
FOcorr = corr(temp4);
%figure();
%subplot(1,2,1);
%imagesc(FOcorr)

[a,b] = pca(temp4,'NumComponents',2);
figure();scatter(b(:,1),b(:,2),'filled')

[A1,B1] = nets_hierarchy(FOcorr,FOcorr,1:12,'');
set(gcf,'Position',[337 350 477 388]);
% manually set ordering from Diego's paper:
new_order = [8,12,3,2,7,1,5,6,10,4,9,11];
print([config.figdir,'3AfMRIClustersFromDiegosPaper'],'-dpng');
%%
% b(:,1) is now each subject's strength of clustering
figure();scatter(b(:,1),b(:,2),'filled')

NoFamilyInfo = [108525, 116322, 146331, 168240, 256540, 280941, 524135, 604537, 657659, 867468];
vars = dlmread([config.fmri_metastates, 'scripts900/vars.txt'],' ');
%Drop = false(N,1);
%for n = NoFamilyInfo, Drop(vars(:,1)==n) = true; end
%vars = vars(~Drop,:);
%twins = twins(~Drop,~Drop);

subject_IDs_fmri = vars(:,1);

% load subject IDs for MEG:
infoMEG= HCP_getparticipantinfo();
subjectIDs_MEG = infoMEG(:,1);

for i=1:length(subjectIDs_MEG)
    has_fmri(i) = any(subject_IDs_fmri==subjectIDs_MEG(i));
    if has_fmri(i)
        subj_mapping(i) = find(subject_IDs_fmri==subjectIDs_MEG(i));
    end
end

%% start by plotting with ANOVA in clusters (b>0)

rng(2)
gmmfit = fitgmdist(b,2);
prob_fit = posterior(gmmfit,b);
prob_fit = prob_fit(:,2);
figure();
scatter(b(prob_fit>0.5,1),b(prob_fit>0.5,2),'filled');hold on;


scatter(b(prob_fit<=0.5,1),b(prob_fit<=0.5,2),'filled');hold on;

clustermember = b(:,1)>0;%prob_fit<=0.5;
print([config.figdir,'3BClustersFromFMRIFO'],'-dpng');

% or try Diego's clustering on other correlation mat:
FOcorr2 = corr(temp4');
%figure();
%subplot(1,2,1);
%imagesc(FOcorr)


[A1,B1] = nets_hierarchy(FOcorr2,FOcorr2,[1:820],'');

% manually extract cluster division:
clustermember = ones(820,1);
clustermember(A1(1:find(A1(1,:)==681))) = 2;


clustermember = prob_fit<=0.5; % above doesn't work!
%
grouplabels = {'Cognitive','Sensorimotor'};
[pval,anovatab] = anova1(1./hmm_2ndlevel.cyctime_mu,clustermember(subj_mapping)+1,'off');
figure('Position',[8 486 1433 312]);
subplot(1,5,1);
boxplot(1./hmm_2ndlevel.cyctime_mu,clustermember(subj_mapping)+1);
set(gca,'XTick',[1:2],'XTickLabel',grouplabels);
plot4paper('fMRI Group','MEG Cycle rate');
title(['p=',num2str(pval)]);

grouplabels = {'Cognitive','Sensorimotor'};
[pval,anovatab] = anova1(hmm_1stlevel.cycle_metrics.rotational_momentum,clustermember(subj_mapping)+1,'off');
subplot(1,5,2);
boxplot(hmm_1stlevel.cycle_metrics.rotational_momentum,clustermember(subj_mapping)+1);
set(gca,'XTick',[1:2],'XTickLabel',grouplabels);
plot4paper('fMRI Group','MEG Cycle strength');
title(['p=',num2str(pval)]);

for i=1:3
    [pval_FO2(i),anovatab_FO2{i}]=anova1(hmm_2ndlevel.FO(:,i),clustermember(subj_mapping)+1,'off');
    subplot(1,5,i+2);
    boxplot(hmm_2ndlevel.FO(:,i),clustermember(subj_mapping)+1);
    plot4paper('fMRI Group',['MEG Metastate ',int2str(i),' FO']);
    title(['p=',num2str(pval_FO2(i))]);
    set(gca,'XTick',[1:2],'XTickLabel',grouplabels);

end
print([config.figdir,'3C_2ndlevelcorrelations'],'-dpng');

%% compare to 1st level hmm FO:
figure('Position',[8 62 1433 736]);

for i=1:12
    [pval_FO1(i),anovatab_FO1{i}]=anova1(hmm_1stlevel.FO(:,i),clustermember(subj_mapping)+1,'off');
    subplot(3,4,i);
    boxplot(hmm_1stlevel.FO(:,i),clustermember(subj_mapping)+1);
    xlabel('fMRI Group'), ylabel(['MEG State ',int2str(i),' FO']);
    title(['p=',num2str(pval_FO1(i))]);
    set(gca,'XTick',[1:2],'XTickLabel',grouplabels);
end

print([config.figdir,'3D_1stlevelcorrelations'],'-dpng');
%%

figure();
subplot(2,2,1);bar(mean(temp4(clustermember==0,new_order)));
title('FMRI HMM states')
subplot(2,2,2);bar(mean(temp4(clustermember==1,new_order)));
title('FMRI HMM states')
subplot(2,2,3);bar(mean(hmm_1stlevel.FO(clustermember(subj_mapping)==0,:)));
title('MEG HMM states')
subplot(2,2,4);bar(mean(hmm_1stlevel.FO(clustermember(subj_mapping)==1,:)));
title('MEG HMM states')

%% and plot big FO matrix:

FOmat = [temp4(subj_mapping,new_order),hmm_1stlevel.FO];
figure();
imagesc(corr(FOmat));

%% actually cluster points rather than arbitrary zero threshold:

figure();
subplot(2,2,1);bar(mean(temp4(clustermember==0,:)));
subplot(2,2,2);bar(mean(temp4(clustermember==1,:)));
subplot(2,2,3);bar(mean(hmm_1stlevel.FO(clustermember(subj_mapping)==0,:)));
subplot(2,2,4);bar(mean(hmm_1stlevel.FO(clustermember(subj_mapping)==1,:)));

%% Can I predict sequence lengths from fMRI trans prob matrix?

for iSj=1:820
    [~,vpath] = max(temp3(:,iSj,:),[],3);
    for k1=1:12
        for k2=1:12
            Pmat(iSj,k1,k2) = sum(vpath(1:end-1)==k1 & vpath(2:end)==k2);
        end
    end
    LT(iSj) = nanmean(cellfun(@nanmean,getStateLifeTimes(vpath,1200*[1 1 1 1])));
    for k1=1:12
        Pmat(iSj,k1,:) = sum(Pmat(iSj,k1,:));
    end
end

%% 

[~,b_pc] = pca(Pmat(subj_mapping,:),'NumComponents',10);
[B,~,~,~,stats] =regress(log(hmm_2ndlevel.cyctime_mu),[ones(length(subj_mapping),1),b_pc]);
%standard_decoding(hmm_2ndlevel.cyctime_mu,Pmat(subj_mapping,:),

stats(3)

%% what do fMRI sequential patterns look like?
temp_VP = reshape(temp.PathVP(:,new_order),4*1200,820,12);
clear vpath hmmT
for iSj=1:820
    [~,vpath{iSj}] = max(temp_VP(:,iSj,:),[],3);
    hmmT{iSj} = [1200,1200,1200,1200];
end
K = 12;
[FO,pvals,t_intervals] = computeLongTermAsymmetry(vpath,hmmT,K);

bonf_ncomparisons = K.^2-K;

mean_direction = squeeze(nanmean(FO(:,:,1,:)-FO(:,:,2,:),4));
mean_assym = squeeze(nanmean((FO(:,:,1,:)-FO(:,:,2,:))./nanmean(FO,3),4));
% figure();subplot(121)
% imagesc(pvals<(0.05/(bonf_ncomparisons)))
% subplot(122)
% imagesc(pvals<(0.05))
hmm_fmri = [];
hmm_fmri.FO_intervals = FO;
hmm_fmri.FO_assym = squeeze((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3));

%%
% this script determines the optimal state ordering for a circular plot; it
% then determines whether such a sequentially organised network could arise
% by chance by random shuffles of the rows of the transmat

optimalseqfile = [config.hmmfolder,'bestseq',int2str(whichstudy),'_fmri.mat'];
if ~isfile(optimalseqfile)
    bestsequencemetrics = optimiseSequentialPattern(FO);
    save(optimalseqfile,'bestsequencemetrics');
else
    load(optimalseqfile);
end
bestseq = bestsequencemetrics{1};
% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
    disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));

rotational_momentum = imag(sum(sum(angleplot.*hmm_fmri.FO_assym)));
hmm_fmri.rotational_momentum = squeeze(rotational_momentum);


%% plot as circular diagram:

cyclicalstateplot(bestseq,mean_direction,pvals<0.0000001*(0.05/bonf_ncomparisons));

print([config.figdir,'1A_Cyclicalpattern_fmri'],'-dpng');
