% TaskSequenceAnalysis

% do the sequential patterns found for the HCP resting state data replicate
% to the inferred task structure?

clear all;
whichstudy = 5; % this flags the HCP task data
config = getStudyDetails(whichstudy);

% other preliminary setup for plotting etc:
color_scheme = set1_cols();

%% Extract baseline RSN activity from resting state scans:

load(fullfile(config.hmmfolder,config.hmmfilename),'hmm','T_all','Gamma');

% get each subject's baseline FO from resting state:
t0 = 0;
FO_sj = zeros(length(T_all),12);
for iSj=1:length(T_all)
    tsj= sum(T_all{iSj}-14);
    FO_sj(iSj,:) = mean(Gamma((t0+1):(t0+tsj),:));
    t0 = t0 + tsj;
end
FO_sj = squeeze(mean(reshape(FO_sj,[3,config.nSj,12])));

%% Task 1: Working memory evoked response:
% extract file structure:
if ~isfile(config.wrkmemfilelist)
    fnames = dir(config.wrkmemdir);
    mat_files_orth = {};
    validfiles = false(length(fnames),1);
    for i=1:length(fnames)
        mat_files_orth{i} = [fnames(i).folder,'/',fnames(i).name];
        validfiles(i) = contains(mat_files_orth{i},'TIM_orth.mat');
    end
    mat_files_orth = mat_files_orth(validfiles);
    save(config.wrkmemfilelist,'mat_files_orth');
else
    load(config.wrkmemfilelist,'mat_files_orth');
end

% go through subjects, loading evoked responses:
% subjdata = readtable('/Users/chiggins/data/HCPAnalysis/behav/unrestricted_aquinn501_4_7_2017_9_4_13.csv');
% subjdata_detailed = readtable('/Users/chiggins/data/HCPAnalysis/behav/vars.txt');
% headers = readtable('/Users/chiggins/data/HCPAnalysis/behav/column_headers.txt');
temp = readtable('/ohba/pi/mwoolrich/datasets/HCP_CH_2022/HCPAnalysis/behav');
subj_ids = [];
for i=1:size(temp,1)
    subj_id = str2num(temp{i,1}{1}(7:12));
    if ~isempty(subj_id) && length(intersect(subj_id,subj_ids))==0
        subj_ids = [subj_ids;subj_id];
    end
end

hastaskdata = false(config.nSj,1);
embeddedlength = 14;

t = 7/250:1/250:4;
t = t(1:946);
t = t-1.5; % 1.5 seconds is image onset
t_ofinterest = t>-0.5 & t<1;
gam_evoked_sj = zeros(946,hmm.K,length(FO_sj));
for i=1:length(FO_sj)
    taskfiles = false(config.nSj,1);
    for i2=1:length(mat_files_orth)
        taskfiles(i2) = contains(mat_files_orth{i2},num2str(subj_ids(i)));
    end
    hastaskdata(i) = any(taskfiles);
    if hastaskdata(i)
        to_analyse = find(taskfiles);
        gamtask = [];
        for i2=1:length(to_analyse)
            temp = load(mat_files_orth{to_analyse(i2)});
            if length(unique(temp.T))>1
                error('Epochs not uniformly sized')
            end
            gamtask = cat(2,gamtask,reshape(temp.Gamma_task,(temp.T(1)-embeddedlength),length(temp.T),hmm.K));
        end
        for k=1:hmm.K
            gam_evoked_sj(:,k,i) = mean(gamtask(:,:,k) - FO_sj(i,k),2);
        end
        temp = reshape(gamtask(t_ofinterest,:,:),sum(t_ofinterest)*size(gamtask,2),hmm.K);
        [~,vpath_task{i}] = max(temp,[],2);
        T_all_task{i} = sum(t_ofinterest) * ones(size(gamtask,2),1);
    end
end
vpath_task = vpath_task(hastaskdata);
T_all_task = T_all_task(hastaskdata);
gam_evoked_sj = gam_evoked_sj(:,:,hastaskdata);
 % this the epoch length of interest
gam_evoked_sj = gam_evoked_sj(t_ofinterest,:,:);
t = t(t_ofinterest);

%% 

plotEvokedStateDist(gam_evoked_sj,t)
print([config.figdir,'FigA_WrkMemEvokedResponse'],'-dpng');

%% and check for cyclical patterns:
[FO_task,pvals_task] = computeLongTermAsymmetry(vpath_task,T_all_task,hmm.K);
bonf_ncomparisons = hmm.K.^2-hmm.K;

mean_direction = squeeze(nanmean(FO_task(:,:,1,:)-FO_task(:,:,2,:),4));
mean_assym = squeeze(mean((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3),4));

optimalseqfile = ['/ohba/pi/mwoolrich/datasets/HCP_CH_2022/ve_output_rest/','bestseq',int2str(3),'.mat'];
load(optimalseqfile);
bestseq = bestsequencemetrics{2};
cyclicalstateplot(bestseq,mean_direction,pvals_task<(0.05/bonf_ncomparisons));
gcf()
print([config.figdir,'FigB_WrkMemSequencePlot'],'-dpng');


%% Task 2: StoryMath task
% extract file structure:
if 1;~isfile(config.storymfilelist)
    fnames = dir(config.storymdir);
    mat_files_orth = {};
    validfiles = false(length(fnames),1);
    for i=1:length(fnames)
        mat_files_orth{i} = [fnames(i).folder,'/',fnames(i).name];
        validfiles(i) = contains(mat_files_orth{i},'TEV_orth.mat');
    end
    mat_files_orth = mat_files_orth(validfiles);
    save(config.storymfilelist,'mat_files_orth');
else
    load(config.storymfilelist,'mat_files_orth');
end

% go through subjects, loading evoked responses:
% subjdata = readtable('/Users/chiggins/data/HCPAnalysis/behav/unrestricted_aquinn501_4_7_2017_9_4_13.csv');
% subjdata_detailed = readtable('/Users/chiggins/data/HCPAnalysis/behav/vars.txt');
% headers = readtable('/Users/chiggins/data/HCPAnalysis/behav/column_headers.txt');
temp = readtable('/ohba/pi/mwoolrich/datasets/HCP_CH_2022/HCPAnalysis/behav/MEGfnames.csv');
subj_ids = [];
for i=1:size(temp,1)
    subj_id = str2num(temp{i,1}{1}(7:12));
    if ~isempty(subj_id) && length(intersect(subj_id,subj_ids))==0
        subj_ids = [subj_ids;subj_id];
    end
end

hastaskdata = false(config.nSj,1);
embeddedlength = 14;


gam_evoked_sj = zeros(1306,hmm.K,length(FO_sj));
t = 7/250:1/250:8;
t = t(1:1306);
%t = t-1.5; % 1.5 seconds is image onset
t_ofinterest = t>0.75 & t<2.75;
for i=1:length(FO_sj)
    taskfiles = false(config.nSj,1);
    for i2=1:length(mat_files_orth)
        taskfiles(i2) = contains(mat_files_orth{i2},num2str(subj_ids(i))) && contains(mat_files_orth{i2},'_8-');
    end
    hastaskdata(i) = any(taskfiles);
    if hastaskdata(i)
        to_analyse = find(taskfiles);
        gamtask = [];
        for i2=1:length(to_analyse)
            temp = load(mat_files_orth{to_analyse(i2)});
            if length(unique(temp.T))>1
                error('Epochs not uniformly sized')
            end
            gamtask = cat(2,gamtask,reshape(temp.Gamma_task,(temp.T(1)-embeddedlength),length(temp.T),hmm.K));
        end
        for k=1:hmm.K
            gam_evoked_sj(:,k,i) = mean(gamtask(:,:,k) - FO_sj(i,k),2);
        end
        temp = reshape(gamtask(t_ofinterest,:,:),sum(t_ofinterest)*size(gamtask,2),hmm.K);
        [~,vpath_task{i}] = max(temp,[],2);
        T_all_task{i} = sum(t_ofinterest) * ones(size(gamtask,2),1);
    end
end
gam_evoked_sj = gam_evoked_sj(:,:,hastaskdata);
 % this the epoch length of interest
gam_evoked_sj = gam_evoked_sj(t_ofinterest,:,:);
t = t(t_ofinterest);
vpath_task = vpath_task(hastaskdata);
T_all_task = T_all_task(hastaskdata);
plotEvokedStateDist(gam_evoked_sj,t)
xlim([min(t),max(t)]);
print([config.figdir,'FigC_StoryMEvokedResponse1'],'-dpng');

%% and check sequential evoked patterns:
[FO_task,pvals_task] = computeLongTermAsymmetry(vpath_task,T_all_task,hmm.K);
bonf_ncomparisons = hmm.K.^2-hmm.K;

mean_direction = squeeze(nanmean(FO_task(:,:,1,:)-FO_task(:,:,2,:),4));
mean_assym = squeeze(mean((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3),4));


optimalseqfile = ['/ohba/pi/mwoolrich/datasets/HCP_CH_2022/ve_output_rest/','bestseq',int2str(3),'.mat'];
load(optimalseqfile);
bestseq = bestsequencemetrics{2};
cyclicalstateplot(bestseq,mean_direction,pvals_task<(0.05/bonf_ncomparisons));
gcf;
print([config.figdir,'FigD_WrkMemSequencePlot'],'-dpng');

%% Task 3: StoryMath task 2
% extract file structure:
if 1;~isfile(config.storymfilelist)
    fnames = dir(config.storymdir);
    mat_files_orth = {};
    validfiles = false(length(fnames),1);
    for i=1:length(fnames)
        mat_files_orth{i} = [fnames(i).folder,'/',fnames(i).name];
        validfiles(i) = contains(mat_files_orth{i},'TEV_orth.mat');
    end
    mat_files_orth = mat_files_orth(validfiles);
    save(config.storymfilelist,'mat_files_orth');
else
    load(config.storymfilelist,'mat_files_orth');
end

% go through subjects, loading evoked responses:
% subjdata = readtable('/Users/chiggins/data/HCPAnalysis/behav/unrestricted_aquinn501_4_7_2017_9_4_13.csv');
% subjdata_detailed = readtable('/Users/chiggins/data/HCPAnalysis/behav/vars.txt');
% headers = readtable('/Users/chiggins/data/HCPAnalysis/behav/column_headers.txt');
temp = readtable('/ohba/pi/mwoolrich/datasets/HCP_CH_2022/HCPAnalysis/behav/MEGfnames.csv');
subj_ids = [];
for i=1:size(temp,1)
    subj_id = str2num(temp{i,1}{1}(7:12));
    if ~isempty(subj_id) && length(intersect(subj_id,subj_ids))==0
        subj_ids = [subj_ids;subj_id];
    end
end

hastaskdata = false(config.nSj,1);
embeddedlength = 14;


gam_evoked_sj = zeros(706,hmm.K,length(FO_sj));
t = 7/250:1/250:8;
t = t(1:706);
%t = t-1.5; % 1.5 seconds is image onset
t_ofinterest = t>0.5 & t<2.5;
for i=1:length(FO_sj)
    taskfiles = false(config.nSj,1);
    for i2=1:length(mat_files_orth)
        taskfiles(i2) = contains(mat_files_orth{i2},num2str(subj_ids(i))) && contains(mat_files_orth{i2},'_9-');
    end
    hastaskdata(i) = any(taskfiles);
    if hastaskdata(i)
        to_analyse = find(taskfiles);
        gamtask = [];
        for i2=1:length(to_analyse)
            temp = load(mat_files_orth{to_analyse(i2)});
            if length(unique(temp.T))>1
                error('Epochs not uniformly sized')
            end
            gamtask = cat(2,gamtask,reshape(temp.Gamma_task,(temp.T(1)-embeddedlength),length(temp.T),hmm.K));
        end
        for k=1:hmm.K
            gam_evoked_sj(:,k,i) = mean(gamtask(:,:,k) - FO_sj(i,k),2);
        end
        temp = reshape(gamtask(t_ofinterest,:,:),sum(t_ofinterest)*size(gamtask,2),hmm.K);
        [~,vpath_task{i}] = max(temp,[],2);
        T_all_task{i} = sum(t_ofinterest) * ones(size(gamtask,2),1);
    end
end
gam_evoked_sj = gam_evoked_sj(:,:,hastaskdata);
 % this the epoch length of interest
gam_evoked_sj = gam_evoked_sj(t_ofinterest,:,:);
t = t(t_ofinterest);
vpath_task = vpath_task(hastaskdata);
T_all_task = T_all_task(hastaskdata);
plotEvokedStateDist(gam_evoked_sj,t)
xlim([min(t),max(t)]);
gcf;
print([config.figdir,'FigE_StoryMEvokedResponse2'],'-dpng');

%% and check sequential evoked patterns:
[FO_task,pvals_task] = computeLongTermAsymmetry(vpath_task,T_all_task,hmm.K);
bonf_ncomparisons = hmm.K.^2-hmm.K;

mean_direction = squeeze(nanmean(FO_task(:,:,1,:)-FO_task(:,:,2,:),4));
mean_assym = squeeze(mean((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3),4));


optimalseqfile = ['/ohba/pi/mwoolrich/datasets/HCP_CH_2022/ve_output_rest/','bestseq',int2str(3),'.mat'];
load(optimalseqfile);
bestseq = bestsequencemetrics{2};
cyclicalstateplot(bestseq,mean_direction,pvals_task<(0.05/bonf_ncomparisons));
gcf;
print([config.figdir,'FigF_WrkMemSequencePlot'],'-dpng');

%% 