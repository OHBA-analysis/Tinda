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

% also load in the best sequence and angleplot
optimalseqfile = [config.hmmfolder,'bestseq',int2str(3),'.mat'];
load(optimalseqfile);
bestseq = bestsequencemetrics{2};

% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
    disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));


%% Task 1: Working memory evoked response:
hmm_1stlevel=[];

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

hastaskdata = false(config.nSj,1);
embeddedlength = 14;

% go through subjects, loading evoked responses:
subjdata = readtable([config.participantcovariates, 'unrestricted_aquinn501_4_7_2017_9_4_13.csv']);
subjdata_detailed = readtable([config.participantcovariates, 'vars.txt']);
temp = readtable([config.participantcovariates, 'MEGfnames.csv']);
subj_ids = [];
for i=1:size(temp,1)
    subj_id = str2num(temp{i,1}{1}(7:12));
    if ~isempty(subj_id) && length(intersect(subj_id,subj_ids))==0
        subj_ids = [subj_ids;subj_id];
    end
end
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
subjdata_detailed=subjdata_detailed(hastaskdata==1,:);

t = 7/250:1/250:4;
t = t(1:946);
t = t-1.5; % 1.5 seconds is image onset
t_ofinterest = t>-0.5 & t<1;
idx_t0 = find(t_ofinterest);
idx_t0 = idx_t0(nearest(t,0));
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
        trialonset{i}=[];
        for i2=1:length(to_analyse)
            temp = load(mat_files_orth{to_analyse(i2)});
            if length(unique(temp.T))>1
                error('Epochs not uniformly sized')
            end
            gamtask = cat(2,gamtask,reshape(temp.Gamma_task,(temp.T(1)-embeddedlength),length(temp.T),hmm.K));
            trialonset{i} = cat(2, trialonset{i}, reshape(zeros(length(temp.Gamma_task),1),(temp.T(1)-embeddedlength),length(temp.T)));
        end
        trialonset{i}(idx_t0,:)=1;
        trialonset{i} = reshape(trialonset{i}, [], 1);
        
        for k=1:hmm.K
            gam_evoked_sj(:,k,i) = mean(gamtask(:,:,k) - FO_sj(i,k),2);
        end
        temp = reshape(gamtask(t_ofinterest,:,:),sum(t_ofinterest)*size(gamtask,2),hmm.K);
        [~,vpath_task{i}] = max(temp,[],2);
        [~,vpath_task_long{i}] = max(reshape(gamtask, [], hmm.K),[],2);
%         vpath_evoked{i} = reshape(vpath_task{i}, sum(t_ofinterest), size(gamtask,2));
        T_all_task{i} = sum(t_ofinterest) * ones(size(gamtask,2),1);
    end
end
trialonset=trialonset(hastaskdata);
vpath_task_long = vpath_task_long(hastaskdata);
vpath_task = vpath_task(hastaskdata);
T_all_task = T_all_task(hastaskdata);
gam_evoked_sj = gam_evoked_sj(:,:,hastaskdata);
 % this the epoch length of interest
gam_evoked_sj = gam_evoked_sj(t_ofinterest,:,:);
t = t(t_ofinterest);


opts = [];
opts.K = 12;
opts.Fs = config.sample_rate;
opts.dropstates=0;
LTmerged = nan(79, 12);
LTmedian = nan(79, 12);
LTstd = nan(79, 12);
FracOcc = nan(79, 12);
ITmerged = nan(79, 12);
ITmedian = nan(79, 12);
ITstd = nan(79, 12);
ix = find(hastaskdata);
for subnum=1:length(vpath_task)
  subnum
  LT = getStateLifeTimes(vpath_task{subnum},T_all_task{subnum},opts);
  LTmerged(ix(subnum),:) = cellfun(@mean,LT);
  LTmedian(ix(subnum),:) = cellfun(@median,LT);
  LTstd(ix(subnum),:) = cellfun(@std,LT);
  FracOcc(ix(subnum),:) = getFractionalOccupancy(vpath_task{subnum},sum(T_all_task{subnum}),opts);
  IT = getStateIntervalTimes(vpath_task{subnum},T_all_task{subnum},opts);
  ITmerged(ix(subnum),:) = cellfun(@mean,IT);
  ITmedian(ix(subnum),:) = cellfun(@median,IT);
  ITstd(ix(subnum),:) = cellfun(@std,IT);
end
  
hmm_1stlevel.WrkMem = [];
hmm_1stlevel.WrkMem.FO = FracOcc;
hmm_1stlevel.WrkMem.LT_mu = LTmerged;
hmm_1stlevel.WrkMem.LT_med = LTmedian;
hmm_1stlevel.WrkMem.LT_std = LTstd;
hmm_1stlevel.WrkMem.IT_mu = ITmerged;
hmm_1stlevel.WrkMem.IT_med = ITmedian;
hmm_1stlevel.WrkMem.IT_std = ITstd;

%% 

plotEvokedStateDist(gam_evoked_sj,t)
print([config.figdir,'FigA_WrkMemEvokedResponse'],'-dpng');

%% and check for cyclical patterns:
% as well as the original Tinda analysis, run it in GLM form after
% subtracting the evoked responde from the task vpath.
[FO_task,pvals_task,~,FO_task_residual, pvals_task_residual] = computeLongTermAsymmetry(vpath_task,T_all_task,hmm.K,[],1);
bonf_ncomparisons = hmm.K.^2-hmm.K;

mean_direction = squeeze(nanmean(FO_task(:,:,1,:)-FO_task(:,:,2,:),4));
mean_assym_raw = squeeze(mean((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3),4));

hmm_1stlevel.WrkMem.raw.FO_intervals = nan(12,12,2,79);
hmm_1stlevel.WrkMem.raw.FO_intervals(:,:,:,hastaskdata) = FO_task;
hmm_1stlevel.WrkMem.raw.FO_assym = nan(12,12,79);
hmm_1stlevel.WrkMem.raw.FO_assym(:,:,hastaskdata) = squeeze((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3));

rotational_momentum = imag(nansum(nansum(angleplot.*hmm_1stlevel.WrkMem.raw.FO_assym)));
hmm_1stlevel.WrkMem.raw.rotational_momentum = squeeze(rotational_momentum);

cyclicalstateplot(bestseq,mean_direction,pvals_task<(0.05/bonf_ncomparisons));
gcf()
print([config.figdir,'FigB2_WrkMemSequencePlot_raw'],'-dpng');

% get the metrics for the GLM residuals (i.e. accounting for evoked vpath)
mean_direction_residual = squeeze(nanmean(FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:),4));
mean_assym_raw = squeeze(mean((FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:))./mean(FO_task_residual,3),4));

hmm_1stlevel.WrkMem.FO_intervals = nan(12,12,2,79);
hmm_1stlevel.WrkMem.FO_intervals(:,:,:,hastaskdata) = FO_task_residual;
hmm_1stlevel.WrkMem.FO_assym = nan(12,12,79);
hmm_1stlevel.WrkMem.FO_assym(:,:,hastaskdata) = squeeze((FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:))./mean(FO_task_residual,3));

rotational_momentum_residual = imag(nansum(nansum(angleplot.*hmm_1stlevel.WrkMem.FO_assym)));
hmm_1stlevel.WrkMem.rotational_momentum = squeeze(rotational_momentum_residual);

cyclicalstateplot(bestseq,mean_direction_residual,pvals_task_residual<(0.05/bonf_ncomparisons));
gcf()
print([config.figdir,'FigB_WrkMemSequencePlot'],'-dpng');

% save metrics for later analysis:
savedir = [config.figdir, 'HMMsummarymetrics.mat'];
if exist(config.metricfile)
    save(savedir,'hmm_1stlevel','-append');
else
    save(savedir,'hmm_1stlevel')
end


% for i=1:length(vpath_evoked)
%     vpathcircle{i}= getCircleVpath(vpath_evoked{i}, bestseq);
%     
%     vpathcircle_task_long{i} = getCircleVpath(vpath_task_long{i}, bestseq);
%     
%     % do Fourier analysis on the circle vpath and epoch
%     [tmp, circlefreq_evoked, circletime_evoked] = fourieranalysis_circleVpath(vpathcircle_task_long{i}, trialonset{i});
%     circlepow_evoked(i,:,:) = squeeze(nanmean(abs(tmp).^2));
%     circlespctrm_evoked(i,:,:) = squeeze(nanmean(spctrm));
% end
% 
% figure; plot(circlefreq_evoked, nanmean(circlepow_evoked,3)), hold on, plot(circlefreq_evoked, nanmean(nanmean(circlepow_evoked,3)), 'k', 'LineWidth', 2) 
% title({'Circle vpath PSD per subject, WorkMem task', sprintf('trial length %g sec (%g Hz)', round(unique(diff(onsetidx))/250,4), round(1/(unique(diff(onsetidx))/250),2))});
% ylabel('PSD'), xlabel('Frequency (Hz)')
% print([config.figdir,'WrkMemCircleFreq'],'-dpng');

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
temp = readtable([config.participantcovariates, 'MEGfnames.csv']);
subj_ids = [];
for i=1:size(temp,1)
    subj_id = str2num(temp{i,1}{1}(7:12));
    if ~isempty(subj_id) && length(intersect(subj_id,subj_ids))==0
        subj_ids = [subj_ids;subj_id];
    end
end

hastaskdata = false(config.nSj,1);
embeddedlength = 14;

nsamples=1306;
gam_evoked_sj = zeros(nsamples,hmm.K,length(FO_sj));
t = 7/250:1/250:8;
t = t(1:1306);
t_ofinterest = t>0.75 & t<2.75;
idx_t0 = find(t_ofinterest);
idx_t0 = idx_t0(nearest(t,0));
for i=1:length(FO_sj)
    taskfiles = false(config.nSj,1);
    for i2=1:length(mat_files_orth)
        taskfiles(i2) = contains(mat_files_orth{i2},num2str(subj_ids(i))) && contains(mat_files_orth{i2},'_8-');
    end
    hastaskdata(i) = any(taskfiles);
    if hastaskdata(i)
        to_analyse = find(taskfiles);
        gamtask = [];
        trialonset{i}=[];
        for i2=1:length(to_analyse)
            temp = load(mat_files_orth{to_analyse(i2)});
            if length(unique(temp.T))>1
                error('Epochs not uniformly sized')
            end
            gamtask = cat(2,gamtask,reshape(temp.Gamma_task,(temp.T(1)-embeddedlength),length(temp.T),hmm.K));
            trialonset{i} = cat(2, trialonset{i}, reshape(zeros(length(temp.Gamma_task),1),(temp.T(1)-embeddedlength),length(temp.T)));
        end
        trialonset{i}(idx_t0,:)=1;
        trialonset{i} = reshape(trialonset{i}, [], 1);
        for k=1:hmm.K
            gam_evoked_sj(:,k,i) = mean(gamtask(:,:,k) - FO_sj(i,k),2);
        end
        temp = reshape(gamtask(t_ofinterest,:,:),sum(t_ofinterest)*size(gamtask,2),hmm.K);
        [~,vpath_task{i}] = max(temp,[],2);
        [~,vpath_task_long{i}] = max(reshape(gamtask, [], hmm.K),[],2);
        T_all_task{i} = sum(t_ofinterest) * ones(size(gamtask,2),1);
    end
end
trialonset = trialonset(hastaskdata);
gam_evoked_sj = gam_evoked_sj(:,:,hastaskdata);
 % this the epoch length of interest
gam_evoked_sj = gam_evoked_sj(t_ofinterest,:,:);
t = t(t_ofinterest);
vpath_task = vpath_task(hastaskdata);
vpath_task_long = vpath_task_long(hastaskdata);
T_all_task = T_all_task(hastaskdata);

opts = [];
opts.K = 12;
opts.Fs = config.sample_rate;
opts.dropstates=0;
LTmerged = nan(79, 12);
LTmedian = nan(79, 12);
LTstd = nan(79, 12);
FracOcc = nan(79, 12);
ITmerged = nan(79, 12);
ITmedian = nan(79, 12);
ITstd = nan(79, 12);
ix = find(hastaskdata);
for subnum=1:length(vpath_task)
  subnum
  LT = getStateLifeTimes(vpath_task{subnum},T_all_task{subnum},opts);
  LTmerged(ix(subnum),:) = cellfun(@mean,LT);
  LTmedian(ix(subnum),:) = cellfun(@median,LT);
  LTstd(ix(subnum),:) = cellfun(@std,LT);
  FracOcc(ix(subnum),:) = getFractionalOccupancy(vpath_task{subnum},sum(T_all_task{subnum}),opts);
  IT = getStateIntervalTimes(vpath_task{subnum},T_all_task{subnum},opts);
  ITmerged(ix(subnum),:) = cellfun(@mean,IT);
  ITmedian(ix(subnum),:) = cellfun(@median,IT);
  ITstd(ix(subnum),:) = cellfun(@std,IT);
end

hmm_1stlevel.StoryM1 = [];
hmm_1stlevel.StoryM1.FO = FracOcc;
hmm_1stlevel.StoryM1.LT_mu = LTmerged;
hmm_1stlevel.StoryM1.LT_med = LTmedian;
hmm_1stlevel.StoryM1.LT_std = LTstd;
hmm_1stlevel.StoryM1.IT_mu = ITmerged;
hmm_1stlevel.StoryM1.IT_med = ITmedian;
hmm_1stlevel.StoryM1.IT_std = ITstd;
%%

plotEvokedStateDist(gam_evoked_sj,t)
xlim([min(t),max(t)]);
print([config.figdir,'FigC_StoryMEvokedResponse1'],'-dpng');

%% and check sequential evoked patterns:
[FO_task,pvals_task,~,FO_task_residual, pvals_task_residual] = computeLongTermAsymmetry(vpath_task,T_all_task,hmm.K,[],1);
bonf_ncomparisons = hmm.K.^2-hmm.K;

mean_direction = squeeze(nanmean(FO_task(:,:,1,:)-FO_task(:,:,2,:),4));
mean_assym = squeeze(mean((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3),4));

hmm_1stlevel.StoryM1.raw.FO_intervals = nan(12,12,2,79);
hmm_1stlevel.StoryM1.raw.FO_intervals(:,:,:,hastaskdata) = FO_task;
hmm_1stlevel.StoryM1.raw.FO_assym = nan(12,12,79);
hmm_1stlevel.StoryM1.raw.FO_assym(:,:,hastaskdata) = squeeze((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3));

rotational_momentum = imag(nansum(nansum(angleplot.*hmm_1stlevel.StoryM1.raw.FO_assym)));
hmm_1stlevel.StoryM1.raw.rotational_momentum = squeeze(rotational_momentum);

cyclicalstateplot(bestseq,mean_direction,pvals_task<(0.05/bonf_ncomparisons));
gcf;
print([config.figdir,'FigD_StoryMSequencePlot1_raw'],'-dpng');

% correcting for evoked vpath:
mean_direction_residual = squeeze(nanmean(FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:),4));
mean_assym_raw = squeeze(mean((FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:))./mean(FO_task_residual,3),4));

hmm_1stlevel.StoryM1.FO_intervals = nan(12,12,2,79);
hmm_1stlevel.StoryM1.FO_intervals(:,:,:,hastaskdata) = FO_task_residual;
hmm_1stlevel.StoryM1.FO_assym = nan(12,12,79);
hmm_1stlevel.StoryM1.FO_assym(:,:,hastaskdata) = squeeze((FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:))./mean(FO_task_residual,3));

rotational_momentum_residual = imag(nansum(nansum(angleplot.*hmm_1stlevel.StoryM1.FO_assym)));
hmm_1stlevel.StoryM1.rotational_momentum = squeeze(rotational_momentum_residual);

cyclicalstateplot(bestseq,mean_direction_residual,pvals_task_residual<(0.05/bonf_ncomparisons));
gcf()
print([config.figdir,'FigD_StoryMSequencePlot1'],'-dpng');

% save metrics for later analysis:
savedir = [config.figdir, 'HMMsummarymetrics.mat']
if exist(config.metricfile)
    save(savedir,'hmm_1stlevel','-append');
else
    save(savedir,'hmm_1stlevel')
end

% % Circle Vpath analysis
% for i=1:length(vpath_evoked)
%     vpathcircle{i}= getCircleVpath(vpath_evoked{i}, bestseq);
%     
%     vpathcircle_task_long{i} = getCircleVpath(vpath_task_long{i}, bestseq);
%     
%     % do Fourier analysis on the circle vpath and epoch
%     [tmp, circlefreq_evoked, circletime_evoked] = fourieranalysis_circleVpath(vpathcircle_task_long{i}, trialonset{i});
%     circlepow_evoked(i,:,:) = squeeze(nanmean(abs(tmp).^2));
%     circlespctrm_evoked(i,:,:) = squeeze(nanmean(spctrm));
% end
% figure; plot(circlefreq_evoked, nanmean(circlepow_evoked,3)), hold on, plot(circlefreq_evoked, nanmean(nanmean(circlepow_evoked,3)), 'k', 'LineWidth', 2) 
% title({'Circle vpath PSD per subject, StoryMath task', sprintf('trial length %g sec (%g Hz)', round(nsamples/250,4), round(1/(nsamples/250),2))});
% ylabel('PSD'), xlabel('Frequency (Hz)')
% print([config.figdir,'StoryMCircleFreq'],'-dpng');

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
temp = readtable([config.participantcovariates, 'MEGfnames.csv']);
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
        gamtask_sj{i} = gamtask;
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
gamtask_sj = gamtask_sj(hastaskdata);
for k=1:length(gamtask_sj)
  gamtask_sj{k} = gamtask_sj{k}(t_ofinterest,:,:);
end
t = t(t_ofinterest);
vpath_task = vpath_task(hastaskdata);
T_all_task = T_all_task(hastaskdata);

opts = [];
opts.K = 12;
opts.Fs = config.sample_rate;
opts.dropstates=0;
LTmerged = nan(79, 12);
LTmedian = nan(79, 12);
LTstd = nan(79, 12);
FracOcc = nan(79, 12);
ITmerged = nan(79, 12);
ITmedian = nan(79, 12);
ITstd = nan(79, 12);
ix = find(hastaskdata);
for subnum=1:length(vpath_task)
  subnum
  LT = getStateLifeTimes(vpath_task{subnum},T_all_task{subnum},opts);
  LTmerged(ix(subnum),:) = cellfun(@mean,LT);
  LTmedian(ix(subnum),:) = cellfun(@median,LT);
  LTstd(ix(subnum),:) = cellfun(@std,LT);
  FracOcc(ix(subnum),:) = getFractionalOccupancy(vpath_task{subnum},sum(T_all_task{subnum}),opts);
  IT = getStateIntervalTimes(vpath_task{subnum},T_all_task{subnum},opts);
  ITmerged(ix(subnum),:) = cellfun(@mean,IT);
  ITmedian(ix(subnum),:) = cellfun(@median,IT);
  ITstd(ix(subnum),:) = cellfun(@std,IT);
end

hmm_1stlevel.StoryM2 = [];
hmm_1stlevel.StoryM2.FO = FracOcc;
hmm_1stlevel.StoryM2.LT_mu = LTmerged;
hmm_1stlevel.StoryM2.LT_med = LTmedian;
hmm_1stlevel.StoryM2.LT_std = LTstd;
hmm_1stlevel.StoryM2.IT_mu = ITmerged;
hmm_1stlevel.StoryM2.IT_med = ITmedian;
hmm_1stlevel.StoryM2.IT_std = ITstd;

%%
plotEvokedStateDist(gam_evoked_sj,t)
xlim([min(t),max(t)]);
gcf;
print([config.figdir,'FigE_StoryMEvokedResponse2'],'-dpng');

%% and check sequential evoked patterns:
[FO_task,pvals_task,~,FO_task_residual, pvals_task_residual] = computeLongTermAsymmetry(vpath_task,T_all_task,hmm.K,[],1);

bonf_ncomparisons = hmm.K.^2-hmm.K;

mean_direction = squeeze(nanmean(FO_task(:,:,1,:)-FO_task(:,:,2,:),4));
mean_assym = squeeze(mean((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3),4));

hmm_1stlevel.StoryM2.raw.FO_intervals = nan(12,12,2,79);
hmm_1stlevel.StoryM2.raw.FO_intervals(:,:,:,hastaskdata) = FO_task;
hmm_1stlevel.StoryM2.raw.FO_assym = nan(12,12,79);
hmm_1stlevel.StoryM2.raw.FO_assym(:,:,hastaskdata) = squeeze((FO_task(:,:,1,:)-FO_task(:,:,2,:))./mean(FO_task,3));

rotational_momentum = imag(nansum(nansum(angleplot.*hmm_1stlevel.StoryM2.raw.FO_assym)));
hmm_1stlevel.StoryM2.raw.rotational_momentum = squeeze(rotational_momentum);

cyclicalstateplot(bestseq,mean_direction,pvals_task<(0.05/bonf_ncomparisons));
gcf;
print([config.figdir,'FigF_StoryMSequencePlot2_raw'],'-dpng');

% correcting for evoked vpath:
mean_direction_residual = squeeze(nanmean(FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:),4));
mean_assym_raw = squeeze(mean((FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:))./mean(FO_task_residual,3),4));

hmm_1stlevel.StoryM2.FO_intervals = nan(12,12,2,79);
hmm_1stlevel.StoryM2.FO_intervals(:,:,:,hastaskdata) = FO_task_residual;
hmm_1stlevel.StoryM2.FO_assym = nan(12,12,79);
hmm_1stlevel.StoryM2.FO_assym(:,:,hastaskdata) = squeeze((FO_task_residual(:,:,1,:)-FO_task_residual(:,:,2,:))./mean(FO_task_residual,3));

rotational_momentum_residual = imag(nansum(nansum(angleplot.*hmm_1stlevel.StoryM2.FO_assym)));
hmm_1stlevel.StoryM2.rotational_momentum = squeeze(rotational_momentum_residual);

cyclicalstateplot(bestseq,mean_direction_residual,pvals_task_residual<(0.05/bonf_ncomparisons));
gcf()
print([config.figdir,'FigF_StoryMSequencePlot2'],'-dpng');

% save metrics for later analysis:
savedir = [config.figdir, 'HMMsummarymetrics.mat'];
if exist(config.metricfile)
    save(savedir,'hmm_1stlevel','-append');
else
    save(savedir,'hmm_1stlevel')
end
