% script called to load viterbi paths inferred and hmm objects and run
% post-hoc sequence analysis:
if ~exist('whichstudy','var')
    whichstudy = 4; % 1 denotes the hmm model run in Higgins2020_neuron
end
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
    if ~isfolder(mat_files_orth{1})
        for k=1:length(mat_files_orth)
            tmp1 = fileparts(config.matfilelist);
            [~, tmp2] = fileparts(mat_files_orth{k});
            mat_files_orth{k} = fullfile(tmp1, tmp2);
        end
    end
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
    fprintf(['\nProcessing subj ',int2str(subnum)]);
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
    LTmedian(subnum,:) = cellfun(@median,LT);
    LTstd(subnum,:) = cellfun(@std,LT);
    FracOcc(subnum,:) = getFractionalOccupancy(vpath{subnum},sum(hmmT{subnum}),opts);
    IT = getStateIntervalTimes(vpath{subnum},hmmT{subnum},opts);
    ITmerged(subnum,:) = cellfun(@mean,IT);
    ITmedian(subnum,:) = cellfun(@median,IT);
    ITstd(subnum,:) = cellfun(@std,IT);
end

% save key metrics to look at covariates later:
hmm_1stlevel = [];
hmm_1stlevel.FO = FracOcc;
hmm_1stlevel.LT_mu = LTmerged;
hmm_1stlevel.LT_med = LTmedian;
hmm_1stlevel.LT_std = LTstd;
hmm_1stlevel.IT_mu = ITmerged;
hmm_1stlevel.IT_med = ITmedian;
hmm_1stlevel.IT_std = ITstd;

%% summary plots:
fontsize = 18;
figure;subplot(111);
distributionPlot(FracOcc,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}});
set(gca,'YLim',[0 1.1*max(FracOcc(:))],'FontSize',fontsize)
title('Fractional Occupancy');plot4paper('RSN-State','Proportion');grid on;
print([config.figdir '0_temporalstats_FO'],'-depsc')

figure;subplot(111);
distributionPlot(LTmerged ./ config.sample_rate * 1000,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Life Times');plot4paper('RSN-State','Time (ms)');grid on;
YL = 1.1*max(LTmerged(:))./ config.sample_rate * 1000;
set(gca,'YLim',[0 YL],'FontSize',fontsize,'FontSize',fontsize);
print([config.figdir '0_temporalstats_LT'],'-depsc')

figure;subplot(111);
distributionPlot(log10(ITmerged ./ config.sample_rate),'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Interval Times');plot4paper('RSN-State','Time (secs)');grid on
YL(2) =1.5* max(mean(log10(ITmerged ./ config.sample_rate)));
YL(1) = min(squash(log10(ITmerged ./ config.sample_rate)));
set(gca,'YLim',YL,'FontSize',fontsize)
set(gca,'YTick',log10([0.05,0.1,0.5,1,5,10]))
y_labels = get(gca,'YTickLabel');
for i=1:length(y_labels)
    y_labels{i}=num2str(10.^(str2num(y_labels{i})),1);
end
set(gca,'YTickLabels',y_labels);
print([config.figdir '0_temporalstats_IT_logscale'],'-depsc')

figure;subplot(111);
distributionPlot(ITmerged ./ config.sample_rate,'showMM',2,'color',{color_scheme{1:size(FracOcc,2)}})
title('Interval Times');plot4paper('RSN-State','Time (secs)');grid on
YL(2) =1.5* max(mean((ITmerged ./ config.sample_rate)));
YL(1) = 0;
set(gca,'YLim',YL,'FontSize',fontsize)
print([config.figdir '0_temporalstats_IT'],'-depsc')

%% and actually compute long term assymetry:

[FO,pvals,t_intervals] = computeLongTermAsymmetry(vpath,hmmT,K);

bonf_ncomparisons = K.^2-K;

mean_direction = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
mean_assym = squeeze(mean((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3),4));
% figure();subplot(121)
% imagesc(pvals<(0.05/(bonf_ncomparisons)))
% subplot(122)
% imagesc(pvals<(0.05))

hmm_1stlevel.FO_intervals = FO;
hmm_1stlevel.FO_assym = squeeze((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3));

%%
% this script determines the optimal state ordering for a circular plot; it
% then determines whether such a sequentially organised network could arise
% by chance by random shuffles of the rows of the transmat

optimalseqfile = [config.hmmfolder,'bestseq',int2str(whichstudy),'.mat'];
if ~isfile(optimalseqfile)
    bestsequencemetrics = optimiseSequentialPattern(FO);
    save(optimalseqfile,'bestsequencemetrics');
else
    load(optimalseqfile);
end
bestseq = bestsequencemetrics{2};
% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
    disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));

rotational_momentum = imag(sum(sum(angleplot.*hmm_1stlevel.FO_assym)));
hmm_1stlevel.rotational_momentum = squeeze(rotational_momentum);





%% plot as circular diagram:

if whichstudy<4
    cyclicalstateplot(bestseq,mean_direction,pvals<(0.05/bonf_ncomparisons));
else
    cyclicalstateplot(bestseq,mean_direction,pvals<0.0000001*(0.05/bonf_ncomparisons));
end
print([config.figdir,'1A_Cyclicalpattern'],'-dpng');

%% or plot as multiple individual state plots:
if whichstudy<4
    cyclicalstateplot_perstate(bestseq,mean_direction,pvals<(0.05/bonf_ncomparisons));
else
    cyclicalstateplot_perstate(bestseq,mean_direction,pvals<0.0000001*(0.05/bonf_ncomparisons));
end
print([config.figdir,'1BStatePathways'],'-dpng');

% also print for legend:
figure('Position', [440 579 114 219]);
quiver(0,0,1,0,'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.8);hold on;
quiver(0,1,1,0,'Color',[0 0 0]+0.8,'LineWidth',2,'MaxHeadSize',0.8);hold on;
axis off;        
print([config.figdir,'1BStatePathways_legend'],'-dpng');

%% analyse by quartiles
percentiles = 0:20:100;
figure('Position',[7 189 1434 609]);
for ip=1:5
    [FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));
    mean_direction_p = squeeze(mean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
    subplot(2,5,ip);
    cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(0.05/bonf_ncomparisons));
    subplot(2,5,5+ip);
    ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
    distributionPlot(ITmerged,'showMM',2,'color',{color_scheme{1:size(FO,2)}});
    set(gca,'YLim',[0 1.1*max(ITmerged(:))],'FontSize',12)
    title('Interval times');plot4paper('RSN-State','Interval time (sec)');grid on;
end
print([config.figdir,'1CCyclicalpatterns_percentiled'],'-dpng');
%hmm_1stlevel.FO_quartile = FO_p;
%% analyse by intervals >2 heartbeats long

% assume HR are greater than 50BPM:
maxHR = (60/50*2);
percentile = [maxHR*config.sample_rate,NaN]; % shortest permissible intervals:
figure('Position',[619 146 381 652]);
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentile);
mean_direction_p = squeeze(nanmean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
subplot(2,1,1);
cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(0.05/bonf_ncomparisons));
subplot(2,1,2);
ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
distributionPlot(ITmerged,'showMM',2,'color',{color_scheme{1:size(FO,2)}});
set(gca,'YLim',[0 1.1*max(ITmerged(:))],'FontSize',12)
title('Interval times');plot4paper('RSN-State','Interval time (sec)');grid on;
print([config.figdir,'1DCyclicalpatterns_greaterthan2sec'],'-dpng');

% save metrics:
hmm_1stlevel.FO_cardiaccontrol = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO,3));

%% analyse intervals <half a respiratory cycle 

% assume resp cycle of 16BPM;
minRR = 60/16*0.5;
percentile = [NaN,minRR*config.sample_rate]; % shortest permissible intervals:
figure('Position',[619 146 381 652]);
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentile);
mean_direction_p = squeeze(nanmean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
subplot(2,1,1);
cyclicalstatesubplot(bestseq,mean_direction_p,pvals_p<(0.05/bonf_ncomparisons));
subplot(2,1,2);
ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
distributionPlot(ITmerged,'showMM',2,'color',{color_scheme{1:size(FO,2)}});
set(gca,'YLim',[0 1.1*max(ITmerged(:))],'FontSize',12)
title('Interval times');plot4paper('RSN-State','Interval time (sec)');grid on;
print([config.figdir,'1ECyclicalpatterns_lessthan2.5sec'],'-dpng');

hmm_1stlevel.FO_respcontrol = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO,3));


% save metrics for later analysis:
if exist(config.metricfile)
    save(config.metricfile,'hmm_1stlevel','-append');
else
    save(config.metricfile,'hmm_1stlevel')
end

%% load wavelet PSDs for each state:
if whichstudy==3
    % for HCP need to recompute run indices (each subject has multiple runs)
    run_inds = zeros(size(hmm.statepath));
    t_offset = 0;
    for i=1:length(hmmT)
        t_length = sum(hmmT{i}) - length(hmmT{1})*(length(hmm.train.embeddedlags)-1);
        run_inds(t_offset + [1:t_length]) = i;
        t_offset = t_offset + t_length;
    end
    [P,f] = loadHMMspectra(config,whichstudy,hmm,run_inds);
elseif whichstudy==4
    % compute FO per subj:
    for i=1:length(vpath)
        for k=1:K
            FO_subj(i,k) = mean(vpath{i}==k);
        end
    end
    [P,f] = loadHMMspectra(config,whichstudy,hmm,[],FO_subj);
else
    [P,f] = loadHMMspectra(config,whichstudy,hmm,hmm.subj_inds);
    
end


% %% temp sanity check: load random subject spectra:
% i = randi(234);
% temp = load(['/Volumes/CamsHD2/HCP_CH/ve_output_rest/WTspect/hmm_spectrawt_sj',int2str(i),'.mat']);
% for i=1:12
%     P(i,:,:,:) = temp.fitmt_subj.state(i).psd;
% end
%% plot average state spectra vs freq (averaged over all parcels):
diagselect = find(eye(config.parc.n_parcels));
Pmean = mean(abs(P(:,:,diagselect)),3);
figure();
ls = {'-','--'};
for i=1:hmm.K
    plot(f,Pmean(i,:),'Color',color_scheme{i},'LineStyle',ls{1+(i>6)},'LineWidth',2);
    hold on;
    leglabels{i} = ['State ',int2str(i)];
end
plot4paper('Freq (Hz)','PSD');
legend(leglabels)
print([config.figdir,'2APSDperstate'],'-dpng');

%% fit NNMF over space and find modes:
P2 = reshape(abs(P(:,:,diagselect)),[12*size(Pmean,2),config.parc.n_parcels]);
%f_select = 1:size(psd_m,1);
%freq_bands = linspace(1,45,length(f_select));
maxP = 6;
nnmffile = [config.hmmfolder,'nnmf_spatialfit.mat'];
if isfile(nnmffile)
    load(nnmffile)
else
    for iter=1:10
        fprintf(['Run ',int2str(iter),'\n'])  
        [atemp btemp]=nnmf_mww(P2,maxP,'replicates',500,'algorithm','als');
        res = mean(mean((P2 - (atemp*btemp)).^2));
        if iter==1 || res<lastres
            a=atemp;
            b=btemp;
            lastres = res;
        end
    end
    clear atemp btemp
    a = reshape(a,[12,size(P2,1)/12,maxP]);
    save(nnmffile,'a','b','maxP');
end

%%
% find lowest power state and position at top of wheel:
[~,lowpower] = min(sum(Pmean,2));
ind = find(bestseq==lowpower);
bestseq_LP = [bestseq(ind:end),bestseq(1:(ind-1))];

% or manually reorder:
%bestseq_LP = [12     9    11     8     6     4     3     2     1     5     7    10];
if whichstudy==1
    bestorder = [5,4,2,6,3,1];
elseif whichstudy==3
    bestorder = [6,2,3,5,1,4];
elseif whichstudy==4
    bestorder = [1:6];
end
figure();
plotCyclicalTimeFreqPattern(bestseq_LP,config.parc,a(:,:,bestorder),b(bestorder,:),f,config)
print([config.figdir,'2BPSD_modes'],'-dpng');

%% refit earlier spatial modes to camcan data:

if whichstudy==4
    MEGUKfolder = '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/CanonicalRS/250Hz/hmm_1to45hz/'
    nnmffile_MEGUK = [MEGUKfolder,'nnmf_spatialfit.mat'];
    load(nnmffile_MEGUK)
    
    opt = statset('maxiter',1);    
    [b_new,a_new]=nnmf_mww(P2',maxP,'replicates',500,'algorithm','als','w0',b','options',opt);
    b_new = b_new';
    a_new = a_new';
    a_new = reshape(a_new,[12,size(P2,1)/12,maxP]);
    
    figure();
    bestorder = [1:6];
    plotCyclicalTimeFreqPattern(bestseq_LP,config.parc,a_new(:,:,bestorder),b_new(bestorder,:),f, config)
    print([config.figdir,'2CPSD_modes_samefit'],'-dpng');

end

%% load age data from camcan and check for relationship with FO assymetry:
if whichstudy==4
    info = camcan_getparticipantinfo(config);
    
    % check how age predicts each FO assymetry measure:
    assym = (FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3);
    for k1=1:12
        for k2=1:12
            Y = permute(assym(k1,k2,1,:),[4,1,2,3]);
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
    temp = (permute(assym,[4,1,2,3]));
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


%% Circular Vpath analysis

W=125;
W_half = (W-1)/2;
for i=1:length(vpath)
    vpathcircle{i}= getCircleVpath(vpath{i}, bestseq);
    tmp = zeros((W-1)/2,1);
    vpathcircle_window{i} = tmp;
    for k=W_half+1:length(vpathcircle{i})-W_half
        vpathcircle_window{i} = [vpathcircle_window{i}; mean(vpathcircle{i}(k-W_half:k+W_half))];
    end
    vpathcircle_window{i} = [vpathcircle_window{i}; tmp];
    
   
    % do Fourier analysis on the circle vpath
    [tmp, circlefreq, circletime] = fourieranalysis_circleVpath(vpathcircle{i}, []);
    circlepow(i,:,:) = squeeze(nanmean(abs(tmp).^2));
    circlespctrm(i,:,:) = squeeze(nanmean(tmp));
    
    [tmp, circlefreq, circletime] = fourieranalysis_circleVpath(vpathcircle_window{i}, []);
    circlepow_window(i,:,:) = squeeze(nanmean(abs(tmp).^2));
    circlespctrm_window(i,:,:) = squeeze(nanmean(tmp));
end

figure; plot(circlefreq, nanmean(circlepow,3)), hold on, plot(circlefreq, nanmean(nanmean(circlepow,3)), 'k', 'LineWidth', 2) 
title('Circle vpath PSD per subject, CamCan Rest');
plot4paper('Frequenzy (Hz)', 'PSD');
print([config.figdir,'1F1_CircleFreq'],'-dpng');

figure; plot(circlefreq, nanmean(circlepow_window,3)), hold on, plot(circlefreq, nanmean(nanmean(circlepow_window,3)), 'k', 'LineWidth', 2) 
title('Windowed circle vpath PSD per subject, CamCan Rest');
plot4paper('Frequenzy (Hz)', 'PSD');
print([config.figdir,'1F2_CircleFreq_windowed'],'-dpng');

