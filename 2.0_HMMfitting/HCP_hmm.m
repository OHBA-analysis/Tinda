% Note: on shared drive this analysis is saved in directory HCP_CH_2022


% this the script to run on the October 2021 preprocced data:

addpath(genpath('C:\Users\chiggins\Google Drive\MATLAB\8.0 Sequential RSNs'));

% dirname = 'C:\Users\chiggins\My Data\HCP_sails\ve-preproc\';
dirname = 'G:\HCP_CH\ve_output_rest\'
dirname_matfiles = 'G:\HCP_CH\ve_output_rest_matfiles\'
if ~isdir(dirname_matfiles)
    mkdir(dirname_matfiles)
end
filenames = dir(dirname);
hdf5files = [];
for i=1:length(filenames)
    if endsWith(filenames(i).name,'hdf5')
        hdf5files = [hdf5files,i];
    end
end
filenames = filenames(hdf5files);
clear fnamepart T_all
%% orthogonalise data and store as mat files:
clear mat_files T_all
for i=1:length(filenames)
    fnamepart{i} = filenames(i).name;
    fname = [dirname,fnamepart{i}];
    X = h5read(fname,'/data');
    X = ROInets.remove_source_leakage(X', 'symmetric');
    X = normalise(X);
    T = length(X);
    T_all{i} = length(X);
    %X = [X;Xtemp];
    SR(i) = h5read(fname,'/sample_rate');
    fnamemat = [dirname_matfiles,fnamepart{i}(1:end-8),'.mat'];
    %if ~isfile(fnamemat)
        save(fnamemat,'X','T');
    %end
    mat_files{i} = fnamemat;
end
SR = unique(SR);
if length(SR)>1
    error('Data not at same sample rate');
end

% sign flipping without sourceleakage correction:

options_signflip = [];
options_signflip.maxlag = 0; 
options_signflip.verbose = 1;
options_signflip.nbatch = 10;
tic;
flips = findflip(mat_files,T_all(1:50),options_signflip);
%flipdata(mat_files,T_all,flips);
timerec(1) = toc;
%%
for i=1:length(filenames)
    fnamepart{i} = filenames(i).name;
    fname = [dirname,fnamepart{i}];
    X = h5read(fname,'/data');
    X = ROInets.remove_source_leakage(X', 'symmetric');
    X = normalise(X');
    T = length(X);
    T_all{i} = length(X);
    %X = [X;Xtemp];
    SR(i) = h5read(fname,'/sample_rate');
    fnamemat = [dirname_matfiles,fnamepart{i}(1:end-8),'_orth.mat'];
    %if ~isfile(fnamemat)
        save(fnamemat,'X','T');
    %end
    mat_files_orth{i} = fnamemat;
end
SR = unique(SR);
if length(SR)>1
    error('Data not at same sample rate');
end

% sign flipping:
options_signflip = [];
options_signflip.maxlag = 3; 
options_signflip.verbose = 1;
options_signflip.nbatch = 10;
options_signflip.maxcyc = 500;
nSj = length(filenames);
nbatch = 15;
for ibatch=1:floor(nSj/nbatch)
    if ibatch==1
        n_subj_this_batch = [1:nbatch];
    else
        n_subj_this_batch = [1,[1:nbatch]+(ibatch-1)*nbatch];
        if any(n_subj_this_batch>nSj)
            n_subj_this_batch(n_subj_this_batch>nSj) = [];
        end
    end
    matfiles_thisbatch = mat_files_orth(n_subj_this_batch);
    T_this_batch = T_all(n_subj_this_batch);
    flipstemp = findflip(matfiles_thisbatch,T_this_batch,options_signflip);
    % set so first row identical, and not flipped:
    if flipstemp(1,ich)
        flipstemp(:,ich) = ~flipstemp(:,ich);
    end
    flipdata(matfiles_thisbatch(2:end),T_this_batch(2:end),flipstemp(2:end,:),[],1);
end

%% Fit the HMM:

Nruns = 5;
embeddedlag = 7; K = 12; Hz = SR; ndim = 78;

options = struct();
options.order = 0;
options.zeromean = 1;
options.covtype = 'full';
options.embeddedlags = -embeddedlag:embeddedlag;
options.pca = ndim*2;
options.K = K;
options.Fs = Hz;
options.verbose = 1;
options.onpower = 0; 
options.standardise = 1;
options.standardise_pc = true;
options.inittype = 'HMM-MAR';
options.cyc = 100;
options.initcyc = 10;
options.initrep = 3;

% stochastic options
options.BIGNinitbatch = 5;
options.BIGNbatch = 5;
options.BIGtol = 1e-7;
options.BIGcyc = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay = 5;
options.BIGforgetrate = 0.7;
options.BIGbase_weights = 0.9;
for irun=1:Nruns

    outputfile = [dirname 'hmm_analysis',int2str(irun),'.mat'];
    fprintf(['Running HMM:']);
    [hmm, Gamma, ~, vpath,~,~,fehist] = hmmmar(mat_files_orth,T_all,options);
    %save(outputfile,'hmm','Gamma','vpath')
    hmm.gamma = Gamma;
    disttoplot = plotMDS_states(hmm);
    [~,new_state_ordering] = sort(disttoplot(:,1));
    new_state_ordering = flipud(new_state_ordering);
    if any(new_state_ordering(2:end) ~= new_state_ordering(1:end-1)+1)
        hmm = hmm_permutestates(hmm,new_state_ordering);
        %hmmfile = [hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
        %save(hmmfile,'new_state_ordering','-append');
        disttoplot = disttoplot(new_state_ordering,:);
        vpath2 = zeros(size(vpath));
        for i=1:hmm.K
            vpath2(vpath==new_state_ordering(i))=i;
        end
    end
    Gamma = hmm.gamma;
    vpath = vpath2;
    hmm = rmfield(hmm,'gamma');
    save(outputfile,'hmm','Gamma','vpath','T_all','fehist')


end

%% find lowest free energy model:
for irun=1:Nruns

    outputfile = [dirname 'hmm_analysis',int2str(irun),'.mat'];
    load(outputfile,'fehist');
    free_en_rec(irun) = fehist(end);
end
[~,best_model] = min(free_en_rec)


%% spectral estimation using wavelets:

outputfile = [dirname 'hmm_analysis',int2str(best_model),'.mat'];
load(outputfile);

dirname_wt = [dirname,'WTspect\']
mkdir(dirname_wt);

N = length(mat_files_orth);
options_mt = struct('Fs',SR); % Sampling rate - for the 25subj it is 300
options_mt.fpass = [1 45];  % band of frequency you're interested in
options_mt.p = 0; %0.01; % interval of confidence  
options_mt.to_do = [1 0]; % turn off pdc
options_mt.order = 0;
options_mt.embeddedlags = -7:7;

% average - NO takes too long to load all data; can just average over
% subject PSDs below
% fitmt = hmmspectramt(mat_files,T_all,Gamma,options_mt);

% per subject
%fitmt_subj = cell(N,1);
d = length(options_mt.embeddedlags) - 1; 
acc = 0; for n=1:N
    fprintf(['\n RUNNING WAVELET ON SUBJECT: ',int2str(n),'\n']);
    load(mat_files_orth{n});
    gamma = Gamma(acc + (1:(sum(T_all{n})-length(T_all{n})*d)),:);
    acc = acc + size(gamma,1);
    fitmt_subj = hmmspectrawt(X,T_all{n},gamma,options_mt);
    fitmt_subj.state = rmfield(fitmt_subj.state,'ipsd');
    fitmt_subj.state = rmfield(fitmt_subj.state,'pcoh');
    fitmt_subj.state = rmfield(fitmt_subj.state,'phase');
    disp(['Subject ' num2str(n)])
    outputfilemt = [dirname_wt 'hmm_spectrawt_sj', int2str(n),'.mat'];
    save(outputfilemt,'fitmt_subj')
end

%%
psd_m = zeros(88,78,12);
inds = find(triu(ones(78),1));
coh_m = zeros(88,length(inds),12);
for n=1:N
    outputfilemt = [dirname_wt 'hmm_spectrawt_sj', int2str(n),'.mat'];
    load(outputfilemt,'fitmt_subj')
    for k=1:12
        psd_m(:,:,k) = psd_m(:,:,k) + abs(fitmt_subj.state(k).psd(:,find(eye(78))));
        coh_m(:,:,k) = coh_m(:,:,k) + fitmt_subj.state(k).coh(:,inds);
    end
end
f = fitmt_subj.state(1).f;

figure();
subplot(1,2,1);
plot(f,squeeze(mean(psd_m(:,:,1:6),2)),'LineWidth',2); hold on;
plot(f,squeeze(mean(psd_m(:,:,7:12),2)),'LineWidth',2,'LineStyle','--')
subplot(1,2,2);
plot(f,squeeze(mean(coh_m(:,:,1:6),2)),'LineWidth',2);hold on;
plot(f,squeeze(mean(coh_m(:,:,7:12),2)),'LineWidth',2,'LineStyle','--');hold on;
legend('1','2','3','4','5','6','7','8','9','10','11','12')

%% NOW: refit above hmm model to task data:
addpath(genpath('C:\Users\chiggins\Google Drive\MATLAB\8.0 Sequential RSNs'));

% dirname = 'C:\Users\chiggins\My Data\HCP_sails\ve-preproc\';
dirname = 'G:\HCP_CH\ve_output_Wrkmem\'
dirname_matfiles = 'G:\HCP_CH\ve_output_Wrkmem_matfiles\'
if ~isdir(dirname_matfiles)
    mkdir(dirname_matfiles)
end
filenames = dir(dirname);
hdf5files = [];
for i=1:length(filenames)
    if endsWith(filenames(i).name,'hdf5') && ~endsWith(filenames(i).name,'trialinfo.mathdf5')
        hdf5files = [hdf5files,i];
    end
end
filenames = filenames(hdf5files);

%%
clear fnamepart T_all
for i=1:length(filenames)
    fnamepart{i} = filenames(i).name;
    fname = [dirname,fnamepart{i}];
    X = h5read(fname,'/data');
    X = ROInets.remove_source_leakage(X, 'symmetric');
    X = normalise(X');
    %T = length(X);
    T = h5read(fname,'/T_downsampled');
    T_all{i} = T;
    %X = [X;Xtemp];
    SR(i) = h5read(fname,'/sample_rate');
    fnamemat = [dirname_matfiles,fnamepart{i}(1:end-8),'_orth.mat'];
    %if ~isfile(fnamemat)
    save(fnamemat,'X','T');
    %end
    mat_files_orth{i} = fnamemat;
end
SR = unique(SR);
if length(SR)>1
    error('Data not at same sample rate');
end

%% load and fit to data:

outputfile = ['G:\HCP_CH\ve_output_rest\' 'hmm_analysis',int2str(best_model),'.mat'];
load(outputfile,'hmm') % careful not to also overwrite T_all!

clear T_all;
for i=1:length(filenames)
    fnamemat = [dirname_matfiles,fnamepart{i}(1:end-8),'_orth.mat'];
    load(fnamemat,'T')
    T_all{i} = T;
end

for i=1:length(filenames)
    fprintf(['\nSession ',int2str(i),' of ',int2str(length(filenames))]);
    fname = [dirname,fnamepart{i}];
    T = h5read(fname,'/T_downsampled');
    Gamma_task = hmmdecode(mat_files_orth(i),T_all(i),hmm,0); % 0 for state time courses
    save(mat_files_orth{i},'Gamma_task','-append');
end

%% quick evoked response plot:
filenames = dir('G:\HCP_CH\ve_output_Wrkmem_matfiles');
mat_files_orth = {};
for i=1:length(filenames)
    if endsWith(filenames(i).name,'orth.mat');
        mat_files_orth = {mat_files_orth{:},[filenames(i).folder,'/',filenames(i).name]};
    end
end
% check all trials same length:
TIMtask = false(length(mat_files_orth),1);
%T_cat = [];
for i=1:length(mat_files_orth)
    if endsWith(mat_files_orth{i},'TIM_orth.mat')  % note 'TIM' refers to image onset, as oposed to TRESP which is onset of button press
        TIMtask(i) = true;
        %T_cat = [T_cat;unique(T_all{i})];
    end
end
%figure();bar(unique(T_cat))
timfiles = mat_files_orth(TIMtask);
%% actual evoked response analysis - quick and dirty, not correcting subject baselines:

load(timfiles{1},'T');
epochlength = unique(T)-14;
for i=1:length(timfiles)
    load(timfiles{i},'Gamma_task','T');
    ntr = length(T);
    epochs = mean(reshape(Gamma_task,epochlength,ntr,12),2);
    evokedresponse(:,:,i) = squeeze(epochs);
end
figure();
meanresponse = mean(evokedresponse,3);
ls = {'-','--'};
for k=1:12
    plot(meanresponse(:,k),'LineWidth',2,'LineStyle',ls{mod(k,2)+1});
    hold on;
end
ylim([0,0.3]);
legend('1','2','3','4','5','6','7','8','9','10','11','12');
%% and to Story / Math task:

% dirname = 'C:\Users\chiggins\My Data\HCP_sails\ve-preproc\';
dirname = 'G:\HCP_CH\ve_output_StoryM\'
dirname_matfiles = 'G:\HCP_CH\ve_output_StoryM_matfiles\'
if ~isdir(dirname_matfiles)
    mkdir(dirname_matfiles)
end
filenames = dir(dirname);
hdf5files = [];
for i=1:length(filenames)
    if endsWith(filenames(i).name,'hdf5') && ~endsWith(filenames(i).name,'trialinfo.mathdf5')
        hdf5files = [hdf5files,i];
    end
end
filenames = filenames(hdf5files);

clear fnamepart T_all
for i=1:length(filenames)
    fnamepart{i} = filenames(i).name;
    fname = [dirname,fnamepart{i}];
    X = h5read(fname,'/data');
    X = ROInets.remove_source_leakage(X, 'symmetric');
    X = normalise(X');
    %T = length(X);
    T = h5read(fname,'/T_downsampled');
    T_all{i} = T;
    %X = [X;Xtemp];
    SR(i) = h5read(fname,'/sample_rate');
    fnamemat = [dirname_matfiles,fnamepart{i}(1:end-8),'_orth.mat'];
    %if ~isfile(fnamemat)
    save(fnamemat,'X','T');
    %end
    mat_files_orth{i} = fnamemat;
end
SR = unique(SR);
if length(SR)>1
    error('Data not at same sample rate');
end

% load and fit to data:

outputfile = ['G:\HCP_CH\ve_output_rest\' 'hmm_analysis',int2str(best_model),'.mat'];
load(outputfile,'hmm') % careful not to also overwrite T_all!

clear T_all;
for i=1:length(filenames)
    fnamemat = [dirname_matfiles,fnamepart{i}(1:end-8),'_orth.mat'];
    load(fnamemat,'T')
    T_all{i} = T;
end
%%
failedsessions = true(length(filenames),1);
for i=1:length(filenames)
    fprintf(['\nSession ',int2str(i),' of ',int2str(length(filenames))]);
    try
        Gamma_task = hmmdecode(mat_files_orth(i),T_all(i),hmm,0); % 0 for state time courses
        save(mat_files_orth{i},'Gamma_task','-append');
        failedsessions(i) = false;
    catch
        failedsessions(i) = true;
    end
end
%% sanity check:
for i=1:length(filenames)
    temp = load(mat_files_orth{i})
	szes(i,1) = size(temp.X,1);
    szes(i,2) = sum(temp.T);
    szes(i,3) = sum(T_all{i});
end

%% quick and dirty epoching of sametime epochs:
%filenames = dir(
sess_to_epoch = zeros(length(mat_files_orth),1);
for i=1:length(mat_files_orth)
    if endsWith(mat_files_orth{i},'TEV_orth.mat')
        sess_to_epoch(i) = unique(T_all{i});
    end
end

epochs1 = find(sess_to_epoch==720);
epochlength = 720-14;
clear evokedresponse
for i=1:length(epochs1)
    load(mat_files_orth{epochs1(i)},'Gamma_task','T');
    ntr = length(T);
    epochs = mean(reshape(Gamma_task,epochlength,ntr,12),2);
    evokedresponse(:,:,i) = squeeze(epochs);
end
figure();
subplot(1,2,1);
meanresponse = mean(evokedresponse,3);
ls = {'-','--'};
for k=1:12
    plot(meanresponse(:,k),'LineWidth',2,'LineStyle',ls{1+(k>6)});
    hold on;
end
ylim([0,0.3]);
legend('1','2','3','4','5','6','7','8','9','10','11','12');

epochs1 = find(sess_to_epoch==1320);
epochlength = 1320-14;
clear evokedresponse
for i=1:length(epochs1)
    load(mat_files_orth{epochs1(i)},'Gamma_task','T');
    ntr = length(T);
    epochs = mean(reshape(Gamma_task,epochlength,ntr,12),2);
    evokedresponse(:,:,i) = squeeze(epochs);
end
subplot(1,2,2);
meanresponse = mean(evokedresponse,3);
ls = {'-','--'};
for k=1:12
    plot(meanresponse(:,k),'LineWidth',2,'LineStyle',ls{1+(k>6)});
    hold on;
end
ylim([0,0.3]);
legend('1','2','3','4','5','6','7','8','9','10','11','12');
suptitle('TEV epochs');
%%
sess_to_epoch = zeros(length(mat_files_orth),1);
for i=1:length(mat_files_orth)
    if endsWith(mat_files_orth{i},'TRESP_orth.mat')
        sess_to_epoch(i) = unique(T_all{i});
    end
end
epochs1 = find(sess_to_epoch==720);
epochlength = 720-14;
clear evokedresponse
for i=1:length(epochs1)
    load(mat_files_orth{epochs1(i)},'Gamma_task','T');
    ntr = length(T);
    epochs = mean(reshape(Gamma_task,epochlength,ntr,12),2);
    evokedresponse(:,:,i) = squeeze(epochs);
end
figure();
meanresponse = mean(evokedresponse,3);
ls = {'-','--'};
for k=1:12
    plot(meanresponse(:,k),'LineWidth',2,'LineStyle',ls{1+(k>6)});
    hold on;
end
ylim([0,0.3]);
legend('1','2','3','4','5','6','7','8','9','10','11','12');
title('TRESP epochs')