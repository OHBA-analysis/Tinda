%% script to run analyses related to figure 2 of Neuron2020 Paper
for whichstudy=1:2
K=12;

%% general setup for this study data:
bestmodel=5;
wd = '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/';
if whichstudy+1==1
    CanonicalRSN = false;
else
    CanonicalRSN = true; % canonical refers to RSNs trained on nottingham data and refitted
end
session_name{1} = 'CanonicalRS/';
session_name{2} = 'Study1/';
session_name{3} ='Study2/';


% Define colours to use in state plots
color_scheme = set1_cols();

% Define sample rate
sample_rate = 250;
    
if CanonicalRSN
    template_string = [int2str(bestmodel),'usingtemplate'];
else
    template_string = [int2str(bestmodel),''];
end

studydir = [wd,session_name{whichstudy+1}];


hmmdir = [studydir,'hmm_1to45hz/'];
hmmfile = [hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed13_K',int2str(K),'_big1_dyn_modelhmm.mat'];
load(hmmfile);    % no need to permute here as will do below 
prepdatafile = [hmmdir,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed13.mat'];
load(prepdatafile,'hmmT','subj_inds');
scan_T = cell2mat(hmmT);

parc_file = ['fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
parc = parcellation(parc_file);
for k=1:hmm.K;statelabels{k}={'RSN-State ',int2str(k)};end

%% SECTION 1: Analyse Transition matrix and order states:
    
% Plot MDS map of state network:
disttoplot = plotMDS_states(hmm);
[~,new_state_ordering] = sort(disttoplot(:,1));
load('/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study1/coherence_state_ordering')
if any(new_state_ordering(2:end) ~= new_state_ordering(1:end-1)+1)
    hmm = hmm_permutestates(hmm,new_state_ordering);
    hmmfile = [hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed13_K',int2str(K),'_big1_dyn_modelhmm.mat'];
    save(hmmfile,'new_state_ordering','-append');
    disttoplot = disttoplot(new_state_ordering,:);
end
Gamma = hmm.gamma;

%% SECTION 2: STATE REPLAY ALIGNMENT
if whichstudy+1 > 1
    datadir = [studydir,'bfnew_1to45hz/'];
    
    % load in alignment data:
    [~,~,triggerpoints,goodsamples] = getSubjectMasks(datadir);
    
    % load replay times:
    replayScores=[];
    if whichstudy+1==2
        load([wd,'GenericReplayData/STUDYI_ReplayOnset/StrReplayOnset'],'ToRall'); % the replay scores for first resting state
        replayScores(1:2:(2*21),:) = ToRall;
        load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
        replayScores(2:2:(2*21),:) = ToRall;
    else
        load([wd,'GenericReplayData/STUDYII_ReplayOnset/STUDYII_PreplayOnset'],'ToRall'); % the replay scores for first resting state
        replayScores(1:2:(2*22),:) = ToRall;
        load([wd,'GenericReplayData/STUDYII_ReplayOnset/STUDYII_ReplayOnset'],'ToRall');
        replayScores(2:2:(2*22),:) = ToRall;
    end
   
    pthreshold = 0.05/hmm.K; %display threshold
    emphasis_states = 1:6;
    emphasis_string= 'SigEmphasis';
    
    t_window = sample_rate/2;
    [~,betas_replay,betas_replay_norm] = plotReplayFig1(Gamma,replayScores,t_window,goodsamples,triggerpoints,scan_T,pthreshold,1,color_scheme,emphasis_states);
    YL = ylim();
    if whichstudy+1==2;YL(2) = YL(2)+0.02;ylim(YL);end
    line([0,0],YL,'LineWidth',2,'Color','red','LineStyle','--');
end
rep_scores{whichstudy} = replayScores;
betas{whichstudy} = betas_replay;
betas_norm{whichstudy}=betas_replay_norm;

  scan_T = cell2mat(hmmT);
  scan_T_studies{whichstudy}=scan_T;
  R = [[1;1+cumsum(scan_T(1:end-1))'],cumsum(scan_T(1:end))'];
  if ~all(R(2:end,1)==R(1:end-1,2)+1)
    % correct any uncorrected offsets in R:
    R(2:end,1) = R(1:end-1,2)+1;
  end
nSes=length(hmm.data_files);
  for iSes=1:nSes
    iSes
    % convert Gamma to padded timecourse
    Gam_long = NaN(length(goodsamples{iSes}),K);
    Gam_long(goodsamples{iSes},:) = Gamma(R(iSes,1):R(iSes,2),:);
    Gam_long = Gam_long(triggerpoints{iSes},:);
    
    G{whichstudy}{iSes} = Gam_long - nanmean(Gam_long);
    
    temp1 = zeros(length(goodsamples{iSes}),1);
    temp1(goodsamples{iSes}) = hmm.statepath(R(iSes,1):R(iSes,2));
    temp1 = temp1(triggerpoints{iSes},:);
    vpath_studies{whichstudy}{iSes} = temp1;
    T_studies{whichstudy}{iSes} = [length(vpath_studies{whichstudy}{iSes})];
  end
end