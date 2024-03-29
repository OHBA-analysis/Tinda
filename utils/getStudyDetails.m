function config = getStudyDetails(whichstudy)
% this script returns the pathnames for all models that have been run
config = [];

% generic output directories:
if isfolder('/Volumes/T5_OHBA/')
  config.figdir = ['/Volumes/T5_OHBA/Projects/Tinda/Study',int2str(whichstudy),'/'];
else
  config.figdir = ['/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study',int2str(whichstudy),'/figures/'];
  config.basedir = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/';
end
  config.resultsdir = ['/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study',int2str(whichstudy),'/'];


if whichstudy==1
  % this is the model run on Higgins2020_Neuron (ie on the MEG UK
  % partnership data of 55 subjects)
  config.nSj = 55;
  config.hmmfolder = '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/CanonicalRS/250Hz/hmm_1to45hz/';
  config.hmmfilename = 'hmm5_parc_giles_symmetric__pcdim80_voxelwise_embed14_K12_big1_dyn_modelhmm.mat';
  config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
  config.sample_rate = 250;
  config.prepdatafile = [config.hmmfolder,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
  config.metricfile = [config.resultsdir, 'HMMsummarymetrics.mat'];
  config.reordering_states = 'coherence'; % can be 'replay' (originally)
elseif whichstudy==2
  % this is the model run on the same MEGUK partnership data with a 4Hz
  % HPF to eliminate any slow frequency waves that might be the sole
  % cause of the cyclical patterns observed
  config.nSj = 55;
  config.hmmfolder = '/ohba/pi/mwoolrich/datasets/Neuron2020Analysis/CanonicalRS/250Hz/hmm_1to45hz/';
  config.hmmfilename = 'hmm1highfreqonly_parc_giles_symmetric__pcdim80_voxelwise_embed14_K12_big1_dyn_modelhighfreqonly.mat';
  config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
  config.sample_rate = 250;
  config.prepdatafile = [config.hmmfolder,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
  config.metricfile = [config.hmmfolder,'HMMsummarymetrics.mat'];
  config.reordering_states = 'coherence'; % can be 'replay' (originally), 'study1matched'
  
elseif whichstudy==3
  % this is the best of five models run on HCP data
  config.nSj = 237/3;
  config.nSes=3;
  config.hmmfolder = '/ohba/pi/mwoolrich/datasets/HCP_CH_2022/ve_output_rest/';
  config.hmmfilename = 'hmm_analysis2.mat';
  config.participantcovariates = '/ohba/pi/mwoolrich/datasets/HCP_CH_2022/HCPAnalysis/behav/';
  config.parc = parcellation('aal_cortical_merged_8mm_stacked.nii.gz');
  config.fmri_metastates = '/ohba/pi/mwoolrich/mvanes/analysis/HCP/Diegov/';
  config.metricfile = [config.resultsdir,'HMMsummarymetrics.mat'];
  % note the parcellation needs to be reordered so labels match data:
  %     reordering = readtable('/Users/chiggins/Documents/MATLAB/rs-mvar-scripts/aal_cortical_details.csv');
  %     reordering = table2array(reordering(:,5))+1; % the +1 converts from python 0 indexing
  %     reordervals = reordering;
  %     for i=1:78
  %         reordervals(i) = find(reordering==i);
  %     end
  %     newweightmask = config.parc.weight_mask;
  %     for i=1:78
  %         newweightmask(:,:,:,(i)) = config.parc.weight_mask(:,:,:,reordervals(i));
  %     end
  %     config.parc.weight_mask = newweightmask;
  config.sample_rate = 240;
  %    config.Poiss_dir
  config.reordering_states = 'study1matched'; %
  
elseif whichstudy==4 % this the CamCan model fit:
    basedir='/ohba/pi/mwoolrich/datasets/CamCan_2021';
  config.nSj = 600;
  config.hmmfolder = fullfile(basedir, 'HMM/');
  config.hmmfilename = 'hmm_analysis1.mat';
  config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
  config.sample_rate = 250;
  %config.prepdatafile = [config.hmmfolder,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
  config.matfilelist = fullfile(basedir, 'HMM/matfiles/filelist.mat');
  config.participantcovariates = fullfile(basedir, 'ParticipantCovariates/');
  config.participantfile = fullfile(basedir, 'ParticipantCovariates/participants.tsv');
  config.secondlevelmodelfile = fullfile(basedir, 'HMM/secondLevelHMM_Poiss_window17_K3.mat');
  config.Poiss_dir = fullfile(basedir, 'HMM/Poissdata_125_overlappingWindows/');
  config.metricfile = [config.resultsdir,'HMMsummarymetrics.mat'];
  config.reordering_states = 'coherence'; 
elseif whichstudy==5
  config.nSj = 237/3;
  config.hmmfolder = '/ohba/pi/mwoolrich/datasets/HCP_CH_2022/ve_output_rest/';
  config.hmmfilename = 'hmm_analysis2.mat';
  config.parc = parcellation('aal_cortical_merged_8mm_stacked.nii.gz');
  config.participantcovariates = '/ohba/pi/mwoolrich/datasets/HCP_CH_2022/HCPAnalysis/behav/';
  config.wrkmemdir = '/ohba/pi/mwoolrich/datasets/HCP_CH_2022/ve_output_Wrkmem_matfiles/';
  config.wrkmemfilelist = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study5/ve_output_Wrkmem_matfiles/filelist.mat';
  config.storymdir = '/ohba/pi/mwoolrich/datasets/HCP_CH_2022/ve_output_StoryM_matfiles/';
  config.storymfilelist = '/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study5/ve_output_StoryM_matfiles/filelist.mat';
  config.sample_rate=240;
  config.metricfile = [config.resultsdir, 'HMMsummarymetrics'];%[config.hmmfolder,'HMMsummarymetrics.mat'];
  config.reordering_states = 'study1matched'; % can be 'replay' (originally)
elseif whichstudy==6 % this is Chet's CamCan model fit:
  basedir='/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study6/';
  config.resultsdir=basedir;
  config.figdir = [config.resultsdir, 'figures/'];
  config.nSj = 612;
  config.hmmfolder = basedir;
  config.hmmfilename = [];%'hmm_analysis.mat';
  config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
  config.sample_rate = 250;
  %config.prepdatafile = [config.hmmfolder,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
  config.matfilelist = fullfile(basedir, 'Poissdata_125_overlappingWindows/filelist.mat');
  config.participantcovariates = basedir;
  config.participantfile = fullfile(basedir, 'camcan_all.csv');
  config.secondlevelmodelfile = fullfile(basedir, 'secondLevelHMM_stoch_Poiss_window125_K3_overlappingWindows.mat');
  config.Poiss_dir = fullfile(basedir, 'Poissdata_16_overlappingWindows/');
  config.metricfile = [config.resultsdir,'HMMsummarymetrics.mat'];
  config.reordering_states = 'study1matched'; 
elseif whichstudy==7 % Python WakeHen fit
    basedir='/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study7/';
  config.nSj = 19;
  config.nSes = 6;
  config.sample_rate = 250;
  config.hmmfolder = basedir;
  config.hmmfilename = [];
  config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
  config.sample_rate = 250;
  config.prepdatafile = [];
  config.metricfile = [config.hmmfolder,'HMMsummarymetrics.mat'];
  config.reordering_states = 'study1matched'; % can be 'replay' (originally), 'study1matched'
  config.Poiss_dir = fullfile(basedir, 'Poissdata_16_overlappingWindows/');
  elseif whichstudy==8 % Python CamCan fit, refitted on WakeHen 
    basedir='/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study8/';
  config.nSj = 19;
  config.nSes = 6;
  config.sample_rate = 250;
  config.hmmfolder = basedir;
  config.hmmfilename = [];
  config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
  config.sample_rate = 250;
  config.prepdatafile = [];
  config.metricfile = [config.hmmfolder,'HMMsummarymetrics.mat'];
  config.reordering_states = 'study1matched'; % can be 'replay' (originally), 'study1matched'
end

if isfolder('/Applications')
  config.workbenchdir = '/Applications/workbench/bin_macosx64';
else
  config.workbenchdir = '/ohba/pi/mwoolrich/mvanes/software/workbench/bin_linux64';
end
if ~isfolder(config.figdir)
  mkdir(config.figdir)
end

config.metricfile_gauss = [config.hmmfolder,'HMM_gauss_summarymetrics.mat'];
end