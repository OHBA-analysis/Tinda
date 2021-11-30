function config = getStudyDetails(whichstudy)
% this script returns the pathnames for all models that have been run
config = [];
if whichstudy==1
    % this is the model run on Higgins2020_Neuron (ie on the MEG UK
    % partnership data of 55 subjects)
    config.nSj = 55;
    config.hmmfolder = '/Users/chiggins/data/Neuron2020/CanonicalRS/250Hz/hmm_1to45hz/';
    config.hmmfilename = 'hmm5_parc_giles_symmetric__pcdim80_voxelwise_embed14_K12_big1_dyn_modelhmm.mat';
    config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
    config.sample_rate = 250;
    config.prepdatafile = [config.hmmfolder,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
elseif whichstudy==2
    % this is the model run on the same MEGUK partnership data with a 4Hz
    % HPF to eliminate any slow frequency waves that might be the sole
    % cause of the cyclical patterns observed
    config.nSj = 55;
    config.hmmfolder = '/Users/chiggins/data/Neuron2020/CanonicalRS/250Hz/hmm_1to45hz/';
    config.hmmfilename = 'hmm1highfreqonly_parc_giles_symmetric__pcdim80_voxelwise_embed14_K12_big1_dyn_modelhighfreqonly.mat';
    config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
    config.sample_rate = 250;
    config.prepdatafile = [config.hmmfolder,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
    
elseif whichstudy==3 
    % this is the best of five models run on HCP data
    config.nSj = 237/3;
    config.hmmfolder = '/Volumes/CamsHD2/HCP_CH/ve_output_rest/';
    config.hmmfilename = 'hmm_analysis2.mat';
    config.parc = parcellation('aal_cortical_merged_8mm_stacked.nii.gz');
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
    
elseif whichstudy==4 % this the CamCan model fit:
    config.nSj = 600;
    config.hmmfolder = '/Volumes/CamsHD2/CamCan_2021/HMM/';
    config.hmmfilename = 'hmm_analysis1.mat';
    config.parc = parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
    config.sample_rate = 250;
    %config.prepdatafile = [config.hmmfolder,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
    config.matfilelist = '/Volumes/CamsHD2/CamCan_2021/HMM/matfiles/filelist.mat';
    config.participantfile = '/Users/chiggins/Documents/CamCan/participants.tsv';
    config.secondlevelmodelfile = '/Volumes/CamsHD2/CamCan_2021/HMM/secondLevelHMM_Poiss_window17_K3.mat';
    config.Poiss_dir = '/Volumes/CamsHD2/CamCan_2021/HMM/Poissdata_125_overlappingWindows/';
elseif whichstudy==5 % this refers to the HCP task epochs
    config.nSj = 237/3;
    config.hmmfolder = '/Volumes/CamsHD2/HCP_CH/ve_output_rest/';
    config.hmmfilename = 'hmm_analysis2.mat';
    config.parc = parcellation('aal_cortical_merged_8mm_stacked.nii.gz');
    config.wrkmemdir = '/Volumes/CamsHD2/HCP_CH/ve_output_Wrkmem_matfiles/';
    config.wrkmemfilelist = '/Volumes/CamsHD2/HCP_CH/ve_output_Wrkmem_matfiles/filelist.mat';
    config.storymdir = '/Volumes/CamsHD2/HCP_CH/ve_output_StoryM_matfiles/';
    config.storymfilelist = '/Volumes/CamsHD2/HCP_CH/ve_output_StoryM_matfiles/filelist.mat';
end
% generic output directories:
config.figdir = ['/Users/chiggins/Google Drive/Doctoral Files/3.0 Experimental Work/8.0SequentialPatterns/Study',int2str(whichstudy),'/'];
if ~isfolder(config.figdir)
    mkdir(config.figdir)
end
config.metricfile = [config.hmmfolder,'HMMsummarymetrics.mat'];
end