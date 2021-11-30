
% to run on maxfiltered, movement corrected data:
%configfile = '/home/chiggins/Documents/MATLAB/camcan_full_movecomptransdef.yaml'

% to run on maxfiltered, movement corrected (but not translation data:
configfile = '/home/chiggins/Documents/MATLAB/camcan_full_movecomp.yaml'

% to run on maxfiltered, no movecomp data
%configfile = '/home/chiggins/Documents/MATLAB/camcan_full_nomovecomp.yaml'


sessions_to_run = [1:639];

camcan_preproc_CHedit(configfile,sessions_to_run)

%% now more manual code to implement sign flipping and orthogonalisation and save to a dedicated folder:

% this to be run on completion of run_camcan preproc script; sign flips and
% orthogonalises data for use with eg TIDE models

MEG_RS_preprocceddir = '/ohba/pi/mwoolrich/datasets/CamCan_2021/MEGRS_preprocessed'
outputdir = '/ohba/pi/mwoolrich/datasets/CamCan_2021/MEGRS_signflippedOrthogonalised/';
if ~isdir(outputdir)
    mkdir(outputdir)
end
if ~exist([outputdir,'filelist.mat'])
    filelist = dir([MEG_RS_preprocceddir,'/dmf*.mat']);
    for i=1:length(filelist)
        fnames{i} = fullfile(MEG_RS_preprocceddir,filelist(i).name);
        D = spm_eeg_load(fnames{i});
        temp = montage(D,'getnumber');
        preproccompleted(i) = temp==5;
        fprintf(['\n subj',int2str(i)])
    end
    fnames = fnames(preproccompleted);


    addpath(genpath('/home/chiggins/Documents/MATLAB/Higgins2020_Neuron/'))

    % first, establish a good template subject:
    S = [];
    S.concat= [];
    S.concat.protocol = 'none';
    S.concat.embed.do = 1;
    S.concat.embed.num_embeddings = 1;
    S.concat.embed.rectify = false;
    S.concat.whiten = 1;
    S.concat.normalisation='voxelwise';
    statenetmats_cov_preflipped = hmm_full_global_cov(fnames,S)
    
    % find optimal template subject:
    modes={'none','abs'};
    diag_offset=15;
    metric_global=zeros(length(statenetmats_cov_preflipped),length(statenetmats_cov_preflipped),length(modes));

    for mm=1:length(modes)
        for subj=1:length(statenetmats_cov_preflipped)
           for subj2=1:length(statenetmats_cov_preflipped)
                if subj2~=subj
                    metric_global(subj,subj2,mm)=matrix_distance_metric(statenetmats_cov_preflipped{subj}.global.netmat_full, statenetmats_cov_preflipped{subj2}.global.netmat_full,diag_offset,modes{mm},[]);
                end
           end
        end
    end

    tmp=sum(metric_global(:,:,2),2);

    [~, template_subj]=max(tmp);
    save([outputdir,'filelist.mat'],'template_subj','fnames');
else
    load([outputdir,'filelist.mat'],'template_subj','fnames');
end

%% now start signflipping in batches:
S=[];
S.roinets_protocol='innovations_mar';
S.innovations_mar_order = 3;            
S.num_iters=1000;
S.prefix='sfold_';
S.num_embeddings = 4;
S.outputdir = outputdir; % directory to save signflipped files to
nbatch = 20;
tic
for ibatch = 18:ceil(length(fnames)/nbatch)
    subjs_to_do = [1:nbatch] + (ibatch-1)*nbatch;
    if any(subjs_to_do)>length(fnames)
        subjs_to_do = subjs_to_do(1):length(fnames);
    end
    if ~ismember(template_subj,subjs_to_do)
        S.Ds = {fnames{subjs_to_do},fnames{template_subj}};
        S.subj_template=length(subjs_to_do)+1;
             
    else
        S.Ds = {fnames{subjs_to_do}};
        S.subj_template=find(subjs_to_do==template_subj);
    end
    [ signflipped_files_out{ibatch}, sign_flip_res{ibatch} ] = find_sign_flips( S );
end
t = toc

sign_flip_results.signflipped_files=signflipped_files_out;
sign_flip_results.results = sign_flip_res;
save([outputdir 'sign_flip_results' ],'-struct','sign_flip_results','-v7.3');

%% check and reload filenames:

signflippedfiles = dir(outputdir);
for i=1:length(signflippedfiles)
    spmfiles(i) = endsWith(signflippedfiles(i).name,'.dat');
end
spmfiles = find(spmfiles);
for i=1:length(spmfiles)
    SFfiles{i} = fullfile(outputdir,signflippedfiles(spmfiles(i)).name);
end


%% Orthogonalise to correct source leakage:
for i=1:length(SFfiles)
    fprintf(['\n Orthogonalising subject ',int2str(i)])
    D = spm_eeg_load(SFfiles{i});
    X = D(:,:);
    X = ROInets.remove_source_leakage(X,'symmetric');
    X = normalise(X,2);
    D(:,:,1) = X;
    D.save();
end