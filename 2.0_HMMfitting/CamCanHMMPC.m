% CamCam_HMM
if isfolder('/Volumes/T5_OHBA/')
  basedir = '/Volumes/T5_OHBA/Projects/Tinda';
else
  basedir = '/ohba/pi/mwoolrich/datasets/CamCan_2021';
end
camcandir_matfiles = fullfile(basedir, 'HMM/matfiles/');
if ~exist([camcandir_matfiles,'filelist.mat'])
    % save as matfiles:
    camcandir = fullfile(basedir,'MEGRS_signflippedOrthogonalised');
    files = dir(camcandir);
    clear to_keep;
    for i=1:length(files)
        to_keep(i) = endsWith(files(i).name,'meg.mat');
    end
    to_keep = find(to_keep);
    clear spm_files mat_files_orth
    for i=1:length(to_keep)
        spm_files{i} = [camcandir,files(to_keep(i)).name];
    end
    if ~isdir(camcandir_matfiles)
        mkdir(camcandir_matfiles);
    end
    clear T_all;
    for i=1:length(spm_files)
        mat_files_orth{i} = [camcandir_matfiles,files(to_keep(i)).name];
        D = spm_eeg_load(spm_files{i});
        gs = good_samples(D);
        X = D(:,:)';
        X = X(gs,:);
        Ton = diff([0,gs,0])==1;
        Toff = diff([0,gs,0])==-1;
        T = find(Toff)-find(Ton);
        T_all{i} = T;
        if sum(T)~=size(X,1)
            error('Size mismatch!');
        end
        save(mat_files_orth{i},'X','T');
        fprintf(['\nSj ',int2str(i),' done.']);
    end

    % and save overall info:
    save([camcandir_matfiles,'filelist'],'mat_files_orth','T_all');
else
    load([camcandir_matfiles,'filelist'],'mat_files_orth','T_all');
    for k=1:length(mat_files_orth)
      mat_files_orth{k} = replace(mat_files_orth{k}, [basedir, '/HMM/matfiles/'], camcandir_matfiles);
    end
end
%% and run HMM:
dirnameHMM = fullfile(basedir, 'HMM/');
% for i=1:length(mat_files_orth)
%     mat_files_orth{i}(1)='E';
% end
Nruns = 1;
embeddedlag = 7; K = 12; Hz = 250; ndim = 38;

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
options.BIGdelay = 1;
options.BIGforgetrate = 0.9;
options.BIGbase_weights = 0.9;
options.BIGcomputeGamma = false; % runs out of memory otherwise

for irun=1:Nruns

    outputfile = [dirnameHMM 'hmm_analysis',int2str(irun),'.mat'];
    fprintf(['Running HMM:']);
    [hmm, ~, ~, ~,~,~,fehist] = hmmmar(mat_files_orth,T_all,options);
    %save(outputfile,'hmm','Gamma','vpath')
%     hmm.gamma = Gamma;
    disttoplot = plotMDS_states(hmm);
    [~,new_state_ordering] = sort(disttoplot(:,1));
    new_state_ordering = flipud(new_state_ordering);
    if any(new_state_ordering(2:end) ~= new_state_ordering(1:end-1)+1)
        hmm = hmm_permutestates(hmm,new_state_ordering);
%         disttoplot = disttoplot(new_state_ordering,:);
%         vpath2 = zeros(size(vpath));
%         for i=1:hmm.K
%             vpath2(vpath==new_state_ordering(i))=i;
%         end
     end
%     Gamma = hmm.gamma;
%     vpath = vpath2;
%     hmm = rmfield(hmm,'gamma');
    save(outputfile,'hmm','T_all','fehist')

    % fit to each data file iteratively:
    for i=1:length(mat_files_orth)
        fprintf(['\nDecoding subj: ',int2str(i)]);
        vpath = hmmdecode(mat_files_orth(i),T_all(i),hmm,1);
        save(mat_files_orth{i},'vpath','-append');
    end
end

