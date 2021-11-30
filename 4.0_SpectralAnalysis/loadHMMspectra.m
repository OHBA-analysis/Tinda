function [P,f] = loadHMMspectra(config,whichstudy,hmm,run_inds,FO_sj)
% this function loads the subject-average state spectra associaed with the
% model fit details specified by the config file.
waveletfolder = [config.hmmfolder,'WTspect/'];
averagespectfile = [waveletfolder,'averagespect.mat'];
if isfile(averagespectfile)
    load(averagespectfile)
    fprintf(['Average state spectra loaded from file: ',averagespectfile]);
else
    if ~isdir(waveletfolder)
        mkdir(waveletfolder)
    end
    fprintf(['Average state spectra not found; recomputing from subject files now...']);
    % compute:
    if whichstudy==1
        % note that waveletoption has the following correspondence:
        %   1   should not be used (had an error in computation)
        %   2   is weighted wavelet psd
        %   3   is identical to above
        %   4   is computed with hamming window gamma weighting
        waveletoptions = 3; % note here that 1 is wrong; 2 is 
        tf_morlet_factor = 3;
        wtfilename = [config.hmmfolder,'hmm5_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(hmm.K),'_big1_dyn_modelhmm_store/',...
            'state_netmats_WT3_',int2str(tf_morlet_factor),'sess_2_vn0_soft_global0.mat'];
        wt = load(wtfilename);
        % note these need to be reordered:
        load([config.hmmfolder,config.hmmfilename],'new_state_ordering');
        psd_temp = wt.psd;
        for k=1:length(new_state_ordering)
            wt.psd(:,k,:,:,:) = psd_temp(:,new_state_ordering(k),:,:,:);
        end
        clear psd_temp;
        % and compute:
        P = zeros(60,config.parc.n_parcels,config.parc.n_parcels,hmm.K);
        Pweighted = zeros(60,config.parc.n_parcels,config.parc.n_parcels,hmm.K);
        FO_sj = zeros(config.nSj,12);
        for i=1:config.nSj
            % load that subject's fractional occupancy:
            FO_sj(i,:) = mean(hmm.gamma(run_inds==i,:));
            for k=1:hmm.K
                Pweighted(:,:,:,k) = Pweighted(:,:,:,k) + FO_sj(i,k)*permute(wt.psd(i,k,:,:,:),[3,4,5,1,2]);
                P(:,:,:,k) = P(:,:,:,k) + permute(wt.psd(i,k,:,:,:),[3,4,5,1,2]);
            end
        end
        Pweighted = permute(Pweighted,[4,1,2,3]);
        P = permute(P,[4,1,2,3]);
        for k=1:12
            Pweighted(k,:,:,:) = Pweighted(k,:,:,:) ./ sum(FO_sj(:,k));
        end
        P = P./i;
        f = wt.f;
    elseif whichstudy==3 || whichstudy ==4
        if whichstudy==3
            nF = 88;
            nSj = config.nSj*3;
            nch = config.parc.n_parcels;
        else
            nF = 88;
            nSj = 300;
            nch = config.parc.n_parcels;
        end
        P = zeros(nF,nch,nch,12);
        Pweighted = zeros(nF,nch,nch,12);
        FO_sj = zeros(nSj,12);
        for i=1:nSj % load all subject's runs and average:
            if mod(i,20)==0
                fprintf(['\nLoading spectra for subj',int2str(i),' of ',int2str(nSj)]);
            end
            load([waveletfolder,'hmm_spectrawt_sj',int2str(i),'.mat']);
            % load that subject's fractional occupancy:
            if nargin<5
                FO_sj(i,:) = mean(hmm.gamma(run_inds==i,:));
            end
            for k=1:12
                Pweighted(:,:,:,k) = Pweighted(:,:,:,k) + FO_sj(i,k)*fitmt_subj.state(k).psd;
                P(:,:,:,k) = P(:,:,:,k) + fitmt_subj.state(k).psd;
            end
        end
        Pweighted = permute(Pweighted,[4,1,2,3]);
        P = permute(P,[4,1,2,3]);
        for k=1:12
            Pweighted(k,:,:,:) = Pweighted(k,:,:,:) ./ sum(FO_sj(:,k));
        end
        P = P./i;
        f = fitmt_subj.state(1).f;
    end
    save(averagespectfile,'P','Pweighted','f');
end

end